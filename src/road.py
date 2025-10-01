import torch
import FV_schemes as fv

class Road:

    b = None # Length of the road. Positive integer.
    dx = None # Lenth of each control-volume.
    L = None # The standard unit for road lenghts. 50m is standard.
    Vmax = None # Speed limit
    gamma = None
    scheme = None
    pad = None # Number of ghost cells on each side of the road
    limiter = None
    left = False
    right = False 
    queue_length = None
    N = None # Number of grid points a unit segment is divided into.
    scheme = None
    id = ""
    boundary_fnc = None #WHAT IS THIS
    initial = None
    CFL = 0.5
    max_dens = 1
    periodic = False
    right_flux = None
    left_flux = None
    left_pos = None
    right_pos = None


    def __init__(self, b, L, N, scheme, id, limiter, Vmax, initial, boundary_fnc, left_pos, right_pos):
        self.b = b 
        self.L = L 
        self.N = N
        self.scheme = scheme
        self.id = id
        self.limiter = limiter
        self.Vmax = Vmax
        self.initial = initial
        self.boundary_fnc = boundary_fnc
        self.left_pos = left_pos
        self.right_pos = right_pos


        # Schemes: 1 - Rusanov, 2 - Lax-Wendroff, 3 - SSP_RK (2.order in time and space), 4 - Euler (2. order in space, 1 order in time)

        match scheme:
            case 0 | 1 | 2:
                self.pad = 1
                self.scheme = scheme
            case 3 | 4:
                self.pad = 2
                self.scheme = scheme
            case _:
                self.pad = 2
                self.scheme = 3

        """
        There might be something wrong with the scaling here.
        THERE IS SOMETHING FISHY GOING ON HERE WITH THE LENGTHS OF STUFF AND SHIT


        self.N_full = b * (N + 2*self.pad)
        self.dx = torch.tensor(1 / self.N_full) #Dimensionless.
        j = torch.arange(0, self.N_full, 1)
        x = (j + 1/2) * self.dx
        self.rho = self.initial(x)
        """


        self.N_full = b*N + 2*self.pad
        self.dx = torch.tensor(L / N) #Dimensionless.
        j = torch.arange(0, b*N, 1)
        
        x = (j + 1/2) * self.dx

        self.rho = torch.zeros(self.N_full)
        self.rho[self.pad:-self.pad] = self.initial(x)
        self.rho[:self.pad] = self.rho[self.pad]
        self.rho[-self.pad:] = self.rho[-self.pad-1]

        #print("IC =", self.initial(x))

        # From here is fine...
        self.left_boundary = self.rho[:self.pad].clone()
        self.right_boundary = self.rho[-self.pad:].clone()

        self.left_flux = torch.tensor(0.0)
        self.right_flux = torch.tensor(0.0)

        self.calculate_gamma()

    def calculate_gamma(self):
        # L*b is the total length of the road
        # Vmax is initialized in km/h, but needs to be converted in m/s, hence divided by 3.6
        self.gamma = self.Vmax/3.6 

    """   
    def calculate_gamma(self):
        # L*b is the total length of the road
        # Vmax is initialized in km/h, but needs to be converted in m/s, hence divided by 3.6
        self.gamma = self.Vmax/(self.L * self.b * 3.6) 
    """

    def max_dt(self):
        #OBS: Here max_dens is used.....
        max_flux = torch.max(torch.abs(fv.d_flux(self.rho, self.gamma)))

        return self.CFL * self.dx / (self.max_dens * max_flux)


    def solve_internally(self, dt):
        # Solves one step of the FV methods. Everything needs to be vecotrized.
        # Does not solve for the obundaries

        match self.scheme:
            case 0:
                pass

            case 1:
                # Rusanov
                # How global is the dt variable? Refers to the dt of network. But that variable is outside this class.

                #Why is this not a class variable...? Because it can be computed by the rho function..... And we do not reuse it.
                F = fv.Rusanov_Flux(self.rho, self.gamma)
                self.rho[self.pad:-self.pad] -= dt/self.dx * (F[self.pad:] - F[:-self.pad])
            
            case 2:
                pass

            #AND SO ON, case 3....
    
    def apply_bc(self, dt, t):
        
        if self.periodic:
            match self.pad:
                case 1:
                    # In condition
                    self.left_boundary[0] = self.rho[-2].clone()

                    # Out condition
                    self.right_boundary[-1] = self.rho[1].clone()
                case 2:
                    # In condition
                    self.left_boundary[0] = self.rho[-4].clone()
                    self.left_boundary[1] = self.rho[-3].clone()

                    # Out condition
                    self.right_boundary[-1] = self.rho[3].clone()
                    self.right_boundary[-2] = self.rho[2].clone()
        
        else:

            if not self.right:
                # Neumann condition? The ghost cells are copies of the boundary cell.
                # No need for match here, just use the pad variable!!!
                match self.pad:
                    case 1:
                        self.right_boundary[-1] = self.rho[-2].clone()
                    case 2:
                        self.right_boundary[-1] = self.rho[-3].clone()
                        self.right_boundary[-2] = self.rho[-3].clone()

            if not self.left:
                # Set some influx to left boundary
                if self.boundary_fnc is None:
                    raise ValueError(f"No boundary function specified for road {self.id}!")
                
                else:
                    f_in = self.boundary_fnc(t)

                    ##############
                    # TODO: queue is calculated and updated here. Also perform CFL computations here.
                    #############

                    new_dt = dt.clone()
                    return  new_dt
        return dt


    def update_boundary_cells(self, dt):
        if not self.periodic:
            self.update_left_boundary(self.left_flux, dt)
            self.update_right_boundary(self.right_flux, dt)
        self.left_flux = torch.tensor(0.0)
        self.right_flux = torch.tensor(0.0)
    
    def update_right_flux(self, incoming_flux):

        self.right_flux = incoming_flux
        dt = torch.tensor(1)  # placeholder, e.g. one cell length
        return dt

    def update_left_flux(self, outgoing_flux):

        self.left_flux = outgoing_flux
        dt = torch.tensor(1)
        return dt

    def update_right_boundary(self, incoming_flux, dt):
        """
        Update the single right ghost cell (pad=1).
        """
        left = self.rho[-2]   # last interior cell
        right = self.rho[-1]  # right ghost cell

        F = torch.zeros(2)  # pad+1 = 2

        # Flux at interface (last interior ↔ ghost)
        F[0] = fv.Rusanov_Flux_2(left.clone(), right.clone(), self.gamma)
        # Boundary flux from junction
        F[1] = incoming_flux.clone()

        # FV update for the ghost cell
        self.right_boundary = self.right_boundary - dt / self.dx * (F[1:] - F[:-1])


    def update_left_boundary(self, outgoing_flux, dt):
            """
            Update the single left ghost cell (pad=1).
            """
            left = self.rho[0]   # left ghost cell
            right = self.rho[1]  # first interior cell

            F = torch.zeros(2)  # pad+1 = 2

            # Boundary flux from junction
            F[0] = outgoing_flux.clone()
            # Flux at interface (ghost ↔ first interior)
            F[1] = fv.Rusanov_Flux_2(left.clone(), right.clone(), self.gamma)

            # FV update for the ghost cell
            self.left_boundary = self.left_boundary - dt / self.dx * (F[1:] - F[:-1])


    def update_boundaries(self):
        # Set the boundary cells of rho equal to the artificial left boundary and right boundary
        self.rho[:self.pad] = self.left_boundary.clone()
        self.rho[-self.pad:] = self.right_boundary.clone()



    # HMMM. HERE WE SUDDENLY USE MAX_DENSE. WICH PROBABLY MAKES SENSE, BUT THEN THE WHOLE ANALYSIS FALLS APPART......
    # ALSO, THIS IS NOT THE TRUE DEMAND CURVE WE DEFINED EARLIER. WEHERE IS THE SPLIT OF THE FUNCTION COMPUTED.
    def demand(self):
        # clone needed?
        return self.max_dens * fv.D(self.rho[-1].clone(), self.gamma)

    def supply(self):
        # Clone needed?
        return self.max_dens * fv.S(self.rho[0].clone(), self.gamma)