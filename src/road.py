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
    max_dense = 1
    periodic = False


    def __init__(self, b, L, N, scheme, id, limiter, Vmax, initial):
        self.b = b 
        self.L = L 
        self.N = N
        self.scheme = scheme
        self.id = id
        self.limiter = limiter
        self.Vmax = Vmax
        self.initial = initial



        """
        There might be something wrong with the scaling here.
        THERE IS SOMETHING FISHY GOING ON HERE WITH THE LENGTHS OF STUFF AND SHIT
        """
        self.N_full
        self.dx = torch.tensor(1 / self.N_full) #Dimensionless.
        self.N_full = b * (N + 2*self.pad)
        j = torch.arange(0, self.N_full, 1)
        x = (j + 1/2) * self.dx
        self.rho = self.initial(x)


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

        
    def calculate_gamma(self):
        # L*b is the total length of the road
        # Vmax is initialized in km/h, but needs to be converted in m/s, hence divided by 3.6
        self.gamma = self.Vmax/(self.L * self.b * 3.6) 
    
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
                F = fv.Rusanov_Flux(self.rho, self.dx, dt, self.gamma)
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
        pass


    def update_boundaries(self):
        # Set the boundary cells of rho equal to the artificial left boundary and right boundary
        self.rho[:self.pad] = self.left_boundary.clone()
        self.rho[-self.pad:] = self.right_boundary.clone()
