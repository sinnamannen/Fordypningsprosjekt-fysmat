import torch
import traffic_light as tl

class Junction:

    roads = None
    entering = None
    leaving = None
    distribution = None
    supply = None
    m = None
    n = None
    #trafficlights = None


    def __init__(self, roads, entering, leaving, distribution):

        ########
        # TODO: make code to check everything is initialized correctly
        #######

        self.roads = roads
        self.entering = entering

        for i in self.entering:
            self.roads[i].right = True

        self.leaving = leaving

        for j in self.leaving:
            self.roads[j].left = True
        
        self.n = len(self.entering)  # number of incoming roads
        self.m = len(self.leaving)   # number of outgoing roads


        self.distribution = torch.tensor(distribution)

        # Priority needs to be implemented at some point
        self.priority = [1/len(entering)] * len(entering)
        
        #self.trafficlights = trafficlights
        #self.coupled_trafficlights = coupled_trafficlights
        self.road_in = [self.roads[i] for i in self.entering]
        self.road_out = [self.roads[i] for i in self.leaving]

    

    def calculate_demand(self, active):
        demands = torch.zeros((self.n, self.m))
        for i, road in enumerate(self.road_in):
            demand = road.demand()
            demands[i,:] = active[i,:] * self.distribution[i] * demand
        return demands




    def calculate_activation(self, t):
        """
        Activation funciton: Creates the matrix showing how much the flux is scaled by from each incoming road to outgoing road.
        Scaled by some constant in [0,1].
        """

        active = torch.ones((self.n, self.m))

        return active
    


    def divide_flux(self, t):

        activation = self.calculate_activation(t)

        demand = self.calculate_demand(activation)

        capacities = [road.supply() for road in self.road_out]

        fluxes = torch.zeros((self.n, self.m))

                # This can likely be done without for loop...
        for j in range(self.m):
            # Update the flux from all roads into road j
            sum_influx = torch.sum(demand[:,j])
            if sum_influx > capacities[j]:
                # If the sum of the fluxes is larger than the capacity, scale down the fluxes
                fluxes[:,j] = demand[:,j] * capacities[j] / sum_influx
            else:
                fluxes[:,j] = demand[:,j]

        # For loops below can likely be removed
        # LIST OF TENSORS
        fluxes_in = [0]*self.n
        fluxes_out = [0]*self.m

        for i in range(self.n):
            # Scaling flux back to correspond to maximum density equal to 1
            fluxes_in[i] = torch.sum(fluxes[i])


        for j in range(self.m):
            # Scaling flux back to correspond to maximum density equal to 1
            fluxes_out[j] = torch.sum(fluxes[:,j])
        
        return fluxes_in, fluxes_out

    def apply_bc(self, dt, t):
        
        fluxes_in, fluxes_out = self.divide_flux(t)

        min_dt = dt + torch.tensor(1)
        for i, road in enumerate(self.road_in):
            min_dt_ = road.update_right_flux(fluxes_in[i] / road.max_dens)
            min_dt = torch.min(min_dt_, min_dt)
        
        for j, road in enumerate(self.road_out):
            min_dt_ = road.update_left_flux(fluxes_out[j] / road.max_dens)
            min_dt = torch.min(min_dt_, min_dt)

        min_dt = dt
        return min_dt






