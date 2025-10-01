import torch

class RoadNetwork:
    def __init__(self, roads, junctions, T, skip_iterations = 1, store_densities=True):
        """
        A network of roads with junctions.

        Args:
            roads (list[Road]): List of Road objects.
            junctions (list[Junction]): List of Junction objects.
            T (float): Final simulation time.
            store_densities (bool): Whether to store rho snapshots.
        """
        self.roads = roads
        self.junctions = junctions
        self.T = T
        self.store_densities = store_densities
        self.skip_iterations = skip_iterations

    def solve(self):
        t = 0.0
        iteration = 0

        # Store true initial condition at t=0.0
        history = {
            i: {0.0: road.rho[road.pad:-road.pad].clone()}
            for i, road in enumerate(self.roads)
        }

        while t < self.T:
            # Step 1: timestep
            dt = min([road.max_dt() for road in self.roads])

            # Step 2: junction fluxes
            for junction in self.junctions:
                junction.apply_bc(dt, t)

            # Step 3: boundaries for open roads
            for road in self.roads:
                road.apply_bc(dt, t)

            # Step 4: update boundary cells
            for road in self.roads:
                road.update_boundary_cells(dt)
                road.update_boundaries()

            # Step 5: solve inside roads
            for road in self.roads:
                road.solve_internally(dt)

            # Advance time
            t += float(dt)

            # Step 6: store updated state
            if self.store_densities and iteration % self.skip_iterations == 0:
                for i, road in enumerate(self.roads):
                    history[i][t] = road.rho[road.pad:-road.pad].clone()
            
            iteration += 1

        return history
