import torch
import traffic_light as tl




class Junctino:


    roads = None
    entering = None
    leaving = None
    distribution = None
    #trafficlights = None


    def __init__(self, roads, entering, leaving, distribution, ):

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


        self.distribution = torch.tensor(distribution)

        # Priority needs to be implemented at some point
        self.priority = [1/len(entering)] * len(entering)
        
        #self.trafficlights = trafficlights
        #self.coupled_trafficlights = coupled_trafficlights
        self.road_in = [self.roads[i] for i in self.entering]
        self.road_out = [self.roads[i] for i in self.leaving]




