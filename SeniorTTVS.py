# Sven Mallach (2024)

class SeniorTTVS:

    def __init__(self, phi, bCC, nodes, arcs, dirs, times, vehicles):
        self.__phi = phi
        self.__breakCostCoefficient = bCC
        self.__nodes = nodes
        self.__deadheadArcs = arcs
        self.__directions = dirs
        self.__timeHorizon = times
        self.__fleet = vehicles

    def getPhi(self):
        return self.__phi

    def getBreakCostCoefficient(self):
        return self.__breakCostCoefficient

    def getDeadheadArcs(self):
        return self.__deadheadArcs

    def getDirections(self):
        return self.__directions

    def getTimeHorizon(self):
        return self.__timeHorizon

    def getFleet(self):
        return self.__fleet

    def getNodes(self):
        return self.__nodes

    def getNodeByID(self, nr):
        for node in self.__nodes:
            if node.getID() == nr:
                return node
        return None
