# Sven Mallach (2024)

class Direction:

    def __init__(self, line, sN, eN, typ):
        self.__line = line
        self.__startNode = sN
        self.__endNode = eN
        self.__type = typ
        self.__maxHeadways = [] # maximum headways per time window
        self.__trips = []

    def addHeadway(self, hw):
        self.__maxHeadways.append(hw)

    def addTrip(self, trip):
        self.__trips.append(trip)

    def getLine(self):
        return self.__line

    def getStartNode(self):
        return self.__startNode

    def getEndNode(self):
        return self.__endNode

    def getType(self):
        return self.__type

    #tw is a time window index starting from zero
    def getMaxHeadway(self, tw):
        return self.__maxHeadways[tw]

    def getTrips(self):
        return self.__trips
