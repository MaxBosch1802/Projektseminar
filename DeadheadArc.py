# Sven Mallach (2024)

class DeadheadArc:

    def __init__(self, ID, tN, typ, l):
        self.__id = ID
        self.__terminalNode = tN
        self.__type = typ
        self.__length = l
        self.__travelTimes = []

    def addTravelTime(self, time):
        self.__travelTimes.append(time)

    def getID(self):
        return self.__id

    def getTerminalNode(self):
        return self.__terminalNode

    def getType(self):
        return self.__type

    def getLength(self):
        return self.__length

    # tw is the index of the time window of interest (starting from zero)
    def getTravelTime(self, tw):
        return self.__travelTimes[tw]
