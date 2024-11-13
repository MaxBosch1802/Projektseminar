# Sven Mallach (2024)

class Node:

    def __init__(self, ID, bC, sCC, fCC):
        self.__id = ID # the depot has ID 0
        self.__breakCapacity = bC
        self.__slowChargeCapacity = sCC
        self.__fastChargeCapacity = fCC
        self.__stoppingTimes = []

    def addStoppingTime(self, min, max):
        element = min, max
        self.__stoppingTimes.append(element)

    def getID(self):
        return self.__id

    def getBreakCapacity(self):
        return self.__breakCapacity

    def getslowChargeCapacity(self):
        return self.__slowChargeCapacity

    def getfastChargeCapacity(self):
        return self.__fastChargeCapacity

    #returns a tuple (delta^tw_min, delta^tw_max) where tw is a time window index starting from zero
    def getMinMaxStoppingTimes(self, tw):
        return self.__stoppingTimes[tw]
