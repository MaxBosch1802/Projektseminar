# Sven Mallach (2024)

class Trip:

    def __init__(self, ID, sT, eT, mSAT, ell, initfin):
        self.__id = ID
        self.__startTime = sT
        self.__endTime = eT
        self.__mainStopArrivalTime = mSAT
        self.__length = ell
        self.__initialFinal = initfin

    def getID(self):
        return self.__id

    def getStartTime(self):
        return self.__startTime

    def getEndTime(self):
        return self.__endTime

    def getLength(self):
        return self.__length

    def getMainStopArrivalTime(self):
        return self.__mainStopArrivalTime

    def getInitialFinal(self):
        return self.__initialFinal


