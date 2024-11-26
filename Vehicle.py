# Sven Mallach (2024)

class Vehicle:

    def __init__(self, ID, typ, uc, pioc):
        self.__id = ID
        self.__type = typ
        self.__usageCost = uc
        self.__pullInOutCost = pioc

    def getID(self):
        return self.__id

    def getType(self):
        return self.__type

    def getUsageCost(self):
        return self.__usageCost

    def getPullInOutCost(self):
        return self.__pullInOutCost

class CombustionVehicle(Vehicle):

    def __init__(self, ID, typ, uc, pioc, ec):
        Vehicle.__init__(self, ID, typ, uc, pioc)
        self.__emissionCoefficient = ec

    def getEmissionCoefficient(self):
        return self.__emissionCoefficient

class ElectricVehicle(Vehicle):

    def __init__(self, ID, typ, uc, pioc, auto, maxCT, minCT):
        Vehicle.__init__(self, ID, typ, uc, pioc)
        self.__autonomy = auto
        self.__maxChargingTime = maxCT
        self.__minChargingTime = minCT

    def getMaxChargingTime(self):
        return self.__maxChargingTime

    def getMinChargingTime(self):
        return self.__minChargingTime

    def getAutonomy(self):
        return self.__autonomy

