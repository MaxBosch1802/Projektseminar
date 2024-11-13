import pulp

from SeniorTTVS import SeniorTTVS

class STTVS_Solve:

    def __init__(self, sttvs):
        self.__sttvs = sttvs
        self.__model = pulp.LpProblem("STTVS", pulp.LpMinimize)

    def generateVariables(self):

        tH = self.__sttvs.getTimeHorizon()
        for i in range(1, len(tH)):
            print("Time Window: [" + str(tH[i-1]) + "," + str(tH[i]) + "]")

        dirs = self.__sttvs.getDirections()

        for d in dirs:
            for i in range(1, len(tH)):
                hw = d.getMaxHeadway(i-1) # maximum headway of time window i-1
                print(d.getType() + "-direction of line " + str(d.getLine()) + " has maximum headway " + str(hw) + " for the " + str(i-1) + "-th time Window")

            trips = d.getTrips()
#            for t in trips:
#                print("Trip " + str(t.getID()) + " has start time " + str(t.getStartTime()) + " and end time " + str(t.getEndTime()) + ".")

        print("TODO")

    def generateConstraints(self):
        print("TODO")

    def solve(self):
        self.__model.solve()

    def writeLPFile(self, filename):
        self.__model.writeLP(filename)
