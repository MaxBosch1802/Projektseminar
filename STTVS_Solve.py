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
            for t in trips:
                print("Trip " + str(t.getID()) + " has start time " + str(t.getStartTime()) + " and end time " + str(t.getEndTime()) + ".")

        self.__x = pulp.LpVariable.dicts("x",(t.getID() for t in trips),cat=pulp.LpBinary)   #x_t Variable = 1 if trip t is used = 0 if not
        self.__z = pulp.LpVariable.dicts("z",((t.getID(),v.getID()) for t in trips for v in self.__sttvs.getFleet()),cat=pulp.LpBinary) # z_tv Variable = 1 if Trip t is assinged to Vehicle V else = 0
        self.__y = pulp.LpVariable.dicts("y",(v.getID() for v in self.__sttvs.getFleet()),cat=pulp.LpBinary)
        print("TODO")


    def generateConstraints(self):
        tH = self.__sttvs.getTimeHorizon()
        dirs = self.__sttvs.getDirections()
        
        for d in dirs:
            trips = d.getTrips()

        #Constraint (1): Each Direction must have exaclty on initial trip
        for d in dirs:
            self.__model += pulp.lpSum(self.__x[t.getID()] for t in trips if t.getInitialFinal() == "initial") == 1
        #Constraint (2): Each Direction must have exaclty on final trip
        for d in dirs:
            self.__model += pulp.lpSum(self.__x[t.getID()] for t in trips if t.getInitialFinal() == "final") == 1 

        #Constraint (4): Link x_t and z_tv variables
        for d in dirs:
            for t in trips:
                self.__model += pulp.lpSum(self.__z[t.getID(), v.getID()] for v in self.__sttvs.getFleet()) == self.__x[t.getID()] 
        #Constraint (8): Link z_tv and y_v variables
        for v in self.__sttvs.getFleet():
            for d in dirs:
                self.__model += pulp.lpSum(self.__z[t.getID(),v.getID()] for t in trips) <= len(trips) * self.__y[v.getID()]
        
        print("TODO")

    def generateObjectiveFunction(self):
        # Simple objective function: Minimize the sum of x_t, y_v, and z_tv variables
        obj_func = pulp.lpSum(self.__x.values()) + pulp.lpSum(self.__y.values()) + pulp.lpSum(self.__z.values())
        self.__model += obj_func

        print("Simple objective function generated.")            

    def solve(self):
        self.__model.solve()
        # Print the status of the solution
        print("Status:", pulp.LpStatus[self.__model.status])

        # Print the optimal objective value (total cost)
        print("Total Cost:", pulp.value(self.__model.objective))

        # Print the chosen variables and their values
        print("Chosen variables:")
        for var in self.__model.variables():
            if var.varValue > 0:
                print(var.name, "=", var.varValue)

    def writeLPFile(self, filename):
        self.__model.writeLP(filename)

