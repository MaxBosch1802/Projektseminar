import pulp
import math
from SeniorTTVS import SeniorTTVS

class STTVS_Solve:

    def __init__(self, sttvs):
        self.sttvs = sttvs
        self.model = pulp.LpProblem("STTVS", pulp.LpMinimize)

        self.xvars = {}  
        self.yvars = {}  
        self.zvars = {} 

    def generateVariables(self):

        dirs = self.sttvs.getDirections()

        # Define xvars
        for d in dirs:
            for t in d.getTrips():
                self.xvars[t.getID()] = pulp.LpVariable('x_' + str(t.getID()), lowBound=0, upBound=1, cat='Binary')

        # Define yvars 
        for v in range(len(self.sttvs.getFleet())):  
            self.yvars[v] = pulp.LpVariable('y_' + str(v), lowBound=0, upBound=1, cat='Binary')

        # Define zvars
        for d in dirs:
            for v in range(len(self.sttvs.getFleet())):  
                for t in d.getTrips():
                    self.zvars[(t.getID(), v)] = pulp.LpVariable('z_' + str(t.getID()) + '_' + str(v), lowBound=0, upBound=1, cat='Binary')


    def generateConstraints(self):
        timeH = self.sttvs.getTimeHorizon()
        dirs = self.sttvs.getDirections()

        # Constraint (1): Exactly one initial trip per direction
        for d in dirs:
            lhs = []
            for t in d.getTrips():
                if t.getInitialFinal() == "initial":  # Check for trips marked as "initial"
                    lhs.append((self.xvars[t.getID()], 1.0))
            con = pulp.LpConstraint(e=pulp.LpAffineExpression(lhs),sense=pulp.LpConstraintEQ,rhs=1)
            self.model += con

        #Constraint (2): Each Direction must have exaclty on final trip
        for d in dirs:
            lhs = []
            for t in d.getTrips():
                if t.getInitialFinal() == "final":  # Check for trips marked as "final"
                    lhs.append((self.xvars[t.getID()], 1.0))
            con = pulp.LpConstraint(e=pulp.LpAffineExpression(lhs),sense=pulp.LpConstraintEQ,rhs=1)
            self.model += con
        

        # Constraint (3): 
        for d in dirs:
            for s in d.getTrips():
                if s.getInitialFinal() != "final":
                    msat = s.getMainStopArrivalTime()
                    twidx = -1

                    for i in range(1, len(timeH)):
                        if msat > timeH[i - 1] and msat <= timeH[i]:
                            twidx = i

                    IJmax = d.getMaxHeadway(twidx - 1)
                    lhs = []
                    lhs.append((self.xvars[s.getID()], -1.0))

                    for t in d.getTrips():
                        if t.getInitialFinal() != "initial":
                            # if (0 < t.getMainStopArrivalTime() - s.getMainStopArrivalTime() < IJmax and  # Condition from constraint (3)
                            #     (t.getTimeWindowIndex() == s.getTimeWindowIndex() or 
                            #     t.getTimeWindowIndex() == s.getTimeWindowIndex() + 1)):  # Time window condition
                            #     lhs.append((self.xvars[t.getID()], 1.0)) 
                            if (0 < t.getMainStopArrivalTime() - s.getMainStopArrivalTime() < IJmax):
                                lhs.append((self.xvars[t.getID()], 1.0)) 
                    con = pulp.LpConstraint(e=pulp.LpAffineExpression(lhs),sense=pulp.LpConstraintGE,rhs=0)
                    self.model += con
        
        # Constraint (4): Linking x and z variables
        for d in dirs:
            for t in d.getTrips():
                lhs = []
                lhs.append((self.xvars[t.getID()], -1.0))  # Coefficient of -1 for x_t

                for v in range(len(self.sttvs.getFleet())):  # Iterate over all vehicles
                    lhs.append((self.zvars[(t.getID(), v)], 1.0))  # Coefficient of +1 for z_tv

                con = pulp.LpConstraint(e=pulp.LpAffineExpression(lhs),sense=pulp.LpConstraintEQ, rhs=0)    
                self.model += con

        # #Constraint (8): Link z_tv and y_v variables
        for d in dirs:
            trips = d.getTrips()
            for i in range(len(trips) - 1):  # Iterate up to the second-to-last trip
                for j in range(i + 1, len(trips)):  # Iterate from i+1 to the end
                    if trips[j] in self.getIncompatibleSuccessorTrips(d, trips[i]) or trips[i] in self.getIncompatibleSuccessorTrips(d, trips[j]):
                        for v in range(len(self.sttvs.getFleet())):
                            con = pulp.LpConstraint(e=pulp.LpAffineExpression([(self.zvars[(trips[i].getID(), v)], 1.0), (self.zvars[(trips[j].getID(), v)], 1.0)]),sense=pulp.LpConstraintLE,rhs=1)                                                                                           
                            self.model += con
        print("TODO")

    def time2window(self, time):
        timeH = self.sttvs.getTimeHorizon()
        twidx = -1
        for i in range(1, len(timeH)):
            if (time >= timeH[i-1] and time < timeH[i]):
                twidx = i - 1
        return twidx

    def getIncompatibleSuccessorTrips(self, d, s):
        dirs = self.sttvs.getDirections()
        timeH = self.sttvs.getTimeHorizon()

        msat = s.getMainStopArrivalTime()
        twidx = self.time2window(msat) 

        minmaxstop = self.sttvs.getNodeByID(d.getEndNode()).getMinMaxStoppingTimes(twidx)
        enat = s.getEndTime()
        twidx_en = self.time2window(enat)
        depotminstop = self.sttvs.getNodeByID(0).getMinMaxStoppingTimes(twidx_en)[0]

        pullintime = 0
        for arc in self.sttvs.getDeadheadArcs():
            if arc.getTerminalNode() == d.getEndNode() and arc.getType() == "in":
                pullintime = arc.getTravelTime(twidx_en)
                break

        incompatible_s = []
        for e in dirs:
            if d.getEndNode() == e.getStartNode():
                for t in e.getTrips():
                    if ((t.getStartTime() - s.getEndTime() >= 0 and t.getStartTime() - s.getEndTime() < minmaxstop[0]) or
                            (t.getStartTime() - s.getEndTime() > minmaxstop[1])):
                        incompatible_s.append(t)
            else:
                if e == d:
                    continue
                for t in e.getTrips():
                    pulloutime = 0
                    for arc in self.sttvs.getDeadheadArcs():
                        if arc.getTerminalNode() == e.getStartNode() and arc.getType() == "out":
                            pullouttime = arc.getTravelTime(twidx_en)
                            break

                    if (0 <= t.getStartTime() - s.getEndTime() < pullintime + depotminstop + pullouttime):
                        incompatible_s.append(t)

        return incompatible_s

    def generateObjectiveFunction(self):
        # Objective function: Minimize the sum of chosen trips and vehicles

        obj = []

        # Add x variables (trips) to the objective function
        for d in self.sttvs.getDirections():
            for t in d.getTrips():
                obj.append((self.xvars[t.getID()], 1.0))  # Assuming each trip has a cost/value of 1

        # Add y variables (vehicles) to the objective function
        for v in range(len(self.sttvs.getFleet())):
            obj.append((self.yvars[v], 1.0))  # Assuming each vehicle has a cost/value of 1

        self.model += pulp.LpAffineExpression(obj)  # Set the objective function  

    def solve(self):
        self.model.solve()
        # Print the status of the model
        print("Status:", pulp.LpStatus[self.model.status])
        print("Total Costs:", pulp.value(self.model.objective)) 

        # Print the chosen x, y, and z variables
        # print("Chosen trips (x):")
        # for d in self.sttvs.getDirections():
        #     for t in d.getTrips():
        #         if self.xvars[t.getID()].varValue == 1:
        #             print(f"  - Trip {t.getID()} (Direction: {d.getLine()}_{d.getType()})")

        # print("Chosen vehicles (y):")
        # for v in range(len(self.sttvs.getFleet())):
        #     if self.yvars[v].varValue == 1:
        #         print(f"  - Vehicle {v}")

        # print("Trip-vehicle assignments (z):")
        # for d in self.sttvs.getDirections():
        #     for t in d.getTrips():
        #         for v in range(len(self.sttvs.getFleet())):
        #             if self.zvars[(t.getID(), v)].varValue == 1:
        #                 print(f"  - Trip {t.getID()} assigned to Vehicle {v}")

    def writeLPFile(self, filename):
        self.__model.writeLP(filename)
        
