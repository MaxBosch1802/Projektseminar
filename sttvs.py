#!/usr/bin/python3

# Sven Mallach (2024)

import sys
import os
import json

from Node import Node
from DeadheadArc import DeadheadArc
from Direction import Direction
from Trip import Trip
from Vehicle import CombustionVehicle, ElectricVehicle
from SeniorTTVS import SeniorTTVS
from STTVS_Solve import STTVS_Solve

if __name__ == '__main__':

    if len(sys.argv) < 2:
        print("Missing argument(s). Usage: ./sttvs.py <json instance file>")
        exit(0)

    filename = sys.argv[1]
    with open(filename, "r") as f_in:
        data = json.load(f_in)

        nodes = []
        arcs = []
        dirs = []

        for nodeset in data["nodes"]:
            nodedict = nodeset["node"]
            nodeIDstr = nodedict["nodeName"]
            nodeID = None
            if ("dep" in nodeIDstr):
                nodeID = 0
            else:
                nodeID = int(nodeIDstr[4:])

            nodeBC = int(nodedict["breakCapacity"])
            nodeSCC = int(nodedict["slowChargeCapacity"])
            nodeFCC = int(nodedict["fastChargeCapacity"])

            thenode = Node(nodeID, nodeBC, nodeSCC, nodeFCC)

            nodeBTs = nodedict["breakingTimes"]
            for sT in nodeBTs:
                minST = int (sT["stoppingTime"]["minStoppingTime"])
                maxST = int (sT["stoppingTime"]["maxStoppingTime"])
                thenode.addStoppingTime(minST, maxST)

            nodes.append(thenode)

        for arcset in data["deadheadArcs"]:
            arcdict = arcset["deadheadArc"]
            arcID = int (arcdict["deadheadArcCode"])
            arcTN = int (arcdict["terminalNode"][4:])
            arcType = "in" if ("In" in arcdict["deadheadType"]) else "out"
            arcLen = float (arcdict["arcLength"])

            thearc = DeadheadArc(arcID, arcTN, arcType, arcLen)

            arcTTs = arcdict["travelTimes"]
            for tT in arcTTs:
                thearc.addTravelTime(int(tT))

            arcs.append(thearc)

        for dirset in data["directions"]:
            dirdict = dirset["direction"]
            dirLine = int (dirdict["lineName"][4:])
            dirSN = int (dirdict["startNode"][4:])
            dirEN = int (dirdict["endNode"][4:])
            dirType = "in" if ("in" in dirdict["directionType"]) else "out"

            thedir = Direction(dirLine, dirSN, dirEN, dirType)

            dirHws = dirdict["headways"]
            for hW in dirHws:
                thedir.addHeadway(int(hW["headway"]["maxHeadway"]))


            dirtrips = dirdict["trips"]
            for dT in dirtrips:
                tripset = dT["trip"]
                tripID = int(tripset["tripId"])
                tripST = int(tripset["startTime"])
                tripET = int(tripset["endTime"])
                tripMSAT = int(tripset["mainStopArrivalTime"])
                tripLen = float(tripset["lengthTrip"])
                tripIF = tripset["isInitialFinalTT"]

                thetrip = Trip(tripID, tripST, tripET, tripMSAT, tripLen, tripIF)

                thedir.addTrip(thetrip)

            dirs.append(thedir)

        times = data["timeHorizon"]

        fleet = []
        vehicleid = 0
        for vehicleset1 in data["fleet"]["vehicleList"]:
            vehicleset2 = vehicleset1["vehicleType"]
            if (vehicleset2["vehicleTypeName"] == "ICE"):
                ucost = int(vehicleset2["usageCost"])
                pioc = float(vehicleset2["pullInOutCost"])
                ec = float(vehicleset2["iceInfo"]["emissionCoefficient"])

                cv = CombustionVehicle(vehicleid, "ICE", ucost, pioc, ec)
                fleet.append(cv)
                vehicleid = vehicleid + 1
            else:
                ucost = int(vehicleset2["usageCost"])
                pioc = float(vehicleset2["pullInOutCost"])
                enum = int(vehicleset2["electricInfo"]["numberVehicle"])
                auto = int(vehicleset2["electricInfo"]["vehicleAutonomy"])
                minCT = int(vehicleset2["electricInfo"]["maxChargingTime"])
                maxCT = int(vehicleset2["electricInfo"]["minChargingTime"])

                ev = ElectricVehicle(vehicleid, "electric", ucost, pioc, enum, auto, minCT, maxCT)
                fleet.append(ev)
                vehicleid = vehicleid + 1
        phi = float(data["fleet"]["phi"])
        bCC = float(data["globalCost"]["breakCostCoefficient"])

        problem = SeniorTTVS(phi, bCC, nodes, arcs, dirs, times, fleet)

        solver = STTVS_Solve(problem)

        solver.generateVariables()
        solver.generateConstraints()

