from basefunctions import *
from vmtkfunctions import *
from vtk import * 
from math import *
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
import csv

pointsWithIds_sph = []
with open('particles_init_SPH.csv', 'rb') as csvfile:
	spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
	rowNumber = 0
	for row in spamreader:
		if(rowNumber > 1):
			tuppleSavingTheRelation = (int(row[0]),[1000 * float(row[1]), 1000 * float(row[2]), 1000 * float(row[3])])
			pointsWithIds_sph.append(tuppleSavingTheRelation)
		rowNumber += 1
		print row[0],row[1]

deformedPointsWithIds_sph = []
with open('particles_end_SPH.csv', 'rb') as csvfile:
	spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
	rowNumber = 0
	for row in spamreader:
		if(rowNumber > 1):
			tuppleSavingTheRelation = (int(row[0]),[1000 * float(row[1]), 1000 * float(row[2]), 1000 * float(row[3])])
			deformedPointsWithIds_sph.append(tuppleSavingTheRelation)
		rowNumber += 1
		print row[0],row[1]

pointsWithIds_FEM = []
with open('particles_init_FEM.csv', 'rb') as csvfile:
	spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
	rowNumber = 0
	for row in spamreader:
		if(rowNumber > 1):
			tuppleSavingTheRelation = (int(row[0]),[float(row[1]), float(row[2]), float(row[3])])
			pointsWithIds_FEM.append(tuppleSavingTheRelation)
		rowNumber += 1
		print row[0],row[1]

deformedPointsWithIds_FEM = []
with open('particles_end_FEM.csv', 'rb') as csvfile:
	spamreader = csv.reader(csvfile, delimiter=',', quotechar='|')
	rowNumber = 0
	for row in spamreader:
		if(rowNumber > 1):
			tuppleSavingTheRelation = (int(row[0]),[float(row[1]), float(row[2]), float(row[3])])
			deformedPointsWithIds_FEM.append(tuppleSavingTheRelation)
		rowNumber += 1
		print row[0],row[1]

def distanceBetweenPoints(p1, p2):
	return sqrt( (p1[0] - p2[0]) * (p1[0] - p2[0])
+ (p1[1] - p2[1]) * (p1[1] - p2[1])
+ (p1[2] - p2[2]) * (p1[2] - p2[2]))



def findIdWithMinDistance(p1, listOfPointsWithIds, numPoints):
	#for listElement
	minDist = 100000000
	minPoint = [0,0,0]
	minId = 0
	listOfTupleWithInfo = []
	for i in range(len(listOfPointsWithIds)):	
		p2 = listOfPointsWithIds[i][1]
		minId = listOfPointsWithIds[i][0]
		minPoint = listOfPointsWithIds[i][1]
		dist = distanceBetweenPoints(p1, p2)
		tupleWithInfo = (dist, minPoint, minId)
		listOfTupleWithInfo.append(tupleWithInfo)
#	print minId, minDist, p1, minPoint
#	print "\n"
#	print listOfTupleWithInfo
	sorted_by_distance = sorted(listOfTupleWithInfo, key=lambda tup: tup[0])
	listToReturn = []
	for i in range(numPoints):
		listToReturn.append(sorted_by_distance[i])
	return listToReturn

def addingSubPlotSPHJustOnePoint(listOfPointsWithIds, listOfDeformedPoints):
	listOfFoundIdsInOrder = []
	listOfFoundIdsInOrderSecondPlot = []
	deformedPointsZHighResolution = []
	listOfPointsSecondPlot = []
	listOfPointsX = []
	listOfPointsZ = []
	xSecondPlot = []
	resolutionOfSecondPlot = 15

	for i in range(resolutionOfSecondPlot):
		arg = -pi + (-acos(5. / 20.) + pi) * (i + 1) / resolutionOfSecondPlot
		listOfPointsSecondPlot.append([8.5 * sin(arg), 0, 18.5 * cos(arg)])	
		listOfPointsX.append(8.5*sin(arg))
		listOfPointsZ.append(18.5*cos(arg))

	for i in range(resolutionOfSecondPlot):		
		closestpoints_id_sph = findIdWithMinDistance(listOfPointsSecondPlot[i], listOfPointsWithIds, 3)
		listOfFoundIdsInOrderSecondPlot.append(closestpoints_id_sph)

	for idsIndex in range(resolutionOfSecondPlot):
#		print listOfDeformedPoints[listOfFoundIdsInOrderSecondPlot[i]], listOfDeformedPoints[listOfFoundIdsInOrderSecondPlot[i]][1][2]
		print "\n"
		closestId1 = listOfFoundIdsInOrderSecondPlot[idsIndex][0][2]
		closestId2 = listOfFoundIdsInOrderSecondPlot[idsIndex][1][2]
		closestId3 = listOfFoundIdsInOrderSecondPlot[idsIndex][2][2]

		dist1 = listOfFoundIdsInOrderSecondPlot[idsIndex][0][0]
		dist2 = listOfFoundIdsInOrderSecondPlot[idsIndex][1][0]
		dist3 = listOfFoundIdsInOrderSecondPlot[idsIndex][2][0]
		sumOfDists = dist1+dist2+dist3
		dist1 = dist1 / sumOfDists#change this
		dist2 = dist2 / sumOfDists
		dist3 = dist3 / sumOfDists

		zClosestId1 = 0
		zClosestId2 = 0
		for i in range(len(listOfPointsWithIds)):
#			print listOfPointsWithIds[i][0], listOfFoundIdsInOrderSecondPlot[idsIndex]
			if(listOfPointsWithIds[i][0] == closestId1):
				print listOfPointsWithIds[i][0], listOfFoundIdsInOrderSecondPlot[idsIndex]
				zClosestId1 = listOfDeformedPoints[i][1][2]
			if(listOfPointsWithIds[i][0] == closestId2):
				print listOfPointsWithIds[i][0], listOfFoundIdsInOrderSecondPlot[idsIndex]
				zClosestId2 = listOfDeformedPoints[i][1][2]
			if(listOfPointsWithIds[i][0] == closestId3):
				zClosestId3 = listOfDeformedPoints[i][1][2]
		deformedPointsZHighResolution.append((zClosestId1 * dist1))

	for idsIndex in range(resolutionOfSecondPlot):
#		print listOfDeformedPoints[listOfFoundIdsInOrderSecondPlot[i]], listOfDeformedPoints[listOfFoundIdsInOrderSecondPlot[i]][1][2]
		closestId1 = listOfFoundIdsInOrderSecondPlot[idsIndex][0][2]
		closestId2 = listOfFoundIdsInOrderSecondPlot[idsIndex][1][2]
		closestId3 = listOfFoundIdsInOrderSecondPlot[idsIndex][2][2]

		dist1 = listOfFoundIdsInOrderSecondPlot[idsIndex][0][0]
		dist2 = listOfFoundIdsInOrderSecondPlot[idsIndex][1][0]
		dist3 = listOfFoundIdsInOrderSecondPlot[idsIndex][2][0]
		sumOfDists = dist1+dist2+dist3
		dist1 = dist1 / sumOfDists
		dist2 = dist2 / sumOfDists
		dist3 = dist3 / sumOfDists
		for i in range(len(listOfPointsWithIds)):
			if(listOfPointsWithIds[i][0] == closestId1):
				xClosestId1 = listOfDeformedPoints[i][1][0]
			if(listOfPointsWithIds[i][0] == closestId2):
				xClosestId2 = listOfDeformedPoints[i][1][0]
			if(listOfPointsWithIds[i][0] == closestId3):
				xClosestId3 = listOfDeformedPoints[i][1][0]
		xSecondPlot.append((xClosestId1 * dist1 ))

#	for i in range(resolutionOfSecondPlot):
#		xSecondPlot.append(listOfDeformedPoints[listOfFoundIdsInOrderSecondPlot[i]][1][0])
	
	plt.plot(listOfPointsX, listOfPointsZ, 'r')
	plt.plot(listOfPointsX, listOfPointsZ, 'o', ms = 10)

	plt.show()
	plt.plot(xSecondPlot, deformedPointsZHighResolution, 'r')
	plt.plot(xSecondPlot, deformedPointsZHighResolution, 'o', ms = 10)

def addingSubPlotSPH(listOfPointsWithIds, listOfDeformedPoints):
	listOfFoundIdsInOrder = []
	listOfFoundIdsInOrderSecondPlot = []
	deformedPointsZHighResolution = []
	undeformedPointsZHighResolution = []
	listOfPointsSecondPlot = []
	listOfPointsX = []
	listOfPointsZ = []
	xSecondPlot = []
	xSecondPlotUndeformed = []
	resolutionOfSecondPlot = 15

	for i in range(resolutionOfSecondPlot):
		arg = -pi + (-acos(5. / 20.) + pi) * (i + 1) / resolutionOfSecondPlot
		listOfPointsSecondPlot.append([8.5 * sin(arg), 0, 18.5 * cos(arg)])	
		listOfPointsX.append(8.5*sin(arg))
		listOfPointsZ.append(18.5*cos(arg))

	for i in range(resolutionOfSecondPlot):		
		closestpoints_id_sph = findIdWithMinDistance(listOfPointsSecondPlot[i], listOfPointsWithIds, 3)
		listOfFoundIdsInOrderSecondPlot.append(closestpoints_id_sph)

	for idsIndex in range(resolutionOfSecondPlot):
#		print listOfDeformedPoints[listOfFoundIdsInOrderSecondPlot[i]], listOfDeformedPoints[listOfFoundIdsInOrderSecondPlot[i]][1][2]
		print "\n"
		closestId1 = listOfFoundIdsInOrderSecondPlot[idsIndex][0][2]
		closestId2 = listOfFoundIdsInOrderSecondPlot[idsIndex][1][2]
		closestId3 = listOfFoundIdsInOrderSecondPlot[idsIndex][2][2]

		dist1 = 1. / listOfFoundIdsInOrderSecondPlot[idsIndex][0][0]
		dist2 = 1. / listOfFoundIdsInOrderSecondPlot[idsIndex][1][0]
		dist3 = 1. / listOfFoundIdsInOrderSecondPlot[idsIndex][2][0]
		sumOfDists = dist1+dist2+dist3
		dist1 = dist1 / sumOfDists#change this
		dist2 = dist2 / sumOfDists
		dist3 = dist3 / sumOfDists

		zClosestId1 = 0
		zClosestId2 = 0
		for i in range(len(listOfPointsWithIds)):
#			print listOfPointsWithIds[i][0], listOfFoundIdsInOrderSecondPlot[idsIndex]
			if(listOfPointsWithIds[i][0] == closestId1):
				print listOfPointsWithIds[i][0], listOfFoundIdsInOrderSecondPlot[idsIndex]
				zClosestId1 = listOfDeformedPoints[i][1][2]
				zClosestId1Undeformed = listOfPointsWithIds[i][1][2]
			if(listOfPointsWithIds[i][0] == closestId2):
				print listOfPointsWithIds[i][0], listOfFoundIdsInOrderSecondPlot[idsIndex]
				zClosestId2 = listOfDeformedPoints[i][1][2]
				zClosestId2Undeformed = listOfPointsWithIds[i][1][2]
			if(listOfPointsWithIds[i][0] == closestId3):
				zClosestId3 = listOfDeformedPoints[i][1][2]
				zClosestId3Undeformed = listOfPointsWithIds[i][1][2]
		deformedPointsZHighResolution.append((zClosestId1 * dist1 + zClosestId2 * dist2 + zClosestId3 * dist3))
		undeformedPointsZHighResolution.append((zClosestId1Undeformed * dist1 + zClosestId2Undeformed * dist2 + zClosestId3Undeformed * dist3))
	powerOfInverseDistances = 3
	for idsIndex in range(resolutionOfSecondPlot):
#		print listOfDeformedPoints[listOfFoundIdsInOrderSecondPlot[i]], listOfDeformedPoints[listOfFoundIdsInOrderSecondPlot[i]][1][2]
		closestId1 = listOfFoundIdsInOrderSecondPlot[idsIndex][0][2]
		closestId2 = listOfFoundIdsInOrderSecondPlot[idsIndex][1][2]
		closestId3 = listOfFoundIdsInOrderSecondPlot[idsIndex][2][2]

		dist1 = pow(1. / listOfFoundIdsInOrderSecondPlot[idsIndex][0][0], powerOfInverseDistances)
		dist2 = pow(1. / listOfFoundIdsInOrderSecondPlot[idsIndex][1][0], powerOfInverseDistances)
		dist3 = pow(1. / listOfFoundIdsInOrderSecondPlot[idsIndex][2][0], powerOfInverseDistances)
		sumOfDists = dist1+dist2+dist3
		dist1 = dist1 / sumOfDists
		dist2 = dist2 / sumOfDists
		dist3 = dist3 / sumOfDists
		for i in range(len(listOfPointsWithIds)):
			if(listOfPointsWithIds[i][0] == closestId1):
				xClosestId1 = listOfDeformedPoints[i][1][0]
				xClosestId1Undeformed = listOfPointsWithIds[i][1][0]				
			if(listOfPointsWithIds[i][0] == closestId2):
				xClosestId2 = listOfDeformedPoints[i][1][0]
				xClosestId2Undeformed = listOfPointsWithIds[i][1][0]				
			if(listOfPointsWithIds[i][0] == closestId3):
				xClosestId3 = listOfDeformedPoints[i][1][0]
				xClosestId3Undeformed = listOfPointsWithIds[i][1][0]
		xSecondPlot.append((xClosestId1 * dist1 + xClosestId2 * dist2 + xClosestId3 * dist3))
		xSecondPlotUndeformed.append((xClosestId1Undeformed * dist1 + xClosestId2Undeformed * dist2 + xClosestId3Undeformed * dist3))

#	for i in range(resolutionOfSecondPlot):
#		xSecondPlot.append(listOfDeformedPoints[listOfFoundIdsInOrderSecondPlot[i]][1][0])
	
	l2, = plt.plot(listOfPointsX, listOfPointsZ, linestyle = '-', color = 'gray', ms = 10, linewidth = 3)
#	plt.plot(listOfPointsX, listOfPointsZ, 'o', ms = 10)

#	plt.show()
	l3, = plt.plot(xSecondPlot, deformedPointsZHighResolution, '-bo', ms = 10, linewidth = 3)
#	plt.plot(xSecondPlot, deformedPointsZHighResolution, 'ro', ms = 10)
	return l2, l3
#	plt.plot(xSecondPlotUndeformed, undeformedPointsZHighResolution, 'blue')
#	plt.plot(xSecondPlotUndeformed, undeformedPointsZHighResolution, 'o', ms = 10)

def addingSubPlotFEM(listOfPointsWithIds, listOfDeformedPoints):
	listOfFoundIdsInOrder = []
	listOfFoundIdsInOrderSecondPlot = []
	deformedPointsZHighResolution = []
	undeformedPointsZHighResolution = []
	listOfPointsSecondPlot = []
	listOfPointsX = []
	listOfPointsZ = []
	xSecondPlot = []
	xSecondPlotUndeformed = []
	resolutionOfSecondPlot = 15

	for i in range(resolutionOfSecondPlot):
		arg = -pi + (-acos(5. / 20.) + pi) * (i + 1) / resolutionOfSecondPlot
		listOfPointsSecondPlot.append([8.5 * sin(arg), 0, 18.5 * cos(arg)])	
		listOfPointsX.append(8.5*sin(arg))
		listOfPointsZ.append(18.5*cos(arg))

	for i in range(resolutionOfSecondPlot):		
		closestpoints_id_sph = findIdWithMinDistance(listOfPointsSecondPlot[i], listOfPointsWithIds, 3)
		listOfFoundIdsInOrderSecondPlot.append(closestpoints_id_sph)

	for idsIndex in range(resolutionOfSecondPlot):
#		print listOfDeformedPoints[listOfFoundIdsInOrderSecondPlot[i]], listOfDeformedPoints[listOfFoundIdsInOrderSecondPlot[i]][1][2]
		print "\n"
		closestId1 = listOfFoundIdsInOrderSecondPlot[idsIndex][0][2]
		closestId2 = listOfFoundIdsInOrderSecondPlot[idsIndex][1][2]
		closestId3 = listOfFoundIdsInOrderSecondPlot[idsIndex][2][2]

		dist1 = 1. / listOfFoundIdsInOrderSecondPlot[idsIndex][0][0]
		dist2 = 1. / listOfFoundIdsInOrderSecondPlot[idsIndex][1][0]
		dist3 = 1. / listOfFoundIdsInOrderSecondPlot[idsIndex][2][0]
		sumOfDists = dist1+dist2+dist3
		dist1 = dist1 / sumOfDists#change this
		dist2 = dist2 / sumOfDists
		dist3 = dist3 / sumOfDists

		zClosestId1 = 0
		zClosestId2 = 0
		for i in range(len(listOfPointsWithIds)):
#			print listOfPointsWithIds[i][0], listOfFoundIdsInOrderSecondPlot[idsIndex]
			if(listOfPointsWithIds[i][0] == closestId1):
				print listOfPointsWithIds[i][0], listOfFoundIdsInOrderSecondPlot[idsIndex]
				zClosestId1 = listOfDeformedPoints[i][1][2]
				zClosestId1Undeformed = listOfPointsWithIds[i][1][2]
			if(listOfPointsWithIds[i][0] == closestId2):
				print listOfPointsWithIds[i][0], listOfFoundIdsInOrderSecondPlot[idsIndex]
				zClosestId2 = listOfDeformedPoints[i][1][2]
				zClosestId2Undeformed = listOfPointsWithIds[i][1][2]
			if(listOfPointsWithIds[i][0] == closestId3):
				zClosestId3 = listOfDeformedPoints[i][1][2]
				zClosestId3Undeformed = listOfPointsWithIds[i][1][2]
		deformedPointsZHighResolution.append((zClosestId1 * dist1 + zClosestId2 * dist2 + zClosestId3 * dist3))
		undeformedPointsZHighResolution.append((zClosestId1Undeformed * dist1 + zClosestId2Undeformed * dist2 + zClosestId3Undeformed * dist3))
	powerOfInverseDistances = 3
	for idsIndex in range(resolutionOfSecondPlot):
#		print listOfDeformedPoints[listOfFoundIdsInOrderSecondPlot[i]], listOfDeformedPoints[listOfFoundIdsInOrderSecondPlot[i]][1][2]
		closestId1 = listOfFoundIdsInOrderSecondPlot[idsIndex][0][2]
		closestId2 = listOfFoundIdsInOrderSecondPlot[idsIndex][1][2]
		closestId3 = listOfFoundIdsInOrderSecondPlot[idsIndex][2][2]

		dist1 = pow(1. / listOfFoundIdsInOrderSecondPlot[idsIndex][0][0], powerOfInverseDistances)
		dist2 = pow(1. / listOfFoundIdsInOrderSecondPlot[idsIndex][1][0], powerOfInverseDistances)
		dist3 = pow(1. / listOfFoundIdsInOrderSecondPlot[idsIndex][2][0], powerOfInverseDistances)
		sumOfDists = dist1+dist2+dist3
		dist1 = dist1 / sumOfDists
		dist2 = dist2 / sumOfDists
		dist3 = dist3 / sumOfDists
		for i in range(len(listOfPointsWithIds)):
			if(listOfPointsWithIds[i][0] == closestId1):
				xClosestId1 = listOfDeformedPoints[i][1][0]
				xClosestId1Undeformed = listOfPointsWithIds[i][1][0]				
			if(listOfPointsWithIds[i][0] == closestId2):
				xClosestId2 = listOfDeformedPoints[i][1][0]
				xClosestId2Undeformed = listOfPointsWithIds[i][1][0]				
			if(listOfPointsWithIds[i][0] == closestId3):
				xClosestId3 = listOfDeformedPoints[i][1][0]
				xClosestId3Undeformed = listOfPointsWithIds[i][1][0]
		xSecondPlot.append((xClosestId1 * dist1 + xClosestId2 * dist2 + xClosestId3 * dist3))
		xSecondPlotUndeformed.append((xClosestId1Undeformed * dist1 + xClosestId2Undeformed * dist2 + xClosestId3Undeformed * dist3))

#	for i in range(resolutionOfSecondPlot):
#		xSecondPlot.append(listOfDeformedPoints[listOfFoundIdsInOrderSecondPlot[i]][1][0])
	
	l2, = plt.plot(listOfPointsX, listOfPointsZ, linestyle = '-', color = 'gray', ms = 10, linewidth = 3)
#	plt.plot(listOfPointsX, listOfPointsZ, 'o', ms = 10)

#	plt.show()
	l3, = plt.plot(xSecondPlot, deformedPointsZHighResolution, '--y^', ms = 10, linewidth = 3)
#	plt.plot(xSecondPlot, deformedPointsZHighResolution, 'y^', ms = 10)
	return l2, l3


# ply2vtk('/home/neubias/Desktop/478.ply', 'home/neubias/Desktop/478.vtk')
#m = readvtk('/home/eric/dataForSimulations/sander314-mechbench-f541b4f0e0b3/problem 1/CardioMechanics/139587/undeformedPolyData.vtk')
#********** We just read the surface **********
def addingSubPlotsOfDeformation(undeformedFileNameFolder):
	undeformed = "/undeformed.vtk"
	undeformedFileName = undeformedFileNameFolder + undeformed
	deformed = "/deformed.vtk"
	deformedFileName = undeformedFileNameFolder + deformed

	reader = vtkUnstructuredGridReader();
	reader.SetFileName(undeformedFileName);
	reader.Update();
	undeformedMesh = reader.GetOutput()

	readerDeformed = vtkUnstructuredGridReader();
	readerDeformed.SetFileName(deformedFileName);
	readerDeformed.Update();
	deformedMesh = readerDeformed.GetOutput()

	locator = vtk.vtkPointLocator()
	locator.SetDataSet(undeformedMesh)
	locator.BuildLocator()

	listOfFoundIdsInOrder = []
	listOfFoundIdsInOrderSecondPlot = []
	deformedPointsZHighResolution = []
	listOfPointsSecondPlot = []
	xSecondPlot = []
	resolutionOfSecondPlot = 15
	for i in range(resolutionOfSecondPlot):
		arg = -pi + (-acos(5. / 20.) + pi) * (i + 1) / resolutionOfSecondPlot
		listOfPointsSecondPlot.append([8.5 * sin(arg), 0, 18.5 * cos(arg)])	
		#xSecondPlot.append(i * 10./resolutionOfSecondPlot);

	for i in range(resolutionOfSecondPlot):
		closestpoint_id = locator.FindClosestPoint(listOfPointsSecondPlot[i])
		listOfFoundIdsInOrderSecondPlot.append(closestpoint_id)

	for i in range(resolutionOfSecondPlot):
		deformedPointsZHighResolution.append(deformedMesh.GetPoint(listOfFoundIdsInOrderSecondPlot[i])[2])

	for i in range(resolutionOfSecondPlot):
		xSecondPlot.append(deformedMesh.GetPoint(listOfFoundIdsInOrderSecondPlot[i])[0])

	print len(xSecondPlot)
	print len(deformedPointsZHighResolution)

	
	l1, = plt.plot(xSecondPlot, deformedPointsZHighResolution)
	plt.plot(xSecondPlot, deformedPointsZHighResolution, 'o', ms = 10)

	return l1
#	plt.set_title(undeformedFileNameFolder)
#	plt.show()

#for strains SPH:
def addingSubPlotSPHStrain(listOfPointsWithIds, listOfDeformedPoints):
	listOfFoundIdsInOrder = []
	listOfFoundIdsInOrderSecondPlot = []
	deformedPointsZHighResolution = []
	listOfPointsSecondPlot = []
	xSecondPlot = []

	listOfPointsX = []
	listOfPointsZ = []
	xSecondPlot = []
	resolutionOfSecondPlot = 9

	for i in range(resolutionOfSecondPlot):
		arg = -pi + (-acos(5. / 20.) + pi) * (i + 1) / resolutionOfSecondPlot
		listOfPointsSecondPlot.append([8.5 * sin(arg), 0, 18.5 * cos(arg)])	
		listOfPointsX.append(8.5*sin(arg))
		listOfPointsZ.append(18.5*cos(arg))

	listOfPointsLongitudinalEndo = []
	listOfPointsLongitudinalMid = []
	listOfPointsLongitudinalEpi = []
	listOfPointsCircumferentialEndo = []
	listOfPointsCircumferentialMid = []
	listOfPointsCircumferentialEpi = []

	n_u = 10
	u1 = -pi
	tEndo = 0.1
	tMid = 0.5
	tEpi = 0.9
	rsEndo = 7
	rlEndo = 17
	rsEpi = 10
	rlEpi = 20
	rsMid = rsEndo + 0.5 * (rsEpi - rsEndo)
	rlMid = rlEndo + 0.5 * (rlEpi - rlEndo)
	
#	for i in range(resolutionOfSecondPlot):
#		arg = -pi + (-acos(5. / 20.) + pi) * (i + 1) / resolutionOfSecondPlot
#		listOfPointsSecondPlot.append([8.5 * sin(arg), 0, 18.5 * cos(arg)])	
#		listOfPointsX.append(8.5*sin(arg))
#		listOfPointsZ.append(18.5*cos(arg))

	print rsMid, rlMid

	for i in range(11):
		v = 0
		u2 = -acos(5. / (17 + 3 * tEndo))
		u = u1 + (u2 - u1) / n_u * (i + 1) * 0.95
		x = rsEndo * sin(u)
		y = 0
		z = rlEndo * cos(u)
		listOfPointsLongitudinalEndo.append([x, y, z])
		print x,y,z

		u2 = -acos(5. / (17 + 3 * tMid))
		u = u1 + (u2 - u1) / n_u * (i + 1) * 0.95
		x = rsMid * sin(u)
		y = 0
		z = rlMid * cos(u)
		listOfPointsLongitudinalMid.append([x, y, z])

		u2 = -acos(5. / (17 + 3 * tEpi))
		u = u1 + (u2 - u1) / n_u * (i + 1) * 0.95
		x = rsEpi * sin(u)
		y = 0
		z = rlEpi * cos(u)
		listOfPointsLongitudinalEpi.append([x, y, z])

		v = math.pi / 10.
		u2 = -acos(5. / (17 + 3 * tEndo))
		u = u1 + (u2 - u1) / n_u * (i + 1) * 0.95
		x = rsEndo * sin(u) * cos(v)
		y = rsEndo * sin(u) * sin(v)
		z = rlEndo * cos(u)
		listOfPointsCircumferentialEndo.append([x, y, z])

		u2 = -acos(5. / (17 + 3 * tMid))
		u = u1 + (u2 - u1) / n_u * (i + 1) * 0.95
		x = rsMid * sin(u) * cos(v)
		y = rsMid * sin(u) * sin(v)
		z = rlMid * cos(u)
		listOfPointsCircumferentialMid.append([x, y, z])

		u2 = -acos(5. / (17 + 3 * tEpi))
		u = u1 + (u2 - u1) / n_u * (i + 1) * 0.95
		x = rsEpi * sin(u) * cos(v)
		y = rsEpi * sin(u) * sin(v)
		z = rlEpi * cos(u)
		listOfPointsCircumferentialEpi.append([x, y, z])
	listOfFoundIdsLongitudinalEndo = []
	listOfFoundIdsLongitudinalMid = []
	listOfFoundIdsLongitudinalEpi = []
	listOfFoundIdsCircumferentialEndo = []
	listOfFoundIdsCircumferentialMid = []
	listOfFoundIdsCircumferentialEpi = []
	numOfPointsToBeConsideredCloseToTheCenterLine = 30
	exponentValueOfInverseDistance = 1
	for i in range(11):
		closestpoint_id_sph = findIdWithMinDistance(listOfPointsLongitudinalEndo[i], listOfPointsWithIds, numOfPointsToBeConsideredCloseToTheCenterLine)
		listOfFoundIdsLongitudinalEndo.append(closestpoint_id_sph)
		closestpoint_id_sph = findIdWithMinDistance(listOfPointsLongitudinalMid[i], listOfPointsWithIds, numOfPointsToBeConsideredCloseToTheCenterLine)
		listOfFoundIdsLongitudinalMid.append(closestpoint_id_sph)
		closestpoint_id_sph = findIdWithMinDistance(listOfPointsLongitudinalEpi[i], listOfPointsWithIds, numOfPointsToBeConsideredCloseToTheCenterLine)
		listOfFoundIdsLongitudinalEpi.append(closestpoint_id_sph)
		closestpoint_id_sph = findIdWithMinDistance(listOfPointsCircumferentialEndo[i], listOfPointsWithIds, numOfPointsToBeConsideredCloseToTheCenterLine)
		listOfFoundIdsCircumferentialEndo.append(closestpoint_id_sph)
		closestpoint_id_sph = findIdWithMinDistance(listOfPointsCircumferentialMid[i], listOfPointsWithIds, numOfPointsToBeConsideredCloseToTheCenterLine)
		listOfFoundIdsCircumferentialMid.append(closestpoint_id_sph)
		closestpoint_id_sph = findIdWithMinDistance(listOfPointsCircumferentialEpi[i], listOfPointsWithIds, numOfPointsToBeConsideredCloseToTheCenterLine)
		listOfFoundIdsCircumferentialEpi.append(closestpoint_id_sph)

	print "STRAINS"
	circumferentialStrainEndo = []
	longitudinalStrainEndo = []
	transmuralStrainEndo = []
	circumferentialStrainMid = []
	longitudinalStrainMid = []
	transmuralStrainMid = []
	circumferentialStrainEpi = []
	longitudinalStrainEpi = []
	transmuralStrainEpi = []
	listOfStrainsLongEndo = []
	listOfAnalyticalComparedStrainsLongEndo = []
	listOfStrainsLongMid = []
	listOfStrainsLongEpi = []
	listOfStrainsCircEndo = []
	listOfStrainsCircMid = []
	listOfStrainsCircEpi = []
	listOfStrainsTransEndo = []
	listOfStrainsTransMid = []
	listOfStrainsTransEpi = []

	listOfStrainsForEachInterpolationPoint = []
	listOfAllStrains = []


	listUndeformedLongEndoX = []
	listDeformedLongEndoX = []
	listUndeformedLongEndoZ = []
	listNumericalUndeformedLongEndoX = []
	listNumericalUndeformedLongEndoZ = []
	listDeformedLongEndoZ = []

	for idsIndex in range(10):
		undeformedPointLongEndo = [0, 0, 0]
		undeformedPointLongMid = [0, 0, 0]
		undeformedPointLongEpi = [0, 0, 0]
		undeformedPointCircumfEndo = [0, 0, 0]
		undeformedPointCircumfMid = [0, 0, 0]
		undeformedPointCircumfEpi = [0, 0, 0]

		deformedPointLongEndo = [0, 0, 0]
		deformedPointLongMid = [0, 0, 0]
		deformedPointLongEpi = [0, 0, 0]
		deformedPointCircumfEndo = [0, 0, 0]
		deformedPointCircumfMid = [0, 0, 0]
		deformedPointCircumfEpi = [0, 0, 0]

		numeratorUndeformedLongEndo  = [0, 0, 0]
		denominatorUndeformedLongEndo = 0
		numeratorDeformedLongEndo  = [0, 0, 0]
		denominatorDeformedLongEndo = 0
		numeratorUndeformedLongEndoNext  = [0, 0, 0]
		denominatorUndeformedLongEndoNext = 0
		numeratorDeformedLongEndoNext  = [0, 0, 0]
		denominatorDeformedLongEndoNext = 0

		numeratorUndeformedLongMid  = [0, 0, 0]
		denominatorUndeformedLongMid = 0
		numeratorDeformedLongMid  = [0, 0, 0]
		denominatorDeformedLongMid = 0
		numeratorUndeformedLongMidNext  = [0, 0, 0]
		denominatorUndeformedLongMidNext = 0
		numeratorDeformedLongMidNext  = [0, 0, 0]
		denominatorDeformedLongMidNext = 0

		numeratorUndeformedLongEpi  = [0, 0, 0]
		denominatorUndeformedLongEpi = 0
		numeratorDeformedLongEpi  = [0, 0, 0]
		denominatorDeformedLongEpi = 0
		numeratorUndeformedLongEpiNext  = [0, 0, 0]
		denominatorUndeformedLongEpiNext = 0
		numeratorDeformedLongEpiNext  = [0, 0, 0]
		denominatorDeformedLongEpiNext = 0

		numeratorUndeformedCircEndo  = [0, 0, 0]
		denominatorUndeformedCircEndo = 0
		numeratorDeformedCircEndo  = [0, 0, 0]
		denominatorDeformedCircEndo = 0

		numeratorUndeformedCircMid  = [0, 0, 0]
		denominatorUndeformedCircMid = 0
		numeratorDeformedCircMid  = [0, 0, 0]
		denominatorDeformedCircMid = 0

		numeratorUndeformedCircEpi  = [0, 0, 0]
		denominatorUndeformedCircEpi = 0
		numeratorDeformedCircEpi  = [0, 0, 0]
		denominatorDeformedCircEpi = 0
		for nPoints in range(numOfPointsToBeConsideredCloseToTheCenterLine):
			closestIdLongEndo = listOfFoundIdsLongitudinalEndo[idsIndex][nPoints][2]
			closestIdLongEndoNext = listOfFoundIdsLongitudinalEndo[idsIndex + 1][nPoints][2]
			closestIdLongMid = listOfFoundIdsLongitudinalMid[idsIndex][nPoints][2]
			closestIdLongMidNext = listOfFoundIdsLongitudinalMid[idsIndex + 1][nPoints][2]
			closestIdLongEpi = listOfFoundIdsLongitudinalEpi[idsIndex][nPoints][2]
			closestIdLongEpiNext = listOfFoundIdsLongitudinalEpi[idsIndex + 1][nPoints][2]

			closestIdCircumferentialEndo = listOfFoundIdsCircumferentialEndo[idsIndex][nPoints][2]
			closestIdCircumferentialMid = listOfFoundIdsCircumferentialMid[idsIndex][nPoints][2]
			closestIdCircumferentialEpi = listOfFoundIdsCircumferentialEpi[idsIndex][nPoints][2]
	
			invDistLongEndo = pow(1. / listOfFoundIdsLongitudinalEndo[idsIndex][nPoints][0], exponentValueOfInverseDistance)
			invDistLongEndoNext = pow(1. / listOfFoundIdsLongitudinalEndo[idsIndex + 1][nPoints][0], exponentValueOfInverseDistance)
			invDistLongMid = pow(1. / listOfFoundIdsLongitudinalMid[idsIndex][nPoints][0], exponentValueOfInverseDistance)
			invDistLongMidNext = pow(1. / listOfFoundIdsLongitudinalMid[idsIndex + 1][nPoints][0], exponentValueOfInverseDistance)
			invDistLongEpi = pow(1. / listOfFoundIdsLongitudinalEpi[idsIndex][nPoints][0], exponentValueOfInverseDistance)
			invDistLongEpiNext = pow(1. / listOfFoundIdsLongitudinalEpi[idsIndex + 1][nPoints][0], exponentValueOfInverseDistance)

			invDistCircumferentialEndo = pow(1. / listOfFoundIdsCircumferentialEndo[idsIndex][nPoints][0], exponentValueOfInverseDistance)
			invDistCircumferentialMid = pow(1. / listOfFoundIdsCircumferentialMid[idsIndex][nPoints][0], exponentValueOfInverseDistance)
			invDistCircumferentialEpi = pow(1. / listOfFoundIdsCircumferentialEpi[idsIndex][nPoints][0], exponentValueOfInverseDistance)

#			print invDistLongEndo, invDistCircumferentialEndo
#			invDistLongEndo = 1. 
#			invDistLongEndoNext = 1. 
#			invDistLongMid = 1. 
#			invDistLongMidNext = 1. 
#			invDistLongEpi = 1. 
#			invDistLongEpiNext = 1. 
#			invDistCircumferentialEndo = 1. 
#			invDistCircumferentialMid = 1. 
#			invDistCircumferentialEpi = 1.

			for i in range(len(listOfPointsWithIds)):
				if((listOfPointsWithIds[i][0] == closestIdLongEndo)):
					undeformedLongEndo = listOfPointsWithIds[i][1]
					deformedLongEndo = listOfDeformedPoints[i][1]
				if((listOfPointsWithIds[i][0] == closestIdLongEndoNext)):
					undeformedLongEndoNext = listOfPointsWithIds[i][1]
					deformedLongEndoNext = listOfDeformedPoints[i][1]
				if((listOfPointsWithIds[i][0] == closestIdLongMid)):
					undeformedLongMid = listOfPointsWithIds[i][1]
					deformedLongMid = listOfDeformedPoints[i][1]
				if((listOfPointsWithIds[i][0] == closestIdLongMidNext)):
					undeformedLongMidNext = listOfPointsWithIds[i][1]
					deformedLongMidNext = listOfDeformedPoints[i][1]
				if((listOfPointsWithIds[i][0] == closestIdLongEpi)):
					undeformedLongEpi = listOfPointsWithIds[i][1]
					deformedLongEpi = listOfDeformedPoints[i][1]
				if((listOfPointsWithIds[i][0] == closestIdLongEpiNext)):
					undeformedLongEpiNext = listOfPointsWithIds[i][1]
					deformedLongEpiNext = listOfDeformedPoints[i][1]
				if((listOfPointsWithIds[i][0] == closestIdCircumferentialEndo)):
					undeformedCircumferentialEndo = listOfPointsWithIds[i][1]
					deformedCircumferentialEndo = listOfDeformedPoints[i][1]
				if((listOfPointsWithIds[i][0] == closestIdCircumferentialMid)):
					undeformedCircumferentialMid = listOfPointsWithIds[i][1]
					deformedCircumferentialMid = listOfDeformedPoints[i][1]
				if((listOfPointsWithIds[i][0] == closestIdCircumferentialEpi)):
					undeformedCircumferentialEpi = listOfPointsWithIds[i][1]
					deformedCircumferentialEpi = listOfDeformedPoints[i][1]

			numeratorUndeformedLongEndo = [numeratorUndeformedLongEndo[index] + undeformedLongEndo[index] * invDistLongEndo for index in range(3)]
			denominatorUndeformedLongEndo = denominatorUndeformedLongEndo + invDistLongEndo
			numeratorDeformedLongEndo = [numeratorDeformedLongEndo[index] + deformedLongEndo[index] * invDistLongEndo for index in range(3)]
			denominatorDeformedLongEndo = denominatorDeformedLongEndo + invDistLongEndo
			numeratorUndeformedLongEndoNext = [numeratorUndeformedLongEndoNext[index] + undeformedLongEndoNext[index] * invDistLongEndoNext for index in range(3)]
			denominatorUndeformedLongEndoNext = denominatorUndeformedLongEndoNext + invDistLongEndoNext
			numeratorDeformedLongEndoNext = [numeratorDeformedLongEndoNext[index] + deformedLongEndoNext[index] * invDistLongEndoNext for index in range(3)]
			denominatorDeformedLongEndoNext = denominatorDeformedLongEndoNext + invDistLongEndoNext

			numeratorUndeformedLongMid = [numeratorUndeformedLongMid[index] + undeformedLongMid[index] * invDistLongMid for index in range(3)]
			denominatorUndeformedLongMid = denominatorUndeformedLongMid + invDistLongMid
			numeratorDeformedLongMid = [numeratorDeformedLongMid[index] + deformedLongMid[index] * invDistLongMid for index in range(3)]
			denominatorDeformedLongMid = denominatorDeformedLongMid + invDistLongMid
			numeratorUndeformedLongMidNext = [numeratorUndeformedLongMidNext[index] + undeformedLongMidNext[index] * invDistLongMidNext for index in range(3)]
			denominatorUndeformedLongMidNext = denominatorUndeformedLongMidNext + invDistLongMidNext
			numeratorDeformedLongMidNext = [numeratorDeformedLongMidNext[index] + deformedLongMidNext[index] * invDistLongMidNext for index in range(3)]
			denominatorDeformedLongMidNext = denominatorDeformedLongMidNext + invDistLongMidNext

			numeratorUndeformedLongEpi = [numeratorUndeformedLongEpi[index] + undeformedLongEpi[index] * invDistLongEpi for index in range(3)]
			denominatorUndeformedLongEpi = denominatorUndeformedLongEpi + invDistLongEpi
			numeratorDeformedLongEpi = [numeratorDeformedLongEpi[index] + deformedLongEpi[index] * invDistLongEpi for index in range(3)]
			denominatorDeformedLongEpi = denominatorDeformedLongEpi + invDistLongEpi
			numeratorUndeformedLongEpiNext = [numeratorUndeformedLongEpiNext[index] + undeformedLongEpiNext[index] * invDistLongEpiNext for index in range(3)]
			denominatorUndeformedLongEpiNext = denominatorUndeformedLongEpiNext + invDistLongEpiNext
			numeratorDeformedLongEpiNext = [numeratorDeformedLongEpiNext[index] + deformedLongEpiNext[index] * invDistLongEpiNext for index in range(3)]
			denominatorDeformedLongEpiNext = denominatorDeformedLongEpiNext + invDistLongEpiNext

			numeratorUndeformedCircEndo = [numeratorUndeformedCircEndo[index] + undeformedCircumferentialEndo[index] * invDistCircumferentialEndo for index in range(3)]
			denominatorUndeformedCircEndo = denominatorUndeformedCircEndo + invDistCircumferentialEndo
			numeratorDeformedCircEndo = [numeratorDeformedCircEndo[index] + deformedCircumferentialEndo[index] * invDistCircumferentialEndo for index in range(3)]
			denominatorDeformedCircEndo = denominatorDeformedCircEndo + invDistCircumferentialEndo

			numeratorUndeformedCircMid = [numeratorUndeformedCircMid[index] + undeformedCircumferentialMid[index] * invDistCircumferentialMid for index in range(3)]
			denominatorUndeformedCircMid = denominatorUndeformedCircMid + invDistCircumferentialMid
			numeratorDeformedCircMid = [numeratorDeformedCircMid[index] + deformedCircumferentialMid[index] * invDistCircumferentialMid for index in range(3)]
			denominatorDeformedCircMid = denominatorDeformedCircMid + invDistCircumferentialMid

			numeratorUndeformedCircEpi = [numeratorUndeformedCircEpi[index] + undeformedCircumferentialEpi[index] * invDistCircumferentialEpi for index in range(3)]
			denominatorUndeformedCircEpi = denominatorUndeformedCircEpi + invDistCircumferentialEpi
			numeratorDeformedCircEpi = [numeratorDeformedCircEpi[index] + deformedCircumferentialEpi[index] * invDistCircumferentialEpi for index in range(3)]
			denominatorDeformedCircEpi = denominatorDeformedCircEpi + invDistCircumferentialEpi


		undeformedLongEndo = [x / denominatorUndeformedLongEndo for x in numeratorUndeformedLongEndo]
		deformedLongEndo = [x / denominatorDeformedLongEndo for x in numeratorDeformedLongEndo]
		undeformedLongEndoNext = [x / denominatorUndeformedLongEndoNext for x in numeratorUndeformedLongEndoNext]
		deformedLongEndoNext = [x / denominatorDeformedLongEndoNext for x in numeratorDeformedLongEndoNext]

		undeformedLongMid = [x / denominatorUndeformedLongMid for x in numeratorUndeformedLongMid]
#		print undeformedLongMid
#		print "next \n"
		deformedLongMid = [x / denominatorDeformedLongMid for x in numeratorDeformedLongMid]
		undeformedLongMidNext = [x / denominatorUndeformedLongMidNext for x in numeratorUndeformedLongMidNext]
		deformedLongMidNext = [x / denominatorDeformedLongMidNext for x in numeratorDeformedLongMidNext]

		undeformedLongEpi = [x / denominatorUndeformedLongEpi for x in numeratorUndeformedLongEpi]
		deformedLongEpi = [x / denominatorDeformedLongEpi for x in numeratorDeformedLongEpi]
		undeformedLongEpiNext = [x / denominatorUndeformedLongEpiNext for x in numeratorUndeformedLongEpiNext]
		deformedLongEpiNext = [x / denominatorDeformedLongEpiNext for x in numeratorDeformedLongEpiNext]

		undeformedCircumferentialEndo = [x / denominatorUndeformedCircEndo for x in numeratorUndeformedCircEndo]
		deformedCircumferentialEndo = [x / denominatorDeformedCircEndo for x in numeratorDeformedCircEndo]

		undeformedCircumferentialMid = [x / denominatorUndeformedCircMid for x in numeratorUndeformedCircMid]
		deformedCircumferentialMid = [x / denominatorDeformedCircMid for x in numeratorDeformedCircMid]

		undeformedCircumferentialEpi = [x / denominatorUndeformedCircEpi for x in numeratorUndeformedCircEpi]
		deformedCircumferentialEpi = [x / denominatorDeformedCircEpi for x in numeratorDeformedCircEpi]

		#interpolating all values of the particles in the point cloud next to the point of interest by weighing them inversally proportional to the distance to the point of interest.
		#using analytical values stored in:
		#listOfPointsLongitudinalEndo = []
		#listOfPointsLongitudinalMid = []
		#listOfPointsLongitudinalEpi = []
		#listOfPointsCircumferentialEndo = []
		#listOfPointsCircumferentialMid = []
		#listOfPointsCircumferentialEpi = []

#		LEndo_0 = distanceBetweenPoints(listOfPointsLongitudinalEndo[idsIndex + 1], listOfPointsLongitudinalEndo[idsIndex])
#		LMid_0 = distanceBetweenPoints(listOfPointsLongitudinalMid[idsIndex + 1], listOfPointsLongitudinalMid[idsIndex])
#		LEpi_0 = distanceBetweenPoints(listOfPointsLongitudinalEpi[idsIndex + 1], listOfPointsLongitudinalEpi[idsIndex])

#		LEndo_f = distanceBetweenPoints(deformedLongEndoNext, deformedLongEndo)
#		LMid_f = distanceBetweenPoints(deformedLongMidNext, deformedLongMid)
#		LEpi_f = distanceBetweenPoints(deformedLongEpiNext, deformedLongEpi)

#		CEndo_0 = distanceBetweenPoints(listOfPointsCircumferentialEndo[idsIndex], listOfPointsLongitudinalEndo[idsIndex])
#		CMid_0 = distanceBetweenPoints(listOfPointsCircumferentialMid[idsIndex], listOfPointsLongitudinalMid[idsIndex])
#		CEpi_0 = distanceBetweenPoints(listOfPointsCircumferentialEpi[idsIndex], listOfPointsLongitudinalEpi[idsIndex])

#		CEndo_f = distanceBetweenPoints(deformedCircumferentialEndo, deformedLongEndo)
#		CMid_f = distanceBetweenPoints(deformedCircumferentialMid, deformedLongMid)
#		CEpi_f = distanceBetweenPoints(deformedCircumferentialEpi, deformedLongEpi)
#For transmural strains, we use pairs of neighbouring endocardium midwall, midwall epicardium and endocardium epicardium points to calculate radial strain at endocardium, epicardium and midwall.
#		TEndo_0 = distanceBetweenPoints(listOfPointsLongitudinalEndo[idsIndex], listOfPointsLongitudinalMid[idsIndex])
#		TMid_0 = distanceBetweenPoints(listOfPointsLongitudinalMid[idsIndex], listOfPointsLongitudinalEpi[idsIndex])
#		TEpi_0 = distanceBetweenPoints(listOfPointsLongitudinalEndo[idsIndex], listOfPointsLongitudinalEpi[idsIndex])

#		TEndo_f = distanceBetweenPoints(deformedLongEndo, deformedLongMid)
#		TMid_f = distanceBetweenPoints(deformedLongMid, deformedLongEpi)
#		TEpi_f = distanceBetweenPoints(deformedLongEndo, deformedLongEpi)

#		numericalStrainLongEndo = (LEndo_f - LEndo_0) / LEndo_0  * 100
#		print "analytical undeformed point", LEndo_0, "analytical strain ", (LEndo_f - LEndo_0) / LEndo_0  * 100, " %"
		LEndo_0 = distanceBetweenPoints(undeformedLongEndoNext, undeformedLongEndo)
#		print "numerical undeformed point", LEndo_0, "numerical strain ", (LEndo_f - LEndo_0) / LEndo_0  * 100, " %"
		LMid_0 = distanceBetweenPoints(undeformedLongMidNext, undeformedLongMid)
		LEpi_0 = distanceBetweenPoints(undeformedLongEpiNext, undeformedLongEpi)

		LEndo_f = distanceBetweenPoints(deformedLongEndoNext, deformedLongEndo)
		LMid_f = distanceBetweenPoints(deformedLongMidNext, deformedLongMid)
		LEpi_f = distanceBetweenPoints(deformedLongEpiNext, deformedLongEpi)

		CEndo_0 = distanceBetweenPoints(undeformedCircumferentialEndo, undeformedLongEndo)
		CMid_0 = distanceBetweenPoints(undeformedCircumferentialMid, undeformedLongMid)
		CEpi_0 = distanceBetweenPoints(undeformedCircumferentialEpi, undeformedLongEpi)

		CEndo_f = distanceBetweenPoints(deformedCircumferentialEndo, deformedLongEndo)
		CMid_f = distanceBetweenPoints(deformedCircumferentialMid, deformedLongMid)
		CEpi_f = distanceBetweenPoints(deformedCircumferentialEpi, deformedLongEpi)
#For transmural strains, we use pairs of neighbouring endocardium midwall, midwall epicardium and endocardium epicardium points to calculate radial strain at endocardium, epicardium and midwall.
		TEndo_0 = distanceBetweenPoints(undeformedLongEndo, undeformedLongMid)
		TMid_0 = distanceBetweenPoints(undeformedLongMid, undeformedLongEpi)
		TEpi_0 = distanceBetweenPoints(undeformedLongEndo, undeformedLongEpi)

		TEndo_f = distanceBetweenPoints(deformedLongEndo, deformedLongMid)
		TMid_f = distanceBetweenPoints(deformedLongMid, deformedLongEpi)
		TEpi_f = distanceBetweenPoints(deformedLongEndo, deformedLongEpi)

#		print undeformedLongEndo
		strainLongEndo = (LEndo_f - LEndo_0) / LEndo_0  * 100
		strainLongMid = (LMid_f - LMid_0) / LMid_0 * 100
		strainLongEpi = (LEpi_f - LEpi_0) / LEpi_0 * 100
		strainCircEndo = (CEndo_f - CEndo_0) / CEndo_0  * 100
		strainCircMid = (CMid_f - CMid_0) / CMid_0 * 100
		strainCircEpi = (CEpi_f - CEpi_0) / CEpi_0 * 100
		strainTransEndo = (TEndo_f - TEndo_0) / TEndo_0  * 100
		strainTransMid = (TMid_f - TMid_0) / TMid_0 * 100
		strainTransEpi = (TEpi_f - TEpi_0) / TEpi_0 * 100

		print LMid_f, LMid_0, strainLongMid, "%", undeformedLongMid, undeformedLongMidNext
		listUndeformedLongEndoX.append(undeformedLongMidNext[0])
		listUndeformedLongEndoZ.append(undeformedLongMidNext[2])
		listDeformedLongEndoX.append(deformedLongMidNext[0])
		listDeformedLongEndoZ.append(deformedLongMidNext[2])

		listNumericalUndeformedLongEndoX.append(undeformedLongEndo[0])
		listNumericalUndeformedLongEndoZ.append(undeformedLongEndo[2])

		listOfStrainsLongEndo.append(strainLongEndo)
#		listOfAnalyticalComparedStrainsLongEndo.append(numericalStrainLongEndo)
		listOfStrainsLongMid.append(strainLongMid)
		listOfStrainsLongEpi.append(strainLongEpi)
		listOfStrainsCircEndo.append(strainCircEndo)
		listOfStrainsCircMid.append(strainCircMid)
		listOfStrainsCircEpi.append(strainCircEpi)
		listOfStrainsTransEndo.append(strainTransEndo)
		listOfStrainsTransMid.append(strainTransMid)
		listOfStrainsTransEpi.append(strainTransEpi)

	listOfAllStrains.append(listOfStrainsLongEndo)
	listOfAllStrains.append(listOfStrainsLongMid)
	listOfAllStrains.append(listOfStrainsLongEpi)
	listOfAllStrains.append(listOfStrainsCircEndo)
	listOfAllStrains.append(listOfStrainsCircMid)
	listOfAllStrains.append(listOfStrainsCircEpi)
	listOfAllStrains.append(listOfStrainsTransEndo)
	listOfAllStrains.append(listOfStrainsTransMid)
	listOfAllStrains.append(listOfStrainsTransEpi)
	return listOfAllStrains
#	print len(listOfPointsLongitudinalMid)
#	listOfPointsLongitudinalMid.pop()
#	x = listOfPointsLongitudinalMid[:len(listOfPointsLongitudinalMid)-1]
#	if len(listOfPointsLongitudinalMid) > 0:
#        	del listOfPointsLongitudinalMid[-1]
#	print len(listOfPointsLongitudinalMid)

#	labelsPoints = ["p1","p2","p3","p4","p5","p6","p7","p8","p9"]
#	x = listDeformedLongEndoX
#	z = listDeformedLongEndoZ
#	plt.plot(x, z, 'ro', ms = 10)
#	plt.plot(x, z)
	# You can specify a rotation for the tick labels in degrees or with keywords.
#	plt.xticks(x, labelsPoints, rotation='vertical')
	# Pad margins so that markers don't get clipped by the axes
#	plt.margins(0.2)
	# Tweak spacing to prevent clipping of tick-labels
#	plt.subplots_adjust(bottom=0.15)
#	plt.show()

#	x = listUndeformedLongEndoX
#	z = listUndeformedLongEndoZ
#	plt.plot(x, z)
#	plt.plot(x, z, 'ro', ms = 10)
	# You can specify a rotation for the tick labels in degrees or with keywords.
#	plt.xticks(x, labelsPoints, rotation='vertical')
	# Pad margins so that markers don't get clipped by the axes
#	plt.margins(0.2)
	# Tweak spacing to prevent clipping of tick-labels
#	plt.subplots_adjust(bottom=0.15)


#	x = listDeformedLongEndoX
#	z = listDeformedLongEndoZ
#	plt.plot(x, z)
#	plt.plot(x, z, 'ro', ms = 10)
	# You can specify a rotation for the tick labels in degrees or with keywords.
#	plt.xticks(x, rotation='vertical')
	# Pad margins so that markers don't get clipped by the axes
#	plt.margins(0.2)
	# Tweak spacing to prevent clipping of tick-labels
#	plt.subplots_adjust(bottom=0.15)

#	x = listUndeformedLongEndoX
#	z = listUndeformedLongEndoZ
#	plt.plot(x, z)
#	plt.plot(x, z, 'ro', ms = 10)
	# You can specify a rotation for the tick labels in degrees or with keywords.
#	plt.xticks(x, rotation='vertical')
	# Pad margins so that markers don't get clipped by the axes
#	plt.margins(0.2)
	# Tweak spacing to prevent clipping of tick-labels
#	plt.subplots_adjust(bottom=0.15)


#for strains FEM:
def addingSubPlots(undeformedFileNameFolder):
	#we need to create the deformedFileName from the undeformedFileNameFolder
#	deformedFileName = undeformedFileNameFolder.join('deformed.vtk')
	undeformed = "/undeformed.vtk"
	undeformedFileName = undeformedFileNameFolder + undeformed
	deformed = "/deformed.vtk"
	deformedFileName = undeformedFileNameFolder + deformed

	reader = vtkUnstructuredGridReader();
	reader.SetFileName(undeformedFileName);
	reader.Update();
	undeformedMesh = reader.GetOutput()


	readerSmall = vtkPolyDataReader();
	#readerSmall.SetFileName('/home/eric/dataForSimulations/sander314-mechbench-f541b4f0e0b3/problem 2/CardioMechanics/139587/undeformedPolyData1.vtk');
	#readerSmall.SetFileName('/home/eric/dataForSimulations/sander314-mechbench-f541b4f0e0b3/problem 2/Simula-FEniCS/8444/undeformedPolyData.vtk');
	#readerSmall.SetFileName('/home/eric/dataForSimulations/sander314-mechbench-f541b4f0e0b3/problem 2/glasgowHeart-IBFE/40x4x4(structure)_120x42x90(fluid)/underofmedPolyData.vtk');
	#readerSmall.Update();
	#m = readerSmall.GetOutput()
	#p = m.GetPoint(1)
	#print(m.GetNumberOfPoints())
	#cl = vmtkcenterlines(m, [0, 0.5, 0.5], [10, 0.5, 0.5])
	#print(cl.GetNumberOfPoints())
	#writer1 = vtkPolyDataWriter()
	#writer1.SetFileName('outputUnstructured.vtk')
	#writer1.SetInput(cl)
	#writer1.Write()

	listOfPoints = []
	for i in range(11):
		listOfPoints.append([i, 0.5, 0.5])

	listOfPointsUp = []
	for i in range(11):
		listOfPointsUp.append([i, 0.5, 1])

	listOfPointsSide = []
	for i in range(11):
		listOfPointsSide.append([i, 1, 0.5])



	# initiate point locator
	locator = vtk.vtkPointLocator()
	locator.SetDataSet(undeformedMesh)
	locator.BuildLocator()

	listOfFoundIdsInOrder = []
	listOfFoundIdsInOrderUp = []
	listOfFoundIdsInOrderSide = []

	for i in range(11):
		closestpoint_id = locator.FindClosestPoint(listOfPoints[i])
		listOfFoundIdsInOrder.append(closestpoint_id)
		print closestpoint_id#we just print the closespoint_id

	for i in range(11):
		closestpoint_id = locator.FindClosestPoint(listOfPointsUp[i])
		listOfFoundIdsInOrderUp.append(closestpoint_id)

	for i in range(11):
		closestpoint_id = locator.FindClosestPoint(listOfPointsSide[i])
		listOfFoundIdsInOrderSide.append(closestpoint_id)

	readerDeformed = vtkUnstructuredGridReader();
	readerDeformed.SetFileName(deformedFileName);
	readerDeformed.Update();
	deformedMesh = readerDeformed.GetOutput()

#	def distanceBetweenPoints(p1, p2):
#		return sqrt( (p1[0] - p2[0]) * (p1[0] - p2[0])
#	+ (p1[1] - p2[1]) * (p1[1] - p2[1])
#	+ (p1[2] - p2[2]) * (p1[2] - p2[2]))



	#arSmall = undeformedMesh.GetPointIds()
	#print arSmall
	print "STRAINS"
	listOfStrainsInX = []
	listOfStrainsInY = []
	listOfStrainsInZ = []

	for i in range(10):
		undeformedPoint = undeformedMesh.GetPoint(listOfFoundIdsInOrder[i])
		undeformedPointNext = undeformedMesh.GetPoint(listOfFoundIdsInOrder[i+1])
		undeformedPointUp = undeformedMesh.GetPoint(listOfFoundIdsInOrderUp[i])
		undeformedPointSide = undeformedMesh.GetPoint(listOfFoundIdsInOrderSide[i])

		deformedPoint = deformedMesh.GetPoint(listOfFoundIdsInOrder[i])
		deformedPointNext = deformedMesh.GetPoint(listOfFoundIdsInOrder[i+1])
		deformedPointUp = deformedMesh.GetPoint(listOfFoundIdsInOrderUp[i])
		deformedPointSide = deformedMesh.GetPoint(listOfFoundIdsInOrderSide[i])

		L0x = distanceBetweenPoints(undeformedPointNext, undeformedPoint)
		L0y = distanceBetweenPoints(undeformedPointSide, undeformedPoint)
		L0z = distanceBetweenPoints(undeformedPointUp, undeformedPoint)

		Lfx = distanceBetweenPoints(deformedPointNext, deformedPoint)
		Lfy = distanceBetweenPoints(deformedPointSide, deformedPoint)
		Lfz = distanceBetweenPoints(deformedPointUp, deformedPoint)

	#	Lfx = sqrt( (deformedPointNext[0] - deformedPoint[0]) * (deformedPointNext[0] - deformedPoint[0])
	#+ (deformedPointNext[1] - deformedPoint[1]) * (deformedPointNext[1] - deformedPoint[1])
	#+ (deformedPointNext[2] - deformedPoint[2]) * (deformedPointNext[2] - deformedPoint[2]))
	#	Lfy = (deformedPointSide[1] - deformedPoint[1])
	#	Lfz = (deformedPointUp[2] - deformedPoint[2])

	#	deformedPointSide[2] - undeformedPoint[2]

		#print L0x, " and deformed ", Lfx, Lfx-L0x
#		if(L0x != 0):
		strainInX = (Lfx - L0x) * 100
#		else:
#			strainInX = 0
#		if(L0y != 0):
		strainInY = (Lfy - L0y) * 100
#		else:
#			strainInY = 0
#		if(L0z != 0):
		strainInZ = (Lfz - L0z) * 100
#		else:
#			strainInZ = 0
		print "p",i,"\t",strainInX,"%",strainInY,"%",strainInZ,"%"
		listOfStrainsInX.append(strainInX)
		listOfStrainsInY.append(strainInY)
		listOfStrainsInZ.append(strainInZ)


	print "Deformed Location of Line (x, 0.5, 0.5)"
	deformedPointsX = []
	deformedPointsZ = []
	for i in range(11):
		deformedPoint = deformedMesh.GetPoint(listOfFoundIdsInOrder[i])
		deformedPointsX.append(deformedPoint[0])
		deformedPointsZ.append(deformedPoint[2])
#		print deformedPoint[0], deformedPoint[2]

	x = [1, 2, 3, 4, 5, 6, 7, 8, 9]
	labelsPoints = ["p1","p2","p3","p4","p5","p6","p7","p8","p9"]
	plt.plot(x, listOfStrainsInX, 'ro', ms = 10)
	plt.plot(x, listOfStrainsInX)
	# You can specify a rotation for the tick labels in degrees or with keywords.
	plt.xticks(x, labelsPoints, rotation='vertical')
	# Pad margins so that markers don't get clipped by the axes
	plt.margins(0.2)
	# Tweak spacing to prevent clipping of tick-labels
	plt.subplots_adjust(bottom=0.15)
#	plt.show()




def addingSubPlotFEMStrain(undeformedFileNameFolder):
	undeformed = "/undeformed.vtk"
	undeformedFileName = undeformedFileNameFolder + undeformed
	deformed = "/deformed.vtk"
	deformedFileName = undeformedFileNameFolder + deformed

	reader = vtkUnstructuredGridReader();
	reader.SetFileName(undeformedFileName);
	reader.Update();
	undeformedMesh = reader.GetOutput()
	readerSmall = vtkPolyDataReader();
	readerDeformed = vtkUnstructuredGridReader();
	readerDeformed.SetFileName(deformedFileName);
	readerDeformed.Update();
	deformedMesh = readerDeformed.GetOutput()


	listOfFoundIdsInOrder = []
	listOfFoundIdsInOrderSecondPlot = []
	deformedPointsZHighResolution = []
	listOfPointsSecondPlot = []
	xSecondPlot = []

	resolutionOfSecondPlot = 9

	listOfPointsLongitudinalEndo = []
	listOfPointsLongitudinalMid = []
	listOfPointsLongitudinalEpi = []
	listOfPointsCircumferentialEndo = []
	listOfPointsCircumferentialMid = []
	listOfPointsCircumferentialEpi = []

	n_u = 10
	u1 = -pi
	tEndo = 0.1
	tMid = 0.5
	tEpi = 0.9
	rsEndo = 7
	rlEndo = 17
	rsEpi = 10
	rlEpi = 20
	rsMid = rsEndo + 0.5 * (rsEpi - rsEndo)
	rlMid = rlEndo + 0.5 * (rlEpi - rlEndo)
	
	print rsMid, rlMid

	for i in range(11):
		v = 0
		u2 = -acos(5. / (17 + 3 * tEndo))
		u = u1 + (u2 - u1) / n_u * (i + 1) * 0.95
		x = rsEndo * sin(u)
		y = 0
		z = rlEndo * cos(u)
		listOfPointsLongitudinalEndo.append([x, y, z])
		print x,y,z

		u2 = -acos(5. / (17 + 3 * tMid))
		u = u1 + (u2 - u1) / n_u * (i + 1) * 0.95
		x = rsMid * sin(u)
		y = 0
		z = rlMid * cos(u)
		listOfPointsLongitudinalMid.append([x, y, z])

		u2 = -acos(5. / (17 + 3 * tEpi))
		u = u1 + (u2 - u1) / n_u * (i + 1) * 0.95
		x = rsEpi * sin(u)
		y = 0
		z = rlEpi * cos(u)
		listOfPointsLongitudinalEpi.append([x, y, z])

		v = math.pi / 10.
		u2 = -acos(5. / (17 + 3 * tEndo))
		u = u1 + (u2 - u1) / n_u * (i + 1) * 0.95
		x = rsEndo * sin(u) * cos(v)
		y = rsEndo * sin(u) * sin(v)
		z = rlEndo * cos(u)
		listOfPointsCircumferentialEndo.append([x, y, z])

		u2 = -acos(5. / (17 + 3 * tMid))
		u = u1 + (u2 - u1) / n_u * (i + 1) * 0.95
		x = rsMid * sin(u) * cos(v)
		y = rsMid * sin(u) * sin(v)
		z = rlMid * cos(u)
		listOfPointsCircumferentialMid.append([x, y, z])

		u2 = -acos(5. / (17 + 3 * tEpi))
		u = u1 + (u2 - u1) / n_u * (i + 1) * 0.95
		x = rsEpi * sin(u) * cos(v)
		y = rsEpi * sin(u) * sin(v)
		z = rlEpi * cos(u)
		listOfPointsCircumferentialEpi.append([x, y, z])
	listOfFoundIdsLongitudinalEndo = []
	listOfFoundIdsLongitudinalMid = []
	listOfFoundIdsLongitudinalEpi = []
	listOfFoundIdsCircumferentialEndo = []
	listOfFoundIdsCircumferentialMid = []
	listOfFoundIdsCircumferentialEpi = []

	locator = vtk.vtkPointLocator()
	locator.SetDataSet(undeformedMesh)
	locator.BuildLocator()
	print "strains_start"
	for i in range(11):
		closestpoint_id = locator.FindClosestPoint(listOfPointsLongitudinalEndo[i])
		listOfFoundIdsLongitudinalEndo.append(closestpoint_id)
		print closestpoint_id
		closestpoint_id = locator.FindClosestPoint(listOfPointsLongitudinalMid[i])
		listOfFoundIdsLongitudinalMid.append(closestpoint_id)
		closestpoint_id = locator.FindClosestPoint(listOfPointsLongitudinalEpi[i])
		listOfFoundIdsLongitudinalEpi.append(closestpoint_id)
		closestpoint_id = locator.FindClosestPoint(listOfPointsCircumferentialEndo[i])
		listOfFoundIdsCircumferentialEndo.append(closestpoint_id)
		closestpoint_id = locator.FindClosestPoint(listOfPointsCircumferentialMid[i])
		listOfFoundIdsCircumferentialMid.append(closestpoint_id)
		closestpoint_id = locator.FindClosestPoint(listOfPointsCircumferentialEpi[i])
		listOfFoundIdsCircumferentialEpi.append(closestpoint_id)
	print "STRAINS"
	circumferentialStrainEndo = []
	longitudinalStrainEndo = []
	transmuralStrainEndo = []
	circumferentialStrainMid = []
	longitudinalStrainMid = []
	transmuralStrainMid = []
	circumferentialStrainEpi = []
	longitudinalStrainEpi = []
	transmuralStrainEpi = []
	listOfStrainsLongEndo = []
	listOfAnalyticalComparedStrainsLongEndo = []
	listOfStrainsLongMid = []
	listOfStrainsLongEpi = []
	listOfStrainsCircEndo = []
	listOfStrainsCircMid = []
	listOfStrainsCircEpi = []
	listOfStrainsTransEndo = []
	listOfStrainsTransMid = []
	listOfStrainsTransEpi = []

	listOfStrainsForEachInterpolationPoint = []
	listOfAllStrains = []
	

	for i in range(10):
		undeformedLongEndo = undeformedMesh.GetPoint(listOfFoundIdsLongitudinalEndo[i])
		deformedLongEndo = deformedMesh.GetPoint(listOfFoundIdsLongitudinalEndo[i])
		undeformedLongEndoNext = undeformedMesh.GetPoint(listOfFoundIdsLongitudinalEndo[i + 1])
		deformedLongEndoNext = deformedMesh.GetPoint(listOfFoundIdsLongitudinalEndo[i + 1])

		undeformedLongMid = undeformedMesh.GetPoint(listOfFoundIdsLongitudinalMid[i])
		deformedLongMid = deformedMesh.GetPoint(listOfFoundIdsLongitudinalMid[i])
		undeformedLongMidNext = undeformedMesh.GetPoint(listOfFoundIdsLongitudinalMid[i + 1])
		deformedLongMidNext = deformedMesh.GetPoint(listOfFoundIdsLongitudinalMid[i + 1])

		undeformedLongEpi = undeformedMesh.GetPoint(listOfFoundIdsLongitudinalEpi[i])
		deformedLongEpi = deformedMesh.GetPoint(listOfFoundIdsLongitudinalEpi[i])
		undeformedLongEpiNext = undeformedMesh.GetPoint(listOfFoundIdsLongitudinalEpi[i + 1])
		deformedLongEpiNext = deformedMesh.GetPoint(listOfFoundIdsLongitudinalEpi[i + 1])

		undeformedCircumferentialEndo = undeformedMesh.GetPoint(listOfFoundIdsCircumferentialEndo[i])
		deformedCircumferentialEndo = deformedMesh.GetPoint(listOfFoundIdsCircumferentialEndo[i])
		undeformedCircumferentialMid = undeformedMesh.GetPoint(listOfFoundIdsCircumferentialMid[i])
		deformedCircumferentialMid = deformedMesh.GetPoint(listOfFoundIdsCircumferentialMid[i])
		undeformedCircumferentialEpi = undeformedMesh.GetPoint(listOfFoundIdsCircumferentialEpi[i])
		deformedCircumferentialEpi = deformedMesh.GetPoint(listOfFoundIdsCircumferentialEpi[i])

		LEndo_0 = distanceBetweenPoints(undeformedLongEndoNext, undeformedLongEndo)
		LMid_0 = distanceBetweenPoints(undeformedLongMidNext, undeformedLongMid)
		LEpi_0 = distanceBetweenPoints(undeformedLongEpiNext, undeformedLongEpi)

		LEndo_f = distanceBetweenPoints(deformedLongEndoNext, deformedLongEndo)
		LMid_f = distanceBetweenPoints(deformedLongMidNext, deformedLongMid)
		LEpi_f = distanceBetweenPoints(deformedLongEpiNext, deformedLongEpi)

		CEndo_0 = distanceBetweenPoints(undeformedCircumferentialEndo, undeformedLongEndo)
		CMid_0 = distanceBetweenPoints(undeformedCircumferentialMid, undeformedLongMid)
		CEpi_0 = distanceBetweenPoints(undeformedCircumferentialEpi, undeformedLongEpi)

		CEndo_f = distanceBetweenPoints(deformedCircumferentialEndo, deformedLongEndo)
		CMid_f = distanceBetweenPoints(deformedCircumferentialMid, deformedLongMid)
		CEpi_f = distanceBetweenPoints(deformedCircumferentialEpi, deformedLongEpi)
		TEndo_0 = distanceBetweenPoints(undeformedLongEndo, undeformedLongMid)
		TMid_0 = distanceBetweenPoints(undeformedLongMid, undeformedLongEpi)
		TEpi_0 = distanceBetweenPoints(undeformedLongEndo, undeformedLongEpi)

		TEndo_f = distanceBetweenPoints(deformedLongEndo, deformedLongMid)
		TMid_f = distanceBetweenPoints(deformedLongMid, deformedLongEpi)
		TEpi_f = distanceBetweenPoints(deformedLongEndo, deformedLongEpi)

		strainLongEndo = (LEndo_f - LEndo_0) / LEndo_0  * 100
		strainLongMid = (LMid_f - LMid_0) / LMid_0 * 100
		strainLongEpi = (LEpi_f - LEpi_0) / LEpi_0 * 100
		print undeformedCircumferentialEndo, undeformedLongEndo
		strainCircEndo = (CEndo_f - CEndo_0) / CEndo_0  * 100
		
		strainCircMid = (CMid_f - CMid_0) / CMid_0 * 100
		strainCircEpi = (CEpi_f - CEpi_0) / CEpi_0 * 100
		strainTransEndo = (TEndo_f - TEndo_0) / TEndo_0  * 100
		strainTransMid = (TMid_f - TMid_0) / TMid_0 * 100
		strainTransEpi = (TEpi_f - TEpi_0) / TEpi_0 * 100

		print LMid_f, LMid_0, strainLongMid, "%", undeformedLongMid, undeformedLongMidNext

		listOfStrainsLongEndo.append(strainLongEndo)
#		listOfAnalyticalComparedStrainsLongEndo.append(numericalStrainLongEndo)
		listOfStrainsLongMid.append(strainLongMid)
		listOfStrainsLongEpi.append(strainLongEpi)
		listOfStrainsCircEndo.append(strainCircEndo)
		listOfStrainsCircMid.append(strainCircMid)
		listOfStrainsCircEpi.append(strainCircEpi)
		listOfStrainsTransEndo.append(strainTransEndo)
		listOfStrainsTransMid.append(strainTransMid)
		listOfStrainsTransEpi.append(strainTransEpi)

	listOfAllStrains.append(listOfStrainsLongEndo)
	listOfAllStrains.append(listOfStrainsLongMid)
	listOfAllStrains.append(listOfStrainsLongEpi)
	listOfAllStrains.append(listOfStrainsCircEndo)
	listOfAllStrains.append(listOfStrainsCircMid)
	listOfAllStrains.append(listOfStrainsCircEpi)
	listOfAllStrains.append(listOfStrainsTransEndo)
	listOfAllStrains.append(listOfStrainsTransMid)
	listOfAllStrains.append(listOfStrainsTransEpi)

	return listOfAllStrains

#	x = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
#	labelsPoints = ["p0","p1","p2","p3","p4","p5","p6","p7","p8","p9"]
#	plt.plot(x, listOfStrainsLongEndo, 'bo')
#	plt.plot(x, listOfStrainsLongEndo)
#	plt.xticks(x, labelsPoints, rotation='vertical')
#	plt.margins(0.2)
#	plt.subplots_adjust(bottom=0.15)

#	plt.show()



#	x = [1, 2, 3, 4, 5, 6, 7, 8, 9]
#	labelsPoints = ["p1","p2","p3","p4","p5","p6","p7","p8","p9"]
#	plt.plot(x, listOfStrainsLongMid, 'ro', ms = 10)
#	plt.plot(x, listOfStrainsLongMid)
	# You can specify a rotation for the tick labels in degrees or with keywords.
#	plt.xticks(x, labelsPoints, rotation='vertical')
	# Pad margins so that markers don't get clipped by the axes
#	plt.margins(0.2)
	# Tweak spacing to prevent clipping of tick-labels
#	plt.subplots_adjust(bottom=0.15)

#	plt.show()
#	plt.plot(x, listOfStrainsLongEpi, 'bo')
#	plt.plot(x, listOfStrainsLongEpi)
	# You can specify a rotation for the tick labels in degrees or with keywords.
#	plt.xticks(x, labelsPoints, rotation='vertical')
	# Pad margins so that markers don't get clipped by the axes
#	plt.margins(0.2)
	# Tweak spacing to prevent clipping of tick-labels
#	plt.subplots_adjust(bottom=0.15)
#	plt.show()

#	plt.show()
#	plt.plot(x, listOfStrainsCircEndo, 'ro', ms = 10)
#	plt.plot(x, listOfStrainsCircEndo)
#	# You can specify a rotation for the tick labels in degrees or with keywords.
#	plt.xticks(x, labelsPoints, rotation='vertical')
	# Pad margins so that markers don't get clipped by the axes
#	plt.margins(0.2)
	# Tweak spacing to prevent clipping of tick-labels
#	plt.subplots_adjust(bottom=0.15)
#	plt.show()

#	plt.plot(x, listOfStrainsCircMid, 'ro', ms = 10)
#	plt.plot(x, listOfStrainsCircMid)
	# You can specify a rotation for the tick labels in degrees or with keywords.
#	plt.xticks(x, labelsPoints, rotation='vertical')
	# Pad margins so that markers don't get clipped by the axes
#	plt.margins(0.2)
	# Tweak spacing to prevent clipping of tick-labels
#	plt.subplots_adjust(bottom=0.15)

#	plt.show()
#	labelsPoints = ["p1","p2","p3","p4","p5","p6","p7","p8","p9"]
#	plt.plot(x, listOfStrainsCircEpi, 'ro', ms = 10)
#	plt.plot(x, listOfStrainsCircEpi)
	# You can specify a rotation for the tick labels in degrees or with keywords.
#	plt.xticks(x, labelsPoints, rotation='vertical')
	# Pad margins so that markers don't get clipped by the axes
#	plt.margins(0.2)
	# Tweak spacing to prevent clipping of tick-labels
#	plt.subplots_adjust(bottom=0.15)
#	plt.show()

#	labelsPoints = ["p1","p2","p3","p4","p5","p6","p7","p8","p9"]
#	plt.plot(x, listOfStrainsTransEndo, 'ro', ms = 10)
#	plt.plot(x, listOfStrainsTransEndo)
	# You can specify a rotation for the tick labels in degrees or with keywords.
#	plt.xticks(x, labelsPoints, rotation='vertical')
	# Pad margins so that markers don't get clipped by the axes
#	plt.margins(0.2)
	# Tweak spacing to prevent clipping of tick-labels
#	plt.subplots_adjust(bottom=0.15)
#	plt.show()

#	labelsPoints = ["p1","p2","p3","p4","p5","p6","p7","p8","p9"]
#	plt.plot(x, listOfStrainsTransMid, 'ro', ms = 10)
#	plt.plot(x, listOfStrainsTransMid)
	# You can specify a rotation for the tick labels in degrees or with keywords.
#	plt.xticks(x, labelsPoints, rotation='vertical')
	# Pad margins so that markers don't get clipped by the axes
#	plt.margins(0.2)
	# Tweak spacing to prevent clipping of tick-labels
#	plt.subplots_adjust(bottom=0.15)
#	plt.show()

#	labelsPoints = ["p1","p2","p3","p4","p5","p6","p7","p8","p9"]
#	plt.plot(x, listOfStrainsTransEpi, 'ro', ms = 10)
#	plt.plot(x, listOfStrainsTransEpi)
	# You can specify a rotation for the tick labels in degrees or with keywords.
#	plt.xticks(x, labelsPoints, rotation='vertical')
	# Pad margins so that markers don't get clipped by the axes
#	plt.margins(0.2)
	# Tweak spacing to prevent clipping of tick-labels
#	plt.subplots_adjust(bottom=0.15)
#	plt.show()



plot_lines = []					
for root, directories, filenames in os.walk('/home/eric/scriptsToHelpSPHCreateGeometryInteraceETC/centerline/meshesWithMoreDegreesOfFreedomProb2'):
	for filename in filenames:
		if(filename == 'undeformed.vtk'):
			undeformedFileNameFolder = os.path.join(root)
			print undeformedFileNameFolder
#			l1 = addingSubPlotsOfDeformation(undeformedFileNameFolder)
			l2,l1 = addingSubPlotFEM(pointsWithIds_FEM, deformedPointsWithIds_FEM)
			l2,l3 = addingSubPlotSPH(pointsWithIds_sph, deformedPointsWithIds_sph)

			plot_lines.append([l1,l2,l3])
			legend1=plt.legend(plot_lines[0],["FEM","undef","SPH"],loc=1)
			plt.gca().add_artist(legend1)
			plt.xlabel('x (mm)')
			plt.ylabel('z (mm)')
			plt.savefig("centerlineComparison", bbox_inches='tight', dpi = 300)
			plt.clf()
#			plt.show()
#plt.show()
#addingSubPlotSPHJustOnePoint(pointsWithIds_sph, deformedPointsWithIds_sph)

#plt.show()



plt.rc('xtick', labelsize=20) 
plt.rc('ytick', labelsize=20) 


#addingSubPlotSPHStrain(pointsWithIds_sph, deformedPointsWithIds_sph)
listOfStrainsFEM = []
listOfStrainsSPH = []
listOfLabels = ['LongEndo','LongMid','LongEpi','CircEndo','CircMid','CircEpi','TransEndo','TransMid','TransEpi']
for root, directories, filenames in os.walk('/home/eric/scriptsToHelpSPHCreateGeometryInteraceETC/centerline/meshesWithMoreDegreesOfFreedomProb2'):
	for filename in filenames:
		if(filename == 'undeformed.vtk'):
			undeformedFileNameFolder = os.path.join(root)
			print undeformedFileNameFolder
			listOfStrainsFEM = addingSubPlotSPHStrain(pointsWithIds_FEM, deformedPointsWithIds_FEM)
			listOfStrainsSPH = addingSubPlotSPHStrain(pointsWithIds_sph, deformedPointsWithIds_sph)
			for graphicToPlot in range(9):
				plot_lines = []
				x = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
				labelsPoints = ["p0","p1","p2","p3","p4","p5","p6","p7","p8","p9"]
				l1, = plt.plot(x, listOfStrainsFEM[graphicToPlot], '--y^', ms = 10, linewidth = 3)
#				l1, = plt.plot(x, listOfStrainsFEM[graphicToPlot], 'b--', linewidth = 3)
				plt.xticks(x, labelsPoints, rotation='vertical')
				plt.margins(0.2)
				plt.subplots_adjust(bottom=0.15)

#				plt.plot(x, listOfStrainsSPH[graphicToPlot], 'ro', ms = 10)
				l2, = plt.plot(x, listOfStrainsSPH[graphicToPlot], '-bo', ms = 10, linewidth = 3)
				plt.xticks(x, labelsPoints, rotation='vertical')
				plt.margins(0.2)
				plt.subplots_adjust(bottom=0.15)

				plot_lines.append([l1,l2])
				if(graphicToPlot < 3):
					legend1=plt.legend(plot_lines[0],["FEM","SPH"],loc=4)
				else:
					legend1=plt.legend(plot_lines[0],["FEM","SPH"],loc=1)
				plt.gca().add_artist(legend1)

				plt.savefig(listOfLabels[graphicToPlot], bbox_inches='tight', dpi = 300)
				plt.clf()
#				plt.show()

