from basefunctions import *
from vmtkfunctions import *
from vtk import *
from math import *
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
import csv

# Previous to this, two simulations of a deforming beam should have been performed using FEM and SPH.
# This file reads the input files of the SPH simulation (particles_init_SPH.csv particles_end_SPH.csv) corresponding to the initial and final configuration
# It has the same input for the FEM simulation (particles_init_FEM.csv and particles_end_FEM.csv).
# The centerline of both simulations is computed and plotted in a figured. Strains in x,y,z is computed as well.


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
	sorted_by_distance = sorted(listOfTupleWithInfo, key=lambda tup: tup[0])
	listToReturn = []
	for i in range(numPoints):
		listToReturn.append(sorted_by_distance[i])
	return listToReturn

#********** We just read the surface **********
def addingSubPlotsOfDeformation(undeformedFileNameFolder):
	undeformed = "/undeformed.vtk"
	undeformedFileName = undeformedFileNameFolder + undeformed
	deformed = "/deformed.vtk"
	deformedFileName = undeformedFileNameFolder + deformed

	reader = vtkUnstructuredGridReader()
	reader.SetFileName(undeformedFileName)
	reader.Update()
	undeformedMesh = reader.GetOutput()

	readerDeformed = vtkUnstructuredGridReader()
	readerDeformed.SetFileName(deformedFileName)
	readerDeformed.Update()
	deformedMesh = readerDeformed.GetOutput()

	locator = vtk.vtkPointLocator()
	locator.SetDataSet(undeformedMesh)
	locator.BuildLocator()

	listOfFoundIdsInOrder = []
	listOfFoundIdsInOrderSecondPlot = []
	deformedPointsZHighResolution = []
	listOfPointsSecondPlot = []
	xSecondPlot = []
	resolutionOfSecondPlot = 10
	for i in range(resolutionOfSecondPlot):
		listOfPointsSecondPlot.append([i * 10./resolutionOfSecondPlot, 0.5, 0.5])
		xSecondPlot.append(i * 10./resolutionOfSecondPlot);

	for i in range(resolutionOfSecondPlot):
		closestpoint_id = locator.FindClosestPoint(listOfPointsSecondPlot[i])
		listOfFoundIdsInOrderSecondPlot.append(closestpoint_id)

	for i in range(resolutionOfSecondPlot):
		deformedPointsZHighResolution.append(deformedMesh.GetPoint(listOfFoundIdsInOrderSecondPlot[i])[2])

	print len(xSecondPlot)
	print len(deformedPointsZHighResolution)
	plt.plot(xSecondPlot, deformedPointsZHighResolution)

#for strains SPH:
def addingSubPlotSPHStrain(listOfPointsWithIds, listOfDeformedPoints):
	listOfDeformedPointsX = []
	listOfDeformedPointsZ = []
	listOfPoints = []
	for i in range(11):
		listOfPoints.append([i, 0.5, 0.5])

	listOfPointsUp = []
	for i in range(11):
		listOfPointsUp.append([i, 0.5, 0.9])

	listOfPointsSide = []
	for i in range(11):
		listOfPointsSide.append([i, 0.9, 0.5])

	listOfFoundIdsInOrder = []
	listOfFoundIdsInOrderUp = []
	listOfFoundIdsInOrderSide = []
	numOfPointsToInterpolate = 200
	exponentOfInterpolation = 1

	for i in range(11):
		closestpoint_id_sph = findIdWithMinDistance(listOfPoints[i], listOfPointsWithIds, numOfPointsToInterpolate)
		listOfFoundIdsInOrder.append(closestpoint_id_sph)

	for i in range(11):
		closestpoint_id_sph = findIdWithMinDistance(listOfPointsUp[i], listOfPointsWithIds, numOfPointsToInterpolate)
		listOfFoundIdsInOrderUp.append(closestpoint_id_sph)

	for i in range(11):
		closestpoint_id_sph = findIdWithMinDistance(listOfPointsSide[i], listOfPointsWithIds, numOfPointsToInterpolate)
		listOfFoundIdsInOrderSide.append(closestpoint_id_sph)

	print "STRAINS"
	listOfStrainsInX = []
	listOfStrainsInY = []
	listOfStrainsInZ = []
	listOfAllStrains = []

	for idsIndex in range(10):

		numeratorUndeformedPoint = [0,0,0]
		numeratorDeformedPoint = [0,0,0]
		numeratorUndeformedPointNext = [0, 0, 0]
		numeratorDeformedPointNext = [0, 0, 0]
		numeratorUndeformedPointUp = [0,0,0]
		numeratorDeformedPointUp = [0,0,0]
		numeratorUndeformedPointSide = [0,0,0]
		numeratorDeformedPointSide = [0,0,0]
		sumOfDists = 0
		sumOfDistsNext = 0
		sumOfDistsUp = 0
		sumOfDistsSide = 0
		for nPoints in range(numOfPointsToInterpolate):
			closestId = listOfFoundIdsInOrder[idsIndex][nPoints][2]
			closestIdNext = listOfFoundIdsInOrder[idsIndex + 1][nPoints][2]
			closestIdUp = listOfFoundIdsInOrderUp[idsIndex][nPoints][2]
			closestIdSide = listOfFoundIdsInOrderSide[idsIndex][nPoints][2]

			inverseDist = pow(1. / max(listOfFoundIdsInOrder[idsIndex][nPoints][0], 1e-3), exponentOfInterpolation)
			inverseDistNext = pow(1. / max(listOfFoundIdsInOrder[idsIndex + 1][nPoints][0], 1e-3), exponentOfInterpolation)
			inverseDistUp = pow(1. / max(listOfFoundIdsInOrderUp[idsIndex][nPoints][0], 1e-3), exponentOfInterpolation)
			inverseDistSide = pow(1. / max(listOfFoundIdsInOrderSide[idsIndex][nPoints][0], 1e-3), exponentOfInterpolation)

			sumOfDists = sumOfDists + inverseDist
			sumOfDistsNext = sumOfDistsNext + inverseDistNext
			sumOfDistsUp = sumOfDistsUp + inverseDistUp
			sumOfDistsSide = sumOfDistsSide + inverseDistSide

			for i in range(len(listOfPointsWithIds)):
				if(listOfPointsWithIds[i][0] == closestId):
					undeformedPoint = listOfPointsWithIds[i][1]
					deformedPoint = listOfDeformedPoints[i][1]
					numeratorUndeformedPoint = [numeratorUndeformedPoint[index] + undeformedPoint[index] * inverseDist for index in range(3)]
					numeratorDeformedPoint = [numeratorDeformedPoint[index] + deformedPoint[index] * inverseDist for index in range(3)]
				if (listOfPointsWithIds[i][0] == closestIdNext):
					undeformedPointNext = listOfPointsWithIds[i][1]
					deformedPointNext = listOfDeformedPoints[i][1]
					numeratorUndeformedPointNext = [numeratorUndeformedPointNext[index] + undeformedPointNext[index] * inverseDistNext for index in range(3)]
					numeratorDeformedPointNext = [numeratorDeformedPointNext[index] + deformedPointNext[index] * inverseDistNext for index in range(3)]
				if(listOfPointsWithIds[i][0] == closestIdUp):
					undeformedPointUp = listOfPointsWithIds[i][1]
					deformedPointUp = listOfDeformedPoints[i][1]
					numeratorUndeformedPointUp = [numeratorUndeformedPointUp[index] + undeformedPointUp[index] * inverseDistUp for index in range(3)]
					numeratorDeformedPointUp = [numeratorDeformedPointUp[index] + deformedPointUp[index] * inverseDistUp for index in range(3)]
				if(listOfPointsWithIds[i][0] == closestIdSide):
					undeformedPointSide = listOfPointsWithIds[i][1]
					deformedPointSide = listOfDeformedPoints[i][1]
					numeratorUndeformedPointSide = [numeratorUndeformedPointSide[index] + undeformedPointSide[index] * inverseDistSide for index in range(3)]
					numeratorDeformedPointSide = [numeratorDeformedPointSide[index] + deformedPointSide[index] * inverseDistSide for index in range(3)]

		undeformedPoint = [x / sumOfDists for x in numeratorUndeformedPoint]
		deformedPoint = [x / sumOfDists for x in numeratorDeformedPoint]
		undeformedPointNext = [x / sumOfDistsNext for x in numeratorUndeformedPointNext]
		deformedPointNext = [x / sumOfDistsNext for x in numeratorDeformedPointNext]
		undeformedPointUp = [x / sumOfDistsUp for x in numeratorUndeformedPointUp]
		deformedPointUp = [x / sumOfDistsUp for x in numeratorDeformedPointUp]
		undeformedPointSide = [x / sumOfDistsSide for x in numeratorUndeformedPointSide]
		deformedPointSide = [x / sumOfDistsSide for x in numeratorDeformedPointSide]

		listOfDeformedPointsX.append(deformedPoint[0])
		listOfDeformedPointsZ.append(deformedPoint[2])

		L0x = distanceBetweenPoints(undeformedPointNext, undeformedPoint)
		L0y = distanceBetweenPoints(undeformedPointSide, undeformedPoint)
		L0z = distanceBetweenPoints(undeformedPointUp, undeformedPoint)

		Lfx = distanceBetweenPoints(deformedPointNext, deformedPoint)
		Lfy = distanceBetweenPoints(deformedPointSide, deformedPoint)
		Lfz = distanceBetweenPoints(deformedPointUp, deformedPoint)

		strainInX = (Lfx - L0x) / L0x * 100
		strainInY = (Lfy - L0y) / L0y * 100
		strainInZ = (Lfz - L0z) / L0z * 100

		print "p",i,"\t",strainInX,"%",strainInY,"%",strainInZ,"%"

		listOfStrainsInX.append(strainInX)
		listOfStrainsInY.append(strainInY)
		listOfStrainsInZ.append(strainInZ)


	listOfAllStrains.append(listOfStrainsInX)
	listOfAllStrains.append(listOfStrainsInY)
	listOfAllStrains.append(listOfStrainsInZ)
	
	return listOfAllStrains

#for strains FEM:
def addingSubPlotFEMStrain(listOfPointsWithIds, listOfDeformedPoints):
	listOfDeformedPointsX = []
	listOfDeformedPointsZ = []
	listOfPoints = []
	for i in range(11):
		listOfPoints.append([i, 0.5, 0.5])

	listOfPointsUp = []
	for i in range(11):
		listOfPointsUp.append([i, 0.5, 0.9])

	listOfPointsSide = []
	for i in range(11):
		listOfPointsSide.append([i, 0.9, 0.5])

	listOfFoundIdsInOrder = []
	listOfFoundIdsInOrderUp = []
	listOfFoundIdsInOrderSide = []
	numOfPointsToInterpolate = 1
	exponentOfInterpolation = 0

	for i in range(11):
		closestpoint_id_sph = findIdWithMinDistance(listOfPoints[i], listOfPointsWithIds, numOfPointsToInterpolate)
		listOfFoundIdsInOrder.append(closestpoint_id_sph)

	for i in range(11):
		closestpoint_id_sph = findIdWithMinDistance(listOfPointsUp[i], listOfPointsWithIds, numOfPointsToInterpolate)
		listOfFoundIdsInOrderUp.append(closestpoint_id_sph)

	for i in range(11):
		closestpoint_id_sph = findIdWithMinDistance(listOfPointsSide[i], listOfPointsWithIds, numOfPointsToInterpolate)
		listOfFoundIdsInOrderSide.append(closestpoint_id_sph)

	print "STRAINS"
	listOfStrainsInX = []
	listOfStrainsInY = []
	listOfStrainsInZ = []
	listOfAllStrains = []

	for idsIndex in range(10):

		numeratorUndeformedPoint = [0,0,0]
		numeratorDeformedPoint = [0,0,0]
		numeratorUndeformedPointNext = [0, 0, 0]
		numeratorDeformedPointNext = [0, 0, 0]
		numeratorUndeformedPointUp = [0,0,0]
		numeratorDeformedPointUp = [0,0,0]
		numeratorUndeformedPointSide = [0,0,0]
		numeratorDeformedPointSide = [0,0,0]
		sumOfDists = 0
		sumOfDistsNext = 0
		sumOfDistsUp = 0
		sumOfDistsSide = 0
		for nPoints in range(numOfPointsToInterpolate):
			closestId = listOfFoundIdsInOrder[idsIndex][nPoints][2]
			closestIdNext = listOfFoundIdsInOrder[idsIndex + 1][nPoints][2]
			closestIdUp = listOfFoundIdsInOrderUp[idsIndex][nPoints][2]
			closestIdSide = listOfFoundIdsInOrderSide[idsIndex][nPoints][2]

			inverseDist = pow(1. / max(listOfFoundIdsInOrder[idsIndex][nPoints][0], 1e-3), exponentOfInterpolation)
			inverseDistNext = pow(1. / max(listOfFoundIdsInOrder[idsIndex + 1][nPoints][0], 1e-3), exponentOfInterpolation)
			inverseDistUp = pow(1. / max(listOfFoundIdsInOrderUp[idsIndex][nPoints][0], 1e-3), exponentOfInterpolation)
			inverseDistSide = pow(1. / max(listOfFoundIdsInOrderSide[idsIndex][nPoints][0], 1e-3), exponentOfInterpolation)

			sumOfDists = sumOfDists + inverseDist
			sumOfDistsNext = sumOfDistsNext + inverseDistNext
			sumOfDistsUp = sumOfDistsUp + inverseDistUp
			sumOfDistsSide = sumOfDistsSide + inverseDistSide

			for i in range(len(listOfPointsWithIds)):
				if(listOfPointsWithIds[i][0] == closestId):
					undeformedPoint = listOfPointsWithIds[i][1]
					deformedPoint = listOfDeformedPoints[i][1]
					numeratorUndeformedPoint = [numeratorUndeformedPoint[index] + undeformedPoint[index] * inverseDist for index in range(3)]
					numeratorDeformedPoint = [numeratorDeformedPoint[index] + deformedPoint[index] * inverseDist for index in range(3)]
				if (listOfPointsWithIds[i][0] == closestIdNext):
					undeformedPointNext = listOfPointsWithIds[i][1]
					deformedPointNext = listOfDeformedPoints[i][1]
					numeratorUndeformedPointNext = [numeratorUndeformedPointNext[index] + undeformedPointNext[index] * inverseDistNext for index in range(3)]
					numeratorDeformedPointNext = [numeratorDeformedPointNext[index] + deformedPointNext[index] * inverseDistNext for index in range(3)]
				if(listOfPointsWithIds[i][0] == closestIdUp):
					undeformedPointUp = listOfPointsWithIds[i][1]
					deformedPointUp = listOfDeformedPoints[i][1]
					numeratorUndeformedPointUp = [numeratorUndeformedPointUp[index] + undeformedPointUp[index] * inverseDistUp for index in range(3)]
					numeratorDeformedPointUp = [numeratorDeformedPointUp[index] + deformedPointUp[index] * inverseDistUp for index in range(3)]
				if(listOfPointsWithIds[i][0] == closestIdSide):
					undeformedPointSide = listOfPointsWithIds[i][1]
					deformedPointSide = listOfDeformedPoints[i][1]
					numeratorUndeformedPointSide = [numeratorUndeformedPointSide[index] + undeformedPointSide[index] * inverseDistSide for index in range(3)]
					numeratorDeformedPointSide = [numeratorDeformedPointSide[index] + deformedPointSide[index] * inverseDistSide for index in range(3)]

		undeformedPoint = [x / sumOfDists for x in numeratorUndeformedPoint]
		deformedPoint = [x / sumOfDists for x in numeratorDeformedPoint]
		undeformedPointNext = [x / sumOfDistsNext for x in numeratorUndeformedPointNext]
		deformedPointNext = [x / sumOfDistsNext for x in numeratorDeformedPointNext]
		undeformedPointUp = [x / sumOfDistsUp for x in numeratorUndeformedPointUp]
		deformedPointUp = [x / sumOfDistsUp for x in numeratorDeformedPointUp]
		undeformedPointSide = [x / sumOfDistsSide for x in numeratorUndeformedPointSide]
		deformedPointSide = [x / sumOfDistsSide for x in numeratorDeformedPointSide]

		listOfDeformedPointsX.append(deformedPoint[0])
		listOfDeformedPointsZ.append(deformedPoint[2])

		L0x = distanceBetweenPoints(undeformedPointNext, undeformedPoint)
		L0y = distanceBetweenPoints(undeformedPointSide, undeformedPoint)
		L0z = distanceBetweenPoints(undeformedPointUp, undeformedPoint)

		Lfx = distanceBetweenPoints(deformedPointNext, deformedPoint)
		Lfy = distanceBetweenPoints(deformedPointSide, deformedPoint)
		Lfz = distanceBetweenPoints(deformedPointUp, deformedPoint)

		strainInX = (Lfx - L0x) / L0x * 100
		strainInY = (Lfy - L0y) / L0y * 100
		strainInZ = (Lfz - L0z) / L0z * 100

		print "p",i,"\t",strainInX,"%",strainInY,"%",strainInZ,"%"

		listOfStrainsInX.append(strainInX)
		listOfStrainsInY.append(strainInY)
		listOfStrainsInZ.append(strainInZ)


	listOfAllStrains.append(listOfStrainsInX)
	listOfAllStrains.append(listOfStrainsInY)
	listOfAllStrains.append(listOfStrainsInZ)
	
	return listOfAllStrains


def addingSubPlots(undeformedFileNameFolder):
	#we need to create the deformedFileName from the undeformedFileNameFolder
#	deformedFileName = undeformedFileNameFolder.join('deformed.vtk')
	undeformed = "/undeformed.vtk"
	undeformedFileName = undeformedFileNameFolder + undeformed
	deformed = "/deformed.vtk"
	deformedFileName = undeformedFileNameFolder + deformed

	reader = vtkUnstructuredGridReader()
	reader.SetFileName(undeformedFileName)
	reader.Update()
	undeformedMesh = reader.GetOutput()
	listOfUndeformedPointsX = []
	listOfUndeformedPointsZ = []

#	readerSmall = vtkPolyDataReader();
	#readerSmall.SetFileName('/home/eric/Downloads/sander314-mechbench-f541b4f0e0b3/problem 1/CardioMechanics/139587/undeformedPolyData1.vtk');
	#readerSmall.SetFileName('/home/eric/Downloads/sander314-mechbench-f541b4f0e0b3/problem 1/Simula-FEniCS/8444/undeformedPolyData.vtk');
#	readerSmall.SetFileName('/home/eric/Downloads/sander314-mechbench-f541b4f0e0b3/problem 1/glasgowHeart-IBFE/40x4x4(structure)_120x42x90(fluid)/underofmedPolyData.vtk');
#	readerSmall.Update();
#	m = readerSmall.GetOutput()

#	p = m.GetPoint(1)
#	print(m.GetNumberOfPoints())

#	cl = vmtkcenterlines(m, [0, 0.5, 0.5], [10, 0.5, 0.5])
#	print(cl.GetNumberOfPoints())

#	writer1 = vtkPolyDataWriter()
#	writer1.SetFileName('outputUnstructured.vtk')
#	writer1.SetInput(cl)
#	writer1.Write()

	listOfPoints = []
	for i in range(11):
		listOfPoints.append([i, 0.5, 0.5])

	listOfPointsUp = []
	for i in range(11):
		listOfPointsUp.append([i, 0.5, 0.9])

	listOfPointsSide = []
	for i in range(11):
		listOfPointsSide.append([i, 0.9, 0.5])



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

	readerDeformed = vtkUnstructuredGridReader()
	readerDeformed.SetFileName(deformedFileName)
	readerDeformed.Update()
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

		listOfUndeformedPointsX.append(deformedPoint[0])
		listOfUndeformedPointsZ.append(deformedPoint[2])

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
		print deformedPoint[0], deformedPoint[2]
	plt.plot(listOfUndeformedPointsX, listOfUndeformedPointsZ, 'ro', ms = 10)
	plt.plot(listOfUndeformedPointsX, listOfUndeformedPointsZ)
	# You can specify a rotation for the tick labels in degrees or with keywords.
	# Pad margins so that markers don't get clipped by the axes
	plt.margins(0.2)
	# Tweak spacing to prevent clipping of tick-labels
	plt.subplots_adjust(bottom=0.15)

def centerLineSPH(listOfPointsWithIds, listOfDeformedPoints):
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
		listOfPointsSecondPlot.append([i * 10./resolutionOfSecondPlot, 0.5, 0.5])
		listOfPointsZ.append(0.5);

	for i in range(resolutionOfSecondPlot):		
		closestpoints_id_sph = findIdWithMinDistance(listOfPointsSecondPlot[i], listOfPointsWithIds, 3)
		listOfFoundIdsInOrderSecondPlot.append(closestpoints_id_sph)
	powerOfInverseDistances = 3
	for idsIndex in range(resolutionOfSecondPlot):
#		print listOfDeformedPoints[listOfFoundIdsInOrderSecondPlot[i]], listOfDeformedPoints[listOfFoundIdsInOrderSecondPlot[i]][1][2]
		print "\n"
		closestId1 = listOfFoundIdsInOrderSecondPlot[idsIndex][0][2]
		closestId2 = listOfFoundIdsInOrderSecondPlot[idsIndex][1][2]
		closestId3 = listOfFoundIdsInOrderSecondPlot[idsIndex][2][2]

		dist1 = pow(1. / max(listOfFoundIdsInOrderSecondPlot[idsIndex][0][0],0.001), powerOfInverseDistances)
		dist2 = pow(1. / max(listOfFoundIdsInOrderSecondPlot[idsIndex][0][0],0.001), powerOfInverseDistances)
		dist3 = pow(1. / max(listOfFoundIdsInOrderSecondPlot[idsIndex][0][0],0.001), powerOfInverseDistances)
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

	for idsIndex in range(resolutionOfSecondPlot):
#		print listOfDeformedPoints[listOfFoundIdsInOrderSecondPlot[i]], listOfDeformedPoints[listOfFoundIdsInOrderSecondPlot[i]][1][2]
		closestId1 = listOfFoundIdsInOrderSecondPlot[idsIndex][0][2]
		closestId2 = listOfFoundIdsInOrderSecondPlot[idsIndex][1][2]
		closestId3 = listOfFoundIdsInOrderSecondPlot[idsIndex][2][2]

		dist1 = pow(1. / max(listOfFoundIdsInOrderSecondPlot[idsIndex][0][0],0.001), powerOfInverseDistances)
		dist2 = pow(1. / max(listOfFoundIdsInOrderSecondPlot[idsIndex][0][0],0.001), powerOfInverseDistances)
		dist3 = pow(1. / max(listOfFoundIdsInOrderSecondPlot[idsIndex][0][0],0.001), powerOfInverseDistances)
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
	
	l2, = plt.plot(xSecondPlot, listOfPointsZ, linestyle = '-', color = 'gray', linewidth = 3)
#	plt.plot(xSecondPlot, listOfPointsZ, 'o')

#	plt.show()
	l3, = plt.plot(xSecondPlot, deformedPointsZHighResolution, '-bo', linewidth = 3, ms = 10)
#	plt.plot(xSecondPlot, deformedPointsZHighResolution, 'ro', ms = 10)
	return l3

def centerLineFEM(listOfPointsWithIds, listOfDeformedPoints):
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
		listOfPointsSecondPlot.append([i * 10./resolutionOfSecondPlot, 0.5, 0.5])
#		xSecondPlot.append(i * 10./resolutionOfSecondPlot);


	for i in range(resolutionOfSecondPlot):		
		closestpoints_id_sph = findIdWithMinDistance(listOfPointsSecondPlot[i], listOfPointsWithIds, 3)
		listOfFoundIdsInOrderSecondPlot.append(closestpoints_id_sph)
	powerOfInverseDistances = 2
	for idsIndex in range(resolutionOfSecondPlot):
#		print listOfDeformedPoints[listOfFoundIdsInOrderSecondPlot[i]], listOfDeformedPoints[listOfFoundIdsInOrderSecondPlot[i]][1][2]
		print "\n"
		closestId1 = listOfFoundIdsInOrderSecondPlot[idsIndex][0][2]
		closestId2 = listOfFoundIdsInOrderSecondPlot[idsIndex][1][2]
		closestId3 = listOfFoundIdsInOrderSecondPlot[idsIndex][2][2]

		dist1 = pow(1. / max(listOfFoundIdsInOrderSecondPlot[idsIndex][0][0],0.001), powerOfInverseDistances)
		dist2 = pow(1. / max(listOfFoundIdsInOrderSecondPlot[idsIndex][0][0],0.001), powerOfInverseDistances)
		dist3 = pow(1. / max(listOfFoundIdsInOrderSecondPlot[idsIndex][0][0],0.001), powerOfInverseDistances)
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

	for idsIndex in range(resolutionOfSecondPlot):
#		print listOfDeformedPoints[listOfFoundIdsInOrderSecondPlot[i]], listOfDeformedPoints[listOfFoundIdsInOrderSecondPlot[i]][1][2]
		closestId1 = listOfFoundIdsInOrderSecondPlot[idsIndex][0][2]
		closestId2 = listOfFoundIdsInOrderSecondPlot[idsIndex][1][2]
		closestId3 = listOfFoundIdsInOrderSecondPlot[idsIndex][2][2]

		dist1 = pow(1. / max(listOfFoundIdsInOrderSecondPlot[idsIndex][0][0],0.001), powerOfInverseDistances)
		dist2 = pow(1. / max(listOfFoundIdsInOrderSecondPlot[idsIndex][0][0],0.001), powerOfInverseDistances)
		dist3 = pow(1. / max(listOfFoundIdsInOrderSecondPlot[idsIndex][0][0],0.001), powerOfInverseDistances)
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
	
	l2, = plt.plot(listOfPointsX, listOfPointsZ, linestyle = '-', color = 'gray', linewidth = 3)
#	plt.plot(listOfPointsX, listOfPointsZ, 'o')

#	plt.show()
	l3, = plt.plot(xSecondPlot, deformedPointsZHighResolution, '--y^', linewidth = 3, ms = 10)
#	plt.plot(xSecondPlot, deformedPointsZHighResolution, 'y^', ms = 10)
#	plt.show()
	return l3



plot_lines = []					
for root, directories, filenames in os.walk('/home/eric/scriptsToHelpSPHCreateGeometryInteraceETC/centerline/meshesWithMoreDegreesOfFreedomProb2'):
	for filename in filenames:
		if(filename == 'undeformed.vtk'):
			undeformedFileNameFolder = os.path.join(root)
			print undeformedFileNameFolder
			
			# calls the function centerLineFEM and centerLineSPH to compute both centerlines. The output is an iterative plot
			# stored in l1 and l3 to obtain the final figure for the journal
			l1 = centerLineFEM(pointsWithIds_FEM, deformedPointsWithIds_FEM)
			l3 = centerLineSPH(pointsWithIds_sph, deformedPointsWithIds_sph)

			plot_lines.append([l1,l3])
			legend1=plt.legend(plot_lines[0],["FEM","SPH"],loc=4)
			plt.gca().add_artist(legend1)
			plt.ylabel('z (mm)')
			plt.xlabel('x (mm)')
			plt.savefig("centerlineComparison", bbox_inches='tight', dpi = 300)

			plt.clf()
plt.rc('xtick', labelsize=20) 
plt.rc('ytick', labelsize=20) 

# it calls functions addingSubPlotSPHStrain and addingSubPlotFEMStrain to obtain the 
# strains for the bar for SPH and FEM
listOfStrainsSPH = addingSubPlotSPHStrain(pointsWithIds_sph, deformedPointsWithIds_sph)
listOfStrainsFEM = addingSubPlotFEMStrain(pointsWithIds_FEM, deformedPointsWithIds_FEM)
listOfLabels = ['x-axis','y-axis','z-axis']
for graphicToPlot in range(3):
	plot_lines = []
	x = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
	labelsPoints = ["p0","p1","p2","p3","p4","p5","p6","p7","p8","p9"]
	l1, = plt.plot(x, listOfStrainsFEM[graphicToPlot], '--y^', ms = 10, linewidth = 3)
#	l1, = plt.plot(x, listOfStrainsFEM[graphicToPlot], '--', linewidth=3)
	plt.xticks(x, labelsPoints, rotation='vertical')
	plt.margins(0.2)
	plt.subplots_adjust(bottom=0.15)

	l2, = plt.plot(x, listOfStrainsSPH[graphicToPlot], '-bo', ms = 10, linewidth = 3)
#	l2, = plt.plot(x, listOfStrainsSPH[graphicToPlot], linewidth=3)
	plt.xticks(x, labelsPoints, rotation='vertical')
	plt.margins(0.2)
	plt.subplots_adjust(bottom=0.15)

	plot_lines.append([l1,l2])
	if(graphicToPlot < 1):
		legend1=plt.legend(plot_lines[0],["FEM","SPH"],loc=4)
	else:
		legend1=plt.legend(plot_lines[0],["FEM","SPH"],loc=1)

	plt.gca().add_artist(legend1)

	plt.savefig(listOfLabels[graphicToPlot], bbox_inches='tight', dpi = 300)
	plt.clf()

