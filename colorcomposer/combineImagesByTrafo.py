import numpy as np
from numpy.linalg import solve
from scipy.misc import imsave
import scipy.odr as odr
from matplotlib import pyplot
from munkres import Munkres
import copy

def getRandomElements(list, numberElements):
	
	if len(list) < numberElements:
		print 'list does not contain enough elements'
		return 0
	randomElements = []
	indices = []
	
	while 1:
		index = np.random.randint(len(list))
		addItem = True
		for element in randomElements:
			if (list[index] == element).all():
				addItem = False
				break
		if addItem:
			randomElements.append(list[index])
			indices.append(index)
			if len(randomElements)==numberElements:
				break
		
	return np.array(randomElements), indices

def performRansac(s,d,dims,pointsNeeded, numberIterations, distanceTolerated):
	indicesSCollection = []
	indicesDCollection = []
	additionalCandidatesS = []
	additionalCandidatesD = []
	trafoCollection = []
	matchedPoints = []
	
	for i in range(numberIterations):
		sSet, indicesS = getRandomElements(s,pointsNeeded)
		dSet, indicesD = getRandomElements(d,pointsNeeded)
		sSet, dSet = alignCandidates(sSet, dSet)
		indicesSCollection.append(indicesS)
		indicesDCollection.append(indicesD)
		trafo = affineMatrix2DFromCorrespondingPoints(sSet, dSet, dims)
		trafoCollection.append(trafo)
		sTransformed = np.dot(trafo, np.vstack([s.T,np.ones(len(s))])).T[:,:2]
		distanceMatrix = getDistanceMatrix(sTransformed, d)
		minimumRow = np.min(distanceMatrix, 1)		
		matchedPoints.append(np.sum(np.where(minimumRow < distanceTolerated,1,0)))
		row, col =np.where(distanceMatrix < distanceTolerated)
		additionalCandidatesD.append(row)
		additionalCandidatesS.append(col)
		if i%500 == 0:
			print '(%i/%i)'%(i,numberIterations)
		
	return matchedPoints, additionalCandidatesS, additionalCandidatesD, indicesSCollection, indicesDCollection, trafoCollection

def preselectPoints(s,d,maxDist):
	snew = []
	dnew = []
	for i in range(len(s)):
		inRange = False
		for j in range(len(d)):
			if ((s[i][0]-d[j][0])**2 + (s[i][1]-d[j][1])**2) < maxDist**2:
				inRange = True
		if inRange:
			snew.append(s[i])		

	for j in range(len(d)):
		inRange = False
		for i in range(len(snew)):
			if ((snew[i][0]-d[j][0])**2 + (snew[i][1]-d[j][1])**2) < maxDist**2:
				inRange = True
		if inRange:
			dnew.append(d[j])
		 
	return np.array(snew), np.array(dnew)

def doRansac(s,d,dims):
	pointsNeeded = 3
	numberIterations = 2000
	numberRestarts = 50
	distanceTolerated = (dims[0]+dims[1])/400  #This sets the distance for the check how much more beads are overlapping 
	toleranceForShearing = 10.04				#set this variable to a small value, to suppress shearing
	maxDist = (dims[0]+dims[1])/2./1			#Beads that have no nearer neighbor of the other color than maxDist are not considered
	print len(s), len(d)
	s,d = preselectPoints(s,d,maxDist)
	print len(s), len(d)
	
	sc=copy.deepcopy(s)		# I don't know why but affineMatrix2DFromCorrespondingPoints changes the values given (looks like a call by reference)
	dc=copy.deepcopy(d)
	
	if len(s) < pointsNeeded or len(d)< pointsNeeded:
		print 'To few points found'
		return
	print 'begin Ransac'
	for i in range(numberRestarts):
		matchedPoints, additionalCandidatesS, additionalCandidatesD, indicesSCollection, indicesDCollection, trafo = performRansac(sc,dc,dims,pointsNeeded, numberIterations, distanceTolerated)
		maximalMatchedPoints = np.max(matchedPoints)
		print maximalMatchedPoints
		if maximalMatchedPoints > np.min(matchedPoints):	
			indexWithMaximalAlignment = np.where(matchedPoints==maximalMatchedPoints)[0][0]
			#breaks if the trafo has no shear
			print trafo[indexWithMaximalAlignment]
			if (1 + toleranceForShearing>(trafo[indexWithMaximalAlignment][0,0]**2+trafo[indexWithMaximalAlignment][0,1]**2) and (trafo[indexWithMaximalAlignment][0,0]**2+trafo[indexWithMaximalAlignment][0,1]**2)>.1 - toleranceForShearing and 1 + toleranceForShearing>(trafo[indexWithMaximalAlignment][1,0]**2+trafo[indexWithMaximalAlignment][1,1]**2) and (trafo[indexWithMaximalAlignment][1,0]**2+trafo[indexWithMaximalAlignment][1,1]**2)>0.1 - toleranceForShearing):
				break
			print 'not enough loops in Ransac min=max . Start Iteration again (%i/%i)' %(i, numberRestarts)
		else:
			doBreak = False
			for j in range(len(trafo)):				
				if (1 + toleranceForShearing>(trafo[j][0,0]**2+trafo[j][0,1]**2) and (trafo[j][0,0]**2+trafo[j][0,1]**2)>.1 - toleranceForShearing and 1 + toleranceForShearing>(trafo[j][1,0]**2+trafo[j][1,1]**2) and (trafo[j][1,0]**2+trafo[j][1,1]**2)>0.1 - toleranceForShearing):
					indexWithMaximalAlignment = j
					doBreak = True
					break
			if doBreak:
				break
			print 'not enough loops in Ransac shear. Start Iteration again (%i/%i)' %(i, numberRestarts)
				#print 'For the next try correlation between near beads will be considered'
				#matchedPoints, additionalCandidatesS, additionalCandidatesD, indicesSCollection, indicesDCollection = tryHungarianApproach(sc,dc, dims)
	print maximalMatchedPoints
	finalIndicesS = []
	finalIndicesD = []
	'''for i in range(len(matchedPoints)):
		if matchedPoints[i] == maximalMatchedPoints:
			#finalIndicesS.append(indicesSCollection[i])
			#finalIndicesD.append(indicesDCollection[i])
			finalIndicesS.append(additionalCandidatesS[i])
			finalIndicesD.append(additionalCandidatesD[i])'''
	
	indexWithMaximalAlignment = np.where(matchedPoints==maximalMatchedPoints)[0][0]
	
	finalIndicesS = additionalCandidatesS[indexWithMaximalAlignment]
	finalIndicesD = additionalCandidatesD[indexWithMaximalAlignment]
			
	finalIndicesS = (np.reshape(np.array(finalIndicesS),-1))
	finalIndicesD = (np.reshape(np.array(finalIndicesD),-1))
	
	finalTrafo = affineMatrix2DFromCorrespondingPoints(s[finalIndicesS], d[finalIndicesD], dims)
	
	print 'Found %i, Points for transformation' %len(finalIndicesS)
	print finalTrafo
	print finalTrafo[0,0]**2+finalTrafo[0,1]**2,finalTrafo[1,0]**2+finalTrafo[1,1]**2
	return finalTrafo,s[finalIndicesS],d[finalIndicesD]

def getDistanceMatrix(s,d):
	numberBeadsRed = d.shape[0]
	numberBeadsGreen = s.shape[0]
	matrix = np.zeros((numberBeadsRed, numberBeadsGreen))		#matrix is the Distance matrix between all red and green points
	for i in range(numberBeadsRed):
		for j in range(numberBeadsGreen):
			matrix[i,j] = np.sqrt((d[i,0]-s[j,0])**2+(d[i,1]-s[j,1])**2)
			
	return matrix

def alignCandidates(s,d):
	matrix = getDistanceMatrix(s,d)
	numberBeadsRed = d.shape[0]
	numberBeadsGreen = s.shape[0]
	
	m = Munkres()
	indexes = m.compute(matrix)
	
	alignedd = []
	aligneds = []
	
	total = 0
	for row, column in indexes:
	    value = matrix[row][column]
	    total += value

	    alignedd.append(d[row])
	    aligneds.append(s[column])
	    indicesDCollection = row
	    indicesSCollection = column
	   
	#trafo = affineMatrix2DFromCorrespondingPoints(s[indicesSCollection], d[indicesDCollection], dims)
	return np.array(aligneds), np.array(alignedd)

def tryHungarianApproach(s,d,dims):
	indexS, indexD = doRansac(s,d,dims)
	getRandomElements(s,3)
	matrix = getDistanceMatrix(s,d)

	print matrix
	##############################################################################
	#There are problems if there are for example 2 green but just one red point near each other
	#This is a very rough ansatz
	minimumRow = np.min(matrix, 1)
	meanMinRow = np.mean(minimumRow)	#gets the mean of the smalles values of each row, this is used to detect points far away from any partner
	
	minimumCol = np.min(matrix, 0)		
	meanMinCol = np.mean(minimumCol)	#same but for the green points
	
	row,col = np.where(matrix<2*meanMinRow)	
	d = d[np.unique(row)]					#red points that are far away from any green point are deleted
	row,col = np.where(matrix<2*meanMinCol)
	s = s[np.unique(col)]					#green points that are far away from any green point are deleted
	
	######################################################
	#Distance matrix is calculated again, but without the outliers
	matrix = getDistanceMatrix(s,d)
	numberBeadsRed = d.shape[0]
	numberBeadsGreen = s.shape[0]
	
	if numberBeadsRed>numberBeadsGreen:
		matrix = np.transpose(matrix)
		transposed = True						#munkres can't handle mxn matrices with n>m
	else:
		transposed = False
	
	m = Munkres()
	indexes = m.compute(matrix)
	
	alignedd = []
	aligneds = []
	
	total = 0
	for row, column in indexes:
	    value = matrix[row][column]
	    total += value
	    #print '(%d, %d) -> %d' % (row, column, value)
	    #print 'total cost: %d' % total
	    if transposed == False:
		    alignedd.append(d[row])
		    aligneds.append(s[column])
		    indicesDCollection = row
		    indicesSCollection = column
	    if transposed == True:
	    	alignedd.append(d[column])
	    	aligneds.append(s[row])
	    	indicesDCollection = column
	    	indicesSCollection = row
	    	
	trafo = affineMatrix2DFromCorrespondingPoints(s[indicesSCollection], d[indicesDCollection], dims)
	sTransformed = np.dot(trafo, np.vstack([s.T,np.ones(len(s))])).T[:,:2]
	distanceMatrix = getDistanceMatrix(sTransformed, d)
	minimumRow = np.min(distanceMatrix, 1)		
	matchedPoints=[0,(np.sum(np.where(minimumRow < distanceTolerated,1,0)))]
	
	additionalCandidatesD = indicesDCollection
	additionalCandidatesS = indicesSCollection
	matchedPoints = [0,10]				#Just fake numbers, but because min != max, the programm will go on
    	   
	#aligneds=np.array(aligneds)
	#alignedd=np.array(alignedd)
	
	#return aligneds, alignedd
	return matchedPoints, additionalCandidatesS, additionalCandidatesD, indicesSCollection, indicesDCollection

def affineMatrix2DFromCorrespondingPoints(s, d, dims):
	n = len(s)
	if len(s) < 3:
		print "at least three points required"
		return
	
	deltaX=dims[0]/2.
	deltaY=dims[1]/2.
	#deltaX=np.mean(d[:,0])
	#deltaY=np.mean(d[:,1])
	#deltaX=0
	#deltaY=0
	
	s[:,0]=s[:,0]-deltaX
	s[:,1]=s[:,1]-deltaY
	d[:,0]=d[:,0]-deltaX
	d[:,1]=d[:,1]-deltaY
	
	m  = np.zeros((3,3))
	
	s = np.hstack([s,np.ones((n,1))]) # add third column with ones
	d = np.hstack([d,np.ones((n,1))]) # add third column with ones
	rx = np.sum(d[:,0:1]*s,axis=0)
	ry = np.sum(d[:,1:2]*s,axis=0)

	for i in range(n):
		c = s[i]
		m += np.outer(c,c)

	solx = solve(m, rx)
	soly = solve(m, ry)
	row3 = np.array([0,0,1])
	
	A = np.vstack([solx,soly,row3])
	R=np.matrix([[1, 0, deltaX],[0, 1, deltaY],[0, 0, 1]])
	Ap = np.dot(R,np.dot(A,R.I))
	
	#print "matrix:",[solx,soly,row3]
	a=np.vstack([solx,soly])
	
	error=(d[:,0:2]).T-np.dot(a,s.T)
	#print error,"=",(d[:,0:2]).T,"-",np.dot(a,s.T)
	#print "error",error
	
	#dA=-np.dot(error,(np.dot(s,np.linalg.pinv(np.dot(s.T,s)))))
	#print np.linalg.pinv(s).shape, error.shape
	dA=-np.dot(error,np.linalg.pinv(s).T)
	
	#print dA
	#print np.dot([solx,soly,row3],s.T)
	
	#return np.vstack([solx,soly,row3])
	return np.array([[Ap[0,0],Ap[0,1],Ap[0,2]],[Ap[1,0],Ap[1,1],Ap[1,2]],[Ap[2,0],Ap[2,1],Ap[2,2]]])


def linFun(B,x):
  return (B[0]*x+B[1])
 
def calculateRegression(s,d):
	
	#pyplot.scatter(s[:,0],d[:,0])
	#pyplot.show()
	#pyplot.scatter(s[:,1],d[:,1])
	#pyplot.show()
	xbead=d[:,0]
	xsbead=s[:,0]
	
	ybead=d[:,1]
	ysbead=s[:,1]
	
	linearModel = odr.Model(linFun)
	myDataX = odr.Data(xbead,xsbead)
	myodrX = odr.ODR(myDataX,linearModel,beta0=[1,1])
	myoutputX = myodrX.run()
	
	myDataY = odr.Data(ybead,ysbead)
	myodrY = odr.ODR(myDataY,linearModel,beta0=[1,1])
	myoutputY = myodrY.run()
	#myoutput.pprint()
	
	parametersX = myoutputX.res_var
	parametersY = myoutputY.res_var
	
	return [parametersX, parametersY]
	
# if run standalone
if __name__ == "__main__":
	import matplotlib.pyplot as plt
	import coords
	file1 = 'data/HeLa5_10000z_3_polyL_actinmEOS2-steve.txt'
	file1_landmarks = 'data/HeLa5_10000z_3_polyL_actinmEOS2-steve.txt_beads.txt'
	file2 = 'data/HeLa5_10000z_3_polyL_ER647.txt'
	file2_landmarks = 'data/HeLa5_10000z_3_polyL_ER647.txt_beads.txt'

	thrash, landmarks1 = coords.readfile(file1_landmarks)
	thrash, landmarks2 = coords.readfile(file2_landmarks)
	landmarks1 = landmarks1[:,:2]
	landmarks2 = landmarks2[:,:2]
	trafo = affineMatrix2DFromCorrespondingPoints(landmarks2, landmarks1)
	print "Transformation: \n", trafo

	dims, cc1 = coords.readfile(file1)
	dims, cc2 = coords.readfile(file2)
	cc2_ones = np.hstack([cc2[:,:2], np.ones((len(cc2),1))]) # add column with ones for affine trafo
	cc2_transformed = np.dot(cc2_ones, trafo.T)[:,:2]
	cc2_transformed = np.hstack([cc2_transformed,cc2[:,2:]]) # add intensity information again

	im1 = coords.coords2Image(dims, cc1)
	im2 = coords.coords2Image(dims, cc2_transformed)

	im_color = np.dstack([im1/np.max(im1),im2/np.max(im2),np.zeros(im1.shape)])
	scipy.misc.imsave(file1+"_processed.png", im_color)
	plt.imshow(im_color)
	plt.plot(landmarks1[:,0]*8, landmarks1[:,1]*8, 'wo', alpha=0.4, scalex=False, scaley=False)
	plt.show()
