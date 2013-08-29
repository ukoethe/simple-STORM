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

    return np.array(randomElements), np.array(indices)

def getGoodPartnerIndex(indicesS, distMat, d):
    # This function tries to find the best red beads for the randomly selected green beads, nearer beads get higher possibilities to be chosen
    indicesD = []
    for i in range(len(indicesS)):
        actuallColumn = distMat[:,indicesS[i]]
        maxActuallColumn = np.max(actuallColumn)
        actuallColumn = np.where(actuallColumn == 0, 0.01, actuallColumn)
        sumCol = np.sum(actuallColumn)
        actuallColumnProb = sumCol / actuallColumn**1.8 #at this point the sensitivity for picking the nearest point can be influenced
        randomValue = np.random.random() * np.sum(actuallColumnProb)
        for j in range(len(actuallColumnProb)):
            if (np.sum(actuallColumnProb[0:j]) < randomValue < (np.sum(actuallColumnProb[0:j+1]))):
                if j in indicesD:
                    while 1: #this loop must break because len(d) >= len(indicesS)
                        nearestMatch = np.where(actuallColumn == np.min(actuallColumn))[0][0]
                        if not nearestMatch in indicesD:
                            indicesD.append(nearestMatch)
                            break
                        else:
                            actuallColumn[nearestMatch] = maxActuallColumn+1
                else:
                    indicesD.append(j)
                break
    indicesD = np.array(indicesD, dtype = int)
    return d[indicesD], indicesD

def performRansac(s,d,dims,pointsNeeded, numberIterations, distanceTolerated):
    indicesSCollection = []
    indicesDCollection = []
    additionalCandidatesS = []
    additionalCandidatesD = []
    trafoCollection = []
    matchedPoints = []
    distanceMatrixSD = getDistanceMatrix(s,d)

    for i in range(numberIterations):
        sSet, indicesS = getRandomElements(s,pointsNeeded)
        dSet, indicesD = getGoodPartnerIndex(indicesS, distanceMatrixSD, d)
        #dSet, indicesD = getRandomElements(d,pointsNeeded)
        #sSet, dSet = alignCandidates(sSet, dSet)
        indicesSCollection.append(indicesS)
        indicesDCollection.append(indicesD)
        trafo = affineMatrix2DFromCorrespondingPoints(sSet, dSet, dims)
        #trafo = affineMatrix2DFromCorrespondingPointsUseTotalLeastSquares(sSet, dSet, dims)
        trafoCollection.append(trafo)
        sTransformed = np.dot(trafo, np.vstack([s.T,np.ones(len(s))])).T[:,:2]
        distanceMatrix = getDistanceMatrix(sTransformed, d)
        minimumRow = np.min(distanceMatrix, 1)
        matchedPoints.append(np.sum(np.where(minimumRow < distanceTolerated,1,0)))  #all red bead that are closer than distanceTolerated to a green bead are added
        row, col =np.where(distanceMatrix < distanceTolerated)
        additionalCandidatesD.append(row)
        additionalCandidatesS.append(col)
        #if i%500 == 0:
        #	print '(%i/%i)'%(i,numberIterations)

    return matchedPoints, additionalCandidatesS, additionalCandidatesD, indicesSCollection, indicesDCollection, trafoCollection

def preselectPoints(s,d,maxDist):
    # If points are far away from any other point of the opposite color, they are not considered
    snew = []
    dnew = []
    for i in range(len(s)):
        for j in range(len(d)):
            if ((s[i][0]-d[j][0])**2 + (s[i][1]-d[j][1])**2) < maxDist**2:
                snew.append(s[i])
                break

    for j in range(len(d)):
        inRange = False
        for i in range(len(snew)):
            if ((snew[i][0]-d[j][0])**2 + (snew[i][1]-d[j][1])**2) < maxDist**2:
                dnew.append(d[j])
                break

    return np.array(snew), np.array(dnew)

def doRansac2(s,d,dims):
    #pointsNeeded = 5
    pointsNeededAtLeast = 3
    numberIterations = 500
    numberRestarts = 100
    distanceTolerated = (dims[0]+dims[1])/500.  #This sets the distance for the check of how much more beads are overlapping
    toleranceForShearing = 0.1                #set this variable to a small value, to suppress shearing
    maxDist = (dims[0]+dims[1])/2./7                #Beads that have no nearer neighbor of the other color than maxDist are not considered
    print len(s), len(d)
    s,d = preselectPoints(s,d,maxDist) #skipps all points that are far off any other point of the other color
    print len(s), len(d)
    pointsNeeded = min([len(s),len(d)])
    sc=copy.deepcopy(s)        #affineMatrix2DFromCorrespondingPoints changes the values given
    dc=copy.deepcopy(d)

    if min([len(s),len(d)]) < pointsNeededAtLeast:
        print 'To few points found, there might have been not enough pairs of matching points'
        return
    
    indicesToCheck = []
    listTrafos = []
    listAdditionalCandS = []
    listAdditionalCandD = []
    print 'begin Ransac'
    for i in range(numberRestarts):
        matchedPoints, additionalCandidatesS, additionalCandidatesD, indicesSCollection, indicesDCollection, trafo = performRansac(sc,dc,dims,pointsNeeded, numberIterations, distanceTolerated)
        maximalMatchedPoints = np.max(matchedPoints)
        print maximalMatchedPoints
        
        if maximalMatchedPoints<pointsNeeded:
            pointsNeeded -=1
            continue
        if maximalMatchedPoints > np.min(matchedPoints) or maximalMatchedPoints ==np.min([len(s),len(d)]):    # more combined points are found than points that are necessary for the calculation of the transformation
            indicesToCheck.append(np.where(np.asarray(matchedPoints)>=pointsNeededAtLeast))
            for i in indicesToCheck[0][0]:
                if len(additionalCandidatesS[i])>0:
                    listTrafos.append(trafo[i])
                    listAdditionalCandS.append(additionalCandidatesS[i])
                    listAdditionalCandD.append(additionalCandidatesD[i])
            #breaks if the trafo has no shear
            
            #if there is just translation and rotation trafo[0,0] should be cos(phi), trafo[0,1] - sin(phi)..., therefore on can check if there is shearing
            
            
                    
            if pointsNeeded >pointsNeededAtLeast:
                pointsNeeded -= 1
                print 'Start Iteration again with decreased number of points for transformation'
            else:
                break
        else:
            for j in range(len(trafo)): # just the points used to calculate the transformation lay above each other, either there are just 3 pairs of beads, or the correct transformation has not been found jet    
                if len(additionalCandidatesS[j])>0:
                    listTrafos.append(trafo[j])
                    listAdditionalCandS.append(additionalCandidatesS[j])
                    listAdditionalCandD.append(additionalCandidatesD[j])
#             if exitRans:
#                 break
            print 'Start Iteration again (%i/%i), reason: too few points found' %(i, numberRestarts)
    

    rmsList = []
    for currTrafo, currAdditionalS, currAdditionalD in zip(listTrafos,listAdditionalCandS,listAdditionalCandD):
        if (np.abs((currTrafo[0,0]**2+currTrafo[0,1]**2)-1) < toleranceForShearing and np.abs((currTrafo[1,0]**2+currTrafo[1,1]**2)-1) <  toleranceForShearing):
            tt = np.dot(currTrafo, np.vstack([s[currAdditionalS].T,np.ones(len(currAdditionalS))])).T[:,:2] # fg transformed in bg coords
            rss = np.sum((tt-d[currAdditionalD])**2)
            rms = np.sqrt(rss/len(currAdditionalD))
            rmsList.append(rms)
            print rms, rss
        else:
            rmsList.append(999999)
#     for i in (bestIndices):
#         trafo = affineMatrix2DFromCorrespondingPoints(s[additionalCandidatesS[i]], d[additionalCandidatesD[i]], dims)
#         trafoList.append(trafo)
#         tt = np.dot(trafo, np.vstack([s[additionalCandidatesS[i]].T,np.ones(len(additionalCandidatesS[i]))])).T[:,:2] # fg transformed in bg coords
#         rss = np.sum((tt-d[additionalCandidatesD[i]])**2)
#         rms = np.sqrt(rss/len(additionalCandidatesD[i]))
#         rmsList.append(rms)
#         print rms, rss
    bestIdx = np.where(rmsList == np.min(rmsList))[0][0]
    print bestIdx
    return listTrafos[bestIdx],s[listAdditionalCandS[bestIdx]],d[listAdditionalCandD[bestIdx]]

def doRansac(s,d,dims):
    #pointsNeeded = 5
    pointsNeededAtLeast = 3
    numberIterations = 1000
    numberRestarts = 100
    distanceTolerated = (dims[0]+dims[1])/400.  #This sets the distance for the check of how much more beads are overlapping
    toleranceForShearing = 0.1				#set this variable to a small value, to suppress shearing
    maxDist = (dims[0]+dims[1])/2./7				#Beads that have no nearer neighbor of the other color than maxDist are not considered
    print len(s), len(d)
    s,d = preselectPoints(s,d,maxDist) #skipps all points that are far off any other point of the other color
    print len(s), len(d)
    pointsNeeded = min([len(s),len(d)])
    sc=copy.deepcopy(s)		#affineMatrix2DFromCorrespondingPoints changes the values given
    dc=copy.deepcopy(d)

    if min([len(s),len(d)]) < pointsNeededAtLeast:
        print 'To few points found, there might have been not enough pairs of matching points'
        return

    print 'begin Ransac'
    for i in range(numberRestarts):
        matchedPoints, additionalCandidatesS, additionalCandidatesD, indicesSCollection, indicesDCollection, trafo = performRansac(sc,dc,dims,pointsNeeded, numberIterations, distanceTolerated)
        maximalMatchedPoints = np.max(matchedPoints)
        print maximalMatchedPoints
        bestIndices = []
        if maximalMatchedPoints<=pointsNeeded:
            pointsNeeded -=1
            continue
        if maximalMatchedPoints > np.min(matchedPoints) or maximalMatchedPoints ==np.min([len(s),len(d)]):	# more combined points are found than points that are necessary for the calculation of the transformation
            indicesWithMaximalAlignment = np.where(matchedPoints==maximalMatchedPoints)[0]
            #breaks if the trafo has no shear
            
            #if there is just translation and rotation trafo[0,0] should be cos(phi), trafo[0,1] - sin(phi)..., therefore on can check if there is shearing
            
            for k in range(len(indicesWithMaximalAlignment)):
                bestTrafo = trafo[indicesWithMaximalAlignment[k]]
                if (np.abs((bestTrafo[0,0]**2+bestTrafo[0,1]**2)-1) < toleranceForShearing and np.abs((bestTrafo[1,0]**2+bestTrafo[1,1]**2)-1) <  toleranceForShearing):
                    bestIndices.append(indicesWithMaximalAlignment[k])
                    
            if len(bestIndices) == 0:
                print 'Start Iteration again (%i/%i), reason: minimum < maximum but shearing' %(i, numberRestarts)
            else:
                break
        else:
            for j in range(len(trafo)): # just the points used to calculate the transformation lay above each other, either there are just 3 pairs of beads, or the correct transformation has not been found jet
                if (toleranceForShearing > np.abs((trafo[j][0,0]**2+trafo[j][0,1]**2)-1) and np.abs((trafo[j][1,0]**2+trafo[j][1,1]**2)-1) <toleranceForShearing):
                    bestIndices.append(j)
#                     exitRans = True
                    break
#             if exitRans:
#                 break
            print 'Start Iteration again (%i/%i), reason: too few points found' %(i, numberRestarts)
    
    rmsList = []
    trafoList = []
    for i in (bestIndices):
        trafo = affineMatrix2DFromCorrespondingPoints(s[additionalCandidatesS[i]], d[additionalCandidatesD[i]], dims)
        trafoList.append(trafo)
        tt = np.dot(trafo, np.vstack([s[additionalCandidatesS[i]].T,np.ones(len(additionalCandidatesS[i]))])).T[:,:2] # fg transformed in bg coords
        rss = np.sum((tt-d[additionalCandidatesD[i]])**2)
        rms = np.sqrt(rss/len(additionalCandidatesD[i]))
        rmsList.append(rms)
        print rms, rss
    bestIdx = bestIndices[np.where(rmsList == np.min(rmsList))[0][0]]
    bestIdxSelection = np.where(rmsList == np.min(rmsList))[0][0]
    return trafoList[bestIdxSelection],s[additionalCandidatesS[bestIdx]],d[additionalCandidatesD[bestIdx]]

#     print maximalMatchedPoints
#     finalIndicesS = []
#     finalIndicesD = []
# 
#     indexWithMaximalAlignment = np.where(matchedPoints==maximalMatchedPoints)[0][0]
# 
#     finalIndicesS = additionalCandidatesS[indexWithMaximalAlignment]
#     finalIndicesD = additionalCandidatesD[indexWithMaximalAlignment]
# 
#     finalIndicesS = (np.reshape(np.array(finalIndicesS),-1))
#     finalIndicesD = (np.reshape(np.array(finalIndicesD),-1))
# 
#     finalTrafo = trafo[indexWithMaximalAlignment]
#     #finalTrafo, delR, delPhi = affineMatrix2DFromCorrespondingPoints(s[finalIndicesS], d[finalIndicesD], dims)
# 
#     print 'Found %i, Points for transformation' %len(finalIndicesS)
#     print finalTrafo
#     print finalTrafo[0,0]**2+finalTrafo[0,1]**2,finalTrafo[1,0]**2+finalTrafo[1,1]**2
#     #return finalTrafo,s[finalIndicesS],d[finalIndicesD], delR, delPhi
#     return finalTrafo,s[finalIndicesS],d[finalIndicesD]

def getDistanceMatrix(s,d):
    numberBeadsRed = d.shape[0]
    numberBeadsGreen = s.shape[0]
    #matrix is the distance matrix between all red and green points
    matrix = np.zeros((numberBeadsRed, numberBeadsGreen))
    for i in range(numberBeadsRed):
        for j in range(numberBeadsGreen):
            matrix[i,j] = np.sqrt((d[i,0]-s[j,0])**2+(d[i,1]-s[j,1])**2)

    return matrix

# def alignCandidates(s,d):
#     matrix = getDistanceMatrix(s,d)
#     numberBeadsRed = d.shape[0]
#     numberBeadsGreen = s.shape[0]
# 
#     m = Munkres()
#     indexes = m.compute(matrix)
# 
#     alignedd = []
#     aligneds = []
# 
#     total = 0
#     for row, column in indexes:
#         value = matrix[row][column]
#         total += value
# 
#         alignedd.append(d[row])
#         aligneds.append(s[column])
# 
#     return np.array(aligneds), np.array(alignedd)


def affineMatrix2DFromCorrespondingPoints(s, d, dims):
    #d = layer 0
    n = len(s)
    if len(s) < 3:
        print "at least three points required"
        return

    #dims = [0,0]
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
    #d.T*s = A * s.T*s
    #A = D.T*S*(S.T*S).I
    row3 = np.array([0,0,1])

    A = np.vstack([solx,soly,row3])
    R=np.matrix([[1, 0, deltaX],[0, 1, deltaY],[0, 0, 1]])
    Ap = np.dot(R,np.dot(A,R.I))

    return np.array([[Ap[0,0],Ap[0,1],Ap[0,2]],[Ap[1,0],Ap[1,1],Ap[1,2]],[Ap[2,0],Ap[2,1],Ap[2,2]]])

def affineMatrix2DFromCorrespondingPointsUseTotalLeastSquares(s, d, dims):
    B0 = np.hstack([s[:,0:2],d[:,0:1]])
    B0 = B0 - np.mean(B0,0)
    U0,sig0,Vt0 = np.linalg.svd(B0)

    B1 = np.hstack([s[:,0:2],d[:,1:2]])
    B1 = B1 - np.mean(B1,0)
    U1,sig1,Vt1 = np.linalg.svd(B1)


    A2 = np.zeros([3,3])

    index1 = np.where(sig1 == np.min(sig1))[0][0]
    index0 = np.where(sig0 == np.min(sig0))[0][0]

    facest0 = -Vt0[index0,0]/(Vt0[index0,2]*np.cos(np.arcsin(-Vt0[index0,1]/Vt0[index0,2])))
    phiest0 = np.arccos(-Vt0[index0,0]/(Vt0[index0,2]*facest0))

    phiest1 = np.arcsin(Vt1[index1,0]/(Vt1[index1,2]))
    facest1 = -Vt1[index1,1]/(Vt1[index1,2]*np.cos(phiest1))

    facest = (facest0 + facest1)/2
    phiest = (phiest0 + phiest1)/2


    A2[0,0] = facest * np.cos(phiest)

    A2[0,1] = np.sin(phiest)
    A2[1,0] = -np.sin(phiest)
    A2[1,1] = facest * np.cos(phiest)

    A2[0,2] = np.mean(d[:,0]) - (A2[0,0] * np.mean(s,0)[0] + A2[0,1] * np.mean(s,0)[1])
    A2[1,2] = np.mean(d[:,1]) - (A2[1,0] * np.mean(s,0)[0] + A2[1,1] * np.mean(s,0)[1])
    A2[2,2] = 1


    return A2

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
    import scipy
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
