import numpy as np
from scipy import integrate
import time
from matplotlib import pyplot
import coords
#import vigra
import scipy.stats
import math
from scipy.stats import spearmanr
from PIL import Image
import matplotlib.pyplot as plot

def contributionFunc(x,G,gamma, gamma_minus_one):
    return(max(G, min(G + 1, x*np.tan(gamma))) - min(G + 1, max(G, x*np.tan(gamma_minus_one))))

def calcContribution(point, gammas, j): #point i, list of all angles, j = angle contribution is calculated for
    return(integrate.quad(contributionFunc,point[0],point[0]+1,(point[1],gammas[j],gammas[j-1])))

def coloc(dataR, dataG):
    #returns heatmap for colocalization based on the angle in the red- green value plot
    #regions with low intensity are filterd out

    if dataR.shape != dataG.shape:
        print('images must have same shape')
        return 0

    tol=0.02
    dataB = np.zeros(dataR.shape)
    dataB[...,0]=np.tan(dataR[...,2]/dataG[...,1])
    dataB[...,0]=np.where((dataB[...,0]-np.pi/2.)**2>tol,0,dataB[...,0])
    maskG=np.where(dataG[...,1]<np.mean(dataG[...,1]),0,dataG[...,1])
    maskR=np.where(dataR[...,2]<np.mean(dataR[...,2]),0,dataR[...,2])
    from matplotlib import pyplot
    #pyplot.matshow(maskR)
    #pyplot.show()
    #pyplot.matshow(maskG)
    #pyplot.show()

    dataB[...,0]=dataB[...,0]*maskG*maskR
    dataB=dataB*255/np.max(dataB)
    dataB = np.array(dataB, dtype=np.uint8)
    print(np.mean(dataB[...,0]))


    plot.matshow(dataB[...,0])
    plot.show()

    return dataB
    #return dataB+dataR, dataG


def Colocdetection(dataR, dataG):
    #similar to the function coloc, but the angle is much wider, and can be set individually for the upper and lower bound
    threshR = 0.3*np.max(dataR)
    threshG = 0.3*np.max(dataG)
    combImg = dataR+dataG
    index = np.where(combImg[...,2]>threshR)
    indexres = np.where(combImg[index[0],index[1],1]>threshG)

    dataB = np.zeros((dataR.shape))

    for i in range(len(indexres[0])):
        x=index[0][indexres[0][i]]
        y=index[1][indexres[0][i]]
        dataB[x,y,0]=np.min([dataR[x,y,2],dataG[x,y,1]])

    dataB=dataB*255/np.max(dataB[x,y,0])
    dataB = np.array(dataB, dtype=np.uint8)
    #return dataB,np.zeros((dataR.shape))
    return dataB+dataR, dataG


def calcPearsonCorrelationCoeff(dataR, dataG):
    dataR=dataR[...,2]
    dataG=dataG[...,1]

    res=np.corrcoef(dataR.flatten(), dataG.flatten())[0,1]

    if np.isnan(res):
        res=0

    return res

def calcOverlapCoeff(dataR, dataG):

    dataR= dataR[...,2].astype(np.float64)
    dataG= dataG[...,1].astype(np.float64)


    numerator=np.sum((dataR*dataG).flatten())
    a=np.sum((dataR*dataR).flatten())
    b=np.sum((dataG*dataG).flatten())
    denominator=np.sqrt(a*b)
    #denominator=np.sqrt(sum(dataR*dataR)*sum(dataG*dataG))

    return numerator/denominator

def calcMandersColocalizationCoeffs(dataR, dataG):
    TG=np.max(dataG)*.1
    TR=np.max(dataR)*.1

    indR = np.where(dataG[:,:,1]>TG)
    indG = np.where(dataR[:,:,2]>TR)
    MR = np.sum(dataR[indR[0],indR[1],2])/float(np.sum((dataR[...,2]).flatten()))
    MG = np.sum(dataG[indG[0],indG[1],1])/float(np.sum((dataG[...,1]).flatten()))

    return MR, MG

def cropROI(data,ROI):
    [xmin,xmax,ymin,ymax]=ROI
    return data[xmin:xmax,ymin:ymax]

def gaussianRect(data, xmean,ymean,varx,vary):
    varx=int(varx)
    vary=int(vary)
    resultMatrix=np.zeros((6*varx+1,6*vary+1))
    for x in range(-3*varx,3*varx+1):
        for y in range(-3*vary,3*vary+1):
            resultMatrix[x,y] = data[xmean+x,ymean+y]*1/np.sqrt(2*np.pi*varx*vary)*np.exp(np.sqrt((x)**2+(y)**2)/(varx*vary))


    return resultMatrix

def createHeatmap2(dataR,dataG):
    dataRm = np.array(dataR[...,2],dtype=np.float32)
    dataRm = np.array(vigra.gaussianSmoothing(dataRm, 3),dtype=np.uint8)
    dataGm = np.array(dataG[...,1],dtype=np.float32)
    dataGm = np.array(vigra.gaussianSmoothing(dataGm, 3),dtype=np.uint8)

    dataRR = dataR[:,:,2] * dataR[:,:,2]
    dataGG = dataG[:,:,1] * dataG[:,:,1]

    dataRR = np.array(dataRR,dtype = np.float32)
    dataGG = np.array(dataGG,dtype = np.float32)

    dataRRm = np.array(vigra.gaussianSmoothing(dataRR, 3),dtype=np.uint8)
    dataGGm = np.array(vigra.gaussianSmoothing(dataGG, 3),dtype=np.uint8)

    dataRG = dataR[:,:,2] * dataG[:,:,1]
    dataRG = np.array(dataRG,dtype=np.float32)
    dataRGm = np.array(vigra.gaussianSmoothing(dataRG, 3),dtype=np.uint8)

    suspectedBackgroundR = scipy.stats.scoreatpercentile(dataR[...,2].flatten(), 50)
    suspectedBackgroundG = scipy.stats.scoreatpercentile(dataG[...,1].flatten(), 50)

    dataRwichtung = np.where(dataR[:,:,2]>suspectedBackgroundR,1,0)
    dataGwichtung = np.where(dataG[:,:,1]>suspectedBackgroundG,1,0)

    dataB = np.zeros((dataR.shape))
    dataB[:,:,0]=(dataRGm-dataRm*dataGm)*dataRwichtung*dataGwichtung


    Var1=np.sqrt(dataRRm-dataRm**2)
    Var2=np.sqrt(dataGGm-dataGm**2)

    dataRVarwichtung = np.where(Var1!=0,1,0)
    dataGVarwichtung = np.where(Var2!=0,1,0)

    Var1=np.where(Var1==0,10,Var1)      #dataRVarwichtung cutts them out anyway
    Var2=np.where(Var2==0,10,Var2)

    dataB[:,:,0]=dataB[:,:,0]/(Var1*Var2)*dataRVarwichtung*dataGVarwichtung
    #mmx = scipy.stats.scoreatpercentile(dataB[dataB>0].flatten(), 90)
    #if mmx.max() > 0:
    #    dataB[dataB[:,:,0]>mmx] = [mmx,0,0,0] # crop maximum at above percentile
    dataB = dataB *255./np.max(dataB)
    dataB = np.array(dataB, dtype=np.uint8)
    dataR=dataR+ dataB
    return dataB
    #return dataB,np.zeros(dataR.shape, dtype=np.uint8)

def createHeatmap(dataR, dataG):
    np.seterr(invalid='ignore')
    numberBinsX=200
    numberBinsY=200
    dims = dataR.shape
    width=int(dims[0]/numberBinsX)
    height=int(dims[1]/numberBinsY)
    ValueHeatmap = np.zeros((numberBinsX,numberBinsY))
    N1=5
    N2=5
    dataB=np.zeros(dims)
    for shift in range(N1):
        for shift2 in range(N2):
            for i in range(numberBinsX-1):
                for j in range(numberBinsY-1):
                    [xmin, xmax,ymin,ymax] = [int(i*width+shift),int((i+1)*width+shift), int(j*height+shift2), int((j+1)*height+shift2)]
                    temp=calcPearsonCorrelationCoeff(dataR[xmin:xmax,ymin:ymax],dataG[xmin:xmax,ymin:ymax])

                    if not np.isnan(temp) and temp > 0:
                        ValueHeatmap[i,j]=temp
                        #dataR[xmin:xmax,ymin:ymax,0]=dataR[xmin:xmax,ymin:ymax,0]+temp*255./(N1*N2)
                        dataB[xmin:xmax,ymin:ymax,0]=dataB[xmin:xmax,ymin:ymax,0]+temp*np.mean(dataR[xmin:xmax,ymin:ymax,2])/(N1*N2)
                        #dataG[xmin:xmax,ymin:ymax,0]=dataG[xmin:xmax,ymin:ymax,0]+temp*255


    #dataB = vigra.RGBImage(np.array(dataB,dtype=np.float32))
    #dataB = np.array(vigra.gaussianSmoothing(dataB, 1),dtype=np.uint8)
    dataB = dataB *255./np.max(dataB)
    dataB = np.array(dataB, dtype=np.uint8)
    dataR=dataR+ dataB
    return [dataR,dataG]

def createColocHeatmap(dataR, dataG):
    np.seterr(invalid='ignore')
    borderWidth = 10
    dataB = np.zeros(dataR.shape)
    dataRc = np.zeros(dataR.shape[0:2])
    dataGc = np.zeros_like(dataRc)
    dataRc[borderWidth:-borderWidth,borderWidth:-borderWidth]=dataR[borderWidth:-borderWidth,borderWidth:-borderWidth,2]
    dataGc[borderWidth:-borderWidth,borderWidth:-borderWidth]=dataG[borderWidth:-borderWidth,borderWidth:-borderWidth,1]


    indicesRc = np.array(np.where(dataRc>0))
    indicesGc = np.array(np.where(dataGc>0))

    for i in range(indicesRc.shape[1]):
        uli = indicesRc[0,i] + borderWidth #upper limit
        lli = indicesRc[0,i] - borderWidth #lower limit
        ulj = indicesRc[1,i] + borderWidth
        llj = indicesRc[1,i] - borderWidth
        res=np.corrcoef(dataRc[lli:uli,llj:ulj].flatten(), dataGc[lli:uli,llj:ulj].flatten())[0,1]
        if np.isnan(res):
            res=0
        dataB[indicesRc[0,i],indicesRc[1,i],0] = res
    for i in range(indicesGc.shape[1]):
        uli = indicesGc[0,i] + borderWidth #upper limit
        lli = indicesGc[0,i] - borderWidth #lower limit
        ulj = indicesGc[1,i] + borderWidth
        llj = indicesGc[1,i] - borderWidth
        res=np.corrcoef(dataRc[lli:uli,llj:ulj].flatten(), dataGc[lli:uli,llj:ulj].flatten())[0,1]
        if np.isnan(res):
            res=0
        #print res
        dataB[indicesGc[0,i],indicesGc[1,i],0] = res

 #   for i in range(borderWidth, dataR.shape[0] - borderWidth):
 #       print i
 #       for j in range(borderWidth, dataR.shape[1] - borderWidth):
 #           uli = i + borderWidth #upper limit
 #           lli = i - borderWidth #lower limit
 #           ulj = j + borderWidth
 #           llj = j - borderWidth
 #           res=np.corrcoef(dataR[lli:uli,llj:ulj].flatten(), dataG[lli:uli,llj:ulj].flatten())[0,1]
 #           if np.isnan(res):
 #               res=0
 #           dataB[i,j,0] = res

    factor = np.mean(np.where(dataR>0))/np.max(np.abs(dataB))
    dataB = np.abs(dataB) * factor

    import matplotlib.pyplot as plot
    plot.matshow(dataB[...,0])
    plot.show()
    return dataB

def getAngleHistogram(dataR, dataG, dims, factor=1):
    import pylab
    K=64+64
    '''dataR=dataR[...,2]
    dataG=dataG[...,1]
    a=dataR[np.where((dataR>10) & (dataG>10))].flatten().reshape(-1,1)[::5]
    b=dataG[np.where((dataR>10) & (dataG>10))].flatten().reshape(-1,1)[::5]
    if len(a) > 0 and len(b) > 0:
        c = np.zeros_like(a)*255
        colors=np.hstack([a,b,c])/255.0
        print a,b
        print
        print
        print np.max(a)
        print np.min(a)
        print np.max(b)
        print np.min(b)
        pylab.scatter(a.flatten(),b.flatten(),color=colors)
        pyplot.xlim((-2,256))
        pyplot.ylim((-2,256))
        pyplot.xlabel('Intensity red channel')
        pyplot.ylabel('Intensity green channel')
        pylab.show()

    angles = np.linspace(0,90,K)'''

    mergedPoints = np.array(mergePoints3(dataR,dataG,factor))
    '''histo = np.zeros(K)

    for j in range(K):
        print j
        start = time.time()
        matrixc = np.zeros((256,256))+2
        for i in range(mergedPoints.shape[0]):
            if matrixc[mergedPoints[i,2],mergedPoints[i,3]] == 2:
                cij,_ = calcContribution([mergedPoints[i,2],mergedPoints[i,3]],angles,j+1)
                matrixc[mergedPoints[i,2],mergedPoints[i,3]] = cij
            else:
                cij = matrixc[mergedPoints[i,2],mergedPoints[i,3]]
            #cij,_ = calcContribution([mergedPoints[i,2],mergedPoints[i,3]],angles,j+1)      #use of j+1, because j is the upper bound

            histo[j] += max(mergedPoints[i,2],mergedPoints[i,3])*cij


        print time.time()-start

    for i in range(K):
        print "Angle: %3.2f, Number entries: %i" %(angles[i]*180/np.pi, histo[i])

    phi=np.linspace(0,90,128)
    pyplot.plot(phi,histo)
    pyplot.show()
    pyplot.scatter(mergedPoints[:,2],mergedPoints[:,3])
    '''
    return 0


def CBC_ROI2(tree, loc, locs, r, R):
    idxR = set(tree.query_ball_point(loc[:2],R))
    idxr = set(tree.query_ball_point(loc[:2],r))
    roi = locs[list(idxR-idxr)]
    L = roi.shape[0]
    if L == 0:
        roi = np.zeros([1,5])
    return roi, L

def CBC_ROI3(tree, loc, r, R):
    idxR = set(tree.query_ball_point(loc[:2],R))
    idxr = set(tree.query_ball_point(loc[:2],r))
    return len(idxR-idxr)

def CBC(locsAo,locsBo, dims, factor, pixelsize = 100):
    import numpy as np
    import scipy.stats
    import scipy as sp
    from scipy import spatial
    locsA = np.copy(locsAo)
    locsB = np.copy(locsBo)
    locsA[:,:2] = locsA[:,:2] * pixelsize
    locsB[:,:2] = locsB[:,:2] * pixelsize
    R = 100#int(5*np.sqrt(pixelsize))
    Int = 5
    thresholdLA = 20
    thresholdLB = 20
    Weight = 1
    import time
    st = time.time()

    listcorr = []
#     locsA[:,:2]*=pixelsize
#     locsB[:,:2]*=pixelsize
    treeA = spatial.cKDTree(locsA[:,:2],10)
    treeB = spatial.cKDTree(locsB[:,:2],10)

    for i in range (0,(len(locsA))):
        print(i, len(locsA))
        # selection

        loc=locsA[i,:]
        LA = CBC_ROI3(treeA, loc, 0 ,R)
        LB = CBC_ROI3(treeB, loc, 0 ,R)

        # correlation
        if (LA>thresholdLA) and (LB>thresholdLB):
            TraceA=np.zeros([Int,2])
            TraceA[:,0]=range (0,R-np.mod(R,Int),(R/Int))
            TraceB=np.zeros([Int,2])
            TraceB[:,0]=range (0,R-np.mod(R,Int),(R/Int))
            for n in range (1,Int):
                m=n*(R/float(Int))
                o=0
                la=CBC_ROI3(treeA,loc,0,m)#CBC_ROI(roiA,o,m)
                lb=CBC_ROI3(treeB,loc,0,m)#CBC_ROI(roiB,o,m)
                TraceA[n,1]=(la/(np.pi*m*m))*(np.pi*R*R/(LA))
                TraceB[n,1]=(lb/(np.pi*m*m))*(np.pi*R*R/(LB))

            Corr=sp.stats.spearmanr(TraceA[:,1], TraceB[:,1])
            if np.isnan(Corr[0]):
                Corr = [0,0]
            #Weighting
            #Dist=np.exp(-Weight*(min(roiB[:,5]))/R)
            NN = 0
            counter = 0
            while NN == 0:
                NN = treeB.query(loc[:2],2+counter)[0][counter]
                counter += 1
            Dist=np.exp(-Weight*(NN)/R)
            Coloc=Corr[0]*Dist
        else:
           Coloc=0
        #saving
        listcorr.append(Coloc)
        #out_file.write(str(loc[0])+' '+str(loc[1])+' '+str(loc[2])+' '+str(int(loc[3]))+' '+str(Coloc)+"\n")

    dataB = np.zeros([dims[1]*factor, dims[0]*factor,4])
    for i in range(len(listcorr)):

        if listcorr[i]>0:
            try: 
                dataB[int(factor * locsA[i,1]/pixelsize),int(factor * locsA[i,0]/pixelsize),0] = listcorr[i]*255
            except IndexError:
                tmp = 1
        else:
            try:
                dataB[int(factor * locsA[i,1]/pixelsize),int(factor * locsA[i,0]/pixelsize),1] = (-listcorr[i])*255
            except IndexError:
                tmp = 1
    dataB = np.array(dataB, dtype=np.uint8)

    print("Total time for colocalization detection:",time.time() - st)
    plot.matshow(dataB[...,0])
    plot.show()

    return dataB
