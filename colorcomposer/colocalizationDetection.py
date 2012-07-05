import numpy as np
from scipy import integrate
import time
from matplotlib import pyplot
import coords
import vigra
import scipy.stats




def Angles(K):
    gammas = np.zeros(K+1)
    for i in range(K+1):
        if i <= K/2.:
            gammas[i] = np.arctan(i/(float(K)/2.))
        else:
            gammas[i] = 90/180.*np.pi - gammas[K-i]
    
    return gammas

def contributionFunc(x,G,gamma, gamma_minus_one):
    return(max(G, min(G + 1, x*np.tan(gamma))) - min(G + 1, max(G, x*np.tan(gamma_minus_one))))

def calcContribution(point, gammas, j): #point i, list of all angles, j = angle contribution is calculated for
    return(integrate.quad(contributionFunc,point[0],point[0]+1,(point[1],gammas[j],gammas[j-1])))

def Colocdetection(dataR, dataG):
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
    
def mergePoints2(dataR, dataG, factor):
    image=dataR+dataG
    mergedPoints =[]
    thresh = 25
    ind = np.where(image!=0)
    counter=0
    while counter < ind[0].shape[0]-1:

        if image[ind[0][counter],ind[1][counter],2] >= thresh and image[ind[0][counter],ind[1][counter],1] >= thresh:
            mergedPoints.append((ind[0][counter], ind[1][counter], image[ind[0][counter],ind[1][counter],2],image[ind[0][counter],ind[1][counter],1]))
            if ind[0][counter] == ind[0][counter+1] and ind[1][counter] == ind[1][counter+1]:   #if one pixel has 2 colors one index must be skipped, otherwise the pixel would be counted twice
                counter +=1
        counter +=1
    
    if not(mergedPoints[-1][0] == ind[0][-1] and mergedPoints[-1][1] == ind[1][-1]):    #if the last point hast NOT two channels it will be added
        if image[ind[0][-1],ind[1][-1],2] > thresh and image[ind[0][-1],ind[1][-1],1] > thresh:
            mergedPoints.append((ind[0][-1], ind[1][-1], image[ind[0][-1],ind[1][-1],2],image[ind[0][-1],ind[1][-1],1]))
    return mergedPoints


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
    TG=10
    TR=10
    
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


def createHeatmap3(dataR,dataG):
    '''rand=5
    for x in range(3*rand,dataR.shape[0]-3*rand):
        for y in range(3*rand,dataR.shape[1]-3*rand):
            dataRG=gaussianRect(dataR[...,2],x,y,rand,rand)
            dataGG=gaussianRect(dataG[...,2],x,y,rand,rand)
            temp = np.corrcoef(dataRG.flatten(), dataGG.flatten())[0,1]
            if np.isnan(temp):
                temp=0
            dataR[x,y,0]=temp
    dataR[...,0]=dataR[...,0]*255/np.max(dataR[...,0])'''
    return dataR, dataG
            
            
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
    return dataR, dataG
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

gammas = Angles(128)
print gammas.shape
#print calcContribution([0,255],gammas,1)

def coordinateBasedColocalization(dataR, dataG):
    
    originaldataR = np.array(dataR)
    originaldataG = np.array(dataG)
    maximumR = np.max(dataR)
    dataR[:,:,2]=np.where(dataR[:,:,2]>=maximumR*0.1,255,0)
    Samsung-Galaxy-Note
    
    dataRnp = vigra.RGBImage(np.array(dataR, dtype = np.float32))    
    labelsR = np.zeros((dataR.shape[0], dataR.shape[1]))   
    labelsR = (vigra.analysis.labelImageWithBackground(dataRnp[...,2], 8,background_value=0))
    
    maximumG = np.max(dataG)
    dataG[:,:,1]=np.where(dataG[:,:,1]>=maximumG*0.1,255,0)
    
    dataGnp = vigra.RGBImage(np.array(dataG, dtype = np.float32))    
    labelsG = np.zeros((dataR.shape[0], dataR.shape[1]))   
    labelsG = (vigra.analysis.labelImageWithBackground(dataGnp[...,1], 8,background_value=0))
    
    '''pyplot.matshow(dataR[...,2],1)
    pyplot.matshow(dataG[...,1],2)
    pyplot.matshow(labelsR,3)
    pyplot.matshow(labelsG,4)
    pyplot.show()
    '''
    
    numberClassesR = np.max(labelsR)
    numberClassesG = np.max(labelsG)
    
    classesR = []
    classesR.append([])
    for i in range(1,numberClassesR + 1):
        coordinatesOfCurrentClass = np.where(labelsR == i,1,0)
        classesR.append(coordinatesOfCurrentClass)
        
    classesG = []
    classesG.append([])
    for i in range(1,numberClassesG + 1):
        coordinatesOfCurrentClass = np.where(labelsG == i,1,0)
        classesG.append(coordinatesOfCurrentClass)
        
    
    DAMatrix = np.zeros(dataR.shape)
    DBMatrix = np.zeros(dataR.shape)
      
    radius = 10

    for i in range(radius,dataR.shape[0]-radius):             
        for j in range(radius,dataR.shape[1]-radius):
            print i,j
            actuallClass = labelsR[i,j]
            if actuallClass != 0:
                Naa = 0
                for x in range(i-radius,i+radius):
                    for y in range(j-radius, j+radius):
                        if (x-i)**2+(y-j)**2<=radius**2:
                            Naa = Naa + classesR[actuallClass][x,y]
            
                if Naa!=0:
                    Rmax=0
                    coords = np.where(classesR[actuallClass] == 1)
                    for n in range(len(coords[0])):
                        tempMax = np.sqrt((coords[0][n]-i)**2+(coords[1][n]-j)**2)
                        if tempMax>Rmax:
                            Rmax=tempMax
                            
                    Naamax = 0
                    Rmax = int(Rmax)
                    for x in range(i-Rmax,i+Rmax):
                        for y in range(j-radius, j+radius):
                            if (x-i)**2+(y-j)**2<=radius**2:
                                Naamax = Naamax + classesR[actuallClass][x,y]
                    Samsung-Galaxy-Note
                    
                    #######Achtung sollte alle classen durchlaufen!!!!
                    Rminb = int(100)
                    coords = np.where(classesG[actuallClass] == 1)
                    for n in range(len(coords[0])):
                        tempMin = np.sqrt((coords[0][n]-i)**2+(coords[1][n]-j)**2)
                        if tempMin<Rminb:
                            Rminb=tempMin
                    if Rminb<2:
                        print i,j            
                    
                    
                    DAMatrix[i,j,2] = Naa / (np.pi*radius**2)* np.pi*Rmax**2/Naamax*np.exp(-Rminb)
                else:
                    DAMatrix[i,j,2] = 0
                
                
                Nab = 0
    
                for x in range(i-radius,i+radius):
                    for y in range(j-radius, j+radius):
                        if (x-i)**2+(y-j)**2<=radius**2:
                            for k in range(1,numberClassesG+1):
                                Nab = Nab + classesG[k][x,y]
                
                
                if Nab!=0:                              #if there is no class B in the radius this block is skipped
                    Rmax=0
                    coords = np.where(classesG[actuallClass] == 1)
                    for n in Samsung-Galaxy-Noterange(len(coords[0])):
                        tempMax = np.sqrt((coords[0][n]-i)**2+(coords[1][n]-j)**2)
                        if tempMax>Rmax:
                            Rmax=tempMax
                            
                    Rmaxb=Rmax
                            
                    Nabmax = 0
                    Rmax=int(Rmax)
                    for x in range(i-Rmax,i+Rmax):
                        for y in range(j-radius, j+radius):
                            if (x-i)**2+(y-j)**2<=radius**2:
                                for k in range(1,numberClassesG+1):
                                    Nabmax = Nabmax + classesG[k][x,y]
                    DBMatrix[i,j,1] = Nab / (np.pi*radius**2)* np.pi*Rmax**2/Nabmax
                else:
                    DBMatrix[i,j,1] = 0   #if there is no class B in the radius
                
                dataR[:,:,2]=np.where(dataR[:,:,2]>=maximumR*0.1,255,0)
            else:
                DAMatrix[i,j,2] = 0
                DBMatrix[i,j,1] = 0
            
            
    dataR=np.where(dataR<0,0,dataR)
    dataG=np.where(dataG<0,0,dataG)
    
    dataR=DAMatrix*255./np.max(DAMatrix)
    dataG=DBMatrix*255./np.max(DBMatrix)
    
    dataRm = np.array(dataR[...,2],dtype=np.float32)
    dataRm = np.array(vigra.gaussianSmoothing(dataRm, 3),dtype=np.uint8)
    dataGm = np.array(dataG[...,1],dtype=np.float32)
    dataGm = np.array(vigra.gaussianSmoothing(dataGm, 3),dtype=np.uint8)
    print dataRm.max()
    dataRG = dataR[:,:,2] * dataG[:,:,1]
    dataRG = np.array(dataRG,dtype=np.float32)
    dataRGm = np.array(vigra.gaussianSmoothing(dataRG, 3),dtype=np.uint8)
    
    dataRwichtung = np.where(dataR[:,:,2]>np.max(dataR)/10.,dataR[:,:,2],0)
    dataGwichtung = np.where(dataG[:,:,1]>np.max(dataG)/10.,dataG[:,:,1],0)
    
    dataB = np.zeros((dataR.shape))
    dataB[:,:,0]=(dataRGm-dataRm*dataGm)*dataRwichtung*dataGwichtung
    
    
    dataB = dataB *255./np.max(dataB)
    
    dataB = np.array(dataB, dtype=np.uint8)
    dataR=originaldataR+ dataB
    dataR = np.array(dataR,dtype=np.uint8)
    dataG = np.array(originaldataG,dtype=np.uint8)
    #return dataB, np.zeros(originaldataR.shape)
    return dataR, dataG
        

    
       