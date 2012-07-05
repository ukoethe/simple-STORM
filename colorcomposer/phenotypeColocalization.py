import os
import Image
import numpy as np
import mahotas
import mahotas.surf
import glob,sys,os
import vigra

def calcPearsonCorrelationCoeff(dataR, dataG):
    '''
    
    meanR=np.mean(dataR)
    meanG=np.mean(dataG)
    
    dataRn= dataR-meanR
    dataGn= dataG-meanG
    np.cov(dataR,dataG)
    numerator=np.sum((dataGn*dataRn).flatten())
    #numerator = np.sum(np.cov(dataR,dataG).flatten())
    #denominator=np.sqrt(np.var(dataR)*np.var(dataG))
    denominator=np.sqrt(sum((dataRn*dataRn).flatten())*sum((dataGn*dataGn).flatten()))
    '''
    
    return np.corrcoef(dataR.flatten(), dataG.flatten())[0,1]

def calcOverlapCoeff(dataR, dataG):

    
    numerator=np.sum((dataR*dataG).flatten())
    a=np.sum((dataR*dataR).flatten())
    b=np.sum((dataG*dataG).flatten())
    denominator=np.sqrt(a*b)
    #denominator=np.sqrt(sum(dataR*dataR)*sum(dataG*dataG))
    
    return numerator/denominator

def calcMandersColocalizationCoeffs(dataR, dataG):
    TG=10
    TR=10
    
    indR = np.where(dataG[:,:]>TG)
    indG = np.where(dataR[:,:]>TR)
    MR = np.sum(dataR[indR[0],indR[1]])/float(np.sum((dataR).flatten()))
    MG = np.sum(dataG[indG[0],indG[1]])/float(np.sum((dataG).flatten()))
    
    return MR, MG


classes={'allstars': 1,'brefeldinA': 2,'control': 3, 'gbf1': 4, 'nocodazole': 5}

print 'hihihi'
#folders = ["/home/herrmannsdoerfer/master/workspace/segmentation/Set1/" , "/home/herrmannsdoerfer/master/workspace/segmentation/Set2/" ,"/home/herrmannsdoerfer/master/workspace/segmentation/Set3/"] 
folders =glob.glob("/home/herrmannsdoerfer/Desktop/SampleCellsUnedited/*")

list_filenames=[]
list_labels=[]
for folder in folders:
    _,marker=folder.split('-')
    label=classes[marker]
    folder=folder+'/confocal'
    files=sorted(glob.glob(os.path.join(folder,'*.tif')))
    list_filenames+=files
    list_labels+=[label]*len(files)

print len(list_filenames),' ',len(list_labels)
count=0

list_PearsonCoeff = []
list_OverlapCoeff = []
list_MandersCoeff = []

list_raw=[]
for k in range(len(list_filenames)/2):
    imgR=vigra.impex.readImage(list_filenames[count]).swapaxes(0,1).view(np.ndarray).astype(np.float32).squeeze()
    imgG=vigra.impex.readImage(list_filenames[count+1]).swapaxes(0,1).view(np.ndarray).astype(np.float32).squeeze()
    list_raw.append((imgR,imgG,list_labels[count]))
    count+=2
    list_PearsonCoeff.append(calcPearsonCorrelationCoeff(imgR,imgG))
    list_OverlapCoeff.append(calcOverlapCoeff(imgR,imgG))
    list_MandersCoeff.append((calcMandersColocalizationCoeffs(imgR,imgG)))

mean_PearsonCoeff = [[],[],[],[],[]]
mean_OverlapCoeff = [[],[],[],[],[]]
mean_MandersCoeff = [[],[],[],[],[]]
for i in range(len(list_PearsonCoeff)):
    aktlabel = list_labels[i*2]-1
    mean_PearsonCoeff[aktlabel].append(list_PearsonCoeff[i])
    mean_OverlapCoeff[aktlabel].append(list_OverlapCoeff[i])
    mean_MandersCoeff[aktlabel].append(list_MandersCoeff[i])
    
print 'hi'
for i in range(5):
    print 'PearsonCoeff class %i, mean: %3.3f, variance: %3.6f' %(i+1,np.mean(mean_PearsonCoeff[i]),np.var(mean_PearsonCoeff[i]))
print
print
for i in range(5):
    print 'OverlapCoeff class %i, mean: %3.3f, variance: %3.6f' %(i+1,np.mean(mean_OverlapCoeff[i]), np.var(mean_OverlapCoeff[i]))

print
print

for i in range(5):
    npMandersCoeff=np.array(mean_MandersCoeff[i])
    print 'MandersCoeff class %i, means: %3.3f %3.3f, variances: %3.6f %3.6f' %(i+1,np.mean(npMandersCoeff[:][0]),np.mean(npMandersCoeff[:][1]),np.var(npMandersCoeff[:][0]),np.var(npMandersCoeff[:][1]))


folders =glob.glob("/home/herrmannsdoerfer/Desktop/SampleCellsUnedited/*")

list_filenames=[]
list_labels=[]
for folder in folders:
    _,marker=folder.split('-')
    label=classes[marker]
    folder=folder+'/storm'
    files=sorted(glob.glob(os.path.join(folder,'*.tif')))
    list_filenames+=files
    list_labels+=[label]*len(files)

print len(list_filenames),' ',len(list_labels)
count=0

list_PearsonCoeff = []
list_OverlapCoeff = []
list_MandersCoeff = []

import Image

list_raw=[]
for k in range(len(list_filenames)):
    imgR=vigra.impex.readImage(list_filenames[count], index = 0).swapaxes(0,1).view(np.ndarray).astype(np.float32)
    imgG=vigra.impex.readImage(list_filenames[count], index = 1).swapaxes(0,1).view(np.ndarray).astype(np.float32)
    list_raw.append((imgR,imgG,list_labels[count]))
    count +=1
    list_PearsonCoeff.append(calcPearsonCorrelationCoeff(imgR,imgG))
    list_OverlapCoeff.append(calcOverlapCoeff(imgR,imgG))
    list_MandersCoeff.append((calcMandersColocalizationCoeffs(imgR,imgG)))

mean_PearsonCoeff = [[],[],[],[],[]]
mean_OverlapCoeff = [[],[],[],[],[]]
mean_MandersCoeff = [[],[],[],[],[]]
for i in range(len(list_PearsonCoeff)):
    aktlabel = list_labels[i]-1
    mean_PearsonCoeff[aktlabel].append(list_PearsonCoeff[i])
    mean_OverlapCoeff[aktlabel].append(list_OverlapCoeff[i])
    mean_MandersCoeff[aktlabel].append(list_MandersCoeff[i])
    
print 'hi'
for i in range(5):
    print 'PearsonCoeff class %i, mean: %3.3f, variance: %3.6f' %(i+1,np.mean(mean_PearsonCoeff[i]),np.var(mean_PearsonCoeff[i]))
print
print
for i in range(5):
    print 'OverlapCoeff class %i, mean: %3.3f, variance: %3.6f' %(i+1,np.mean(mean_OverlapCoeff[i]), np.var(mean_OverlapCoeff[i]))

print
print
for i in range(5):
    npMandersCoeff=np.array(mean_MandersCoeff[i])
    print 'MandersCoeff class %i, means: %3.3f %3.3f, variances: %3.6f %3.6f' %(i+1,np.mean(npMandersCoeff[:,0]),np.mean(npMandersCoeff[:,1]),np.var(npMandersCoeff[:,0]),np.var(npMandersCoeff[:,1]))