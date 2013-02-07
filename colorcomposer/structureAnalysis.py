import PyQt4
from PyQt4 import QtGui, QtCore
import sys
import os.path
import csv
import worker
import vigra
import numpy as np
import scipy.spatial.distance as scsd
#import ui_mainwindowstructureanalysis
from ui_mainwindowstructureanalysis import *
import matplotlib.pyplot as plt
import coords
import Image
from os.path import basename
import colocalizationDetection
import copy


class ImageData:
    def __init__(self):
        self.npimage = []      # numpy-images
        self.coords = []      # coordinate lists
        self.color = []        # color of the images
        self.dims   = []
        self.labels = []
        self.numberclasses = 0
        self.numberElementsOfEachClass = []
        self.meanOfEachClass = []
        self.meanx = []
        self.meany = []
        self.maximalValue = 0
        self.listOfCoordinatesOfEachClass = []
        self.cov_matrix = []
    
    def findClasses(self, factor):
        self.numberclasses = np.max(self.labels)
        self.numberElementsOfEachClass = np.zeros((self.numberclasses + 1))
        self.meanx = np.zeros((self.numberclasses + 1))
        self.meany = np.zeros((self.numberclasses + 1))
        for i in range(self.labels.shape[0]):
            for j in range(self.labels.shape[1]):
                #print "[",i,",",j,"]=",self.labels[i,j]
                self.numberElementsOfEachClass[self.labels[i,j]] += 1
                self.meanx[self.labels[i,j]] += j
                self.meany[self.labels[i,j]] += i
        False
        for i in range(1, self.numberclasses + 1):
            self.meanx[i] = self.meanx[i] / float(self.numberElementsOfEachClass[i]) / float(factor)
            self.meany[i] = self.meany[i] / float(self.numberElementsOfEachClass[i]) / float(factor)
        
    def findCovarianceMatrices(self, factor):
        for i in range(0, self.numberclasses + 1):
            self.listOfCoordinatesOfEachClass.append([])
        
        for i in range(self.labels.shape[0]):
            for j in range(self.labels.shape[1]):
                self.listOfCoordinatesOfEachClass[self.labels[i,j]].append([i/float(factor),j/float(factor)])
        
        self.cov_matrix.append([])
        for i in range(1, self.numberclasses + 1):
            print i
            self.cov_matrix.append(np.cov(self.listOfCoordinatesOfEachClass[i],y=None,rowvar=False))
            print self.cov_matrix
    
    def cutOff(self, number, factor):
        
        self.numberElementsOfEachClass = np.zeros((self.numberclasses + 1))

        for i in range(self.labels.shape[0]):
            for j in range(self.labels.shape[1]):
                print i,j
                self.numberElementsOfEachClass[self.labels[i,j]] += 1
        
        toSmallCluster = []
        for i in range(1, self.numberclasses+1):
            if number >= self.numberElementsOfEachClass[i]:
                toSmallCluster.append(i)
        
        for i in toSmallCluster:
            for j in range(int(self.numberElementsOfEachClass[i])):
                print i,j
                self.npimage[self.listOfCoordinatesOfEachClass[i][j][0]*factor,self.listOfCoordinatesOfEachClass[i][j][1]*factor,:] = 0 
                
              
    def calculateEverything(self, number, factor, indexPic, Self):
        #Find classes and mean
        self.numberclasses = np.max(self.labels)
        self.numberElementsOfEachClass = np.zeros((self.numberclasses + 1))
        self.meanx = np.zeros((self.numberclasses + 1))
        self.meany = np.zeros((self.numberclasses + 1))
        for i in range(self.labels.shape[0]):
            for j in range(self.labels.shape[1]):
                #print "[",i,",",j,"]=",self.labels[i,j]
                self.numberElementsOfEachClass[self.labels[i,j]] += 1
                self.meanx[self.labels[i,j]] += j
                self.meany[self.labels[i,j]] += i
        
        for i in range(1, self.numberclasses + 1):
            self.meanx[i] = self.meanx[i] / float(self.numberElementsOfEachClass[i]) / float(factor)
            self.meany[i] = self.meany[i] / float(self.numberElementsOfEachClass[i]) / float(factor)
        Self.main_window.statusBar().showMessage('Classes found')
        
        #Get list of coordinates
        for i in range(0, self.numberclasses + 1):
            self.listOfCoordinatesOfEachClass.append([])
        
        for i in range(self.labels.shape[0]):
            for j in range(self.labels.shape[1]):
                self.listOfCoordinatesOfEachClass[self.labels[i,j]].append([i/float(factor),j/float(factor)])
        Self.main_window.statusBar().showMessage('Coordinates listed')
        
        #cut off to small blobs
        toSmallCluster = []
        for i in range(1, self.numberclasses+1):
            if number >= self.numberElementsOfEachClass[i]:
                toSmallCluster.append(i)
        
        for i in toSmallCluster:
            for j in range(int(self.numberElementsOfEachClass[i])):
                print i,j
                self.npimage[self.listOfCoordinatesOfEachClass[i][j][0]*factor,self.listOfCoordinatesOfEachClass[i][j][1]*factor,:] = 0 
        Self.main_window.statusBar().showMessage('Small blobs deleted')
        
        if toSmallCluster != []:
            imag = self.npimage
            imag = vigra.RGBImage(np.array(imag, dtype = np.float32))         
            self.labels = np.zeros((imag.shape[0], imag.shape[1]))
            self.labels = (vigra.analysis.labelImageWithBackground(imag[:,:,2-indexPic], 8,background_value=0))
            self.findClasses(factor)
             #Get list of coordinates
            for i in range(0, self.numberclasses + 1):
                self.listOfCoordinatesOfEachClass.append([])
            
            for i in range(self.labels.shape[0]):
                for j in range(self.labels.shape[1]):
                    self.listOfCoordinatesOfEachClass[self.labels[i,j]].append([i/float(factor),j/float(factor)])
            Self.main_window.statusBar().showMessage('Coordinates listed')
            Self.main_window.statusBar().showMessage('Connected Components recalculated')
                
        #Find cov matrix
        self.cov_matrix.append([])
        for i in range(1, self.numberclasses + 1):
            print i
            self.cov_matrix.append(np.cov(self.listOfCoordinatesOfEachClass[i],y=None,rowvar=False))
            print self.cov_matrix
        Self.main_window.statusBar().showMessage('Covariance Matrix calculated')    
        
class CursorGraphicsScene(QtGui.QGraphicsScene):
    '''Add a cursor (ellipseItem) to the QGraphicsView with appropriate 
    MouseEvent-Handling'''
    def __init__(self):
        self.lastCursorPosition=[0,0]
        QtGui.QGraphicsScene.__init__(self)
        self.cursorRadius = 10.
        pen = QtGui.QPen(QtGui.QColor("white"))
        self.cursor = self.addEllipse(0,0,2.*self.cursorRadius,2.*self.cursorRadius, pen)
        self.cursor.setZValue(10.) # always foreground

        
    def mouseReleaseEvent(self, event):
        QtGui.QGraphicsScene.mouseReleaseEvent(self, event) # forward to individual items
        self.cursor.setPos(event.scenePos()-QtCore.QPointF(self.cursorRadius,self.cursorRadius)) # move cursor
        self.emit(QtCore.SIGNAL("cursorMoved(float, float, float)"), event.scenePos().x(), event.scenePos().y(), self.cursorRadius)
        self.lastCursorPosition=[event.scenePos().x(),event.scenePos().y()]
        
    
    def getCursorPosition(self):
        return self.lastCursorPosition
        
    def setCursorRadius(self, radius):
        self.cursorRadius = radius
        rect = self.cursor.rect()
        rect.setHeight(2.*radius)
        rect.setWidth(2.*radius)
        self.cursor.setRect(rect)
        

class structureAnalysis(QtCore.QObject):
    def __init__(self):
        QtCore.QObject.__init__(self)
        self.main_window = QtGui.QMainWindow()
        self.ui = Ui_MainWindowStructureAnalysis()
        self.ui.setupUi(self.main_window)
        self.main_window.show()

        self.Images  = []
        
        self.m_npimages = []      # numpy-images
        self.m_coords = []      # coordinate lists
        self.m_colors = []        # list of colors for the images
        self.m_factor = 8.         # factor for scaling coordinates to an image
        self.m_dims   = []
        self.m_labels = []
        self.m_qimage = QtGui.QImage()
        self.scene = CursorGraphicsScene()
        self.ui.graphicsView.setScene(self.scene)
        self.pixmap = self.scene.addPixmap(QtGui.QPixmap())
        self.pixmap.setScale(1/self.m_factor)
        self.m_numberImages = 0
        self.m_colornames = ["red", "green", "blue"]
        self.m_DistanceMatrix = []
        self.m_lastPath = "."
        self.m_lastPathSave = "."
        
        self.m_colocalizationHeatmap = []
        self.m_showcolocalizationHeatmap = True
        
        self.m_copyImageRed = []
        self.m_copyImageGreen = []
        self.m_showImages = True
        
        self.sigma = 1
        self.threshold = 50
        
        self.allcolors = [(1,0,0),(0,1,0),(0,0,1),(1,1,0),(0,1,1)] 
        self.connectSignals()
        
    def connectSignals(self):
        QtCore.QObject.connect(self.ui.actionColocalization, QtCore.SIGNAL("triggered()"), self.colocalization)
        QtCore.QObject.connect(self.ui.actionToggle_colocalization_heatmap, QtCore.SIGNAL("triggered()"), self.toggleColocalizationHeatmap)
        QtCore.QObject.connect(self.ui.actionToggle_Images, QtCore.SIGNAL("triggered()"), self.toggleImages)
        self.ui.setSigmaSlider.valueChanged.connect(self.setSigmaSliderChanged)
        self.ui.setSigmaSpinBox.valueChanged.connect(self.setSigmaSpinBoxChanged)
        self.ui.thresholdSlider.valueChanged.connect(self.thresholdSliderChanged)
        self.ui.thresholdSpinBox.valueChanged.connect(self.thresholdSpinBoxChanged)
        self.ui.doSmoothingButton.clicked.connect(self.doSmoothing)
        self.ui.thresholdButton.clicked.connect(self.doThreshold)
        self.ui.actionOpen_File.triggered.connect(self.getFilename)
        self.ui.actionSave_File.triggered.connect(self.saveFile)
        self.ui.actionSave_Features.triggered.connect(self.saveFeatures)
        self.ui.actionReset.triggered.connect(self.reset)
        self.ui.actionShow_covariance_matrix.triggered.connect(self.showCovar)
        self.ui.actionCalculate_distance_matrix.triggered.connect(self.calcDistanceMatrix)
        self.ui.actionShow_Heatmatrix.triggered.connect(self.showHeatmatrix)
        self.ui.addFileButton.clicked.connect(self.getFilename)
        self.ui.zoomInButton.clicked.connect(self.zoomInButton_clicked)
        self.ui.zoomOutButton.clicked.connect(self.zoomOutButton_clicked)
        self.ui.findConnectedComponentsButton.clicked.connect(self.findConnectedComponents)
        self.ui.cutOffButton.clicked.connect(self.cutOff)
        QtCore.QObject.connect(self.scene, QtCore.SIGNAL("cursorMoved(float, float, float)"), self.displayStats)
        
        if 0:
            self.addFile('/home/herrmannsdoerfer/Desktop/SampleCells/01-control/confocal/Pos075_S001C0 - AveProj.tif')
            self.addFile('/home/herrmannsdoerfer/Desktop/SampleCells/01-control/confocal/Pos075_S001C1 - AveProj.tif')
            
        
        if 0:
            self.addFile('/home/herrmannsdoerfer/Desktop/SampleCells/01-control/confocal/Pos081_S001C0 - AveProj.tif')
            self.addFile('/home/herrmannsdoerfer/Desktop/SampleCells/01-control/confocal/Pos081_S001C1 - AveProj.tif')
            imgR, imgG = colocalizationDetection.coloc(self.Images[0].npimage,self.Images[1].npimage)
            self.Images[0].npimage =imgR
            self.Images[1].npimage =imgG
            self.recalculateResult()
        
        if 0:
            self.Images.append(ImageData())
            self.Images.append(ImageData())
            imgR=np.zeros((3810,2960,4))
            imgG=np.zeros((3810,2960,4))
            self.m_numberImages=2
            imgR[...,2]=vigra.impex.readImage('/home/herrmannsdoerfer/Desktop/SampleCellsUnedited/01-control/storm/Pos11_2_aligned.tif', index = 0).swapaxes(0,1).view(np.ndarray).astype(np.float32)[:,:,0]
            imgG[...,1]=vigra.impex.readImage('/home/herrmannsdoerfer/Desktop/SampleCellsUnedited/01-control/storm/Pos11_2_aligned.tif', index = 1).swapaxes(0,1).view(np.ndarray).astype(np.float32)[:,:,0]
            self.Images[0].npimage =imgR
            self.Images[1].npimage =imgG
            self.showHeatmatrix()
        
        if 0:
            self.addFile('/home/herrmannsdoerfer/master/workspace/PythonPrograms/Colocalization/ColocBlockRed.txt')
            self.addFile('/home/herrmannsdoerfer/master/workspace/PythonPrograms/Colocalization/ColocBlockGreen.txt')
            self.showHeatmatrix()
            
        if 0:
            self.addFile('/home/herrmannsdoerfer/master/workspace/PythonPrograms/Colocalization/ColocBlockRed.txt')
            self.addFile('/home/herrmannsdoerfer/master/workspace/PythonPrograms/Colocalization/ColocBlockGreen.txt')
            #self.showHeatmatrix()
            imgR, imgG = colocalizationDetection.Colocdetection(self.Images[0].npimage,self.Images[1].npimage)
            self.Images[0].npimage =imgR
            self.Images[1].npimage =imgG
            self.recalculateResult()
              
    def showHeatmatrix(self):
        #imgR,imgG = colocalizationDetection.createHeatmap2(self.Images[0].npimage,self.Images[1].npimage)
        #imgR,imgG = colocalizationDetection.coordinateBasedColocalization(self.Images[0].npimage,self.Images[1].npimage)
        imgR, imgG = colocalizationDetection.Colocdetection(self.Images[0].npimage,self.Images[1].npimage)
        self.Images[0].npimage =imgR
        self.Images[1].npimage =imgG
        self.recalculateResult()
        
    def setSigmaSliderChanged(self, value):
        self.ui.setSigmaSpinBox.setValue(value/10.)
        self.sigma = value/10.
        
    def setSigmaSpinBoxChanged(self, value): 
        self.ui.setSigmaSlider.setValue(value*10)
        self.sigma = value
        
    def thresholdSliderChanged(self, value):
        self.ui.thresholdSpinBox.setValue(value/10.)
        self.threshold = value/10.
        
    def thresholdSpinBoxChanged(self, value): 
        self.ui.thresholdSlider.setValue(value*10)
        self.threshold = value
    
    def colocalization(self):
        imgR = self.Images[0].npimage
        imgG = self.Images[1].npimage
        self.m_colocalizationHeatmap = colocalizationDetection.coloc(imgR,imgG)
        #Histogram = colocalizationDetection.getAngleHistogram(imgR, imgG, self.m_dims)
        pearsonCoeff = colocalizationDetection.calcPearsonCorrelationCoeff(imgR, imgG)
        MR, MG = colocalizationDetection.calcMandersColocalizationCoeffs(imgR, imgG)
        overlap = colocalizationDetection.calcOverlapCoeff(imgR, imgG)
        self.recalculateResult()    
        
    def toggleColocalizationHeatmap(self):
        if self.m_colocalizationHeatmap != []:
            if self.m_showcolocalizationHeatmap:
                self.Images[0].npimage[...,0] = self.m_colocalizationHeatmap[...,0]
            else:
                self.Images[0].npimage[...,0] = np.zeros(self.Images[0].npimage[:,:,2].shape)
            self.m_showcolocalizationHeatmap = not self.m_showcolocalizationHeatmap
            self.recalculateResult()    
    
    def toggleImages(self):
        if self.m_copyImageRed == []:
            self.m_copyImageRed = copy.deepcopy(self.Images[0].npimage)
            self.m_copyImageGreen = copy.deepcopy(self.Images[1].npimage)
            
        if self.m_showImages:
            self.Images[0].npimage = self.m_copyImageRed
            self.Images[1].npimage = self.m_copyImageGreen
        else:
            self.m_copyImageRed = copy.deepcopy(self.Images[0].npimage)
            self.Images[0].npimage[...,2] = np.zeros(self.Images[0].npimage[...,2].shape)
            self.m_copyImageGreen = copy.deepcopy(self.Images[1].npimage)
            self.Images[1].npimage[...,1] = np.zeros(self.Images[0].npimage[...,1].shape)
        
        self.m_showImages = not self.m_showImages
        
        self.recalculateResult()
                
    def doSmoothing(self):
        #if len(self.m_npimages):
            #for i, imag in enumerate(self.m_npimages):
        if self.ui.selectChannelcomboBox.currentIndex() == 0:
            for i in range(self.m_numberImages):
                imag = self.Images[i].npimage
                imag = vigra.RGBImage(np.array(imag,dtype=np.float32))
                self.Images[i].npimage = np.array(vigra.gaussianSmoothing(imag, self.sigma),dtype=np.uint8)
                #self.m_npimages[i] = np.array(vigra.gaussianSmoothing(imag, self.sigma),dtype=np.uint8)
                
        else:
            i = self.ui.selectChannelcomboBox.currentIndex()-1
            imag = self.Images[i].npimagenp.var
            imag = vigra.RGBImage(np.array(imag,dtype=np.float32))
            self.Images[i].npimage = np.array(vigra.gaussianSmoothing(imag, self.sigma),dtype=np.uint8)
    
        self.recalculateResult()
        self.main_window.statusBar().showMessage('Smoothing done')
        
    def doThreshold(self):
        if self.ui.selectChannelcomboBox.currentIndex() == 0:
            for i in range(self.m_numberImages):
                maximum = np.max(self.Images[i].npimage[:,:,:])
                imag = self.Images[i].npimage
                imag[:,:,2-i]=np.where(imag[:,:,2-i]>=maximum*self.threshold/100,255,0)
                self.Images[i].npimage=imag

        else:
            i = self.ui.selectChannelcomboBox.currentIndex() - 1
            imag = self.Images[i].npimage
            imag[:,:,2-i]=np.where(imag[:,:,2-i]>=255*self.threshold/100,255,0)
            self.Images[i].npimage=imag
        
        self.recalculateResult()
        self.main_window.statusBar().showMessage('Thresholding done')
    
    def getFilename(self):
        filename = QtGui.QFileDialog.getOpenFileName(self.main_window,"Open Image", self.m_lastPath, filter="Coordinate Files (*.txt *.tif)");
        self.addFile(filename)
    
    def addFile(self, filename):
        
        temp = filename.split("/")[-1]
        self.m_lastPath = filename.split(temp)[0]
        if filename:
            extension = filename.split(".")[-1]
            self.Images.append(ImageData())
            self.addImage(filename,self.allcolors[self.ui.fileWidget.count()%len(self.allcolors)],extension) 
            self.ui.fileWidget.addItem(filename)
            self.recalculateResult()
                
            
            
    def recalculateResult(self):
        #if len(self.m_npimages) == 0:
        self.ui.selectChannelcomboBox.clear()
        self.ui.selectChannelcomboBox.addItem("all colors")
        if self.m_numberImages == 0:    
            self.pixmap.setPixmap(QtGui.QPixmap())
            return
        
        #height, width = self.m_npimages[0].shape[0], self.m_npimages[0].shape[1]
        height, width = self.Images[0].npimage.shape[0], self.Images[0].npimage.shape[1]
        self.m_qimage = QtGui.QImage(width, height, QtGui.QImage.Format_RGB32)
        self.m_qimage.fill(0)
        painter = QtGui.QPainter(self.m_qimage);
        painter.setCompositionMode(QtGui.QPainter.CompositionMode_Plus);
        #for i, np_img in enumerate(self.m_npimages):
        for i in range(self.m_numberImages):
            self.ui.selectChannelcomboBox.addItem(self.m_colornames[i])
            np_img = self.Images[i].npimage
            qimg_i = QtGui.QImage(np_img,np_img.shape[1],np_img.shape[0],QtGui.QImage.Format_RGB32)  
            painter.drawImage(0, 0, qimg_i)
        self.pixmap.setPixmap(QtGui.QPixmap.fromImage(self.m_qimage))
        self.ui.graphicsView.fitInView(self.pixmap, QtCore.Qt.KeepAspectRatio)

     
    def reset(self):
        self.Images = []
        self.m_numberImages = 0
        self.ui.fileWidget.clear()
        self.recalculateResult()  
                      
    def addImage(self, filename, color=(1,0,0), extension = "txt"):
        '''load a file and add it to the list of QImages'''
        if (extension == "txt"):
            npimg, cc, dims = worker.loadImage(filename, color, self.m_factor)
            self.scene.setCursorRadius(dims[0]/50.)
        if (extension == "tif"):
            im = Image.open(str(filename))
            imarray = np.array(im)
            npimg, cc, dims = coords.Image2coords(imarray, color)
            self.scene.setCursorRadius(dims[0]/50./self.m_factor)
            #self.loadtif(filename, color, self.m_factor)
        self.Images[self.m_numberImages].npimage = npimg
        self.Images[self.m_numberImages].coords = cc
        self.Images[self.m_numberImages].color = color
        self.Images[self.m_numberImages].dims = dims
        print 'image loaded'
        self.m_numberImages = self.m_numberImages + 1
        
     
        
           
    def zoomOutButton_clicked(self):
        self.zoomGraphicsView(1/1.2)

    def zoomInButton_clicked(self):
        self.zoomGraphicsView(1.2)
        
    def zoomGraphicsView(self, factor):
        self.ui.graphicsView.scale(factor,factor)

    def findConnectedComponents(self):

       # for i, imag in enumerate(self.m_npimages):
        for i in range(self.m_numberImages):
            imag = self.Images[i].npimage
            imag = vigra.RGBImage(np.array(imag, dtype = np.float32))
            #self.m_npimages[i] = np.array(vigra.analysis.labelImage(imag, 4), dtype=np.uint8)
            #self.m_labels[i] = np.zeros((imag.shape[0], imag.shape[1]))
            self.Images[i].labels = np.zeros((imag.shape[0], imag.shape[1]))
            #self.m_labels[i] = vigra.analysis.labelImage(imag[:,:,2-i], 4)
            self.Images[i].labels = (vigra.analysis.labelImageWithBackground(imag[:,:,2-i], 8,background_value=0))
            self.Images[i].findClasses(self.m_factor)            
            self.recalculateResult()
        self.main_window.statusBar().showMessage('Connected components found')
                
    def displayStats(self, y, x, radius):   #the changed and y position compensates that the first component is in a matrix "y"
        number = ""
        mean   = "" 
        row1cov= ""
        row2cov= ""
        blobIndex = ""
        for i in range(self.m_numberImages):
            if self.Images[i].cov_matrix != []:
                selectedClass = self.Images[i].labels[int(x*self.m_factor),int(y*self.m_factor)]
                if np.max(self.Images[i].npimage[int(x*self.m_factor),int(y*self.m_factor),:]) == 0:
                    number += "  bg  "
                    mean   += "  bg  "
                    row1cov += " bg bg "
                    row2cov += " bg bg "
                    blobIndex += " bg "
                else:
                    number += " %3.0f" % self.Images[i].numberElementsOfEachClass[int(np.max(selectedClass))]
                    mean   += "(%3.2f, %3.2f)" % (self.Images[i].meanx[selectedClass] ,self.Images[i].meany[selectedClass])
                    blobIndex += " %3i" % selectedClass
                    if self.Images[i].cov_matrix[selectedClass] != []:
                        row1cov += " [%3.2f, %3.2f" %(self.Images[i].cov_matrix[selectedClass][0,0], self.Images[i].cov_matrix[selectedClass][0,1])
                        row2cov += "  %3.2f, %3.2f]" %(self.Images[i].cov_matrix[selectedClass][1,0], self.Images[i].cov_matrix[selectedClass][1,1])
                
        self.ui.sizelabel.setText(number)
        self.ui.meanlabel.setText(mean)
        self.ui.cov1label.setText(row1cov)
        self.ui.cov2label.setText(row2cov)
        self.ui.blobIndexlabel.setText(blobIndex)
        
    def countDetections(cc,x,y,radius):
        for i, lab in enumerate(self.m_labels):
            labelNumber = lab[x,y]
     
    def showCovar(self):
        for i in range(self.m_numberImages):
            self.Images[i].findCovarianceMatrices(self.m_factor) 
        self.main_window.statusBar().showMessage('Covariance matrix computed')      
    
    def calcDistanceMatrix(self):
        xvec = []
        yvec = []
        for i in range(self.m_numberImages):
            xvec = np.hstack([xvec, self.Images[i].meanx[1:]])  #the first entry is a dummy to correlate the detected classes directly with the entries of the vectors
            yvec = np.hstack([yvec, self.Images[i].meany[1:]])  #the mean of the first class (class 1) is mean[1], therefore the mean[0] entry is skipped
        X = np.vstack([xvec, yvec])
        X = np.array(X)
        X = X.T
        print X
        self.m_DistanceMatrix = scsd.squareform(scsd.pdist(X))
        print self.m_DistanceMatrix.shape
        print self.m_DistanceMatrix
        self.main_window.statusBar().showMessage('Distance matrix computed')
      
    def saveFile(self):
        listItem = self.ui.fileWidget.item(0)
        filename = listItem.data(QtCore.Qt.DisplayRole).toString()
        name = QtCore.QFileInfo(filename).fileName()
        newfilename = QtGui.QFileDialog.getSaveFileName(self.main_window, "Save Coordinate List %s" % name, filter="Coordinate List (*.txt)")
        if not newfilename.contains('.'):
            newfilename += ".txt" # add default extension
        
        print newfilename
        fp = open(newfilename, 'wb')
        csvwriter = csv.writer(fp, delimiter=' ', lineterminator='\n')
        csvwriter.writerow([int(self.Images[0].dims[0]), int(self.Images[0].dims[1])])
        for i in range(int(self.Images[0].dims[0]*self.m_factor)):
            for j in range(int(self.Images[0].dims[1]*self.m_factor)):
                if np.max(self.Images[0].npimage[i,j,:]) != 0:
                    temp = [i/float(self.m_factor),j/float(self.m_factor),np.max(self.Images[0].npimage[i,j,:])]
                    csvwriter.writerow(temp)
        fp.close()
        
    def saveFeatures(self):
        newfilename = QtGui.QFileDialog.getSaveFileName(self.main_window, "Save Coordinate" ,self.m_lastPathSave ,filter="Coordinate List (*.txt)")
        temp = newfilename.split("/")[-1]
        self.m_lastPathSave = newfilename.split(temp)[0]
        if not newfilename.contains('.'):
            newfilename += ".txt" # add default extension
        fp = open(newfilename, 'wb')
        csvwriter = csv.writer(fp, delimiter=' ', lineterminator='\n')
        
        for i in range(1,self.Images[0].numberclasses + 1):
            print i
            csvwriter.writerow([self.Images[0].meanx[i], self.Images[0].meany[i], self.Images[0].numberElementsOfEachClass[i], self.Images[0].cov_matrix[i][0,0],self.Images[0].cov_matrix[i][0,1],self.Images[0].cov_matrix[i][1,0],self.Images[0].cov_matrix[i][1,1]])
  
        fp.close()  
        
    def cutOff(self):
        if self.ui.selectChannelcomboBox.currentIndex() == 0:
            for i in range(self.m_numberImages):
                self.Images[i].calculateEverything(self.ui.numberCutOffspinBox.value(), self.m_factor, i, self)
        else:
            i = self.ui.selectChannelcomboBox.currentIndex() - 1
            self.Images[i].calculateEverything(self.ui.numberCutOffspinBox.value())
        self.recalculateResult()
        
if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    m_app = structureAnalysis()
    sys.exit(app.exec_())