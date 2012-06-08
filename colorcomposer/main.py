 #!/usr/bin/env python

import PyQt4
import sys
from ui_mainwindow import *
from PyQt4 import QtGui, QtCore
import worker
import coords
import numpy as np
import beadLocalisationVariance as beadLocalisation
import transformcontroller as tc
import colocalizationDetection

# This class represents a bead Graphics within
# the visual editor environment.
class BeadCircle(QtGui.QGraphicsObject):
	def __init__(self, x , y , color, layerindex, r=5.0):
		QtGui.QGraphicsObject.__init__(self)
		self.setFlags(QtGui.QGraphicsItem.ItemIsSelectable)
		#~ self.setRect( -r , -r , 2.*r, 2.*r )
		self.setPos(x,y)
		self.radius = r
		self.boundingRect = QtCore.QRectF()
		self.boundingRect.setCoords( -r, -r, r, r)
		self.color = QtGui.QColor(color)
		self.color.setAlphaF(0.75)
		self.selectedcolor = QtGui.QColor(color)

		self.defaultpen = QtGui.QPen(self.color)
		self.defaultpen.setWidth(1.0)
		self.markpen = QtGui.QPen(self.selectedcolor)
		self.markpen.setWidth(1.0)

		self.active = 0
		self.layerindex = layerindex
		self.setZValue(layerindex)

	def paint(self,painter,option,pwidget):
		if self.active:
			painter.setPen( self.markpen )
			painter.drawEllipse(QtCore.QPointF(0.,0.), self.radius/2., self.radius/2.)
		else:
			painter.setPen( self.defaultpen )
			painter.drawEllipse( self.boundingRect )

	def boundingRect(self):
		return self.boundingRect

	def mouseDoubleClickEvent( self, event):
		self.active = not self.active
		self.emit(QtCore.SIGNAL("beadSelected(float,float,int,bool)"), \
			self.x(),self.y(),self.layerindex,self.active)
		self.update() #redraw
		
	def BigBeadCircle(self):
		print 'updated'
		self.update()

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
		
	def setCursorRadius(self, radius):
		self.cursorRadius = radius
		rect = self.cursor.rect()
		rect.setHeight(2.*radius)
		rect.setWidth(2.*radius)
		self.cursor.setRect(rect)
	
	def getCursorPosition(self):
		return self.lastCursorPosition
		
	def cursorRadius(self):
		return self.cursorRadius

class MyColorcomposerApp(QtCore.QObject):
	def __init__(self):
		QtCore.QObject.__init__(self)
		self.main_window = QtGui.QMainWindow()
		self.ui = Ui_MainWindow()
		self.ui.setupUi(self.main_window)
		self.main_window.show()
		self.transformcontroller = tc.TransformController(self)
		
		self.m_npimages = []  	# numpy-images
		self.m_coords = []      # coordinate lists
		self.m_colors = []    	# list of colors for the images
		self.m_factor = 4.   	# factor for scaling coordinates to an image
		self.m_dims   = []
		self.m_qimage = QtGui.QImage()
		self.scene = CursorGraphicsScene()
		self.ui.graphicsView.setScene(self.scene)
		self.pixmap = self.scene.addPixmap(QtGui.QPixmap())
		self.pixmap.setScale(1/self.m_factor)

		self.connectSignals()
		
		self.cutoff=.95
		self.maxVariance=1
		self.maxDistance=3

	def connectSignals(self):
		QtCore.QObject.connect(self.ui.zoomOutButton, QtCore.SIGNAL("clicked()"), self.zoomOutButton_clicked)
		QtCore.QObject.connect(self.ui.zoomInButton, QtCore.SIGNAL("clicked()"), self.zoomInButton_clicked)
		QtCore.QObject.connect(self.ui.addFileButton, QtCore.SIGNAL("clicked()"	), self.addFileButton_clicked)
		QtCore.QObject.connect(self.ui.doTransformationButton, QtCore.SIGNAL("clicked()"), self.doTransformationButton_clicked)
		QtCore.QObject.connect(self.ui.actionAuto_detect_beads, QtCore.SIGNAL("triggered()"), self.autoDetectBeads_triggered)
		QtCore.QObject.connect(self.ui.actionClear_all, QtCore.SIGNAL("triggered()"), self.clearAll_triggered)
		QtCore.QObject.connect(self.ui.actionColocalization, QtCore.SIGNAL("triggered()"), self.Colocalization)
		QtCore.QObject.connect(self.ui.actionAbout, QtCore.SIGNAL("triggered()"), self.about_triggered)
		QtCore.QObject.connect(self.ui.actionExport_Composed_image, QtCore.SIGNAL("triggered()"), self.exportComposed_triggered)
		QtCore.QObject.connect(self.ui.actionExport_transformed_coordinates, QtCore.SIGNAL("triggered()"), self.exportTransformedCoordinates_triggered)
		#QtCore.QObject.connect(self.ui.cursorRadiusSpinBox, QtCore.SIGNAL("valueChanged(double)"), self.scene.setCursorRadius)
		#QtCore.QObject.connect(self.ui.cursorRadiusSlider, QtCore.SIGNAL("valueChanged(int)"), self.SliderChanged)
		QtCore.QObject.connect(self.scene, QtCore.SIGNAL("cursorMoved(float, float, float)"), self.displayStats)
		
		self.ui.deleteBeadButton.clicked.connect(self.deleteBead)
		self.ui.addRedBeadButton.clicked.connect(self.addRedBead)
		self.ui.addGreenBeadButton.clicked.connect(self.addGreenBead)
		
		self.ui.cursorRadiusSlider.valueChanged.connect(self.RadiusSliderChanged)
		self.ui.cursorRadiusSpinBox.valueChanged.connect(self.RadiusSpinBoxChanged)
		self.ui.PercentOfOccurrencySlider.valueChanged.connect(self.PercentOfOccurrencySliderChanged)
		self.ui.PercentOfOccurrencySpinBox.valueChanged.connect(self.PercentOfOccurrencySpinBoxChanged)
		self.ui.maximalVarianceSlider.valueChanged.connect(self.maximalVarianceSliderChanged)
		self.ui.maximalVarianceSpinBox.valueChanged.connect(self.maximalVarianceSpinBoxChanged)
		self.ui.maximalDistanceSlider.valueChanged.connect(self.maximalDistanceSliderChanged)
		self.ui.maximalDistanceSpinBox.valueChanged.connect(self.maximalDistanceSpinBoxChanged)

		
	def zoomGraphicsView(self, factor):
		self.ui.graphicsView.scale(factor,factor)

	def addImage(self, filename, color=(1,1,0)):
		'''load a file and add it to the list of QImages'''
		npimg, cc, dims = worker.loadImage(filename, color, self.m_factor)
		self.m_npimages.append(npimg)
		self.m_coords.append(cc)
		self.m_colors.append(color)
		self.m_dims.append(dims)
		if len(self.m_colors)==2:
			imgR = self.m_npimages[0]
			imgG = self.m_npimages[1]
			#Histogram = colocalizationDetection.getAngleHistogram(self.m_coords[0], self.m_coords[1], self.m_dims)
			pearsonCoeff = colocalizationDetection.calcPearsonCorrelationCoeff(imgR, imgG)
			MR, MG = colocalizationDetection.calcMandersColocalizationCoeffs(imgR, imgG)
			overlap = colocalizationDetection.calcOverlapCoeff(imgR, imgG)
		
			message='Pearson = %.3f '%pearsonCoeff+' Manders Green = %.3f '%MG+' Manders Red = %.3f '%MR+' overalp coef= %.3f '%overlap
			self.main_window.statusBar().showMessage(message)
		
		
	def recalculateResult(self):
		if len(self.m_npimages) == 0:
			self.pixmap.setPixmap(QtGui.QPixmap())
			return
		height, width = self.m_npimages[0].shape[0], self.m_npimages[0].shape[1]
		self.m_qimage = QtGui.QImage(width, height, QtGui.QImage.Format_RGB32)
		self.m_qimage.fill(0)
		painter = QtGui.QPainter(self.m_qimage);
		painter.setCompositionMode(QtGui.QPainter.CompositionMode_Plus);
		for i, np_img in enumerate(self.m_npimages):
			trafo = self.transformcontroller.getTransform(i)
			f = self.m_factor
			qimg_i = QtGui.QImage(np_img,np_img.shape[1],np_img.shape[0],QtGui.QImage.Format_RGB32)
			qimg_i = qimg_i.transformed(trafo, QtCore.Qt.SmoothTransformation)
			trueMatrix = QtGui.QImage.trueMatrix(trafo, np_img.shape[1], np_img.shape[0]) # Qt does an additional shift, that we have to compensate
			painter.drawImage(f*trafo.dx()-trueMatrix.dx(),f*trafo.dy()-trueMatrix.dy(), qimg_i)
		self.pixmap.setPixmap(QtGui.QPixmap.fromImage(self.m_qimage))
		self.ui.graphicsView.fitInView(self.pixmap, QtCore.Qt.KeepAspectRatio)
	
	def Colocalization(self):
		imgR = self.m_npimages[0]
		imgG = self.m_npimages[1]
		#Histogram = colocalizationDetection.getAngleHistogram(self.m_coords[0], self.m_coords[1], self.m_dims)
		pearsonCoeff = colocalizationDetection.calcPearsonCorrelationCoeff(imgR, imgG)
		print 'Pearson Coefficient: ',pearsonCoeff
		overlap = colocalizationDetection.calcOverlapCoeff(imgR, imgG)
		print 'Overlap coefficient: ',overlap
		MR, MG = colocalizationDetection.calcMandersColocalizationCoeffs(imgR, imgG)
		print MR, MG
		
		Histogram = colocalizationDetection.getAngleHistogram(imgR, imgG, self.m_dims)
		
		
	def autoMarkBeads(self):
		for i in xrange(self.ui.fileWidget.count()):
			listItem = self.ui.fileWidget.item(i)
			filename = listItem.data(QtCore.Qt.DisplayRole).toString()
			beadPositions = beadLocalisation.detectBeadsFromFile(filename,self.cutoff,self.maxDistance,self.maxVariance)	
			for b in beadPositions:
				color = QtGui.QColor.fromRgbF(*self.m_colors[i])
				bead = BeadCircle(b[0],b[1],color,i)
				bead.connect(bead, QtCore.SIGNAL("beadSelected(float,float,int,bool)"), self, QtCore.SIGNAL("setBeadSelected(float,float,int,bool)"))
				bead.connect(self, QtCore.SIGNAL("hideBeads()"), bead.hide)
				self.scene.addItem(bead)
			self.emit(QtCore.SIGNAL("unselectAllBeads()"))

	def getCircledBead(self, layer):						#returns the bead (layer=0 for red, layer=1 for green) in the rectangular selected by single mouseclick
		listItem = self.ui.fileWidget.item(layer)
		filename = listItem.data(QtCore.Qt.DisplayRole).toString()
		dims, cc = coords.readfile(filename)
		rect=CursorGraphicsScene.getCursorPosition(self.scene)
		maxDist=CursorGraphicsScene.cursorRadius(self.scene)
		roi=[rect[0]-maxDist, rect[0]+maxDist, rect[1]-maxDist, rect[1]+maxDist]
		candidate = coords.cropROI(cc, roi)[:,:4]

		if len(candidate) < self.cutoff*dims[2]:
			print 'Number of points below percentage of occurrency 2'
			bead='false'
		else:
			mm,stddev,intensitiy = beadLocalisation.beadVariance(candidate)		
			color = QtGui.QColor.fromRgbF(*self.m_colors[layer])
			bead = BeadCircle(mm[0],mm[1],color,layer)
		if len(candidate) > dims[2]: 
			print ' WARNING: Number of points below percentage of occurrency because len(candidate) =%d > dims[2] = %d'%(len(candidate),dims[2])
			#bead='false'
		return bead
	
	def saveImage(self, filename):
		if not self.m_qimage.save(filename):
			QtGui.QMessageBox.critical(self.main_window, "Error saving file", "Could not save image as %s" %filename)
		return

	########################################
	# SIGNAL handlers                      #
	########################################
	def zoomOutButton_clicked(self):
		self.zoomGraphicsView(1/1.2)

	def zoomInButton_clicked(self):
		self.zoomGraphicsView(1.2)

	def addFileButton_clicked(self):
		filename = QtGui.QFileDialog.getOpenFileName(self.main_window,
			"Open Image", filter="Coordinate Files (*.txt)");
		colors=[(1,0,0),(0,1,0),(0,0,1),(1,1,0),(0,1,1)] # FIXME
		if filename:
			self.addImage(filename,colors[self.ui.fileWidget.count()%len(colors)])
			self.ui.fileWidget.addItem(filename)
			self.recalculateResult()

	def doTransformationButton_clicked(self):
		self.transformcontroller.calculateTransform([self.m_npimages[0].shape[0]/self.m_factor, self.m_npimages[0].shape[1]/self.m_factor])
		self.emit(QtCore.SIGNAL("hideBeads()"))
		self.recalculateResult()

	def autoDetectBeads_triggered(self):
		print "auto-detecting beads..."
		self.autoMarkBeads()

	def clearAll_triggered(self):
		self.emit(QtCore.SIGNAL("hideBeads()")) # we don't need the landmarks any longer
		self.m_npimages = []
		self.m_colors = []
		self.m_coords = []
		self.ui.fileWidget.clear()
		self.recalculateResult()
		
	
	def about_triggered(self):
		version = "0.1 betacoord2im"
		copyright = "(c) Joachim Schleicher <J.Schleicher@stud.uni-heidelberg.de>"
		license = "This is free software; There is NO warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE"
		QtGui.QMessageBox.about(self.main_window, "About Colorcomposer", 
		"ColorComposer Version: " + version + "\n" + copyright + "\n" + license)

	def exportComposed_triggered(self):
		filename = QtGui.QFileDialog.getSaveFileName(self.main_window, "Save Image", filter="Image (*.png *.jpg *.bmp)")
		if not filename.contains('.'):
			filename += ".png" # add default extension
		self.saveImage(filename)

	def exportTransformedCoordinates_triggered(self):
		'''export coordinate file list of tranformed points for all layers'''
		for i in xrange(1,self.ui.fileWidget.count()):
			listItem = self.ui.fileWidget.item(i)
			filename = listItem.data(QtCore.Qt.DisplayRole).toString()
			name = QtCore.QFileInfo(filename).fileName() # without path
			newfilename = QtGui.QFileDialog.getSaveFileName(self.main_window, "Save Coordinate List %s" % name, filter="Coordinate List (*.txt)")
			if not newfilename.contains('.'):
				newfilename += ".txt" # add default extension
			dimensions, points = coords.readfile(filename)
			p_transformed = self.transformcontroller.doTransform(points[:,0:2], i)
			points[:,0:2] = p_transformed
			coords.writefile(newfilename, dimensions[0:2], points)

	def displayStats(self, x, y, radius):
		number = "number:   "
		intensity = "intensity: "
		numberFrame = ""
		for i in range(len(self.m_coords)):
			num, ii = worker.countDetections(self.m_coords[i],x,y,radius)
			number += " %3.0f" % num
			intensity += " %3.0f" % ii
			numberFrame += "  %i" % self.m_dims[i][2]
		self.ui.statsIntensity.setText(intensity)
		self.ui.statsNumber.setText(number)
		self.ui.statsNumberOfFrames.setText(numberFrame)
		
	def deleteBead(self):			
		deltax=1.5
		deltay=1.5
		for obj in self.scene.items():		
			lay=int(obj.zValue())				#lay is the value of the zValue 0=red, 1=green
			try:								#if the x and y means of the object are the same as the means of the new created bead, that is returned from getCirceldBead, obj will be removed
				#if obj.pos().x()==self.getCircledBead(lay).pos().x() and obj.pos().y()==self.getCircledBead(lay).pos().y():
				if self.getCircledBead(lay).pos().x()-deltax<obj.pos().x()<self.getCircledBead(lay).pos().x()+deltax and self.getCircledBead(lay).pos().y()-deltay<obj.pos().y()<self.getCircledBead(lay).pos().y()+deltay:
					self.scene.removeItem(obj)
			except AttributeError:
				print ' '
	
	
	def addRedBead(self):
		bead=self.getCircledBead(0) #returns a BeadCircle Object or the string 'false' if no bead could be found at the given location
		if bead!='false':
			bead.connect(bead, QtCore.SIGNAL("beadSelected(float,float,int,bool)"), self, QtCore.SIGNAL("setBeadSelected(float,float,int,bool)"))
			bead.connect(self, QtCore.SIGNAL("hideBeads()"), bead.hide)
			self.scene.addItem(bead)

	def addGreenBead(self):		
		bead=self.getCircledBead(1) #returns a BeadCircle Object or the string 'false' if no bead could be found at the given location
		if bead!='false':
			bead.connect(bead, QtCore.SIGNAL("beadSelected(float,float,int,bool)"), self, QtCore.SIGNAL("setBeadSelected(float,float,int,bool)"))
			bead.connect(self, QtCore.SIGNAL("hideBeads()"), bead.hide)
			self.scene.addItem(bead)
	
	def RadiusSliderChanged(self,value):
		self.scene.setCursorRadius(value/10)
		self.ui.cursorRadiusSpinBox.setValue(value/10)
		
		
	def RadiusSpinBoxChanged(self,value):
		self.ui.cursorRadiusSlider.setValue(value*10)
		self.scene.setCursorRadius(value)
		
	def PercentOfOccurrencySliderChanged(self,value):
		self.ui.PercentOfOccurrencySpinBox.setValue(value/10.)	
		self.cutoff=value/10/100.	#/10 is the correction needed because of the slider /100 is because cutoff needs to be between 0 and 1
	
		
	def PercentOfOccurrencySpinBoxChanged(self,value):
		self.ui.PercentOfOccurrencySlider.setValue(value*10)
		self.cutoff=value/100.		#cutoff needs to be between 0 and 1
		
	def maximalVarianceSliderChanged(self,value):
		self.ui.maximalVarianceSpinBox.setValue(value/10.)	
		self.maxVariance=value/10.	
	
		
	def maximalVarianceSpinBoxChanged(self,value):
		self.ui.maximalVarianceSlider.setValue(value*10)
		self.maxVariance=value		
		
	def maximalDistanceSliderChanged(self,value):
		self.ui.maximalDistanceSpinBox.setValue(value/10.)	
		self.maxDistance=value/10.	
	
		
	def maximalDistanceSpinBoxChanged(self,value):
		self.ui.maximalDistanceSlider.setValue(value*10)
		self.maxDistance=value		
	
	
		
if __name__ == "__main__":
	app = QtGui.QApplication(sys.argv)
	m_app = MyColorcomposerApp()
	sys.exit(app.exec_())
