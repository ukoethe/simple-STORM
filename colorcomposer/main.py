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

# This class represents a bead Graphics within
# the visual editor environment.
class BeadCircle(QtGui.QGraphicsObject):
	def __init__(self, x , y , color, layerindex, r=2.0):
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

class CursorGraphicsScene(QtGui.QGraphicsScene):
	'''Add a cursor (ellipseItem) to the QGraphicsView with appropriate 
	MouseEvent-Handling'''
	def __init__(self):
		QtGui.QGraphicsScene.__init__(self)
		self.cursorRadius = 10.
		pen = QtGui.QPen(QtGui.QColor("white"))
		self.cursor = self.addEllipse(0,0,2.*self.cursorRadius,2.*self.cursorRadius, pen)
		self.cursor.setZValue(10.) # always foreground

	def mouseReleaseEvent(self, event):
		QtGui.QGraphicsScene.mouseReleaseEvent(self, event) # forward to individual items
		self.cursor.setPos(event.scenePos()-QtCore.QPointF(self.cursorRadius,self.cursorRadius)) # move cursor
		self.emit(QtCore.SIGNAL("cursorMoved(float, float, float)"), event.scenePos().x(), event.scenePos().y(), self.cursorRadius)
		
	def setCursorRadius(self, radius):
		self.cursorRadius = radius
		rect = self.cursor.rect()
		rect.setHeight(2.*radius)
		rect.setWidth(2.*radius)
		self.cursor.setRect(rect)

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
		self.m_factor = 8.   	# factor for scaling coordinates to an image
		self.m_qimage = QtGui.QImage()
		self.scene = CursorGraphicsScene()
		self.ui.graphicsView.setScene(self.scene)
		self.pixmap = self.scene.addPixmap(QtGui.QPixmap())
		self.pixmap.setScale(1/self.m_factor)

		self.connectSignals()

	def connectSignals(self):
		QtCore.QObject.connect(self.ui.zoomOutButton, QtCore.SIGNAL("clicked()"), self.zoomOutButton_clicked)
		QtCore.QObject.connect(self.ui.zoomInButton, QtCore.SIGNAL("clicked()"), self.zoomInButton_clicked)
		QtCore.QObject.connect(self.ui.addFileButton, QtCore.SIGNAL("clicked()"), self.addFileButton_clicked)
		QtCore.QObject.connect(self.ui.doTransformationButton, QtCore.SIGNAL("clicked()"), self.doTransformationButton_clicked)
		QtCore.QObject.connect(self.ui.actionAuto_detect_beads, QtCore.SIGNAL("triggered()"), self.autoDetectBeads_triggered)
		QtCore.QObject.connect(self.ui.actionClear_all, QtCore.SIGNAL("triggered()"), self.clearAll_triggered)
		QtCore.QObject.connect(self.ui.actionAbout, QtCore.SIGNAL("triggered()"), self.about_triggered)
		QtCore.QObject.connect(self.ui.actionExport_Composed_image, QtCore.SIGNAL("triggered()"), self.exportComposed_triggered)
		QtCore.QObject.connect(self.ui.actionExport_transformed_coordinates, QtCore.SIGNAL("triggered()"), self.exportTransformedCoordinates_triggered)
		QtCore.QObject.connect(self.ui.cursorRadiusSpinBox, QtCore.SIGNAL("valueChanged(double)"), self.scene.setCursorRadius)
		QtCore.QObject.connect(self.scene, QtCore.SIGNAL("cursorMoved(float, float, float)"), self.displayStats)

	def zoomGraphicsView(self, factor):
		self.ui.graphicsView.scale(factor,factor)

	def addImage(self, filename, color=(1,1,0)):
		'''load a file and add it to the list of QImages'''
		npimg, cc = worker.loadImage(filename, color, self.m_factor)
		self.m_npimages.append(npimg)
		self.m_coords.append(cc)
		self.m_colors.append(color)
		
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

	def autoMarkBeads(self):
		for i in xrange(self.ui.fileWidget.count()):
			listItem = self.ui.fileWidget.item(i)
			filename = listItem.data(QtCore.Qt.DisplayRole).toString()
			beadPositions = beadLocalisation.detectBeadsFromFile(filename)
			for b in beadPositions:
				color = QtGui.QColor.fromRgbF(*self.m_colors[i])
				bead = BeadCircle(b[0],b[1],color,i)
				bead.connect(bead, QtCore.SIGNAL("beadSelected(float,float,int,bool)"), self, QtCore.SIGNAL("setBeadSelected(float,float,int,bool)"))
				bead.connect(self, QtCore.SIGNAL("hideBeads()"), bead.hide)
				self.scene.addItem(bead)
			self.emit(QtCore.SIGNAL("unselectAllBeads()"))

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
		self.transformcontroller.calculateTransform()
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
		version = "0.1 beta"
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
		for i in range(len(self.m_coords)):
			num, ii = worker.countDetections(self.m_coords[i],x,y,radius)
			number += " %3.0f" % num
			intensity += " %3.0f" % ii
		self.ui.statsIntensity.setText(intensity)
		self.ui.statsNumber.setText(number)

if __name__ == "__main__":
	app = QtGui.QApplication(sys.argv)
	m_app = MyColorcomposerApp()
	sys.exit(app.exec_())
