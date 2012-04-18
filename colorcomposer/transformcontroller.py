import numpy as np
from PyQt4 import QtGui, QtCore
import combineImagesByTrafo as calcTrafo

class Landmark:
	'''model class that holds the data of a single bead'''
	def __init__(self, pos, layerindex):
		self.active = 1
		self.pos = np.asarray(pos) #center
		self.layerindex = layerindex #belongs to i-th image

	def __eq__(self, other):
		return (self.pos==other.pos).all() and self.layerindex==other.layerindex

	def __str__(self):
		return "(%.2f, %.2f) %i" %(self.pos[0], self.pos[1], self.layerindex)

	def __repr__(self):
		return "Landmark( [%f, %f], %i )" % (self.pos[0], self.pos[1], self.layerindex)

class TransformController(QtCore.QObject):
	def __init__(self, view):
		QtCore.QObject.__init__(self)
		self.m_landmarksList = []
		self.m_transformationList = []
		self.m_imageList = []
		self.m_combinedImage = 0
		self.m_transform = {}
		self.mainview = view
		self.connectSignals()
		
	def connectSignals(self):
		QtCore.QObject.connect(self.mainview, QtCore.SIGNAL("setBeadSelected(float,float,int,bool)"), self.setBeadSelected)
		QtCore.QObject.connect(self.mainview, QtCore.SIGNAL("unselectAllBeads()"), self.unselectAllBeads)
		

	def getResultImage(self):
		return self.combinedImage

	def setBeadSelected(self, x,y,layerindex,active=True):
		'''Called by appropriate signal'''
		newlandmark = Landmark([x,y],layerindex)
		if active:
			self.m_landmarksList.append(newlandmark)
		else:
			self.m_landmarksList.remove(newlandmark)
		return

	def unselectAllBeads(self):
		self.m_landmarksList = []

	def calculateTransform(self):
		'''Calculates the transformation matrix based on the active beads
		and updates the color image accordingly'''
		print "doTransform executed."
		print self.m_landmarksList
		bg = self.getCoordinates(layer=0)
		layer = 0
		while True: # for all layers
			layer += 1
			fg = self.getCoordinates(layer)
			if len(fg) == 0:
				break # probably last layer reached
			if len(fg) != len(bg):
				print "Different number of foreground points in layer %i and background points in layer 0" % layer
				continue
			print "transforming layer %i using %i beads as landmarks" % (layer, len(fg))
			self.m_transform[layer] = calcTrafo.affineMatrix2DFromCorrespondingPoints(fg, bg)
			tt = np.dot(self.m_transform[layer], np.vstack([fg.T,np.ones(len(fg))])).T[:,:2] # fg transformed in bg coords
			rss = np.sum((tt-bg)**2)
			rms = np.sqrt(rss/len(fg))
			print "Residual sum of squares (RSS): %fpx^2,     RMS: %fpx" % (rss, rms)
		return

	def getCoordinates(self, layer):
		cc = []
		for l in self.m_landmarksList:
			if l.active and l.layerindex == layer:
				cc.append(l.pos)
		return np.asarray(cc)

	def getTransform(self, i):
		'''return the transformation matrix for homogenous coordinates that
		transforms layer i into coordinates of layer 0'''
		if self.m_transform.has_key(i):
			trafo = self.m_transform[i]
		else:
			trafo = np.diag([1,1,1]) # identity
		# convert to QMatrix
		matrix = QtGui.QMatrix(trafo[0,0], trafo[1,0], trafo[0,1], trafo[1,1], trafo[0,2], trafo[1,2])
		return matrix

	def doTransform(self, points, layer):
		'''transform the points with the previously calculated transform
		from the coordinates at layer l to the coordinate system of the
		background layer (layer 0)'''
		return np.dot(self.m_transform[layer], np.vstack([points.T,np.ones(len(points))])).T[:,:2] # fg transformed in bg coords

