import numpy as np
from PyQt4 import QtGui, QtCore
import combineImagesByTrafo as calcTrafo
import coords
import copy

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
        self.m_errorMatrix = []
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
            print "bead added"
        else:
            self.m_landmarksList.remove(newlandmark)
            print "bead removed"
        return

    def unselectAllBeads(self):
        self.m_landmarksList = []
        print "all bead unselected"

    def calculateTransform(self,dims):
        '''Calculates the transformation matrix based on the active beads
        and updates the color image accordingly'''
        print "doTransform executed."
        print self.m_landmarksList
        bg = self.getCoordinates(layer=0)
        layer = 0
        heatMatrix = [np.zeros([dims[1],dims[0]])]
        while True: # for all layers
            layer += 1
            fg = self.getCoordinates(layer)
            if len(fg) == 0:
                break # probably last layer reached
            #fg,bg = calcTrafo.distanceMatrix(fg,bg,dims)
            #if len(fg) != len(bg):
            #	print "Different number of foreground points in layer %i and background points in layer 0" % layer
            #	continue

            print "transforming layer %i using %i beads as landmarks" % (layer, len(fg))

            fgc = copy.deepcopy(fg)
            bgc = copy.deepcopy(bg)
            self.m_transform[layer], fg, bg = calcTrafo.doRansac2(fgc, bgc, dims)
            #self.m_transform[layer] = calcTrafo.affineMatrix2DFromCorrespondingPoints(fg, bg, dims)
            tt = np.dot(self.m_transform[layer], np.vstack([fg.T,np.ones(len(fg))])).T[:,:2] # fg transformed in bg coords
            rss = np.sum((tt-bg)**2)
            rms = np.sqrt(rss/len(fg))
            print "Residual sum of squares (RSS): %fpx^2,     RMS: %fpx" % (rss, rms)

            self.m_errorMatrix.append(self.getHeatmatrix(dims, fg, bg,rms))
            bg = bgc
        return self.m_errorMatrix

    def getHeatmatrix(self, dims, fg, bg,rms):


        dims=[int(dims[1]),int(dims[0])]	#now the first dimension is assosiated with right
        deltaX=np.mean(bg[:,0])
        deltaY=np.mean(bg[:,1])

        n = len(fg)

        fg[:,0]=fg[:,0]-deltaX
        fg[:,1]=fg[:,1]-deltaY
        bg[:,0]=bg[:,0]-deltaX
        bg[:,1]=bg[:,1]-deltaY

        #paramX,paramY = calcTrafo.calculateRegression(fg, bg)
        Heatmatrix = np.zeros((dims[1],dims[0]))
        contributionX = np.zeros([dims[0],dims[1]])
        contributionY = np.zeros([dims[0],dims[1]])
        student_coeff95 = np.ones(2000) * 2
        student_coeff95[:19] = [12.7, 4.3, 3.18, 2.77, 2.57, 2.447, 2.365, 2.30, 2.25, 2.28,2.22,2.18,2.16,2.145,2.13,2.12,2.11,2.10,2.09]

        constantX = student_coeff95[n-2]*rms/np.sqrt(2)
        constantY = student_coeff95[n-2]*rms/np.sqrt(2)
        #The error for x and y direction is calculated independently, (confidence interval)
        '''
        for i in range(int(dims[0])):
            contributionX[i] = constantX * np.sqrt(1/len(bg)+((i-deltaX)-np.mean(bg[:,0]))**2/np.sum((bg[:,0]-np.mean(bg[:,0]))**2))

        for j in range(int(dims[1])):
            contributionY[j] = constantY * np.sqrt(1/len(bg)+((j-deltaY)-np.mean(bg[:,1]))**2/np.sum((bg[:,1]-np.mean(bg[:,1]))**2))
        '''
        MX_inv = (np.matrix(fg).T * np.matrix(fg)).I

        for i in range(int(dims[0])):
            for j in range(int(dims[1])):
                x0 = np.matrix([i-deltaX, j- deltaY]).T
                #contributionX[i,j] = constantX**2 * (x0.T * MX_inv * x0)
                contributionX[i,j] = rms**2 * (x0.T * MX_inv * x0)
                #contributionY[i,j] = constantY**2 * (x0.T * MX_inv * x0)
                contributionY[i,j] = rms**2 * (x0.T * MX_inv * x0)

        '''
        contribMatX=np.ones((dims[1],dims[0]))*contributionX
        contribMatY=np.tile(contributionY,(int(dims[0]),1)).T
        Heatmatrix[:,:,0] = np.sqrt(contribMatX+contribMatY)
        '''
        #contribMatX=np.ones((dims[1],dims[0]))*contributionX**2
        #contribMatY=np.tile(contributionY**2,(int(dims[0]),1)).T
        Heatmatrix = np.sqrt(contributionX+contributionY).T

        from matplotlib import pyplot
        pyplot.matshow(Heatmatrix)
        pyplot.colorbar()
        pyplot.gray()
        pyplot.show()

        '''deltaX = dims[0]/2
        deltaY = dims[1]/2

        Heatmatrix = np.ones([dims[0], dims[1],4])

        for i in range(int(dims[0])):
            for j in range(int(dims[1])):
                Heatmatrix[i,j,0] = delR + np.abs(np.sqrt((i-deltaX)**2+(j-deltaY)**2) * delPhi/(2*np.pi))


        from matplotlib import pyplot
        pyplot.matshow(Heatmatrix[...,0])
        pyplot.colorbar()
        pyplot.gray()
        pyplot.show()
        Heatmatrix=Heatmatrix.T'''
        return Heatmatrix

    def getCoordinates(self, layer):
        cc = []
        for l in self.m_landmarksList:
            if l.active and l.layerindex == layer:
                cc.append(l.pos)
        return np.asarray(cc)

    def getTransform(self, i):
        '''return the transformation matrix for homogeneous coordinates that
        transforms layer i into coordinates of layer 0'''
        if self.m_transform.has_key(i):
            trafo = self.m_transform[i]
        else:
            trafo = np.diag([1,1,1]) # identity
        # convert to QMatrix
        matrix = QtGui.QMatrix(trafo[0,0], trafo[1,0], trafo[0,1], trafo[1,1], trafo[0,2], trafo[1,2])
        return matrix, trafo

    def doTransform(self, points, layer):
        '''transform the points with the previously calculated transform
        from the coordinates at layer l to the coordinate system of the
        background layer (layer 0)'''

        return np.dot(self.m_transform[layer], np.vstack([points.T,np.ones(len(points))])).T[:,:2] # fg transformed in bg coords
        
    def addTransformationError(self, cc, factor, layer):
        for i in range(cc.shape[0]):
            #print i, cc.shape[0], self.m_errorMatrix[layer].shape
            cc[i,5] += self.m_errorMatrix[layer][cc[i,1],cc[i,0]]
