import coords
import numpy as np


def loadImage(filename, color=(1.,1.,0.), factor=4):
	print "reading file %s" % filename
	dims, cc = coords.readfile(filename)
	img = coords.coords2Image(dims, cc, factor=factor)
	colorimg = np.zeros((img.shape[0],img.shape[1],4),dtype=np.uint8)
	mx = np.max(img)
	colorimg[:,:,0] = color[2]*img*(255./mx) # blue = LSB
	colorimg[:,:,1] = color[1]*img*(255./mx) # green
	colorimg[:,:,2] = color[0]*img*(255./mx) # red

	return colorimg, cc

def countDetections(cc,x,y,radius):
	roi = x-radius, x+radius, y-radius, y+radius
	cc = coords.cropROI(cc, roi)
	number = 0
	intensity = 0.
	center = np.array([x,y])
	cc[:,4] = np.sqrt(np.sum((cc[:,0:2]-center)**2,axis=1))
	idxs = cc[:,4]<radius
	if(idxs!=[]):
		cc = cc[idxs]
		return len(cc), np.sum(cc[:,3])
	else: # empty set
		return 0, 0.
