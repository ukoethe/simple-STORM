import coords
import numpy as np
import os.path
from scipy import signal

def loadImage(filename, color=(1.,1.,0.), factor=4):
	print "reading file %s" % filename
	
	dims, cc = coords.readfile(filename)
	img = coords.coords2Image(dims, cc, factor=factor)
	colorimg = np.zeros((img.shape[0],img.shape[1],4),dtype=np.uint8)
	mx = np.max(img)
	colorimg[:,:,0] = color[2]*img*(255./mx) # blue = LSB
	colorimg[:,:,1] = color[1]*img*(255./mx) # green
	colorimg[:,:,2] = color[0]*img*(255./mx) # red

	return colorimg, cc, dims

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

def smooth_image_according_to_heatmatrix(img, heatmatrix, factor):
	filter_width = 2
	heatmatrix = np.array(heatmatrix[1])
	dim0, dim1, dim2 = img.shape
	res_img = np.zeros_like(img)
	counter = 0
	for i in range(filter_width, int(dim0 - filter_width)):
		#print 'smooth', i
		for j in range(filter_width, int(dim1 - filter_width)):
			if img[i,j,1] != 0:
				counter = counter + 1
				Gauss1d = np.matrix(signal.get_window(('gaussian',heatmatrix[i/factor,j/factor]/3*factor),2*filter_width+1))
				Gauss2d = Gauss1d.T * Gauss1d/(2*np.pi*(heatmatrix[i/factor,j/factor]/3*factor)**2)
				#print np.sum(Gauss2d)
				res_img[i-filter_width:i+filter_width+1, j-filter_width:j+filter_width+1,1] += Gauss2d*img[i,j,1]
				'''for k0 in range(-filter_width, filter_width+1):
					for k1 in range(-filter_width, filter_width+1):
						res_img[i+k0,j+k1,1] = res_img[i+k0,j+k1,1] + Gauss2d[k0 + filter_width, k1 + filter_width] * img[i,j,1]'''
	print np.mean(res_img[...,1])
	print 'counter alt', counter
	return res_img

def smooth_image_according_to_heatmatrix_new(img, heatmatrix, factor):
	filter_width = 4
	heatmatrix = np.array(heatmatrix[1])
	heatmatrix[np.where(heatmatrix == 0)] = 0.001 # this doesn't change much but prevents gaussian from becoming singular
	dim0, dim1, dim2 = img.shape
	res_img = np.zeros_like(img)
	pos0, pos1 = np.where(img[filter_width:-filter_width,filter_width:-filter_width,1] > 0)
	pos0 = pos0 + filter_width #indices have to be corrected, because they are searched in a smaller array (with borders of width filter_width)
	pos1 = pos1 + filter_width 
	for k in range(len(pos0)):
		i = pos0[k]
		j = pos1[k]
		Gauss1d = np.matrix(signal.get_window(('gaussian',heatmatrix[i/factor,j/factor]/3*factor),2*filter_width+1))
		Gauss2d = Gauss1d.T * Gauss1d/(np.sum(Gauss1d)**2)
		#print i,j, dim0,dim1
		res_img[i-filter_width:i+filter_width+1, j-filter_width:j+filter_width+1,1] += Gauss2d*img[i,j,1]

	return res_img