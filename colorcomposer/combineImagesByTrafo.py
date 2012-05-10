import numpy as np
from numpy.linalg import solve
from scipy.misc import imsave

def affineMatrix2DFromCorrespondingPoints(s, d):
	n = len(s)
	assert len(s) == len(d) and s.shape[1]==2 and d.shape[1]==2
	if len(s) < 3:
		print "at least three points required"
		return
	
		
	m  = np.zeros((3,3))
	
	s = np.hstack([s,np.ones((n,1))]) # add third column with ones
	d = np.hstack([d,np.ones((n,1))]) # add third column with ones
	rx = np.sum(d[:,0:1]*s,axis=0)
	ry = np.sum(d[:,1:2]*s,axis=0)

	for i in range(n):
		c = s[i]
		m += np.outer(c,c)

	solx = solve(m, rx)
	soly = solve(m, ry)
	row3 = np.array([0,0,1])
	
	print [solx,soly,row3]
	a=np.vstack([solx,soly])
	error=(d[:,0:2]).T-np.dot(a,s.T)
	
	print error
	s[:,0]=s[:,0]-sum(s[:,0])/len(s)
	s[:,1]=s[:,1]-sum(s[:,1])/len(s)
	#dA=-np.dot(error,(np.dot(s,np.linalg.pinv(np.dot(s.T,s)))))
	dA=-np.dot(error,(np.dot(s,np.linalg.pinv(np.dot(s.T,s)))))
	
	print dA
	print np.dot(dA,s.T)
	
	return np.vstack([solx,soly,row3])

# if run standalone
if __name__ == "__main__":
	import matplotlib.pyplot as plt
	import coords
	file1 = 'data/HeLa5_10000z_3_polyL_actinmEOS2-steve.txt'
	file1_landmarks = 'data/HeLa5_10000z_3_polyL_actinmEOS2-steve.txt_beads.txt'
	file2 = 'data/HeLa5_10000z_3_polyL_ER647.txt'
	file2_landmarks = 'data/HeLa5_10000z_3_polyL_ER647.txt_beads.txt'

	thrash, landmarks1 = coords.readfile(file1_landmarks)
	thrash, landmarks2 = coords.readfile(file2_landmarks)
	landmarks1 = landmarks1[:,:2]
	landmarks2 = landmarks2[:,:2]
	trafo = affineMatrix2DFromCorrespondingPoints(landmarks2, landmarks1)
	print "Transformation: \n", trafo

	dims, cc1 = coords.readfile(file1)
	dims, cc2 = coords.readfile(file2)
	cc2_ones = np.hstack([cc2[:,:2], np.ones((len(cc2),1))]) # add column with ones for affine trafo
	cc2_transformed = np.dot(cc2_ones, trafo.T)[:,:2]
	cc2_transformed = np.hstack([cc2_transformed,cc2[:,2:]]) # add intensity information again

	im1 = coords.coords2Image(dims, cc1)
	im2 = coords.coords2Image(dims, cc2_transformed)

	im_color = np.dstack([im1/np.max(im1),im2/np.max(im2),np.zeros(im1.shape)])
	scipy.misc.imsave(file1+"_processed.png", im_color)
	plt.imshow(im_color)
	plt.plot(landmarks1[:,0]*8, landmarks1[:,1]*8, 'wo', alpha=0.4, scalex=False, scaley=False)
	plt.show()
