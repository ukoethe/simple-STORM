import numpy as np
from numpy.linalg import solve
from scipy.misc import imsave

def affineMatrix2DFromCorrespondingPoints(s, d, dims):
	n = len(s)
	assert len(s) == len(d) and s.shape[1]==2 and d.shape[1]==2
	if len(s) < 3:
		print "at least three points required"
		return
	
	s[:,0]=s[:,0]-100#-dims[0]
	s[:,1]=s[:,1]-100#-dims[1]
	d[:,0]=d[:,0]-100#-dims[0]
	d[:,1]=d[:,1]-100#-dims[1]
	
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
	
	A = np.vstack([solx,soly,row3])
	R=np.matrix([[1, 0, 100],[0, 1, 100],[0, 0, 1]])
	Ap = np.dot(R,np.dot(A,R.I))
	
	print "matrix:",[solx,soly,row3]
	a=np.vstack([solx,soly])
	error=(d[:,0:2]).T-np.dot(a,s.T)
	#print error,"=",(d[:,0:2]).T,"-",np.dot(a,s.T)
	print "error",error
	s[:,0]=s[:,0]
	s[:,1]=s[:,1]
	
	#dA=-np.dot(error,(np.dot(s,np.linalg.pinv(np.dot(s.T,s)))))
	#print np.linalg.pinv(s).shape, error.shape
	dA=-np.dot(error,np.linalg.pinv(s).T)
	
	#print dA
	print np.dot([solx,soly,row3],s.T)
	
	#return np.vstack([solx,soly,row3])
	return np.array([[Ap[0,0],Ap[0,1],Ap[0,2]],[Ap[1,0],Ap[1,1],Ap[1,2]],[Ap[2,0],Ap[2,1],Ap[2,2]]])
	
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
