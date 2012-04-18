import numpy as np

import coords as coord2im
import sys

def dist(a, b):
	a1 = np.array(a)
	b1 = np.array(b)
	return np.linalg.norm(a1[:2]-b1[:2])

def nearestNeighbor(p, beads, maxDist):
	x,y,intensity = p
	nearestIdx = 0 #initial guess
	nearestDist = 32768 #large number
	for i in range(len(beads)):
		if not x-maxDist <= beads[i][0] <= x+maxDist:
			continue
		if not y-maxDist <= beads[i][1] <= y+maxDist:
			continue
		#check nearest neighbor
		if dist(beads[i], p)<nearestDist:
			nearestIdx = i
			nearestDist = dist(beads[i], p)
	return nearestDist, nearestIdx

def beadVariance(positions, mean=None):
	pos = positions[:,0:2]
	intensities = positions[:,3]
	if mean is None:
		mean=np.mean(pos, 0)
	dists = np.sqrt(np.sum((pos-mean)**2, axis=1))
	stddev = np.sqrt(np.mean(dists**2))
	intensity = np.mean(intensities)
	return mean, stddev, intensity

def detectBeads(dims, cc, maxDist=2, maxStdDev=0.3):
	numFrames = cc[-1,2]+1

	beads = []
	singleConsidered=0
	skipped = 0 #counter
	scatterplotData = []
	meanData = []
	# sort beads in 'bins'
	intensities = cc[:,3]
	#~ allCandidates = cc[intensities > 2*np.mean(intensities)] # only high intensities
	allCandidates = cc # all spots
	initialCandiates = allCandidates[allCandidates[:,2]==0] # from frame 0
	for x, y, frame, intensity in initialCandiates[:,0:4]:
		nearestDist, nearestIdx = nearestNeighbor((x,y,intensity), beads, maxDist)
		if nearestDist > maxDist:
			beads.append([x,y,intensity])
		else:
			print "removing beads too close together: ", (x,y,intensity), beads[nearestIdx]
			beads.remove(beads[nearestIdx])

	for beadCandidate in beads: # over all candidates
		x,y,intensity = beadCandidate
		roi = [x-maxDist, x+maxDist, y-maxDist, y+maxDist]
		candidate = coord2im.cropROI(cc, roi)[:,:4]
		singleConsidered+=len(candidate)
		if len(candidate) > numFrames or len(candidate) < 0.950*numFrames: # allowing few percent not detected
			skipped += 1
			continue
		mm,stddev,intensitiy = beadVariance(candidate)
		if stddev < maxStdDev:
			print "mean: ", mm, "variance of dist: ", stddev, "@intensity: ", intensity, "#", len(candidate)
			scatterplotData.append( (intensity, stddev))
			meanData.append(mm)
		else:
			print "IGNORED (variance too large): mean: ", mm, "variance of dist: ", stddev, "@intensity: ", intensity, "#", len(candidate)
	print singleConsidered, "single detections considered."
	print skipped, "initial candidates skipped."
	return meanData, scatterplotData

def detectBeadsFromFile(filename, maxDist=2, maxStdDev=0.3):
	dims, cc = coord2im.readfile(filename)
	import time
	start = time.time()
	beadMeans, beadScatter = detectBeads(dims, cc, maxDist, maxStdDev)
	print "calculated in %f seconds" % (time.time()-start)
	return beadMeans

if __name__ == "__main__":
	if len(sys.argv != 2):
		print "Usage: %s coordsfile.txt" % sys.argv[0]
		sys.exit(1)
	filename = sys.argv[1]
	dims, cc = coord2im.readfile(filename)
	beadMeans, beadScatter = detectBeads(dims, cc)

	print "number of beads: ", len(beadMeans)
	#~ print singleConsidered, "single detections considered."
	#~ print skipped, "beads skipped."
	print "mean variance: %f" % (np.mean(beadScatter, axis=0)[1])

	#plot:
	import matplotlib.pyplot as plt
	xx,yy=np.hsplit(np.array(beadScatter),np.array([1]))
	plt.plot(xx,yy,'rx')
	plt.xlabel('Beat Intensity')
	plt.ylabel('stddev of detected position')
	plt.title('STORM: localisation precision at different beat intensities')
	plt.ylabel('stddev of detected position')
	plt.savefig(filename+"_beadVariance.png")

	plt.figure()
	xx,yy=np.hsplit(np.array(beadMeans),np.array([1]))
	plt.plot(xx,yy,'ro')
	plt.title(filename)
	plt.xlabel('x')
	plt.ylabel('y')
	plt.savefig(filename+"_beads.png")
	plt.show()

	if len(beadMeans) > 0:
		meanData_np = np.array(beadMeans)
		outData = np.hstack([meanData_np, np.arange(len(meanData_np)).reshape(-1,1)]) # add index column
		coord2im.writefile(filename+"_beads.txt", dims, outData)
