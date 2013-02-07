import numpy as np
import vigra
import coords as coord2im
import sys
import time

cutoff=95

def getCandidates(dims,cc,frames):
	#finds good candidates without to much redundancy
	for i in range(cc.shape[0]):
		if cc[i,2]>frames:
			numberPoints=i
			break
	tempArr = cc[0:numberPoints,:]
	tempArr.tolist().sort()
	cc=np.array(tempArr)
	candidates =[]
	candidates.append(cc[0,:])
	counter=0
	limit=5
	addNewPoint=True
	for x,y in cc[1:,:2]:
		counter +=1
		for i in range(len(candidates)):	
			if ((candidates[i][0]-x)**2+ (candidates[i][1]-y)**2)<limit:
				addNewPoint = False
		if addNewPoint:
			candidates.append(cc[counter])
		addNewPoint = True
	print len(cc), len(candidates)
	return np.array(candidates)

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


def detectBeads(dims, cc, cutoff,maxDist=2, maxStdDev=1.3):
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
	if dims[2]>=50:
		initialCandidates = getCandidates(dims, cc, 50)
	else:
		initialCandidates = getCandidates(dims, cc, dims[2])
		
	#initialCandidates = allCandidates[allCandidates[:,2]<100] # from frame 0 to 49
	counter=0
	'''for x, y, frame, intensity in initialCandidates[:,0:4]:
		counter=counter+1
		print counter
		nearestDist, nearestIdx = nearestNeighbor((x,y,intensity), beads, maxDist)
		if nearestDist > 2*maxDist:
			beads.append([x,y,intensity])
		else:
			print "removing beads too close together: ", (x,y,intensity), beads[nearestIdx]
			beads.remove(beads[nearestIdx])'''
	
	for x, y, frame, intensity in initialCandidates[:,0:4]:
		counter=counter+1
		#print counter
		start2=time.time()
		nearestDist, nearestIdx = nearestNeighbor((x,y,intensity), beads, maxDist)
		print 'nearestDist: %1.2f, Bead: (%3.2f, %3.2f)' %(nearestDist, x,y)
		if nearestDist > 2*maxDist:
			beads.append([x,y,intensity])
			#print 'greater', time.time()-start2
		else:
			#print "merging beads too close together: ", (x,y,intensity), beads[nearestIdx]
			beads.append([(x+beads[nearestIdx][0])/2., (y+beads[nearestIdx][1])/2.,(intensity+beads[nearestIdx][2])/2.])
			beads.remove(beads[nearestIdx])
			#print 'smaller', time.time()-start2


	for beadCandidate in beads: # over all candidates
		x,y,intensity = beadCandidate
		roi = [x-maxDist, x+maxDist, y-maxDist, y+maxDist]
		candidate = coord2im.cropROI(cc, roi)[:,:4]
		singleConsidered+=len(candidate)
		#if len(candidate) > numFrames or len(candidate) < cutoff*numFrames: # allowing few percent not detected
		if len(candidate) < cutoff*numFrames: # allowing few percent not detected
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
	
	meanData = deleteSimilarBeads(meanData, maxDist, cc)

	
	return meanData, scatterplotData

def deleteSimilarBeads(mm, maxDist, cc):
	#merges Beads that are too close
	cleanmm = []
	for i in range(len(mm)):
		dist=1000
		if len(cleanmm) == 0:
			cleanmm.append(mm[i])
		else:
			for j in range(len(cleanmm)):
				currdist = (cleanmm[j][0]-mm[i][0])**2+(cleanmm[j][1]-mm[i][1])**2
				if currdist < dist:
					dist = currdist
					index = j
			if dist < 10:
				x = (cleanmm[index][0] + mm[i][0])/2
				y = (cleanmm[index][1] + mm[i][1])/2
				roi = [x-2*maxDist, x+2*maxDist, y-2*maxDist, y+2*maxDist]
				candidate = coord2im.cropROI(cc, roi)[:,:4]
				mm2,stddev,intensitiy = beadVariance(candidate)
				cleanmm[index][0] = mm2[0]
				cleanmm[index][1] = mm2[1]
			else: 
				cleanmm.append(mm[i])
				
	return np.array(cleanmm)

def detectBeadsFromFile(filename, cutoff,  maxDist=2,maxStdDev=1.1):
	dims, cc = coord2im.readfile(filename)
	import time
	start = time.time()
	beadMeans, beadScatter = detectBeads(dims, cc, cutoff,maxDist, maxStdDev)
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
