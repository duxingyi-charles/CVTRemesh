
import sys
import os
import re  #regexp
import numpy as np


def mesh_stats(filename):
	"calc and save mesh statistics"
	#check
	fext = filename.split('.')[-1]
	if (fext != "obj"):
		return False

	#read obj
	with open(filename) as f:
		vert = []
		face = []
		for line in f.readlines():
			data = re.split(' +', line)
			if (data[0] == "v"):
				vert.append(np.array([float(data[1]), float(data[2]), float(data[3])]))
			elif (data[0] == "f"):
				if (data[1].find("//")!=-1):
					for idx in range(1,4):
						data[idx] = data[idx][:data[idx].find("//")]
				face.append([int(data[1]), int(data[2]), int(data[3])])

	#mesh stats
	area = []
	quality = []
	tri_angle = []
	neighbor = [set() for _ in xrange(len(vert))]
	num_obtuse = 0
	num_lt_thirty = 0

	for f in face:
		v1 = vert[f[0]-1]
		v2 = vert[f[1]-1]
		v3 = vert[f[2]-1]
		e1 = v3 - v2
		e2 = v1 - v3
		e3 = v2 - v1
		# Vertex Neighbor
		neighbor[f[0]-1].add(f[1]-1)
		neighbor[f[0]-1].add(f[2]-1)
		neighbor[f[1]-1].add(f[0]-1)
		neighbor[f[1]-1].add(f[2]-1)
		neighbor[f[2]-1].add(f[0]-1)
		neighbor[f[2]-1].add(f[1]-1)
		# Area and Quality
		S = 0.5 * np.linalg.norm(np.cross(e1, e2))
		length = np.array([np.linalg.norm(e1), np.linalg.norm(e2), np.linalg.norm(e3)])
		p = 0.5 * np.sum(length)
		h  = length.max()
		Q = 6 * S / (np.sqrt(3) * p * h)
		area.append(S)
		quality.append(Q)
		# Angle
		dot = np.array([np.vdot(-e1,e2), np.vdot(-e2,e3), np.vdot(-e3,e1)])
		mullen = np.array([length[0]*length[1], length[1]*length[2], length[2]*length[0]])
		##
		div = []
		for i in xrange(3):
			if abs(dot[i]) > abs(mullen[i]):
				break;
			else:
				div.append(dot[i] / mullen[i])
		div = np.array(div)
		if div.size != 3:
			print("Warning: invalid face")
			continue
		##
		ang = np.arccos(dot/mullen)
		ang = ang * (180 / np.pi)
		if (ang.max() > 90):
			num_obtuse += 1
		if (ang.min() < 30):
			num_lt_thirty += 1
		tri_angle.append(np.array(ang))


	area = np.array(area)
	quality = np.array(quality)
	tri_angle = np.array(tri_angle)
	valence = [len(x) for x in neighbor]
	valence = np.array(valence)

	# Quality
	minQ = quality.min()
	avgQ = quality.mean()

	# Area
	area = area / area.mean()
	minS = area.min()
	maxS = area.max()
	stdS = area.std()

	# Angle
	max_angle = tri_angle.max(axis=1)
	min_angle = tri_angle.min(axis=1)

	minA = min_angle.min()
	maxA = max_angle.max()
	elem_avg_minA = min_angle.mean()
	elem_avg_maxA = max_angle.mean()
	stdA = tri_angle.std()

	lt_thirty_rate = num_lt_thirty * 100.0 / len(face)
	obtuse_rate = num_obtuse * 100.0 / len(face)

	# Singularity
	hist = np.histogram(valence, 13, (-0.5, 12.5))[0]
	v567_rate = (hist[5] + hist[6] + hist[7]) * 100.0 / len(vert)
	v6_rate = hist[6] * 100.0 / len(vert)
	#test
	print(hist)

	# read log file for num_iter, time, final energy and gradient
	cut_pos = filename.rfind("rvd")
	logfile = filename[:cut_pos] + filename[cut_pos+4:]
	logfile = logfile[:logfile.rfind(".")] + ".log"

	with open(logfile) as f:
		lines = f.readlines()
		timeline = lines[1].strip('\n')
		times = timeline.split('\t')
		lastline = lines[-2].strip('\n')
		lastIter = lastline.split('\t')


	# Save
	filename = filename[:-3] + 'txt'
	with open(filename, "w") as f:
		f.write("|V|, %-5i, |F|, %-5i, |t_obt|, %-5i, minQ, %-6f, avgQ, %-6f, minS, %-6f, maxS, %-6f, stdS, %-6f, minA, %-6f, maxA, %-6f, elemMinA, %-6f, elemMaxA, %-6f, stdA, %-6f, <30, %-6f, >90, %-6f, V567, %-6f , v6, %-6f ,"
		% (len(vert), len(face), num_obtuse, minQ, avgQ, minS, maxS, stdS, minA, maxA, elem_avg_minA, elem_avg_maxA, stdA, lt_thirty_rate, obtuse_rate, v567_rate, v6_rate))
		f.write("Tcvt, %s, Tfcvt, %s, niter, %s, Ecvt, %s, Gcvt, %s, Edir, %s, Gdir, %s, Gtot, %s, " 
			% (times[2], times[4], lastIter[0], lastIter[1], lastIter[2], lastIter[3], lastIter[4], lastIter[5]))
		f.write("Singularity, "+str(hist) + ",")
	return True


# Main Program

infile = sys.argv[1]
if (os.path.isdir(infile)):
	for f in os.listdir(infile):
		f = infile + '\\' + f
		if (os.path.isfile(f) and f.split('.')[-1] == 'obj'):
			print("processing %s" % f)
			mesh_stats(f)
elif (os.path.isfile(infile) and infile.split('.')[-1] == 'obj'):
	mesh_stats(infile)
else:
	print 'input is neither obj file or directory'





