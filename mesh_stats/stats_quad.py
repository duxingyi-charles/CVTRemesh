import sys
import os
import re #regexp
import numpy as np

def load_obj(filename, vertices, faces):
	'''
	load obj file
	vertices: list of vertex coordinates (numpy array)
	faces: list of list of vertex index
	'''
	fext = filename.split('.')[-1]
	if fext!='obj':
		print ("Error: %s is not an obj file.") % filename
		return False

	#
	vertices[:] = []
	faces[:] = []
	with open(filename) as f:
		for line in f.readlines():
			words = re.split(' +', line)
			if words[0]=='v':
				vertices.append(np.array([float(words[1]), float(words[2]), float(words[3])]))
			elif words[0]=='f':
				nf = len(faces)
				faces.append([])
				for w in words[1:]:
					faces[nf].append(int(w.split('/')[0]) - 1)


def calc_tri_quality(t):
	assert len(t) == 3

	edges = [t[(i+1)%3]-t[i] for i in xrange(3)]
	lengths = [np.linalg.norm(e) for e in edges]
	length2 = [l*l for l in lengths]
	mullen = lengths[0] * lengths[1] * lengths[2]
	sumlen = lengths[0] + lengths[1] + lengths[2]
	sumlen2 = length2[0] + length2[1] + length2[2]

	tri_area = 0.5 * np.linalg.norm(np.cross(edges[0], edges[1]))
	sin = [2*tri_area/(lengths[i]*lengths[(i+1)%3]) for i in xrange(3)]
	sin2 = [s*s for s in sin]

	return (16 * tri_area * tri_area) / (mullen * sumlen)
	#return 4 * np.sqrt(3) * tri_area / sumlen2
	#return np.array(sin2).min()

def calc_quad_quality(q):
	'''
	q: list of 4 vert coordinate of the quad face
	return the quality of the quad face
	'''
	assert len(q) == 4
	scale = (np.sqrt(2) + 1) / 2
	#scale = 2 * np.sqrt(3) / 3
	#scale = 1.5
	tri = [[q[i], q[(i+1)%4], q[(i+2)%4]] for i in range(4)]
	triN = [np.cross(t[1]-t[0], t[2]-t[0]) for t in tri]
	triN = [n/np.linalg.norm(n) for n in triN]
	mulNorm = [np.inner(triN[i], triN[j]) for i in range(4) for j in range(4) if i<j]
	planess = 1 - (2/np.pi) * np.arccos(np.array(mulNorm).min())

	triQ = [calc_tri_quality(t) for t in tri]
	min_triQ = np.array(triQ).min()

	return scale * min_triQ * planess




def calc_quad_angles(q):
	'''
	q: list of 4 vert coordinates of the quad face
	return a list of 4 angles of the quad face
	'''
	assert len(q) == 4
	edges = [q[(i+1)%4]-q[i] for i in xrange(4)]
	length = [np.linalg.norm(e) for e in edges]
	dot = [np.vdot(-edges[(i+3)%4], edges[i]) for i in xrange(4)]
	mullen = [length[(i+3)%4]*length[i] for i in xrange(4)]
	ang = np.arccos(np.array(dot)/np.array(mullen))
	ang = ang * (180 / np.pi)
	return ang




def mesh_stats(filename):
	"calc mesh statistics"
	vertices = []
	faces = []
	load_obj(filename, vertices, faces)

	nvert = len(vertices)
	nface = len(faces)

	num_edges = []
	neighbor = [set() for _ in xrange(nvert)]

	quad_quality = []
	quad_angle = []
	quads = []


	for f in faces:
		nedge = len(f)
		# num_edges
		num_edges.append(nedge)
		# neighbor
		for i in xrange(nedge):
			vi = f[i]
			vp = f[(i-1+nedge) % nedge]
			vn = f[(i+1) % nedge]
			neighbor[vi].add(vp)
			neighbor[vi].add(vn)
		# quad quality and quad angle
		if nedge==4:
			quads.append(f[:])
			q = []
			for i in range(4):
				q.append(vertices[f[i]])
			quad_quality.append(calc_quad_quality(q))
			quad_angle.append(calc_quad_angles(q))

	quad_quality = np.array(quad_quality)
	quad_angle = np.array(quad_angle)
	valence = [len(x) for x in neighbor]
	valence = np.array(valence)

	# Quality
	minQ = quad_quality.min()
	avgQ = quad_quality.mean()
	argminQ = quad_quality.argmin()

	# Angle
	max_angle = quad_angle.max(axis=1)
	min_angle = quad_angle.min(axis=1)

	minA = min_angle.min()
	maxA = max_angle.max()
	argminA = min_angle.argmin()
	argmaxA = max_angle.argmax()
	elem_avg_minA = min_angle.mean()
	elem_avg_maxA = max_angle.mean()
	stdA = quad_angle.std()

	# vert singularity histogram
	vs_hist = np.histogram(valence, 9, (-0.5, 8.5))[0]
	vs_hist = vs_hist * (100.0/nvert)
	v4_rate = vs_hist[4]

	# face singularity histogram
	fs_hist = np.histogram(num_edges, 9, (-0.5, 8.5))[0]
	fs_hist = fs_hist * (100.0/nface) 
	f4_rate = fs_hist[4]

	# quad quality histogram
	qq_hist = np.histogram(quad_quality, 50, (0.0, 1.0))[0]
	qq_hist = qq_hist * (100.0/len(quad_quality))

	# quad angle histogram
	qa_hist = np.histogram(quad_angle, 90, (0.0, 180.0))[0]
	assert (quad_angle.size == 4 * len(quad_quality))
	qa_hist= qa_hist * (100.0/(4*len(quad_angle)))


	# Save
	filename = filename[:-4] + '.txt'
	with open(filename, "w") as f:
		f.write("|V|, %i, |F|, %i, minQ, %f, avgQ, %f, minA, %f, elemMinA, %f, maxA, %f, elemMaxA, %f, stdA, %f, v4, %f, f4, %f, \n"
		% (nvert, nface, minQ, avgQ, minA, elem_avg_minA, maxA, elem_avg_maxA, stdA, v4_rate, f4_rate))
		f.write("minQ Quad, " + ','.join([str(int(x)) for x in quads[argminQ]]) + ",\n")
		f.write("minA Quad, " + ','.join([str(int(x)) for x in quads[argminA]]) + ",\n")
		f.write("maxA Quad, " + ','.join([str(int(x)) for x in quads[argmaxA]]) + ",\n")
		f.write("Vert Sing, " + ','.join([str(x) for x in vs_hist]) + ",\n")
		f.write("Face Sing, " + ','.join([str(x) for x in fs_hist]) + ",\n")
		f.write("Quality Hist, " + ','.join([str(x) for x in qq_hist]) + ",\n")
		f.write("Angle Hist, " + ','.join([str(x) for x in qa_hist]) + ",\n")
	# Quad Angles
	ang_file = filename[:-4] + '_ang.csv'
	with open(ang_file, "w") as f:
		f.write(','.join([str(float(x)) for y in quad_angle for x in y]))
	# Quad Quality
	qual_file = filename[:-4] + '_qual.csv'
	with open(qual_file, "w") as f:
		f.write(','.join([str(float(x)) for x in quad_quality]))
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



