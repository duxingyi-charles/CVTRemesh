import sys
import re # regexp


def obj2pts(objfile):
	fext = objfile.split('.')[-1]
	if fext!='obj':
		print 'Wrong Input: not obj file!'
	vertices = []
	with open(objfile) as f:
		for line in f.readlines():
			words = re.split(' +', line)
			if words[0] == 'v':
				vertices.append([float(words[1]), float(words[2]), float(words[3])])
	# write pts file
	ptsfile = objfile[:-3] + 'pts'
	with open(ptsfile, "w") as f:
		f.write("%i\n" % len(vertices))
		for v in vertices:
			f.write("%f  %f  %f\n" % (v[0], v[1], v[2]))


# Main
infile = sys.argv[1]
obj2pts(infile)