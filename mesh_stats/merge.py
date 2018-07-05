import sys
import os

header = []
header.append("No.")
data = []

indir = sys.argv[1]
num=0
if (os.path.isdir(indir)):
	for f in os.listdir(indir):
		f = indir + "\\" + f
		if(os.path.isfile(f) and f.split('.')[-1] == 'txt'):
			with open(f) as file:
				data.append([])
				data[-1].append(f.split('_')[-1].split('.')[0])
				words = file.readlines()[0].split(',')
				if (num==0):
					for i in range(len(words)/2):
						header.append(words[2*i])
				for i in range(len(words)/2):
					data[-1].append(words[2*i+1])
				num += 1

	# write
	filename = indir + "\\" + "results.log"
	with open(filename, 'w') as f:
		for h in header:
			f.write(h+',')
		f.write("\n")
		for item in data:
			for d in item:
				f.write(d+',')
			f.write("\n")




