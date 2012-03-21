#This script takes the output of replicate runs of Structure and formats
#them for CLUMPP. The first argument is an output name. Files may be entered at the command
#line with the wildcard expansion (on UNIX systems). 

import sys


if len(sys.argv) < 2:
	print "Usage: python structure2clumpp.py outputfilename inputfile1 inputfile2 ..."
filenames = sys.argv[2:]
outfilename = sys.argv[1]
outfile = open(outfilename,'w')

def is_line_Qmatrix(line):
	if line.find("%Miss") > 0:
		return 2
	
	if len(line.split()) > 0:
		try:
			int(line.rstrip().split()[0])
			return 1
		except ValueError:
			return 0
	else:
		return 0


for file in filenames:
	run = open(file,'U')
	print file
	Qmatrix = False
	idNum = 1
	for line in run.readlines():
		
		Q = is_line_Qmatrix(line)
		if Qmatrix:
			if Q == 1:
				line = line.rstrip().split()
				line[1] = str(idNum)
				line = "\t".join(line)+'\n'
				outfile.write(line)
				idNum +=1
			elif Q == 0:
				Qmatrix = False
		if Q == 2:
			Qmatrix = True
		
	run.close()

outfile.close()