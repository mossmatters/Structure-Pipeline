error_message = """This script will parse a set of structure results and create a file like the
 "Simulation Summary" that the front-end program produces. List the input files on the
 command line, wildcard expansion should work on Unix-like systems.
 Before the file names, indicate the Max K:
 
 python simsum.py 10 output1_f output2_f output3_f
 
 or
 
 python simsum.py 10 *f"""

def main():
	import sys
	import re
	
	numMatch = re.compile('[-+]?[0-9]*\.?[0-9]+')
	if len(sys.argv) < 2:
		print error_message
		return 0
	else:
		maxK = int(sys.argv[1])
		filelist = sys.argv[2:]
	filestem = filelist[0].split("_")[0]
	outputfile = open(filestem+'.simsum','w')
	header = ['K','LnPr','MeanPr','VarPr']
	
	for k in range(maxK):
		fsthead = 'Fst' + str(k+1)
		header.append(fsthead)
	header.append('\n')
	
	outputfile.write("\t".join(header))
	for file in filelist:
		sim = open(file,'U')
		simsum = []
		for line in sim:
			if line.find('populations assumed') > 0:
				simsum.append(re.findall(numMatch,line)[0])
			elif line.find('Estimated Ln Prob of Data') > -1:
				simsum.append(re.findall(numMatch,line)[0])
			elif line.find('Mean value of ln likelihood') > -1:
				simsum.append(re.findall(numMatch,line)[0])
			elif line.find('Variance of ln likelihood') > -1:
				simsum.append(re.findall(numMatch,line)[0])
			elif line.find('Mean value of Fst') > -1:
				simsum.append(re.findall(numMatch,line)[1])
		simsum.append('\n')
		outputfile.write("\t".join(simsum))
		sim.close()
		
	outputfile.close()
	
if __name__ == '__main__': main()