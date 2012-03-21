# Last updated: 11/4/2010

#input: tab-delimited text file of GENALEX format, with a header denoting the locus of each column:

# Nloci	Nindiv	Npop	NPop1	NPop2	NPop3
# dataid			Pop1name	Pop2name	Pop3name
# ID	Pop	1	1	1	2	2	3	3	3
# MJ1	A	124	125	126	89	90	250	240	230

#So it will be expecting two columns of ID data followed by the loci.

#Supports the following program files:

#Structure
#PopDist
#Multilocus
#Structurama
#Prabclus
#AWclust
#Migrate-n

#Supports multiple ploidy levels, as determined by the header row.
#NOTE: not all programs support ploidy > 2!!
#If an allele is double-digits, will convert to three digit number in Prabclus/GenoDive
	#Missing data will be converted to 000 in Prabclus/GenoDive, ? in Structurama and Migrate
	
#At the command line, need to specify only the input file name and the output file formats. 

import sys

class metadata:
	"""loads data from various file sources into a meta-format, which can be called by other functions for output"""

	def __init__(self,filetype):
		self.inTitle = 0
		self.IDs = []
		self.Pops = []
		self.genotypes = []
		self.filetype = filetype
		
	def loadIDs(self, ID):
		self.IDs.append(ID)
		
	def loadPops(self, Pop):
		self.Pops.append(Pop)
	
	def loadgenotypes(self, genotype):
		self.genotypes.append(genotype)
		
	def filestats(self,NumGenotypes,numLoci,numPops,ploidy=0,loci_names=None,analysis_ID=None):
		self.ploidy = ploidy
		self.numGenotypes = NumGenotypes
		self.numLoci = numLoci
		self.loci_names = loci_names
		self.analysis_ID = analysis_ID
		self.numPops = numPops

###############COMMON FUNCTIONS############################

def getploidy(header):
	loci_names = []
	column_count = []

	for column in xrange(len(header)):
		if header[column] not in loci_names:
			loci_names.append(header[column])
	for locus in loci_names:
		column_count.append(header.count(locus))
		
	ploidy = max(column_count)


	return ploidy, column_count,loci_names

def get_genotype(individual,column_count):
	genotype = []
	column_position = 0
	
	for step in range(len(column_count)):
		num_alleles = column_count[step]
		locus = individual[column_position:column_position+num_alleles]
		genotype.append(locus)
		column_position += num_alleles	
	return genotype
	

def missing_replace(alleles,missing_character):
	"""Input is list of alleles. Missing data (0 in GenAlEx) will be replaced with program-specific character"""
	for allele in range(len(alleles)):
		if alleles[allele] == '0':
			alleles[allele]= missing_character
	return alleles
	
def doubledigit_fix(alleles):
	"""Fixes doubledigit alleles to triple digits for formats like genepop."""
	for allele in range(len(alleles)):
		if len(alleles[allele]) < 3:
			alleles[allele] = "0%s" % (alleles[allele])
	return alleles

###################INPUT FUNCTIONS##########################

def genalex_in(inputfilename):
	input = metadata('GENALEX')

	infile = open(inputfilename,'U')
	rawdata = infile.readlines()
	
	header1=rawdata.pop(0).rstrip().split('\t')
	numLoci = int(header1[0])
	numGenotypes = int(header1[1])
	numPops = int(header1[2])
	
	header2=rawdata.pop(0).rstrip().split('\t')
	
	header3=rawdata.pop(0).rstrip().split('\t')
	ploidy,col_count,loci_names = getploidy(header3[2:])
	
	input.filestats(numGenotypes,numLoci,numPops,ploidy=ploidy,loci_names=loci_names,analysis_ID=header2[0])
	
	if len(rawdata) != numGenotypes:
		print "%d genotypes found. Expected %d genotypes from file header. Check and re-try!" % (len(rawdata), numGenotypes)
		return 0
	
	for row in xrange(len(rawdata)):
		rowlen = len(rawdata[row].rstrip().split('\t'))
		if rowlen-2 != numLoci*ploidy:
			print "Header indicates ploidy = %d, Number of Loci = %d. Expecting %d columns, found %d columns for individual %d. Check and retry!" %(ploidy, numLoci, ploidy*numLoci+2, rowlen, row)   
			return 0
			
	for row in xrange(len(rawdata)):
		individual = rawdata[row].rstrip().split('\t')
		input.loadIDs(individual.pop(0))
		input.loadPops(individual.pop(0))
		genotype = get_genotype(individual,col_count)
		input.loadgenotypes(genotype)
	
	print "GenAlEx input complete!"
	return input

#######################OUTPUT FUNCTIONS######################


def structure_out(inputdata,outputfileroot):
	IDs = inputdata.IDs
	numGenotypes = inputdata.numGenotypes
	genotypes =inputdata.genotypes
	loci_names = inputdata.loci_names
	analysis_ID = inputdata.analysis_ID
	
	outputfilename = "%s.structure" %(outputfileroot)
	outputfile = open(outputfilename,'w')
	
	locistring = "\t".join(loci_names) + "\n"
	outputfile.write(locistring)
	
	if len(inputdata.Pops) > 0:
		numPops = inputdata.numPops
		Pops_names = inputdata.Pops
		Pops = []
		namestonums = {}
		for name in range(len(Pops_names)):
			if Pops_names[name] in namestonums:
				Pops.append(str(namestonums[Pops_names[name]]))
			else: 
				namestonums[Pops_names[name]] = len(namestonums) + 1
				Pops.append(str(namestonums[Pops_names[name]]))
	
	for genotype in xrange(numGenotypes):
		ID = IDs[genotype]
		Pop = Pops[genotype]
		individual = genotypes[genotype]
		
		loci = []
		for locus in range(len(individual)):
			alleles = [x for x in individual[locus]]
			
			alleles_out = "\t".join(alleles)
			loci.append(alleles_out)
		ind_out = "%s\t%s\t%s\n" %(ID,Pop,"\t".join(loci))
		outputfile.write(ind_out)

			
	print "Structure output complete!"
	outputfile.close()
	return

def structurama_out(inputdata,outputfileroot):
	missing = "?"
	numLoci = inputdata.numLoci
	numGenotypes = inputdata.numGenotypes
	IDs = inputdata.IDs
	genotypes =inputdata.genotypes
	outputfilename = "%s.structurama" %(outputfileroot)
	outputfile = open(outputfilename,'w')
	
	
	Header = """#NEXUS
	
	begin data;
		dimensions nind=%d nloci=%d;
		info\n""" %(numGenotypes,numLoci)
	outputfile.write(Header)
	
	
	for genotype in xrange(numGenotypes):
		ID = IDs[genotype]
		individual = genotypes[genotype]
		
		loci = []
		for locus in range(len(individual)):
			alleles = [x for x in individual[locus]]
			alleles = missing_replace(alleles,missing)
			
			alleles_out = "(%s)" %(",".join(alleles))
			loci.append(alleles_out)
		
		if genotype == numGenotypes-1:
			ind_out = "\t%s\t%s;\n" %(ID," ".join(loci))
		else:
			ind_out = "\t%s\t%s,\n" %(ID," ".join(loci))
		outputfile.write(ind_out)
		
	Footer = "end;"
	
	outputfile.write(Footer)

			
	print "Structurama output complete! Note: UNIX line endings."
	outputfile.close()
	return


def popdist_out(inputdata,outputfileroot):
	missing = "000"
	numLoci = inputdata.numLoci
	numGenotypes = inputdata.numGenotypes
	IDs = inputdata.IDs
	genotypes =inputdata.genotypes
	loci_names = inputdata.loci_names
	analysis_ID = inputdata.analysis_ID
	
	outputfilename = "%s.popdist" %(outputfileroot)
	outputfile = open(outputfilename,'w')
	
	
	Header = analysis_ID + "\n"
	outputfile.write(Header)
	for locus in loci_names:
		outputfile.write(locus+"\n")
	
	for genotype in xrange(numGenotypes):
		ID = IDs[genotype]
		individual = genotypes[genotype]
		
		loci = []
		for locus in range(len(individual)):
			alleles = [x for x in individual[locus]]
			alleles = missing_replace(alleles,missing)
			alleles = doubledigit_fix(alleles)
			
			alleles_out = "".join(alleles)
			loci.append(alleles_out)
		ind_out = "%s,\t%s\n" %(ID,"\t".join(loci))
		outputfile.write('POP\n')
		outputfile.write(ind_out)

			
	print "Popdist output complete!"
	outputfile.close()
	return
def multilocus_out(inputdata,outputfileroot):
	missing = "?"
	numLoci = inputdata.numLoci
	numGenotypes = inputdata.numGenotypes
	IDs = inputdata.IDs
	genotypes =inputdata.genotypes
	loci_names = inputdata.loci_names
	ploidy = inputdata.ploidy
	
	outputfilename = "%s.multilocus" %(outputfileroot)
	outputfile = open(outputfilename,'w')
	
	
	for genotype in xrange(numGenotypes):
		ID = IDs[genotype]
		individual = genotypes[genotype]
		
		loci = []
		for locus in range(len(individual)):
			alleles = [x for x in individual[locus]]
			alleles = missing_replace(alleles,missing)
			alleles_out = "/".join(alleles)
			loci.append(alleles_out)
		ind_out = "%s\n" %("\t".join(loci))
		outputfile.write(ind_out)

			
	print "Multilocus output complete!"
	outputfile.close()

	return
	
	return
def prabclus_out(inputdata,outputfileroot):
	missing = "000"
	numLoci = inputdata.numLoci
	numGenotypes = inputdata.numGenotypes
	IDs = inputdata.IDs
	genotypes =inputdata.genotypes
	loci_names = inputdata.loci_names
	ploidy = inputdata.ploidy
	
	outputfilename = "%s.prabclus" %(outputfileroot)
	outputfile = open(outputfilename,'w')
	
	
	for genotype in xrange(numGenotypes):
		ID = IDs[genotype]
		individual = genotypes[genotype]
		
		loci = []
		for locus in range(len(individual)):
			alleles = [x for x in individual[locus]]
			alleles = missing_replace(alleles,missing)
			alleles = doubledigit_fix(alleles)
			if ploidy == 1:
				alleles_out = "".join(alleles+alleles)
			else:
				alleles_out = "".join(alleles)
			loci.append(alleles_out)
		ind_out = "%s,\t%s\n" %(ID,"\t".join(loci))
		outputfile.write(ind_out)

			
	print "Prabclus output complete!"
	outputfile.close()

	return
def genepop_out(inputdata,outputfileroot):
	return


def awclust_out(inputdata,outputfileroot):
	"""Converts the file into allele-based presence absence for use in the R module AWclust"""
	missing_data = '-1'
	ID_List = inputdata.IDs
	dataset = inputdata.genotypes
	ploidy = inputdata.ploidy
	outputfilename = "%s.awclust" %(outputfileroot)
	numLoci = inputdata.numLoci
	alleles = []
	for x in range(numLoci):
		alleles.append([])
	
	for line in dataset:
		for locus in range(numLoci):
			for allele in range(ploidy):
				if int(line[locus][allele]) != 0:
					if int(line[locus][allele])  not in alleles[locus]:
						alleles[locus].append(int(line[locus][allele]))
						
	AWclust = []
	for individual in range(len(ID_List)):
		AWclust.append([])
		for locus in range(numLoci):
			for allele in alleles[locus]:
				if int(dataset[individual][locus][0]) == 0:
					allelecount = missing_data
				else:
					allelecount = dataset[individual][locus].count(str(allele))
				AWclust[individual].append(str(allelecount))
				
	AWclust_flip = zip(*AWclust)

	out = open(outputfilename,'w')
	out.write(" ".join(ID_List)+'\n')
	for line in range(len(AWclust_flip)):
		out.write(" ".join(AWclust_flip[line])+'\n')
	out.close()
	
	print "AWClust output complete!"
	
	return

def migrate_out(inputdata,outputfileroot):
	"""Converts the file for migrate-n. Populations are split as in input file."""
	missing = "?"
	numLoci = inputdata.numLoci
	numGenotypes = inputdata.numGenotypes
	IDs = inputdata.IDs
	genotypes =inputdata.genotypes
	#loci_names = inputdata.loci_names
	ploidy = inputdata.ploidy
	pops = inputdata.Pops
	
	popNames = list(set(pops))
	popSplit = []
	
	header = " %d %d . %s\n" %(len(popNames),numLoci,outputfileroot)
	outputfilename = "%s.mig" %(outputfileroot)
	outputfile = open(outputfilename,'w')
	outputfile.write(header)
	
	for x in range(len(popNames)):
		popname = popNames[x]
		popSplit.append([])
		for id_index in [id_index for id_index,id_pop in enumerate(pops) if id_pop == popNames[x]]:
			popSplit[x].append(id_index)
		
		pop_size = len(popSplit[x])
		
		pop_header = "%d %s\n" %(pop_size,popname)
		outputfile.write(pop_header)
		
		for genotype in xrange(pop_size):
			record_index = popSplit[x][genotype]
			ID = IDs[record_index]
			if len(ID) < 10:
				ID = ID + '_'*(10-len(ID))
			elif len(ID) > 10:
				ID = ID[0:9]
			
			individual = genotypes[record_index]
			
			loci = []
			for locus in range(len(individual)):
				alleles = [y for y in individual[locus]]
				alleles = missing_replace(alleles,missing)
				if ploidy == 1:
					alleles_out = [y + ".?" for y in alleles][0]
				else:
					alleles_out = ".".join(alleles)
				loci.append(alleles_out)
			ind_out = "%s %s\n" %(ID," ".join(loci))
			outputfile.write(ind_out)
			
	print "Migrate output complete"
	outputfile.close()
	return

programs = {
				"STRUCTURE":structure_out,
				"STRUCTURAMA":structurama_out,
				"POPDIST":popdist_out,
				"MULTILOCUS":multilocus_out,
				"PRABCLUS":prabclus_out,
				"GENEPOP":genepop_out,
				"AWCLUST":awclust_out,
				"MIGRATE":migrate_out
				
				}
def main():
	
	if len(sys.argv) < 3:
		print "\n\nProper usage: \npython popgen_parse.py inputfilename outputfiletype1 outputfiletype2..."
		print "\n\nSupported filetypes:"
		for item in sorted(programs.iterkeys()):
			print item
		return 0
	
	inputfilename = sys.argv[1]
	
	try:
		open(inputfilename)
	except IOError:
		print "File %s not found! Check spelling/captialization and try again!" % (inputfilename)
		return 0
	
	if "." in inputfilename:
		inputroot = inputfilename.split(".")[0]
	else:
		inputroot = inputfilename
	
	outputformats = [x.upper() for x in sys.argv[2:]]
	
	inputdata = genalex_in(inputfilename)	
	
	if inputdata:
		for format in outputformats:
			programs[format](inputdata,inputroot)
	
		
if __name__ == '__main__': main()