# Structure Plot in R, by Matthew G. Johnson. September 1, 2011.
#
# This script will plot the familiar Structure barplot (each bar is an individual).
# Each bar will be an individual with colors corresponding to the inferred cluster.
# The input file is to be taken directly from CLUMPP-- a whitespace-delimited text file
# with NO HEADER. The putative (original) populations should be in column 4, and the Q-matrix
# should start in column 6. 

# You should specify your original population names, and colors for clusters.

# The only thing below the dashed line that might need adjustment is the text line, if your
# original populations have names that are too big.

#You need to change the following lines:
filename = "/Users/mehmattski/Desktop/Sphagnum/Guwass/Guwass_cluster/guwass_calym_k2.outfile"
k = 2
barcolors = c('gray80','gray50','red','green','blue','orange','pink')
poplabels = c('S. triseriporosum','S. guwassanense')

#----------------------------------------------
#You shouldn't need to change these lines
mydata = read.table(filename,header=FALSE)
klast = k + 5
par(mar=c(4,0,0,0))
bp = barplot(t(mydata[,6:klast]),
	col=barcolors,
	space = 0,
	axes=F)
#Automatically grab beginning and end of each population for plotting of lines and labels
popsizes = tapply(mydata$V1,mydata$V4,length)
numpops = length(popsizes)
poploc = 0
popbegins = rep(NA,k)
popends = rep(NA,k)
for(x in 1:numpops){
	popbegins[x] = poploc
	popends[x] = poploc + popsizes[x]
	poploc = poploc + popsizes[x]
}
popmidpoints = (popbegins+popends)/2

#Puts a dark line between each orignial population
abline(v=c(0,popends),lwd=3,xpd=T)
#Label the populations underneath, at the midpoint of each population.
for(x in 1:numpops){
	mtext(poplabels[x],side=1,at=popmidpoints[x],padj=2,cex=1.5,font=4)
}
