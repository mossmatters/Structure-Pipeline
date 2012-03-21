# This script takes as input a file of replicate STRUCTURE 
# runs for many values of K. For each run, the value of K and the 
# posterior probability "Ln(P|D)" column should be labeled "LnPr" in a tab-delimited file.
# Delta-K is calculated by the method of Evanno et al 2005 Mol. Ecol. 14, 2611-20:
# DeltaK = abs(L(K+1) - 2*L(K) + L(K-1)) / sd(L(K))

inputfile = "/Users/mehmattski/Desktop/Sphagnum/Guwass/Guwass_cluster/100K1M.simsum"
simsum = read.delim(inputfile)

maxK = max(simsum$K)

Kmeans = tapply(simsum$LnPr,simsum$K,mean)
Kvars = tapply(simsum$LnPr,simsum$K,var)
Ksds = Kvars^0.5

deltaK = rep(NA,maxK)

for(x in 2:maxK){deltaK[x] = abs(Kmeans[x+1] - 2*Kmeans[x] + Kmeans[x-1])/Ksds[x]}
plot(deltaK,type = 'b',xlab='K',cex.axis=0.75)

#Use this to get the plot of likelihood versus K
#plot(LnPr~K,data=simsum,ylab='L(P|D)')