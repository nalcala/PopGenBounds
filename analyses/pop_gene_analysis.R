library(PopGenBounds)
library(tibble)


# Let's first create an allele frequency matrix with subpopulations as rows and alleles as columns:

popa = matrix(c(0.4,0.4,0.2,rep(0,3*5),
                rep(0,3),0.4,0.4,0.2,rep(0,3*4),
                rep(0,3*2),0.4,0.4,0.2,rep(0,3*3),
                rep(0,3*3),0.4,0.4,0.2,rep(0,3*2),
                rep(0,3*4),0.4,0.4,0.2,rep(0,3),
                rep(0,3*5),0.4,0.4,0.2),ncol=18,byrow = T)
# The Diff function computes the differentiation statistics along with their bounds given the most frequent allele M
Diff(popa)

# We can plot the allele frequency table with function ggfreqtable:
ggfreqtable(popa)

#And plot the value within the bounds of FST, G'ST, and D with function ggbounds:
ggbounds(M=Diff(popa)$value[1],FST=Diff(popa)$value[2],GpST = Diff(popa)$value[3],D=Diff(popa)$value[4],K=nrow(popa))

# From vignette: Yellow Toad populations

Toad_microsat_freq = PopGenBounds::Toad_microsat_freq

Diff(Toad_microsat_freq[[2]])
# Plotting the allele frequency tables
ggfreqtable(Toad_microsat_freq[[2]])


# Computing the statistics
#For all subpopulations together (K=47)

Diff_toad = lapply(Toad_microsat_freq,Diff)

# For pairs of subpopulations (K=2)
Diff_toad_pairs = lapply(Toad_microsat_freq,function(f){sapply(1:46, function(i){ sapply((i+1):47, function(j) Diff(f[c(i,j),])$value)})} )

# Plotting the data
#For K=47
toad.tib = tibble(M=unlist( sapply(Diff_toad,function(x){x[1,2]})),
                  FST=unlist( sapply(Diff_toad,function(x){x[2,2]})),
                  GpST=unlist( sapply(Diff_toad,function(x){x[3,2]})),
                  D=unlist( sapply(Diff_toad,function(x){x[4,2]})))

ggtoad = ggbounds(M=toad.tib$M,FST=toad.tib$FST,GpST=toad.tib$GpST,D=toad.tib$D,K=47)

ggtoad[[1]]+ ggtoad[[2]]+ ggtoad[[3]]


