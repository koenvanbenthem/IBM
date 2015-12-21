# ############### Genetic determinisms  ##############%
# dominance<-1 # for additive effects only, must be 0
# overdominance<-0 # non-null values generate overdominance
# 
# nbLoci<-10 #number of loci controling the trait phenotype
# nbAlleles<-10 #number of existing alleles per loci
# 
# SDZ<-1#standard deviation of locus effect on size
# SDR<-2#standard deviation of locus effect on reproductive success

FunctionGvalues<-function(nbLoci=10,nbAlleles=10,dominance=0.5,overdominance=0,SDeffects=1,SDalleles=1)
{
  gvalues<-array(data=NA,dim=c(nbAlleles,nbAlleles,nbLoci),dimnames=list(paste("A",1: nbAlleles,sep=""),paste("A",1: nbAlleles,sep=""),paste("L",1:nbLoci,sep=""))) # Initialising a matrix that will contain the genotypic effects on the/a trait
  for(L in 1:nbLoci)
  {
    # Setting the effects for the homozygotes [all loci]
    effect<-abs(rnorm(n=1,mean=0,sd=SDeffects))# alter the locus importance in a realistic way (many small-effect loci, few major loci)
    diag(gvalues[,,L])<-2*rnorm(n=dim(gvalues)[1],mean=0,sd=effect*SDalleles)
    # Setting the effects for the heterozygotes
    for(A in 1:(nbAlleles-1))# loop for off-diagonal = heterozygotes (additive and dominance effects)
    {
      for (D in (A+1):nbAlleles)
      {
        d<-dominance*runif(n=1,min=-0.5-overdominance,max=0.5+overdominance)
        gvalues[A,D,L]<-(0.5-d)*gvalues[A,A,L]+(0.5+d)*gvalues[D,D,L] # mean of additive effects + dominance, over diagonal
        gvalues[D,A,L]<-(0.5-d)*gvalues[A,A,L]+(0.5+d)*gvalues[D,D,L] # the same below diagonal    
      }
    }
  }
  return(gvalues)
}

