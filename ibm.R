##### In this file we will try to come up with some things ####
###############################################################%

########################%
####### Type of program
## Written in R
## Object Oriented
#######################%

#######################%
# Main object "Leprechaun": an individual of our species of interest, the Irish Leprechaun, small in size, but with great powers. To make it more French for Timothee, we assume that their favourite food is camembert. Leprechaun is not very choosy, and mates completely random.
##########################################%

# The object contains the following values

# Static [these numbers do not change after initialisation]
#  V ID (unique identifier number) [integer]
#  V pID (two numbers referring to the parents of the Leprechaun, if none: NA) [vector of two integers]
#  V Year of birth (timestep at which the individual was born)	[integer]
#  (- Genome (?) (two vectors of length N coding for both chromosomes of N loci in the genome.) [two vectors of N integers]
#  - Heritable phenotypic trait value of interest (z) (e.g. birth weight) (changes/constant through life depending on trait)
#  - (Possibly: breeding value A))
#  - Rather simulate physically independent loci, otherwise we need to simulate recombination on the chromosomes. Over short time periods, 100 recombinations fragments (ie independent loci) sounds realistic.
#  - Explicit coding of traits by many independent diploid loci. T/Users/gauravbaruha/Documents/New Folder With Items copy/Zurich PhD Research/IBM/MS THESIS DATA/demographics.Rhe simplest model: z = mean + sum_loci(a1_locus + a2_locus) + environment. One can add explicit dominance and epistasis, as well as interactions with environment. 
#  - A possibility is to draw the a of the different alleles from a N(0,V). Each locus can have a different V, and thus a different importance. 
#  - A large number of loci (>20) will give easily patterns expected from quantitative genetics. We can draw randomly the number of loci per trait. 
#  - My main concern at the moment is the initialisation of the genetic diversity: it will be hard to avoid a fast decline of diversity at the beginning. One possibility is using neutral expectations of diversity (n-coalescent or Ewens distribution)
#  - We probably do not need mutations if we consider a population over no more than some tens of generations.
#  V Sex (M/F). Could be genetically determined by on locus, thus allowing random fluctuations of sex ratio and thus population structure.

# Dynamic [these numbers do change after initialisation]
#  V alive (boolean, true/false) [boolean]
#  V age (a) (timesteps since birth) [integer]
#  - stage (possibly instead of/in addition to age) (juvenile, adult, etc)
#  V size (x)

####################%
###### Relations that need to be defined between size (x), age (a), trait (z) and vital rates
# (- Possibly include population density d in functions)
# V Survival(a,x,z,d) (logistic function) 
# - Growth (a,x,z,d) (either transition probability to next stage or absolute growth)
# - Reproduction probability: p_repr(a,x,z,d) (logistic function)
# - Number of offspring: n_offspring(a,x,z,d) (poisson distribution)
# - Offspring size x distribution: x_offspring(x_mother,z_mother,x_father,z_father) (not required if working with stages) (prability density function)
# - Offspring trait z distribution: z_offspring(x_mother,z_mother,x_father,z_father,A_mother,A_father) + rnorm(0,V_E) (prability density function)
# - For survival and reproduction, we can consider functions of the form W=exp(-((z-Optimum)^2)/(2*width)), that is, gaussian stabilizing selection. We could then let the Optimum shift. 

####################
####### Environmental aspects
# - Extra variation in z due to unexplained factors (i.e. V_E)) (assumed constant over time and constant accross all individuals)
# - Changes in the environment affecting survival(x,z), growth(x,z), p_repr(x,z), n_offspring(x,z), f_offspring(x_mother,z_mother) (i.e. changing selection)

########################
########## 'Settings' of ancestral population from which simulation can start:
# - n start individuals
# - Start trait z distribution
# - Start age (/stage) distribution
# - Assign sexes to individuals
# - (Possibly: start a values)

########## Perform in each time step the following actions over all alive individuals at t=0
## In following order:
# survival(x,z)
### For those who survive:
# growth(x,z)
# p_repr(x,z)
### For those who reproduce:
# Random mating between all reproductive males and females
# n_offspring(x,z)
# x_offspring(x_mother,z_mother,x_father,z_father)
# z_offspring(x_mother,z_mother,x_father,z_father,a_mother,a_father,m) + rnorm(0,V_E)
# Random sex assigned to offspring
# Offspring added to population
### End of timestep

####################################################%
################ GLOBAL VARIABLES AND COUNTERS #####
####################################################%

################ Random seed #######################

set.seed(12)

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
# ############### Genetic determinisms  ##############%
# dominance<-1 # for additive effects only, must be 0
# overdominance<-0 # non-null values generate overdominance
# 
# nbLoci<-10 #number of loci controling the trait phenotype
# nbAlleles<-10 #number of existing alleles per loci
# 
# SDZ<-1#standard deviation of locus effect on size
# SDR<-2#standard deviation of locus effect on reproductive success

############### Genetic determinism Z
gvaluesZ<-FunctionGvalues(nbLoci = 10,nbAlleles = 10,dominance = 1,overdominance = 0,SDeffects =1,SDalleles = 1)
############### Genetic determinism H 
gvaluesR<-FunctionGvalues(nbLoci = 10,nbAlleles = 10,dominance = 1,overdominance = 0,SDeffects = 1,SDalleles = 2)
mean(apply(gvaluesZ,1,FUN = function(x){apply(X = x,MARGIN = 1,FUN = var)}))
mean(apply(gvaluesR,1,FUN = function(x){apply(X = x,MARGIN = 1,FUN = var)}))
save(gvaluesZ=gvaluesZ,gvaluesR=gvaluesR,file = "GeneticMatrices")

#################################################################################################################
################# Just to define parameters when playing with the elements of the main function #################
#################################################################################################################
# DIR= "Data/Simple"
# SurvivalSelection1=0.1#linear coefficient on a logit scale for Survival ~ ... + size +size^2
# SurvivalSelection2=0#quadratic coefficient on a logit scale for Survival ~ ... + size + size^2; negative value=balancing selection
# fertilitySelection1=0.1#linear coefficient on a log scale for reproduction ~ ... + size + size^2
# fertilitySelection2=0#quadratic coefficient on a log scale for reproduction ~ ... + size + size^2; negative value=balancing selection
# camembertSelection=0.1#increases the number of offspring, multiplied by Number Of Camemberts at the power 1/3 (because its a volume, uhuh)
# SurvivalPenaltyForRepro=0.2 #trade-off repro survival
# MeanCamembert=20000 # expected number of camemberts available per year
# SDCamembert=10 # standard deviation of the annual camembert availability per year
# CamembertTrend=-200 # change per year in expected number of camembert after the year 10
# MeanBirthSize=10 # hmmm... guess what, this is the mean size at birth
# PlasticityBirthSize=1 # 
# MeanGrowth=1.245 # expected growth factor at birth (this decreases with age to converge to 1. Also, there is some variance around this expectation, that can make individuals shrink)
# SexualMaturity=1 # age a which reproduction becomes possible
# PlasticityReproduction=0 # plasticity in individual fixed quality for reproduction. At the moment there is none, the trait is alreay plastic enough and we are not really interested in it
# MaternalEffect=0.1 # proportion of the mother size inherited (plastically) to the offspring size
# MeanRepro=0.1 # mean annual reproductive success
# StudyLength=50 # total number of years in the study. The first years are likely weird as we start from a single juvenile cohort, the 10 first years could be discarded
# InitialPopSize=200 # guess
# sizeSurvivalInflection=0 # hmmm, no idea, I guess we will delete that
# CamembertthresholdGrowth=20 # minimum number of camemberts necessary for growth??? TO CLARIFY
# GeneticMatrices="GeneticMatrices" # file where genetic matrices are stored

#################################################################################
#################################################################################

LeprechaunSimul<-function(DIR= "Data/Simple",SurvivalSelection1=0.5,SurvivalSelection2=0,fertilitySelection1=0.1,fertilitySelection2=0,
                          camembertSelection=0.1,SurvivalPenaltyForRepro=0.2,MeanCamembert=4000,SDCamembert=10,CamembertTrend=-5,
                          MeanBirthSize=10,MeanGrowth=1.245,SexualMaturity=1, replica=1,
                          PlasticityBirthSize=1,PlasticityReproduction=0,MaternalEffect=1,MeanRepro=0.1,
                          StudyLength=2000,InitialPopSize=70,sizeSurvivalInflection=0,CamembertthresholdGrowth=30,GeneticMatrices="GeneticMatrices"){
   
  load(file = GeneticMatrices)
  nbLociZ<-dim(gvaluesZ)[3]#number of loci targeting Z
  nbAllelesZ<-dim(gvaluesZ)[2]#number of existing alleles per locus targeting Z
  nbLociR<-dim(gvaluesR)[3]#number of loci targeting Reproduction (through unobserved traits)
  nbAllelesR<-dim(gvaluesR)[2]#number of existing alleles per locus targeting Reproduction
    
  
  #folder<-"Simple"
  converter<-paste(DIR,"/conv.csv",sep="")
  #dir.create(file.path("Data", folder))
  cat("filename,SS1,SS2,FS1,FS2,CAMSEL,SDCAM",file=converter,append=FALSE)
  
  
  
  filename<-paste(DIR,"/pop_SS1_",SurvivalSelection1*10,"_SS2_",SurvivalSelection2*100,"_FS1_",fertilitySelection1*10,
                  "_FS2_",fertilitySelection2*100,"_CAMSEL_",camembertSelection*10,"_SDCAM_",SDCamembert,"_replica_",replica, ".csv",sep="")
  cat("\n",filename,"\t",SurvivalSelection1,"\t",SurvivalSelection2,"\t",fertilitySelection1,"\t",fertilitySelection2,"\t",camembertSelection,"\t",SDCamembert,"\t", replica,file=converter,append=TRUE)
  
  
#filename1<-data.frame(DIR,"/pop_SS1_",SurvivalSelection1*10,"_SS2_",SurvivalSelection2*100,"_FS1_",fertilitySelection1*10, "_FS2_",fertilitySelection2*100,"_CAMSEL_",camembertSelection*10,"_SDCAM_",SDCamembert)


#filename<-data.frame(DIR,YR,pop[[i]]@ID,pop[[i]]@size,pop[[i]]@bvs,pop[[i]]@FixRepro,pop[[i]]@bvr,pop[[i]]@camemberts,pop[[i]]@sex,pop[[i]]@ARS,pop[[i]]@age,pop[[i]]@pID[1],pop[[i]]@pID[2])


#a<-data.frame(filenmame,ID,size,camemberts,age,sex,"\t(",object@pID[1],",",object@pID[2],")\t",object@Birth,"\t",object@alive )  
################ Counter for the IDs ###############%
CID<-as.integer(1) 
################ Counter for the current year ######%
YR<-0 
################ Base population parameters ########


####################################################%
############### Definition of the class ############
####################################################%

setClass(
	Class="Leprechaun",
	representation=representation(
		ID = "integer",
		pID = "integer",
		age = "integer",
		Birth = "integer",
		alive = "logical",
		size = "numeric",
    FixRepro = "numeric",
    camemberts = "integer",
		sex = "character",
		DNAZ = "matrix",
		DNAR = "matrix",
    bvs = "numeric",
    bvr = "numeric",
    ARS = "integer"
	)
)

############### Definition of the basic methods (for printing to the screen and initialisation)######
###############################################################################################%
setMethod("show","Leprechaun",
	function(object){
		cat(object@ID,"\t",object@size,"\t",object@camemberts,"\t",object@age,"\t",object@sex,"\t(",object@pID[1],",",object@pID[2],")\t",object@Birth,"\t",object@alive,"\n",sep="")
	}
)

setMethod("initialize","Leprechaun",function(.Object,parent1,parent2,sibs){#parent1 is the mother
  .Object@DNAZ<-matrix(NA,nrow=2,ncol=nbLociZ)
  .Object@DNAR<-matrix(NA,nrow=2,ncol=nbLociR)
	if(missing(parent1)){
		parent1<-NA
		.Object@DNAZ[1,]<-floor(runif(nbLociZ,min=1,max=nbAllelesZ+1))
		.Object@DNAR[1,]<-floor(runif(nbLociR,min=1,max=nbAllelesR+1))
	}else{
			#weight1<-pop[[parent1]]@size
			.Object@DNAZ[1,]<-pop[[parent1]]@DNAZ[cbind(floor(runif(n=nbLociZ, min=1, max=3)), 1:nbAllelesZ)]
			.Object@DNAR[1,]<-pop[[parent1]]@DNAR[cbind(floor(runif(n=nbLociR, min=1, max=3)), 1:nbAllelesR)]
	}
	if(missing(parent2)){
		parent2<-NA
		.Object@DNAZ[2,]<-floor(runif(n = nbLociZ, min=1, max=nbAllelesZ+1))
		.Object@DNAR[2,]<-floor(runif(n = nbLociR, min=1, max=nbAllelesR+1))
	}else{
	  .Object@DNAZ[2,]<-pop[[parent2]]@DNAZ[cbind(floor(runif(n=nbAllelesZ, min=1, max=3)), 1:nbAllelesZ)]
	  .Object@DNAR[2,]<-pop[[parent2]]@DNAR[cbind(floor(runif(n=nbAllelesR ,min=1, max=3)), 1:nbAllelesR)]
	}
  
  if(missing(sibs)){
    sibs <- 1
  }
  
	.Object@age<-as.integer(0)
	.Object@ID<-CID
	.Object@pID<-c(as.integer(parent1),as.integer(parent2))
	.Object@Birth<-as.integer(YR)
	.Object@alive<-TRUE
  BreedingValueSize<-0
  for (Locus in 1:nbLociZ)#take the mean of genetic values
    {
      BreedingValueSize<-BreedingValueSize+(gvaluesZ[ .Object@DNAZ[1,Locus], .Object@DNAZ[2,Locus], Locus]/nbLociZ)
    }
  .Object@bvs<-BreedingValueSize
	size<-MeanBirthSize*(1/(0.875+0.125*sibs))+BreedingValueSize
  
  if(!is.na(parent1)){
  size<-size+MaternalEffect*pop[[parent1]]@size
  }
  .Object@size<-abs(rnorm(n=1,mean=size,sd=PlasticityBirthSize)) # sd plasticity birth size
  
  
	BreedingValueReproduction<-0
	for (Locus in 1:nbLociR)#take the mean of genetic values
	{
	  BreedingValueReproduction<-BreedingValueReproduction+(gvaluesR[ .Object@DNAR[1,Locus], .Object@DNAR[2,Locus], Locus]/nbLociR)
	}
  .Object@bvr<-BreedingValueReproduction
  .Object@FixRepro<-rnorm(n=1,mean=.Object@bvr,sd=PlasticityReproduction) #plasticity reproductive quality
  
  .Object@camemberts<-as.integer(0)
  
  .Object@ARS<-as.integer(0)#annual reproductive success
  
	if(runif(1)>0.5){.Object@sex<-'F'}else{.Object@sex<-'M'}

	CID<<-as.integer(CID+1)
  
	return(.Object)
})

################### Definition of more biologically relevant methods (e.g. survival)######
####################################################################################%

##### Implementing the famous bathtub, ages 1 to 20#####
bathtub<-function(age,shift=6,BaseMortality=0.1){
  p<-BaseMortality*exp(-(age-shift)/4)+(-1+exp((age-shift)*log(2)/(30-shift)))
  p[p>1]<-1
  return(p)
}
#curve(expr =bathtub(x),0,20,add = F,ylim=c(0,1))
#prod((1-bathtub(0:15)))
#prod(sizeSurvival(0:10,size = 10,camemberts = 100,ARS =1))
# Incorporate the effect of size to the survival function
sizeSurvival<-function(age,size,camemberts,ARS){
  #sizedeviation<-size-(MeanBirthSize*prod(Growth[1:(age+1)]))
  p<-bathtub(age)
  if(p<1)#because size does not prevent animals of maximal age to die out
    {
      plogit<-log(p/(1-p))
      Philogit<-plogit-SurvivalSelection1*(size-sizeSurvivalInflection)-SurvivalSelection2*(size-sizeSurvivalInflection)^2+ARS*SurvivalPenaltyForRepro
      p<-exp(Philogit)/(1+exp(Philogit))
      if(camemberts<100)
        {
          p<-1-(1-p)*(camemberts/100)^(1/2)
        }
    }
  return(1-p)
}

# Applying the bathtub in a surival function
setGeneric("Surv",function(Object){standardGeneric("Surv")})

setMethod("Surv","Leprechaun",function(Object){
	
	if(runif(1)>sizeSurvival(Object@age,Object@size,Object@camemberts,Object@ARS)){
		Object@alive<-FALSE
		DEAD<<-c(DEAD,Object@ID)
	}
	
	return(Object)
	
})

# Simple function, simply adds 1 to the age
setGeneric("Age",function(Object){standardGeneric("Age")})

setMethod("Age","Leprechaun",function(Object){
	Object@age<-as.integer(Object@age+1)
	return(Object)
})

# Growth // Sizes change, this fact is known by many - if not all
setGeneric("Grow",function(Object){standardGeneric("Grow")})

setMethod("Grow","Leprechaun",function(Object){
  MeanGrowthAge<-(MeanGrowth+Object@age-1)/Object@age #maximal growth converge to 1 with increasing age
  SdGrowthAge<-(MeanGrowthAge-1)/2
  
  GrowthC<- 1+ (MeanGrowthAge-1)*(2/(1+exp(-0.1*(Object@camemberts-CamembertthresholdGrowth)))-1)#gives the camemebert specific expected growth
  
	Object@size<-Object@size*rnorm(1,mean=GrowthC,sd = SdGrowthAge)# but the animal can always shrink
	
	  if(Object@size > GrowthC )
	
	return(Object)
})

# Retrieving the ID of an individual
setGeneric("IDretrieve",function(Object){standardGeneric("IDretrieve")})

setMethod("IDretrieve","Leprechaun",function(Object){
  return(Object@ID)
})

# Retrieving the size of an individual
setGeneric("Size",function(Object){standardGeneric("Size")})

setMethod("Size","Leprechaun",function(Object){
  return(Object@size)
})

# Retrieving the reproductive quality of an individual
setGeneric("FixRepro",function(Object){standardGeneric("FixRepro")})

setMethod("FixRepro","Leprechaun",function(Object){
  return(Object@FixRepro)
})

# Retrieving the sex of an individual
setGeneric("Sex",function(Object){standardGeneric("Sex")})

setMethod("Sex","Leprechaun",function(Object){
  return(Object@sex)
})

# Calculating the number of offspring for a females
setGeneric("Num_off",function(Object){standardGeneric("Num_off")})

setMethod("Num_off","Leprechaun",function(Object){  
  lambda<-exp(log(MeanRepro)+Object@FixRepro+fertilitySelection1*(Object@size-MeanBirthSize)/10+fertilitySelection2*((Object@size-MeanBirthSize)/10)^2
                +camembertSelection*((Object@camemberts)^(1/3))/10)
  repro<-rpois(n=1,lambda=lambda)+1
  
  logitV<-Object@age-SexualMaturity
  p<-1/(1+exp(-logitV))
  Object@ARS<-as.integer(repro)*rbinom(1,size = 1,prob = p)#*disaster[YR]
	return(Object)
})


# Function for camembert attributions 
setGeneric("Food",function(Object,Camams){standardGeneric("Food")})

setMethod("Food","Leprechaun",function(Object,Camams){
  Object@camemberts<-as.integer(Camams)
  return(Object)
})

############### Creating an initial population with 100 individuals#######
pop<-c(new("Leprechaun"))
for(i in 2:InitialPopSize){
	pop<-c(pop,new("Leprechaun"))
}
############### Function for printing values########
print_info<-function(YR,ALIVE,CID,camembert){
	cat("\nAt the beginning of year:",YR,"\nThere are:",length(ALIVE),"Leprechauns\n-----------------\n")
	cat("ALIVE:",ALIVE,"\n")
	cat("Current ID:",CID,"\n")
  	cat("\n-----------------\nThe camembert production is:",camembert,"\n")
	return()
}

############### List of living individuals [their indices], this will save time later, because dead individuals are not looped over
ALIVE<-1:length(pop)

cat("t,ID,z,bvs,FixRepro,bvr,C,s,ARS,age,p1,p2,phi",file=filename,append=FALSE)

#disaster<-as.integer(rep(c(1,1,1),times = 1+StudyLength/3))


############### The start of time
for(YR in 1:StudyLength){
  
  ExpectedCamembert<-ifelse(test = YR>200, yes = MeanCamembert+(YR-200)*CamembertTrend, no = MeanCamembert)#after the year 10, availability decreases 
  if(ExpectedCamembert<0){
    ExpectedCamembert<-1
  }

  camembert<-abs(round(rnorm(n=1,mean=ExpectedCamembert,sd=ExpectedCamembert/SDCamembert),digits=0)) # resources for year YR
  
  #print_info(YR,ALIVE,CID,camembert)
	#### Competition for resources
  sizes<-unlist(lapply(pop[ALIVE], function(x) Size(x)))
  
  
  HunterQualities <- runif(n = length(ALIVE),min = 0.5,max = 1)  +sizes/20  #+ (rnorm(length(ALIVE), mean=1,sd=0.2))
  
  if(length(HunterQualities)==1){HunterQualities=1}
  
  #podium<-table(factor(sample(as.character(ALIVE),size=camembert,replace=T,prob=HunterQualities),levels=ALIVE))
  #camams<-as.numeric(podium[match(ALIVE,names(podium))])
  camams<-round(camembert*HunterQualities/sum(HunterQualities),digits = 0)
  #pop[ALIVE]<-lapply(1:length(ALIVE), function(x) Food(pop[[ALIVE[x]]],camams[x]))
  pop[ALIVE]<-lapply(1:length(ALIVE), function(x) Food(pop[[ALIVE[x]]],camams[x]))
  
  #### Survival
  DEAD<-c()
	pop[ALIVE]<-lapply(pop[ALIVE],Surv)
	ALIVE<-ALIVE[!(ALIVE %in% DEAD)]
	
  if(length(ALIVE)==0){
    for(i in DEAD){
      cat("\n",YR,",",pop[[i]]@ID,",",pop[[i]]@size,",",pop[[i]]@bvs,",",pop[[i]]@FixRepro,",",pop[[i]]@bvr,",",pop[[i]]@camemberts,",",pop[[i]]@sex,",",pop[[i]]@ARS,",",pop[[i]]@age,",",pop[[i]]@pID[1],",",pop[[i]]@pID[2],",",0,file=filename,append=TRUE)
    }
    break
  }
	#### Age+1 and growth
	pop[ALIVE]<-lapply(pop[ALIVE],Age)
	pop[ALIVE]<-lapply(pop[ALIVE],Grow)
	
	#### Reproduction  ### Not the easiest part (but necessary if we want to allow the population to grow)
	
	##########
	### Part dedicated to retrieving the indices of all living males and of all living females
	##########
	males<-lapply(pop,Sex)=="M" # Determine which individuals are males -- logicaly all other individuals should be females... However, this includes dead ones... Simply a list of T,T,F,F,T,F,F,....
	females<-which(!males) # Get the indices of the non-males (that is females..)
	males<-which(males)   # Get the indices of the males
	females<-intersect(females,ALIVE) # Retrieve the indices of the living(!) females
	males <-intersect(males,ALIVE) # Retrieve the indises of the living males
	#cat(females)
	
	##########
	# Part dedicated to breeding.. 
	##########
	
	# We take a female based approach: we determine for each females 
	from<-CID
	pop[females]<-lapply(pop[females],FUN = Num_off)
	for(i in females){
		Noffs<-pop[[i]]@ARS
		if(Noffs>0 & length(males>0)){
			#Determine the father
			fat<-sample(males,1)
			for(j in 1:Noffs){
				#cat(j,"\n")
				# Create the offspring
				pop<-c(pop,new("Leprechaun",parent1=i,parent2=fat,sibs=Noffs))
			}
		}		
	}
  
	if(from!=CID){
	ALIVE<-c(ALIVE,(from):(CID-1))
	}
	### Everything should be written to a dataframe, to make sure we have all the values for ever and ever
  for(i in ALIVE){
    cat("\n",YR,",",pop[[i]]@ID,",",pop[[i]]@size,",",pop[[i]]@bvs,",",pop[[i]]@FixRepro,",",pop[[i]]@bvr,",",pop[[i]]@camemberts,",",pop[[i]]@sex,",",pop[[i]]@ARS,",",pop[[i]]@age,",",pop[[i]]@pID[1],",",pop[[i]]@pID[2],",",1,file=filename,append=TRUE)
  }
  for(i in DEAD){
    cat("\n",YR,",",pop[[i]]@ID,",",pop[[i]]@size,",",pop[[i]]@bvs,",",pop[[i]]@FixRepro,",",pop[[i]]@bvr,",",pop[[i]]@camemberts,",",pop[[i]]@sex,",",pop[[i]]@ARS,",",pop[[i]]@age,",",pop[[i]]@pID[1],",",pop[[i]]@pID[2],",",0,file=filename,append=TRUE)
  }
}
 return(filename)
}
#end of LeprechaunSimul

# THE COMMAND FOR RUNNING THE CODE: LeprechaunSimul(StudyLength=50,DIR="Data")


