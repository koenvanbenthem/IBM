########################%
####### Type of program
## Written in R
## Object Oriented
#######################%

LeprechaunSimul<-function(DIR= "Data/Simple",datafile="dataset",SurvivalSelection1=0.5,SurvivalSelection2=0,FertilitySelection1=0.1,FertilitySelection2=0,
                          FoodSelection=0.1,SurvivalPenaltyForRepro=0.2,MeanFood=4000,SDFood=10,FoodTrend=-5,
                          MeanBirthSize=10,MeanGrowth=1.245,SexualMaturity=1, Replica=1,
                          PlasticityBirthSize=1,PlasticityReproduction=0,MaternalEffect=1,MeanRepro=0.1,
                          StudyLength=2000,InitialPopSize=70,SizeSurvivalInflection=0,FoodThresholdGrowth=30,GeneticMatrices="GeneticMatrices"){
   
  load(file = "GeneticMatrices.RData")
  nbLociZ<-dim(gvaluesZ)[3]#number of loci targeting Z
  nbAllelesZ<-dim(gvaluesZ)[2]#number of existing alleles per locus targeting Z
  nbLociR<-dim(gvaluesR)[3]#number of loci targeting Reproduction (through unobserved traits)
  nbAllelesR<-dim(gvaluesR)[2]#number of existing alleles per locus targeting Reproduction

  source("class.R",local=TRUE) # Loads in the class and function definitions - this is done per run to adapt to changing sizes of genomes

################ Counter for the IDs ###############%
CID<-as.integer(1) 
################ Counter for the current year ######%
YR<-0 
################ Base population parameters ########

############### Creating an initial population with 100 individuals#######
pop<-c(new("Leprechaun"))
for(i in 2:InitialPopSize){
	pop<-c(pop,new("Leprechaun"))
}
############### Function for printing values########
print_info<-function(YR,ALIVE,CID,food){
	cat("\nAt the beginning of year:",YR,"\nThere are:",length(ALIVE),"Leprechauns\n-----------------\n")
	cat("ALIVE:",ALIVE,"\n")
	cat("Current ID:",CID,"\n")
  	cat("\n-----------------\nThe food production is:",food,"\n")
	return()
}

filename<-paste(DIR,"/",datafile,".csv",sep="")

############### List of living individuals [their indices], this will save time later, because dead individuals are not looped over
ALIVE<-1:length(pop)

cat("t,ID,z,bvs,FixRepro,bvr,food,s,ARS,age,hetero,p1,p2,phi",file=filename,append=FALSE)

############### The start of time
for(YR in 1:StudyLength){
  
  ExpectedFood<-ifelse(test = YR>200, yes = MeanFood+(YR-200)*FoodTrend, no = MeanFood)#after the year 10, availability decreases 
  if(ExpectedFood<0){
    ExpectedFood<-1
  }

  food<-abs(round(rnorm(n=1,mean=ExpectedFood,sd=SDFood/ExpectedFood),digits=0)) # resources for year YR
  
  #print_info(YR,ALIVE,CID,camembert)
	#### Competition for resources
  sizes<-unlist(lapply(pop[ALIVE], function(x) Size(x)))
  
  
  HunterQualities <- runif(n = length(ALIVE),min = 0.5,max = 1)  +sizes/20  #+ (rnorm(length(ALIVE), mean=1,sd=0.2))
  
  if(length(HunterQualities)==1){HunterQualities=1}
  
  portions<-round(food*HunterQualities/sum(HunterQualities),digits = 0)
  pop[ALIVE]<-lapply(1:length(ALIVE), function(x) Food(pop[[ALIVE[x]]],portions[x]))
  
  #### Survival
  DEAD<-c()
	pop[ALIVE]<-lapply(pop[ALIVE],Surv)
	ALIVE<-ALIVE[!(ALIVE %in% DEAD)]
	
  if(length(ALIVE)==0){
    for(i in DEAD){
      cat("\n",YR,",",pop[[i]]@ID,",",pop[[i]]@size,",",pop[[i]]@bvs,",",pop[[i]]@FixRepro,",",pop[[i]]@bvr,",",pop[[i]]@camemberts,",",pop[[i]]@sex,",",pop[[i]]@ARS,",",pop[[i]]@age,",",pop[[i]]@hetero,",",pop[[i]]@pID[1],",",pop[[i]]@pID[2],",",0,file=filename,append=TRUE)
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
    cat("\n",YR,",",pop[[i]]@ID,",",pop[[i]]@size,",",pop[[i]]@bvs,",",pop[[i]]@FixRepro,",",pop[[i]]@bvr,",",pop[[i]]@food,",",pop[[i]]@sex,",",pop[[i]]@ARS,",",pop[[i]]@age,",",pop[[i]]@hetero,",",pop[[i]]@pID[1],",",pop[[i]]@pID[2],",",1,file=filename,append=TRUE)
  }
  for(i in DEAD){
    cat("\n",YR,",",pop[[i]]@ID,",",pop[[i]]@size,",",pop[[i]]@bvs,",",pop[[i]]@FixRepro,",",pop[[i]]@bvr,",",pop[[i]]@food,",",pop[[i]]@sex,",",pop[[i]]@ARS,",",pop[[i]]@age,",",pop[[i]]@hetero,",",pop[[i]]@pID[1],",",pop[[i]]@pID[2],",",0,file=filename,append=TRUE)
  }
}
 return(filename)
}
#end of LeprechaunSimul

# THE COMMAND FOR RUNNING THE CODE: LeprechaunSimul(StudyLength=50,DIR="Data")


