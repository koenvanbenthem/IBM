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
    food = "integer",
    sex = "character",
    DNAZ = "matrix",
    DNAR = "matrix",
    bvs = "numeric",
    bvr = "numeric",
    ARS = "integer"
  )
)

setGeneric("Food",function(Object,food){standardGeneric("Food")})
setGeneric("Num_off",function(Object){standardGeneric("Num_off")})
setGeneric("Age",function(Object){standardGeneric("Age")})
setGeneric("Grow",function(Object){standardGeneric("Grow")})
setGeneric("Surv",function(Object){standardGeneric("Surv")})

setGeneric("IDretrieve",function(Object){standardGeneric("IDretrieve")})
setGeneric("Size",function(Object){standardGeneric("Size")})
setGeneric("FixRepro",function(Object){standardGeneric("FixRepro")})
setGeneric("Sex",function(Object){standardGeneric("Sex")})

############### Definition of the basic methods (for printing to the screen and initialisation)######
###############################################################################################%
setMethod("show","Leprechaun",
          function(object){
            cat(object@ID,"\t",object@size,"\t",object@food,"\t",object@age,"\t",object@sex,"\t(",object@pID[1],",",object@pID[2],")\t",object@Birth,"\t",object@alive,"\n",sep="")
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
  
  .Object@food<-as.integer(0)
  
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

# Incorporate the effect of size to the survival function
sizeSurvival<-function(age,size,food,ARS){
  #sizedeviation<-size-(MeanBirthSize*prod(Growth[1:(age+1)]))
  p<-bathtub(age)
  if(p<1)#because size does not prevent animals of maximal age to die out
  {
    plogit<-log(p/(1-p))
    Philogit<-plogit-SurvivalSelection1*(size-SizeSurvivalInflection)-SurvivalSelection2*(size-SizeSurvivalInflection)^2+ARS*SurvivalPenaltyForRepro
    p<-exp(Philogit)/(1+exp(Philogit))
    if(food<100)
    {
      p<-1-(1-p)*(food/100)^(1/2)
    }
  }
  return(1-p)
}

# Applying the bathtub in a surival function
setMethod("Surv","Leprechaun",function(Object){
  
  if(runif(1)>sizeSurvival(Object@age,Object@size,Object@food,Object@ARS)){
    Object@alive<-FALSE
    DEAD<<-c(DEAD,Object@ID)
  }
  
  return(Object)
  
})

# Simple function, simply adds 1 to the age
setMethod("Age","Leprechaun",function(Object){
  Object@age<-as.integer(Object@age+1)
  return(Object)
})

# Growth // Sizes change, this fact is known by many - if not all
setMethod("Grow","Leprechaun",function(Object){
  MeanGrowthAge<-(MeanGrowth+Object@age-1)/Object@age #maximal growth converge to 1 with increasing age
  SdGrowthAge<-(MeanGrowthAge-1)/2
  
  GrowthC<- 1+ (MeanGrowthAge-1)*(2/(1+exp(-0.1*(Object@food-FoodThresholdGrowth)))-1)#gives the camemebert specific expected growth
  
  Object@size<-Object@size*rnorm(1,mean=GrowthC,sd = SdGrowthAge)# but the animal can always shrink
  
  if(Object@size > GrowthC )
    
    return(Object)
})

# Retrieving the ID of an individual
setMethod("IDretrieve","Leprechaun",function(Object){
  return(Object@ID)
})

# Retrieving the size of an individual
setMethod("Size","Leprechaun",function(Object){
  return(Object@size)
})

# Retrieving the reproductive quality of an individualsetGeneric("Surv",function(Object){standardGeneric("Surv")})
setMethod("FixRepro","Leprechaun",function(Object){
  return(Object@FixRepro)
})

# Retrieving the sex of an individual
setMethod("Sex","Leprechaun",function(Object){
  return(Object@sex)
})

# Calculating the number of offspring for a females
setMethod("Num_off","Leprechaun",function(Object){  
  lambda<-exp(log(MeanRepro)+Object@FixRepro+FertilitySelection1*(Object@size-MeanBirthSize)/10+FertilitySelection2*((Object@size-MeanBirthSize)/10)^2
              +FoodSelection*((Object@food)^(1/3))/10)
  repro<-rpois(n=1,lambda=lambda)+1
  
  logitV<-Object@age-SexualMaturity
  p<-1/(1+exp(-logitV))
  Object@ARS<-as.integer(repro)*rbinom(1,size = 1,prob = p)#*disaster[YR]
  return(Object)
})


# Function for camembert attributions
setMethod("Food","Leprechaun",function(Object,food){
  Object@food<-as.integer(food)
  return(Object)
})