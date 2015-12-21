### Set your own working directory
### ------------------------------
#setwd("/Users/koen/Documents/R_scripts/IBM/IBM/") # koen
setwd("/Users/gauravbaruha/Desktop/PhD Research/Zurich PhD Research/IBM/GITHUB_IBM/IBM") # gaurav

### Run

# 1. Set seed
set.seed(12)

# 2. Set genotype map and save
source("genetics.R")

gvaluesZ<-FunctionGvalues(nbLoci = 10,nbAlleles = 10,dominance = 1,overdominance = 0,SDeffects =1,SDalleles = 1)
gvaluesR<-FunctionGvalues(nbLoci = 10,nbAlleles = 10,dominance = 1,overdominance = 0,SDeffects = 1,SDalleles = 2)

save(gvaluesZ=gvaluesZ,gvaluesR=gvaluesR,file = "GeneticMatrices")

# 3. Load in class defintions
source("ibm.R")

# 4. Run the simulation
LeprechaunSimul(DIR= "Data", StudyLength=5)

# Other parameters:
# SurvivalSelection1=0.5,SurvivalSelection2=0,fertilitySelection1=0.1,fertilitySelection2=0,
# camembertSelection=0.1,SurvivalPenaltyForRepro=0.2,MeanCamembert=4000,SDCamembert=10,CamembertTrend=-5,
# MeanBirthSize=10,MeanGrowth=1.245,SexualMaturity=1, replica=1,
# PlasticityBirthSize=1,PlasticityReproduction=0,MaternalEffect=1,MeanRepro=0.1,
# StudyLength=2000,InitialPopSize=70,sizeSurvivalInflection=0,CamembertthresholdGrowth=30,GeneticMatrices="GeneticMatrices")