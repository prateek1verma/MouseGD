rm(list = ls())
gc()
save.image()  # overwrites .RData with empty workspace
# 
# 
# setwd("~/Documents/GitHub/MouseGD/MGDrivEmouse2")
# pkg <- "MGDrivEmouse2"
# remove.packages("MGDrivEmouse2")
# # Remove from all libraries on your machine
# for (lib in .libPaths()) {
#   if (pkg %in% rownames(installed.packages(lib = lib))) {
#     message("Removing from: ", lib)
#     try(remove.packages(pkg, lib = lib), silent = TRUE)
#   }
# }
# 
# pkgbuild::clean_dll("~/Documents/GitHub/MouseGD/MGDrivEmouse2")
# system("rm -rf /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/MGDrivEmouse2/R/MGDrivEmouse2.rdb")
# system("rm -rf /Library/Frameworks/R.framework/Versions/4.5-arm64/Resources/library/MGDrivEmouse2")
# 
# setwd("~/Documents/GitHub/MouseGD/MGDrivEmouse2")
# devtools::document("~/Documents/GitHub/MouseGD/MGDrivEmouse2")
# devtools::install("~/Documents/GitHub/MouseGD/MGDrivEmouse2")

# Then load it
library(MGDrivEmouse2)
library(BBmisc)


####################
# Output Folder
####################
outFolder <- "outFolder"
dir.create(path = outFolder)

# # Check if the folder exists before deleting
# if (dir.exists(outFolder)) {
#   unlink(outFolder, recursive = TRUE, force = TRUE)
#   cat("Folder deleted:", outFolder, "\n")
# } else {
#   cat("Folder does not exist:", outFolder, "\n")
# }

####################
# Simulation Parameters
####################
# days to run the simulation
tMax <- 7475

# number of Monte Carlo iterations
nRep <- 2


# each Monte Carlo iteration gets its own folder
folderNames <- file.path(outFolder,
                         formatC(x = 1:nRep, width = 3, format = "d", flag = "0"))

# biological parameters
# bioParameters <- list(betaK= 6, litters = 7.5, tGest=19, tNursing=23, tAdo=37, muAI = (1/690), muJI = 0, muN = 0, muG = 0, theta = 22.4)

# biological parameters (including density-independent mortality rates)
# bioParameters <- list(betaK  = 6,
#                       litters = 7.5,
#                       tGest = 19,
#                       tNursing = 23,
#                       tAdo = 37,
#                       muAd = 1/690,
#                       muAI = 1/690,
#                       muJI = 0.005, # Assuming 85% survival over 37 days
#                       muN  = 0.012, # Assuming 75% survival over 23 days
#                       muG  = 0.0028, # Assuming 95% survival over 19 days
#                      theta = 22.4)

bioParameters <- list(betaK  = 6,
                      litters = 7.5,
                      tGest = 19,
                      tNursing = 23,
                      tAdo = 37,
                      muAd = 1/690,
                      muAI = 1/690,
                      muJI = 0.00, # Assuming 85% survival over 37 days
                      muN  = 0.00, # Assuming 75% survival over 23 days
                      muG  = 0.00, # Assuming 95% survival over 19 days
                      theta = 22.4)

sitesNumber <- 2 # number of patches
cap <- c(10000,40000) # adult carrying capacity for small patch and big patch respectively

# auxiliary function
triDiag <- function(upper, lower){
  
  # return matrix
  retMat <- matrix(data = 0, nrow = length(upper) + 1, ncol = length(upper) + 1)
  
  # set index values for upper/lower triangles
  indx <- 1:length(upper)
  
  # set forward/backward migration using matrix access
  retMat[cbind(indx+1,indx)] <- lower
  retMat[cbind(indx,indx+1)] <- upper
  
  # set stay probs
  diag(x = retMat) <- 1-rowSums(x = retMat)
  
  return(retMat)
}

# fill movement matrix
# Remember, rows need to sum to 1.
moveMat <- triDiag(upper = rep.int(x = 0.2, times = sitesNumber-1),
                   lower = rep.int(x = 0.05, times = sitesNumber-1))

# batch migration is disabled by setting the probability to 0
batchMigration <- basicBatchMigration(batchProbs=0,
                                      sexProbs=c(.5,.5),
                                      numPatches=sitesNumber)

####################
# Basic Inheritance pattern
####################
# CRISPR-Cas9 Sox9 Homing Gene Drive
# devtools::source_url(
#   "https://raw.githubusercontent.com/eabrown2378/MGDrivE/master/MGDrivE/R/SoxHomingDrive.R"
# )
cube <- SoxHomingDrive(cM = 0.999, cF = 0.999, chM = 0.99, chF = 0.99, crM = (1/3), crF = (1/3))


cube$phi <- c(0,0,0,0,0,0,0,0,0,0, # 'mWW', 'mWH', 'mWR', 'mWB', 'mHH', 'mHR', 'mHB', 'mRR', 'mRB', 'mBB',
              1,0,1,1,0,0,0,1,1,1) # 'fWW', 'fWH', 'fWR', 'fWB', 'fHH', 'fHR', 'fHB', 'fRR', 'fRB', 'fBB'
# adjusts for cube using sex-specific genotypes
# and specifies that females with 'H' allele develop with male morphology

cube$s <- c(1,1,1,0.1,1,1,0.1,1,0.1,0, #'mWW', 'mWH', 'mWR', 'mWB', 'mHH', 'mHR', 'mHB', 'mRR', 'mRB', 'mBB',
            1,0,1,0.1,0,0,0,1,0.1,0)   #'fWW', 'fWH', 'fWR', 'fWB', 'fHH', 'fHR', 'fHB', 'fRR', 'fRB', 'fBB'
# fitness reduction assumed for each genotype


####################
# Setup releases and batch migration
####################
# set up the empty release vector
#  MGDrivE pulls things out by name
patchReleases <- replicate(n=sitesNumber,
                           expr={list(maleReleases=NULL,femaleReleases=NULL,
                                      eggReleases=NULL,matedFemaleReleases=NULL)},
                           simplify=FALSE)

# choose release parameters

releasesParameters_SP <- list(releasesStart=2000,
                              releasesNumber=50,
                              releasesInterval=30,
                              releaseProportion=500)

releasesParameters_BP <- list(releasesStart=2000,
                              releasesNumber=50,
                              releasesInterval=30,
                              releaseProportion=500)

# specify toxicological parameters

exposure  = 1.26 # estimated average dose of toxicant in mice
expTime = 7 # length of toxicity experiment used to generate exposure-response curve

# daily average toxicant-induced mortality in mice
# calculated from exposure-response equation and the
# length of the corresponding toxicity test:

toxCurve = (1/(1+exp(-3.2*(log(exposure)-log(1.774)))))/expTime 

# the time interval at which the toxicant is
# present in the system, assuming no chemical degradation,
# assuming chemical CAN be removed:

toxTime = c((releasesParameters_BP$releasesStart
             #+ releasesParameters_BP$releasesNumber*releasesParameters_BP$releasesInterval
),(releasesParameters_BP$releasesStart +
     #releasesParameters_BP$releasesNumber*releasesParameters_BP$releasesInterval
     + 56)) 

# generate small patch release vector
ReleasesVector_SP <- generateReleaseVector(driveCube=cube,
                                           releasesParameters=releasesParameters_SP)

# generate big patch release vector
ReleasesVector_BP <- generateReleaseVector(driveCube=cube,
                                           releasesParameters=releasesParameters_BP)

# put releases into the proper place in the release list
#  This specifies the releases for the small patch (SP) and big patch (BP) respectively 
patchReleases[[1]]$maleReleases <- ReleasesVector_SP


patchReleases[[2]]$maleReleases <- ReleasesVector_BP



####################
# Combine parameters and run!
####################
# setup parameters for the network. This builds a list of parameters required for
# every population in the network.
netPar <- parameterizeMGDrivE(runID = 1,
                              simTime = tMax,
                              nPatch = sitesNumber,
                              beta = bioParameters$betaK,
                              litters = bioParameters$litters,
                              tGest = bioParameters$tGest,
                              tNursing = bioParameters$tNursing,
                              tAdo = bioParameters$tAdo,
                              muAd = bioParameters$muAd,
                              muAI = bioParameters$muAI,
                              muJI = bioParameters$muJI,
                              muN = bioParameters$muN,
                              muG = bioParameters$muG,
                              k = cap,
                              theta = bioParameters$theta,
                              inheritanceCube = cube,
                              AdPopRatio_M = matrix(c(1,1,0,0),2,2, dimnames = list(NULL,c("mWW","fWW"))),
                              AdPopRatio_F = matrix(c(0,0,1,1),2,2, dimnames = list(NULL,c("mWW","fWW"))))

# set MGDrivE to run stochastic
setupMGDrivE(stochasticityON = TRUE, verbose = FALSE)

# build network prior to run
MGDrivESim <- Network$new(params=netPar,
                          driveCube=cube,
                          patchReleases=patchReleases,
                          migrationMale=moveMat,
                          migrationFemale=moveMat,
                          migrationBatch=batchMigration,
                          directory=folderNames,
                          verbose = TRUE,
                          toxMort = toxCurve,
                          toxInt = toxTime)
# run simulation
MGDrivESim$multRun(verbose = TRUE)


####################
# Post Analysis
####################
# First, split output by patch
# Second, aggregate females by their mate choice
for(i in 1:nRep){
  splitOutput(readDir = folderNames[i], remFile = TRUE, verbose = FALSE)
  aggregateFemales(readDir = folderNames[i], genotypes = cube$genotypesID,
                   remFile = TRUE, verbose = TRUE)
}

# plot output of first run to see effect
plotMGDrivESingle(readDir=folderNames[1],totalPop = TRUE,lwd=3.5,alpha=1)

# # plot all repetitions together
# png(paste(outFolder, ".png", sep = ""), width = 1000, height = 750)
# plotMGDrivEMult(readDir=outFolder,lwd=0.35,alpha=0.75)
# dev.off()
# 
# out_path = as.character(paste(getwd(),outFolder, sep = "/"))
# 
# 
# timeSeq = seq(2365,tMax, by = 365) # select time points each year, assuming no leap years, starting
# # 5 years after initial genedrive/pesticide deployment
# 
# ext_table = MGDcpt(tPoints = timeSeq) # returns [[1]] a dataframe for the percent of simulations
# # where each genotype (and the total population) reached zero
# # and [[2]] a list of data frames for the total abundance and
# # genotype abundances for each simulation (1 frame per simulation)
# 
# ext_table[[1]] # percent of simulations that went to zero (for conditional probability table in Bayes Net)
# 
# write.csv(x = ext_table[[1]], file = paste(outFolder, "_zeroProbs", ".csv", sep = "")) #write to csv
# 
# totCol = ncol(ext_table[[1]]) # column number corresponding to the total population of mice (not genotype-specific)
# 
# zero = ext_table[[1]][,totCol]
# names(zero) = NULL
# 
# sparse = (abs((Reduce(`+`, lapply(ext_table[[2]], `>`, 0))/length(ext_table[[2]])*100)
#               - (Reduce(`+`, lapply(ext_table[[2]], `>`, 1000))/length(ext_table[[2]])*100)))[,totCol]
# names(sparse) = NULL
# 
# low = (abs((Reduce(`+`, lapply(ext_table[[2]], `>`, 1000))/length(ext_table[[2]])*100)
#            - (Reduce(`+`, lapply(ext_table[[2]], `>`, 10000))/length(ext_table[[2]])*100)))[,totCol]
# names(low) = NULL
# 
# med = (abs((Reduce(`+`, lapply(ext_table[[2]], `>`, 10000))/length(ext_table[[2]])*100)
#            - (Reduce(`+`, lapply(ext_table[[2]], `>`, 25000))/length(ext_table[[2]])*100)))[,totCol]
# names(med) = NULL
# 
# high = (Reduce(`+`, lapply(ext_table[[2]], `>`, 25000))/length(ext_table[[2]])*100)[,totCol]
# names(high) = NULL
# 
# # create table of population probability distribution (rows = year; columns = percent of simulations
# # within population threshold)
# totDistribution = data.frame(zero, sparse, low, med, high)
# 
# # write to csv (use this file as part of the conditional probability table in the Bayesian network)
# write.csv(x = totDistribution, file = paste(outFolder, "_popDist", ".csv", sep = ""))
# 
