pkgname <- "MGDrivE"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
base::assign(".ExTimings", "MGDrivE-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('MGDrivE')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("Network")
### * Network

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: Network
### Title: Network Class Definition
### Aliases: Network
### Keywords: R6 class

### ** Examples

 ## Not run: 
##D  # There are no simple examples for this, so looking at the vignettes would be
##D  #  most useful.
##D 
##D  # Complete manual with examples, but none explored in depth.
##D  vignette("MGDrivE-Examples", package = "MGDrivE")
##D 
##D  # One example, explored in great detail. This is probably more helpful.
##D  vignette("MGDrivE-Run", package = "MGDrivE")
##D 
##D  
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("Network", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("aggregateFemales")
### * aggregateFemales

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: aggregateFemales
### Title: Aggregate Female Output by Genotype
### Aliases: aggregateFemales

### ** Examples

## Not run: 
##D # This example assumes user has already run MGDrivE and generated output.
##D #  This also assumes that the user has already split output by patch.
##D # See vignette for complete example.
##D 
##D # set read/write directory
##D fPath <- "path/to/data/containing/folder"
##D 
##D # Need genotypes from the cube run in the simulation
##D #  This is dependent on the simulation run
##D #  Using Mendelian cube for this example
##D cube <- cubeMendelian()
##D 
##D # no return value from function
##D aggregateFemales(readDir= fPath, writeDir = NULL, genotypes = cube$genotypesID,
##D                  remFile = TRUE)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("aggregateFemales", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("aggregateOutput")
### * aggregateOutput

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: aggregateOutput
### Title: Aggregate Output Over Landscape
### Aliases: aggregateOutput

### ** Examples

## Not run: 
##D # This assumes user has run MGDrivE and output is in fPath.
##D #  See vignette for examples on how to run MGDrivE
##D 
##D # read/write dirs
##D fPath <- "folder/containing/output"
##D oPath <- "folder/to/write/stuff"
##D 
##D # first, split output by patch and aggregate females by mate genotype
##D # remember, cube is for example and changes with simulation
##D #  landscape aggregation will work if females are not aggregated, but it's slower
##D cube <- cubeMendelian()
##D 
##D splitOutput(readDir = fPath, writeDir = NULL, remFile = TRUE)
##D aggregateFemales(readDir= fPath, writeDi = NULL, genotypes = cube$genotypesID,
##D                  remFile = TRUE)
##D 
##D # aggregate mosquitoes over entire landscape
##D #  no return value
##D aggregateOutput(readDir = fPath, writeDir = NULL)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("aggregateOutput", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("basicBatchMigration")
### * basicBatchMigration

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: basicBatchMigration
### Title: Make List of Batch Migration Parameters
### Aliases: basicBatchMigration

### ** Examples

# to setup for 3 patches
batchMigration = basicBatchMigration(batchProbs = 1e-5, sexProbs = c(0.1, 0.01), numPatches = 3)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("basicBatchMigration", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("basicRepeatedReleases")
### * basicRepeatedReleases

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: basicRepeatedReleases
### Title: Make List of Modified Mosquito Releases
### Aliases: basicRepeatedReleases

### ** Examples

## Not run: 
##D # Setup for 3 patches but only release in the first with a defined release
##D #  schedule, for the cube cubeHomingDrive:
##D 
##D patchReleases = replicate(n = 3, expr = {
##D   list(maleReleases = NULL, femaleReleases = NULL, eggReleases = NULL, matedFemaleReleases = NULL)
##D },simplify = FALSE)
##D 
##D patchReleases[[1]]$femaleReleases = MGDrivE::basicRepeatedReleases(releaseStart = 5,
##D                                                           releaseEnd = 30,
##D                                                           releaseInterval = 5,
##D                                                           releaseMatrix = matrix(c(5,100),1,2))
##D 
##D patchReleases[[1]]$maleReleases = MGDrivE::basicRepeatedReleases(releaseStart = 50,
##D                                                         releaseEnd = 60,
##D                                                         releaseInterval = 1,
##D                                                         releaseMatrix = matrix(c(5,100),1,2))
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("basicRepeatedReleases", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("calcCos")
### * calcCos

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: calcCos
### Title: Calculate Geodesic Distance - Cosine Method
### Aliases: calcCos

### ** Examples

# two-column matrix with latitude/longitude, in degrees
latLong = cbind(runif(n = 5, min = 0, max = 90),
                runif(n = 5, min = 0, max = 180))

# cosine distance formula
distMat = calcCos(latLongs = latLong)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("calcCos", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("calcExpKernel")
### * calcExpKernel

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: calcExpKernel
### Title: Calculate Exponential Stochastic Matrix
### Aliases: calcExpKernel

### ** Examples

# setup distance matrix
# two-column matrix with latitude/longitude, in degrees
latLong = cbind(runif(n = 5, min = 0, max = 90),
                runif(n = 5, min = 0, max = 180))

# Vincenty Ellipsoid  distance formula
distMat = calcVinEll(latLongs = latLong)

# calculate exponential distribution over distances
#  rate is just for example
kernMat = calcExpKernel(distMat = distMat, rate = 10)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("calcExpKernel", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("calcGammaKernel")
### * calcGammaKernel

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: calcGammaKernel
### Title: Calculate Gamma Stochastic Matrix
### Aliases: calcGammaKernel

### ** Examples

# setup distance matrix
# two-column matrix with latitude/longitude, in degrees
latLong = cbind(runif(n = 5, min = 0, max = 90),
                runif(n = 5, min = 0, max = 180))

# Vincenty Ellipsoid  distance formula
distMat = calcVinEll(latLongs = latLong)

# calculate gamma distribution over distances
#  shape and rate are just for example
kernMat = calcGammaKernel(distMat = distMat, shape = 1, rate = 1)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("calcGammaKernel", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("calcHaversine")
### * calcHaversine

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: calcHaversine
### Title: Calculate Geodesic Distance - Haversine Method
### Aliases: calcHaversine

### ** Examples

# two-column matrix with latitude/longitude, in degrees
latLong = cbind(runif(n = 5, min = 0, max = 90),
                runif(n = 5, min = 0, max = 180))

# Haversine distance formula
distMat = calcHaversine(latLongs = latLong)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("calcHaversine", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("calcHurdleExpKernel")
### * calcHurdleExpKernel

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: calcHurdleExpKernel
### Title: Calculate Zero-inflated Exponential Stochastic Matrix
### Aliases: calcHurdleExpKernel

### ** Examples

# setup distance matrix
# two-column matrix with latitude/longitude, in degrees
latLong = cbind(runif(n = 5, min = 0, max = 90),
                runif(n = 5, min = 0, max = 180))

# Vincenty Ellipsoid  distance formula
distMat = calcVinEll(latLongs = latLong)

# calculate hurdle exponential distribution over distances
#  rate and point mass are just for example
kernMat = calcHurdleExpKernel(distMat = distMat, rate = 1/1e6, p0 = 0.1)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("calcHurdleExpKernel", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("calcLognormalKernel")
### * calcLognormalKernel

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: calcLognormalKernel
### Title: Calculate Lognormal Stochastic Matrix
### Aliases: calcLognormalKernel

### ** Examples

# setup distance matrix
# two-column matrix with latitude/longitude, in degrees
latLong = cbind(runif(n = 5, min = 0, max = 90),
                runif(n = 5, min = 0, max = 180))

# Vincenty Ellipsoid  distance formula
distMat = calcVinEll(latLongs = latLong)

# calculate lognormal distribution over distances
#  mean and standard deviation are just for example
kernMat = calcLognormalKernel(distMat = distMat, meanlog = 100, sdlog = 10)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("calcLognormalKernel", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("calcOmega")
### * calcOmega

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: calcOmega
### Title: Solve for Omega (additional genotype-specific mortality)
### Aliases: calcOmega

### ** Examples

# reduce lifespan by 10%
#  Example mu is an average for Aedes
newOmega <- calcOmega(mu = 0.11, lifespanReduction = 0.90)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("calcOmega", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("calcQuantiles")
### * calcQuantiles

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: calcQuantiles
### Title: Summary Statistics for Stochastic MGDrivE
### Aliases: calcQuantiles

### ** Examples

## Not run: 
##D # This function assumes network$multRun() has been performed, or several
##D #  network$oneRun() have been performed and all of the data has been split
##D #  and aggregated.
##D 
##D # read/write paths
##D fPath <- "path/to/folder/ofFolders/with/data"
##D oPath <- "my/path/output"
##D 
##D # here, only calculate mean, no quantiles
##D #  no return value
##D calcQuantiles(readDir = fPath, writeDir = oPath, mean = TRUE,
##D               quantiles = NULL)
##D 
##D # here, calculate 2.5% and 97.5% quantiles
##D calcQuantiles(readDir = fPath, writeDir = oPath, mean = FALSE,
##D               quantiles = c(0.025, 0.975))
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("calcQuantiles", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("calcVinEll")
### * calcVinEll

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: calcVinEll
### Title: Calculate Geodesic Distance - Vincenty Ellipsoid Method
### Aliases: calcVinEll

### ** Examples

# two-column matrix with latitude/longitude, in degrees
latLong = cbind(runif(n = 5, min = 0, max = 90),
                runif(n = 5, min = 0, max = 180))

# Vincenty Ellipsoid  distance formula
distMat = calcVinEll(latLongs = latLong)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("calcVinEll", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("calcVinSph")
### * calcVinSph

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: calcVinSph
### Title: Calculate Geodesic Distance - Vincenty Sphere Method
### Aliases: calcVinSph

### ** Examples

# two-column matrix with latitude/longitude, in degrees
latLong = cbind(runif(n = 5, min = 0, max = 90),
                runif(n = 5, min = 0, max = 180))

# Vincenty Sphere  distance formula
distMat = calcVinSph(latLongs = latLong)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("calcVinSph", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("calcZeroInflation")
### * calcZeroInflation

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: calcZeroInflation
### Title: Calculates the zero-inflation part of a hurdle exponential
###   kernel.
### Aliases: calcZeroInflation

### ** Examples

# setup distance matrix
# two-column matrix with latitude/longitude, in degrees
latLong = cbind(runif(n = 5, min = 0, max = 90),
                runif(n = 5, min = 0, max = 180))

# Vincenty Ellipsoid  distance formula
distMat = calcVinEll(latLongs = latLong)

# get hurdle height
# Lets assume 80% stay probs and adult mortality of 0.1
hHeight <- calcZeroInflation(stayThroughLifespanProbability = 0.80,
                             adultMortality = 0.1)

# calculate hurdle exponential distribution over distances
kernMat = calcHurdleExpKernel(distMat = distMat, rate = 10, p0 = hHeight)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("calcZeroInflation", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("cube2csv")
### * cube2csv

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: cube2csv
### Title: Export a Cube to .csv
### Aliases: cube2csv

### ** Examples

## Not run: 
##D # output directory
##D oPath <- "path/to/write/output"
##D 
##D # setup inheritance cube for export, using Mendelian as the example
##D cube <- cubeMendelian()
##D 
##D # write out
##D cube2csv(cube = cube, directory = oPath, digits = 3)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("cube2csv", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("eraseDirectory")
### * eraseDirectory

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: eraseDirectory
### Title: Erase all files in a directory
### Aliases: eraseDirectory

### ** Examples

## Not run: 
##D # Path to directory, can tilde expand
##D myPath <- "~/path/to/write/output"
##D 
##D # Erase directory
##D #  No return value
##D eraseDirectory(directory = myPath)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("eraseDirectory", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("generateReleaseVector")
### * generateReleaseVector

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: generateReleaseVector
### Title: Make List of Modified Mosquito Releases
### Aliases: generateReleaseVector

### ** Examples

# setup a drive cube, using Mendelian as the example
cube <- cubeMendelian()

# setup release parameter list
#  releasesStart is the time of first release
#  releasesNumber is the number of releases
#  releasesInterval is the number of days between releases
#  releaseProportion is the number of mosquitoes released
relParams <- list(releasesStart = 25, releasesNumber = 1,
                  releasesInterval = 0, releaseProportion = 10)

# generate male releases
mRelVec <- generateReleaseVector(driveCube = cube,
                                 releasesParameters = relParams)

# generate mated female releases
fRelVec <- generateReleaseVector(driveCube = cube,
                                 releasesParameters = relParams,
                                 nameGenotypes = list(c("AA","AA", 10),
                                                      c("AA","aa", 10)))




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("generateReleaseVector", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("parameterizeMGDrivE")
### * parameterizeMGDrivE

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: parameterizeMGDrivE
### Title: parameterizeMGDrivE
### Aliases: parameterizeMGDrivE

### ** Examples

# using default parameters for 2 patches
#  using different population sizes for patches
simPars <- parameterizeMGDrivE(nPatch = 2, simTime = 365,
                               AdPopEQ = c(100,200), inheritanceCube = cubeMendelian())




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("parameterizeMGDrivE", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plotMGDrivEMult")
### * plotMGDrivEMult

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plotMGDrivEMult
### Title: Plot
### Aliases: plotMGDrivEMult

### ** Examples

## Not run: 
##D # Requires the user to have run MGDrivE, logically stochastic, analyzed
##D #  the data, and stored it in the directory shown below.
##D # See vignette for complete example
##D 
##D # Folder where single run is stored
##D fPath <- "path/to/data/containing/folder"
##D 
##D # plot output to see effect
##D plotMGDrivEMult(readDir=fPath,totalPop = TRUE,lwd=3.5,alpha=1)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plotMGDrivEMult", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("plotMGDrivESingle")
### * plotMGDrivESingle

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: plotMGDrivESingle
### Title: Plot
### Aliases: plotMGDrivESingle

### ** Examples

## Not run: 
##D # Requires the user to have run MGDrivE, deterministic or stochastic, analyzed
##D #  the data, and stored it in the directory shown below.
##D # See vignette for complete example
##D 
##D # Folder where single run is stored
##D fPath <- "path/to/data/containing/folder"
##D 
##D # plot output to see effect
##D plotMGDrivESingle(readDir=fPath,totalPop = TRUE,lwd=3.5,alpha=1)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("plotMGDrivESingle", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("retrieveOutput")
### * retrieveOutput

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: retrieveOutput
### Title: Retrieve Output
### Aliases: retrieveOutput

### ** Examples

## Not run: 
##D # Example assumes user has run and analyzed MGDrivE.
##D #  See vignette for examples of how to do that.
##D 
##D # set read directory
##D fPath <- "path/to/split/aggregated/output"
##D 
##D # read in data as nested lists
##D dataList <- retrieveOutput(readDir = fPath)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("retrieveOutput", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("setupMGDrivE")
### * setupMGDrivE

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: setupMGDrivE
### Title: Setup MGDrivE
### Aliases: setupMGDrivE

### ** Examples

# run deterministic MGDrivE
setupMGDrivE(stochasticityON = FALSE)

# run stochastic MGDrivE
setupMGDrivE(stochasticityON = TRUE)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("setupMGDrivE", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("splitOutput")
### * splitOutput

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: splitOutput
### Title: Split Output by Patch
### Aliases: splitOutput

### ** Examples

## Not run: 
##D # This example assumes user has already run MGDrivE and generated output.
##D #  If that's untree, see vignette for complete example
##D fPath <- "path/to/data/containing/folder"
##D oPath <- "path/to/write/output"
##D 
##D # split data by patch, keep original files
##D #  no return value
##D splitOutput(readDir = fPath, writeDir = oPath, remFile = FALSE)
##D 
##D # Alternatively, remove the original files and write new ones in their place
##D fPath <- "path/to/data/containing/folder"
##D 
##D splitOutput(readDir = fPath, writeDir = NULL, remFile = TRUE)
## End(Not run)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("splitOutput", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
