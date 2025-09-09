###############################################################################
#       _   __     __                      __
#      / | / /__  / /__      ______  _____/ /__
#     /  |/ / _ \/ __/ | /| / / __ \/ ___/ //_/
#    / /|  /  __/ /_ | |/ |/ / /_/ / /  / ,<
#   /_/ |_/\___/\__/ |__/|__/\____/_/  /_/|_|
#
#   parameterizeMGDrivE
#   Original Code by: Marshall Lab
#   jared_bennett@berkeley.edu
#   December 2019
#   MODIFIED BY: ETHAN A. BROWN (JUL 12 2020)
#   ebrown23@nd.edu
###############################################################################

#' parameterizeMGDrivE
#'
#' Generate parameters for simulation on a \code{\link{Network}}.
#' Parameters include: average generation time \eqn{g}, population growth rate \eqn{R_{m}},
#' juvenile mortality \eqn{\mu_{Ju}}, and juvenile survival \eqn{\theta_{Ju}}, which
#' are shared between patches and calculated by \code{\link{calcAverageGenerationTime}},
#' \code{\link{calcPopulationGrowthRate}}, and \code{\link{calcLarvalStageMortalityRate}}. \cr
#' Patch-specific parameters \eqn{\alpha} and \eqn{L_{eq}}
#' are calculated for each patch by \code{\link{calcDensityDependentDeathRate}}
#' and \code{\link{calcLarvalPopEquilibrium}}.
#'
#' @param runID Begin counting runs with this set of parameters from this value
#' @param nPatch Number of \code{\link{Patch}}
#' @param simTime Maximum time to run simulation
#' @param sampTime Times to sample, used as tNow %% sampTime, default is every day
#' @param moveVar Variance of stochastic movement (not used in diffusion model of migration).
#' It affects the concentration of probability in the Dirchlet simplex, small values
#' lead to high variance and large values lead to low variance.
#' @param tGest Length of gestation stage
#' @param tNursing Length of nursing stage
#' @param tAdo Length of adolescent stage
#' @param beta Female litter size of wild-type
#' @param muAd Wild-type daily adult mortality (1/muAd is average wild-type lifespan)
#' @param muAI Daily adult mortality rate used for derived quantities
#' @param muJI Juvenile mortality rate
#' @param muN Nursing-stage mortality rate
#' @param muG Gestation-stage mortality rate
#' @param k Single number or vector of adult population carrying capacity at equilibrium
#' (single number implies all patches have the same population, but populations for
#' individual patches can be specified with a vector)
#' @param theta Shape parameter for mediating the carrying capacity for adult mice
#' @param AdPopRatio_F May be empty; if not, a vector gives the wildtype gene frequencies
#' among adult females at the beginning of simulation or a matrix provides different
#' initial frequencies for each patch (every row is a different patch, must have nrow = nPatch)
#' @param AdPopRatio_M May be empty; if not, a vector gives the wildtype gene frequencies
#' among adult males at the beginning of simulation or a matrix provides different
#' initial frequencies for each patch (every row is a different patch, must have nrow = nPatch)
#' @param inheritanceCube Inheritance cube to check/set population ratios at the beginning of the simulation
#' @param litters average number of litters per female mouse per year
#'
#' @examples
#' # using default parameters for 2 patches
#' #  using different population sizes for patches
#' simPars <- parameterizeMGDrivE(nPatch = 2, simTime = 365,
#'                                k = c(100,200), inheritanceCube = cubeMendelian())
#'
#' @export
parameterizeMGDrivE <- function(
  runID = 1L,
  nPatch,
  simTime,
  sampTime = 1L,
  moveVar = 1000L,
  tGest = 19L,
  tNursing = 23L,
  tAdo = 37L,
  beta = 6,
  litters = 7.5,
  muAd = I(1/690),
  muAI = I(1/690),
  muJI = 0,
  muN = 0,
  muG = 0,
  k,
  theta = 22.4,
  AdPopRatio_F,
  AdPopRatio_M,
  inheritanceCube
){

  # check required parameters


  if(any(missing(nPatch),missing(simTime),missing(k),missing(inheritanceCube))){
    stop("nPatch, simTime, k, and inheritanceCube must be provided by the user.")
  }

  # make empty parameter list
  pars = list()

  # fill list
  pars$nPatch = nPatch
  pars$simTime = simTime
  pars$sampTime = sampTime
  pars$moveVar = moveVar
  pars$runID = runID

  # biological parameters
  pars$timeJu = c("G"=tGest, "N"=tNursing, "A"=tAdo)
  pars$timeAd = 1/muAI-sum(pars$timeJu)
  pars$beta = beta # assumes average of 7.5 pups per litter. The decimal is appropriate because it is used to
  # calculate a rate parameter (lambda) in a Poisson distribution, so the number of pups will always be an integer
                        #
  pars$litters = litters # assumes female mice have a median number of litters per year (assume no leap year)
                               # that determines the probability that each female will birth a litter on a given day

  # initial parameters
  pars$muAd = muAd
  pars$muAI = muAI
  pars$muJI = muJI
  pars$muN = muN
  pars$muG = muG
  pars$theta = theta
  pars$k = k
  if(length(pars$k) == 1){
    pars$k = rep.int(x = pars$k, times = nPatch)
  } else if(length(pars$k)!=nPatch){
    stop("length of k vector must be 1 or nPatch (number of patches)")
  }



  # derived parameters
  pars$g = calcAverageGenerationTime(pars$timeJu,pars$timeAd)

  pars$thetaAd = 1-pars$muAI


  # setup female initial pop ratio
  if(missing(AdPopRatio_F)){
    # default behaviour - this way nothing needs to be specified

    # now, two options. This is not sex based, so only 1 wild-type
    #  or, sex based cube, we have an X and a Y, females are not Y
    if(length(inheritanceCube$wildType) == 1){
      # only 1 wild-type, so no problem
      pars$AdPopRatio_F <- matrix(data = 1, nrow = nPatch, ncol = length(inheritanceCube$wildType),
                                  dimnames = list(NULL,inheritanceCube$wildType))
    } else if(length(inheritanceCube$wildType) == 2){
      # one female and one male wild-type.
      # get index of female wild-type
      whichGeno <- grep(pattern = "Y", x = inheritanceCube$wildType, fixed = TRUE, invert = TRUE)
      # setup default matrix
      pars$AdPopRatio_F <- matrix(data = 0, nrow = nPatch, ncol = length(inheritanceCube$wildType),
                                  dimnames = list(NULL,inheritanceCube$wildType))
      # set all to female genotype
      pars$AdPopRatio_F[ ,whichGeno] <- 1

    } else {
      stop("Default AdPopRatio_F only handles 1 or 2 wild-type genotypes.\n
           Please provide a matrix specifying female initial genotypes for every patch.")
    }

  } else if(is.null(dim(AdPopRatio_F)) ){
    # behaviour of user supplied 1 patch worth of weights.
    # has to be supplied as a vector, the matrix stuff makes this way too difficult
    if(abs(sum(AdPopRatio_F) - 1) > sqrt(.Machine$double.eps) ) stop('AdPopRatio_F must sum to 1')
    # check that columns names are in inheritance cube
    if(is.null(names(AdPopRatio_F)) || !all(names(AdPopRatio_F) %in% inheritanceCube$genotypesID)) {
      stop("Names for AdPopRatio_F must be specified as one of the genotypesID in the inheritance cube.")
    }
    # set all patches equal
    pars$AdPopRatio_F <- matrix(data = AdPopRatio_F, nrow = nPatch, ncol = length(AdPopRatio_F),
                                byrow = TRUE, dimnames = list(NULL,names(AdPopRatio_F)) )

  } else if(dim(AdPopRatio_F)[1] == nPatch){
    # behaviour if user supplies a matrix of probabilities.
    # each row is a different patch
    # check that all patches sum to 1
    if(any(abs(rowSums(AdPopRatio_F) - 1) > sqrt(.Machine$double.eps)) ) stop('Each row of AdPopRatio_F must sum to 1')
    # check that columns names are in inheritance cube
    if(is.null(colnames(AdPopRatio_F)) || !all(colnames(AdPopRatio_F) %in% inheritanceCube$genotypesID)) {
      stop("Column names for AdPopRatio_F must be specified as one of the genotypesID in the inheritance cube.")
    }
    # store
    pars$AdPopRatio_F <- AdPopRatio_F

  } else {
    stop("AdPopRatio_F has been miss specified.\n
         Left blank - default, all populations are the same and begin as wild-type individuals\n
         Vector - a named vector that sums to one. All populations will be the same\n
         Matrix - an nPatch by nGenotype matrix with column names and all rows sum to 1.
         Specifies each population individually.")
  }


  # setup male initial pop ratio
  if(missing(AdPopRatio_M)){
    # default behaviour - this way nothing needs to be specified

    # now, two options. This is not sex based, so only 1 wild-type
    #  or, sex based cube, we have an X and a Y, females are not Y
    if(length(inheritanceCube$wildType) == 1){
      # only 1 wild-type, so no problem
      pars$AdPopRatio_M <- matrix(data = 1, nrow = nPatch, ncol = length(inheritanceCube$wildType),
                                  dimnames = list(NULL,inheritanceCube$wildType))
    } else if(length(inheritanceCube$wildType) == 2){
      # one female and one male wild-type.
      # get index of female wild-type
      whichGeno <- grep(pattern = "Y", x = inheritanceCube$wildType, fixed = TRUE)
      # setup default matrix
      pars$AdPopRatio_M <- matrix(data = 1, nrow = nPatch, ncol = length(inheritanceCube$wildType),
                                  dimnames = list(NULL,inheritanceCube$wildType))
      # set all to female genotype
      pars$AdPopRatio_M[ ,whichGeno] <- 1

    } else {
      stop("Default AdPopRatio_M only handles 1 or 2 wild-type genotypes.\n
           Please provide a matrix specifying female initial genotypes for every patch.")
    }

  } else if(is.null(dim(AdPopRatio_M)) ){
    # behaviour of user supplied 1 patch worth of weights.
    # has to be supplied as a vector, the matrix stuff makes this way too difficult
    if(abs(sum(AdPopRatio_M) - 1) > sqrt(.Machine$double.eps) ) stop('AdPopRatio_M must sum to 1')
    # check that columns names are in inheritance cube
    if(is.null(names(AdPopRatio_M)) || !all(names(AdPopRatio_M) %in% inheritanceCube$genotypesID)) {
      stop("Names for AdPopRatio_M must be specified as one of the genotypesID in the inheritance cube.")
    }
    # set all patches equal
    pars$AdPopRatio_M <- matrix(data = AdPopRatio_M, nrow = nPatch, ncol = length(AdPopRatio_M),
                                byrow = TRUE, dimnames = list(NULL,names(AdPopRatio_M)) )

  } else if(dim(AdPopRatio_M)[1] == nPatch){
    # behaviour if user supplies a matrix of probabilities.
    # each row is a different patch
    # check that all patches sum to 1
    if(any(abs(rowSums(AdPopRatio_M) - 1) > sqrt(.Machine$double.eps)) ) stop('Each row of AdPopRatio_M must sum to 1')
    # check that columns names are in inheritance cube
    if(is.null(colnames(AdPopRatio_M)) || !all(colnames(AdPopRatio_M) %in% inheritanceCube$genotypesID)) {
      stop("Column names for AdPopRatio_M must be specified as one of the genotypesID in the inheritance cube.")
    }
    # store
    pars$AdPopRatio_M <- AdPopRatio_M

  } else {
    stop("AdPopRatio_M has been miss specified.\n
         Left blank - default, all populations are the same and begin as wild-type individuals\n
         Vector - a named vector that sums to one. All populations will be the same\n
         Matrix - an nPatch by nGenotype matrix with column names and all rows sum to 1.
         Specifies each population individually.")
  }

  # check the list
  invisible(Map(f = check, pars))


  # if pass the check, return the parameter vector
  return(pars)
}

########################################################################
# Equations and Equilibrium Parameters for parameterizeMGDrivE()
########################################################################

# check for positive parameter values
check <- function(x){
  if(is.numeric(x)||is.integer(x)){
    if(any(x < 0)){
      stop("only nonnegative parameter values allowed")
    }
  }
}


#' Calculate Average Generation Time
#'
#' Calculate \eqn{g}, average generation time, given by: \deqn{g=T_e+T_l+T_p+\frac{1}{\mu_{AI}}}
#'
#' @param stagesDuration Vector of lengths of juvenile stages, \eqn{T_{g}, T_{n}, T_{a}}
#' @param tAd Vector of lengths of juvenile stages, \eqn{T_{a}}
#'
calcAverageGenerationTime <- function(stagesDuration, tAd){
  return(sum(c(stagesDuration,tAd)))
}


