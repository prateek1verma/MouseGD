###############################################################################
#        ____        __       __
#       / __ \____ _/ /______/ /_
#      / /_/ / __ `/ __/ ___/ __ \
#     / ____/ /_/ / /_/ /__/ / / /
#    /_/    \__,_/\__/\___/_/ /_/
#
#   Patch Class MOUSE Population Simulation
#   Original Code by: Marshall Lab
#   jared_bennett@berkeley.edu
#   December 2019
#   MODIFIED BY: ETHAN A. BROWN (JUL 12 2020)
#   ebrown23@nd.edu
###############################################################################
###############################################################################
# Daily Simulation
###############################################################################

#' Daily Population Dynamics for a Patch
#'
#' Run population dynamics (including migration) for this patch. \cr
#' Performed in this order, see the following for each function: \cr
#' Adult Death: \code{\link{oneDay_adultDeath_deterministic_Patch}} or \code{\link{oneDay_adultDeath_stochastic_Patch}} \cr
#' Pupa Death/Maturation: \code{\link{oneDay_adoDM_deterministic_Patch}} or \code{\link{oneDay_adoDM_stochastic_Patch}} \cr
#' Larva Death/Maturation: \code{\link{oneDay_nursingDM_deterministic_Patch}} or \code{\link{oneDay_nursingDM_stochastic_Patch}} \cr
#' Gestating pup Death/Maturation: \code{\link{oneDay_gestDM_deterministic_Patch}} or \code{\link{oneDay_gestDM_stochastic_Patch}} \cr
#' Maturation: \code{\link{oneDay_maturation_deterministic_Patch}} or \code{\link{oneDay_maturation_stochastic_Patch}} \cr
#' Releases: \code{\link{oneDay_releases_Patch}} \cr
#' Mating: \code{\link{oneDay_mating_deterministic_Patch}} or \code{\link{oneDay_mating_stochastic_Patch}} \cr
#' Conceive pups: \code{\link{oneDay_conceive_deterministic_Patch}} or \code{\link{oneDay_conceive_stochastic_Patch}} \cr
#' Release Gestating Pups: \code{\link{oneDay_gestReleases_Patch}} \cr
#' Migration: \code{\link{oneDay_migrationOut_deterministic_Patch}} or \code{\link{oneDay_migrationOut_stochastic_Patch}} \cr
#'
oneDay_PopDynamics_Patch <- function(){

  ##################
  #Death/Maturation#
  ##################
  self$oneDay_adultD()
  self$oneDay_adoDM()
  self$oneDay_nursingDM()
  self$oneDay_gestDM()

  ##########
  #Maturation#
  ##########
  self$oneDay_maturation()

  ##########
  #Releases#
  ##########
  self$oneDay_releases()

  ########
  #Mating#
  ########
  self$oneDay_mating()

  ##########
  #Conceive Pups#
  ##########
  self$oneDay_layPups()
  self$oneDay_releasePups()

  ###########
  #Migration#
  ###########
  self$oneDay_migrationOut()

}


###############################################################################
# Death/Maturation
###############################################################################
##########
# Adult
##########
#' Deterministc Adult Survival
#'
#' Daily adult survival is calculated according to \deqn{\overline{\overline{Af_{[t-1]}}} * (1-\mu_{ad}) * \overline{\omega_{m/f}}},
#' where \eqn{\mu_{ad}} corresponds to adult mortality rate and \eqn{\overline{\omega_{m/f}}}
#' corresponds to genotype-specific male/female mortality effects.
#'
#' Toxcity is also included in this adult mortality process for a specified time interval
#' within the simulation. Daily survival is caluclated by by \eqn{1-\mu_{tox}}. At this time,
#' there is no parameter for genotype-specific toxicant susceptibility
#'

oneDay_adultDeath_deterministic_Patch <- function(){

#'
#'  density dependent function for adults mediated by carrying capacity
#'  function \eqn{k/{A_{t-1} + k}} applied to each genotype in the adult
#'  population vectors (male and female density-dependence calculated
#'  separately)
#'
  density <- ((private$NetworkPointer$get_k(ix = private$patchID)/2)/
               (private$NetworkPointer$get_k(ix = private$patchID)  + sum(private$popMale) + sum(private$popUnmated))) ^
    (1/private$NetworkPointer$get_theta())

  lifeM <- density * (1-private$NetworkPointer$get_muAd())

  # probability of survival
  probHolderM = private$NetworkPointer$get_omega() * lifeM

  # males that live through the day
  private$popMale[] = private$popMale * probHolderM



  lifeF <- density * (1-private$NetworkPointer$get_muAd())

  # probability of survival
  probHolderF = private$NetworkPointer$get_omega() * lifeF

  private$popUnmated[] = private$popUnmated * probHolderF

  # females that live through the day
  # this works becase female genotypes are rows, and R applies column-wise, thereby
  #  properly applying the vector to each female genotype in order
  private$popFemale[] = private$popUnmated %o% normalise(private$popMale) #make sure the female matrix keeps up

  tNow <- private$NetworkPointer$get_tNow()
  interval = sort(private$NetworkPointer$get_toxInt())
  dayToxSurv = 1 - private$NetworkPointer$get_toxMort()

  if (dayToxSurv < 1 & tNow >= interval[1] & tNow <= interval[2]){

    private$popMale[] = private$popMale * dayToxSurv
    private$popUnmated[] = private$popUnmated * dayToxSurv

    private$popFemale[] = private$popUnmated %o% normalise(private$popMale) #make sure the female matrix keeps up
  }


}

#' Stochastic Adult Survival
#'
#' Daily adult survival is sampled from a binomial distribution where survival
#' probability is given by \deqn{1-{\mu_{ad}} * \overline{\omega_m/f}}.
#' \eqn{\mu_{ad}} corresponds to adult mortality rate and \eqn{\overline{\omega_m/f}}
#' corresponds to genotype-specific mortality effects.
#'
#' Toxcity is also included in this adult mortality process for a specified time interval
#' within the simulation. Survival is sampled from a binomial distribution with survival
#' probability given by \eqn{1-\mu_{tox}}. At this time, there is no parameter for
#' genotype-specific toxicant susceptibility
#'


oneDay_adultDeath_stochastic_Patch <- function(){

  #'
  #'  density dependent function for adults mediated by carrying capacity
  #'  function \eqn{k/{A_{t-1} + k}} applied to each genotype in the adult
  #'  population vectors (male and female density-dependence calculated
  #'  separately)
  #'

  density <- ((private$NetworkPointer$get_k(ix = private$patchID)/2)/
                (private$NetworkPointer$get_k(ix = private$patchID) + sum(private$popMale) + (sum(private$popUnmated)))) ^
    (1/private$NetworkPointer$get_theta())

  lifeM <- (1-private$NetworkPointer$get_muAd()) * density

  # probability of survival
  probHolderM = private$NetworkPointer$get_omega() * lifeM

  # males that live through the day
  private$popMale[] <- rbinom(n = private$NetworkPointer$get_genotypesN(),
                              size = round(private$popMale),
                              prob = probHolderM)

  lifeF <- (1-private$NetworkPointer$get_muAd()) * density

  # probability of survival
  probHolderF = private$NetworkPointer$get_omega() * lifeF

  # females that live through the day
  # this also works because of how R applies things and fills matrices
  private$popUnmated[] = rbinom(n = private$NetworkPointer$get_genotypesN(),
                               size = round(private$popUnmated),
                               prob = probHolderF)

  private$popFemale[] = private$popUnmated %o% normalise(private$popMale) #make sure the female matrix keeps up

  #aggregates female genotype into their total abundance to sub for popUnmated instead of resetting after one day

  tNow <- private$NetworkPointer$get_tNow()
  interval = sort(private$NetworkPointer$get_toxInt())
  dayToxSurv = 1 - private$NetworkPointer$get_toxMort()

  if (dayToxSurv < 1 & tNow >= interval[1] & tNow <= interval[2]){

    private$popMale[] <- rbinom(n = private$NetworkPointer$get_genotypesN(),
                                size = round(private$popMale),
                                prob = dayToxSurv)

    private$popUnmated[] = rbinom(n = private$NetworkPointer$get_genotypesN(),
                                  size = round(private$popUnmated),
                                  prob = dayToxSurv)

    private$popFemale[] = private$popUnmated %o% normalise(private$popMale) #make sure the female matrix keeps up
  }

}


##########
# Pupa
##########
#' Deterministic Pupa Death and Maturation
#'
#' Daily adolescent survival is calculated according to \deqn{\overline{P_{[t-1]}} * {1-\mu_{JI}}},
#' where \eqn{\mu_{JI}} corresponds to adolescent mortality due to infection. \cr
#' See \code{\link{parameterizeMGDrivE}} for how these parameters are derived.
#'
oneDay_adoDM_deterministic_Patch <- function(){

  adoStart <- private$NetworkPointer$get_timeJu(stage = 'G') + private$NetworkPointer$get_timeJu(stage = 'N') + 1
  adoEnd <- private$NetworkPointer$get_timeJu()
  survAdo <- 1 - private$NetworkPointer$get_muJI()

  # Treat last day differently because the adolescent mice start to mature
  #  This does not handle the continuous to discrete time conversion artifact, see
  #  the maturation function for that.
  private$popHolder[] <- private$popJuvenile[ ,adoEnd] * survAdo

  # check if there are other days, then
  # run loop backwards to move populations as we go
  if((adoEnd - adoStart) > 0){
    for(i in (adoEnd-1):adoStart){
      private$popJuvenile[ ,i+1] = private$popJuvenile[ ,i] * survAdo
    } # end loop
  }

}

#' Stochastic Pupa Death and Maturation
#'
#' Daily adolescent survival is sampled from a binomial distribution, where survival
#' probability is given by \deqn{1-\mu_{JI}}. \eqn{\mu_{JI}} corresponds
#' to adolescent mortality due to infection. \cr
#' See \code{\link{parameterizeMGDrivE}} for how these parameters are derived.
#'
oneDay_adoDM_stochastic_Patch <- function(){

  # things to reuse
  nGeno <- private$NetworkPointer$get_genotypesN()
  adoStart <- private$NetworkPointer$get_timeJu(stage = 'G') + private$NetworkPointer$get_timeJu(stage = 'N') + 1
  adoEnd <- private$NetworkPointer$get_timeJu()
  survAdo <- 1 - private$NetworkPointer$get_muJI()

  # Treat last day differently because the adolescent mice start to mature
  #  This does not handle the continuous to discrete time conversion artifact, see
  #  the maturation function for that.
  private$popHolder[] <- rbinom(n = nGeno,
                                size = private$popJuvenile[ ,adoEnd],
                                prob = survAdo)

  # check if there are other days, then
  # run loop backwards to move populations as we go
  if((adoEnd - adoStart) > 0){
    for(i in (adoEnd-1):adoStart){
      private$popJuvenile[ ,i+1] <- rbinom(n = nGeno,
                                          size = private$popJuvenile[ ,i],
                                          prob = survAdo)
    } # end loop
  }

}


##########
# Larva
##########
#' Deterministic Larva Death and Maturation
#'
#' Calculate the number of nursing pups surviving from day to day, given by:
#' \deqn{\overline{L_{[t-1]}} * {1-\mu_{aq}}
#' See \code{\link{parameterizeMGDrivE}} for how these parameters are derived.
#' Maturation has no parameters, so the final day of nursing pups naturally enter the adolescent state.
#'
oneDay_nursingDM_deterministic_Patch <- function(){
  nGeno <- private$NetworkPointer$get_genotypesN()
  nursingStart <- private$NetworkPointer$get_timeJu(stage = 'G') + 1
  nursingEnd <- private$NetworkPointer$get_timeJu(stage = 'G') + private$NetworkPointer$get_timeJu(stage = 'N')


  # run loop backwards to move populations as we go
  for(i in nursingEnd:nursingStart){
    private$popJuvenile[ ,i+1] = private$popJuvenile[ ,i]
  } # end loop

}

#' Stochastic Larva Death and Maturation
#'
#' The daily number of nursing pups surviving is drawn from a binomial distribution, where
#' survival probability is given by \deqn{1-\mu_{aq}}
#' See \code{\link{parameterizeMGDrivE}} for how these parameters are derived.
#' Maturation has no parameters, so the final day of nursing pups naturally enter the adolescent state.
#'
oneDay_nursingDM_stochastic_Patch <- function(){

  # things to reuse
  nGeno <- private$NetworkPointer$get_genotypesN()
  nursingStart <- private$NetworkPointer$get_timeJu(stage = 'G') + 1
  nursingEnd <- private$NetworkPointer$get_timeJu(stage = 'G') + private$NetworkPointer$get_timeJu(stage = 'N')


  # run loop backwards to move populations as we go
  for(i in nursingEnd:nursingStart){
    private$popJuvenile[ ,i+1] = rbinom(n = nGeno,
                                       size = private$popJuvenile[ ,i],
                                       prob = 1)
  } # end loop

}


##########
# Gestating Pup
##########
#' Deterministic Gestating Pup Death and Maturation
#'
#' Daily gestating pup survival is calculated according to \deqn{\overline{E_{[t-1]}} * {1-\mu_{aq}}},
#' where \eqn{\mu_{aq}} corresponds to daily non-density-dependent juvenile mortality.
#' Gestating Pups transition into nursing pups at the end of \eqn{T_e}. \cr
#' See \code{\link{parameterizeMGDrivE}} for how these parameters are derived.
#'
oneDay_gestDM_deterministic_Patch <- function(){

  gestEnd <- private$NetworkPointer$get_timeJu(stage = 'G')

  # run loop backwards to move populations as we go
  for(i in gestEnd:1){
    private$popJuvenile[ ,i+1] = private$popJuvenile[ ,i]
  } # end loop

}

#' Stochastic Gestating Pup Death and Maturation
#'
#' Daily gestating pup survival is sampled from a binomial distribution, where survival
#' probability is given by \eqn{1-\mu_{aq}}. \eqn{\mu_{aq}} corresdponds
#' to daily non-density-dependent juvenile mortality. \cr
#' Gestating Pups transition into nursing pups at the end of \eqn{T_e}. \cr
#' See \code{\link{parameterizeMGDrivE}} for how these parameters are derived.
#'
oneDay_gestDM_stochastic_Patch <- function(){

  # things to reuse
  nGeno <- private$NetworkPointer$get_genotypesN()
  gestEnd <- private$NetworkPointer$get_timeJu(stage = 'G')

  # run loop backwards to move populations as we go
  for(i in gestEnd:1){
    private$popJuvenile[ ,i+1] = rbinom(n = nGeno,
                                       size = private$popJuvenile[ ,i],
                                       prob = 1)
  } # end loop

}


###############################################################################
# Maturation
###############################################################################

#' Deterministc Maturation
#'
#' Pupa first undergo one extra day of survival, calculated as \deqn{\overline{P_{[t-1]}} * {1-\mu_{ad}}}.
#' This is an artifact of the conversion from continuous to discrete time (as mentioned
#' in the original Hancock paper this model is derived from). \cr
#' Then, maturation into adult males is calculated as \deqn{{1-\overline{\phi}} * \overline{P_{[t]}}}
#' and into adult females as \deqn{\overline{\phi} * \overline{P_{[t]}}}
#'
oneDay_maturation_deterministic_Patch <- function(){

  # one extra death to match continuous time math
  #  This is an artifact of being discrete time
  private$popHolder[] = private$popHolder * (1-private$NetworkPointer$get_muAd())

  # get sex ratio for emergence
  phi = private$NetworkPointer$get_phi()

  # perform genotype-specific sex ratio and sex-dependent emergence
  private$popMale[] = private$popMale + private$popHolder * (1-phi) * private$NetworkPointer$get_xiM()

  private$popUnmated[] = private$popUnmated + private$popHolder * phi * private$NetworkPointer$get_xiF()

  private$popFemale[] = private$popUnmated %o% normalise(private$popMale) #make sure the female matrix keeps up

}

#' Stochastic Maturation
#'
#' Pupa first undergo one extra day of survival, calculated as a binomial over
#' \deqn{\overline{P_{[t-1]}} * {1-\mu_{ad}}}.
#' This is an artifact of the conversion from continuous to discrete time (as mentioned
#' in the original Hancock paper this model is derived from). \cr
#' Then, maturation is sampled from a binomial, where \eqn{{1-\overline{\phi}}} is
#' the genotype-specific probability of becoming male, and \eqn{\overline{\phi}}
#' is the genotype-specific of becoming female.
#'
oneDay_maturation_stochastic_Patch <- function(){

  # reuse
  nGeno = private$NetworkPointer$get_genotypesN()

  # there needs to be 1 extra death for discrete time
  private$popHolder[] = rbinom(n = nGeno,
                               size = private$popHolder,
                               prob =  1-private$NetworkPointer$get_muAd())

  # pull male/female distinction
  # only need to pull female, male is just popHolder - female
  private$popPupSex[] = rbinom(n = nGeno,
                             size = private$popHolder,
                             prob = private$NetworkPointer$get_phi() )

  # genotype-specific and sex-dependent emergence
  private$popMale[] = private$popMale + rbinom(n = nGeno,
                                               size = private$popHolder - private$popPupSex,
                                               prob = private$NetworkPointer$get_xiM() )
  private$popUnmated[] = private$popUnmated + rbinom(n = nGeno,
                                                     size = private$popPupSex,
                                                     prob = private$NetworkPointer$get_xiF() )

  private$popFemale[] = private$popUnmated %o% normalise(private$popMale) #make sure the female matrix keeps up

}


###############################################################################
# Releases
###############################################################################

#' Release Male/Female/Mated-Female Mosquitoes in a Patch
#'
#' Based on this patch's release schedule, \code{\link{generateReleaseVector}},
#' this function handles daily releases.
#'
oneDay_releases_Patch <- function(){

  # things to reuse
  tNow <- private$NetworkPointer$get_tNow()

  ##########
  # male
  ##########
  if((length(private$maleReleases) > 0) && (private$maleReleases[[1]]$tRelease <= tNow) ){
    # combine release
    # 1st column is index, 2nd column is amount
    idxRel <- private$maleReleases[[1]]$nRelease[ ,1]
    private$popMale[idxRel] = private$popMale[idxRel] + private$maleReleases[[1]]$nRelease[ ,2]
    # remove finished release
    private$maleReleases[[1]] = NULL
  }

  ##########
  # female
  ##########
  if((length(private$femaleReleases) > 0) && (private$femaleReleases[[1]]$tRelease <= tNow) ){
    # combine unmated females
    # 1st column is index, 2nd column is amount
    idxRel <- private$femaleReleases[[1]]$nRelease[ ,1]
    private$popUnmated[idxRel] = private$popUnmated[idxRel] + private$femaleReleases[[1]]$nRelease[ ,2]
    # remove finished release
    private$femaleReleases[[1]] = NULL
  }

  private$popFemale[] = private$popUnmated %o% normalise(private$popMale) #make sure the female matrix keeps up

  ##########
  # mated female
  ##########
  if((length(private$matedFemaleReleases) > 0) && (private$matedFemaleReleases[[1]]$tRelease <= tNow) ){
    # combine mated females
    # 1st column is female index, 2nd column is male index, 3rd column is amount
    idxRel <- private$matedFemaleReleases[[1]]$nRelease[ ,1:2,drop=FALSE]
    private$popFemale[idxRel] = private$popFemale[idxRel] + private$matedFemaleReleases[[1]]$nRelease[ ,3]
    # remove finished release
    private$matedFemaleReleases[[1]] = NULL


  }

}

#' Determine number of newly conceived mice in a Patch
#'
#' Based on this patch's release schedule, \code{\link{generateReleaseVector}},
#' this function handles daily conception.
#'
oneDay_gestReleases_Patch <- function(){

  # things to reuse
  tNow <- private$NetworkPointer$get_tNow()

  ##########
  # newly gestating
  ##########
  if((length(private$gestReleases) > 0) && (private$gestReleases[[1]]$tRelease <= tNow) ){
    # combine gestating pups
    # 1st column is index, 2nd column is amount
    idxRel <- private$gestReleases[[1]]$nRelease[ ,1]
    private$popJuvenile[idxRel,1] = private$popJuvenile[idxRel,1] + private$gestReleases[[1]]$nRelease[ ,2]
    # remove finished release
    private$gestReleases[[1]] = NULL
  }

}


###############################################################################
# Mating
###############################################################################

#' Deterministc Mating
#'
#' Mating is calculated as the outer product of newly emerging adult females and
#' all-current adult males, modulated by \eqn{\overline{\overline{\eta}}}, the genotype-specific
#' male mating fitness. \eqn{\overline{\overline{\eta}}} corresponds to a number of females calculated from the
#' average number of litters per year \eqn{litters/365} (rows)
#' and male (columns) genotypes, to perform any type of assortative mating. \cr
#'
oneDay_mating_deterministic_Patch <- function(){

  # Check if there are males
  if(sum(private$popMale) > 0){
    # mate females with males normalized by their mating ability
    #  add to current females

    # things to reuse
    nGeno = private$NetworkPointer$get_genotypesN()
    lits = private$NetworkPointer$get_litters()
    private$popMating[] = round(private$popUnmated * lits/365)

    # step through each female genotype to mate.
    for(i in 1:nGeno){
      private$popMatches[i, ] = private$popMating[i] *
        normalise(private$popMale * private$NetworkPointer$get_eta(i))

    private$popMating[] = 0
    }

  }
}





#' Stochastic Mating
#'
#' Mating for each newly emerging adult female genotype is sampled from a multinomial
#' distribution with probabilities equal to the adult male population vector
#' multiplied by \eqn{\overline{\overline{\eta}}}, the genotype-specific
#' male mating fitness. \eqn{\overline{\overline{\eta}}} corresponds a number of females calculated using a
#' binomial distribution, with the probability of mating for each femal equal to \eqn{litters/365} or the probability
#' of a female having a litter on any given day (rows)
#' and male (columns) genotypes, to perform any type of assortative mating. \cr
#'
oneDay_mating_stochastic_Patch <- function(){

  # things to reuse
  nGeno = private$NetworkPointer$get_genotypesN()
  lits = private$NetworkPointer$get_litters()
  # check if there are males
  if(sum(private$popMale) > 0){
    # mating probs for males
    mProb = numeric(length = nGeno)

    # loop over each female genotype, mate with available males
    for(i in 1:nGeno){
      # get mating prob for males for each female genotype
      mProb[] <- normalise(private$popMale * private$NetworkPointer$get_eta(i))

      private$popMating[i] = rbinom(n = 1,
                                    size = private$popUnmated[i],
                                    prob = lits/365)

      # check if there are females to mate, skip if not
      # check if males have a chance of mating these females, if not then skip
      #  This serves 2 purposes:
      #    1 - skip drawing a multinomial if all probs are zero
      #    2 - the multinomial fails if all probs are 0
      if( (private$popUnmated[i] > 0) && (sum(mProb) != 0) ){
        private$popMatches[i, ] = rmultinom(n = 1,
                                            size = private$popMating[i],
                                            prob = mProb)
        private$popMating[] = 0



  }
    }
  }
}





###############################################################################
# Conceive offspring
###############################################################################

#' Deterministc Gestation
#'
#'
#' Calculate the number of pups conceived by female mice following:
#' \deqn{\overline{O{T_x}} = \sum_{j=1}^{n} \Bigg{ \bigg{ {\beta*\overline{s} * \overline{ \overline{Af_{[t]}}}} * {\overline{\overline{\overline{Ih}}} \bigg} * \Lambda \Bigg}^{\top}_{ij}}}
#'
oneDay_conceive_deterministic_Patch <- function(){

  # things to reuse
  nGeno = private$NetworkPointer$get_genotypesN()
  matingFemBetaS = private$popMatches * (private$NetworkPointer$get_beta() * private$NetworkPointer$get_s())
  private$popMatches[] = 0

  #fill offspring cube with parents
  for(slice in 1:nGeno){

    private$popJuvenile[slice, 1] = sum(matingFemBetaS *
                                         private$NetworkPointer$get_drivecubeindex(NULL,NULL,slice) *
                                         private$NetworkPointer$get_tau(NULL,NULL,slice) )

  } # end loop over genotypes

}

#' Stochastic gestation
#'
#' Calculate the number of pups conceived by female mice following:
#' \deqn{\overline{O{T_x}} = \sum_{j=1}^{n} \Bigg{ \bigg{ {\beta*\overline{s} * \overline{ \overline{Af_{[t]}}}} * {\overline{\overline{\overline{Ih}}} \bigg} * \Lambda \Bigg}^{\top}_{ij}}}
#' The deterministic result for number of pups is used as the rate parameter of a Poisson-distributed
#' number of actual pups conceived.
#'
oneDay_conceive_stochastic_Patch <- function(){

  # things to reuse
  nGeno = private$NetworkPointer$get_genotypesN()
  matingFemBetaS = private$popMatches * (private$NetworkPointer$get_beta() * private$NetworkPointer$get_s())
  private$popMatches[] = 0

  #fill offspring cube with parents
  for(slice in 1:nGeno){

    private$popJuvenile[slice, 1] = sum(rpois(n = nGeno*nGeno,
                                             lambda = matingFemBetaS *
                                               private$NetworkPointer$get_drivecubeindex(NULL,NULL,slice) *
                                               private$NetworkPointer$get_tau(NULL,NULL,slice) )
    )

  } # end loop over genotypes

}


