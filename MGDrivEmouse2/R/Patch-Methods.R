###############################################################################
#        ____        __       __
#       / __ \____ _/ /______/ /_
#      / /_/ / __ `/ __/ ___/ __ \
#     / ____/ /_/ / /_/ /__/ / / /
#    /_/    \__,_/\__/\___/_/ /_/#
#
#   Patch Class Implementation
#   Original Code by: Marshall Lab
#   jared_bennett@berkeley.edu
#   December 2019
#   MODIFIED BY: ETHAN A. BROWN (JUL 9 2020)
#   ebrown23@nd.edu
###############################################################################

#' Set Initial Population
#'
#' This hidden function distributes the population at time 0 in the steady-state
#' conformation. This involves splitting adults into male and female.
#'
#' @param k Carrying capacity constant
#' @param adultRatioF Genotype specific ratio for adult females
#' @param adultRatioM Genotype specific ratio for adult males
#' @param muAI daily survival of adult stage
#' @param muJI juvenile infection mortality parameter
#' @param muN nursing mortality parameter
#' @param muG gestation mortality parameter
#' @param timeAd approximate length of adult life stage
#'

set_initialPopulation_Patch <- function(k = k, adultRatioF = adultRatioF, adultRatioM = adultRatioM,
                                        muAI = muAI, muJI = muJI, muN = muN, muG = muG, timeAd = timeAd){
  
  
  ##########
  # set male population breakdown
  ##########
  private$popMale[names(adultRatioM)] = adultRatioM * k/2
  
  ##########
  # set mated female population breakdown
  ##########
  # this isn't exactly correct, it mates all males to all females, ignoring
  # the genotype-specific male mating abilities
  private$popUnmated[names(adultRatioF)] = adultRatioF * k/2
  private$popFemale = private$popUnmated %o% normalise(private$popMale)
  
  ##########
  # set juvenile population breakdown
  ##########
  B <- muAI * k
  
  timeJu <- try(get("timeJu", inherits = TRUE), silent = TRUE)
  if(inherits(timeJu, "try-error")){
    timeJu <- c("G" = 0L, "N" = 0L, "A" = 0L)
  }
  
  tG <- as.integer(timeJu["G"])
  tN <- as.integer(timeJu["N"])
  tA <- as.integer(timeJu["A"])
  
  survA <- 1 - muJI
  survN <- 1 - muN
  survG <- 1 - muG
  
  nCols <- ncol(private$popJuvenile)
  aStart <- tG + tN + 1L
  aEnd <- tG + tN + tA
  nStart <- tG + 1L
  nEnd <- tG + tN
  
  for(g in names(private$popMale)){
    Bg <- B * (adultRatioF[g] + adultRatioM[g]) / 2
    pop_vec <- numeric(nCols)
    
    # adolescent stage
    if(tA > 0){
      pop_vec[aEnd] <- Bg / survA
      if(tA > 1){
        for(i in (aEnd-1L):aStart){
          pop_vec[i] <- pop_vec[i+1L] / survA
        }
      }
    }
    
    # nursing stage
    if(tN > 0){
      pop_vec[nEnd] <- pop_vec[aStart] / survN
      if(tN > 1){
        for(i in (nEnd-1L):nStart){
          pop_vec[i] <- pop_vec[i+1L] / survN
        }
      }
    }
    
    # gestation stage
    if(tG > 0){
      pop_vec[tG] <- pop_vec[nStart] / survG
      if(tG > 1){
        for(i in (tG-1L):1L){
          pop_vec[i] <- pop_vec[i+1L] / survG
        }
      }
    }
    
    private$popJuvenile[g, ] <- pop_vec
  }
  
  private$popJuvenilet0 = private$popJuvenile
  
  
}

#' Set Initial Population Deterministic
#'
#' Calls \code{\link{set_initialPopulation_Patch}} to initialize a steady-state
#' population distribution.
#'
#' @param k Carrying capacity constant
#' @param adultRatioF Genotype specific ratio for adult females
#' @param adultRatioM Genotype specific ratio for adult males
#' @param muAI daily survival of adult stage
#' @param muJI juvenile infection mortality parameter
#' @param muN nursing mortality parameter
#' @param muG gestation mortality parameter
#' @param timeAd approximate length of adult life stage
#'
set_population_deterministic_Patch <- function(k = k, adultRatioF = adultRatioF, adultRatioM = adultRatioM,
                                               muAI = muAI, muJI = muJI, muN = muN, muG = muG, timeAd = timeAd){
  
  self$initialPopulation(k = k, adultRatioF = adultRatioF, adultRatioM = adultRatioM,
                         muAI = muAI, muJI = muJI, muN = muN, muG = muG, timeAd = timeAd)
  
}

#' Set Initial Population Stochastic
#'
#' Calls \code{\link{set_initialPopulation_Patch}} to initialize a steady-state
#' population distribution. Populations are then rounded to integer values.
#'
#' @param k Carrying capacity constant
#' @param adultRatioF Genotype specific ratio for adult females
#' @param adultRatioM Genotype specific ratio for adult males
#' @param muAI daily survival of adult stage
#' @param muJI juvenile infection mortality parameter
#' @param muN nursing mortality parameter
#' @param muG gestation mortality parameter
#' @param timeAd approximate length of adult life stage
#'
set_population_stochastic_Patch <- function(k = k, adultRatioF = adultRatioF, adultRatioM = adultRatioM,
                                            muAI = muAI, muJI = muJI, muN = muN, muG = muG, timeAd = timeAd){
  
  # set initial population
  self$initialPopulation(k = k, adultRatioF = adultRatioF, adultRatioM = adultRatioM,
                         muAI = muAI, muJI = muJI, muN = muN, muG = muG, timeAd = timeAd)
  
  ##########
  # make everything an integer
  ##########
  private$popJuvenile[] <- round(private$popJuvenile)
  private$popUnmated[] <- round(private$popUnmated)
  private$popMale[] <- round(private$popMale)
  private$popFemale[] <- round(private$popFemale)
  
}

#' Reset Patch to Initial Conditions
#'
#' Resets a patch to its initial configuration so that a new one does not have
#' to be created and allocated in the network (for Monte Carlo simulation).
#'
#' @param verbose Chatty? Default is TRUE
#'
reset_Patch <- function(verbose = TRUE){
  
  if(verbose){cat("reset patch ",private$patchID,"\n",sep="")}
  
  # reset population
  private$popJuvenile[] = private$popJuvenilet0
  private$popMale[] = private$popMalet0
  private$popUnmated[] = private$popUnmatedt0
  private$popFemale[] = private$popFemalet0
  private$popFemale[] = private$popFemalet0
  
  
  # Reset Mosquito Releases
  private$gestReleases = private$NetworkPointer$get_patchReleases(private$patchID,"Gest")
  private$maleReleases = private$NetworkPointer$get_patchReleases(private$patchID,"M")
  private$femaleReleases = private$NetworkPointer$get_patchReleases(private$patchID,"F")
  private$matedFemaleReleases = private$NetworkPointer$get_patchReleases(private$patchID,"mF")
  
}

#' Initialize Output from Focal Patch
#'
#' Writes output to the text connections specified in the enclosing \code{\link{Network}}.
#'
oneDay_initOutput_Patch <- function(){
  
  ##########
  # headers
  ##########
  if(private$patchID == 1){
    # males
    writeLines(text = paste0(c("Time","Patch",private$NetworkPointer$get_genotypesID()), collapse = ","),
               con = private$NetworkPointer$get_conADM(), sep = "\n")
    # females
    femaleCrosses = c(t(outer(private$NetworkPointer$get_genotypesID(),private$NetworkPointer$get_genotypesID(),FUN = paste0)))
    writeLines(text = paste0(c("Time","Patch",femaleCrosses), collapse = ","),
               con = private$NetworkPointer$get_conADF(),sep = "\n")
  }
  
  ##########
  # males
  ##########
  writeLines(text = paste0(c(1,private$patchID,private$popMale),collapse = ","),
             con = private$NetworkPointer$get_conADM(), sep = "\n")
  
  ##########
  # females
  ##########
  writeLines(text = paste0(c(1,private$patchID,c(t(private$popFemale))),collapse = ","),
             con = private$NetworkPointer$get_conADF(), sep = "\n")
  
}

#' Write Output from Focal Patch
#'
#' Writes output to the text connections specified in the enclosing \code{\link{Network}}.
#'
oneDay_writeOutput_Patch <- function(){
  
  tNow = private$NetworkPointer$get_tNow()
  
  # write males
  ADMout = paste0(c(tNow,private$patchID,private$popMale),collapse = ",")
  writeLines(text = ADMout,con = private$NetworkPointer$get_conADM(),sep = "\n")
  
  # write females
  ADFout = paste0(c(tNow,private$patchID,c(t(private$popFemale))),collapse = ",")
  writeLines(text = ADFout,con = private$NetworkPointer$get_conADF(),sep = "\n")
  
}
