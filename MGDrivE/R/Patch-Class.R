###############################################################################
#        ____        __       __
#       / __ \____ _/ /______/ /_
#      / /_/ / __ `/ __/ ___/ __ \
#     / ____/ /_/ / /_/ /__/ / / /
#    /_/    \__,_/\__/\___/_/ /_/
#
#   Patch Class Definition
#   Original Code by: Marshall Lab
#   jared_bennett@berkeley.edu
#   December 2019
#   MODIFIED BY: ETHAN A. BROWN (JUL 9 2020)
#   ebrown23@nd.edu
###############################################################################

#' Patch Class Definition
#'
#' A Patch is a single well-mixed population that is the smallest unit of simulation for MGDrivE.
#'
#' @docType class
#' @format An \code{\link{R6Class}} generator object
#' @keywords R6 class
#'
#' @importFrom stats rbinom rmultinom rpois
#'
#' @section **Constructor**:
#'  * patchID: integer ID of this patch
#'  * genotypesID: character vector of genotypes
#'  * timeJu: integer vector of length 3 specifying the length of each juvenile stage
#'  * numPatches: integer, total number of patches in this simulation
#'  * adultEQ: integer, total adult population in this patch for the duration of the simulation
#'  * k: double, carrying capacity parameter, see \code{\link{parameterizeMGDrivE}}
#'  * muAI: additional adult mortality parameter
#'  * muJI: juvenile infection mortality parameter
#'  * muN: nursing mortality parameter
#'  * muG: gestation mortality parameter
#'  * adultRatioF: named double vector, distribution of adult female genotypes, see \code{\link{parameterizeMGDrivE}}
#'  * adultRatioM: named double vector, distribution of adult male genotypes, see \code{\link{parameterizeMGDrivE}}
#'  * gestReleases: gestating pup release schedule for this patch, see \code{\link{basicRepeatedReleases}}
#'  * maleReleases: male release schedule for this patch, see \code{\link{basicRepeatedReleases}}
#'  * femaleReleases: female release schedule for this patch, see \code{\link{basicRepeatedReleases}}
#'  * matedFemaleReleases: mated females release schedule for this patch, see \code{\link{basicRepeatedReleases}}
#'
#' @section **Methods**:
#'  * set_NetworkPointer: see \code{\link{set_NetworkPointer_Patch}}
#'  * get_maleMigration: see \code{\link{get_maleMigration_Patch}}
#'  * get_femaleMigration: see \code{\link{get_femaleMigration_Patch}}
#'  * initialPopulation: see \code{\link{set_initialPopulation_Patch}}
#'  * setPopulation: see \code{\link{set_population_deterministic_Patch}} or \code{\link{set_population_stochastic_Patch}}
#'  * reset: see \code{\link{reset_Patch}}
#'  * oneDay_initOutput: see \code{\link{oneDay_initOutput_Patch}}
#'  * oneDay_writeOutput: see \code{\link{oneDay_writeOutput_Patch}}
#'  * oneDay_migrationOut: see \code{\link{oneDay_migrationOut_deterministic_Patch}} or \code{\link{oneDay_migrationOut_stochastic_Patch}}
#'  * oneDay_migrationIn: see \code{\link{oneDay_migrationIn_Patch}}
#'  * oneDay_PopDynamics: see \code{\link{oneDay_PopDynamics_Patch}}
#'  * oneDay_adultD: see \code{\link{oneDay_adultDeath_deterministic_Patch}} or \code{\link{oneDay_adultDeath_stochastic_Patch}}
#'  * oneDay_adoDM: see \code{\link{oneDay_adoDM_deterministic_Patch}} or \code{\link{oneDay_adoDM_stochastic_Patch}}
#'  * oneDay_nursingDM: see \code{\link{oneDay_nursingDM_deterministic_Patch}} or \code{\link{oneDay_nursingDM_stochastic_Patch}}
#'  * oneDay_gestDM: see \code{\link{oneDay_gestDM_deterministic_Patch}} or \code{\link{oneDay_gestDM_stochastic_Patch}}
#'  * oneDay_maturation: see \code{\link{oneDay_maturation_deterministic_Patch}} or \code{\link{oneDay_maturation_stochastic_Patch}}
#'  * oneDay_releases: see \code{\link{oneDay_releases_Patch}}
#'  * oneDay_releasePups: see \code{\link{oneDay_gestReleases_Patch}}
#'  * oneDay_mating: see \code{\link{oneDay_mating_deterministic_Patch}} or \code{\link{oneDay_mating_stochastic_Patch}}
#'  * oneDay_layPups: see \code{\link{oneDay_oviposit_deterministic_Patch}} or \code{\link{oneDay_oviposit_stochastic_Patch}}
#'
#' @section **Fields**:
#'  * patchID: integer ID of this patch
#'  * popJuvenile: matrix, nGenotype x sum(timeJuvenile), holding all gestating pups, nursing pups, and adolescent pups
#'  * popMale: vector, nGenotype x 1, holds adult males
#'  * popFemale: matrix, nGenotype x nGenotype, holds mated adult females
#'  * popHolder: vector, nGenotype x 1, temporary population storage
#'  * popPupSex: vector, nGenotype x 1, used in stochastic maturation as another temporary population
#'  * popUnmated: vector, nGenotype x 1, holds unmated females
#'  * mMig: matrix, nGenotype x nPatches, holds outbound males for migration, see \code{\link{oneDay_migrationOut_deterministic_Patch}} or \code{\link{oneDay_migrationOut_stochastic_Patch}}
#'  * fMig: array, nGenotype x nGenotype x nPatches, holds outbound females for migration, see \code{\link{oneDay_migrationOut_deterministic_Patch}} or \code{\link{oneDay_migrationOut_stochastic_Patch}}
#'  * popJuvenilet0: matrix, nGenotype x sum(timeJuvenile), holding all gestating pups, nursing pups, and maturation for reset, see \code{\link{reset_Patch}}
#'  * popMalet0: vector, nGenotype x 1, holds adult males for reset see \code{\link{reset_Patch}}
#'  * popFemalet0: matrix, nGenotype x nGenotype, holds mated adult females for reset see \code{\link{reset_Patch}}
#'  * gestReleases: list of gestating pup releases for this patch. See \code{\link{oneDay_gestReleases_Patch}}
#'  * maleReleases: list of adult male releases for this patch. See \code{\link{oneDay_releases_Patch}}
#'  * femaleReleases: list of adult female releases for this patch. See \code{\link{oneDay_releases_Patch}}
#'  * matedFemaleReleases: list of mated adult female releases for this patch. See \code{\link{oneDay_releases_Patch}}
#'  * NetworkPointer: a reference to enclosing \code{\link{Network}}
#'
Patch <- R6::R6Class(classname = "Patch",
            portable = TRUE,
            cloneable = FALSE,
            lock_class = FALSE,
            lock_objects = FALSE,
            class = FALSE,

            # public members
            public = list(

                #################################################
                # Constructor
                #################################################

                initialize = function(patchID, genotypesID, timeJu, timeAd, numPatches,
                                      k, muAd, muAI, muJI, muN, muG,
                                      adultRatioF, adultRatioM,
                                      gestReleases = NULL,
                                      maleReleases = NULL,
                                      femaleReleases = NULL,
                                      matedFemaleReleases = NULL
                                      ){

                  # ID of this patch
                  private$patchID = patchID

                  # initialize objects for simulation. This way, they have dimensions and names
                  nGeno = length(genotypesID)
                  private$popJuvenile = matrix(data = 0,  nrow = nGeno, ncol = sum(timeJu),
                                              dimnames = list(genotypesID, NULL))
                  #private$popAdult = matrix(data = 0,  nrow = nGeno, ncol = 2,
                                              #dimnames = list(genotypesID, c("M","F")))
                  private$popMale = setNames(object = numeric(length = nGeno), nm = genotypesID)
                  private$popFemale = matrix(data = 0, nrow = nGeno, ncol = nGeno,
                                             dimnames = list(genotypesID, genotypesID))
                  private$popHolder = setNames(object = numeric(length = nGeno), nm = genotypesID)
                  private$popPupSex = setNames(object = numeric(length = nGeno), nm = genotypesID)
                  private$popUnmated = setNames(object = numeric(length = nGeno), nm = genotypesID)
                  private$popMating = setNames(object = numeric(length = nGeno), nm = genotypesID)
                  private$popMatches = matrix(data = 0, nrow = nGeno, ncol = nGeno,
                                             dimnames = list(genotypesID, genotypesID))

                  private$mMig = matrix(data=0, nrow=nGeno, ncol=numPatches)
                  private$fMig = array(data = 0, dim=c(nGeno,nGeno,numPatches))

                  private$muAI = muAI
                  private$muJI = muJI
                  private$muN = muN
                  private$muG = muG

                  # set initial population
                  self$setPopulation(k = k,
                                     adultRatioF = adultRatioF,
                                     adultRatioM = adultRatioM,
                                     timeAd = timeAd,
                                     muAd = muAd)

                  # store reset populations
                  private$popJuvenilet0 = private$popJuvenile
                  private$popMalet0 = private$popMale
                  private$popUnmatedt0 = private$popUnmated
                  private$popFemalet0 = private$popFemale

                  # Mosquito Releases
                  private$gestReleases = gestReleases
                  private$maleReleases = maleReleases
                  private$femaleReleases = femaleReleases
                  private$matedFemaleReleases = matedFemaleReleases

                } # end constructor
              ),

            # private members
            private = list(

              patchID = NULL,

              # temporary populations
              popJuvenile = NULL,
              popMale = NULL,
              popFemale = NULL,
              popHolder = NULL,
              popPupSex = NULL, # only used in stochastic maturation function
              popUnmated = NULL,

              muAI = NULL,
              muJI = NULL,
              muN = NULL,
              muG = NULL,

              # migration
              mMig = NULL,
              fMig = NULL,

              # reset populations
              popJuvenilet0 = NULL,
              popMalet0 = NULL,
              popUnmatedt0 = NULL,
              popFemalet0 = NULL,
              popMatingt0 = NULL,
              popMatchest0 = NULL,

              # releases
              gestReleases = NULL,
              maleReleases = NULL,
              femaleReleases = NULL,
              matedFemaleReleases = NULL,

              # pointers
              NetworkPointer = NULL

            ) # end private list
)


###############################################################################
# Getters & Setters
###############################################################################

#' Set Network Pointer
#'
#' Set a reference to the enclosing \code{\link{Network}} object
#'
#' @param NetworkPointer A \code{\link{Network}} object
#'
set_NetworkPointer_Patch <- function(NetworkPointer){private$NetworkPointer = NetworkPointer}

Patch$set(which = "public",name = "set_NetworkPointer",
          value = set_NetworkPointer_Patch,overwrite = TRUE
)

#' Get maleMigration
#'
#' Return outbound males (nGenotypes X nPatch integer matrix)
#'
get_maleMigration_Patch <- function(){return(private$mMig)}

Patch$set(which = "public",name = "get_maleMigration",
          value = get_maleMigration_Patch,overwrite = TRUE
)

#' Get femaleMigration
#'
#' Return outbound females (nGenotypes X nGenotypes X nPatch array)
#'
get_femaleMigration_Patch <- function(){return(private$fMig)}

Patch$set(which = "public",name = "get_femaleMigration",
          value = get_femaleMigration_Patch,overwrite = TRUE
)
