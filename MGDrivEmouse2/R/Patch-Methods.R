###############################################################################
#        ____        __       __
#       / __ \____ _/ /______/ /_
#      / /_/ / __ `/ __/ ___/ __ \
#     / ____/ /_/ / /_/ /__/ / / /
#    /_/    \__,_/\__/\___/_/ /_/
#
#   Patch Class Implementation (Rewritten)
#   Original Code by: Marshall Lab
#   jared_bennett@berkeley.edu
#   December 2019
#
#   Rewritten by: ChatGPT (MouseGD equilibrium-friendly init)
#   Notes:
#     * Keeps external interface compatible with prior code: initialPopulation,
#       set_population_*_Patch, reset_Patch, init/write output.
#     * Adds equilibrium seeding of juvenile pipeline (gestation/nursing/adolescence)
#       with optional density-independent juvenile mortalities.
#     * Balances adult emergence (t+1) to expected adult deaths (t) given the
#       current adult state, K, theta, omega, muAd, and (optionally) toxicant.
#     * Uses NetworkPointer getters for timeJu & parameters; avoids try(get(...)).
#     * Fixes duplicate popFemale reset; rounds only in stochastic init.
###############################################################################

# ---------- Utilities ---------------------------------------------------------

# defensive normalizer (should already exist in codebase)
normalise <- function(x){
  s <- sum(x)
  if(s > 0){ x / s } else { x }
}

# clip helper
.clip_nonneg <- function(x){ ifelse(is.finite(x) & x >= 0, x, 0) }

# ---------- Initial Population -----------------------------------------------

#' Set Initial Population (adults + female mating matrix; juveniles seeded later)
#'
#' @param k Carrying capacity constant (expected adults at equilibrium ~ K)
#' @param adultRatioF Named vector: genotype-specific ratio for adult females
#' @param adultRatioM Named vector: genotype-specific ratio for adult males
#' @param use_eta_at_t0 logical; if TRUE, use genotype-specific male mating
#'        fitness (eta) when building popFemale at t0; otherwise proportional to
#'        male abundance only (backwards compatible). Default FALSE.
#'
set_initialPopulation_Patch <- function(k = k, adultRatioF = adultRatioF, adultRatioM = adultRatioM, 
                                        use_eta_at_t0 = FALSE, ...){
  
  # ---- normalize ratio vectors to be safe ----
  if(sum(adultRatioM) <= 0) stop("adultRatioM must sum to > 0")
  if(sum(adultRatioF) <= 0) stop("adultRatioF must sum to > 0")
  adultRatioM <- adultRatioM / sum(adultRatioM)
  adultRatioF <- adultRatioF / sum(adultRatioF)
  
  # ---- set adult population breakdown (K/2 per sex) ----
  private$popMale[names(adultRatioM)]    <- adultRatioM * (k/2)
  private$popUnmated[names(adultRatioF)] <- adultRatioF * (k/2)
  
  # ---- build female mating matrix popFemale (rows=female genos, cols=male genos) ----
  if(sum(private$popMale) > 0){
    if(isTRUE(use_eta_at_t0)){
      # include eta per female genotype
      nG <- private$NetworkPointer$get_genotypesN()
      pf <- matrix(0, nrow = nG, ncol = nG)
      for(i in 1:nG){
        maleProb_i <- normalise(private$popMale * private$NetworkPointer$get_eta(i))
        pf[i, ] <- private$popUnmated[i] * maleProb_i
      }
      private$popFemale[] <- pf
    } else {
      # proportional to male abundance (legacy behavior)
      private$popFemale[] <- private$popUnmated %o% normalise(private$popMale)
    }
  } else {
    private$popFemale[] <- 0
  }
}

# ---------- Juvenile Equilibrium Seeding -------------------------------------

#' Seed Juveniles to Equilibrium (fills gestation, nursing, adolescence)
#'
#' Uses current adults to compute expected daily offspring by offspring genotype
#' (deterministic expectation), optionally scaled so that expected adult
#' emergence tomorrow matches expected adult deaths today. Fills the juvenile
#' pipeline forward using geometric survival within each stage.
#'
seed_juveniles_equilibrium_Patch <- function(){
  nG <- private$NetworkPointer$get_genotypesN()
  
  # ----- Stage lengths -----
  tG <- private$NetworkPointer$get_timeJu(stage = 'G')
  tN <- private$NetworkPointer$get_timeJu(stage = 'N')
  tTot <- private$NetworkPointer$get_timeJu()         # total columns
  tA <- tTot - (tG + tN)
  if(tA < 0) stop("Inconsistent juvenile times: A < 0")
  
  # ----- Juvenile (density-independent) survival per day -----
  # These getters must exist in your backend; otherwise set to 1 for no mortality.
  sG <- 1 - private$NetworkPointer$get_muG()   # gestation DI mortality
  sN <- 1 - private$NetworkPointer$get_muN()   # nursing DI mortality
  sA <- 1 - private$NetworkPointer$get_muJI()  # adolescence DI mortality
  
  # ----- Adult survival today (as in oneDay_adultDeath_*_Patch) -----
  muAd  <- private$NetworkPointer$get_muAd()
  omega <- private$NetworkPointer$get_omega()     # vector length nG
  K     <- private$NetworkPointer$get_k(ix = private$patchID)
  theta <- private$NetworkPointer$get_theta()
  
  A_m <- sum(private$popMale)
  A_f <- sum(private$popUnmated)
  density <- ((K/2) / (K + A_m + A_f))^(1/theta)
  
  # base adult daily survival per genotype (both sexes share same formula in code)
  pSurv_g <- omega * ((1 - muAd) * density)
  
  # include toxicant window if active at t0
  tNow <- private$NetworkPointer$get_tNow()
  interval <- sort(private$NetworkPointer$get_toxInt())
  dayToxSurv <- 1 - private$NetworkPointer$get_toxMort()
  if(length(interval) == 2 && is.finite(dayToxSurv) && dayToxSurv < 1 &&
     tNow >= interval[1] && tNow <= interval[2]){
    pSurv_g <- pSurv_g * dayToxSurv
  }
  
  # expected adult deaths today
  expDeaths <- sum(private$popMale    * (1 - pSurv_g)) +
    sum(private$popUnmated * (1 - pSurv_g))
  expDeaths <- max(expDeaths, 0)
  
  # ----- Expected births by offspring genotype (deterministic expectation) -----
  lits <- private$NetworkPointer$get_litters()   # litters per year
  beta <- private$NetworkPointer$get_beta()
  sFec <- private$NetworkPointer$get_s()
  
  # expected mated females by genotype (same simplification as t0 female matrix)
  matedFem_i <- private$popUnmated * (lits/365)
  maleProb <- if(sum(private$popMale) > 0) normalise(private$popMale) else rep(0, nG)
  popMatches_exp <- matrix(0, nrow = nG, ncol = nG)
  for(i in 1:nG){
    popMatches_exp[i, ] <- matedFem_i[i] * maleProb
  }
  
  matingFemBetaS <- popMatches_exp * (beta * sFec)
  
  births_vec <- numeric(nG)
  for (g in 1:nG){
    births_vec[g] <- sum(
      matingFemBetaS *
        private$NetworkPointer$get_drivecubeindex(NULL, NULL, g) *
        private$NetworkPointer$get_tau(NULL, NULL, g)
    )
  }
  
  # ----- Scale births so E[emergence tomorrow] = E[deaths today] -----
  phi <- private$NetworkPointer$get_phi()
  xiM <- private$NetworkPointer$get_xiM()
  xiF <- private$NetworkPointer$get_xiF()
  
  emerg_factor <- (sG^tG) * (sN^tN) * (sA^tA) * (1 - muAd) * ((1 - phi) * xiM + phi * xiF)
  denom <- sum(births_vec * emerg_factor)
  alpha <- if (denom > 0) expDeaths / denom else 0
  alpha <- .clip_nonneg(alpha)
  
  births_vec <- alpha * births_vec
  
  # ----- Fill juvenile pipeline forward with geometric survival -----
  private$popJuvenile[,] <- 0
  # gestation columns 1..tG
  if (tG > 0) {
    for (j in 1:tG) {
      private$popJuvenile[, j] <- births_vec * (sG^(j-1))
    }
  }
  # nursing columns tG+1..tG+tN
  if (tN > 0) {
    for (j in 1:tN) {
      private$popJuvenile[, tG + j] <- births_vec * (sG^tG) * (sN^(j-1))
    }
  }
  # adolescence columns tG+tN+1..tTot
  if (tA > 0) {
    for (j in 1:tA) {
      private$popJuvenile[, tG + tN + j] <- births_vec * (sG^tG) * (sN^tN) * (sA^(j-1))
    }
  }
  
  # save t0 copy
  private$popJuvenilet0[] <- private$popJuvenile
}

# ---------- Public Init Wrappers ---------------------------------------------

#' Set Initial Population (Deterministic)
set_population_deterministic_Patch <- function(k = k,
                                               adultRatioF = adultRatioF,
                                               adultRatioM = adultRatioM,
                                               use_eta_at_t0 = FALSE, ...){
  self$initialPopulation(k = k,
                         adultRatioF = adultRatioF,
                         adultRatioM = adultRatioM,
                         use_eta_at_t0 = use_eta_at_t0, ...)
  
  # equilibrium seed juveniles
  self$seedJuvenilesEquilibrium()
}

#' Set Initial Population (Stochastic)
set_population_stochastic_Patch <- function(k = k,
                                            adultRatioF = adultRatioF,
                                            adultRatioM = adultRatioM,
                                            use_eta_at_t0 = FALSE, ...){
  # set initial adults & female matrix
  self$initialPopulation(k = k,
                         adultRatioF = adultRatioF,
                         adultRatioM = adultRatioM,
                         use_eta_at_t0 = use_eta_at_t0, ...)
  
  # seed juveniles using expectation, then round to integers
  self$seedJuvenilesEquilibrium()
  
  private$popJuvenile[] <- round(private$popJuvenile)
  private$popUnmated[]  <- round(private$popUnmated)
  private$popMale[]     <- round(private$popMale)
  private$popFemale[]   <- round(private$popFemale)
}

# ---------- Reset -------------------------------------------------------------

#' Reset Patch to Initial Conditions
#'
#' Restores t0 state and reloads scheduled releases.
reset_Patch <- function(verbose = TRUE){
  if(verbose){cat("reset patch ", private$patchID, "\n", sep = "")}
  
  # reset population
  private$popJuvenile[] <- private$popJuvenilet0
  private$popMale[]     <- private$popMalet0
  private$popUnmated[]  <- private$popUnmatedt0
  private$popFemale[]   <- private$popFemalet0
  
  # Reset Releases
  private$gestReleases        <- private$NetworkPointer$get_patchReleases(private$patchID, "Gest")
  private$maleReleases        <- private$NetworkPointer$get_patchReleases(private$patchID, "M")
  private$femaleReleases      <- private$NetworkPointer$get_patchReleases(private$patchID, "F")
  private$matedFemaleReleases <- private$NetworkPointer$get_patchReleases(private$patchID, "mF")
}

# ---------- Output ------------------------------------------------------------

#' Initialize Output from Focal Patch
oneDay_initOutput_Patch <- function(){
  # headers (once, from patch 1)
  if(private$patchID == 1){
    # males header
    writeLines(text = paste0(c("Time","Patch", private$NetworkPointer$get_genotypesID()), collapse = ","),
               con = private$NetworkPointer$get_conADM(), sep = "\n")
    
    # females header (female x male crosses)
    femaleCrosses <- c(t(outer(private$NetworkPointer$get_genotypesID(),
                               private$NetworkPointer$get_genotypesID(),
                               FUN = paste0)))
    writeLines(text = paste0(c("Time","Patch", femaleCrosses), collapse = ","),
               con = private$NetworkPointer$get_conADF(), sep = "\n")
  }
  
  # write t=1 snapshot
  writeLines(text = paste0(c(1, private$patchID, private$popMale), collapse = ","),
             con = private$NetworkPointer$get_conADM(), sep = "\n")
  
  writeLines(text = paste0(c(1, private$patchID, c(t(private$popFemale))), collapse = ","),
             con = private$NetworkPointer$get_conADF(), sep = "\n")
}

#' Write Output from Focal Patch (each day)
oneDay_writeOutput_Patch <- function(){
  tNow <- private$NetworkPointer$get_tNow()
  
  # males
  ADMout <- paste0(c(tNow, private$patchID, private$popMale), collapse = ",")
  writeLines(text = ADMout, con = private$NetworkPointer$get_conADM(), sep = "\n")
  
  # females
  ADFout <- paste0(c(tNow, private$patchID, c(t(private$popFemale))), collapse = ",")
  writeLines(text = ADFout, con = private$NetworkPointer$get_conADF(), sep = "\n")
}

# ---------- Public Facade Methods (bind to private impl) ----------------------

# These are convenience aliases you can add in the R6 class definition to
# expose the rewritten helpers with your preferred names, e.g.:
#   seedJuvenilesEquilibrium = function() private$seed_juveniles_equilibrium_Patch()
#   initialPopulation        = function(...) private$set_initialPopulation_Patch(...)
# Make sure to wire them in your class constructor.
