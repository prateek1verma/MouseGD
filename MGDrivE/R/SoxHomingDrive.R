
SoxHomingDrive <- function(cM = 1.0, cF = 1.0, chM = 0, crM = 0, chF = 0,
                            crF = 0, eta = NULL, phi = NULL,
                            omega = NULL, xiF = NULL, xiM = NULL, s = NULL){
  
  ## safety cWecks in case someone is dumb
  if(any(c(cM, cF, chM, crM, chF, crF)>1) || any(c(cM, cF,chM, crM, chF, crF)<0)){
    stop("Parameters are rates, they must be between 0 and 1.")
  }
  
  
  ## define matrices
  ## Matrix Dimensions Key: [femaleGenotype,maleGenotype,offspringGenotype]
  gtype <- c('mWW', 'mWH', 'mWR', 'mWB', 'mHH', 'mHR', 'mHB', 'mRR', 'mRB', 'mBB',
            'fWW', 'fWH', 'fWR', 'fWB', 'fHH', 'fHR', 'fHB', 'fRR', 'fRB', 'fBB')
 
  
  size <- length(gtype)
  
  tMatrix <- array(data=0, dim=c(size, size, size), dimnames=list(gtype, gtype, gtype)) #transition matrix
  
  ## fill tMatrix with probabilities
  #('WW', 'WH', 'WR', 'WB', 'HH', 'HR', 'HB', 'RR', 'RB', 'BB')

  tMatrix['fWW','mWW', c('fWW','mWW')] <- 1/2
  
  tMatrix['fWR','mWW', c('fWW', 'fWR','mWW', 'mWR')] <- c( 1, 1)/4
  tMatrix['fWR','mWR', c('fWW', 'fWR', 'fRR','mWW', 'mWR', 'mRR')] <- c( 1/2, 1, 1/2)/4
  
  tMatrix['fWB','mWW', c('fWW', 'fWB','mWW', 'mWB')] <- c( 1, 1)/4
  tMatrix['fWB','mWR', c('fWW', 'fWR', 'fWB', 'fRB','mWW', 'mWR', 'mWB', 'mRB')] <- c( 1, 1, 1, 1)/8
  tMatrix['fWB','mWB', c('fWW', 'fWB', 'fBB','mWW', 'mWB', 'mBB')] <- c( 1/2, 1, 1/2)/4
  
  tMatrix['fHH','mHH', c('fHH','mHH')] <- 1/2
  
  tMatrix['fHR','mHH', c('fHH', 'fHR','mHH', 'mHR')] <- c( 1, 1)/4
  tMatrix['fHR','mHR', c('fHH', 'fHR', 'fRR','mHH', 'mHR', 'mRR')] <- c( 1/2, 1, 1/2)/4
  
  tMatrix['fHB','mHH', c('fHH', 'fHB','mHH', 'mHB')] <- c( 1, 1)/4
  tMatrix['fHB','mHR', c('fHH', 'fHR', 'fHB', 'fRB','mHH', 'mHR', 'mHB', 'mRB')] <- c( 1, 1, 1, 1)/8
  tMatrix['fHB','mHB', c('fHH', 'fHB', 'fBB','mHH', 'mHB', 'mBB')] <- c( 1/2, 1, 1/2)/4
  
  tMatrix['fRR','mWW', c('fWR','mWR')] <- 1/2
  tMatrix['fRR','mWR', c('fWR', 'fRR','mWR', 'mRR')] <- c( 1, 1)/4
  tMatrix['fRR','mWB', c('fWR', 'fRB','mWR', 'mRB')] <- c( 1, 1)/4
  tMatrix['fRR','mHH', c('fHR','mHR')] <- 1/2
  tMatrix['fRR','mHR', c('fHR', 'fRR','mHR', 'mRR')] <- c( 1, 1)/4
  tMatrix['fRR','mHB', c('fHR', 'fRB','mHR', 'mRB')] <- c( 1, 1)/4
  tMatrix['fRR','mRR', c('fRR','mRR')] <- 1/2
  
  tMatrix['fRB','mWW', c('fWR', 'fWB','mWR', 'mWB')] <- c( 1, 1)/4
  tMatrix['fRB','mWR', c('fWR', 'fWB', 'fRR', 'fRB','mWR', 'mWB', 'mRR', 'mRB')] <- c( 1, 1, 1, 1)/8
  tMatrix['fRB','mWB', c('fWR', 'fWB', 'fRB', 'fBB','mWR', 'mWB', 'mRB', 'mBB')] <- c( 1, 1, 1, 1)/8
  tMatrix['fRB','mHH', c('fHR', 'fHB','mHR', 'mHB')] <- c( 1, 1)/4
  tMatrix['fRB','mHR', c('fHR', 'fHB', 'fRR', 'fRB','mHR', 'mHB', 'mRR', 'mRB')] <- c( 1, 1, 1, 1)/8
  tMatrix['fRB','mHB', c('mHR', 'mHB', 'mRB', 'mBB','fHR', 'fHB', 'fRB', 'fBB')] <- c( 1, 1, 1, 1)/8
  tMatrix['fRB','mRR', c('fRR', 'fRB','mRR', 'mRB')] <- c( 1, 1)/4
  tMatrix['fRB','mRB', c('fRR', 'fRB', 'fBB','mRR', 'mRB', 'mBB')] <- c( 1/2, 1, 1/2)/4
  
  tMatrix['fBB','mWW', c('fWB','mWB')] <- 1/2
  tMatrix['fBB','mWR', c('fWB', 'fRB','mWB', 'mRB')] <- c( 1, 1)/4
  tMatrix['fBB','mWB', c('fWB', 'fBB','mWB', 'mBB')] <- c( 1, 1)/4
  tMatrix['fBB','mHH', c('fHB','mHB')] <- 1/2
  tMatrix['fBB','mHR', c('fHB', 'fRB','mHB', 'mRB')] <- c( 1, 1)/4
  tMatrix['fBB','mHB', c('fHB', 'fBB','mHB', 'mBB')] <- c( 1, 1)/4
  tMatrix['fBB','mRR', c('fRB','mRB')] <- 1/2
  tMatrix['fBB','mRB', c('fRB', 'fBB','mRB', 'mBB')] <- c( 1, 1)/4
  tMatrix['fBB','mBB', c('fBB','mBB')] <- 1/2
  
  #####
  
  tMatrix['fWW','mWR', c('fWW', 'fWR','mWW', 'mWR')] <- c( 1, 1)/4
  
  tMatrix['fWW','mWB', c('fWW', 'fWB','mWW', 'mWB')] <- c( 1, 1)/4
  tMatrix['fWR','mWB', c('fWW', 'fWR', 'fWB', 'fRB','mWW', 'mWR', 'mWB', 'mRB')] <- c( 1, 1, 1, 1)/8
  
  tMatrix['fHH','mHR', c('fHH', 'fHR','mHH', 'mHR')] <- c( 1, 1)/4
  
  tMatrix['fHH','mHB', c('fHH', 'fHB','mHH', 'mHB')] <- c( 1, 1)/4
  tMatrix['fHR','mHB', c('fHH', 'fHR', 'fHB', 'fRB','mHH', 'mHR', 'mHB', 'mRB')] <- c( 1, 1, 1, 1)/8
  
  tMatrix['fWW','mRR', c('fWR','mWR')] <- 1/2
  tMatrix['fWR','mRR', c('fWR', 'fRR','mWR', 'mRR')] <- c( 1, 1)/4
  tMatrix['fWB','mRR', c('fWR', 'fRB','mWR', 'mRB')] <- c( 1, 1)/4
  tMatrix['fHH','mRR', c('fHR','mHR')] <- 1/2
  tMatrix['fHR','mRR', c('fHR', 'fRR','mHR', 'mRR')] <- c( 1, 1)/4
  tMatrix['fHB','mRR', c('fHR', 'fRB','mHR', 'mRB')] <- c( 1, 1)/4
  
  tMatrix['fWW','mRB', c('fWR', 'fWB','mWR', 'mWB')] <- c( 1, 1)/4
  tMatrix['fWR','mRB', c('fWR', 'fWB', 'fRR', 'fRB','mWR', 'mWB', 'mRR', 'mRB')] <- c( 1, 1, 1, 1)/8
  tMatrix['fWB','mRB', c('fWR', 'fWB', 'fRB', 'fBB','mWR', 'mWB', 'mRB', 'mBB')] <- c( 1, 1, 1, 1)/8
  tMatrix['fHH','mRB', c('fHR', 'fHB','mHR', 'mHB')] <- c( 1, 1)/4
  tMatrix['fHR','mRB', c('fHR', 'fHB', 'fRR', 'fRB','mHR', 'mHB', 'mRR', 'mRB')] <- c( 1, 1, 1, 1)/8
  tMatrix['fHB','mRB', c('mHR', 'mHB', 'mRB', 'mBB','fHR', 'fHB', 'fRB', 'fBB')] <- c( 1, 1, 1, 1)/8
  tMatrix['fRR','mRB', c('fRR', 'fRB','mRR', 'mRB')] <- c( 1, 1)/4
  
  tMatrix['fWW','mBB', c('fWB','mWB')] <- 1/2
  tMatrix['fWR','mBB', c('fWB', 'fRB','mWB', 'mRB')] <- c( 1, 1)/4
  tMatrix['fWB','mBB', c('fWB', 'fBB','mWB', 'mBB')] <- c( 1, 1)/4
  tMatrix['fHH','mBB', c('fHB','mHB')] <- 1/2
  tMatrix['fHR','mBB', c('fHB', 'fRB','mHB', 'mRB')] <- c( 1, 1)/4
  tMatrix['fHB','mBB', c('fHB', 'fBB','mHB', 'mBB')] <- c( 1, 1)/4
  tMatrix['fRR','mBB', c('fRB','mRB')] <- 1/2
  tMatrix['fRB','mBB', c('fRB', 'fBB','mRB', 'mBB')] <- c( 1, 1)/4
  
  ## set the other half of the matrix that is symmetric
  # Boolean matrix for subsetting, used several times
  boolMat <- upper.tri(x = tMatrix[ , ,1], diag = FALSE)
  # loop over depth, set upper triangle
  for(z in 1:size){tMatrix[ , ,z][boolMat] <- t(tMatrix[ , ,z])[boolMat]}
  
  #####
  
  tMatrix['fWW','mHH', c('fWH','mWH')] <- 1/2
  tMatrix['fWR','mHH', c('fWH', 'fHR','mWH', 'mHR')] <- c( 1, 1)/4
  tMatrix['fWB','mHH', c('fWH', 'fHB','mWH', 'mHB')] <- c( 1, 1)/4
  
  tMatrix['fWW','mHR', c('fWH', 'fWR','mWH', 'mWR')] <- c( 1, 1)/4
  tMatrix['fWR','mHR', c('fWH', 'fWR', 'fHR', 'fRR','mWH', 'mWR', 'mHR', 'mRR')] <- c( 1, 1, 1, 1)/8
  tMatrix['fWB','mHR', c('fWH', 'fWR', 'fHB', 'fRB','mWH', 'mWR', 'mHB', 'mRB')] <- c( 1, 1, 1, 1)/8
  
  tMatrix['fWW','mHB', c('fWH', 'fWB','mWH', 'mWB')] <- c( 1, 1)/4
  tMatrix['fWR','mHB', c('fWH', 'fWB', 'fHR', 'fRB','mWH', 'mWB', 'mHR', 'mRB')] <- c( 1, 1, 1, 1)/8
  tMatrix['fWB','mHB', c('fWH', 'fWB', 'fHB', 'fBB','mWH', 'mWB', 'mHB', 'mBB')] <- c( 1, 1, 1, 1)/8
  
  
  tMatrix['fHH','mWW', c('fWH','mWH')] <- 1/2
  tMatrix['fHH','mWR', c('fWH', 'fHR','mWH', 'mHR')] <- c( 1, 1)/4
  tMatrix['fHH','mWB', c('fWH', 'fHB','mWH', 'mHB')] <- c( 1, 1)/4
  
  tMatrix['fHR','mWW', c('fWH', 'fWR','mWH', 'mWR')] <- c( 1, 1)/4
  tMatrix['fHR','mWR', c('fWH', 'fWR', 'fHR', 'fRR','mWH', 'mWR', 'mHR', 'mRR')] <- c( 1, 1, 1, 1)/8
  tMatrix['fHR','mWB', c('fWH', 'fWR', 'fHB', 'fRB','mWH', 'mWR', 'mHB', 'mRB')] <- c( 1, 1, 1, 1)/8
  
  tMatrix['fHB','mWW', c('fWH', 'fWB','mWH', 'mWB')] <- c( 1, 1)/4
  tMatrix['fHB','mWR', c('fWH', 'fWB', 'fHR', 'fRB','mWH', 'mWB', 'mHR', 'mRB')] <- c( 1, 1, 1, 1)/8
  tMatrix['fHB','mWB', c('fWH', 'fWB', 'fHB', 'fBB','mWH', 'mWB', 'mHB', 'mBB')] <- c( 1, 1, 1, 1)/8
  
  
  # hetrozygous homing
  tMatrix['fWW','mWH', c('fWW', 'fWH', 'fWR', 'fWB','mWW', 'mWH', 'mWR', 'mWB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/4
  tMatrix['fWR','mWH', c('fWW', 'fWH', 'fWR', 'fWB',
                       'fHR', 'fRR', 'fRB','mWW', 'mWH', 'mWR', 'mWB',
                       'mHR', 'mRR', 'mRB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM + 1-cM, cM*(1-chM)*(1-crM),
                                               1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/8
  tMatrix['fWB','mWH', c('fWW', 'fWH', 'fWR', 'fWB',
                       'fHB', 'fRB', 'fBB','mWW', 'mWH', 'mWR', 'mWB',
                       'mHB', 'mRB', 'mBB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM) + 1-cM,
                                               1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/8
  tMatrix['fHH','mWH', c('fWH', 'fHH',
                       'fHR', 'fHB','mWH', 'mHH',
                       'mHR', 'mHB')] <- c((1-cM), 1+cM*chM,
                                         cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/4
  tMatrix['fHR','mWH', c('fWH', 'fHH', 'fHR', 'fHB',
                       'fWR', 'fRR', 'fRB','mWH', 'mHH', 'mHR', 'mHB',
                       'mWR', 'mRR', 'mRB')] <- c((1-cM), 1+cM*chM, cM*(1-chM)*crM + 1+cM*chM,
                                               cM*(1-chM)*(1-crM), (1-cM),
                                               cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/8
  tMatrix['fHB','mWH', c('fWH', 'fHH', 'fHR', 'fHB',
                       'fWB', 'fRB', 'fBB','mWH', 'mHH', 'mHR', 'mHB',
                       'mWB', 'mRB', 'mBB')] <- c((1-cM), 1+cM*chM, cM*(1-chM)*crM,
                                                  cM*(1-chM)*(1-crM) + 1+cM*chM, (1-cM),
                                                  cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/8
  tMatrix['fRR','mWH', c('fWR', 'fHR', 'fRR', 'fRB','mWR', 'mHR', 'mRR', 'mRB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/4
  tMatrix['fRB','mWH', c('fWR', 'fHR', 'fRR', 'fRB',
                       'fWB', 'fHB', 'fBB','mWR', 'mHR', 'mRR', 'mRB',
                       'mWB', 'mHB', 'mBB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM) + cM*(1-chM)*crM,
                                               1-cM, 1+cM*chM, cM*(1-chM)*(1-crM))/8
  tMatrix['fBB','mWH', c('fWB', 'fHB', 'fRB', 'fBB','mWB', 'mHB', 'mRB', 'mBB')] <- c(1-cM, 1+cM*chM, cM*(1-chM)*crM, cM*(1-chM)*(1-crM))/4
  tMatrix['fWH','mWH', ] <- c((1-cF)*(1-cM), (1+cF*chF)*(1-cM) + (1-cF)*(1+cM*chM),
                            cF*(1-chF)*crF*(1-cM) + (1-cF)*cM*(1-chM)*crM,
                            cF*(1-chF)*(1-crF)*(1-cM) + (1-cF)*cM*(1-chM)*(1-crM),
                            (1+cF*chF)*(1+cM*chM),
                            cF*(1-chF)*crF*(1+cM*chM) + (1+cF*chF)*cM*(1-chM)*crM,
                            cF*(1-chF)*(1-crF)*(1+cM*chM) + (1+cF*chF)*cM*(1-chM)*(1-crM),
                            cF*(1-chF)*crF*cM*(1-chM)*crM,
                            cF*(1-chF)*(1-crF)*cM*(1-chM)*crM + cF*(1-chF)*crF*cM*(1-chM)*(1-crM),
                            cF*(1-chF)*(1-crF)*cM*(1-chM)*(1-crM))/8
  tMatrix['fWH','mWW', c('fWW', 'fWH', 'fWR', 'fWB','mWW', 'mWH', 'mWR', 'mWB')] <- c(1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF))/4
  tMatrix['fWH','mWR', c('fWW', 'fWH', 'fWR', 'fWB',
                         'fHR', 'fRR', 'fRB','mWW', 'mWH', 'mWR', 'mWB',
                         'mHR', 'mRR', 'mRB')] <- c(1-cF, 1+cF*chF, cF*(1-chF)*crF + 1-cF, cF*(1-chF)*(1-crF),
                                                    1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF))/8
  tMatrix['fWH','mWB', c('fWW', 'fWH', 'fWR', 'fWB',
                         'fHB', 'fRB', 'fBB','mWW', 'mWH', 'mWR', 'mWB',
                         'mHB', 'mRB', 'mBB')] <- c(1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF) + 1-cF,
                                                    1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF))/8
  tMatrix['fWH','mHH', c('fWH', 'fHH',
                         'fHR', 'fHB','mWH', 'mHH',
                         'mHR', 'mHB')] <- c((1-cF), 1+cF*chF,
                                             cF*(1-chF)*crF, cF*(1-chF)*(1-crF))/4
  tMatrix['fWH','mHR', c('fWH', 'fHH', 'fHR', 'fHB',
                         'fWR', 'fRR', 'fRB','mWH', 'mHH', 'mHR', 'mHB',
                         'mWR', 'mRR', 'mRB')] <- c((1-cF), 1+cF*chF, cF*(1-chF)*crF + 1+cF*chF,
                                                    cF*(1-chF)*(1-crF), (1-cF),
                                                    cF*(1-chF)*crF, cF*(1-chF)*(1-crF))/8
  tMatrix['fWH','mHB', c('fWH', 'fHH', 'fHR', 'fHB',
                         'fWB', 'fRB', 'fBB','mWH', 'mHH', 'mHR', 'mHB',
                         'mWB', 'mRB', 'mBB')] <- c((1-cF), 1+cF*chF, cF*(1-chF)*crF,
                                                    cF*(1-chF)*(1-crF) + 1+cF*chF, (1-cF),
                                                    cF*(1-chF)*crF, cF*(1-chF)*(1-crF))/8
  tMatrix['fWH','mRR', c('fWR', 'fHR', 'fRR', 'fRB','mWR', 'mHR', 'mRR', 'mRB')] <- c(1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF))/4
  tMatrix['fWH','mRB', c('fWR', 'fHR', 'fRR', 'fRB',
                         'fWB', 'fHB', 'fBB','mWR', 'mHR', 'mRR', 'mRB',
                         'mWB', 'mHB', 'mBB')] <- c(1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF) + cF*(1-chF)*crF,
                                                    1-cF, 1+cF*chF, cF*(1-chF)*(1-crF))/8
  tMatrix['fWH','mBB', c('fWB', 'fHB', 'fRB', 'fBB','mWB', 'mHB', 'mRB', 'mBB')] <- c(1-cF, 1+cF*chF, cF*(1-chF)*crF, cF*(1-chF)*(1-crF))/4
  
  
  
  #protection from underflow errors
  tMatrix[tMatrix < .Machine$double.eps] <- 0
  
  ## initialize viability mask. No mother/father-specific death, so use basic mask
  
  viabilityMask <- array(data = 1L, dim = c(size,size,size), dimnames = list(gtype, gtype, gtype))
  
  modifiers = cubeModifiers(gtype, eta = eta, phi = phi, omega = omega, xiF = xiF, xiM = xiM, s = s)
  
  ## put everytWing into a labeled list to return
  return(list(
    ih = tMatrix,
    tau = viabilityMask,
    genotypesID = gtype,
    genotypesN = size,
    wildType = c("mWW","fWW"),
    eta = modifiers$eta,
    phi = modifiers$phi,
    omega = modifiers$omega,
    xiF = modifiers$xiF,
    xiM = modifiers$xiM,
    s = modifiers$s,
    releaseType = "mHH"
  ))
  
  
}



  
  
  
  
  
  
  
  