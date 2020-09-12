#' Function by: Ethan A. Brown (March 4 2020)
#'
#'
#' Consolidates output of MGDrivE into a conditional probability table for use in a Bayesian network
#'
#'
#'
#' @param out_path. complete file path to the output folder (where the simulation directories are located)
#' @param tPoints specific time points that you are interested in (vector)
#' @param cube. inheritance cube object (same one you used to run the model)
#' @param nrep. number of simulations
#' @param sitesNumber. number of patches
#'
#' @export

MGDcpt = function(out_path. = out_path, #complete file path to the output folder (where the simulation directories are located)
                  tPoints = 1L,
                  cube. = cube, #inheritance cube object (same one you use to run the model)
                  nrep. = nRep, #number of simulations
                  sitesNumber. = sitesNumber)#number of patches
{


  #make sure to look at your 'outFolder' and 'out_path' objects. make sure all "/" are present


  list.dirs <- function(path=".", pattern=NULL, all.dirs=FALSE,
                        full.names=FALSE, ignore.case=FALSE) {
    # use full.names=TRUE to pass to file.info
    all <- list.files(path, pattern, all.dirs,
                      full.names=TRUE, recursive=FALSE, ignore.case)
    dirs <- all[file.info(all)$isdir]
    # determine whether to return full names or just dir names
    if(isTRUE(full.names))
      return(dirs)
    else
      return(basename(dirs))
  }

  sims = list.dirs(out_path.)

  bn_csvs = list.files(as.character(paste(out_path.,"/", sims, sep = "")))

  def = getwd()

  setwd(as.character(paste(out_path., sep = "")))

  directory  = paste("./",rep(sims, each = sitesNumber., times = 2),"/", bn_csvs, sep = "")

  frames = list()

  for (i in seq_along(directory)){
    frames[[i]] = read.csv(directory[i])[tPoints,]
  }

  ngtype = cube.[["genotypesN"]]

  named_frames = c()
  for (i in 1:length(frames)){
    named_frames[[i]] = frames[[i]][,2:(ngtype+1)]
    named_frames[[i]]["names"] <- bn_csvs[i]

  }

  run_index = c()
  for (i in 1:nrep.){
    pattern = paste("Run",sims,sep = "")
    run_index[[i]] = grep(pattern[i], named_frames)
  }

  run_index = BBmisc::dapply(run_index, as.vector)
  run_index2 = rep(1:nrep., each = sitesNumber., times = 2)

  op = split(named_frames, run_index2) #split genotype frequencies by run

  op2 = lapply(op, function(x) Reduce(`+`, lapply(x, function(y) y[-ncol(y)])))

  tots = list()
  freqs = list()
  for(i in 1:length(op2)){
    tots[[i]] = apply(op2[[i]][,], MARGIN = 1, sum)
    freqs[[i]] = cbind(op2[[i]],tots[[i]])
    colnames(freqs[[i]])[(1+ngtype)] <- "Total"
  }

  cpt = Reduce(`+`, lapply(freqs, `==`, 0))/length(freqs)*100


  return(list(cpt, freqs))


  setwd(def)
}

########

