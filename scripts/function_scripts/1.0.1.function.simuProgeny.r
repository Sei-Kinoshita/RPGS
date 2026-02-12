# haploArrayP1 <- f4HaploArray[, , sireNameNow, ]
# haploArrayP2 <- f4HaploArray[, , damNameNow, ]
# mrkEffLst <- mrkEffAllTraitLst
# nProg1 <- 1
# nProg2 <- 4
# nProg3 <- 500

###' @param haplo4PArrayEachChrom [array] mrkNoEachChrom x parentsNo x 2(Haplo) x indNo
###' @param rec [vector] recombination rate in 1 chromosome
###'
###### Function to simulate gametes ######
simulateGametes4Parents <- function(haplo4PArrayEachChrom, rec) {
  
  samples <- runif(length(rec), min = 0, max = 1)
  ### Decide cross-over occurs or not at each marker
  crossOver <- ((rec - samples) >= 0)
  ### Count times of cross-over
  selectHaplo <- cumsum(crossOver) %% 2
  selectHaplo <- selectHaplo + 1
  
  ### Select which haplotype to be gamete
  whichHaplo <- (runif(1) < 0.5) + 1
  
  gamHaplo <- haplo4PArrayEachChrom[, , whichHaplo]
  
  ### Change genotype of where cross-over happened
  gamHaplo[selectHaplo == 2, ] <- haplo4PArrayEachChrom[selectHaplo == 2, , 3 - whichHaplo]
  
  return(gamHaplo)
}




###' @param nProg [integer] number of progeny one individual generates
###' @param haploArray [array] haploType of 1 individual, mrkNo x parentsNo x 2(Haplo) (3 dimension)
###' @return progeny number length list of gametes
###' 
##### Function to generate multiple gametes in one individual #####
gametesOneInd4Parents <- function(haploArray, nProg = 1) {
  
  nChrom <- length(mrkNameEachChrLst)
  
  gametesLst <- sapply(X = 1:nProg, FUN = function(progNo) {
    
    gameteOneIndMat <- do.call(what = rbind, 
                               args = sapply(X = 1:nChrom, FUN = function(chrNo) {
                                 
                                 mrkNameEachChr <- mrkNameEachChrLst[[chrNo]]
                                 rrEachChr <- rrEachChrLst[[chrNo]]
                                 
                                 haploEachInd <- haploArray[mrkNameEachChr, , ]
                                 
                                 gameteHaplo <- simulateGametes4Parents(haplo4PArrayEachChrom = haploEachInd, 
                                                                        rec = rrEachChr)
                                 
                                 return(gameteHaplo)
                                 
                               }, simplify = FALSE))
    
    return(gameteOneIndMat)
    
  }, simplify = FALSE)
  names(gametesLst) <- paste0("gamete_", 1:nProg)
  
  return(gametesLst)
}




##### Function of mating #####
###' @param haploArrayP1 [array] haploType of mating candidates (Parent1), mrkNo x parentsNo x 2(Haplo) (3 dimension)
###' @param haploArrayP2 [array] haploType of mating candidates (Parent2), mrkNo x parentsNo x 2(Haplo) (3 dimension)
###' @param nProg [integer] numeber of progeny, same as the gametes of the parents
###' 
mate2Inds4Parents <- function(haploArrayP1, haploArrayP2, nProg = 1) {
  
  ### gametes of each parent
  gameteP1 <- gametesOneInd4Parents(haploArray = haploArrayP1, nProg = nProg)
  gameteP2 <- gametesOneInd4Parents(haploArray = haploArrayP2, nProg = nProg)
  
  ### progeny number length list
  progHaploLst <- lapply(X = 1:nProg, FUN = function(progNo) {
    progHaploArray <- array(data = NA, 
                            dim = c(dim(gameteP1[[progNo]]), 2, 1), 
                            dimnames = c(dimnames(gameteP1[[progNo]]), 
                                         list(paste0("Ploidy", 1:2))))
    progHaploArray[, , 1, 1] <- gameteP1[[progNo]]
    progHaploArray[, , 2, 1] <- gameteP2[[progNo]]
    
    return(progHaploArray)
  })
  allProgHaploArray <- do.call(what = abind::abind, args = progHaploLst)
  dimnames(allProgHaploArray)[[4]] <- paste0("Progeny_", 1:nProg)
  
  return(allProgHaploArray)
}




##### Function of selfing #####
###' @param haploArray [array] haploType of 1 individual, mrkNo x parentsNo x 2(Haplo) (3 dimension)
###' 
selfing4Parents <- function(haploArray, nProg = 1) {
  
  gameteSelfing1 <- gametesOneInd4Parents(haploArray = haploArray, nProg = nProg)
  gameteSelfing2 <- gametesOneInd4Parents(haploArray = haploArray, nProg = nProg)
  
  if (nProg == 1) {
    ### if only 1 progeny, return haplotype array
    progSelfingHaploArray <- array(data = NA, 
                                   dim = c(dim(gameteSelfing1[[1]]), 2), 
                                   dimnames = c(dimnames(gameteSelfing1[[1]]), 
                                                list(paste0("Ploidy", 1:2))))
    progSelfingHaploArray[, , 1] <- gameteSelfing1[[1]]
    progSelfingHaploArray[, , 2] <- gameteSelfing2[[1]]
    
    return(progSelfingHaploArray)
    
  } else {
    
    ### if more than 1 progeny, return list
    progSelfingHaploLst <- lapply(X = 1:nProg, FUN = function(progNo) {
      progSelfingHaploArray <- array(data = NA, 
                                     dim = c(dim(gameteSelfing1[[progNo]]), 2, 1), 
                                     dimnames = c(dimnames(gameteSelfing1[[progNo]]), 
                                                  list(paste0("Ploidy", 1:2))))
      progSelfingHaploArray[, , 1, 1] <- gameteSelfing1[[progNo]]
      progSelfingHaploArray[, , 2, 1] <- gameteSelfing2[[progNo]]
      
      return(progSelfingHaploArray)
    })
    
    multiProgSelfingHaploArray <- do.call(what = abind::abind, 
                                          args = progSelfingHaploLst)
    dimnames(multiProgSelfingHaploArray)[[4]] <- paste0("Progeny_", 1:nProg)
    
    return(multiProgSelfingHaploArray)
  }
}




##### Function of mating chain between 1 mating pair #####
### mating at first generation and SSD for n generations
###' @param haploArrayP1 [array] haploType of mating candidates (Parent1), nMrk x nParents x 2(Haplo) (3 dimension)
###' @param haploArrayP2 [array] haploType of mating candidates (Parent2), nMrk x nParents x 2(Haplo) (3 dimension)
###' @param mrkEffLst [list] marker effects, list of nTrait, each contains matrix (nMrk x nParents - 1)
###' @param nProg1 [integer] number of progeny from F4->F5
###' @param nProg2 [integer] number of progeny from F5->G1(cross)
###' @param nProg3 [integer] number of progeny from G1->G2
###' 
##### Function of mating chain between 1 mating pair #####
matingChain4Parents <- function(haploArrayP1, haploArrayP2, mrkEffLst, 
                                nProg1 = 1, nProg2 = 4, nProg3 = 500) {
  
  ### selfing from F4 to F5
  f5HaploInd1 <- selfing4Parents(haploArray = haploArrayP1, nProg = nProg1)
  f5HaploInd2 <- selfing4Parents(haploArray = haploArrayP2, nProg = nProg1)
  
  ### mating to make F1 from F5
  f1HaploAllProgLst <- lapply(X = 1:nProg1, FUN = function(prog1No) {
    f1HaploLst <- mate2Inds4Parents(haploArrayP1 = f5HaploInd1, 
                                    haploArrayP2 = f5HaploInd2, 
                                    nProg = nProg2)
  })
  
  f2HaploAllProgLst <- lapply(X = 1:nProg1, FUN = function(prog1No) {
    
    f2HaploLst <- lapply(X = 1:nProg2, FUN = function(prog2No) {
      selfing4Parents(haploArray = f1HaploAllProgLst[[prog1No]][, , , prog2No], 
                      nProg = nProg3)
    })
    
    return(f2HaploLst)
  })
  
  f2GvVarLst <- lapply(X = 1:nProg1, FUN = function(prog1No) {
    gvProg2Mat <- do.call(what = rbind, 
                          args = lapply(X = 1:nProg2, FUN = function(prog2No) {
                            
                            
                            f2HaploArrayAllProg <- f2HaploAllProgLst[[prog1No]][[prog2No]]
                            f2GenoMatAllProg <- f2HaploArrayAllProg[, , 1, ] + f2HaploArrayAllProg[, , 2, ]
                            
                            ### calculate genotypic value of each progeny
                            gvF2AllTraitMat <- do.call(what = cbind, 
                                                       args = lapply(X = mrkEffLst, FUN = function(mrkEffMat) {
                                                         apply(X = f2GenoMatAllProg[, -1, ], MARGIN = 3, 
                                                               FUN = function(x) sum(x * mrkEffMat))
                                                       }))
                            
                            return(gvF2AllTraitMat)
                          }))
    
    return(gvProg2Mat)
  })
  
  if (nProg1 == 1) {
    
    return(f2GvVarLst[[1]])
    
  } else {
    
    return(f2GvVarLst)
    
  }
}