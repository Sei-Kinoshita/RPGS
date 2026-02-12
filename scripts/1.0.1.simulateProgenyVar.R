###### 0. Settings ######
##### 0.1. Remove objects #####
rm(list = ls())
gc(reset = TRUE);gc(reset = TRUE)

##### 0.2. Load Packages and Functions #####
source("scripts/function_scripts/1.0.1.function.simuProgeny.r")

##### 0.3. Define Parameters #####
scriptID <- "1.0.1."

nProg1 <- 1
nProg2 <- 4
nProg3 <- 500

nSimu <- 500
nCores <- 30

nPerTop <- 0.01 ### top 0.01 percent individual's GV instead of variance

##### 0.4. Define Directories #####
crsCombDir <- "data/crossComb/"
haploDir <- "data/genotype/"
mrkEffDir <- "data/markerEffects/"
linkageDir <- "data/linkageInfo/"

saveDir <- paste0("result/1.0.simuProgVar/")
if (!dir.exists(saveDir)) {
  dir.create(saveDir, recursive = TRUE)
}


###### 1. Load Data ######
##### 1.1. Haplotype & Genotype Data #####
fileNameHaplo <- paste0(haploDir, "0.2.F4895_Haplotype_4ParentsArray.rds")
f4HaploArray <- readRDS(file = fileNameHaplo)

fileNameGeno <- paste0(haploDir, "0.1.F4895_Genotype.csv")
f4GenoMat <- as.matrix(read.csv(file = fileNameGeno, header = TRUE, row.names = 1))

##### 1.2. Selected Cross Combination #####
fileNameSelCrs <- paste0(crsCombDir, "0.1.F4_crsComb_selected_nComb_1000.rds")
crsCombDf <- readRDS(file = fileNameSelCrs)
### cross with high anthocyanin
isHighCrs <- apply(X = crsCombDf, MARGIN = 1, FUN = function(x) {
  (f4GenoMat[x[1], "v1055"] + f4GenoMat[x[2], "v1055"]) == -2 
})
crsCombDf <- crsCombDf[isHighCrs, ]
nComb <- nrow(crsCombDf)
selIndName <- unique(unlist(crsCombDf))
combName <- apply(X = crsCombDf, MARGIN = 1, FUN = paste0, collapse = " X ")

##### 1.3. Marker Effects Data #####
fileNameMrkEff <- paste0(mrkEffDir, "0.2.F4_allTraits_bindedMrkEffLst.rds")
mrkEffAllTraitLst <- readRDS(file = fileNameMrkEff)[1:2]

nTrait <- length(mrkEffAllTraitLst)
traitNames <- names(mrkEffAllTraitLst)
### GV of F4
genoF4 <- f4HaploArray[, -1, , 1] + f4HaploArray[, -1, , 2] ### remove sekiho
gvF4Mat <- do.call(what = cbind, 
                   args = lapply(X = mrkEffAllTraitLst, FUN = function(oneTraitMat) {
                     apply(X = genoF4, MARGIN = 3, FUN = function(x) sum(x * oneTraitMat))
                   }))

##### 1.4. Linkage Information #####
### marker name in each chromosome
fileNameMrkNameEachChr <- paste0(linkageDir, "0.2.F4_mrkName_eachChromLst.rds")
mrkNameEachChrLst <- readRDS(file = fileNameMrkNameEachChr)
### recombination rate in each chromosome
fileNameRREachChr <- paste0(linkageDir, "0.2.F4_rr_eachChromLst.rds")
rrEachChrLst <- readRDS(file = fileNameRREachChr)




###### 2. Simulate Genetic Variance of Progeny ######
fileNameProgGv <- paste0(saveDir, scriptID, "progGvAllSimuLst_nSimu_", nSimu, 
                         "_nComb_", nComb, ".rds")
if (!file.exists(fileNameProgGv)) {
  
  gvAllCrsLst <- rep(list(NULL), nComb)
  
  for (crsCombNo in 1:nComb) {
    # crsCombNo <- 1
    # tictoc::tic()
    
    cat("------ Cross Combination No.", crsCombNo, "/", nComb, " will be performed ------\n")
    
    sireNameNow <- crsCombDf$sire[crsCombNo]
    damNameNow <- crsCombDf$dam[crsCombNo]
    
    ### simulation
    gvAllSimuLst <- pbmcapply::pbmclapply(X = 1:nSimu, FUN = function(simuNo) {
      matingChain4Parents(haploArrayP1 = f4HaploArray[, , sireNameNow, ],
                          haploArrayP2 = f4HaploArray[, , damNameNow, ],
                          mrkEffLst = mrkEffAllTraitLst, 
                          nProg1 = nProg1, 
                          nProg2 = nProg2, 
                          nProg3 = nProg3)
      
    }, mc.cores = nCores)
    names(gvAllSimuLst) <- 1:nSimu
    
    gvAllCrsLst[[crsCombNo]] <- gvAllSimuLst
    
    gc(reset = TRUE);gc(reset = TRUE)
    # tictoc::toc()
  }
  
  ### save
  saveRDS(object = gvAllCrsLst, file = fileNameProgGv)
  
} else {
  gvAllCrsLst <- readRDS(file = fileNameProgGv)
}

### top 1% individual's GV
fileNameGvTop <- paste0(saveDir, scriptID, "progGv_top_", nPerTop * 100, 
                        "%_nComb_", nComb, ".rds")
if (!file.exists(fileNameGvTop)) {
  topGvAllCrs <- do.call(what = rbind, 
                         args = lapply(X = gvAllCrsLst, FUN = function(oneCrsLst) {
                           
                           allSimuTopGv <- sapply(X = oneCrsLst, FUN = function(oneSimuGvMat) {
                             apply(X = oneSimuGvMat, MARGIN = 2, FUN = function(x) {
                               sort(x, decreasing = TRUE)[nrow(oneSimuGvMat) * nPerTop]
                             })
                           })
                           oneCrsTopGv <- apply(X = allSimuTopGv, MARGIN = 1, 
                                                FUN = mean)
                           
                           return(oneCrsTopGv)
                         }))
  rownames(topGvAllCrs) <- combName
  ### save
  saveRDS(object = topGvAllCrs, file = fileNameGvTop)
} else {
  topGvAllCrs <- readRDS(fileNameGvTop)
}