###### 0. Settings ######
##### 0.1. Remove objects #####
rm(list = ls())
gc(reset = TRUE);gc(reset = TRUE)

##### 0.2. Load Packages and Functions #####
require(stringr)

##### 0.3. Define Parameters #####
scriptID <- "0.1."

popNames <- c("S827", "S840", "S844")
nSelCrs <- 1000 ### number of selected pairs with high gv mean

##### 0.4. Define Directories #####
genoDir <- "data/genotype/"
mrkEffDir <- "data/markerEffects/"

saveDir <- "data/crossComb/"
if (!dir.exists(saveDir)) {
  dir.create(saveDir, recursive = TRUE)
}


###### 1. Load Data ######
### genotype
fileNameHaplo <- paste0(genoDir, "0.2.F4895_Haplotype_4ParentsArray.rds")
haploArrayF4 <- readRDS(file = fileNameHaplo)
genoArrayF4 <- haploArrayF4[, , , 1] + haploArrayF4[, , , 2]
nF4 <- dim(genoArrayF4)[3]
f4Names <- dimnames(genoArrayF4)[[3]]

### marker effects
fileNameMrkEff <- paste0(mrkEffDir, "0.2.F4_allTraits_bindedMrkEffLst.rds")
mrkEffLst <- readRDS(file = fileNameMrkEff)


###### 2. Calculate GV ######
gvF4 <- sapply(X = 1:2, FUN = function(traitNo) {
  apply(X = genoArrayF4[, -1, ], MARGIN = 3, FUN = function(x) {
    sum(x * mrkEffLst[[traitNo]])
  })
})
colnames(gvF4) <- names(mrkEffLst)[1:2]

gvScaledF4 <- apply(X = gvF4, MARGIN = 2, FUN = scale)

selIndexF4 <- apply(X = gvScaledF4, MARGIN = 1, FUN = sum)
names(selIndexF4) <- f4Names

### rank GV
selIndexSorted <- sort(x = selIndexF4, decreasing = TRUE)


###### 3. Make Cross Pair ######
crossGrid <- expand.grid(rep(list(1:nF4), 2))
crsCombNoGrid <- crossGrid[crossGrid[, 1] <= crossGrid[, 2], ]
nComb <- nrow(crsCombNoGrid)
### df with line name
crsCombDf <- data.frame(sire = f4Names[crsCombNoGrid[, 1]], 
           dam = f4Names[crsCombNoGrid[, 2]])
### GV for each cross pair
selIndexCrsMat <- as.matrix(data.frame(sire = selIndexF4[crsCombDf$sire], 
                            dam = selIndexF4[crsCombDf$dam]))
rownames(selIndexCrsMat) <- 1:nComb

gvMeanCrsVec <- apply(X = selIndexCrsMat, MARGIN = 1, FUN = function(x) sum(x) / 2)
gvMeanCrsSorted <- sort(x = gvMeanCrsVec, decreasing = TRUE)

### select cross pairs with high GV mean
selCrsCombDf <- crsCombDf[as.numeric(names(gvMeanCrsSorted))[1:nSelCrs], ]


###### 3. Save ######
fileNameAllCrsComb <- paste0(saveDir, scriptID, "F4_crsComb_all_nComb_", nComb, ".rds")
saveRDS(object = crsCombDf, file = fileNameAllCrsComb)
fileNameSelCrsComb <- paste0(saveDir, scriptID, "F4_crsComb_selected_nComb_", nSelCrs, ".rds")
saveRDS(object = selCrsCombDf, file = fileNameSelCrsComb)
