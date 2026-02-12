###### 0. Settings ######
##### 0.1. Remove objects #####
rm(list = ls())
gc(reset = TRUE);gc(reset = TRUE)

##### 0.2. Load Packages and Functions #####
require(stringr)
require(dplyr)
require(plotly)
require(htmlwidgets)
require(rPref)

##### 0.3. Define Parameters #####
scriptID <- "2.0.1."

nPerTop <- 0.01 ### top 0.01 percent individual's GV instead of variance

traitNames <- c("Perillaldehyde", "RosmarinicAcid")

##### 0.4. Define Directories #####
crsCombDir <- "data/crossComb/"
haploDir <- "data/genotype/"
mrkEffDir <- "data/markerEffects/"
progGvDir <- "result/1.0.simuProgVar/"

saveDir <- paste0("result/2.0.Pareto_progGV_simu/")
if (!dir.exists(saveDir)) {
  dir.create(saveDir, recursive = TRUE)
}


###### 1. Load Data ######
##### 1.1. Genotype #####
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

##### 1.3. GV by top 1% #####
fileNameGvTop <- paste0(progGvDir, "1.0.1.progGv_top_1%_nComb_196.rds")
progGvTop <- readRDS(file = fileNameGvTop)

##### 1.4. Intercept at F4 #####
muPeril27 <- 4.207890
muRos27 <- 23.48521
muPeril40 <- 2.795166
muRos40 <- 25.12886
muMat <- rbind(c(muPeril27, muRos27, 
               c(muPeril40, muRos40)))
rownames(muMat) <- c("S827", "S840")
colnames(muMat) <- traitNames


###### 2. Draw Pareto (rPref) ######
### Inter- or Intra- cross
crsPatternVec <- apply(X = crsCombDf, MARGIN = 1, FUN = function(x) {
  p1Pop <- str_sub(x[1], 1, 4)
  p2Pop <- str_sub(x[2], 1, 4)
  if (p1Pop == p2Pop) {
    return(paste0(p1Pop, "_Within"))
  } else {
    return("Across")
  }
})

plotDf <- data.frame(progGvTop[, traitNames], Cross = crsPatternVec)
plotDf$Cross <- factor(plotDf$Cross)
plotDf$label <- paste0("ID: ", rownames(plotDf))

### GV + intercept
plotDf[, 1:2] <- do.call(what = rbind, 
                         args = lapply(X = 1:nrow(plotDf), FUN = function(indNo) {
                           if (plotDf[indNo, 3] == "S827_Within") {
                             return(plotDf[indNo, 1:2] + muMat[1, ])
                           } else if (plotDf[indNo, 3] == "S840_Within") {
                             return(plotDf[indNo, 1:2] + muMat[2, ])
                           } else {
                             return(plotDf[indNo, 1:2] + (muMat[1, ] + muMat[2, ]) / 2)
                           }
                         }))

### pareto point
paretoPt27 <- psel(df = plotDf[plotDf$Cross == "S827_Within", ], 
                   pref = high(RosmarinicAcid) * high(Perillaldehyde), top_level = 1)
paretoPt40 <- psel(df = plotDf[plotDf$Cross == "S840_Within", ], 
                   pref = high(RosmarinicAcid) * high(Perillaldehyde), top_level = 1)
paretoPtAcross <- psel(df = plotDf[plotDf$Cross == "Across", ], 
                       pref = high(RosmarinicAcid) * high(Perillaldehyde), top_level = 1)

### plot
fileNamePlotPng <- paste0(saveDir, scriptID, "ParetoPoint_progGv_simulation_top_",
                          nPerTop * 100, "_Per_individual.png")
png(filename = fileNamePlotPng, width = 1200, height = 800)
p <- ggplot(data = plotDf, aes(x = Perillaldehyde, y = RosmarinicAcid, text = label)) +
  geom_point(color = "gray", size = 6) +
  geom_point(data = paretoPt27, aes(x = Perillaldehyde, y = RosmarinicAcid),
             color = "skyblue", size = 8) +
  geom_point(data = paretoPt40, aes(x = Perillaldehyde, y = RosmarinicAcid),
             color = "pink", size = 8) +
  geom_point(data = paretoPtAcross, aes(x = Perillaldehyde, y = RosmarinicAcid),
             color = "darkgreen", size = 8) +
  theme_bw() +
  theme(text = element_text(size = 24))
print(p)
dev.off()
pNew <- ggplotly(p)

### save
fileNamePlotly <- paste0(saveDir, scriptID, "ParetoPoint_progGv_simulation_top_", 
                         nPerTop * 100, "_Per_individual.html")
saveWidget(widget = pNew, file = fileNamePlotly, selfcontained = TRUE)

