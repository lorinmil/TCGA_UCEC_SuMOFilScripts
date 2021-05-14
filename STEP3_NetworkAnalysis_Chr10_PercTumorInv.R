library(TCGAbiolinks)
library(SummarizedExperiment)
library(GenomicRanges)
library(stringr)
library(parallel)
library(PMA)

#Must be downloaded from GibHub
library(SuMOFil)


#---Load TCGA Data---#
#Load the data from Step 2 (update location)
load("DATA LOCATION HERE/UCEC_TCGA.Rdata")



#---Drill down to the data for analysis---#

#Pull out the outcome of interest
y <- ucec.clinical[keepIDs, 'pct_tumor_invasion']

#Filter down to just the chromosomes of interest
methy.measure <- methy.measure[methyChroms %in% c('chr10'), ]
trans.measure <- trans.measure[transChroms %in% c('chr10'), ]

#Clean up NA values and items with only one value (e.g. all zeros)
nonEmptyOutcomes <- which(!is.na(y) & y<=100 & y>=0)
y <- y[nonEmptyOutcomes]
trans.measure <- trans.measure[, keepIDs[nonEmptyOutcomes]]
methy.measure <- methy.measure[, keepIDs[nonEmptyOutcomes]]
keepTrans <- apply(trans.measure, 1, function(z) ifelse(sum(is.na(z))>0 || length(unique(z))<=5, FALSE, TRUE))
trans.measure <- trans.measure[keepTrans, ]
keepMethyl <- apply(methy.measure, 1, function(z) ifelse(sum(is.na(z))>0 || length(unique(z))<=5, FALSE, TRUE))
methy.measure <- methy.measure[keepMethyl, ]

#Restructure the matrices
trans.measure <- t(trans.measure)
methy.measure <- t(methy.measure)

#Scale the data
methy.measure_unscaled <- methy.measure
trans.measure_unscaled <- trans.measure
colnames_methy <- colnames(methy.measure)
colnames_trans <- colnames(trans.measure)
methy.measure <- scale(methy.measure)
colnames(methy.measure) <- colnames_methy
trans.measure <- scale(trans.measure)
colnames(trans.measure) <- colnames_trans

#Filter down the gene info
trans.geneInfo <- trans.geneInfo[trans.geneInfo$geneID %in% colnames(trans.measure),]
methy.geneInfo <- methy.geneInfo[methy.geneInfo$geneID %in% colnames(methy.measure),]

#Only keep necessary items (for memory purposes)
rm(list=base::setdiff(ls(), c("trans.measure", "methy.measure", "y", "trans.geneInfo", "methy.geneInfo",
                              "trans.measure_unscaled", "methy.measure_unscaled")))


#(Optional) Save this data, only run this once
save.image("DATA LOCATION HERE/TCGAdata_Chr10_PercTumorInv.Rdata")

#Load the data (if the above has already been ran)
load("DATA LOCATION HERE/TCGAdata_Chr10_PercTumorInv.Rdata")





#----------------------Perform methods before filtering----------------#

#Perform CCA from tibshirani
time1 <- Sys.time()
CCA_PMA_tune <- PMA::CCA.permute(x=trans.measure,
                            z=methy.measure,
                            typex="standard",
                            typez="standard",
                            outcome="quantitative",
                            y=y)
CCA_PMA <- PMA::CCA(x=trans.measure,
                    z=methy.measure,
                    typex="standard",
                    typez="standard",
                    outcome="quantitative",
                    y=y,
                    penaltyx=CCA_PMA_tune$bestpenaltyx,
                    penaltyz=CCA_PMA_tune$bestpenaltyz)
time2 <- Sys.time()
timeB4Filter_CCA_PMA <- time2 - time1
CCA_PMA$runTime <- timeB4Filter_CCA_PMA

#(Optional) Save outputs
saveRDS(CCA_PMA, "DATA LOCATION HERE/CCA_PMA_before.rda")

#Read data (if above has already been saved so don't have to re-execute)
CCA_PMA <- readRDS("DATA LOCATION HERE/CCA_PMA_before.rda")
selectedPMAx_before <- which(CCA_PMA$u!=0)
selectedPMAg_before <- which(CCA_PMA$v!=0)










#---Perform SuMO-Fil---#
time1 <- Sys.time()
filterResults <- SuMOFil(x=trans.measure,
                         g=methy.measure,
                         y=y,
                         numClusters_1=3,
                         numClusters_2=3)
time2 <- Sys.time()
timeSuMOFil <- time2 - time1


#(optional) Save outputs
saveRDS(filterResults, "DATA LOCATION HERE/filterResults_SuMOFil.rda")
filterResults <- readRDS("DATA LOCATION HERE/filterResults_SuMOFil.rda")


#Perform CCA from tibshirani
time1 <- Sys.time()
CCA_PMA_SuMOFiltune <- CCA.permute(x=trans.measure[,-as.numeric(filterResults$removeX_both)],
                            z=methy.measure[,-as.numeric(filterResults$removeG_both)],
                            typex="standard",
                            typez="standard",
                            outcome="quantitative",
                            y=y)
CCA_PMASuMOFil <- PMA::CCA(x=trans.measure[,-as.numeric(filterResults$removeX_both)],
                    z=methy.measure[,-as.numeric(filterResults$removeG_both)],
                    typex="standard",
                    typez="standard",
                    outcome="quantitative",
                    y=y,
                    penaltyx=CCA_PMA_SuMOFiltune$bestpenaltyx,
                    penaltyz=CCA_PMA_SuMOFiltune$bestpenaltyz)
time2 <- Sys.time()
timeSuMOFil_CCA_PMA <- time2 - time1
CCA_PMASuMOFil$runTime <- timeSuMOFil_CCA_PMA

#(Optional) Save outputs
saveRDS(CCA_PMASuMOFil, "DATA LOCATION HERE/chr10_PercTumorInv/CCA_PMASuMOFil.rda")

#Read outputs (only if it wasn't reran)
CCA_PMASuMOFil <- readRDS("DATA LOCATION HERE/CCA_PMASuMOFil.rda")

#What was selected after the filter?
selectedPMAx_SuMOFIL <- (1:ncol(trans.measure))[-as.numeric(filterResults$removeX_both)][which(CCA_PMASuMOFil$u!=0)]
selectedPMAg_SuMOFIL <- (1:ncol(methy.measure))[-as.numeric(filterResults$removeG_both)][which(CCA_PMASuMOFil$v!=0)]
#How many were filtered that were previously selected?
sum(selectedPMAx_before %in% as.numeric(filterResults$removeX_both))
sum(selectedPMAg_before %in% as.numeric(filterResults$removeG_both))
#How many genes maintained selection after filtering?
sum(selectedPMAx_SuMOFIL %in% selectedPMAx_before)
sum(selectedPMAg_SuMOFIL %in% selectedPMAg_before)








#---Perform the mean filtering---#
time1 <- Sys.time()
xMean <- as.numeric(colMeans(trans.measure_unscaled))
gMean <- as.numeric(colMeans(methy.measure_unscaled))
#Cluster the means
xMeanClust <- kmeans(xMean, centers=9)
gMeanClust <- kmeans(gMean, centers=9)
#Filter the X features
clusterMeansX <- xMeanClust$centers
smallestX <- which(clusterMeansX == min(clusterMeansX))
removeX_mean <- rep("Keep", ncol(trans.measure))
removeX_mean[xMeanClust$cluster == smallestX] <- "Remove"
#Filter the G features
clusterMeansG <- gMeanClust$centers
smallestG <- which(clusterMeansG == min(clusterMeansG))
removeG_mean <- rep("Keep", ncol(methy.measure))
removeG_mean[gMeanClust$cluster == smallestG] <- "Remove"
time2 <- Sys.time()
time2 - time1

#How many were removed?
sum(removeX_mean=="Remove")
sum(removeG_mean=="Remove")



#Perform CCA from tibshirani
time1 <- Sys.time()
CCA_PMA_LowMeantune <- CCA.permute(x=trans.measure[,as.numeric(which(removeX_mean=="Keep"))],
                                   z=methy.measure[,as.numeric(which(removeG_mean=="Keep"))],
                                   typex="standard",
                                   typez="standard",
                                   outcome="quantitative",
                                   y=y)
CCA_PMALowMean <- PMA::CCA(x=trans.measure[,as.numeric(which(removeX_mean=="Keep"))],
                           z=methy.measure[,as.numeric(which(removeG_mean=="Keep"))],
                           typex="standard",
                           typez="standard",
                           outcome="quantitative",
                           y=y,
                           penaltyx=CCA_PMA_LowMeantune$bestpenaltyx,
                           penaltyz=CCA_PMA_LowMeantune$bestpenaltyz)
time2 <- Sys.time()
timeLowMean_CCA_PMA <- time2 - time1
CCA_PMALowMean$runTime <- timeLowMean_CCA_PMA

#(Optional) Save outputs
saveRDS(CCA_PMALowMean, "DATA LOCATION HERE/CCA_PMALowMean.rda")

#Read the outputs (only if it hasn't been rerun)
CCA_PMALowMean <- readRDS("DATA LOCATION HERE/CCA_PMALowMean.rda")

#What was selected after the filter?
selectedPMAx_LowMean <- (1:ncol(trans.measure))[as.numeric(which(removeX_mean=="Keep"))][which(CCA_PMALowMean$u!=0)]
selectedPMAg_LowMean <- (1:ncol(methy.measure))[as.numeric(which(removeG_mean=="Keep"))][which(CCA_PMALowMean$v!=0)]
#How many were filtered that were previously selected?
sum(selectedPMAx_before %in% as.numeric(which(removeX_mean!="Keep")))
sum(selectedPMAg_before %in% as.numeric(which(removeG_mean!="Keep")))
#How many genes maintained selection after filtering?
sum(selectedPMAx_LowMean %in% selectedPMAx_before)
sum(selectedPMAg_LowMean %in% selectedPMAg_before)









#---Perform the variance filtering---#
time1 <- Sys.time()
#Fetch the feature variances
xVar <- resample::colVars(trans.measure_unscaled)
gVar <- resample::colVars(methy.measure_unscaled)
#Cluster the variances
xVarClust <- kmeans(xVar, centers=9)
gVarClust <- kmeans(gVar, centers=9)
#Filter the X features
clusterMeansX <- xVarClust$centers
smallestX <- which(clusterMeansX == min(clusterMeansX))
removeX_var <- rep("Keep", ncol(trans.measure))
removeX_var[xVarClust$cluster == smallestX] <- "Remove"
#Filter the G features
clusterMeansG <- gVarClust$centers
smallestG <- which(clusterMeansG == min(clusterMeansG))
removeG_var <- rep("Keep", ncol(methy.measure))
removeG_var[gVarClust$cluster == smallestG] <- "Remove"
time2 <- Sys.time()
time2 - time1

#How many removed?
sum(removeX_var=="Remove")
sum(removeG_var=="Remove")

#Perform CCA from tibshirani
time1 <- Sys.time()
CCA_PMA_LowVartune <- CCA.permute(x=trans.measure[,as.numeric(which(removeX_var=="Keep"))],
                                   z=methy.measure[,as.numeric(which(removeG_var=="Keep"))],
                                   typex="standard",
                                   typez="standard",
                                   outcome="quantitative",
                                   y=y)
CCA_PMALowVar <- PMA::CCA(x=trans.measure[,as.numeric(which(removeX_var=="Keep"))],
                           z=methy.measure[,as.numeric(which(removeG_var=="Keep"))],
                           typex="standard",
                           typez="standard",
                           outcome="quantitative",
                           y=y,
                           penaltyx=CCA_PMA_LowVartune$bestpenaltyx,
                           penaltyz=CCA_PMA_LowVartune$bestpenaltyz)
time2 <- Sys.time()
timeLowVar_CCA_PMA <- time2 - time1
CCA_PMALowVar$runTime <- timeLowVar_CCA_PMA

#(Optional) Save outputs
saveRDS(CCA_PMALowVar, "DATA LOCATION HERE/CCA_PMALowVar")

#Read the outputs (only if it hasn't been rerun)
CCA_PMALowVar <- readRDS("DATA LOCATION HERE/chr10_PercTumorInv/CCA_PMALowVar")

#What was selected after the filter?
selectedPMAx_LowVar <- (1:ncol(trans.measure))[as.numeric(which(removeX_var=="Keep"))][which(CCA_PMALowVar$u!=0)]
selectedPMAg_LowVar <- (1:ncol(methy.measure))[as.numeric(which(removeG_var=="Keep"))][which(CCA_PMALowVar$v!=0)]
#How many were filtered that were previously selected?
sum(selectedPMAx_before %in% as.numeric(which(removeX_var!="Keep")))
sum(selectedPMAg_before %in% as.numeric(which(removeG_var!="Keep")))
#How many genes maintained selection after filtering?
sum(selectedPMAx_LowVar %in% selectedPMAx_before)
sum(selectedPMAg_LowVar %in% selectedPMAg_before)













#-----------------------------------Compile everything into one table on selections-------------------------#

#---X: gene expression-----#
#Collect all features that were selected
selectedX_comb <- sort(unique(c(selectedPMAx_before, selectedPMAx_SuMOFIL, selectedPMAx_LowMean, selectedPMAx_LowVar)))
#Fetch the gene names
X_comb <- data.frame(Idx=selectedX_comb, 
                     geneID=colnames(trans.measure)[selectedX_comb], 
                     geneName=trans.geneInfo$geneName[match(x=colnames(trans.measure)[selectedX_comb], table=trans.geneInfo$geneID)],
                     selected_unfiltered=rep(0, length(selectedX_comb)),
                     selected_SuMOFil=rep(0, length(selectedX_comb)),
                     selected_lowmeans=rep(0, length(selectedX_comb)),
                     selected_lowvar=rep(0, length(selectedX_comb)),
                     filtered_SuMOFil=rep(0, length(selectedX_comb)),
                     filtered_lowmeans=rep(0, length(selectedX_comb)),
                     filtered_lowvar=rep(0, length(selectedX_comb)),
                     stringsAsFactors=F)
#Mark if its been selected by the different methods
X_comb[selectedX_comb %in% selectedPMAx_before, "selected_unfiltered"] <- 1
X_comb[selectedX_comb %in% selectedPMAx_SuMOFIL, "selected_SuMOFil"] <- 1
X_comb[selectedX_comb %in% selectedPMAx_LowMean, "selected_lowmeans"] <- 1
X_comb[selectedX_comb %in% selectedPMAx_LowVar, "selected_lowvar"] <- 1
#Mark if it was filtered by the different methods
X_comb[selectedX_comb %in% as.numeric(filterResults$removeX_both), "filtered_SuMOFil"] <- 1
X_comb[selectedX_comb %in% which(removeX_mean=="Remove"), "filtered_lowmeans"] <- 1
X_comb[selectedX_comb %in% which(removeX_var=="Remove"), "filtered_lowvar"] <- 1

#(Optional) Export to CSV
write.csv(X_comb, "DATA LOCATION HERE/XSelections_PercTumorInv.csv", na="")



#---G: DNA methylation-----#
#Collect all features that were selected
selectedG_comb <- sort(unique(c(selectedPMAg_before, selectedPMAg_SuMOFIL, selectedPMAg_LowMean, selectedPMAg_LowVar)))
#Fetch the gene names
G_comb <- data.frame(Idx=selectedG_comb, 
                     geneID=colnames(methy.measure)[selectedG_comb], 
                     geneSymbol=methy.geneInfo$geneSymbol[match(x=colnames(methy.measure)[selectedG_comb], table=methy.geneInfo$geneID)],
                     transcriptID=methy.geneInfo$transcriptID[match(x=colnames(methy.measure)[selectedG_comb], table=methy.geneInfo$geneID)],
                     selected_unfiltered=rep(0, length(selectedG_comb)),
                     selected_SuMOFil=rep(0, length(selectedG_comb)),
                     selected_lowmeans=rep(0, length(selectedG_comb)),
                     selected_lowvar=rep(0, length(selectedG_comb)),
                     filtered_SuMOFil=rep(0, length(selectedG_comb)),
                     filtered_lowmeans=rep(0, length(selectedG_comb)),
                     filtered_lowvar=rep(0, length(selectedG_comb)),
                     stringsAsFactors=F)
#Mark if its been selected by the different methods
G_comb[selectedG_comb %in% selectedPMAg_before, "selected_unfiltered"] <- 1
G_comb[selectedG_comb %in% selectedPMAg_SuMOFIL, "selected_SuMOFil"] <- 1
G_comb[selectedG_comb %in% selectedPMAg_LowMean, "selected_lowmeans"] <- 1
G_comb[selectedG_comb %in% selectedPMAg_LowVar, "selected_lowvar"] <- 1
#Mark if it was filtered by the different methods
G_comb[selectedG_comb %in% as.numeric(filterResults$removeG_both), "filtered_SuMOFil"] <- 1
G_comb[selectedG_comb %in% which(removeG_mean=="Remove"), "filtered_lowmeans"] <- 1
G_comb[selectedG_comb %in% which(removeG_var=="Remove"), "filtered_lowvar"] <- 1

#(Optional) Export to CSV
write.csv(G_comb, "DATA LOCATION HERE/GSelections_PercTumorInv.csv", na="")








#Hamming distance - X
library(e1071)
e1071::hamming.distance(X_comb$selected_unfiltered, X_comb$selected_SuMOFil)
e1071::hamming.distance(X_comb$selected_unfiltered, X_comb$selected_lowmeans)
e1071::hamming.distance(X_comb$selected_unfiltered, X_comb$selected_lowvar)
e1071::hamming.distance(X_comb$selected_SuMOFil, X_comb$selected_lowmeans)
e1071::hamming.distance(X_comb$selected_SuMOFil, X_comb$selected_lowvar)
e1071::hamming.distance(X_comb$selected_lowmeans, X_comb$selected_lowvar)

#Hamming distance - G
library(e1071)
e1071::hamming.distance(G_comb$selected_unfiltered, G_comb$selected_SuMOFil)
e1071::hamming.distance(G_comb$selected_unfiltered, G_comb$selected_lowmeans)
e1071::hamming.distance(G_comb$selected_unfiltered, G_comb$selected_lowvar)
e1071::hamming.distance(G_comb$selected_SuMOFil, G_comb$selected_lowmeans)
e1071::hamming.distance(G_comb$selected_SuMOFil, G_comb$selected_lowvar)
e1071::hamming.distance(G_comb$selected_lowmeans, G_comb$selected_lowvar)

