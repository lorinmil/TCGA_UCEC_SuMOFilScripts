library(TCGAbiolinks)
library(SummarizedExperiment)



#### clincal ###
query.ucec.clinical <- GDCquery(project = "TCGA-UCEC",data.category = "Clinical",legacy = FALSE)

GDCdownload(query = query.ucec.clinical, method='client', directory = "./examples/TCGA")

#NOTE: HAVE TO EXCLUDE PATIENTS WITH TXT FILES!!! (otherwise it errors out)
fileNames <- query.ucec.clinical$results[[1]]$file_name
query.ucec.clinical$results[[1]] <- query.ucec.clinical$results[[1]][substr(fileNames, nchar(fileNames)-3, nchar(fileNames))!=".txt",]

ucec.clinical <- GDCprepare_clinic(query = query.ucec.clinical,
																	clinical.info = "patient",
																	directory = "./examples/TCGA")

#Save the data
save.image(file = "./examples/Rdata/UCEC/clinical_UCEC.Rdata")
saveRDS(ucec.clinical, "./examples/Rdata/UCEC/clinical.rda")







### transcripts ###
query.ucec.transcripts <- GDCquery(project = "TCGA-UCEC",
                                   data.category = "Transcriptome Profiling",
                                   data.type = "Gene Expression Quantification",
                                   legacy = FALSE,
                                   workflow.type = "HTSeq - FPKM-UQ")

GDCdownload(query = query.ucec.transcripts, directory = "./examples/TCGA")


trans.UCEC <- GDCprepare(query = query.ucec.transcripts,
                         directory = "./examples/TCGA",
                         summarizedExperiment = FALSE)


trans.UCEC.sumy <- GDCprepare(query = query.ucec.transcripts,
                              directory = "./examples/TCGA",
                              summarizedExperiment = TRUE)

trans.measure <- assay(trans.UCEC.sumy)

trans.compound.info <- rowRanges(trans.UCEC.sumy)

trans.patient.info <- colData(trans.UCEC.sumy)


save.image(file = "./examples/Rdata/UCEC/transcripts_UCEC.Rdata")
