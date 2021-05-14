library(TCGAbiolinks)
library(SummarizedExperiment)
library(GenomicRanges)
library(stringr)
library(parallel)

#Load in SuMOFil library from GitHub
library(SuMOFil)

#------------------------------Set up data-------------------------------#
#Load data (update path to where TCGA Rdata files are located from STEP 1)
setwd('./Folder2')

#Read the TCGA files that were constructed earlier
load('./results/examples/Rdata/UCEC/clinical_UCEC.rdata')
load('./results/examples/Rdata/UCEC/dnamethy_UCEC.rdata')
load('./results/examples/Rdata/UCEC/transcripts_UCEC.rdata')

#Pull out gene info
trans.rowRanges <- rowRanges(trans.UCEC.sumy)
methy.rowRanges <- rowRanges(methy.UCEC.sumy)
#Organize gene info into a matrix
trans.geneInfo <- data.frame(geneID = trans.rowRanges$ensembl_gene_id,
                             geneName = trans.rowRanges$external_gene_name,
                             origGeneID = trans.rowRanges$original_ensembl_gene_id)
methy.geneInfo <- data.frame(geneID = methy.rowRanges$Composite.Element.REF,
                             geneSymbol = methy.rowRanges$Gene_Symbol,
                             transcriptID = methy.rowRanges$Transcript_ID)

#Extract the chromosomal information
methyChroms <- as.character(seqnames(methy.compound.info))
transChroms <- as.character(seqnames(trans.UCEC.sumy))


#Remove the unnecessary stuff from the R environment to save space
rm(list=base::setdiff(ls(), c("trans.measure", "methy.measure", "ucec.clinical", "methyChroms", "transChroms", "",
                              "trans.geneInfo", "methy.geneInfo")))

#Extract the common patient IDs
colnames(trans.measure) <- substr(colnames(trans.measure), 9, 12)
colnames(methy.measure) <- substr(colnames(methy.measure), 9, 12)
keepIDs <- colnames(methy.measure)[colnames(methy.measure) %in% colnames(trans.measure)]
rownames(ucec.clinical) <- ucec.clinical$patient_id
ucec.clinical <- ucec.clinical[keepIDs,]



#---Save up to this point---#
#ONLY RUN ONCE! (update location where to put cleaned up TCGA data)
save.image("OUTPUT DATA LOCATION HERE/UCEC_TCGA.Rdata")


