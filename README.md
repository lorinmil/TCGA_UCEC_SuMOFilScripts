# Introduction

This respository contains R code to replicate the real data application found in https://doi.org/10.1101/2020.03.12.985077. These instructions will query the UCEC project data from TCGA; pull out the DNA methylation, gene expression, and clinical datasets; limit data to chromosome 10; clean/scale the data; apply SuMO-Fil, low mean, and low variance filtering; execute supervised SCCA on the unfiltered and filtered data by each of the different techniques.

Note that the results from the filtering and network detections may also be viewed directly from the FinalSelectionOutputs.xlsx file.

# Instructions

Required R packages: TCGAbiolinks, SummarizedExperiment, GenomicRanges, stringr, parallel

## STEP 1: Query TCGA data from database

To make data set RDT1 from extracting data from TCGA:

1, create a folder called Folder2 and paths "Folder2/examples/TCGA/" and "Folder2/results/examples/Rdata/UCEC/".
2, Download all the files to Folder2.
3, In Rstudio, set the working directory to Folder2, and run the following commands:

   source("./STEP1_tcga_ucec prepare clinical.R")    
   source("./STEP1_tcga_ucec prepare dna methylation.R")    
   source("./STEP1_tcga_ucec prepare transcripts.R")
    
   Note that it can take quite a long time from hours to days to download the above database depending on your connection speed.
   The raw database is stored under the path Folder2/examples/TCGA, and it is **64.4 GB**   
   The "GDCdownload" command is just needed once. Once the database is downloaded, that command can be commented out.   
   The usable .Rdata file is stored under the path Folder2/results/examples/Rdata/UCEC/   
   clinical_UCEC.Rdata is 135 KB   
   dnamethy_UCEC.Rdata is **4.29 GB**   
   transcripts_UCEC.Rdata is 392 MB   
   So, arrange sufficient storage before downloading.
   
   ## STEP 2: Clean up Initial TCGA data
   
   Run the code from the STEP2_CleanTCGAData.R file. Make sure to update folder locations.
   
   ## STEP3 3: Conduct Analysis for Percent Tumor Invasion on Chromosome 10  
   
   Run the code from the STEP3_NetworkAnalysis_Chr10_PercTumorInv.R file. Make sure to update folder locations.
   
   ## STEP3 4: Conduct Analysis for Total Pelvic Lymph Node Ratio on Chromosome 10
   
   Run the code from the STEP4_NetworkAnalysis_Chr10_PelvLNR.R file. Make sure to update folder locations.
   
