# =============================================================================
# Script 02: Pipeline
# Project:   TCGA-LIHC DEG Analysis
# Course:    INB 321G: Computational Biology
#            UT Austin, Spring 2025
# Input:     Data/CleanGroup1.csv
#            Data/CleanGroup2.csv
# Output:    Results/Fig1_PCA1.png
#            Results/Fig2_PCA2.png
#            Results/Fig3_VolcanoPlot.png
#            SupplementaryFiles/Samples.csv
#            SupplementaryFiles/Loci.csv
# =============================================================================

# ----- Load Libraries --------------------------------------------------------
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(biomaRt)
library(ggplot2)  
library(ggrepel)

# ----- Set Working Directory -------------------------------------------------
setwd("/Volumes/DONGHOA/Kim's UT 25 BIO Undergrad/4 Spring 2025/INB 321G/project/TCGA-LIHC-DEGanalysis")

# ----- Load Clean Data -------------------------------------------------------
group1 <- read.csv("Data/CleanGroup1.csv", fill=T, header = 1)
group2 <- read.csv("Data/CleanGroup2.csv", fill=T, header = 1)

# ----- Check sample sheet information ----------------------------------------
checkedColumns <- c("Project.ID","Sample.ID","Data.Type")

# Check for the expected column names
if(!all(checkedColumns%in%colnames(group1))){stop("Oh, no! You are missing important column names. Check your input file / code for group 1.")}
if(!all(checkedColumns%in%colnames(group2))){stop("Oh, no! You are missing important column names. Check your input file / code for group 2.")}

# Check for duplicate sample IDs
if(length(unique(group1$Sample.ID))<nrow(group1)){stop("Oh, no! You have some duplicate values in the Sample.ID column of group 1; investigate and correct this.")}
if(length(unique(group2$Sample.ID))<nrow(group2)){stop("Oh, no! You have some duplicate values in the Sample.ID column of group 2; investigate and correct this.")}

# Check for too many sample IDs
if((nrow(group1)+nrow(group2))>250 ){stop("You have too many samples and might run into computer issues. Strongly consider subsetting your sample sheet!!!")
}else if( (nrow(group1)+nrow(group2))>150 ){stop("You have a large number of samples. Consider subsetting your sample sheets.")}

# Check for the data.type needed for DESeq2
if(!all(c(group1$Data.Type,group2$Data.Type)=="Gene Expression Quantification")){
  stop("This pipeline is only built to use Gene Expression Quantification data")
}

# ----- Search for data to download -------------------------------------------
##### Connect the grouping information to the download tool #####

# Create a character vector that determines what projects to search through on the GDC based on the grouping files
projects <- unique(c(group1$Project.ID, group2$Project.ID))

# Do the initial search of the GDC
query1 <- GDCquery(
  project = projects,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  access = "open"
)

# Get a data.frame of the initial search results
samDf <- query1$results[[1]]

# Determine which search results correspond to the samples in the group files 
samDf$inGrp1 <- samDf$sample.submitter_id%in%group1$Sample.ID
samDf$inGrp2 <- samDf$sample.submitter_id%in%group2$Sample.ID
samDf$sampleLogic <- samDf$inGrp1 | samDf$inGrp2

# Explore your search results to make sure you have selected what you wanted to select
View(samDf[order(-samDf$sampleLogic),])

# Save the desired file barcodes (a.k.a. cases) to use to filter another query of the database
barcodes <- samDf$cases[samDf$sampleLogic]

# Search the GDC again, but this time filter the results by the file barcodes
query2 <- GDCquery(
  project = projects,
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  access = "open",
  barcode = barcodes
)

# ----- Check your query results ----------------------------------------------

# Check that the grouping files do not overlap (Should be 0)
sum(group1$Sample.ID%in%group2$Sample.ID)
sum(group2$Sample.ID%in%group1$Sample.ID)

# Check that your second search had any results at all (Should be >0)
samDf2 <- query2$results[[1]]
nrow(samDf2) 
length(group1$Sample.ID)
length(group2$Sample.ID)

# Check that your second search had results overlapping with group files (Should be >=6 and not >>>50)
sum(samDf2$sample.submitter_id%in%group1$Sample.ID)
sum(samDf2$sample.submitter_id%in%group2$Sample.ID)


# ----- Download and prepare data for analysis --------------------------------
# Download the data files (Don't download >1 Gb on Edupod server)
GDCdownload(query = query2, method = "api", files.per.chunk = 10)

# Read the downloaded files and structure the sequencing data into a single object
dds1 <- GDCprepare(query2)

# Add new columns to the sample data for use in downstream data cleaning
dds1$group1 <- dds1$sample_submitter_id%in%group1$Sample.ID
dds1$group2 <- dds1$sample_submitter_id%in%group2$Sample.ID

# Add a factor column to control how DESeq2 will group the samples
dds1$comp <- factor(x = dds1$group2,levels = c(FALSE,TRUE))

# Convert the sequencing data into DESeqDataSet class (a subclass of RangedSummarizedExperiment)
dds1 <- DESeqDataSet(dds1,design = ~comp)

# Visualize your samples' clinical data and ensure they are what you expected
df_dds1 <- as.data.frame(colData(dds1))
View(df_dds1)


# ----- Analyze your count data -----------------------------------------------
##### Detect bad samples based on PCA of normalized counts #####
# Create a function to visualize overall clustering of the samples via PCA
pcaFun <- function(ddsFun){
  # Check arguments of pcaFun prior to analysis
  if(class(ddsFun)=="DESeq2"){
    stop("ddsFun needs to be a DESeq2 class object")
  }
  if(!"comp"%in%colnames(colData(ddsFun))){
    stop("The column 'comp' needs to be the name of a factor column in ddsFun's sample descriptor section")
  }
  
  # Calculate normalized values
  ddsFun <- DESeq2::estimateSizeFactors(ddsFun)
  normCounts <- as.data.frame(counts(ddsFun, normalized = T))
  
  # Calculate the variance of counts across each gene in order to subset data
  # This is needed because only rows with > 0 variance can be used in PCA
  # The apply loops function var() (variance) across rows of normCounts
  perLocusVar   <- apply(normCounts,MARGIN = 1,var) 
  
  # Determine which genes to use based on variance
  pcaLocusLogic <- perLocusVar>median(perLocusVar[perLocusVar>0])
  if(sum(pcaLocusLogic)==0){stop("Your data did not contain loci with variable counts")}
  
  # Subsets and transposes the dataset for analysis with prcomp
  pcaInput <- t(normCounts[which(pcaLocusLogic),])
  
  # Calculate principal components with scaling (to make each locus equal in weight)
  pca <- prcomp(pcaInput,scale. = T)
  
  # Calculate the percent variance explained by each principal component
  pca_var <- pca$sdev^2
  pve <- pca_var / sum(pca_var)
  
  # Plots a scatterplot of the samples visualized via principal component analysis
  plot(
    x = pca$x[,1],y = pca$x[,2],asp = 1,
    xlab = paste0("PC1: ",round(pve[1],3)*100,"%"),
    ylab = paste0("PC2: ",round(pve[2],3)*100,"%"),
    col  = colData(ddsFun)$comp,
    pch  = as.numeric(colData(ddsFun)$comp)+1
  )
  
  # Return a list of the data
  outLs <- list(x=pca$x,pve=pve)
  return(outLs)
}
pca1Out <- pcaFun(ddsFun = dds1)


# Analyze the PCA results to determine undesirable samples / barcodes
pca1Cutoffs <- 150
abline(h = pca1Cutoffs, col="red")
dds1$pca1Logic  <- pca1Out$x[,2] > pca1Cutoffs # for pca2 (y-axis)
badSamples      <- dds1$barcode[dds1$pca1Logic]

# Subset the data to certain samples based on the SummarizedExperiment documentation
dds2 <- dds1[,!dds1$barcode%in%badSamples] 
dim(dds2) 

# After initial sample filtering, repeat PCA to see if additional filtering is needed
pca2Out <- pcaFun(dds2) 


# ----- Detect bad loci based having high frequency of 0 in either grouping ---
# Calculate proportion with zero raw counts per grouping
locusCounts <- as.data.frame(counts(dds2, normalized = F))
locusCounts_grp1At0 <- rowMeans(locusCounts[,dds2$group1]==0)
locusCounts_grp2At0 <- rowMeans(locusCounts[,dds2$group2]==0)

# Visualize the proportion of the samples with 0 read depth per grouping
hist(x = c(locusCounts_grp1At0,locusCounts_grp2At0),100)
maxPct0Cutoff  <- 0.90
abline(v = maxPct0Cutoff,col="red")

# Determine which loci meet the read depth based filter across both groups
locusCounts_logic <- locusCounts_grp1At0 < maxPct0Cutoff & locusCounts_grp2At0 < maxPct0Cutoff

# Subset the data to certain loci based on the SummarizedExperiment documentation
dds <- dds2[locusCounts_logic,]
dim(dds) # 38836 genes across 99 samples


# ----- Check subset count data -----------------------------------------------
# Calculate your total gene and sample count
nrow(dds) #The number of genes should be >>> 1,000.
ncol(dds) #The number of samples 

# Calculate the number of samples per grouping
class(dds$comp) #This needs to be a factor type object
table(dds$comp) #You should have two values, both more than 6. 

# Calculate the minimum sample raw depth (Needs to be above 0)
min(colMeans(assays(dds)[[1]]))

#  Calculate the maximum proportion of samples with 0 counts per group (Needs to be less than 100%)
max(rowMeans(assays(dds)[[1]][,dds$group1]==0))
max(rowMeans(assays(dds)[[1]][,dds$group2]==0))


# ----- Analyze the data with DESeq2 ------------------------------------------
# Repeat the conversion to DESeqDataSet class object to make sure it's still well formatted
dds <- DESeqDataSet(dds,design = ~comp)

# Recalculate  the size factors used in depth normalization after all the previous filtering
dds <- estimateSizeFactors(dds)

# Run the actual differential expression calculations
dds <- DESeq(dds)

# Organize the results of the differential expression calculations into a table
res <- results(dds) 

# Combine the results of differential expression with descriptions of the loci
resOutput <- cbind(as.data.frame(res),as.data.frame(rowRanges(dds)))

# Remove a duplicated column that causes problems downstream otherwise
resOutput <- resOutput[,!duplicated(colnames(resOutput))]

# Plot a simple scatterplot / volcano plot of the analyzed genes
plot(resOutput$log2FoldChange,-log10(resOutput$padj))

# ----- Analyze the data with DESeq2 ------------------------------------------
source("Scripts/DEG_Functions.R")

### nicePCA
out1 <- nicePCA(pca1Out, dds1$group2)
out1_line <- out1 + geom_hline(yintercept = 150, color = "red")
out1_line

ggsave(out1_line, filename = "Results/Fig1_PCA1.jpeg", device = "jpeg")

out2 <- nicePCA(pca2Out, dds2$group2)
out2

ggsave(out2, filename = "Results/Fig2_PCA2.jpeg", device = "jpeg")

### niceVolcano
volcPlot <- niceVolcano(resOutput, 0.05, 1.5)
volcPlot

ggsave(filename = "Results/Fig3_VolcanoPlot.png", volcPlot, dpi = 600, width = 84*2, height = 84, units = "mm")

### resSummary
resSummary(resOutput)


# ----- Supplementary Files ---------------------------------------------------
source("Scripts/SupplementaryDataFiles.R")

## Samples.csv
df_dds <- as.data.frame(colData(dds))
samplesFile(df_dds, "SupplementaryFiles/Samples.csv")
 
## Loci.csv
lociFile(resOutput, "SupplementaryFiles/Loci.csv")

