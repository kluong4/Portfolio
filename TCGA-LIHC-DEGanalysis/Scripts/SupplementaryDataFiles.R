# ------------------------------------------------------------------------------
# Supplementary data code:
#   - Samples.csv: 
#       - The sample barcodes used in your final analysis and how they were 
#         grouped is the minimum information needed in these files.
#   - Loci.csv: 
#       - The ensemble id of the loci used in your final analysis, their gene 
#         name, their log2 fold-change, their p-value, and their FDR-adjusted 
#         p-value is the minimum information needed.
# ------------------------------------------------------------------------------

# Function that creates Samples.csv
samplesFile <- function(data, file_name) {
  
  # extract the desired columns individually
  barcodes         <- data$barcode
  pathologic_stage <- data$ajcc_pathologic_stage
  group <- ifelse(grepl("^Stage II\\b", pathologic_stage), "Group 1", NA)
  group[is.na(group)] <- "Group 2"
  
  sample_type      <- data$sample_type #primary tumor
  tissue_organ <- data$tissue_or_organ_of_origin #liver
  diagnosis <- data$primary_diagnosis #hepta
  
  # Construct desired data frame
  samples <- data.frame(barcodes, pathologic_stage, group, sample_type, tissue_organ, diagnosis)
  
  # Check if order is all the same
  if(sum(samples$barcodes == data$barcode) != 99) {
    stop("order of sample barcodes is not the same")
  }
  
  if(sum(samples$pathologic_stage == data$ajcc_pathologic_stage) != 99) {
    stop("order of the pathologic stage is not the same")
  }
  
  if(sum(samples$sample_type == data$sample_type) != 99) {
    stop("order of the sample type is not the same")
  }
  
  write.csv(samples, file_name)
}

# Function that creates Loci.csv
lociFile <- function(data, file_name) {
  
  # Subset the desired columns
  loci <- data[c("gene_id", "gene_name", "log2FoldChange", "pvalue", "padj")]
  
  # Change column name from gene_id to ensemble_id
  colnames(loci)[1] <- "ensemble_id"
  
  write.csv(loci, file_name, row.names = FALSE)
}