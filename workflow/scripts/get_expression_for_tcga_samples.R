#!/usr/bin/env Rscript
################################################################################
## Retrieve gene expression data from TCGA-BRCA dataset 
## Author: Lourdes Rosano
## Date created: Nov 2023
## R Version: 4.2.3
################################################################################


suppressPackageStartupMessages({
    library(optparse)
    library(TCGAbiolinks)
    library(SummarizedExperiment)
})


##################
## MAIN SCRIPT
##################

# Parse command line arguments
option_list <- list(
  make_option("--tcga_ids", default=NULL, type="character", help="String containing list of TCGA sample ids to retrieve gene expression for, separated by spaces."),
  make_option("--genes", default=NULL, type="character", help="Path to input file containing list of gene ids to retrieve gene expression for. Input file is assumed to contain one id per row (without header)."),
  make_option("--library_id", default=NULL, type="character", help="Library id corresponding to the provided gene set."),
  make_option("--output_dir", default=NULL, type="character", help="Path to the output directory where result files should be written to."),
  make_option("--download_tmp_dir", default=NULL, type="character", help="Path to a temporary directory where GDC data can be downloaded into."),
  make_option("--success_file", default=NULL, type="character", help="Path to the output success file required by the pipeline.")
)
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Split provided input string of TCGA sample ids into separate ids
tcga_sample_ids <- strsplit(opt$tcga_ids, split="\\s+")[[1]]
print(tcga_sample_ids)

# Load input file of gene ids of interest. Assumes no header and one row per id.
gene_ids <- read.table(opt$genes, header=FALSE, stringsAsFactors=FALSE, quote="", comment.char="", check.names=FALSE)

# Sanity check that no duplicated gene symbols were provided
if (sum(duplicated(gene_ids[[1]])) > 0){
  stop("Error! Input file '--genes' contains duplicated gene identifiers!")
  quit(save="no", status=1)
}


# Retrieve gene expression data from the TCGA-BRCA dataset for the provided set of TCGA samples
query <- GDCquery(
    project = "TCGA-BRCA",
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification", 
    workflow.type = "STAR - Counts",
    barcode = tcga_sample_ids
)
GDCdownload(query, method = "api", directory = opt$download_tmp_dir)
rda_file <- file.path(opt$output_dir, "TCGA-BRCA.RData")
data <- GDCprepare(query = query, directory = opt$download_tmp_dir, save = TRUE, save.filename = rda_file, remove.files.prepared = TRUE)

output_files <- c()

# Iterate individual TCGA samples of interest and further process their gene expression data
for (tcga_sample in tcga_sample_ids){
    print(tcga_sample)

    # Retrieve gene expression data for each TCGA sample of interest
    sub_data <- assays(data[,data$barcode == tcga_sample])
    # Merge all available assays for the sample into one gene expression matrix
    merged_data <- do.call(cbind, as.list(sub_data))
    # Assign assay names as columns
    colnames(merged_data) <- names(sub_data)
    # Add info about gene IDs, symbols and types
    df <- as.data.frame(merged_data)
    df$gene_id <- rownames(rowData(data))
    df$gene_name <- rowData(data)$gene_name
    df$gene_type <- rowData(data)$gene_type

    # Write complete gene expression matrix for the current sample
    full_out_file <- file.path(opt$output_dir, paste0(tcga_sample, ".gene_expression.tsv"))
    # Reorder data in preparation for writing to output
    write.table(df[,c(7:9,1:6)], file=full_out_file, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

    # Get subset of expression matrix for the provided genes of interest
    is_found <- (df$gene_name %in% gene_ids[[1]])
    subdf <- df[is_found,,drop=F]

    # Get info about the number of genes that could not be matched
    found_genes <- (gene_ids[[1]] %in% unique(subdf$gene_name))
    genes_not_found <- gene_ids[[1]][!found_genes]
    print(paste0("Input genes not found in TCGA expression matrix: ", length(genes_not_found)))

    # Write subset of expression matrix for the current sample and genes of interest 
    sub_out_file <- file.path(opt$output_dir, paste0(tcga_sample, ".", opt$library_id, ".gene_expression.mapped_genes.tsv")) 
    # Reorder data in preparation for writing to output
    write.table(subdf[,c(7:9,1:6)], file=sub_out_file, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE)

    output_files <- c(output_files, sub_out_file)
}

# Assign TCGA sample names to the already written expression files
names(output_files) <- tcga_sample_ids

# Write success file listing all the filepaths to the files of interest (i.e. subset of expression matrix for samples and genes of interest)
write.table(as.data.frame(output_files), file=opt$success_file, quote=FALSE, sep="\t", row.names=TRUE, col.names=FALSE)
