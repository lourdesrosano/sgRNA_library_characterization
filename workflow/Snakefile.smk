# Master Snakefile for sgRNA library characterization

from snakemake.utils import min_version
import os
import glob

# Define minimum Snakemake version
min_version("7.32.4")

# Include Config file
configfile: "../config/config.yaml"

# Inputs and outputs
IN_DIR = config["input_output"]["input_dir"]
OUT_DIR = config["input_output"]["output_dir"]
LIBRARY_MAP = config["input_output"]["library_map"]
TMP_DIR = config["input_output"]["tmp_dir"]

# Include common helper functions used in the pipeline
include: "rules/misc_snake.smk"
# Include pipeline rules
include: "rules/sgrna_library_characterization.smk"

rule sgRNA_library_characterization:
    input:
        expand(BWA_OUT + "{library}.bam", library = library_ids),
        expand(BWA_OUT + "{library}.bam.flagstat", library = library_ids),
        expand(GENE_MAP_OUT + "{library}.mapping_info.tsv", library = library_ids),
        expand(GENE_MAP_OUT + "{library}.gene_comparison.tsv", library = library_ids),
        expand(GENE_EXPR_OUT + "{library}.get_expression_for_tcga_samples.success.txt", library = library_ids)
    output:
        OUT_DIR + "complete.txt"
    shell:
        "date > {output}"
