# sgRNA_library_characterization

A Snakemake workflow dedicated to the characterization of a sgRNA library.

## Workflow description

This pipeline performs the mapping of sgRNA sequences to the GRCh38 reference genome and evaluates gene expression based on the mapping results. The pipeline is implemented using Snakemake and consists of several steps:
* Input: Multi FASTA file
* Build genome index with BWA
* Map sgRNA library to genome (BWA backtrack)
* Sort aligned data by coordinate
* Report alignment stats with Samtools flagstat
* Extraction of alignment info
* Annotation of alignment and comparison of gene names
* Retrieval of TCGA-BRCA gene expression

![alt text](https://github.com/lourdesrosano/sgRNA_library_characterization/blob/main/images/rulegraph.png?raw=true)


## Installations

### Conda/Mamba

Installation instructions can be found in the [Conda documentation](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html).

### Snakemake

Installation instructions can be found in the [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).


## Before running the pipeline

- **Supported OS:**
The pipeline is intended to be run on a Linux or MacOS X system.

- **Internet connection:**
Some steps of this workflow perform online queries. Please make sure that this is possible on your computing system.


### Required input files

- **Input FASTA:**
The pipeline expects one Multi FASTA file per library, and assumes the input FASTA file contains single-end sequences and that it has already been trimmed to remove any adapters. Note that the filenames must be identical to their corresponding library IDs specified in the library map (see below).

- **Config file:**
  - Input directory  
    Before running the pipeline, file [config/config.yaml](config/config.yaml) needs to be adapted to contain the **full paths to the inputs and outputs** for the intended analysis (provided in the
    first section `input_output`).
  - Resource information  
    In addition to the input and output paths, further resource information must be provided in the section `data_resources`, namely, the reference genome file and the gencode gene annotation file.

- **Library map:**
File must be provided in the config file (`library_map`), and corresponds to a list of the library ids to be analyzed (one row per library id). The library map must contain a column with the header `library_id`. This ID will be used to name files and identify the library throughout the pipeline. See file [config/library_map.tsv](config/library_map.tsv) for a template of the required format for the library map.


- **TCGA sample map:**
File must be provided in the config file (`tcga_samples`), and corresponds to a list of the TCGA sample ids to be analyzed (one row per sample id). The sample map must contain a column with the header `sample`. These IDs will be used to retrieve gene expression data from the TCGA-BRCA dataset. See file [config/tcga_samples.tsv](config/tcga_samples.tsv) for a template of the required format for the TCGA sample map.



### Required resources

These correspond to files that are required for running the pipeline, which are listed in the `data_resources` section of the config.

- **Reference assembly genome GRCh38:**
This file corresponds to `reference_genome` in the config and is **NOT provided** in the pipeline due to its size.
Download the file from the [NCBI FTP Server](https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/) and unzip:
```
> wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
> gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz 
```

- **Gencode annotation file:**
This file corresponds to `gencode_bed` in the config and is directly provided under `resources` ([here](resources/gencode.v44.primary_assembly.annotation.genes.ucsc.bed)).
The available file was generated as follows:
```
# Download file from Gencode server
> wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz
# Unzip file
> gunzip gencode.v44.primary_assembly.annotation.gtf.gz

# Generate BED by extracting and parsing gene information
> awk -F'\t' '{if($3 == "gene"){for(i=1;i<=7;++i){printf $i"\t"}; split($9, myName, "gene_name \""); split(myName[2], myName2, "\""); print myName2[1]}}' gencode.v44.primary_assembly.annotation.gtf | cut -f1,4,5,7,8 > gencode.v44.primary_assembly.annotation.genes.bed 

# Use tool recontig to ensure that contig names are correct and match the reference genome
> recontig convert -b GRCh38 -c gencode2UCSC -f bed gencode.v44.primary_assembly.annotation.genes.bed | egrep -v '^#' | sort -k1,1V -k2,2n -k3,3n > gencode.v44.primary_assembly.annotation.genes.ucsc.bed  
```


### Running sgRNA_library_characterization

- **Software dependencies with Conda:**
This pipeline follows the best practices of the Snakemake workflow manager in providing the software needed to run the pipeline in per-rule conda environments. Those environmnents are specified under directory [envs/](workflow/envs/) in yaml files that are named `{rule_name}.yaml`. The easiest way to install and use the software is by running Snakemake with parameter `--use-conda`. Snakemake will try to find the environments of the yaml files the rules point to, and install them if they are not already available. The directory for installing the conda environments can be specified with the `--conda-prefix` parameter.

- **Software dependencies with Docker/Singularity:**
As an alternative to using Conda, it is also possible to define a (docker) container and execute Snakemake with parameter `--use-singularity`, which will execute the job within a container that is spawned from a given image. Instructions for automatically building an image with Docker are provided in a [Dockerfile](Dockerfile). For more detailed information, please refer to the corresponding [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/snakefiles/deployment.html).

- **Set pipeline location:**
Parameter `--directory` must point to the `workflow` repository of the pipeline (i.e. corresponding location to the copy of the Git repository in the user's machine). This is necessary in order to ensure that the relative paths to the pipeline scripts can be correctly resolved and interpreted by Snakemake.

Example call:

```
snakemake -s workflow/Snakefile.smk  --configfile config/config.yaml  -p -k --use-conda --conda-prefix /path/to/conda_envs/ --cores 1 --directory /path/to/sgRNA_library_characterization/workflow/   
```
Note that if the pipeline is run on a compute cluster with a job scheduling system (e.g. LSF), the commands need to be adjusted accordingly.


### Relevant output

* Directory `bwa`: one `.bam.flagstat` file per library analyzed, listing the number of sgRNA sequence alignments per FLAG type (e.g. mapped, duplicates, primary, secondary, etc.).
* Directory `gene_mapping`: one `.mapping_info.tsv` file per library analyzed, listing the chromosome name, start position, end position and strand of each mapped sgRNA sequence.
* Directory `gene_expression`: one `.gene_expression.mapped_genes.tsv` file per library analyzed and TCGA sample of interest, containing the matrix of gene expression for the given TCGA sample and the gene IDs annotated based on the mapped sgRNA sequences.


### Future development

This pipeline should be extended further in order to include the following steps:

- **Standardization of FASTA gene symbols to avoid mismatches:**

Upon inspection of the annotated genes that could not be matched to those from the FASTA descriptions of the library, it became apparent that many genes in the FASTA descriptions use deprecated symbols. By leveraging file [hgnc_complete_set.txt](https://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/hgnc_complete_set.txt) from the [HGNC archive](https://www.genenames.org/download/archive/), it would be possible to convert the gene symbols from the FASTA descriptions of the library to their most up-to-date version, and in this manner avoid multiple mismatches.


### Adapting/Integrating rules in Snakemake

Snakemake is a Python-based workflow management system for building and executing pipelines. A pipeline is made up of [rules](workflow/rules/sgrna_library_characterization.smk) that represent single steps of the analysis. In a [yaml config file](config/config.yaml) parameters and rule-specific input can be adjusted to a new analysis without changing the rules. In a [master snake file](workflow/Snakefile.smk) the desired end points of the analysis are specified, which in turn defines the steps that should be included and executed in the given pipeline run (i.e. with the input and the desired output defined, Snakemake is able infer all steps that have to be performed in-between). It is important to make sure that the format of the input and output of each rule is compatible with the previous and the subsequent rule. For more detailed information, please have a look at the excellent [Snakemake documentation](https://snakemake.readthedocs.io/en/stable/index.html).
