# Rules for sgRNA library characterization

containerized: "docker://lourdesrosano/sgrna_library_characterization:latest"


# This rule indexes a reference genome using BWA
rule bwa_index:
    input:
        reference = config["data_resources"]["reference_genome"]
    output:
        reference_amb = config["data_resources"]["reference_genome"] + ".amb",
        reference_ann = config["data_resources"]["reference_genome"] + ".ann",
        reference_bwt = config["data_resources"]["reference_genome"] + ".bwt",
        reference_pac = config["data_resources"]["reference_genome"] + ".pac",
        reference_sa = config["data_resources"]["reference_genome"] + ".sa"
    conda:
        "../envs/bwa_v0.7.17.yaml"
    log:
        config["data_resources"]["reference_genome"] + "bwa_index.log"
    shell:
        "bwa index -a bwtsw {input.reference} "
        "2>{log}"


if not "BWA_IN" in globals():
    BWA_IN = IN_DIR
if not "BWA_OUT" in globals():
    BWA_OUT = OUT_DIR + "bwa/"

# This rule aligns short single-end reads using BWA-backtrack
rule bwa_backtrack:
    input:
        fasta = BWA_IN + "{library}.fa",
        reference = config["data_resources"]["reference_genome"],
        reference_amb = config["data_resources"]["reference_genome"] + ".amb",
        reference_ann = config["data_resources"]["reference_genome"] + ".ann",
        reference_bwt = config["data_resources"]["reference_genome"] + ".bwt",
        reference_pac = config["data_resources"]["reference_genome"] + ".pac",
        reference_sa = config["data_resources"]["reference_genome"] + ".sa"
    output:
        sam = BWA_OUT + "{library}.sam"
    conda:
        "../envs/bwa_v0.7.17.yaml"
    params:
        sai = BWA_OUT + "{library}.sai"
    log:
        BWA_OUT + "{library}.bwa_backtrack.log"
    shell:
        "bwa aln {input.reference} {input.fasta} > {params.sai} && "
        "bwa samse {input.reference} {params.sai} {input.fasta} > {output.sam} "
        "2>{log}"


# This rule sorts BAM by coordinate using samtools sort
rule samtools_sort:
    input:
        sam = "{library}.sam"
    output:
        bam = "{library}.bam"
    conda:
        "../envs/samtools_v1.18.yaml"
    log:
        "{library}.samtools_sort.log"
    shell:
        "samtools sort {input.sam} > {output.bam} "
        "2>{log}"


# This rule extracts and reports BAM statistics using samtools flagstat
rule samtools_flagstat:
    input:
        bam = "{library}.bam"
    output:
        flagstat = "{library}.bam.flagstat"
    conda:
        "../envs/samtools_v1.18.yaml"
    log:
        "{library}.samtools_flagstat.log"
    shell:
        "samtools flagstat {input.bam} > {output.flagstat} "
        "2>{log}"


if not "GENE_MAP_IN" in globals():
    GENE_MAP_IN = BWA_OUT
if not "GENE_MAP_OUT" in globals():
    GENE_MAP_OUT = OUT_DIR + "gene_mapping/"

# This rule converts a BAM file into a BED file using bedtools bamtobed
rule bedtools_bamtobed:
    input:
        bam = GENE_MAP_IN + "{library}.bam"
    output:
        bed = GENE_MAP_OUT + "{library}.bed"
    conda:
        "../envs/bedtools_v2.30.0.yaml"
    log:
        GENE_MAP_OUT + "{library}.bedtools_bamtobed.log"
    shell:
        "bedtools bamtobed -i {input.bam} > {output.bed} "
        "2>{log}"


# This rule extracts relevant mapping info from the alignment file in BED format
rule get_mapping_info:
    input:
        bed = "{library}.bed"
    output:
        outfile = "{library}.mapping_info.tsv"
    log:
        "{library}.get_mapping_info.log"
    shell:
        "echo -e \"Mapped_sgRNA_sequence\tChrom\tStart\tEnd\tStrand\" > {output.outfile} && "
        "awk -v FS=\"\t\" -v OFS=\"\t\" '{{print $4,$1,$2,$3,$6}}' {input.bed} >> {output.outfile} "
        "2>{log}"


# This rule prepares the input BED file to be annotated
# (extracts unique regions and gene names from the FASTA descriptions)
rule get_bed_for_annotation:
    input:
        bed = "{library}.bed"
    output:
        bed = "{library}.uniq_gene_regions.bed"
    log:
        "{library}.get_bed_for_annotation.log"
    shell:
        "cut -f1-4 {input.bed} | awk -v FS=\"|\" -v OFS=\"\t\" '$1=$1' | "
	"cut -f1-3,6 | sort | uniq | sort -k1,1V -k2,2n -k3,3n > {output.bed} "
        "2>{log}"


# This rule annotates a BED file with the overlapping genes from Gencode using bedtools intersect
rule annotate_bed_with_genes:
    input:
        bed = "{library}.uniq_gene_regions.bed",
        annotation_bed = config["data_resources"]["gencode_bed"]
    output:
        bed = "{library}.annotated.bed"
    conda:
        "../envs/bedtools_v2.30.0.yaml"
    log:
        "{library}.annotate_bed_with_genes.log"
    shell:
	# Block below is required so that awk command can be interpreted correctly by Snakemake  
        r"""
        bedtools intersect -a {input.bed} -b {input.annotation_bed} -wa -wb > tmp1.txt && 
        bedtools intersect -a {input.bed} -b {input.annotation_bed} -wa -v > tmp2.txt && 
        cat tmp1.txt tmp2.txt | sort -k1,1V -k2,2n -k3,3n | 
        awk '{{sub(/\r/,""); if(NF > 1){{printf $0; for(i = NF; i < 9; ++i){{printf "\t*"}}; printf "\n"}}}}' > {output.bed} && 
        rm tmp1.txt tmp2.txt 
        2>{log}
        """


# This rule reformats an annotated BED file and compares the assigned gene annotations with the original ones
# (genes based on mapping vs. genes based on FASTA descriptions)
rule reformat_and_compare_gene_annotations:
    input:
        bed = "{library}.annotated.bed"
    output:
        tsv = "{library}.gene_comparison.tsv"
    log:
        "{library}.reformat_and_compare_gene_annotations.log"
    conda:
        "../envs/reformat_and_compare_gene_annotations.yaml"
    shell:
        "python scripts/reformat_and_compare_gene_annotations.py "
	"--infile {input.bed} "
	"--outfile {output.tsv} "
        "2>{log}"


# This rule extracts the set of gene annotations to be used for retrieving TCGA-BRCA expression data
rule get_library_gene_set:
    input:
        tsv = "{library}.gene_comparison.tsv"
    output:
        outfile = "{library}.gene_set.txt"
    log:
        "{library}.get_library_gene_set.log"
    shell:
        "cut -f5 {input.tsv} | tr \';\' \'\n\' | sort | uniq | grep -v \'\\*\'> {output.outfile} "
        "2>{log}"


if not "GENE_EXPR_IN" in globals():
    GENE_EXPR_IN = GENE_MAP_OUT
if not "GENE_EXPR_OUT" in globals():
    GENE_EXPR_OUT = OUT_DIR + "gene_expression/"

# This rule extracts the gene expression matrix for the provided set of TCGA-BRCA samples
# and the provided set of gene annotations from rule 'get_library_gene_set'
rule get_expression_for_tcga_samples:
    input:
        infile = GENE_EXPR_IN + "{library}.gene_set.txt"
    output:
        success = GENE_EXPR_OUT + "{library}.get_expression_for_tcga_samples.success.txt"
    params:
        tcga_sample_list = tcga_sample_ids,
        output_dir = GENE_EXPR_OUT,
        tmp_dir = TMP_DIR
    conda:
        "../envs/get_expression_for_tcga_samples.yaml"
    log:
        GENE_EXPR_OUT + "{library}.get_expression_for_tcga_samples.log"
    shell:
        "Rscript scripts/get_expression_for_tcga_samples.R "
        "--tcga_ids '{params.tcga_sample_list}' "
        "--genes {input.infile} "
        "--library_id {wildcards.library} "
        "--output_dir {params.output_dir} "
        "--download_tmp_dir {params.tmp_dir} "
        "--success_file {output.success} "
        "2>{log}"
