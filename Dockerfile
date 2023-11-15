FROM condaforge/mambaforge:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="42309fa21793f4b0f83e1e85ecf1fa4e6c8bcd938599c9178704e19748b9ccc8"

# Step 1: Retrieve conda environments

# Conda environment:
#   source: envs/bedtools_v2.30.0.yaml
#   prefix: /conda-envs/bb168372ce6b183032329fef3ce19beb
#   name: bedtools_v2.30.0
#   channels:
#     - bioconda
#   dependencies:
#     - bedtools=2.30.0
RUN mkdir -p /conda-envs/bb168372ce6b183032329fef3ce19beb
COPY envs/bedtools_v2.30.0.yaml /conda-envs/bb168372ce6b183032329fef3ce19beb/environment.yaml

# Conda environment:
#   source: envs/bwa_v0.7.17.yaml
#   prefix: /conda-envs/e2efcf0503ddb0ec3b71aa4735cf3244
#   name: bwa_v0.7.17
#   channels:
#     - bioconda
#   dependencies:
#     - bwa=0.7.17
RUN mkdir -p /conda-envs/e2efcf0503ddb0ec3b71aa4735cf3244
COPY envs/bwa_v0.7.17.yaml /conda-envs/e2efcf0503ddb0ec3b71aa4735cf3244/environment.yaml

# Conda environment:
#   source: envs/get_expression_for_tcga_samples.yaml
#   prefix: /conda-envs/a45d23bea30009bb1a45ca05aa0266e0
#   name: get_expression_for_tcga_samples
#   channels:
#     - bioconda
#     - conda-forge
#     - r
#   dependencies:
#     - r-base
#     - r-optparse
#     - bioconductor-tcgabiolinks
RUN mkdir -p /conda-envs/a45d23bea30009bb1a45ca05aa0266e0
COPY envs/get_expression_for_tcga_samples.yaml /conda-envs/a45d23bea30009bb1a45ca05aa0266e0/environment.yaml

# Conda environment:
#   source: envs/reformat_and_compare_gene_annotations.yaml
#   prefix: /conda-envs/bd5a8f156ced9db075dc4ef6abb08b2d
#   name: reformat_and_compare_gene_annotations
#   channels:
#     - anaconda
#   dependencies:
#     - python=3.12.0
RUN mkdir -p /conda-envs/bd5a8f156ced9db075dc4ef6abb08b2d
COPY envs/reformat_and_compare_gene_annotations.yaml /conda-envs/bd5a8f156ced9db075dc4ef6abb08b2d/environment.yaml

# Conda environment:
#   source: envs/samtools_v1.18.yaml
#   prefix: /conda-envs/dea57e10724c95e8a78cd9930e5b5f34
#   name: samtools_v1.18
#   channels:
#     - bioconda
#     - conda-forge
#   dependencies:
#     - samtools=1.18
RUN mkdir -p /conda-envs/dea57e10724c95e8a78cd9930e5b5f34
COPY envs/samtools_v1.18.yaml /conda-envs/dea57e10724c95e8a78cd9930e5b5f34/environment.yaml

# Step 2: Generate conda environments

RUN mamba env create --prefix /conda-envs/bb168372ce6b183032329fef3ce19beb --file /conda-envs/bb168372ce6b183032329fef3ce19beb/environment.yaml && \
    mamba env create --prefix /conda-envs/e2efcf0503ddb0ec3b71aa4735cf3244 --file /conda-envs/e2efcf0503ddb0ec3b71aa4735cf3244/environment.yaml && \
    mamba env create --prefix /conda-envs/a45d23bea30009bb1a45ca05aa0266e0 --file /conda-envs/a45d23bea30009bb1a45ca05aa0266e0/environment.yaml && \
    mamba env create --prefix /conda-envs/bd5a8f156ced9db075dc4ef6abb08b2d --file /conda-envs/bd5a8f156ced9db075dc4ef6abb08b2d/environment.yaml && \
    mamba env create --prefix /conda-envs/dea57e10724c95e8a78cd9930e5b5f34 --file /conda-envs/dea57e10724c95e8a78cd9930e5b5f34/environment.yaml && \
    mamba clean --all -y
