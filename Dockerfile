FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="4cd48c6d54ef82c5a5cd171bf8a457087f244e20c0a521c057b04e4533e1ad22"
# Conda environment:
#   source: workflow/envs/amr.yaml
#   prefix: /conda-envs/2fcd4755cf11099b5eb7dd73b434f516
#   name: amr
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - blast
#     - abricate
#   prefix: /home/asus/miniforge3/envs/amr
RUN mkdir -p /conda-envs/2fcd4755cf11099b5eb7dd73b434f516
COPY workflow/envs/amr.yaml /conda-envs/2fcd4755cf11099b5eb7dd73b434f516/environment.yaml
# Conda environment:
#   source: workflow/envs/asm.yaml
#   prefix: /conda-envs/714bb2b69ba3873f4ec1fbc30af65e8c
#   name: asm
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - shovill
#     - quast
#     - flye
#   prefix: /home/asus/miniforge3/envs/asm
RUN mkdir -p /conda-envs/714bb2b69ba3873f4ec1fbc30af65e8c
COPY workflow/envs/asm.yaml /conda-envs/714bb2b69ba3873f4ec1fbc30af65e8c/environment.yaml
# Conda environment:
#   source: workflow/envs/clair3.yaml
#   prefix: /conda-envs/67e7c42b2b6333df5f806f7292b45578
#   name: clair3
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - clair3
#   prefix: /home/asus/miniforge3/envs/clair3
RUN mkdir -p /conda-envs/67e7c42b2b6333df5f806f7292b45578
COPY workflow/envs/clair3.yaml /conda-envs/67e7c42b2b6333df5f806f7292b45578/environment.yaml
# Conda environment:
#   source: workflow/envs/coverage.yaml
#   prefix: /conda-envs/aca386b0e3ef6d2120f6a04a9c353de1
#   name: coverage
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - samtools
#     - mosdepth
RUN mkdir -p /conda-envs/aca386b0e3ef6d2120f6a04a9c353de1
COPY workflow/envs/coverage.yaml /conda-envs/aca386b0e3ef6d2120f6a04a9c353de1/environment.yaml
# Conda environment:
#   source: workflow/envs/ont.yaml
#   prefix: /conda-envs/281e5f875e2437da27eb5a4ba5997474
#   name: ont
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - filtlong
#     - minimap2
#     - samtools
RUN mkdir -p /conda-envs/281e5f875e2437da27eb5a4ba5997474
COPY workflow/envs/ont.yaml /conda-envs/281e5f875e2437da27eb5a4ba5997474/environment.yaml
# Conda environment:
#   source: workflow/envs/phylo.yaml
#   prefix: /conda-envs/ce397c09e0adb4c89ffb0136af8d297e
#   name: phylo
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - snp-dists
#     - snp-sites
#     - iqtree
#     - gubbins
#   prefix: /home/asus/miniforge3/envs/phylo
RUN mkdir -p /conda-envs/ce397c09e0adb4c89ffb0136af8d297e
COPY workflow/envs/phylo.yaml /conda-envs/ce397c09e0adb4c89ffb0136af8d297e/environment.yaml
# Conda environment:
#   source: workflow/envs/plot.yaml
#   prefix: /conda-envs/d474078edf5a49889e5c939b16ca00cb
#   name: plot
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - python=3.11
#     - pandas
#     - matplotlib
#     - pyyaml
#     - biopython
RUN mkdir -p /conda-envs/d474078edf5a49889e5c939b16ca00cb
COPY workflow/envs/plot.yaml /conda-envs/d474078edf5a49889e5c939b16ca00cb/environment.yaml
# Conda environment:
#   source: workflow/envs/qc.yaml
#   prefix: /conda-envs/254629f540476d8a72418932a073a585
#   name: qc
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - fastqc
#     - fastp=1.3.7
#     - multiqc
#     - nanoplot
#     - seqkit
RUN mkdir -p /conda-envs/254629f540476d8a72418932a073a585
COPY workflow/envs/qc.yaml /conda-envs/254629f540476d8a72418932a073a585/environment.yaml
# Conda environment:
#   source: workflow/envs/ref.yaml
#   prefix: /conda-envs/18d789a616a95030492eb91335ae3a6a
#   name: ref
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - ncbi-datasets-cli
#     - unzip
RUN mkdir -p /conda-envs/18d789a616a95030492eb91335ae3a6a
COPY workflow/envs/ref.yaml /conda-envs/18d789a616a95030492eb91335ae3a6a/environment.yaml
# Conda environment:
#   source: workflow/envs/sim.yaml
#   prefix: /conda-envs/d1f3191c60d912ee3e3454836b96d07d
#   name: sim
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - wgsim
RUN mkdir -p /conda-envs/d1f3191c60d912ee3e3454836b96d07d
COPY workflow/envs/sim.yaml /conda-envs/d1f3191c60d912ee3e3454836b96d07d/environment.yaml
# Conda environment:
#   source: workflow/envs/snippy.yaml
#   prefix: /conda-envs/191b380bfc77297b7ea85c837f5c0488
#   name: snippy
#   channels:
#     - conda-forge
#     - bioconda
#   dependencies:
#     - mlst
#     - snippy
#     - snp-dists
#   prefix: /home/asus/miniforge3/envs/snippy
RUN mkdir -p /conda-envs/191b380bfc77297b7ea85c837f5c0488
COPY workflow/envs/snippy.yaml /conda-envs/191b380bfc77297b7ea85c837f5c0488/environment.yaml

RUN conda env create --prefix /conda-envs/2fcd4755cf11099b5eb7dd73b434f516 --file /conda-envs/2fcd4755cf11099b5eb7dd73b434f516/environment.yaml && \
    conda env create --prefix /conda-envs/714bb2b69ba3873f4ec1fbc30af65e8c --file /conda-envs/714bb2b69ba3873f4ec1fbc30af65e8c/environment.yaml && \
    conda env create --prefix /conda-envs/67e7c42b2b6333df5f806f7292b45578 --file /conda-envs/67e7c42b2b6333df5f806f7292b45578/environment.yaml && \
    conda env create --prefix /conda-envs/aca386b0e3ef6d2120f6a04a9c353de1 --file /conda-envs/aca386b0e3ef6d2120f6a04a9c353de1/environment.yaml && \
    conda env create --prefix /conda-envs/281e5f875e2437da27eb5a4ba5997474 --file /conda-envs/281e5f875e2437da27eb5a4ba5997474/environment.yaml && \
    conda env create --prefix /conda-envs/ce397c09e0adb4c89ffb0136af8d297e --file /conda-envs/ce397c09e0adb4c89ffb0136af8d297e/environment.yaml && \
    conda env create --prefix /conda-envs/d474078edf5a49889e5c939b16ca00cb --file /conda-envs/d474078edf5a49889e5c939b16ca00cb/environment.yaml && \
    conda env create --prefix /conda-envs/254629f540476d8a72418932a073a585 --file /conda-envs/254629f540476d8a72418932a073a585/environment.yaml && \
    conda env create --prefix /conda-envs/18d789a616a95030492eb91335ae3a6a --file /conda-envs/18d789a616a95030492eb91335ae3a6a/environment.yaml && \
    conda env create --prefix /conda-envs/d1f3191c60d912ee3e3454836b96d07d --file /conda-envs/d1f3191c60d912ee3e3454836b96d07d/environment.yaml && \
    conda env create --prefix /conda-envs/191b380bfc77297b7ea85c837f5c0488 --file /conda-envs/191b380bfc77297b7ea85c837f5c0488/environment.yaml && \
    conda clean --all -y
