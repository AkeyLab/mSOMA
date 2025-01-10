# Description: Example of running the betabin nextflow pipeline on SLURM
# uses singularity instead of docker, still requires nextflow to be installed
# memory/time/cpu requirements are specified in the SLURM.config
nextflow run \
    -config SLURM.config \
    -with-singularity docker://rbiermanpu/somatic-mutations:dev \
    --sample-sheet sample_sheet_chrM.tsv \
    --output-dir nextflow_output \
    --seq-length 150 \
    mSOMA.nf
