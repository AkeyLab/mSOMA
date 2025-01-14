## Installing mSOMA
mSOMA is a collection of multiple tools, scripts, and programming languages
that have been packaged together with conda, follow the installation instructions
in the main README.md file to install the package.

## Example data

Included in this `examples/` directory is example data created by
using the [ART](https://www.niehs.nih.gov/research/resources/software/biostatistics/art)
tool to simulate RNA-seq data from chrM on GRCh38 in the `synthetic_data/` folder.

The following synthetic data files are included for usage with the `msoma` cli tool:
* chrM_seed1_ART.bam
* chrM_seed1_ART.bam.bai
* chrM.bed
* chrM.fa

## Running mSOMA on the example data after installation

After you have `msoma` installed, you can run the following command to
run both the `msoma count` and `msoma mle` commands on the example data:

```bash
sh run_msoma.sh
```

You can first run `msoma check-dependencies` to ensure that all dependencies are installed before running the `run_msoma.sh` script.

Note that running `msoma` on the syntethic data will result in some NaN P_VALs. This is due to the small size of the data and the fact that the data is synthetic. Running `msoma` on real data should not result in NaN P_VALs.

## Running mSOMA with nextflow with docker

As part of the GitHub Actions CI/CD pipeline, we have included a `nextflow`
script that runs `msoma` on the example data using the `msoma` docker image that
is created earlier in the msoma pipeline.

Docker is used here instead of conda to ensure the GHA pipeline runs
on the current version of the software and all dependencies are installed.
Installing from bioconda would not reflect the current version of the software
under development.

The following files are part of the `nextflow` pipeline:
* `mSOMA.nf`
* `sample_sheet_chrM.tsv`

The `mSOMA.nf` nextflow script can also be used as a starting point
for running `msoma` on your own data, but I'd recommend using nextflow
with conda environments for a more reproducible pipeline.