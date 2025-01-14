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