
mSOMA is a collection of multiple tools, scripts, and programming languages
that have been packaged together to provide a docker image
that should make it easier to get this tool up and running on your system.

## Running the Docker image from Docker Hub

To run the docker image from Docker Hub, you can run the following command:

```docker run -rm rbiermanpu/somatic-mutations:dev```

## Building and running the Docker image locally

If you want to build the docker image locally, you can do so by running
the following command from the root of the repository where the
Dockerfile is located. Note that you need to have docker installed and
the build process will take a while to complete, 20+ minutes, and the image takes approximately 2 GB of space.

```docker build -t local_betabin --target prod .```

Then you can run the docker image by running the following command:
```docker run --rm local_betabin betabin```

which should output the help message of the betabin tool.

You can run the example mitochondrial data by running the following command from within the examples directory (`cd examples`):
```docker run --rm -v $PWD:/betabinomial local_betabin count --fasta synthetic_data/chrM.fa --bed synthetic_data/chrM.bed --seq-length 150 -o test_out.counts.gz synthetic_data/chrM_seed1_ART.bam```

This will run the betabin tool on the example data and output
the results to the current directory with prefix `test_out`.
