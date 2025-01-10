# Build based on bioconductor image since we need to install Bioconductor R packages
FROM bioconductor/bioconductor_docker:RELEASE_3_17 as base

# Change to root for installations
USER root

# Install system dependencies (bioconductor_docker is built on debian so use apt-get)
RUN apt-get update -y && apt-get install -y \
    wget \
    make \
    g++ \
    python3-pip \
    lbzip2 \
    libz-dev

# Install dependencies for betabinomial MLE fitting R script
# Install Bioconductor-based packages
RUN R -e 'BiocManager::install(ask = F)' && R -e 'BiocManager::install(c("survcomp","Biostrings","qvalue", ask=F))'

# For packages that are not available through bioconductor, install using install2.r from littler
RUN apt-get install -y -qq littler \
    && install2.r --error --deps TRUE VGAM \
    && install2.r --error --deps TRUE argparse \
    && install2.r --error --deps TRUE tidyverse \
    && install2.r --error --deps TRUE data.table \
    && install2.r --error --deps TRUE bbmle \
    && install2.r --error --deps TRUE dplyr

# Install samtools
RUN wget -O samtools.tar.bz2 https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2 \
    && tar -xf samtools.tar.bz2 \
    && cd samtools-1.18 \
    && ./configure --without-curses --disable-bz2 --disable-lzma \
    && make \
    && make install \
    && cd ..

# Install bamUtil (requires libStatGen)
RUN wget -O bamUtil.tar.gz https://github.com/statgen/bamUtil/archive/refs/tags/v1.0.15.tar.gz \
     && wget -O libStatGen.tar.gz https://github.com/statgen/libStatGen/archive/refs/tags/v1.0.15.tar.gz \
     && tar -xzf bamUtil.tar.gz \
     && tar -xzf libStatGen.tar.gz \
     && mv libStatGen-1.0.15 libStatGen \
     && cd bamUtil-1.0.15 \
     && make \
     && make install

RUN pip3 install --upgrade pip

# Install MtSOMA python package, copy package files to container
RUN mkdir /msoma
WORKDIR /msoma
ADD pyproject.toml .
ADD src src

# Test build stage
FROM base as test
ADD tests tests
ADD requirements_test.txt .
RUN pip3 install -r requirements_test.txt
RUN pip3 install -e .
ENTRYPOINT ["pytest", "--docker", "--cov=src/", "--cov-report=xml", "--cov-report=html", "tests"]

# Production build stage to run MtSOMA
# Switch back to non-root user for execution
FROM base as prod
RUN pip3 install .
USER 1001
