FROM rocker/r-ver:3.5.1

# From rocker/tidyverse image, but without RStudio!
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
  libxml2-dev \
  libcairo2-dev \
  libsqlite3-dev \
  libssh2-1-dev \
  libcurl4-openssl-dev \
  libssl-dev \
  unixodbc-dev \
  && install2.r --error \
    --deps TRUE \
    tidyverse \
    caTools \
    BiocManager \
    MonoPhy \
    ape \
    future

COPY data/* data/

COPY scripts/* scripts/

RUN Rscript scripts/generate_taxonomy.R family \
    && Rscript scripts/generate_taxonomy.R order \
    && Rscript scripts/generate_monophyly.R family \
    && Rscript scripts/generate_monophyly.R order \
    && Rscript scripts/generate_fossils.R \
    && rm -rf data/* scripts/*

CMD ["R"]
