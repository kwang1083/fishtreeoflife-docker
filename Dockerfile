FROM rocker/r-ver:3.5.2

RUN apt-get update -qq && apt-get -y --no-install-recommends install \
  libxml2-dev \
  libcairo2-dev \
  libsqlite-dev \
  libmariadbd-dev \
  libmariadb-client-lgpl-dev \
  libpq-dev \
  libssh2-1-dev \
  libcurl4-openssl-dev \
  libedit2 \
  libssl-dev \
  libmagick++-dev \
  && install2.r --error \
    --deps TRUE \
    tidyverse \
    dplyr \
    devtools \
    formatR \
    remotes \
    selectr \
    ape \
    future \
    MonoPhy

COPY downloads/* downloads/

COPY scripts/* scripts/

RUN Rscript scripts/generate_taxonomy.R family \
    && Rscript scripts/generate_taxonomy.R order \
    && Rscript scripts/generate_monophyly.R family \
    && Rscript scripts/generate_monophyly.R order \
    && Rscript scripts/generate_fossils.R \
    && rm -rf scripts/*

CMD ["R"]
