FROM rocker/r-ver:3.5.2

RUN apt-get update -qq && apt-get -y --no-install-recommends install \
    libxml2-dev \
    libcairo2-dev \
    libsqlite-dev \
    libssh2-1-dev \
    libcurl4-openssl-dev \
    libedit2 \
    libssl-dev \
    libmagick++-dev \
    curl \
    && rm -rf /tmp/* \
    && apt-get autoremove -y \
    && apt-get autoclean -y \
    && rm -rf /var/lib/apt/lists/*

RUN curl -Lo /usr/local/bin/install2.r https://github.com/eddelbuettel/littler/raw/master/inst/examples/install2.r \
    && chmod +x /usr/local/bin/install2.r

RUN install2.r --error \
    --ncpus -1 \
    tidyverse \
    ape \
    future \
    MonoPhy \
    && rm -rf /tmp/*

COPY scripts/* scripts/
COPY downloads/* downloads/

RUN Rscript scripts/generate_taxonomy.R family \
    && Rscript scripts/generate_taxonomy.R order \
    && Rscript scripts/generate_monophyly.R family \
    && Rscript scripts/generate_monophyly.R order \
    && Rscript scripts/generate_fossils.R \
    && rm -rf scripts/*

CMD ["R"]
