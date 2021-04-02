FROM rocker/tidyverse:4.0.4 AS build

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    libglpk40 \
  && rm -rf /var/lib/apt/lists/*

RUN install2.r --error \
    --ncpus -1 \
    ape \
    future \
    phangorn \
    && rm -rf /tmp/*

COPY downloads/* downloads/

COPY R/lib.R R/

COPY R/generate_taxonomy.R R/
RUN Rscript R/generate_taxonomy.R

COPY R/generate_monophyly.R R/monophy_minimal.R R/
RUN Rscript R/generate_monophyly.R family \
    && Rscript R/generate_monophyly.R order

COPY R/generate_fossils.R R/
RUN Rscript R/generate_fossils.R

FROM alpine:3.13.4 AS files

COPY --from=build assets /assets
COPY --from=build _fossils /_fossils
COPY --from=build _data /_data
COPY --from=build downloads /downloads

CMD ["sh"]
