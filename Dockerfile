FROM rocker/tidyverse:3.6.1 AS build

RUN install2.r --error \
    --ncpus -1 \
    ape \
    future \
    phangorn \
    && rm -rf /tmp/*

COPY downloads/* downloads/

COPY scripts/lib.R scripts/

COPY scripts/generate_taxonomy.R scripts/
RUN Rscript scripts/generate_taxonomy.R

COPY scripts/generate_monophyly.R scripts/monophy_minimal.R scripts/
RUN Rscript scripts/generate_monophyly.R family \
    && Rscript scripts/generate_monophyly.R order

COPY scripts/generate_fossils.R scripts/
RUN Rscript scripts/generate_fossils.R

FROM scratch AS files

COPY --from=build _assets /
COPY --from=build _fossils /
COPY --from=build _data /
COPY --from=build downloads /

CMD ["bash"]
