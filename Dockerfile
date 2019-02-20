FROM rocker/tidyverse:3.5.2

RUN install2.r --error \
    --ncpus -1 \
    ape \
    future \
    phangorn \
    && rm -rf /tmp/*

COPY scripts/* scripts/
COPY downloads/* downloads/

RUN Rscript scripts/generate_taxonomy.R

RUN Rscript scripts/generate_monophyly.R family
    && Rscript scripts/generate_monophyly.R order

RUN Rscript scripts/generate_fossils.R

RUN rm -rf scripts/*

CMD ["bash"]
