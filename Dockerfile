FROM rocker/tidyverse:3.5.2

RUN install2.r --error \
    --ncpus -1 \
    ape \
    future \
    phangorn \
    progress \
    && rm -rf /tmp/*

COPY scripts/* scripts/
COPY downloads/* downloads/

RUN Rscript scripts/generate_taxonomy.R \
    && Rscript scripts/generate_monophyly.R family \
    && Rscript scripts/generate_monophyly.R order \
    && Rscript scripts/generate_fossils.R \
    && rm -rf scripts/*

CMD ["bash"]
