FROM r-base

RUN R -e 'pkg <- c("MonoPhy", "tidyverse", "glue", "future"); install.packages(pkg, repos = "https://cran.rstudio.com", Ncpu = parallel::detectCores()); if (!all(pkg) %in% installed.packages()) q(status = 1, save = "no")'

COPY data/* .
