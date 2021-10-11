#!/usr/bin/env Rscript

library(ape)
library(tidyverse)
library(dplyr)
library(readr)
library(ggplot2)
library(stringr)
library(future)
library(listenv)

width <- 1000 - (30 * 2)
height <- width * 3

slugify <- function(str) {
    str %>% str_replace_all("[^a-zA-Z0-9-]", "-") %>% str_replace_all("-+", "-") %>% str_replace("-$", "") %>% tolower()
}

tax %<-% read_csv("downloads/PFC_taxonomy.csv.xz", col_types = cols(.default = "c"))
wanted_ranks <- c("class", "subclass", "infraclass", "megacohort",
                  "supercohort", "cohort", "subcohort", "infracohort", "section",
                  "subsection", "division", "subdivision", "series", "superorder",
                  "order", "suborder", "infraorder", "family")

output <- list()
generate_rank_data <- function(df) {
    out <- list()
    out$species <- df$genus.species
    out
}

for (rank in wanted_ranks) {
    splat <- split(tax, tax[[rank]])
    for (named_rank in names(splat)) {
        output[[rank]][[named_rank]] <- generate_rank_data(splat[[named_rank]])
    }
}

tree <- read.tree("downloads/actinopt_12k_treePL.tre.xz")
fossil_nodes <- read_csv("downloads/fossil_data.csv")

fossil_nodes$idx <- seq_len(nrow(fossil_nodes))

# Add a new column `node` with the node number of that calibration
fossil_nodes <- group_by(fossil_nodes, group) %>% mutate(node = getMRCA(tree, c(left, right)))

dir.create("assets/img", recursive = TRUE)
png("assets/img/vertical_tree@3x.png", width = width * 3, height = height * 3)
plot(tree, show.tip.label = FALSE, no.margin = TRUE)
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
res <- fossil_nodes %>% mutate(x = lastPP$xx[node], y = lastPP$yy[node],
                               devx = grconvertX(x, to = "device"),
                               devy = grconvertY(y, to = "device"),
                               slug = slugify(fossil)) %>% ungroup()
dev.off()

png("assets/img/vertical_tree@2x.png", width = width * 2, height = height * 2)
plot(tree, show.tip.label = FALSE, no.margin = TRUE)
dev.off()

png("assets/img/vertical_tree@1x.png", width = width, height = height)
plot(tree, show.tip.label = FALSE, no.margin = TRUE)
dev.off()

output_cols <- list()
for(i in 1:nrow(res)) {
    row <- res[i,]
    left <- str_replace(row$left, "_", " ")
    right <- str_replace(row$right, "_", " ")
    for (rank in wanted_ranks) {
        splat <- split(tax, tax[[rank]])
        for (named_rank in names(splat)) {
            species <- output[[rank]][[named_rank]]$species
            found <- FALSE
            if(left %in% species && right %in% species) {
                output_cols[[rank]] = c(output_cols[[rank]], named_rank)
                found <- TRUE
                break
            }
        }
        if(!found) {
            output_cols[[rank]] = c(output_cols[[rank]], NA)
        }
    }
}

res <- cbind(res, as.data.frame(output_cols))

dir.create("_data/", recursive = TRUE)
res %>% write_csv("_data/fossil_data.csv")
res %>% transmute(clade = clade_pretty, fossil, left, right, min, max, locality, authority, age_authority) %>% write_csv("_data/fossil_pretty.csv")

mdpath <- "_fossils/"
dir.create(mdpath, recursive = T)

for (ii in 1:nrow(res)) {
    sink(file.path(mdpath, paste0(res$idx[ii], ".md")))
    cat("---\n")
    cat(paste0("title: ", res$fossil[ii], "\n"))
    cat(paste0("slug: ", res$slug[ii], "\n"))
    cat(paste0("description: Fossil calibration data for ", res$fossil[ii], ", an extinct species of fish. Includes taxonomy authority and locality references, and cross-references to living taxa.\n"))
    cat("\n---\n")
    sink(NULL)
}
