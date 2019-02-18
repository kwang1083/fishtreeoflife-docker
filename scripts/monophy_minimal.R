requireNamespace("phangorn")
# from phytools
getDescendants<-function(tree,node,curr=NULL){
	if(!inherits(tree,"phylo")) stop("tree should be an object of class \"phylo\".")
	if(is.null(curr)) curr<-vector()
	daughters<-tree$edge[which(tree$edge[,1]==node),2]
	curr<-c(curr,daughters)
	if(length(curr)==0&&node<=Ntip(tree)) curr<-node
	w<-which(daughters>Ntip(tree))
	if(length(w)>0) for(i in 1:length(w)) 
		curr<-getDescendants(tree,daughters[w[i]],curr)
	return(curr)
}

# from MonoPhy package
AssessMonophyly <- function (tree, taxonomy = NULL, verbosity = 5, outliercheck = TRUE, outlierlevel = 0.5) {
    if (!is.rooted(tree)) {
        stop("Phylogeny must be rooted!")
    }
    if (is.null(taxonomy)) {
        for (i in 1:length(tree$tip.label)) {
            if (grepl(("_| "), tree$tip.label[i]) == FALSE) {
                stop("Tip labels do not contain underscore separating genus name from species epithet!")
            }
        }
        f <- function(s) strsplit(s, ("_| "))[[1]][1]
        split.taxa <- sapply(tree$tip.label, f)
        taxa <- as.vector(unique(split.taxa))
        taxsets <- c("taxa")
    }
    else {
        if (length(taxonomy[, 1]) != length(tree$tip.label)) {
            stop("Number of rows of taxonomy file is not equal to number of taxa (note: if your table has a header, you should specify header=TRUE when importing)!")
        }
        if (length(taxonomy[1, ]) < 2) {
            stop("Taxonomy file needs at least 2 columns: tip labels and taxonomic group!")
        }
        taxchecktree <- c()
        for (itaxcheck in 1:length(tree$tip.label)) {
            taxchecktree <- c(taxchecktree, tree$tip.label[itaxcheck] %in% 
                              taxonomy[, 1])
        }
        taxintruderstree <- c()
        if ("FALSE" %in% taxchecktree) {
            positionstree <- grep("FALSE", taxchecktree)
            taxintruderstree <- c(taxintruderstree, tree$tip.label[positionstree])
            message(paste("\n"), appendLF = TRUE)
            message(paste("Tip-labels which do not occur in taxonomy file:", 
                          "[", length(taxintruderstree), "/", length(tree$tip.label), 
                          "]", collapse = " "), appendLF = TRUE)
            message(paste(taxintruderstree, collapse = ", "), 
                    appendLF = TRUE)
            message(paste("\n"), appendLF = TRUE)
        }
        taxcheckfile <- c()
        for (itaxcheck2 in 1:length(taxonomy[, 1])) {
            taxcheckfile <- c(taxcheckfile, taxonomy[itaxcheck2, 
                              1] %in% tree$tip.label)
        }
        taxintrudersfile <- c()
        if ("FALSE" %in% taxcheckfile) {
            positionsfile <- grep("FALSE", taxcheckfile)
            taxintrudersfile <- c(taxintrudersfile, as.character(taxonomy[positionsfile, 1]))
            message(paste("Taxon names in file which do not occur in tip-labels of tree:", 
                          "[", length(taxintrudersfile), "/", length(taxonomy[, 1]), "]", collapse = " "), appendLF = TRUE)
            message(paste(taxintrudersfile, collapse = ", "), appendLF = TRUE)
            message(paste("\n"), appendLF = TRUE)
        }
        if ("FALSE" %in% (taxchecktree)) {
            stop("The taxon names of tree and taxonfile do not match (see above)!")
        }
        if ("FALSE" %in% (taxcheckfile)) {
            stop("The taxon names of tree and taxonfile do not match (see above)!")
        }
        taxsets <- list()
        for (jtax in 1:(length(taxonomy[1, ]) - 1)) {
            nametax <- paste("taxa", jtax, sep = "")
            tmp <- as.vector(unique(taxonomy[, (jtax) + 1]))
            taxsets[[nametax]] <- tmp
        }
    }
    finallist <- list()
    for (ifullround in 1:length(taxsets)) {
        if (is.null(taxonomy)) {
            taxa <- taxa
        }
        else {
            taxa <- unlist(taxsets[ifullround])
            taxa <- taxa[!taxa %in% c("unknown", NA)]
        }
        intruder.genus <- list()
        intruder.genus.full <- c()
        intruder.species <- list()
        intruder.species.full <- c()
        intruder.names <- c()
        if (outliercheck == TRUE) {
            outlist.summary <- matrix(NA, nrow = 6, ncol = 3)
            dfheaders <- c("Taxon", "Monophyly", "MRCA", "#Tips", 
                           "Delta-Tips", "#Intruders", "Intruders", "#Outliers", 
                           "Outliers")
            outlist <- matrix(NA, nrow = length(taxa), ncol = 9)
            outlier.species <- list()
            outlier.species.full <- c()
            outlier.names <- c()
        }
        else {
            outlist.summary <- matrix(NA, nrow = 5, ncol = 3)
            dfheaders <- c("Taxon", "Monophyly", "MRCA", "#Tips", 
                           "Delta-Tips", "#Intruders", "Intruders")
            outlist <- matrix(NA, nrow = length(taxa), ncol = 7)
        }
        tip.states.matrix <- matrix(NA, nrow = length(tree$tip.label), 
                                    ncol = 3)
        for (i in 1:length(taxa)) {
            if (is.null(taxonomy)) {
                ancnode <- getMRCA(tree, tip = c(tree$tip.label[c(grep(paste("^", 
                                                                             taxa[i], "_", sep = ""), tree$tip.label))]))
            }
            else {
                subtips <- subset(taxonomy, as.character(taxonomy[, 
                                                         (ifullround + 1)]) == as.character(taxa[i]))
                subtipsnr <- c()
                for (sbts in 1:nrow(subtips)) {
                    sbtname <- subtips[sbts, 1]
                    sbtnr <- which(tree$tip.label == sbtname)
                    subtipsnr <- c(subtipsnr, sbtnr)
                }
                ancnode <- getMRCA(tree, tip = c(subtipsnr))
            }
            if (length(ancnode) == 0) {
                if (outliercheck == TRUE) {
                    outlist[i, ] <- c(taxa[i], "Monotypic", "NA", 
                                      1, "NA", "NA", "", "NA", "")
                }
                else {
                    outlist[i, ] <- c(taxa[i], "Monotypic", "NA", 
                                      1, "NA", "NA", "")
                }
            }
            else {
                anctips <- getDescendants(tree, ancnode)
                ancnames <- tree$tip.label[c(anctips)]
                ancnames <- ancnames[!is.na(ancnames)]
                if (is.null(taxonomy)) {
                    taxtips <- tree$tip.label[c(grep(paste("^", 
                                                           taxa[i], "_", sep = ""), tree$tip.label))]
                }
                else {
                    taxtips <- subtips[, 1]
                }
                if (length(ancnames) == length(taxtips)) {
                    if (outliercheck == TRUE) {
                        outlist[i, ] <- c(taxa[i], "Yes", ancnode, 
                                          length(taxtips), "0", "0", "", "NA", "")
                    }
                    else {
                        outlist[i, ] <- c(taxa[i], "Yes", ancnode, 
                                          length(taxtips), "0", "0", "")
                    }
                }
                else {
                    intruder.tips <- setdiff(ancnames, taxtips)
                    if (is.null(taxonomy)) {
                        f2 <- function(s) strsplit(s, ("_| "))[[1]][1]
                        split.taxa2 <- sapply(intruder.tips, f2)
                        intruder.taxa <- as.vector(unique(split.taxa2))
                    }
                    else {
                        subtaxa <- c()
                        for (j in 1:length(intruder.tips)) {
                            subtaxon <- rbind(subset(taxonomy, taxonomy[, 
                                                     1] == intruder.tips[j]))
                            subtaxa <- rbind(subtaxa, subtaxon)
                        }
                        intruder.taxa <- as.vector(unique(subtaxa[, 
                                                          ifullround + 1]))
                    }
                    outlier.tips <- c()
                    if (outliercheck == TRUE) {
                        if (length(phangorn::Children(tree, ancnode)) > 2) {
                            outlier.tips <- c()
                        }
                        else {
                            tiplevels <- length(taxtips)/length(ancnames)
                            if (tiplevels < outlierlevel) {
                                start.node <- ancnode
                                while (tiplevels < outlierlevel) {
                                    subtaxtips <- c()
                                    subancnames <- c()
                                    parent.node <- start.node
                                    daughter.nodes <- phangorn::Children(tree, parent.node)
                                    if (length(daughter.nodes) > 2) {
                                        anctips1 <- getDescendants(tree, 
                                                                   parent.node)
                                        ancnames1 <- tree$tip.label[c(anctips1)]
                                        subancnames <- ancnames1[!is.na(ancnames1)]
                                        subtaxtips <- intersect(taxtips, 
                                                                subancnames)
                                        start.node <- parent.node
                                        break
                                    }
                                    daughter1 <- daughter.nodes[1]
                                    daughter2 <- daughter.nodes[2]
                                    anctips1 <- getDescendants(tree, daughter1)
                                    ancnames1 <- tree$tip.label[c(anctips1)]
                                    ancnames1 <- ancnames1[!is.na(ancnames1)]
                                    taxtips1 <- intersect(taxtips, ancnames1)
                                    anctips2 <- getDescendants(tree, daughter2)
                                    ancnames2 <- tree$tip.label[c(anctips2)]
                                    ancnames2 <- ancnames2[!is.na(ancnames2)]
                                    taxtips2 <- intersect(taxtips, ancnames2)
                                    nodechoice <- which(c(length(taxtips1), 
                                                          length(taxtips2)) == max(c(length(taxtips1), 
                                                                                     length(taxtips2))))
                                    if (length(nodechoice) > 1) {
                                        nodechoice2 <- which(c((length(taxtips1)/length(anctips1)), 
                                                               (length(taxtips2)/length(anctips2))) == 
                                        max(c((length(taxtips1)/length(anctips1)), 
                                              (length(taxtips2)/length(anctips2)))))
                                        if (length(nodechoice2) > 1) {
                                            subtaxtips <- c(taxtips1, taxtips2)
                                            subancnames <- c(ancnames1, ancnames2)
                                            start.node <- parent.node
                                            break
                                        }
                                        else if (nodechoice2 == 1) {
                                            subtaxtips <- taxtips1
                                            subancnames <- ancnames1
                                            start.node <- daughter1
                                        }
                                        else if (nodechoice2 == 2) {
                                            subtaxtips <- taxtips2
                                            subancnames <- ancnames2
                                            start.node <- daughter2
                                        }
                                    }
                                    else if (nodechoice == 1) {
                                        subtaxtips <- taxtips1
                                        subancnames <- ancnames1
                                        start.node <- daughter1
                                    }
                                    else if (nodechoice == 2) {
                                        subtaxtips <- taxtips2
                                        subancnames <- ancnames2
                                        start.node <- daughter2
                                    }
                                    tiplevels <- length(subtaxtips)/length(subancnames)
                                }
                                if (tiplevels < 1 & length(daughter.nodes) <= 
                                    2) {
                                    EDtaxtips1 <- c()
                                    EDtaxtips2 <- c()
                                    repeat {
                                        EDparent.node <- start.node
                                        if (EDparent.node <= length(tree$tip.label)) {
                                            subtaxtips <- tree$tip.label[EDparent.node]
                                            subancnames <- tree$tip.label[EDparent.node]
                                            start.node <- EDparent.node
                                            break
                                        }
                                        EDdaughter.nodes <- phangorn::Children(tree, 
                                                                     EDparent.node)
                                        EDdaughter1 <- EDdaughter.nodes[1]
                                        EDdaughter2 <- EDdaughter.nodes[2]
                                        EDanctips1 <- getDescendants(tree, 
                                                                     EDdaughter1)
                                        EDancnames1 <- tree$tip.label[c(EDanctips1)]
                                        EDancnames1 <- EDancnames1[!is.na(EDancnames1)]
                                        EDtaxtips1 <- intersect(taxtips, 
                                                                EDancnames1)
                                        EDanctips2 <- getDescendants(tree, 
                                                                     EDdaughter2)
                                        EDancnames2 <- tree$tip.label[c(EDanctips2)]
                                        EDancnames2 <- EDancnames2[!is.na(EDancnames2)]
                                        EDtaxtips2 <- intersect(taxtips, 
                                                                EDancnames2)
                                        if (length(EDtaxtips1) == 0) {
                                            start.node <- EDdaughter2
                                        }
                                        else if (length(EDtaxtips2) == 0) {
                                            start.node <- EDdaughter1
                                        }
                                        else {
                                            subtaxtips <- c(EDtaxtips1, EDtaxtips2)
                                            subancnames <- c(EDancnames1, EDancnames2)
                                            start.node <- EDparent.node
                                            break
                                        }
                                    }
                                }
                                outlier.tips <- setdiff(taxtips, subtaxtips)
                                if (length(outlier.tips) != 0) {
                                    outlier.species <- c(outlier.species, 
                                                         list(Tips = outlier.tips))
                                    outlier.species.full <- c(outlier.species.full, 
                                                              outlier.tips)
                                    outlier.names <- c(outlier.names, taxa[i])
                                }
                                intruder.tips <- setdiff(subancnames, 
                                                         subtaxtips)
                                if (is.null(taxonomy)) {
                                    f2 <- function(s) strsplit(s, ("_| "))[[1]][1]
                                    split.taxa2 <- sapply(intruder.tips, 
                                                          f2)
                                    intruder.taxa <- as.vector(unique(split.taxa2))
                                }
                                else {
                                    subtaxa <- c()
                                    for (j in 1:length(intruder.tips)) {
                                        subtaxon <- rbind(subset(taxonomy, 
                                                                 taxonomy[, 1] == intruder.tips[j]))
                                        subtaxa <- rbind(subtaxa, subtaxon)
                                    }
                                    intruder.taxa <- as.vector(unique(subtaxa[, 
                                                                      ifullround + 1]))
                                }
                            }
                        }
                        if (length(intruder.taxa) != 0) {
                            intruder.genus <- c(intruder.genus, list(Taxa = intruder.taxa))
                            intruder.genus.full <- c(intruder.genus.full, 
                                                     intruder.taxa)
                            intruder.species <- c(intruder.species, 
                                                  list(Tips = intruder.tips))
                            intruder.species.full <- c(intruder.species.full, 
                                                       intruder.tips)
                            intruder.names <- c(intruder.names, taxa[i])
                        }
                        if (outliercheck == TRUE) {
                            if (length(intruder.taxa) <= verbosity) {
                                if (length(outlier.tips) <= verbosity) {
                                    outlist[i, ] <- c(taxa[i], "No", ancnode, 
                                                      length(taxtips), (length(ancnames) - 
                                                                        length(taxtips)), length(intruder.taxa), 
                                                      paste(intruder.taxa, collapse = ", "), 
                                                      length(outlier.tips), paste(outlier.tips, 
                                                                                  collapse = ", "))
                                }
                                else {
                                    outlist[i, ] <- c(taxa[i], "No", ancnode, 
                                                      length(taxtips), (length(ancnames) - 
                                                                        length(taxtips)), length(intruder.taxa), 
                                                      paste(intruder.taxa, collapse = ", "), 
                                                      length(outlier.tips), paste(outlier.tips[1], 
                                                                                  "and", (length(outlier.tips) - 
                                                                                          1), "more.", collapse = ", "))
                                }
                            }
                            else {
                                if (length(outlier.tips) <= verbosity) {
                                    outlist[i, ] <- c(taxa[i], "No", ancnode, 
                                                      length(taxtips), (length(ancnames) - 
                                                                        length(taxtips)), length(intruder.taxa), 
                                                      paste(intruder.taxa[1], "and", (length(intruder.taxa) - 
                                                                                      1), "more.", collapse = ", "), 
                                                      length(outlier.tips), paste(outlier.tips, 
                                                                                  collapse = ", "))
                                }
                                else {
                                    outlist[i, ] <- c(taxa[i], "No", ancnode, 
                                                      length(taxtips), (length(ancnames) - 
                                                                        length(taxtips)), length(intruder.taxa), 
                                                      paste(intruder.taxa[1], "and", (length(intruder.taxa) - 
                                                                                      1), "more.", collapse = ", "), 
                                                      length(outlier.tips), paste(outlier.tips[1], 
                                                                                  "and", (length(outlier.tips) - 
                                                                                          1), "more.", collapse = ", "))
                                }
                            }
                        }
                        else {
                            if (length(intruder.taxa) <= verbosity) {
                                outlist[i, ] <- c(taxa[i], "No", ancnode, 
                                                  length(taxtips), (length(ancnames) - 
                                                                    length(taxtips)), length(intruder.taxa), 
                                                  paste(intruder.taxa, collapse = ", "))
                            }
                            else {
                                outlist[i, ] <- c(taxa[i], "No", ancnode, 
                                                  length(taxtips), (length(ancnames) - 
                                                                    length(taxtips)), length(intruder.taxa), 
                                                  paste(intruder.taxa[1], "and", (length(intruder.taxa) - 
                                                                                  1), "more.", collapse = ", "))
                            }
                        }
                    }
                }
            }
        }
        intruder.genus.all <- unique(intruder.genus.full)
        intruder.species.all <- unique(intruder.species.full)
        outlier.species.all <- c()
        if (outliercheck == TRUE) {
            outlier.species.all <- unique(outlier.species.full)
        }
        outframe <- data.frame(outlist)
        names(outframe) <- dfheaders
        rownames(outframe) <- outframe[, 1]
        outframe[, 1] <- NULL
        names(intruder.genus) <- intruder.names
        names(intruder.species) <- intruder.names
        if (is.null(taxonomy)) {
            tip.states.matrix[, 1] <- tree$tip.label
            f3 <- function(s) strsplit(s, ("_| "))[[1]][1]
            tip.states.matrix[, 2] <- sapply(tip.states.matrix[, 
                                             1], f3)
        }
        else {
            tip.states.matrix[, 1] <- as.vector(taxonomy[, 1])
            tip.states.matrix[, 2] <- as.vector(taxonomy[, ifullround + 
                                                1])
        }
        for (i in 1:length(tree$tip.label)) {
            if (tip.states.matrix[i, 1] %in% intruder.species.all == 
                TRUE) {
                tip.states.matrix[i, 3] <- "Intruder"
            }
        else if (tip.states.matrix[i, 1] %in% outlier.species.all == 
                 TRUE) {
            tip.states.matrix[i, 3] <- "Outlier"
        }
    else if (outframe[tip.states.matrix[i, 2], "Monophyly"] == 
             "Monotypic") {
        tip.states.matrix[i, 3] <- "Monophyletic"
    }
else if (outframe[tip.states.matrix[i, 2], "Monophyly"] == 
         "Yes") {
    tip.states.matrix[i, 3] <- "Monophyletic"
}
            else if (outframe[tip.states.matrix[i, 2], "Monophyly"] == 
                     "No") {
                tip.states.matrix[i, 3] <- "Non-Monophyletic"
            }
        else if (tip.states.matrix[i, 2] == "unknown" | is.na(tip.states.matrix[i, 
                                                              2])) {
            tip.states.matrix[i, 3] <- "unknown"
        }
        }
        tip.states.frame <- as.data.frame(tip.states.matrix)
        colnames(tip.states.frame) <- c("Tip", "Taxon", "Status")
        if (outliercheck == TRUE) {
            names(outlier.species) <- outlier.names
            outlist.summary[, 1] <- c("Total", "Monophyletic", 
                                      "Non-Monophyletic", "Monotypic", "Intruder", 
                                      "Outlier")
        }
        else {
            outlist.summary[, 1] <- c("Total", "Monophyletic", 
                                      "Non-Monophyletic", "Monotypic", "Intruder")
        }
        counttable <- table(outframe[, "Monophyly"])
        countframe <- as.data.frame(counttable)
        rownames(countframe) <- countframe[, 1]
        countframe[, 1] <- NULL
        mono.count <- c()
        nonmono.count <- c()
        for (itaxcount1 in 1:length(taxa)) {
            if (outlist[itaxcount1, 2] == "Yes") {
                mono.count <- c(mono.count, as.numeric(outlist[itaxcount1, 
                                                       4]))
            }
        }
        for (itaxcount2 in 1:length(taxa)) {
            if (outlist[itaxcount2, 2] == "No") {
                nonmono.count <- c(nonmono.count, as.numeric(outlist[itaxcount2, 
                                                             4]))
            }
        }
        if (outliercheck == TRUE) {
            outlist.summary[, 2] <- c(length(taxa), countframe["Yes", 
                                      "Freq"], countframe["No", "Freq"], countframe["Monotypic", 
                                      "Freq"], length(intruder.genus.all), length(outlier.names))
            outlist.summary[, 3] <- c(length(tree$tip.label), 
                                      sum(mono.count), sum(nonmono.count), countframe["Monotypic", 
                                                                                      "Freq"], length(intruder.species.all), length(outlier.species.all))
        }
        else {
            outlist.summary[, 2] <- c(length(taxa), countframe["Yes", 
                                      "Freq"], countframe["No", "Freq"], countframe["Monotypic", 
                                      "Freq"], length(intruder.genus.all))
            outlist.summary[, 3] <- c(length(tree$tip.label), 
                                      sum(mono.count), sum(nonmono.count), countframe["Monotypic", 
                                                                                      "Freq"], length(intruder.species.all))
        }
        outframe.summary <- data.frame(outlist.summary)
        rownames(outframe.summary) <- outframe.summary[, 1]
        outframe.summary[, 1] <- NULL
        colnames(outframe.summary) <- c("Taxa", "Tips")
        if (outliercheck == TRUE) {
            outputlist <- list(IntruderTaxa = intruder.genus, 
                               IntruderTips = intruder.species, OutlierTaxa = outlier.names, 
                               OutlierTips = outlier.species, result = outframe, 
                               summary = outframe.summary, TipStates = tip.states.frame)
        }
        else {
            outputlist <- list(IntruderTaxa = intruder.genus, 
                               IntruderTips = intruder.species, result = outframe, 
                               summary = outframe.summary, TipStates = tip.states.frame)
        }
        if (is.null(taxonomy)) {
            nameout <- "Genera"
        }
        else if (!is.null(taxonomy)) {
            if (colnames(taxonomy)[ifullround] == paste("V", 
                                                        ifullround, sep = "")) {
                nameout <- paste("Taxlevel", ifullround, sep = "_")
            }
        else {
            nameout <- colnames(taxonomy)[ifullround + 1]
        }
        }
        finallist[[nameout]] <- outputlist
    }
    finallist
}

