NMFEM
=====

Generate protein-protein interaction network modules from selected genes, using spin-glass algorithm.

Installation
------------

    install.packages(c("dplyr", "ggplot2", "devtools", "igraph"))
    source("https://bioconductor.org/biocLite.R")
    biocLite(c("topGO", "org.Mm.eg.db", "limma", "DESeq2", "clusterProfiler", "org.Hs.eg.db"))
    devtools::install_github('lanagarmire/NMFEM')

Documentation
-------------

Help for the main functions can be found at:

    ?nmf_subpopulation
    ?spinglass_procedure

Example Code
------------

    library(NMFEM)
    subp <- nmf_subpopulation(toy)
    res <- spinglass_procedure(fpkm, phe, leading_genes, mppi, 'mouse')
