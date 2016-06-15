NMFEM
=====

Generate protein-protein interaction network modules from selected genes, using spin-glass algorithm.

Prerequisites
-------------

Please install the prerequisites before using NMFEM:

    install.packages(c("NMF", "igraph", "foreach", "dplyr"))
    source("https://bioconductor.org/biocLite.R")
    biocLite(c("topGO", "org.Mm.eg.db", "limma", "DESeq2"))

A parallelism implementation is required. For unix-based system,

    install.packages("doMC")

For Windows,

    install.packages("doParallel")


Usage
-----

`source('NMFEM.R')` then call `run_workflow`. This will generate top PPI modules using a list of seed genes provided.

The list of genes can be generated using `nmf` and `extractFeatures` functions provided in the `NMF` R-package.


Example
-------

`example.R` demonstrates a typical usage of NMFEM with seed genes generated from various methods. This script outputs the found modules to `*_modules.txt` files under the working folder.


Main Function
-------------

    run_workflow
      ( expr_matrix_
      , grouping_vec_
      , selected_genes_
      , ppi_edges_file_
      , eg_db_
      , ppi_edges_file_header_ = T
      , verbose_level_         = 1
      , n_top_modules_         = 5
      , n_threads_             = 1
      , min_size_              = 10
      , max_size_              = 100 )

**NOTE:** When running this function, please make sure that the row names of `expr_matrix_`, the values in `selected_genes_`, and the records in `ppi_edges_file_` all belong to the same set of symbols. Gene symbols from different speices may have different capitalizations.

Parameters
----------

`expr_matrix_`:

The expression matrix, each row is a gene and each column is a sample. The row names should
be gene symbols and the column names should be sample ids. Use the correct capticalization
for gene symbols. For mouse genes,only the first letter is capitalized (e.g. Tp53); for human
genes, all letters are capitalized (e.g. TP53).

`grouping_vec_` :

A factor vector, indicating the groupping. Its length should be
the same as the number of columns in `expr_matrix_`. Example:

    factor(c(1,1,1,2,2,2,2,2,2))

`selected_genes_` :

A character vector. The genes that are used as seed genes in the FEM algorithm. Typically
generated from another method such as NMF.

`ppi_edges_file_` :

The file containing the protein-protein interaction (PPI) network file.
The format is examplified in `mppi.txt` for mouse genes, and `hppi.txt` for human genes.

`eg_db_` :

The ID conversion database. Such as `'org.Mm.eg.db'` or `'org.Hs.eg.db'`.

`ppi_edges_file_header_` :

Whether the PPI file has a header row. Default: `True`.

`verbose_level_` :

How much information is printed. Either 0, 1 or 2, default: 1.

`n_top_modules_` :

Number of top modules to pick. Default: 5.

`n_threads_` :

Number of threads to use. Default: 1.

`min_size_` :

Module size lower cut-off. Default: 10.

`max_size_` :

Module size upper cut-off. Default: 100.

Example Code
------------

    library(org.Mm.eg.db)
    library(foreach)
    library(doMC)
    library(topGO)
    library(DESeq2)
    library(limma)
    library(igraph)
    library(dplyr)

    source('NMFEM.R')

    se <- readRDS('epitSE2.RDS')
    rw <- assays(se)$fpkm
    min0 <- min(setdiff(rw, 0))
    rwl <- log(rw/min0 + 1)

    phe <- colData(se)$phe

    leading_genes <- list()
    leading_genes[['NMF']] <- readRDS('leading_genes/leadingGenesNMF.RDS')
    leading_genes[['DESeq2']] <- readRDS('leading_genes/leadingGenesDESeq2.RDS')
    leading_genes[['EdgeR']] <- readRDS('leading_genes/leadingGenesEdgeR.RDS')
    leading_genes[['Monocle']] <- readRDS('leading_genes/leadingGenesMonocle.RDS')
    leading_genes[['B2013']] <- readRDS('leading_genes/leadingGenesB2013.RDS')
    leading_genes[['SCDE']] <- readRDS('leading_genes/leadingGenesSCDE.RDS')
    leading_genes[['SemiNMF']] <- readRDS('leading_genes/leadingGenesSemiNMF.RDS')
    leading_genes[['MAST']] <- readRDS('leading_genes/leadingGenesMAST.RDS')

    for (method_index in 1:length(leading_genes)) {
      final_tb <- run_workflow(rwl, phe, leading_genes[[method_index]], 'mppi.txt', 'org.Mm.eg.db', verbose_level_ = 2, n_threads_ = 1)
      
      write.table(final_tb, paste0(names(leading_genes)[method_index], '_modules.txt'), sep='\t', row.names=F, col.names=T, quote=F)
    }
