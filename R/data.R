#' Human Protein Protein Interaction Network
#'
#' This network is extracted from Human Integrated Protein-Protein Interaction rEference.
#' http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/index.php
#'
#' @format A data frame with two columns \code{g1} and \code{g2}. Each row
#'   indicates a pair of genes that have putative interations.
'hppi'

#' Mouse Protein Protein Interaction Network
#'
#' This network is extracted from the interaction database provided by BioGRID.
#' http://thebiogrid.org/
#'
#' @format A data frame with two columns \code{g1} and \code{g2}. Each row
#'   indicates a pair of genes that have putative interations.
'mppi'

#' Log FPKM table of distal mouse lung epithelial cells from two developmental stages (E14.5 and E16.5)
#'
#' http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE52583
#'
#' @format A matrix with rows being the genes sequenced and columns being the samples.
'log_fpkm'

#' Phenotype vector for distal mouse lung epithelial cells dataset
#'
#' @format A data frame with five variables: \code{year}, \code{sex},
#'   \code{name}, \code{n} and \code{prop} (\code{n} divided by total number
#'   of applicants in that year).
'phe'

#' A toy dataset
#'
#' @format A matrix with rows being the genes sequenced and columns being the samples.
'toy'