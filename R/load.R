library(taigr)
library(magrittr)
library(celllinemapr)
library(tidyverse)
library(useful)

#' Loads omics features from taiga 
#'
#' @return a vector of gene cell line x gene effect for the given gene.
#' 
#' @export
#' 
load_achilles <- function(gene){
  xpr <- load.from.taiga(data.name='depmap-a0ab',data.file='Achilles_gene_effect')
  rownames(xpr) %<>% arxspan.to.ccle()
  colnames(xpr) <- str_split_fixed(string = colnames(xpr),pattern = " ",n = 2)[,1]
  xpr[,gene]
}

#' Loads omics features from taiga 
#'
#' @return a named list containing RNA expression, RPPA, Micro RNA, and Proteomics
#' 
#' @export
#' 
load_omics <- function(){
  ge <- load.from.taiga(data.name='mts011-9d70',data.file='ge')
  rppa <- rppa <- load.from.taiga(data.name='mts011-9d70',data.file='rppa')
  mirna <- load.from.taiga(data.name='mts011-9d70',data.file='mirna')
  prot <- load.from.taiga(data.name='total-proteome--5c50',data.file='normalized_protein_abundance')
  rownames(prot) %<>% arxspan.to.ccle()
  colnames(prot) <- str_c("PROT_",str_split_fixed(string = colnames(prot),pattern = " ",n = 2)[,1])
  list("RNA expression" = ge, "RPPA" = rppa, "Micro RNA" = mirna,"Proteomics" = prot)
}
