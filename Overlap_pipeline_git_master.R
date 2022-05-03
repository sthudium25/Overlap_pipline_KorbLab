## Overlap analaysis pipeline for a folder containing the lists of genes to be overlapped.
## This will produce a folder of pairwise overlap text files, a table of the numbers
## in each of the original lists and overlap lists, as well as statistics, and 
## a folder containing pairwise venn diagrams (pdfs)
library(tidyverse)
library(glue)

## RUN THESE FUNCTIONS BUT DON'T CHANGE
hyperGeomTest <- function(overlap, sample1, sample2, background) {
  sum(dhyper(overlap:sample1, sample2, background - sample2, sample1))
}

calculateOverlaps <- function(gene.lists, backgroundSize = 18000) {
  nms <- combn( names(gene.lists) , 2 , FUN = paste0 , collapse = "_" , simplify = FALSE )
  ll <- combn( gene.lists , 2 , simplify = FALSE )
  out.OLnum <- lapply( ll , function(x) length( intersect( x[[1]] , x[[2]] ) ) )
  out.OLgenes <- lapply( ll , function(x) intersect( x[[1]] , x[[2]] ) ) 
  out.OLgenes = setNames( out.OLgenes , nms )
  
  for (i in 1:length(out.OLgenes)) {
    filename = names(out.OLgenes)[i]
    write_tsv(as_tibble(out.OLgenes[[i]]), file = glue('../geneOverlapLists/{filename}.txt'), col_names = F)
  }
  
  List1 = unlist(lapply(ll, function(x) names(x)[1]))
  List2 = unlist(lapply(ll, function(x) names(x)[2]))
  ith = unlist(lapply(ll, function(x) length(x[[1]])))
  jth = unlist(lapply(ll, function(x) length(x[[2]])))
  overlap = unlist(lapply(out.OLnum, '[', 1))
  
  tmp = cbind(List1, List2, ith, jth, overlap) %>% as.data.frame()
  
  out <- tmp2 %>% 
    mutate(across(ith:overlap, .fns = as.numeric),
           hyper_pval = mapply(hyperGeomTest, overlap, ith, jth, backgroundSize))
  out
}
## set working directory as the folder that contain the genelist files 
setwd("/Users/Sammy/Desktop/UPenn/Projects/ASD_genes/Experiments/68.ReviewerResponseExperiment/recount3_forMouseStudies/OverlapAnalysis_wSignatureLists/genelists")

files <- dir()

gene.lists <- list()
sample.names <- c()

for (i in files) {
  sample = str_remove(i, '\\..*')
  tmp = read_tsv(i, col_names = F) %>% pull(X1)
  sample.names = c(sample.names, sample)
  gene.lists = c(gene.lists, list(tmp))
}

names(gene.lists) = sample.names

## The empty dataframe created below has five columns: The first two name the lists that went into the 
## overlap, and columns 3-5 contain the number of genes in the ith list, jth list, and the number of genes
## shared by the ith and jth lists, respectively.

dir.create(path = "../geneOverlapLists")

OLValuesDf <- calculateOverlaps(gene.lists, 18000)


write_csv(OLValuesDf, path = "../OverlapsResultsTable.csv")

library(VennDiagram)
list1 <- split(OLValuesDf, rownames(OLValuesDf))
dir.create("../vennOutput")
setwd("../vennOutput")
index <- 1
for (list in list1) {
  pdf(file = paste0(list[[1]], "_", list[[2]], ".pdf"))
  x <- draw.pairwise.venn(list[[3]], list[[4]], list[[5]], 
                          category = c(list[[1]], list[[2]]),
                          cat.pos = c(0, 0))
  grid.draw(x)
  dev.off()
  index <- index + 1
}
