## Overlap analaysis pipeline for a folder containing the lists of genes to be overlapped.
## This will produce a folder of pairwise overlap text files, a table of the numbers
## in each of the original lists and overlap lists, as well as statistics, and 
## a folder containing pairwise venn diagrams (pdfs)
library(tidyverse)
library(glue)
library(VennDiagram)

## RUN THESE FUNCTIONS BUT DON'T CHANGE
hyperGeomTest <- function(overlap, sample1, sample2, background) {
  sum(dhyper(overlap:sample1, sample2, background - sample2, sample1))
}

calculateOverlaps <- function(gene.lists, backgroundSize = 18000, out.folder = 'geneOverlapLists') {
  if (!(dir.exists(glue('../{out.folder}')))) {
    dir.create(glue('../{out.folder}'))
  } 
  nms <- combn( names(gene.lists) , 2 , FUN = paste0 , collapse = "_" , simplify = FALSE )
  ll <- combn( gene.lists , 2 , simplify = FALSE )
  out.OLnum <- lapply( ll , function(x) length( intersect( x[[1]] , x[[2]] ) ) )
  out.OLgenes <- lapply( ll , function(x) intersect( x[[1]] , x[[2]] ) ) 
  out.OLgenes = setNames( out.OLgenes , nms )
  
  for (i in 1:length(out.OLgenes)) {
    filename = names(out.OLgenes)[i]
    write_tsv(as_tibble(out.OLgenes[[i]]), file = glue('../{out.folder}/{filename}.txt'), col_names = F)
  }
  
  List1 = unlist(lapply(ll, function(x) names(x)[1]))
  List2 = unlist(lapply(ll, function(x) names(x)[2]))
  ith = unlist(lapply(ll, function(x) length(x[[1]])))
  jth = unlist(lapply(ll, function(x) length(x[[2]])))
  overlap = unlist(lapply(out.OLnum, '[', 1))
  
  tmp = cbind(List1, List2, ith, jth, overlap) %>% as.data.frame()
  
  out <- tmp %>% 
    mutate(across(ith:overlap, .fns = as.numeric),
           hyper_pval = mapply(hyperGeomTest, overlap, ith, jth, backgroundSize),
           backgroundSize = backgroundSize)
  out
}

output.venns.toPDF <-  function(ol.list, analysis.name) {
  venndir = if_else(ol.list[[6]] < 0.05, 
                    paste0("../", analysis.name, '_vennSignif/'), 
                    paste0("../", analysis.name, '_vennInsignif/'))
  pdf(file = file.path(paste0(venndir, ol.list[[1]], "_", ol.list[[2]], ".pdf")))
  x <- draw.pairwise.venn(ol.list[[3]], ol.list[[4]], ol.list[[5]], 
                          category = c(ol.list[[1]], ol.list[[2]]),
                          cat.pos = c(0, 0))
  grid.draw(x)
  dev.off()
}


## set working directory as the folder that contain the genelist files 
setwd("../genelistsLimma/")

files <- dir()
files

gene.lists <- list()
sample.names <- c()

for (i in files) {
  sample = str_remove(i, '\\..*')
  tmp = read_tsv(i, col_names = F) %>% pull(X1)
  sample.names = c(sample.names, sample)
  gene.lists = c(gene.lists, list(tmp))
}

names(gene.lists) = sample.names

analysis.name = 'Limma_mouseNewStudies_OLs'

OLValuesDf <- calculateOverlaps(gene.lists, backgroundSize = 17987, out.folder = analysis.name)

write_csv(OLValuesDf, file = glue("../{analysis.name}_res.csv"))

venns.list <- split(OLValuesDf, rownames(OLValuesDf))
dir.create(glue("../{analysis.name}_vennInsignif"))
dir.create(glue("../{analysis.name}_vennSignif"))


lapply(venns.list, output.venns.toPDF, analysis.name = analysis.name)




