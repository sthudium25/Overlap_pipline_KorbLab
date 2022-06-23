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
           p.value = mapply(hyperGeomTest, overlap, ith, jth, backgroundSize),
           p.adjust = p.adjust(p.value, method = 'bonf'),
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

readFilesToList <- function(files.vec, file.type = c('txt', 'csv', 'xlsx'), col.names = FALSE, sheet=1) {
  stopifnot("Error: at least of of the files listed does not exist." = all(file.exists(files.vec)))
  if (file.type == 'txt') {
    out <- plyr::llply(files.vec, .fun=read_tsv, col_names = col.names) 
  } 
  if (file.type == 'csv') {
    out <- plyr::llply(files.vec, .fun=read_csv, col_names = col.names)
  } 
  if (file.type == 'xlsx') {
    out <- plyr::llply(files.vec, .fun=readxl::read_xlsx, 
                       col_names = col.names, sheet = sheet)
  } 
  names(out) <- str_remove(basename(files.vec), '\\..*')
  out
}

## set working directory as the folder that contain the genelist files 
setwd("path/to/genelists")

# Read in gene lists that you want to overlap
files <- list.files(path = '.', pattern = '.txt', full.names = T)
files
gene.lists <- readFilesToList(files, file.type = 'txt')

# Give this analysis a name
analysis.name = ''

# Compute and output pairwise overlaps + a summary file
OLValuesDf <- calculateOverlaps(gene.lists, backgroundSize = 13739, out.folder = analysis.name)
write_csv(OLValuesDf, file = glue("../{analysis.name}_res.csv"))

# Draw venn diagrams of overlaps
venns.list <- split(OLValuesDf, rownames(OLValuesDf))
dir.create(glue("../{analysis.name}_vennInsignif"))
dir.create(glue("../{analysis.name}_vennSignif"))

lapply(venns.list, output.venns.toPDF, analysis.name = analysis.name)




