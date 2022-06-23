## An alternative Overlap analysis pipeline for two sets of gene lists.
## setwd to a directory with two subfolders, a comparison set and an external/experimental set
## Rather than performing every pairwise overlap, this script will perform each comparator to external overlap
## Just like the original pairwise OL script, this will produce an output folder with overlap list text files,
## a summary csv file with stats, and venn diagrams of the overlaps. 

library(tidyverse)
library(glue)
library(VennDiagram)
library(rlist)
library(tictoc)

## RUN THESE FUNCTIONS BUT DON'T CHANGE
hyperGeomTest <- function(overlap, sample1, sample2, background) {
  sum(dhyper(overlap:sample1, sample2, background - sample2, sample1))
}

comparatorOverlap <- function(comparator, external, 
                              backgroundSize = 18000, analysis.name = 'Overlap', 
                              out.folder = 'geneOverlapLists') {
  
  # Create output folder if it doesn't exist already
  createNewDir(glue('{analysis.name}_outputs'))
  createNewDir(glue('{analysis.name}_outputs/{out.folder}'))

  len.comp = length(comparator)
  len.ext = length(external)
  
  # Calculate overlaps 
  out.OLgenes <- rlist::list.flatten(lapply(comparator, function(x) {
    list.map(external, intersect(.,x))}))
  
  # Write out overlap text files
  for (i in 1:length(out.OLgenes)) {
    filename = names(out.OLgenes)[i]
    write_tsv(as_tibble(out.OLgenes[[i]]), 
              file = glue('{analysis.name}_outputs/{out.folder}/{filename}.txt'), 
              col_names = F)
  }
  
  # Create and return summary overlap dataframe
  exts <- pivot_longer(data = as_tibble(sapply(comparator, function(x) rep(length(x), len.ext))), 
                       cols = everything(), 
                       names_to = 'list1', 
                       values_to = 'list1_length') %>%
    arrange(list1)
  
  list2 = rep(names(external), len.comp)
  list2_length = rep(unlist(lapply(external, length)), len.comp)
  overlap = unlist(lapply(out.OLgenes, length))
  
  tmp = cbind(exts, list2, list2_length, overlap) %>%
    as_tibble() %>%
    relocate(list2, .after = list1)

  out <- tmp %>% 
    mutate(across(list1_length:overlap, .fns = as.numeric),
           p.value = mapply(hyperGeomTest, overlap, list1_length, list2_length, backgroundSize),
           p.adj = p.adjust(p.value, method = 'bonf'),
           backgroundSize = backgroundSize)
  
  write_csv(out, file = glue("{analysis.name}_outputs/{analysis.name}_summary.csv"))
  
  out
}

output.venns.toPDF <-  function(ol.list, analysis.name) {
  venndir = if_else(ol.list[[6]] < 0.05, 
                    paste0(analysis.name, '_outputs/', analysis.name, '_vennSignif/'), 
                    paste0(analysis.name, '_outputs/', analysis.name, '_vennInsignif/'))
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
    out <- plyr::llply(files.vec, read_tsv, col_names = col.names) 
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

createNewDir <- function(dir.name) {
  if (!(dir.exists(dir.name))) {
    dir.create(dir.name)
  } 
}

## set working directory as the folder that contain the genelist files 
# setwd("/Users/Sammy/Desktop/UPenn/Projects/ASD_genes/Experiments/68.ReviewerResponseExperiment/recount3_forMouseStudies/SRP143453_Shank3mut/genelists")

analysis.name = 'test1'

files_ext <- list.files('externalLists/', full.names = T)
files_ext

files_comp <- list.files('testComparatorLists/', full.names = T)
files_comp

ext_lists <- readFilesToList(files_ext, file.type = 'txt')
ext_lists <- lapply(ext_lists, deframe)
comp_lists <- readFilesToList(files_comp, file.type = 'txt')
comp_lists <- lapply(comp_lists, deframe)

OLValuesDf = comparatorOverlap(comparator = comp_lists, external = ext_lists,
                               analysis.name = analysis.name,
                               backgroundSize = 13739)

# Draw venn diagrams of overlaps
venns.list <- split(OLValuesDf, rownames(OLValuesDf))
createNewDir(glue("{analysis.name}_outputs/{analysis.name}_vennInsignif"))
createNewDir(glue("{analysis.name}_outputs/{analysis.name}_vennSignif"))

lapply(venns.list, output.venns.toPDF, analysis.name = analysis.name)




