## Overlap analaysis pipeline for a folder containing the lists of genes to be overlapped.
library(tidyverse)

## set working directory as the folder that contain the genelist files 
setwd("/Users/Sammy/Desktop/UPenn/Projects/DOT1L/Dot1L_H3.3ChIP_analysis/genelists")

files <- dir()

## The empty dataframe created below has five columns: The first two name the lists that went into the 
## overlap, and columns 3-5 contain the number of genes in the ith list, jth list, and the number of genes
## shared by the ith and jth lists, respectively.

## The rows of the dataframe correspond to the number of overlaps performed, which can be calculated
## for a given number of lists, n: (n * (n-1)) / 2
n_files <- length(files)
n_minus_1 <- n_files - 1
listLength <- (n_files * n_minus_1) / 2

OLValuesDf <- data.frame(matrix(nrow = listLength, ncol = 5))
rowIdx <- 1
names(OLValuesDf) <- c("List1", "List2", "ith", "jth", "overlap")
dir.create(path = "../output")
for( i in files ){
  D <- data.frame(matrix(nrow=0,ncol=1))
  d <- read.table(i,header=F,stringsAsFactors=F)
  for( j in files) {
    if (j > i) {
      d2 <- read.table(j, header = F, stringsAsFactors = F)
      overlap <- d %>% intersect(d2)
      OL_df <- paste("D", i, j, sep = "_")
      OL_df<- assign(OL_df, rbind(D,overlap))
      myfile <- file.path("../output", 
                          paste0("OL", "_", i, "_", j, ".txt"))
      write_tsv(OL_df, path = myfile, col_names = F)
      OLValuesDf[rowIdx, ] <- c(i, j, count(d), count(d2), count(overlap))
      rowIdx <-  rowIdx + 1
    }
  }
}

## Change the totalPop to match the background expressed genelist for your experiment
totalPop <- 18000

hyperGeomTest <- function(overlap, sample1, sample2, background) {
  sum(dhyper(overlap:sample1, sample2, background - sample2, sample1))
}

OLValuesDf <- OLValuesDf %>% mutate(hyper_pval = mapply(hyperGeomTest, 
                                                        OLValuesDf$overlap, 
                                                        OLValuesDf$ith, 
                                                        OLValuesDf$jth, 
                                                        totalPop))
OLValuesDf$List1 <- gsub("\\..*", "", OLValuesDf$List1)
OLValuesDf$List2 <- gsub("\\..*", "", OLValuesDf$List2)

write_csv(OLValuesDf, path = "../OverlapsTable_HypergeometricTestPvals.csv")

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
