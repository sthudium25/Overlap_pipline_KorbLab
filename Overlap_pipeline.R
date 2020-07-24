## Overlap analaysis pipeline for a folder containing the lists of genes to be overlapped.
library(tidyverse)
## set working directory as the folder that contain the genelist files 
setwd("/Users/Sammy/Desktop/UPenn/Projects/Overlap_pipline/test_genelists")
files <- dir()

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
    }
  }
}
