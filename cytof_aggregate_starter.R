library(flowCore)
library(tidyverse)

## modify paths to fit your local dir structure
dir_data <- "data/acute_covid_cytof"
files <- list.files(dir_data, pattern=".fcs")

set.seed(0)
drop_channels <- c()
nmax <- 1e5
data_list <- sapply(paste0(dir_data, "/", files), function(path) {
  message(path)
  ff <- read.FCS(path)
  mat <- asinh(exprs(ff)/5) # change this to e.g. 500 for flow data
  
  nice_names <- ff@parameters@data$desc %>% 
    str_split("_") %>% 
    sapply("[[",2) # this is splitting on _ and taking the second piece
  ## change as needed for flow data
  
  colnames(mat) <- nice_names
  
  sel_cells <- sample(nrow(mat), min(nrow(mat), nmax))
  
  return(mat[sel_cells,setdiff(colnames(mat),drop_channels)])
})

data <- do.call(rbind, data_list)
n_each <- sapply(data_list, nrow)
samples <- rep(files, n_each)

save(data, samples, file=paste0(dir_data, "/data.rda"))




