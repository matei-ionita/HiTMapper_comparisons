library(flowCore)
library(tidyverse)
library(wadeTools)

## modify paths to fit your local dir structure
dir_data <- "data/acute_covid_flow"
files <- list.files(dir_data, pattern=".fcs")

set.seed(0)
drop_channels <- c()
nmax <- 1e5 # number of rows subsetted
data_list <- sapply(paste0(dir_data, "/", files), function(path) {
  message(path)
  ff <- read.FCS(path)
  not_fluo <- grep("Time|FSC|SSC", colnames(ff))
  fluo <- colnames(ff)[-not_fluo]
  ff <- swap.names(ff) # this is from WadeTools
  gate <- rectangleGate("Ghost V500"=c(-Inf, bx(10^4))) # because I didn't gate live/dead on Omiq
  ff_gated <- Subset(ff, gate)
  
  # mat <- asinh(exprs(ff)/5) # change this to e.g. 500 for flow data
  mat <- asinh(exprs(ff_gated)/500)
  
  names <- colnames(ff_gated)
  nice_names <- names %>% 
    str_replace("Granzyme B", "Granzyme-B") %>% # change granzyme B name so it doesn't split on the space there
    str_split(" ") %>% # pass even column name and split on every space
    sapply("[", 1) # apply a function to the list and grab 1st item in array for each list item
  colnames(mat) <- nice_names
  
  sel_cells <- sample(nrow(mat), min(nrow(mat), nmax))
  
  return(mat[sel_cells,setdiff(colnames(mat),drop_channels)])
})

data <- do.call(rbind, data_list)
n_each <- sapply(data_list, nrow)
samples <- rep(files, n_each)

save(data, samples, file=paste0(dir_data, "/data.rda"))