# Author - Daniel Smith
# Email - dansmi@ceh.ac.uk

# Script to simulate energy flux of food webs with and without
# parasites in three North American estuarine webs

# This code has been re-written to better effectively manage memory
# Usage on UNIX only reccomended mutlicore CPU which can be modified 
# by using the mc.cores function in mclapply core = 12 and 100GB RAM 
# more cores = more RAM

rm(list = ls())

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(igraph)
library(NetIndices)
library(mice)
library(janitor)
library(fluxweb)
library(parallel)

# Start -------------------------------------------------------------------

# To edit how many cores the process takes over edit the 
# corecount variable at the start of each script e.g:
# corecount = 8 | 16 | 4 | etc...

# This is currently set in each individual script 
# would be nice to have this be set in this base script:
# but clearing memory between each script takes precedence

## Below no longer accurate after some optomisations and reducing simulation
## intensity as had no effect on results (e.g webs simulated 5000 -> 500).

# corecount = 2, uses 16GB RAM
# corecount = 4, uses 32GB RAM

# Going to split the analysis into three portions divided by foodweb, 
# allowing efficient clearing of memory between food webs

set.seed(1234)

# Create File Path if doesn't exist

filepath = file.path("rdata")

ifelse(!dir.exists(filepath), dir.create(filepath), FALSE)

print("Sourcing script for bsq...")
source("code/bsq.R", echo = F)
save.image(file='rdata/bsq.RData')
print("Finished BSQ!")
rm(list = ls())
gc()


print("Sourcing script for csm...")
source("code/csm.R", echo = F)
save.image(file='rdata/csm.RData')
print("Finished CSM!")
rm(list = ls())
gc()


print("Sourcing script for epb...")
source("code/epb.R", echo = F)
save.image(file='rdata/epb.RData')
print("Finished EPB!")
rm(list = ls())
gc()

## Put all results together for plotting / Analysis
source("code/results-assembly.R")

