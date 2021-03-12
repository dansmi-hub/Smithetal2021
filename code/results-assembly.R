# Consolodating the results into flux values and a data frame

library(tidyverse)

# Read in the data --------------------------------------------------------

# Webs with parasites
bp <- read_csv("results/bsq_parasites_flux.csv")
cp <- read_csv("results/csm_parasites_flux.csv")
ep <- read_csv("results/epb_parasites_flux.csv")

# Webs with no parasites
bs <- read_csv("results/bsq_no_para_flux.csv")
cs <- read_csv("results/csm_no_para_flux.csv")
es <- read_csv("results/epb_no_para_flux.csv")

# Bind
pf <-  bind_rows(bp, cp, ep) %>% mutate(simulated = "Has Parasites")
sf <- bind_rows(bs, cs, es) %>% mutate(simulated = "Replaced Parasites")

# Clean
data <- bind_rows(pf, sf) %>% 
  # transform vars
  mutate(lflux = log(flux),
         zlflux = scale(flux)) %>% 
  # clean up interaction strings
  mutate(interaction = str_remove(interaction, pattern = "_bsq"),
         interaction = str_remove(interaction, pattern = "_csm"),
         interaction = str_remove(interaction, pattern = "_epb"),
         interaction = str_remove(interaction, pattern = "sim"))

# Recode factor levels for interaction data
data <- data %>% mutate(interaction = as.factor(interaction))

# make totalflux the reference factor
data$interaction <- fct_relevel(data$interaction, "total", "herb", "detr", "carn", "para")

# two instances where detr had no flux - hinders lognormal modelling
data <- data %>% filter(flux != 0)

# Fri Mar 12 13:23:36 2021 ------------------------------
# Write data to a nice rdat file

write_rds(data, file = "results/fluxresults.rds")





