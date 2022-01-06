library(tidyverse)

og <- bind_rows(
  bsq <- read_csv("results/stats_webs_bsq.csv"),
  csm <- read_csv("results/stats_webs_csm.csv"),
  epb <- read_csv("results/stats_webs_epb.csv")
)

gen <- bind_rows(
  bsq <- read_csv("results/stats_nichemodelwebs_bsq.csv"),
  csm <- read_csv("results/stats_nichemodelwebs_csm.csv"),
  epb <- read_csv("results/stats_nichemodelwebs_epb.csv")
)


og$Web <- c("BSQ", "CSM", "EPB")
og$Parasites <- "Yes"

gen$Web <- c("BSQ", "CSM", "EPB")
gen$Parasites <- "No"

bind_rows(og, gen) %>% write_csv("tables/Web_Statistics.csv")

