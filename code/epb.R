# If you want to run just this web you'll have to set params here:

rm(list = ls())

corecount = 4

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(igraph)
library(NetIndices)
library(mice)
library(janitor)
library(fluxweb)
library(parallel)

# Load EPB ----------------------------------------------------------------

nodes <- read.csv("data/EPBweb_Nodes.csv") %>% 
  clean_names() %>% as_tibble()
links <- read.csv("data/EPBweb_Links.csv") %>% 
  clean_names() %>%  as_tibble()

# Who are pathogens? - need to remove them not required or point of analysis
takemefromweb <- nodes %>% filter(consumer_strategy_stage != "pathogen") %>% pull(species_id_stage_id)

nodes <- nodes %>% filter(consumer_strategy_stage != "pathogen")

links <- links %>% 
  filter(consumer_species_id_stage_id %in% takemefromweb) %>% 
  filter(resource_species_id_stage_id %in% takemefromweb)

# Also two species of parasitic plants should be removed - likely function differnt to animal parasites
# "Salt Marsh Bird's Beak" and "Dodder" - This is unique to epb
noplantparasites_ids <- nodes %>% 
  filter(working_name == "Salt Marsh Bird's Beak" | working_name == "Dodder") %>% 
  pull(species_id_stage_id)

`%nin%` <- negate(`%in%`)

nodes <- nodes %>% filter(species_id_stage_id %nin% noplantparasites_ids)

links <- links %>%
  filter(consumer_species_id_stage_id %nin% noplantparasites_ids) %>% 
  filter(resource_species_id_stage_id %nin% noplantparasites_ids)

## --------------------------- Test If Para Plants Removed

if (nodes$working_name %>% has_element("Salt Marsh Bird's Beak") == T) {
  stop("Parasite plants not taken out of nodes")
}

if (nodes$working_name %>% has_element("Dodder") == T) {
  stop("Parasite plants not taken out of nodes")
}

## --------------------------- Test links for same thing

var <- nodes %>% filter(working_name == "Salt Marsh Bird's Beak" | working_name == "Dodder") %>% pull(species_id_stage_id)

if (any(links$consumer_species_id_stage_id %in% var)) {
  stop("Para plants are in the links")
} 


# Clean -------------------------------------------------------------------

foo <- 
  nodes %>%
  dplyr::select(M = body_size_g,
                N = abundance_no_ha,
                B = biomass_kg_ha,
                organismal_group,
                consumer_strategy_stage,
                species_id_stage_id,
                node_id) %>% 
  mutate(M = M / 1000)

impute_me <- 
  foo %>% 
  mutate(logM = log(M),
         logN = case_when(N > 0 ~ log(N)),
         logB = log(B)) %>% 
  # phylum = na_if(phylum, ""),
  # class = na_if(class, ""),
  # order = na_if(order, "")) %>% 
  dplyr::select(logM, logN, logB)

write_csv(impute_me, "results/epb_original.csv")

# Predictor Matrix --------------------------------------------------------

pred <- make.predictorMatrix(impute_me)

pred[c("logM", "logB"), "logN"] <-  0

meth <- make.method(impute_me)

meth["logN"] <- "~I(logB - logM)"

# Impute ------------------------------------------------------------------

imputed_epb <- parlmice(impute_me, maxit = 50, printFlag = FALSE, meth = meth, predictorMatrix = pred,
                        n.core = corecount, n.imp.core = 25)

plot1 <- plot(imputed_epb)

png("plots/epb_convergence.png", width = 7, height = 5, units = 'in', res = 320)
plot1 # Make plot
dev.off()

complete_epb <- 
  complete(imputed_epb, "long") %>% 
  group_split(.imp)

# write this to a .csv
bind_rows(complete_epb) %>% write_csv("results/epb_imputed.csv")

# Easy dplyr::select -------------------------------------------------------------

df_nameswant <- c("species_id_stage_id",
                  "M",
                  "N",
                  "node_type",
                  "working_name",
                  "organismal_group",
                  "consumer_strategy_stage")

# Clean Imputed -----------------------------------------------------------

flux_epb <-
  complete_epb %>% 
  map(dplyr::select, logM, logN) %>% 
  map(cbind, nodes) %>% 
  map(mutate,
      # name the mutate variables
      M = exp(logM), 
      N = exp(logN)) %>%
  map(dplyr::select, df_nameswant) %>% 
  map(as_tibble)

# Recode Organism Type Factors --------------------------------------------

org_type_epb <- 
  flux_epb[[1]] %>% 
  pull(organismal_group) %>% 
  as.factor() %>%
  fct_recode(
    plant = "vascular plant",
    plant = "microphytobenthos",
    plant = "macroalgae",
    # animals
    animal = "protist",
    animal = "annelid",
    animal = "leech",
    animal = "nemertean",
    animal = "bivalve",
    animal = "snail",
    animal = "mosquito",
    animal = "branchiuran",
    animal = "amphipod",
    animal = "copepod",
    animal = "dipteran",
    animal = "isopod",
    animal = "ostracod",
    # animal = "spider",
    animal = "water boatman",
    animal = "burrowing shrimp",
    animal = "crab",
    animal = "fish",
    animal = "elasmobranch",
    animal = "bird",
    # animal = "mammal",
    # animal = "myxozoan",
    animal = "monogenean",
    animal = "cestode",
    para = "nematode",
    animal = "acanthocephalan",
    animal = "anthozoan",
    animal = "holothurian",
    animal = "phoronid",
    animal = "turbellarian",
    # viruslabel as detritus for now,
    para = "virus",
    para = "trematode"
  )

# Recode Consumer Factors -------------------------------------------------

con_type_epb <- 
  flux_epb[[1]] %>% 
  pull(consumer_strategy_stage) %>% 
  as.factor() %>%
  fct_recode(
    plant = "autotroph",
    detritus = "detritus",
    detritus = "pathogen",
    animal = "predator",
    animal = "detritivore",
    animal = "micropredator",
    animal = "parasitoid",
    para = "macroparasite",
    para = "nonfeeding",
    para = "parasitic castrator",
    para = "trophically transmitted parasite"
  )

feed_type_epb <- 
  flux_epb[[1]] %>% 
  pull(consumer_strategy_stage) %>% 
  as.factor() %>%
  fct_recode(
    autotroph = "autotroph",
    detritus = "detritus",
    detritus = "pathogen",
    carnivore = "predator",
    detritivore = "detritivore",
    carnivore = "micropredator",
    # carnivore = "parasitoid",
    para = "macroparasite",
    para = "nonfeeding",
    para = "parasitic castrator",
    para = "trophically transmitted parasite"
  )

# Construct Matrix ---------------------------------------------------------

mat_epb <- 
  links %>% 
  mutate(resource_id = as.character(resource_species_id_stage_id),
         consumer_id = as.character(consumer_species_id_stage_id)) %>% 
  dplyr::select(resource_id, consumer_id) %>% 
  graph_from_data_frame(directed = T, vertices = flux_epb[[1]]$species_id_stage_id) %>% 
  as_adj() %>% 
  as.matrix()

# Vectors -----------------------------------------------------------------

M_epb <- 
  flux_epb %>% 
  lapply(pull, M)

N_epb <- 
  flux_epb %>% 
  lapply(pull, N)

B_epb <- 
  flux_epb %>% 
  map(mutate,
      B = M * N) %>% 
  map(pull, B)

# Define Metabolic Types --------------------------------------------------

met_types <-  c("ecto_vert", "endo_vert", "invert")


# Define Losses and Efficiencies  -----------------------------------------

losses_epb <- 
  lapply(M_epb, function(x){
    losses = rep(NA, length(x))
    ecto.vert = met_types == "ecto_vert"
    endo.vert = met_types == "endo_vert"
    inv = met_types == "invert"
    losses[ecto.vert] = 18.18 * x[ecto.vert] ^ (-0.29)
    losses[endo.vert] = 19.5 * x[endo.vert] ^ (-0.29)
    losses[inv] = 18.18 * x[inv] ^ (-0.29)
    losses
  })

efficiencies_epb <- 
  lapply(M_epb, function(x){
    efficiencies_epb = rep(NA, length(x[[1]]))
    efficiencies_epb[con_type_epb == "animal"] = 0.906
    efficiencies_epb[con_type_epb == "para"] = 0.906
    efficiencies_epb[con_type_epb == "plant"] = 0.545 
    efficiencies_epb[con_type_epb == "detritus"] = 0.158
    efficiencies_epb
  })


# Fluxweb -----------------------------------------------------------------

epb_fluxes <- 
  lapply(1:length(N_epb), function(x){
    try({fluxing(
      mat = mat_epb,
      biomasses = B_epb[[x]],
      losses = losses_epb[[x]],
      efficiencies = efficiencies_epb[[x]]
    )}, silent = TRUE)
  })

# which simulations worked?
epb_fluxes <- Filter(is.numeric, epb_fluxes)

# pull the ones that did work
# epb_fluxes <- epb_fluxes[1:100]


# Graphs ------------------------------------------------------------------

# network structure:
graph <- 
  links %>%
  mutate(resource = as.character(resource_species_id_stage_id),
         consumer = as.character(consumer_species_id_stage_id)) %>% 
  dplyr::select(resource, consumer) %>% 
  graph_from_data_frame(directed = T, vertices = NULL)

write_graph(graph, "results/epb_graph.txt", format = "edgelist")

# in matrix form
mat <-
  graph %>% 
  get.adjacency() %>% # dcgClass matrix
  as.matrix() # traditional matrix form

TrophInd <-  NetIndices::TrophInd(mat)

# network stats

Stats_matrix <- function(mat){
  S = nrow(mat)
  L = sum(mat)
  
  basal = sum(colSums(mat) == 0)/S
  top   = sum(colSums(t(mat)) == 0)/S
  int   = 1 - basal - top
  gen   = mean(colSums(mat))
  vun   = mean(rowSums(mat))
  gensd = sd(colSums(mat)/(L/S))
  vunsd = sd(rowSums(mat)/(L/S))
  
  return(c("Nodes" = S,
           "Links" = L,
           "Basal" = basal,
           "Intermediate" = int,
           "Top" = top,
           "Generality" = gen,
           "Generality_sd" = gensd,
           "Vulnerability" = vun,
           "Vulnerabilty_sd" = vunsd))
}

simple_stats <-  Stats_matrix(mat)
data.frame(as.list(simple_stats)) %>% write_csv("results/stats_webs_epb.csv")

###########################################################################
# Counting --------------------------------------------------------------
# Now counting overall flux value in webs for EPB with parasites
############################################################################

# Initialise Vectors ------------------------------------------------------

herb_epb = rep(NA, length(epb_fluxes))
carn_epb = rep(NA, length(epb_fluxes))
detr_epb = rep(NA, length(epb_fluxes))
para_epb = rep(NA, length(epb_fluxes))
total_epb = rep(NA, length(epb_fluxes))

# Fill vectors ------------------------------------------------------------

for (x in 1:length(epb_fluxes)) {
  # herb_epb
  herb_epb[[x]] = sum(rowSums(epb_fluxes[[x]][org_type_epb == "plant", ]))
  # carn_epb
  carn_epb[[x]] = sum(rowSums(epb_fluxes[[x]][org_type_epb == "animal", ]))
  # detr_epb
  detr_epb[[x]] = sum(rowSums(epb_fluxes[[x]][org_type_epb == "detritus", ]))
  # para_epb
  para_epb[[x]] = sum(rowSums(epb_fluxes[[x]][org_type_epb == "para", ]))
  # total_epb
  total_epb[[x]] = sum(epb_fluxes[[x]])
}

# Gather Data -------------------------------------------------------------

epb_dat <- 
  cbind(herb_epb, carn_epb, detr_epb, para_epb, total_epb) %>% 
  as_tibble %>% 
  gather(key = interaction, value = flux) %>% 
  mutate(system = "EPB")

# Write Data --------------------------------------------------------------

write_csv(epb_dat, "results/epb_parasites_flux.csv")

###########################################################################
# Moving to Simulating Webs  ----------------------------------------------
# Now simuulating webs without parasites randomly over many iterations
###########################################################################

# Niche Model -------------------------------------------------------------

Niche.model <- function(S, L, N = 1){
  C <- L/S^2
  if(N==1){
    n <- sort(runif(S))
    beta <- (1 - 2 * C) / (2 * C)
    r <- n*(1 - (1 - runif(S))^(1/beta))
    c <- r/2 + runif(S) * (n - r/2)
    web <- matrix(0,S,S)
    min.n <- c-r/2
    max.n <- c+r/2
    for(i in 1:S){
      diet <- c(1:S)[c(which(n>min.n[i]), which(n<max.n[i]))[duplicated(c(which(n>min.n[i]), which(n<max.n[i])))]]
      web[diet,i] <- 1
    }
    dimnames(web) <- list(1:length(web[,1]), 1:length(web[,1]))
  }
  if(N>1){
    web <- list()
    for(j in 1:N){
      n <- sort(runif(S))
      beta <- (1 - 2 * C) / (2 * C)
      r <- n*(1 - (1 - runif(S))^(1/beta))
      c <- r/2 + runif(S) * (n - r/2)
      web[[j]] <- matrix(0,S,S)
      min.n <- c-r/2
      max.n <- c+r/2
      for(i in 1:S){
        diet <- c(1:S)[c(which(n>min.n[i]), which(n<max.n[i]))[duplicated(c(which(n>min.n[i]), which(n<max.n[i])))]]
        web[[j]][diet,i] <- 1
      }
      dimnames(web[[j]]) <- list(1:length(web[[j]][,1]), 1:length(web[[j]][,1]))
    }
  }
  web
}

# Generate Webs with Prior Params -----------------------------------------

S = as.numeric(simple_stats["Nodes"])
L = as.numeric(simple_stats["Links"])

print("Simulating 500 Webs")
niche_epb <- Niche.model(S = S, L = L, N = 500)
print("Done...")

# Get TL of Simulated Webs ------------------------------------------------

niche_epb_TL <- lapply(niche_epb, NetIndices::TrophInd)

TrophInd$species_id_stage_id <- as.numeric(rownames(TrophInd))

suppressMessages(
troph_epb <- 
  flux_epb %>% 
  map(left_join, TrophInd) %>% 
  map(cbind, 
      org_type = org_type_epb,
      con_type = con_type_epb) %>%
  bind_rows()
)

# get fluxes for each node
suppressMessages(
epb_data_para <- lapply(epb_fluxes, rowSums) %>% 
  # now take the names of those generated and the fluxes 
  # put into a dataframe to bind too
  lapply(function(x) {
    test <- data.frame(flux = x, species_id_stage_id = as.numeric(names(x)))
  }) %>% 
  lapply(left_join, troph_epb) %>% 
  bind_rows() %>% 
  distinct(flux, species_id_stage_id, .keep_all = TRUE)
)

write_csv(epb_data_para, "results/epb_raw.csv")

# Cut into Trophic Levels -------------------------------------------------

epb_TL <- 
  troph_epb %>% 
  # distinct
  distinct() %>%  
  # select useful cols
  dplyr::select(org_type, TL) %>% 
  # drop NA TL value from isolated node in web
  drop_na() %>% 
  # bin TL into .5 trophic level intervals
  mutate( 
    # arguments for mutate
    TL_cut = cut_width(TL, .5)) %>%
  # split into trophic cut groups
  group_split(TL_cut) %>% 
  # pull vector of org_type for each group
  lapply(pull, org_type)

# need to drop animals and non-feeding stage parasites from the first trophic level
epb_TL[[1]] <- 
  epb_TL[[1]] %>% 
  enframe() %>% 
  subset(value != "para") %>% 
  subset(value != "animal") %>% 
  droplevels() %>% 
  pull(value)


epb_M <- 
  troph_epb %>% 
  # distinct
  distinct() %>%  
  # remove parasite abundances from Null Model
  filter(org_type != "para") %>%  
  # select useful cols
  dplyr::select(M, TL) %>% 
  # drop NA TL value from isolated node in web
  drop_na() %>% 
  # bin TL into .5 trophic level intervals
  mutate( 
    # arguments for mutate
    TL_cut = cut_width(TL, .5)) %>%
  # split into trophic cut groups
  group_split(TL_cut) %>% 
  # pull vector of org_type for each group
  lapply(pull, M)

epb_N <- 
  troph_epb %>% 
  # distinct
  distinct() %>%  
  # remove parasite abundances from Null Model
  filter(org_type != "para") %>%  
  # select useful cols
  dplyr::select(N, TL) %>% 
  # drop NA TL value from isolated node in web
  drop_na() %>% 
  # bin TL into .5 trophic level intervals
  mutate( 
    # arguments for mutate
    TL_cut = cut_width(TL, .5)) %>%
  # split into trophic cut groups
  group_split(TL_cut) %>% 
  # pull vector of org_type for each group
  lapply(pull, N)

epb_sim_TL <- 
  niche_epb_TL %>% 
  # add organism column, and cut TL's into bins of 0.5
  map(mutate,
      TL_cut = cut_width(TL, .5),
      org_type = NA,
      M = NA,
      N = NA)

# Sampling and matching ---------------------------------------------------

out_epb <- list()

for (i in 1:length(epb_sim_TL)) {
  # gets the list structure of the food web we want to sample into
  x <- epb_sim_TL[[i]]
  y <- group_split(x, TL_cut)
  z <- map(y, pull, org_type)
  # create M list of structure z
  out_epb$m[[i]] <- z
  # create N list of structure z
  out_epb$n[[i]] <- z
  # create org list of structure z
  out_epb$org[[i]] <- z
  # sample the M data
  for (j in 1:length(z)) {
    if (j < length(epb_N)) {
      out_epb$m[[i]][[j]] <- sample(epb_M[[j]], replace = TRUE, length(z[[j]]))
      # now find the index numbers of those sampled
      a <- match(x = out_epb$m[[i]][[j]], table = epb_M[[j]])
      # now take these index numbers and take their respective values for N
      out_epb$n[[i]][[j]] <- epb_N[[j]][a]
      # now give them the matching animal type
      out_epb$org[[i]][[j]] <- epb_TL[[j]][a]
    } else {
      out_epb$m[[i]][[j]] <- sample(epb_M[[length(epb_M)]], replace = TRUE, length(z[[j]]))
      # now find the index numbers of those sampled
      a <- match(x = out_epb$m[[i]][[j]], table = epb_M[[length(epb_M)]])
      # now take these index numbers and take their respective values for N
      out_epb$n[[i]][[j]] <- epb_N[[length(epb_M)]][a]
      # now give them the matching animal type
      out_epb$org[[i]][[j]] <- epb_TL[[length(epb_M)]][a]
    }
  }
}

# assign trophic levels
for (i in 1:length(epb_sim_TL)) {
  out_epb$TL[[i]] <- epb_sim_TL[[i]]$TL
}

# Assign vector for biomass
out_epb$b <- out_epb$m

# This calculates biomass from N and M for each organism in the web
for(i in seq(out_epb$b)){
  for(j in seq(out_epb$b[[i]])){
    out_epb$b[[i]][[j]] <- out_epb$m[[i]][[j]] * out_epb$n[[i]][[j]]
  }
}

# Testing!! ---------------------------------------------------------------

# Unit test to find plants and check simulation sampling
for (j in 1:500) {
  data <- c()
  for (i in seq_along(out_epb$org[[j]])) {
    test <- out_epb$org %>% map(i) %>% unlist
    data[[i]] <- "plant" %in% test
    data <- map_lgl(data, isTRUE)
    # remove first trophic level
    data <- data[-1]
  }
  # predicate to stop
  if (any(isTRUE(data))) {
    stop(paste0("Simulated webs have plants in multiple trophic levels: ", " \n Web with issue is ", i, "/500"))
  }
}

# Unit test to find detritus and check simulation sampling
for (j in 1:500) {
  data <- c()
  for (i in seq_along(out_epb$org[[j]])) {
    test <- out_epb$org %>% map(i) %>% unlist
    data[[i]] <- "detritus" %in% test
    data <- map_lgl(data, isTRUE)
    # remove first trophic level
    data <- data[-1]
  }
  # predicate to stop
  if (any(isTRUE(data))) {
    stop(paste0("Simulated webs have detritus in multiple trophic levels: ", " \n Web with issue is ", i, "/500"))
  }
}

# Unlist ------------------------------------------------------------------

out_epb$org <- lapply(out_epb$org, unlist)
out_epb$m <- lapply(out_epb$m, unlist)
out_epb$n <- lapply(out_epb$n, unlist)
out_epb$b <- lapply(out_epb$b, unlist)


# Losses and Efficiencies -------------------------------------------------

losses_epb_sim <- 
  lapply(out_epb$m, function(x){
    losses = rep(NA, length(x))
    ecto.vert = met_types == "ecto_vert"
    endo.vert = met_types == "endo_vert"
    inv = met_types == "invert"
    losses[ecto.vert] = 18.18 * x[ecto.vert] ^ (-0.29)
    losses[endo.vert] = 19.5 * x[endo.vert] ^ (-0.29)
    losses[inv] = 18.18 * x[inv] ^ (-0.29)
    losses
  })

efficiencies_epb_sim <- 
  lapply(seq(out_epb$org), function(x) {
    efficiencies_epb_sim = rep(NA, length(out_epb$org[[x]]))
    efficiencies_epb_sim[out_epb$org[[x]] == "animal"] = 0.906
    efficiencies_epb_sim[out_epb$org[[x]] == "para"] = 0.906
    efficiencies_epb_sim[out_epb$org[[x]] == "plant"] = 0.545
    efficiencies_epb_sim[out_epb$org[[x]] == "detritus"] = 0.906
    efficiencies_epb_sim
  })


# Flux Calc ---------------------------------------------------------------

epb_fluxes_sim <- 
  lapply(seq(niche_epb), function(x){
    try(fluxing(
      mat = niche_epb[[x]],
      biomasses = out_epb$b[[x]],
      losses = losses_epb_sim[[x]],
      efficiencies = efficiencies_epb_sim[[x]]
    ), silent = TRUE)
  })

winners_epb <- Filter(is.numeric, epb_fluxes_sim)
length(winners_epb) # how many without negative flux value

# Which index values of list worked/failed? -------------------------------

winners_epb_index <- # these ones won
  unlist(lapply(epb_fluxes_sim, function(x){
    is.numeric(x)
  }))

# filter like this
out_epb$m[winners_epb_index]

# Filter index winners ----------------------------------------------------

# need to filter these so we aggregate the right fluxes for each animal type
winners_out_epb <- list()
winners_out_epb$m <- out_epb$m[winners_epb_index]
winners_out_epb$n <- out_epb$n[winners_epb_index]
winners_out_epb$b <- out_epb$b[winners_epb_index]
winners_out_epb$org <- out_epb$org[winners_epb_index]
# Add TL
winners_out_epb$TL <- out_epb$TL[winners_epb_index]
# Add FLuxes for each row?
winners_out_epb$flux <- lapply(epb_fluxes_sim[winners_epb_index], rowSums)

# get fluxes for each node in the simulated webs
epb_raw2 <- transpose(winners_out_epb) %>% 
  lapply(bind_rows) %>% 
  bind_rows(.id = "Iteration") %>% 
  distinct()

write_csv(epb_raw2, "results/epb_raw2.csv")


# Aggregate fluxes --------------------------------------------------------

herb_epbsim = rep(NA, length(winners_epb))
carn_epbsim = rep(NA, length(winners_epb))
detr_epbsim = rep(NA, length(winners_epb))
para_epbsim = rep(NA, length(winners_epb))
total_epbsim = rep(NA, length(winners_epb))

for (i in seq_along(winners_epb)) {
  # print(i)
  tryCatch({
    # herbivory
    herb_epbsim[[i]] = sum(rowSums(winners_epb[[i]][winners_out_epb$org[[i]] == 'plant', ]))
    # carnivory
    carn_epbsim[[i]] = sum(rowSums(winners_epb[[i]][winners_out_epb$org[[i]] == 'animal', ]))
    # para
    para_epbsim[[i]] = sum(rowSums(winners_epb[[i]][winners_out_epb$org[[i]] == 'para', ]))
    # detritivory
    detr_epbsim[[i]] = sum(rowSums(winners_epb[[i]][winners_out_epb$org[[i]] == 'detritus', ]))
    # total 
    total_epbsim[[i]] = sum(winners_epb[[i]])
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

# errors where detritus is only a vector and not array because only one detritus in TL
# Get web stats for those niche model webs that won and we used -----------

niche_stats <- lapply(niche_epb[winners_epb_index], Stats_matrix)
bind_rows(niche_stats) %>% sapply(mean) %>% enframe() %>% pivot_wider() %>% 
  write_csv("results/stats_nichemodelwebs_epb.csv")

# Write niche model graphs ------------------------------------------------

niche_webs <- niche_epb[winners_epb_index] %>% lapply(graph_from_adjacency_matrix)

n <- 1:length(niche_webs)
names(niche_webs) <- sprintf("epb_nicheweb_%03d", n)

# Write these to a directory

for(i in names(niche_webs)){
  write.graph(niche_webs[[i]], file.path("graphs", "epb", paste0(i,".txt")))
}


# Gather and Write -----------------------------------

sim_epb_dat <- 
  cbind(herb_epbsim, carn_epbsim, detr_epbsim, total_epbsim, para_epbsim) %>% 
  as_tibble() %>% 
  gather(key = interaction, value = flux) %>% 
  mutate(system = "EPB")

write_csv(sim_epb_dat, "results/epb_no_para_flux.csv")

save.image("rdata/epb.RData")

print("Done - check results!")
