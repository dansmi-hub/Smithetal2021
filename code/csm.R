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

# Load CSM ----------------------------------------------------------------

nodes <- read.csv("data/CSMweb_Nodes.csv") %>% 
  clean_names() %>% as_tibble()
links <- read.csv("data/CSMweb_Links.csv") %>% 
  clean_names() %>%  as_tibble()

# Who are pathogens? - need to remove them not required or point of analysis
takemefromweb <- nodes %>% filter(consumer_strategy_stage != "pathogen") %>% pull(species_id_stage_id)

nodes <- nodes %>% filter(consumer_strategy_stage != "pathogen")

links <- links %>% 
  filter(consumer_species_id_stage_id %in% takemefromweb) %>% 
  filter(resource_species_id_stage_id %in% takemefromweb)

# Also two species of parasitic plants should be removed - likely function differnt to animal parasites
# "Salt Marsh Bird's Beak" and "Dodder" - This is unique to CSM
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

write_csv(impute_me, "results/csm_original.csv")

# Predictor Matrix --------------------------------------------------------

pred <- make.predictorMatrix(impute_me)

pred[c("logM", "logB"), "logN"] <-  0

meth <- make.method(impute_me)

meth["logN"] <- "~I(logB - logM)"

# Impute ------------------------------------------------------------------

imputed_csm <- parlmice(impute_me, maxit = 50, printFlag = FALSE, meth = meth, predictorMatrix = pred,
                        n.core = corecount, n.imp.core = 25)

plot1 <- plot(imputed_csm)

png("plots/csm_convergence.png", width = 7, height = 5, units = 'in', res = 320)
plot1 # Make plot
dev.off()

complete_csm <- 
  complete(imputed_csm, "long") %>% 
  group_split(.imp)

# write this to a .csv
bind_rows(complete_csm) %>% write_csv("results/csm_imputed.csv")

# Easy dplyr::select -------------------------------------------------------------

df_nameswant <- c("species_id_stage_id",
                  "M",
                  "N",
                  "node_type",
                  "working_name",
                  "organismal_group",
                  "consumer_strategy_stage")

# Clean Imputed -----------------------------------------------------------

flux_csm <-
  complete_csm %>% 
  map(dplyr::select, logM, logN) %>% 
  map(cbind, nodes) %>% 
  map(mutate,
      # name the mutate variables
      M = exp(logM), 
      N = exp(logN)) %>%
  map(dplyr::select, df_nameswant) %>% 
  map(as_tibble)

# Recode Organism Type Factors --------------------------------------------

org_type_csm <- 
  flux_csm[[1]] %>% 
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
    animal = "spider",
    animal = "water boatman",
    animal = "burrowing shrimp",
    animal = "crab",
    animal = "fish",
    animal = "elasmobranch",
    animal = "bird",
    animal = "mammal",
    animal = "myxozoan",
    animal = "monogenean",
    animal = "cestode",
    para = "nematode",
    animal = "acanthocephalan",
    # animal = "anthozoan",
    # animal = "holothurian",
    # animal = "phoronid",
    # animal = "turbellarian",
    # viruslabel as detritus for now,
    detritus = "virus",
    para = "trematode"
  )

# Recode Consumer Factors -------------------------------------------------

con_type_csm <- 
  flux_csm[[1]] %>% 
  pull(consumer_strategy_stage) %>% 
  as.factor() %>%
  fct_recode(
    plant = "autotroph",
    detritus = "detritus",
    detritus = "pathogen",
    animal = "predator",
    animal = "detritivore",
    animal = "micropredator",
    # animal = "parasitoid",
    para = "macroparasite",
    para = "nonfeeding",
    para = "parasitic castrator",
    para = "trophically transmitted parasite"
  )

feed_type_csm <- 
  flux_csm[[1]] %>% 
  pull(consumer_strategy_stage) %>% 
  as.factor() %>%
  fct_recode(
    autotroph = "autotroph",
    detritus = "detritus",
    detritus = "pathogen",
    carnivore = "predator",
    detritivore = "detritivore",
    carnivore = "micropredator",
    # animal = "parasitoid",
    para = "macroparasite",
    para = "nonfeeding",
    para = "parasitic castrator",
    para = "trophically transmitted parasite"
  )

# Construct Matrix ---------------------------------------------------------

mat_csm <- 
  links %>% 
  mutate(resource_id = as.character(resource_species_id_stage_id),
         consumer_id = as.character(consumer_species_id_stage_id)) %>% 
  dplyr::select(resource_id, consumer_id) %>% 
  graph_from_data_frame(directed = T, vertices = flux_csm[[1]]$species_id_stage_id) %>% 
  as_adj() %>% 
  as.matrix()

# Vectors -----------------------------------------------------------------

M_csm <- 
  flux_csm %>% 
  lapply(pull, M)

N_csm <- 
  flux_csm %>% 
  lapply(pull, N)

B_csm <- 
  flux_csm %>% 
  map(mutate,
      B = M * N) %>% 
  map(pull, B)

# Define Metabolic Types --------------------------------------------------

met_types <-  c("ecto_vert", "endo_vert", "invert")


# Define Losses and Efficiencies  -----------------------------------------

losses_csm <- 
  lapply(M_csm, function(x){
    losses = rep(NA, length(x))
    ecto.vert = met_types == "ecto_vert"
    endo.vert = met_types == "endo_vert"
    inv = met_types == "invert"
    losses[ecto.vert] = 18.18 * x[ecto.vert] ^ (-0.29)
    losses[endo.vert] = 19.5 * x[endo.vert] ^ (-0.29)
    losses[inv] = 18.18 * x[inv] ^ (-0.29)
    losses
  })

efficiencies_csm <- 
  lapply(M_csm, function(x){
    efficiencies_csm = rep(NA, length(x[[1]]))
    efficiencies_csm[con_type_csm == "animal"] = 0.906
    efficiencies_csm[con_type_csm == "para"] = 0.906
    efficiencies_csm[con_type_csm == "plant"] = 0.545 
    efficiencies_csm[con_type_csm == "detritus"] = 0.158
    efficiencies_csm
  })


# Fluxweb -----------------------------------------------------------------

csm_fluxes <- 
  lapply(1:length(N_csm), function(x){
    try({fluxing(
      mat = mat_csm,
      biomasses = B_csm[[x]],
      losses = losses_csm[[x]],
      efficiencies = efficiencies_csm[[x]]
    )}, silent = TRUE)
  })

# which simulations worked?
csm_fluxes <- Filter(is.numeric, csm_fluxes)

# pull the ones that did work
# csm_fluxes <- csm_fluxes[1:100]


# Graphs ------------------------------------------------------------------

# network structure:
graph <- 
  links %>%
  mutate(resource = as.character(resource_species_id_stage_id),
         consumer = as.character(consumer_species_id_stage_id)) %>% 
  dplyr::select(resource, consumer) %>% 
  graph_from_data_frame(directed = T, vertices = NULL)

write_graph(graph, "results/csm_graph.txt", format = "edgelist")

# in matrix form
mat <-
  graph %>% 
  get.adjacency() %>% # dcgClass matrix
  as.matrix() # traditional matrix form

## --------------------------- Testing
var <- mat %>% rownames() %>% as.numeric()
if (any(var %in% noplantparasites_ids)) {
  stop("Rownames hava paraplants")
} 

var <- mat %>% colnames() %>% as.numeric()
if (any(var %in% noplantparasites_ids)) {
  stop("Colnames have para plants")
} 
## --------------------------- Testing

TrophInd <- NetIndices::TrophInd(mat)

## --------------------------- Testing TrophInd
var <- TrophInd %>% rownames %>% as.numeric()
if (any(var %in% noplantparasites_ids)) {
  stop("Para plants are in the sampling webs")
} 


## ---------------------------

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
data.frame(as.list(simple_stats)) %>% write_csv("results/stats_webs_csm.csv")

###########################################################################
# Counting --------------------------------------------------------------
# Now counting overall flux value in webs for CSM with parasites
############################################################################

# Initialise Vectors ------------------------------------------------------

herb_csm = rep(NA, length(csm_fluxes))
carn_csm = rep(NA, length(csm_fluxes))
detr_csm = rep(NA, length(csm_fluxes))
para_csm = rep(NA, length(csm_fluxes))
total_csm = rep(NA, length(csm_fluxes))

# Fill vectors ------------------------------------------------------------

for (x in 1:length(csm_fluxes)) {
  # herb_csm
  herb_csm[[x]] = sum(rowSums(csm_fluxes[[x]][org_type_csm == "plant", ]))
  # carn_csm
  carn_csm[[x]] = sum(rowSums(csm_fluxes[[x]][org_type_csm == "animal", ]))
  # detr_csm
  detr_csm[[x]] = sum(rowSums(csm_fluxes[[x]][org_type_csm == "detritus", ]))
  # para_csm
  para_csm[[x]] = sum(rowSums(csm_fluxes[[x]][org_type_csm == "para", ]))
  # total_csm
  total_csm[[x]] = sum(csm_fluxes[[x]])
}

# Gather Data -------------------------------------------------------------

csm_dat <- 
  cbind(herb_csm, carn_csm, detr_csm, para_csm, total_csm) %>% 
  as_tibble %>% 
  gather(key = interaction, value = flux) %>% 
  mutate(system = "CSM")

# Write Data --------------------------------------------------------------

write_csv(csm_dat, "results/csm_parasites_flux.csv")

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
niche_csm <- Niche.model(S = S, L = L, N = 500)
print("Done...")

# Get TL of Simulated Webs ------------------------------------------------

niche_csm_TL <- lapply(niche_csm, NetIndices::TrophInd)

TrophInd$species_id_stage_id <- as.numeric(rownames(TrophInd))

suppressMessages(
troph_csm <- 
  flux_csm %>% 
  map(left_join, TrophInd) %>% 
  map(cbind, 
      org_type = org_type_csm,
      con_type = con_type_csm) %>%
  bind_rows()
)

# get fluxes for each node
suppressMessages(
csm_data_para <- lapply(csm_fluxes, rowSums) %>% 
  # now take the names of those generated and the fluxes 
  # put into a dataframe to bind too
  lapply(function(x) {
    test <- data.frame(flux = x, species_id_stage_id = as.numeric(names(x)))
  }) %>% 
  lapply(left_join, troph_csm) %>% 
  bind_rows() %>% 
  distinct(flux, species_id_stage_id, .keep_all = TRUE)
)

write_csv(csm_data_para, "results/csm_raw.csv")
# Cut into Trophic Levels -------------------------------------------------

csm_TL <- 
  troph_csm %>% 
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
csm_TL[[1]] <- 
  csm_TL[[1]] %>% 
  enframe() %>% 
  subset(value != "para") %>% 
  subset(value != "animal") %>% 
  droplevels() %>% 
  pull(value)


csm_M <- 
  troph_csm %>% 
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

csm_N <- 
  troph_csm %>% 
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

csm_sim_TL <- 
  niche_csm_TL %>% 
  # add organism column, and cut TL's into bins of 0.5
  map(mutate,
      TL_cut = cut_width(TL, .5),
      org_type = NA,
      M = NA,
      N = NA)

# Sampling and matching ---------------------------------------------------

out_csm <- list()

for (i in 1:length(csm_sim_TL)) {
  # gets the list structure of the food web we want to sample into
  x <- csm_sim_TL[[i]]
  y <- group_split(x, TL_cut)
  z <- map(y, pull, org_type)
  # create M list of structure z
  out_csm$m[[i]] <- z
  # create N list of structure z
  out_csm$n[[i]] <- z
  # create org list of structure z
  out_csm$org[[i]] <- z
  # sample the M data
  for (j in 1:length(z)) {
    if (j < length(csm_N)) {
      out_csm$m[[i]][[j]] <- sample(csm_M[[j]], replace = TRUE, length(z[[j]]))
      # now find the index numbers of those sampled
      a <- match(x = out_csm$m[[i]][[j]], table = csm_M[[j]])
      # now take these index numbers and take their respective values for N
      out_csm$n[[i]][[j]] <- csm_N[[j]][a]
      # now give them the matching animal type
      out_csm$org[[i]][[j]] <- csm_TL[[j]][a]
    } else {
      out_csm$m[[i]][[j]] <- sample(csm_M[[length(csm_M)]], replace = TRUE, length(z[[j]]))
      # now find the index numbers of those sampled
      a <- match(x = out_csm$m[[i]][[j]], table = csm_M[[length(csm_M)]])
      # now take these index numbers and take their respective values for N
      out_csm$n[[i]][[j]] <- csm_N[[length(csm_M)]][a]
      # now give them the matching animal type
      out_csm$org[[i]][[j]] <- csm_TL[[length(csm_M)]][a]
    }
  }
}

# assign trophic levels
for (i in 1:length(csm_sim_TL)) {
  out_csm$TL[[i]] <- csm_sim_TL[[i]]$TL
}

# Assign vector for biomass
out_csm$b <- out_csm$m

# This calculates biomass from N and M for each organism in the web
for(i in seq(out_csm$b)){
  for(j in seq(out_csm$b[[i]])){
    out_csm$b[[i]][[j]] <- out_csm$m[[i]][[j]] * out_csm$n[[i]][[j]]
  }
}

# Testing!! ---------------------------------------------------------------

# Unit test to find plants and check simulation sampling
for (j in 1:500) {
  data <- c()
  for (i in seq_along(out_csm$org[[j]])) {
    test <- out_csm$org %>% map(i) %>% unlist
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
  for (i in seq_along(out_csm$org[[j]])) {
    test <- out_csm$org %>% map(i) %>% unlist
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

out_csm$org <- lapply(out_csm$org, unlist)
out_csm$m <- lapply(out_csm$m, unlist)
out_csm$n <- lapply(out_csm$n, unlist)
out_csm$b <- lapply(out_csm$b, unlist)


# Losses and Efficiencies -------------------------------------------------

losses_csm_sim <- 
  lapply(out_csm$m, function(x){
    losses = rep(NA, length(x))
    ecto.vert = met_types == "ecto_vert"
    endo.vert = met_types == "endo_vert"
    inv = met_types == "invert"
    losses[ecto.vert] = 18.18 * x[ecto.vert] ^ (-0.29)
    losses[endo.vert] = 19.5 * x[endo.vert] ^ (-0.29)
    losses[inv] = 18.18 * x[inv] ^ (-0.29)
    losses
  })

efficiencies_csm_sim <- 
  lapply(seq(out_csm$org), function(x) {
    efficiencies_csm_sim = rep(NA, length(out_csm$org[[x]]))
    efficiencies_csm_sim[out_csm$org[[x]] == "animal"] = 0.906
    efficiencies_csm_sim[out_csm$org[[x]] == "para"] = 0.906
    efficiencies_csm_sim[out_csm$org[[x]] == "plant"] = 0.545
    efficiencies_csm_sim[out_csm$org[[x]] == "detritus"] = 0.906
    efficiencies_csm_sim
  })


# Flux Calc ---------------------------------------------------------------

csm_fluxes_sim <- 
  lapply(seq(niche_csm), function(x){
    try(fluxing(
      mat = niche_csm[[x]],
      biomasses = out_csm$b[[x]],
      losses = losses_csm_sim[[x]],
      efficiencies = efficiencies_csm_sim[[x]]
    ), silent = TRUE)
  })

winners_csm <- Filter(is.numeric, csm_fluxes_sim)
length(winners_csm) # how many without negative flux value

# Which index values of list worked/failed? -------------------------------

winners_csm_index <- # these ones won
  unlist(lapply(csm_fluxes_sim, function(x){
    is.numeric(x)
  }))

# filter like this
out_csm$m[winners_csm_index]

# Filter index winners ----------------------------------------------------

# need to filter these so we aggregate the right fluxes for each animal type
winners_out_csm <- list()
winners_out_csm$m <- out_csm$m[winners_csm_index]
winners_out_csm$n <- out_csm$n[winners_csm_index]
winners_out_csm$b <- out_csm$b[winners_csm_index]
winners_out_csm$org <- out_csm$org[winners_csm_index]
# Add TL
winners_out_csm$TL <- out_csm$TL[winners_csm_index]
# Add FLuxes for each row?
winners_out_csm$flux <- lapply(csm_fluxes_sim[winners_csm_index], rowSums)

# get fluxes for each node in the simulated webs
csm_raw2 <- transpose(winners_out_csm) %>% 
  lapply(bind_rows) %>% 
  bind_rows(.id = "Iteration") %>% 
  distinct()

write_csv(csm_raw2, "results/csm_raw2.csv")


# Aggregate fluxes --------------------------------------------------------

herb_csmsim = rep(NA, length(winners_csm))
carn_csmsim = rep(NA, length(winners_csm))
detr_csmsim = rep(NA, length(winners_csm))
para_csmsim = rep(NA, length(winners_csm))
total_csmsim = rep(NA, length(winners_csm))

for (i in seq_along(winners_csm)) {
  # print(i)
  tryCatch({
    # herbivory
    herb_csmsim[[i]] = sum(rowSums(winners_csm[[i]][winners_out_csm$org[[i]] == 'plant', ]))
    # carnivory
    carn_csmsim[[i]] = sum(rowSums(winners_csm[[i]][winners_out_csm$org[[i]] == 'animal', ]))
    # para
    para_csmsim[[i]] = sum(rowSums(winners_csm[[i]][winners_out_csm$org[[i]] == 'para', ]))
    # detritivory
    detr_csmsim[[i]] = sum(rowSums(winners_csm[[i]][winners_out_csm$org[[i]] == 'detritus', ]))
    # total 
    total_csmsim[[i]] = sum(winners_csm[[i]])
  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

# errors where detritus is only a vector and not array because only one detritus in TL


# Get web stats for those niche model webs that won and we used -----------

niche_stats <- lapply(niche_csm[winners_csm_index], Stats_matrix)
bind_rows(niche_stats) %>% sapply(mean) %>% enframe() %>% pivot_wider() %>% 
  write_csv("results/stats_nichemodelwebs_csm.csv")

# Write niche mdoel graphs ------------------------------------------------

niche_webs <- niche_csm[winners_csm_index] %>% lapply(graph_from_adjacency_matrix)

n <- 1:length(niche_webs)
names(niche_webs) <- sprintf("csm_nicheweb_%03d", n)

# Write these to a directory

for(i in names(niche_webs)){
  write.graph(niche_webs[[i]], file.path("graphs", "csm", paste0(i,".txt")))
}

# Gather and Write -----------------------------------

sim_csm_dat <- 
  cbind(herb_csmsim, carn_csmsim, detr_csmsim, total_csmsim, para_csmsim) %>% 
  as_tibble() %>% 
  gather(key = interaction, value = flux) %>% 
  mutate(system = "CSM")

write_csv(sim_csm_dat, "results/csm_no_para_flux.csv")

save.image("rdata/csm.RData")

print("Done!")
