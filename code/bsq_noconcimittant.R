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
library(pbmcapply)

# Load BSQ ----------------------------------------------------------------

nodes <- read.csv("data/BSQweb_Nodes__New.csv") %>% 
  clean_names() %>% as_tibble()
links <- read.csv("data/BSQweb_New_Links.csv") %>% 
  clean_names() %>%  as_tibble()

# Who are pathogens? - need to remove them not required or point of analysis
takemefromweb <- nodes %>% filter(consumer_strategy_stage != "pathogen") %>% pull(species_id_stage_id)

nodes <- nodes %>% filter(consumer_strategy_stage != "pathogen")

links <- links %>% 
  filter(consumer_species_id_stage_id %in% takemefromweb) %>% 
  filter(resource_species_id_stage_id %in% takemefromweb)
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

write_csv(impute_me, "results/bsq_original.csv")

# Predictor Matrix --------------------------------------------------------

pred <- make.predictorMatrix(impute_me)

pred[c("logM", "logB"), "logN"] <-  0

meth <- make.method(impute_me)

meth["logN"] <- "~I(logB - logM)"


# Impute ------------------------------------------------------------------

imputed_bsq <- parlmice(impute_me, maxit = 50, printFlag = FALSE, meth = meth, predictorMatrix = pred,
                        n.core = corecount, n.imp.core = 25)

plot1 <- plot(imputed_bsq)

png("plots/bsq_convergence.png", width = 7, height = 5, units = 'in', res = 320)
plot1 # Make plot
dev.off()


complete_bsq <- 
  complete(imputed_bsq, "long") %>% 
  group_split(.imp)

# write this to a .csv
bind_rows(complete_bsq) %>% write_csv("results/bsq_imputed.csv")

# Easy dplyr::select -------------------------------------------------------------

df_nameswant <- c("species_id_stage_id",
                  "M",
                  "N",
                  "node_type",
                  "working_name",
                  "organismal_group",
                  "consumer_strategy_stage")

# Clean Imputed -----------------------------------------------------------

flux_bsq <-
  complete_bsq %>% 
  map(dplyr::select, logM, logN) %>% 
  map(cbind, nodes) %>% 
  map(mutate,
      # name the mutate variables
      M = exp(logM), 
      N = exp(logN)) %>%
  map(dplyr::select, df_nameswant) %>% 
  map(as_tibble)

# Recode Organism Type Factors --------------------------------------------

org_type_bsq <- 
  flux_bsq[[1]] %>% 
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
    animal = "anthozoan",
    animal = "amphipod",
    animal = "copepod",
    animal = "dipteran",
    animal = "isopod",
    animal = "ostracod",
    animal = "cumacean",
    animal = "water boatman",
    animal = "burrowing shrimp",
    animal = "crab",
    animal = "fish",
    animal = "elasmobranch",
    animal = "bird",
    animal = "leptostracan",
    animal = "ophiuroid",
    # animal = "monogenean",
    animal = "cestode",
    para = "nematode",
    animal = "acanthocephalan",
    animal = "tanaidacean",
    animal = "beetle",
    # animal = "phoronid",
    # animal = "turbellarian",
    # viruslabel as detritus for now,
    # detritus = "virus",
    para = "trematode"
  )

# Recode Consumer Factors -------------------------------------------------

con_type_bsq <- 
  flux_bsq[[1]] %>% 
  pull(consumer_strategy_stage) %>% 
  as.factor() %>%
  fct_recode(
    plant = "autotroph",
    detritus = "detritus",
    # detritus = "pathogen", # Now removed
    animal = "predator",
    animal = "detritivore",
    animal = "micropredator",
    animal = "parasitoid",
    para = "macroparasite",
    para = "nonfeeding",
    para = "parasitic castrator",
    para = "trophically transmitted parasite"
  )

feed_type_bsq <- 
  flux_bsq[[1]] %>% 
  pull(consumer_strategy_stage) %>% 
  as.factor() %>%
  fct_recode(
    autotroph = "autotroph",
    detritus = "detritus",
    # detritus = "pathogen", # Now removed
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

mat_bsq <- 
  links %>% 
  mutate(resource_id = as.character(resource_species_id_stage_id),
         consumer_id = as.character(consumer_species_id_stage_id)) %>% 
  dplyr::select(resource_id, consumer_id) %>% 
  graph_from_data_frame(directed = T, vertices = flux_bsq[[1]]$species_id_stage_id) %>% 
  as_adj() %>% 
  as.matrix()

# Vectors -----------------------------------------------------------------

M_bsq <- 
  flux_bsq %>% 
  lapply(pull, M)

N_bsq <- 
  flux_bsq %>% 
  lapply(pull, N)

B_bsq <- 
  flux_bsq %>% 
  map(mutate,
      B = M * N) %>% 
  map(pull, B)

# Define Metabolic Types --------------------------------------------------

met_types <-  c("ecto_vert", "endo_vert", "invert")


# Define Losses and Efficiencies  -----------------------------------------

losses_bsq <- 
  lapply(M_bsq, function(x){
    losses = rep(NA, length(x))
    ecto.vert = met_types == "ecto_vert"
    endo.vert = met_types == "endo_vert"
    inv = met_types == "invert"
    losses[ecto.vert] = 18.18 * x[ecto.vert] ^ (-0.29)
    losses[endo.vert] = 19.5 * x[endo.vert] ^ (-0.29)
    losses[inv] = 18.18 * x[inv] ^ (-0.29)
    losses
  })

efficiencies_bsq <- 
  lapply(M_bsq, function(x){
    efficiencies_bsq = rep(NA, length(x[[1]]))
    efficiencies_bsq[con_type_bsq == "animal"] = 0.906
    efficiencies_bsq[con_type_bsq == "para"] = 0.906
    efficiencies_bsq[con_type_bsq == "plant"] = 0.545 
    efficiencies_bsq[con_type_bsq == "detritus"] = 0.158
    efficiencies_bsq
  })


# Fluxweb -----------------------------------------------------------------

bsq_fluxes <- 
  lapply(1:length(N_bsq), function(x){
    try({fluxing(
      mat = mat_bsq,
      biomasses = B_bsq[[x]],
      losses = losses_bsq[[x]],
      efficiencies = efficiencies_bsq[[x]]
    )})
  })

# which simulations worked?
bsq_fluxes <- Filter(is.numeric, bsq_fluxes)

# length(bsq_fluxes)
# pull the ones that did work
# bsq_fluxes <- bsq_fluxes[1:100]


# Graphs ------------------------------------------------------------------

# network structure:
graph <- 
  links %>%
  mutate(resource = as.character(resource_species_id_stage_id),
         consumer = as.character(consumer_species_id_stage_id)) %>% 
  dplyr::select(resource, consumer) %>% 
  graph_from_data_frame(directed = T, vertices = NULL)

write_graph(graph, "results/bsq_graph.txt", format = "edgelist")

# in matrix form
mat <-
  graph %>% 
  get.adjacency() %>% # dcgClass matrix
  as.matrix() # traditional matrix form

TrophInd <- TrophInd(mat)

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
data.frame(as.list(simple_stats)) %>% write_csv("results/stats_webs_bsq.csv")

###########################################################################
# Counting --------------------------------------------------------------
# Now counting overall flux value in webs for BSQ with parasites
############################################################################

# Initialise Vectors ------------------------------------------------------

herb_bsq = rep(NA, length(bsq_fluxes))
carn_bsq = rep(NA, length(bsq_fluxes))
detr_bsq = rep(NA, length(bsq_fluxes))
para_bsq = rep(NA, length(bsq_fluxes))
total_bsq = rep(NA, length(bsq_fluxes))

# Fill vectors ------------------------------------------------------------

for (x in 1:length(bsq_fluxes)) {
  # herb_bsq
  herb_bsq[[x]] = sum(rowSums(bsq_fluxes[[x]][org_type_bsq == "plant", ]))
  # carn_bsq
  carn_bsq[[x]] = sum(rowSums(bsq_fluxes[[x]][org_type_bsq == "animal", ]))
  # detr_bsq
  detr_bsq[[x]] = sum(rowSums(bsq_fluxes[[x]][org_type_bsq == "detritus", ]))
  # para_bsq
  para_bsq[[x]] = sum(rowSums(bsq_fluxes[[x]][org_type_bsq == "para", ]))
  # total_bsq
  total_bsq[[x]] = sum(bsq_fluxes[[x]])
}

# Gather Data -------------------------------------------------------------

bsq_dat <- 
  cbind(herb_bsq, carn_bsq, detr_bsq, para_bsq, total_bsq) %>% 
  as_tibble %>% 
  gather(key = interaction, value = flux) %>% 
  mutate(system = "BSQ")

# Write Data --------------------------------------------------------------

write_csv(bsq_dat, "results/bsq_parasites_flux.csv")

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
niche_bsq <- Niche.model(S = S, L = L, N = 500)
print("Done...")

# Get TL of Simulated Webs ------------------------------------------------

print("Get TL of simulated webs")
niche_bsq_TL <- lapply(niche_bsq, NetIndices::TrophInd)
print("Done")

TrophInd$species_id_stage_id <- as.numeric(rownames(TrophInd))

suppressMessages(
  # Suppress the "dplyr "JoiningBy X messge"
  troph_bsq <-
    flux_bsq %>%
    map(left_join, TrophInd) %>%
    map(cbind,
        org_type = org_type_bsq,
        con_type = con_type_bsq) %>%
    bind_rows()
)

# get fluxes for each node
suppressMessages(
  # Suppress the "dplyr "JoiningBy X messge"
  bsq_data_para <- lapply(bsq_fluxes, rowSums) %>%
    # now take the names of those generated and the fluxes
    # put into a dataframe to bind too
    lapply(function(x) {
      test <-
        data.frame(flux = x, species_id_stage_id = as.numeric(names(x)))
    }) %>%
    lapply(left_join, troph_bsq) %>%
    bind_rows() %>%
    # there are duplicates
    distinct(flux, species_id_stage_id, .keep_all = T)
)

write_csv(bsq_data_para, "results/bsq_raw.csv")
# Cut into Trophic Levels -------------------------------------------------

bsq_TL <- 
  troph_bsq %>% 
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
bsq_TL[[1]] <- 
  bsq_TL[[1]] %>% 
  enframe() %>% 
  subset(value != "para") %>% 
  subset(value != "animal") %>% 
  droplevels() %>% 
  pull(value)


bsq_M <- 
  troph_bsq %>% 
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

bsq_N <- 
  troph_bsq %>% 
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

bsq_sim_TL <- 
  niche_bsq_TL %>% 
  # add organism column, and cut TL's into bins of 0.5
  map(mutate,
      TL_cut = cut_width(TL, .5),
      org_type = NA,
      M = NA,
      N = NA)

# Sampling and matching ---------------------------------------------------

out_bsq <- list()

for (i in 1:length(bsq_sim_TL)) {
  # gets the list structure of the food web we want to sample into
  x <- bsq_sim_TL[[i]]
  y <- group_split(x, TL_cut)
  z <- map(y, pull, org_type)
  # create M list of structure z
  out_bsq$m[[i]] <- z
  # create N list of structure z
  out_bsq$n[[i]] <- z
  # create org list of structure z
  out_bsq$org[[i]] <- z
  # sample the M data
  for (j in 1:length(z)) {
    if (j < length(bsq_N)) {
      out_bsq$m[[i]][[j]] <- sample(bsq_M[[j]], replace = TRUE, length(z[[j]]))
      # now find the index numbers of those sampled
      a <- match(x = out_bsq$m[[i]][[j]], table = bsq_M[[j]])
      # now take these index numbers and take their respective values for N
      out_bsq$n[[i]][[j]] <- bsq_N[[j]][a]
      # now give them the matching animal type
      out_bsq$org[[i]][[j]] <- bsq_TL[[j]][a]
    } else {
      out_bsq$m[[i]][[j]] <- sample(bsq_M[[length(bsq_M)]], replace = TRUE, length(z[[j]]))
      # now find the index numbers of those sampled
      a <- match(x = out_bsq$m[[i]][[j]], table = bsq_M[[length(bsq_M)]])
      # now take these index numbers and take their respective values for N
      out_bsq$n[[i]][[j]] <- bsq_N[[length(bsq_M)]][a]
      # now give them the matching animal type
      out_bsq$org[[i]][[j]] <- bsq_TL[[length(bsq_M)]][a]
    }
  }
}

# assign trophic levels
for (i in 1:length(bsq_sim_TL)) {
  out_bsq$TL[[i]] <- bsq_sim_TL[[i]]$TL
}

# Assign vector for biomass
out_bsq$b <- out_bsq$m

# This calculates biomass from N and M for each organism in the web
for(i in seq(out_bsq$b)){
  for(j in seq(out_bsq$b[[i]])){
    out_bsq$b[[i]][[j]] <- out_bsq$m[[i]][[j]] * out_bsq$n[[i]][[j]]
  }
}


# Testing!! ---------------------------------------------------------------

# Unit test to find plants and check simulation sampling
for (j in 1:500) {
  data <- c()
  for (i in seq_along(out_bsq$org[[j]])) {
    test <- out_bsq$org %>% map(i) %>% unlist
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
  for (i in seq_along(out_bsq$org[[j]])) {
    test <- out_bsq$org %>% map(i) %>% unlist
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

out_bsq$org <- lapply(out_bsq$org, unlist)
out_bsq$m <- lapply(out_bsq$m, unlist)
out_bsq$n <- lapply(out_bsq$n, unlist)
out_bsq$b <- lapply(out_bsq$b, unlist)

# Losses and Efficiencies -------------------------------------------------

losses_bsq_sim <- 
  lapply(out_bsq$m, function(x){
    losses = rep(NA, length(x))
    ecto.vert = met_types == "ecto_vert"
    endo.vert = met_types == "endo_vert"
    inv = met_types == "invert"
    losses[ecto.vert] = 18.18 * x[ecto.vert] ^ (-0.29)
    losses[endo.vert] = 19.5 * x[endo.vert] ^ (-0.29)
    losses[inv] = 18.18 * x[inv] ^ (-0.29)
    losses
  })

efficiencies_bsq_sim <- 
  lapply(seq(out_bsq$org), function(x) {
    efficiencies_bsq_sim = rep(NA, length(out_bsq$org[[x]]))
    efficiencies_bsq_sim[out_bsq$org[[x]] == "animal"] = 0.906
    efficiencies_bsq_sim[out_bsq$org[[x]] == "para"] = 0.906
    efficiencies_bsq_sim[out_bsq$org[[x]] == "plant"] = 0.545
    efficiencies_bsq_sim[out_bsq$org[[x]] == "detritus"] = 0.906
    efficiencies_bsq_sim
  })


# Flux Calc ---------------------------------------------------------------

bsq_fluxes_sim <- 
  lapply(seq(niche_bsq), function(x){
    try(fluxing(
      mat = niche_bsq[[x]],
      biomasses = out_bsq$b[[x]],
      losses = losses_bsq_sim[[x]],
      efficiencies = efficiencies_bsq_sim[[x]]
    ), silent = TRUE)
  })

winners_bsq <- Filter(is.numeric, bsq_fluxes_sim)
length(winners_bsq) # how many without negative flux value

# Which index values of list worked/failed? -------------------------------

winners_bsq_index <- # these ones won
  unlist(lapply(bsq_fluxes_sim, function(x){
    is.numeric(x)
  }))

# Filter index winners ----------------------------------------------------

# need to filter these so we aggregate the right fluxes for each animal type
winners_out_bsq <- list()
winners_out_bsq$m <- out_bsq$m[winners_bsq_index]
winners_out_bsq$n <- out_bsq$n[winners_bsq_index]
winners_out_bsq$b <- out_bsq$b[winners_bsq_index]
winners_out_bsq$org <- out_bsq$org[winners_bsq_index]
# Add TL
winners_out_bsq$TL <- out_bsq$TL[winners_bsq_index]
# Add FLuxes for each row?
winners_out_bsq$flux <- lapply(bsq_fluxes_sim[winners_bsq_index], rowSums)

# get fluxes for each node in the simulated webs
bsq_raw2 <- transpose(winners_out_bsq) %>% 
  lapply(bind_rows) %>% 
  bind_rows(.id = "Iteration") %>% 
  distinct()

write_csv(bsq_raw2, "results/bsq_raw2.csv")

# Aggregate fluxes --------------------------------------------------------

herb_bsqsim = rep(NA, length(winners_bsq))
carn_bsqsim = rep(NA, length(winners_bsq))
detr_bsqsim = rep(NA, length(winners_bsq))
para_bsqsim = rep(NA, length(winners_bsq))
total_bsqsim = rep(NA, length(winners_bsq))

for (i in seq_along(winners_bsq)) {
  # print(i)
  tryCatch({
    # herbivory
    herb_bsqsim[[i]] = sum(rowSums(winners_bsq[[i]][winners_out_bsq$org[[i]] == 'plant', ]))
    # carnivory
    carn_bsqsim[[i]] = sum(rowSums(winners_bsq[[i]][winners_out_bsq$org[[i]] == 'animal', ]))
    # para
    para_bsqsim[[i]] = sum(rowSums(winners_bsq[[i]][winners_out_bsq$org[[i]] == 'para', ]))
    # detritivory
    detr_bsqsim[[i]] = sum(rowSums(winners_bsq[[i]][winners_out_bsq$org[[i]] == 'detritus', ]))
    # total 
    total_bsqsim[[i]] = sum(winners_bsq[[i]])
  }, error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
}

# errors where detritus is only a vector and not array because only one detritus in TL


# Get web stats for those niche model webs that won and we used -----------

niche_stats <- lapply(niche_bsq[winners_bsq_index], Stats_matrix)
bind_rows(niche_stats) %>% sapply(mean) %>% enframe() %>% pivot_wider() %>% 
  write_csv("results/stats_nichemodelwebs_bsq.csv")


# Get webs as matrix to run comparative analysis --------------------------

niche_webs <- niche_bsq[winners_bsq_index] %>% lapply(graph_from_adjacency_matrix)

n <- 1:length(niche_webs)
names(niche_webs) <- sprintf("bsq_nicheweb_%03d", n)

# Write these to a directory

for(i in names(niche_webs)){
  write.graph(niche_webs[[i]], file.path("graphs", "bsq", paste0(i,".txt")))
}

# Gather and Write --------------------------------------------------------

sim_bsq_dat <- 
  cbind(herb_bsqsim, carn_bsqsim, detr_bsqsim, total_bsqsim, para_bsqsim) %>% 
  as_tibble() %>% 
  gather(key = interaction, value = flux) %>% 
  mutate(system = "BSQ") %>% 
  drop_na()

write_csv(sim_bsq_dat, "results/bsq_no_para_flux.csv")

