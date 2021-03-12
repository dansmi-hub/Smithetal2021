# STAN model --------------------------------------------------------------

####################################################
## THE EASY STATS PACKAGES see, performance, etc  ##
## should be installed from github like so for    ##
## this script:                                   ##
## remotes::install_github("easystats/easystats") ##

# library(brms) didn't use int he end but commented example included
library(rstanarm)
library(tidyverse)
library(emmeans)
library(ggeffects)
library(performance)
library(effectsize)
library(see)
library(modelr)
library(tidybayes)
library(bayestestR)

# Fri Mar 12 13:26:34 2021 ------------------------------
## Read in the data 

data <- read_rds("results/fluxresults.rds")

## Quick note that web:BSQ, interaction:Total and simulated:HasParasites
## is the baseline intercept category 

# ## BRMS version - less efficient 
# # default scaled priors set here:
# get_prior(formula = lflux ~ simulated * interaction * system, data =data)
# 
# # Model:
# mod1brms <- brm(# All interactions included
#   formula = lflux ~ simulated * interaction * system, 
#   # Normally distributed after log transform
#   family = gaussian(), data = data, cores = 4, sample_prior = TRUE)
#  # priors auto set, student T distributions above


# ## rstanarm --- more efficient sampling and precompiled better priors
mod1 <- stan_glm(# All interactions included
  formula = lflux ~ simulated * interaction * system,
  # Normally distributed after log transform
  family = gaussian(), data = data, cores = 4,
  # priors set with to a scale of 10 - weakly informative
  prior_intercept = normal(0, 1, autoscale = TRUE),
  prior = normal(0, 1, autoscale = TRUE))

# Posterior predictive check. Does dark blue line follow lighter lines?
# Does dark blue dot sit in the middle of a circle of dots?
pp_check(mod1, plotfun = "stat_2d", stat = c("mean", "sd"))
pp_check(mod1)

# BRMS stan output of model
summary(mod1)

# Model analogous to anova
mod2 <- stan_glm(# All interactions included
  formula = lflux ~ simulated + interaction + system, 
  # Normally distributed after log transform
  family = gaussian(), data = data, cores = 4, 
  # priors set with to a scale of 10 - weakly informative
  prior_intercept = normal(0, 10, autoscale = TRUE),
  prior = normal(0, 10, autoscale = TRUE))

# Effect size of predictors if needed by standardising parameter values
effectsize::standardize_parameters(mod2) # Gives effect sizes as anova eta2 would


# What fix eff? If you want to manually add coeffs together
fixef(mod1)

# Check parameter estimates vs intercept 
bayestestR::hdi(mod1) %>% plot()

# Quick predict to check model output
ggeffects::ggpredict(mod1, terms = c("simulated", "system", "interaction")) %>% plot()

# Some extra diagnostics and posterior output draws if needed
post <- posterior_predict(mod1)
loomod1 <- loo(mod1, cores = 4)
plot(loomod1) # Pareto all fine no outliers in model

# launch_shinystan(mod1) # Diagnositics if needed, non needed fit well. 

# Plotting

# Summarised hdci width posterior draws to usabel df
plot_data <- mod1 %>%
  emmeans::emmeans( ~ system + interaction + simulated) %>%
  gather_emmeans_draws() %>%
  median_hdci()

# All draws from posterior to useable df
plot_draws <-  mod1 %>%
  emmeans::emmeans( ~ system + interaction + simulated) %>%
  gather_emmeans_draws() %>% 
  as.data.frame()


##### Flux and System and Interaction Plotting ######

# Final plot - but needs themeing
plot_data %>% 
  ggplot(aes(x = system, y = .value, col = simulated, group = simulated)) +
  # geom_point(data = plot_draws, aes(y = .value, ymin = NULL, ymax = NULL),
  #             position = position_jitterdodge(jitter.width = 0.75), alpha = 0.1) +
  # tidybayes::stat_slabh(data = plot_draws) +
  geom_pointinterval(position = position_dodge(width = .75), 
                     aes(ymin = .lower, ymax = .upper))+
  geom_errorbar(position = position_dodge(width = .75), col = "black", 
                aes(ymin = .lower, ymax = .upper, width = .5)) +
  facet_grid(~ interaction) +
  see::theme_modern() +
  ylim(c(5,20))
