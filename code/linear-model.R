####################################################
## THE EASY STATS PACKAGES see, performance, etc  ##
## should be installed from github like so for    ##
## this script:                                   ##
## remotes::install_github("easystats/easystats") ##


# The linear model apporach and plots
library(tidyverse)
library(parameters)
library(see)
library(performance)
library(ggeffects)

# Read in the data

data <- read_rds("results/fluxresults.rds")

# GLM and ANOVA model -----------------------------------------------------

lm1 <- lm(lflux ~ system+(interaction*simulated), data = data)
lm2 <- lm(lflux ~ (system*interaction*simulated), data = data)

# What are the model parameters?
parameters::model_parameters(lm2)

# Anova check for if interaction term is signif
anova(lm1, lm2)
car::Anova(lm2)

# Parameters
preds <- parameters::model_parameters(lm2)

### Robust parameters CI intervels if needed
# parameters::ci_robust(lm2)

# Eta2 of model
effectsize::eta_squared(lm2)

# Dirty and quick plot to see effects
ggpredict(lm2, terms = c("system", "simulated", "interaction")) %>% plot()


# check the model performance all in one go
performance::check_model(lm2) # Extreme values - distribution for this?

# Plot the model predictions vs the data
data <- cbind(data, predict(lm2, interval = "confidence", level = 0.95))


###### Plotting the data vs model outputs #####

plot_labs <- c(para = "Parasites", carn = "Carnivory", detr = "Detrivory", herb = "Herbivory", total = "Total")


# Basic plotting of mdoel results
data %>% 
  ggplot(aes(x = system, y = lflux, col = simulated)) +
  # geom_point(position = position_jitterdodge(), alpha = 0.05) +
  geom_errorbar(aes(y = fit, ymin = lwr, ymax = upr), position = position_dodge()) +
  facet_grid(~ interaction, labeller = as_labeller(plot_labs))# + theme_Publication()


# Add to the summarised values 
summ_data <- data %>% group_by(system, simulated, interaction) %>% summarise(mean_fit = mean(fit))


# Publication plot
data %>% mutate(interaction = fct_relevel(interaction, c("total", "para", "carn", "herb", "detr"))) %>% 
  ggplot(aes(x = system, y = lflux, col = simulated)) +
  geom_boxplot(aes(x = system, y = lflux, fill = simulated)) +
  # geom_violin(aes(x = system, y = lflux, fill = simulated)) +
  geom_errorbar(aes(y = fit, ymin = lwr, ymax = upr, width = .5, group = simulated), position = position_dodge(width = 0.75), col = "black") +
  geom_point(data = summ_data, aes(y = mean_fit, fill = simulated, group = simulated), position = position_dodge(width = 0.75), pch = 21, col = "black", size = 3) +
  facet_grid(~ interaction, labeller = as_labeller(plot_labs)) + 
  ylab(expression(Log~Flux~(J.yr^{-1}))) +
  xlab("Web") +
  see::theme_lucid() +
  ylim(c(5,20)) + theme(legend.position = "bottom") # + labs(col = "Simulated")

ggsave("plots/Fig2Box.png", width = 297, height = 210, dpi = 300, units = "mm", bg = "white")

