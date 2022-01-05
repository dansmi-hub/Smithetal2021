### Box-Cox Model for extreme values and weird distribution

####################################################
## THE EASY STATS PACKAGES see, performance, etc  ##
## should be installed from github like so for    ##
## this script:                                   ##
## remotes::install_github("easystats/easystats") ##
####################################################

# The linear model apporach and plots
library(tidyverse)
library(parameters)
library(see)
library(performance)
library(ggeffects)
library(MASS)

# Read in the data

data <- read_rds("results/fluxresults.rds")

# BoxCox Transformed LM -----------------------------------------------------

# Which model?
lm1 <- lm(flux ~ system+(interaction*simulated), data = data)
lm2 <- lm(lflux ~ (system*interaction*simulated), data = data)
# Anova check for if interaction term is signif
anova(lm1, lm2)
# Shows lm2 version is needed with 3way interaction

## Andrews figure mentions lambda of 0.06. Check on same page.

bc2 = MASS::boxcox(flux ~ (system*interaction*simulated), data = data)

# Appropriate boxcox here 
(lambda <- bc2$x[which.max(bc2$y)])


## BoxCox LM ------
lm2 <- lm(((flux^lambda-1)/lambda) ~ (system*interaction*simulated), data = data)

# What are the model parameters?
parameters::model_parameters(lm2)


check_normality(lm2) %>% plot(type = "qq")
# plot(lm2)

# ANOVA  ------------------------------------------------------------------

car::Anova(lm2)
capture.output(car::Anova(lm2), file = "tables/anova_results.txt")


# Parameters
preds <- parameters::model_parameters(lm2)
write.table(preds, file = "tables/model_params.txt")

### Robust parameters CI intervels if needed
parameters::ci_robust(lm2)

# Eta2 of model
effectsize::eta_squared(lm2)

##Write to file
write.table(effectsize::eta_squared(lm2), file = "tables/eta-2.txt")

# Dirty and quick plot to see effects
ggpredict(lm2, terms = c("system", "simulated", "interaction")) %>% plot()

# Plot the model predictions vs the data
data <- cbind(data, predict(lm2, interval = "confidence", level = 0.95))


# Original Box plot -------------------------------------------------------


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
data %>% 
  mutate(interaction = fct_relevel(interaction, c("total", "para", "carn", "herb", "detr"))) %>% 
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


# BeeswarmPLot?

# Y scale is not log but label is log. Is this plot the original flux values 


































