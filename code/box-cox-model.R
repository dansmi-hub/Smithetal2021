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

































