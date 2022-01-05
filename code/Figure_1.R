### PLotting imputed density vs orignial Density per web

library(tidyverse)

imputation_data = 
  bind_rows(
  bsq_im = read_csv("results/bsq_imputed.csv") %>% mutate(web = "BSQ"),
  csm_im = read_csv("results/csm_imputed.csv") %>% mutate(web = "CSM"),
  epb_im = read_csv("results/epb_imputed.csv") %>% mutate(web = "EPB"),
)

orignal_data =
  bind_rows(
  bsq_og = read_csv("results/bsq_original.csv") %>% mutate(web = "BSQ"),
  csm_og = read_csv("results/csm_original.csv") %>% mutate(web = "CSM"),
  epb_og = read_csv("results/epb_original.csv") %>% mutate(web = "EPB")
  )


logM = imputation_data %>% 
  ggplot(aes(x = logM, group = .imp)) +
  geom_density(col = "grey40") +
  geom_density(data = original_data, col = "red", inherit.aes = FALSE, aes(x = logM)) +
  see::theme_lucid() +
  theme(legend.position = "bottom") +
  ylab("Density") +
  facet_wrap(~ web)
logN = imputation_data %>% 
  ggplot(aes(x = logN, group = .imp)) +
  geom_density(col = "grey40") +
  geom_density(data = original_data, col = "red", inherit.aes = FALSE, aes(x = logN)) +
  see::theme_lucid() +
  theme(legend.position = "bottom") +
  ylab("Density") +
  facet_wrap(~ web)
logB = imputation_data %>% 
  ggplot(aes(x = logB, group = .imp)) +
  geom_density(col = "grey40") +
  geom_density(data = original_data, col = "red", inherit.aes = FALSE, aes(x = logB)) +
  see::theme_lucid() +
  theme(legend.position = "bottom") +
  ylab("Density") +
  facet_wrap(~ web)


gridded_plot = cowplot::plot_grid(logM, logN, logB, nrow = 3)












