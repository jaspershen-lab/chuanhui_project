library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

load("3_data_analysis/1_data_preparation_lipoprotein/lipoprotein_data.rda")
load("3_data_analysis/2_data_preparation_metabolites/metabolite_data.rda")

dir.create("3_data_analysis/8_correlation_network", showWarnings = FALSE)
setwd("3_data_analysis/8_correlation_network")

library(tidymass)

lipoprotein_data <-
  lipoprotein_data %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(class = "lipoprotein")

metabolite_data <-
  metabolite_data %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(class = "metabolite")

data <-
  rbind(lipoprotein_data, metabolite_data)

cor_data <-
  massstat::cor_mass_dataset(data)

correlation <-
cor_data$correlation %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("from") %>%
  pivot_longer(
    cols = -from,
    names_to = "to",
    values_to = "correlation"
  ) %>% 
  dplyr::filter(from != to) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(edge_id = paste(sort(c(from, to)), collapse = "_")) %>% 
  dplyr::distinct(edge_id, .keep_all = TRUE)


p <-
  cor_data$p %>% 
  as.data.frame() %>% 
  tibble::rownames_to_column("from") %>%
  pivot_longer(
    cols = -from,
    names_to = "to",
    values_to = "p_fdr"
  ) %>% 
  dplyr::filter(from != to) %>% 
  dplyr::rowwise() %>% 
  dplyr::mutate(edge_id = paste(sort(c(from, to)), collapse = "_")) %>% 
  dplyr::distinct(edge_id, .keep_all = TRUE)

cor_data <-
  correlation %>% 
  dplyr::left_join(p, by = c("edge_id", "from", "to"))


sum(cor_data$correlation > 0.5 & cor_data$p_fdr < 0.05)
