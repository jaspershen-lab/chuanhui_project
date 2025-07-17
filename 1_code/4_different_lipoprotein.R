library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

load("3_data_analysis/1_data_preparation_lipoprotein/lipoprotein_data.rda")

dir.create("3_data_analysis/4_different_lipoprotein", showWarnings = FALSE)
setwd("3_data_analysis/4_different_lipoprotein")

library(tidymass)

####PCA use all lipoprotein variables
pca_object <-
lipoprotein_data %>% 
  scale_data(center = TRUE) %>% 
  run_pca()

plot <-
lipoprotein_data %>% 
  pca_score_plot(pca_object = pca_object, color_by = "group") +
  scale_color_manual(values = ra_dm_color) +
  scale_fill_manual(values = ra_dm_color)

ggsave(plot,
       filename = "pca_plot_all_lipoprotein.pdf",
       width = 6,
       height = 5)

###biomarker discovery
dm_sample_id <-
  lipoprotein_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(disease == "DM") %>%
  dplyr::pull(sample_id)

ra_sample_id <-
  lipoprotein_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(disease == "RA") %>%
  dplyr::pull(sample_id)

lipoprotein_data <-
  lipoprotein_data %>%
  mutate_fc(
    control_sample_id = dm_sample_id,
    case_sample_id = ra_sample_id,
    mean_median = "mean",
    return_mass_dataset = TRUE
  ) %>%
  scale_data(center = TRUE) %>%
  mutate_p_value(
    control_sample_id = dm_sample_id,
    case_sample_id = ra_sample_id,
    method = "t.test",
    p_adjust_method = "fdr"
  )

volcano_plot <-
  lipoprotein_data %>%
  volcano_plot(
    fc_column_name = "fc",
    p_value_column_name = "p_value_adjust",
    fc_up_cutoff = 1,
    fc_down_cutoff = 1,
    add_text = TRUE,
    text_from = "full_name"
  )

ggsave(volcano_plot,
       filename = "volcano_plot.pdf",
       width = 9,
       height = 8)

write.csv(lipoprotein_data@variable_info,
          file = "variable_info.csv",
          row.names = FALSE)


####PCA use marker
pca_object <-
  lipoprotein_data %>% 
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(p_value_adjust < 0.05) %>% 
  scale_data(center = TRUE) %>% 
  run_pca()

plot <-
  lipoprotein_data %>% 
  pca_score_plot(pca_object = pca_object, color_by = "group") +
  scale_color_manual(values = ra_dm_color) +
  scale_fill_manual(values = ra_dm_color)
plot

ggsave(plot,
       filename = "pca_plot_lipoprotein_biomarker.pdf",
       width = 6,
       height = 5)
