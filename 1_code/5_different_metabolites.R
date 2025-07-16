library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

load("3_data_analysis/2_data_preparation_metabolites/metabolite_data.rda")

dir.create("3_data_analysis/5_different_metabolites", showWarnings = FALSE)
setwd("3_data_analysis/5_different_metabolites")

library(tidymass)

###biomarker discovery
dm_sample_id <-
  metabolite_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(disease == "DM") %>%
  dplyr::pull(sample_id)

ra_sample_id <-
  metabolite_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(disease == "RA") %>%
  dplyr::pull(sample_id)

metabolite_data <-
  metabolite_data %>%
  mutate_fc(
    control_sample_id = dm_sample_id,
    case_sample_id = ra_sample_id,
    mean_median = "mean",
    return_mass_dataset = TRUE
  ) %>%
  mutate_p_value(
    control_sample_id = dm_sample_id,
    case_sample_id = ra_sample_id,
    method = "t.test",
    p_adjust_method = "fdr"
  )

volcano_plot <-
  metabolite_data %>%
  volcano_plot(
    fc_column_name = "fc",
    p_value_column_name = "p_value_adjust",
    fc_up_cutoff = 1,
    fc_down_cutoff = 1,
    add_text = TRUE,
    text_from = "variable_id"
  )
volcano_plot
ggsave(volcano_plot,
       filename = "volcano_plot.pdf",
       width = 9,
       height = 8)

write.csv(metabolite_data@variable_info,
          file = "variable_info.csv",
          row.names = FALSE)
