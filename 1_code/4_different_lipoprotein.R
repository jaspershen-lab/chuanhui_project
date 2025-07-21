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
plot
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

save(lipoprotein_data,
     file = "lipoprotein_data.rda",
     compress = "xz")

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
volcano_plot
ggsave(volcano_plot,
       filename = "volcano_plot.pdf",
       width = 9,
       height = 8)

library(openxlsx)

openxlsx::write.xlsx(lipoprotein_data@variable_info,
                     file = "variable_info.xlsx",
                     asTable = TRUE)

openxlsx::write.xlsx(
  lipoprotein_data@variable_info %>% dplyr::filter(p_value < 0.05),
  file = "lipoprotein_markers.xlsx",
  asTable = TRUE
)


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


###Heatmap
expression_data <-
  lipoprotein_data %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(p_value_adjust < 0.05) %>%
  scale_data(center = TRUE) %>%
  extract_expression_data()

rownames(expression_data) <-
  lipoprotein_data %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(p_value_adjust < 0.05) %>%
  extract_variable_info() %>%
  pull(full_name)

sample_info <-
  lipoprotein_data %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(p_value_adjust < 0.05) %>%
  extract_sample_info()


variable_info <-
  lipoprotein_data %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::filter(p_value_adjust < 0.05) %>%
  extract_variable_info()

library(ComplexHeatmap)
# Compute log2FC values
log2fc_values <- log2(variable_info$fc)

# Assign colors: red if > 0, blue if < 0
log2fc_colors <- ifelse(log2fc_values > 0, "red", "blue")

plot <-
  Heatmap(
    as.matrix(expression_data),
    row_split = 2,
    name = "Expression",
    show_row_names = TRUE,
    show_column_names = TRUE,
    cluster_columns = TRUE,
    cluster_rows = TRUE,
    row_title = "Lipoprotein variables",
    column_title = "Samples",
    row_title_side = "left",
    column_title_side = "top",
    row_dend_side = "left",
    column_dend_side = "top",
    row_names_gp = gpar(fontsize = 5),
    column_names_gp = gpar(fontsize = 5),
    border = TRUE,
    top_annotation = HeatmapAnnotation(
      group = sample_info$disease,
      col = list(group = ra_dm_color),
      show_legend = TRUE
    ),
    right_annotation =
      rowAnnotation(
        log2FC = anno_barplot(log(variable_info$fc, 2), gp = gpar(fill = log2fc_colors)),
        log10AdjustedP = anno_points(
          -log10(variable_info$p_value_adjust),
          pch = 21,
          size = unit(2, "mm"),
          gp = gpar(fill = "black")
        )
      )
  )

plot <-
  ggplotify::as.ggplot(plot)
plot

ggsave(plot,
       filename = "heatmap_lipoprotein_biomarker.pdf",
       width = 10,
       height = 6)
