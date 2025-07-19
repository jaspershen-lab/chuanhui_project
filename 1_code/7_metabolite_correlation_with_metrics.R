library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

load("3_data_analysis/2_data_preparation_metabolites/metabolite_data.rda")

dir.create("3_data_analysis/5_different_metabolites", showWarnings = FALSE)
setwd("3_data_analysis/5_different_metabolites")

library(tidymass)

expression_data <-
  metabolite_data %>%
  extract_expression_data()


sample_info <-
  metabolite_data %>%
  extract_sample_info()

variable_info <-
  metabolite_data %>%
  extract_variable_info()

metabolite_data_dm <-
  metabolite_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(disease == "DM")

metabolite_data_ra <-
  metabolite_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(disease == "RA")


# metabolite_metrics_correlation_dm <-
#   1:nrow(metabolite_data_dm@variable_info) %>%
#   purrr::map(function(i) {
#     cat(i, " ")
#     test_pwv <-
#       cor.test(
#         as.numeric(metabolite_data_dm@expression_data[i, ]),
#         metabolite_data_dm@sample_info$pwv,
#         method = "spearman"
#       )
#
#     test_FMD <-
#       cor.test(
#         as.numeric(metabolite_data_dm@expression_data[i, ]),
#         metabolite_data_dm@sample_info$FMD,
#         method = "spearman"
#       )
#
#     test_avgcimt <-
#       cor.test(
#         as.numeric(metabolite_data_dm@expression_data[i, ]),
#         metabolite_data_dm@sample_info$avgcimt,
#         method = "spearman"
#       )
#
#     test_maxcimt <-
#       cor.test(
#         as.numeric(metabolite_data_dm@expression_data[i, ]),
#         metabolite_data_dm@sample_info$maxcimt,
#         method = "spearman"
#       )
#
#     test_cavi_avg <-
#       cor.test(
#         as.numeric(metabolite_data_dm@expression_data[i, ]),
#         as.numeric(metabolite_data_dm@sample_info$cavi_avg),
#         method = "spearman"
#       )
#
#     rbind(
#       data.frame(
#         variable_id = metabolite_data_dm@variable_info$variable_id[i],
#         metrics = "pwv",
#         cor = test_pwv$estimate,
#         p_value = test_pwv$p.value
#       ),
#       data.frame(
#         variable_id = metabolite_data_dm@variable_info$variable_id[i],
#         metrics = "FMD",
#         cor = test_FMD$estimate,
#         p_value = test_FMD$p.value
#       ),
#       data.frame(
#         variable_id = metabolite_data_dm@variable_info$variable_id[i],
#         metrics = "avgcimt",
#         cor = test_avgcimt$estimate,
#         p_value = test_avgcimt$p.value
#       ),
#       data.frame(
#         variable_id = metabolite_data_dm@variable_info$variable_id[i],
#         metrics = "maxcimt",
#         cor = test_maxcimt$estimate,
#         p_value = test_maxcimt$p.value
#       ),
#       data.frame(
#         variable_id = metabolite_data_dm@variable_info$variable_id[i],
#         metrics = "cavi_avg",
#         cor = test_cavi_avg$estimate,
#         p_value = test_cavi_avg$p.value
#       )
#     )
#   }) %>%
#   do.call(rbind, .) %>%
#   as.data.frame()
#
# metabolite_metrics_correlation_dm$cor[is.na(metabolite_metrics_correlation_dm$cor)] <- 0
# metabolite_metrics_correlation_dm$p_value[is.na(metabolite_metrics_correlation_dm$p_value)] <- 1
#
# library(plyr)
#
# metabolite_metrics_correlation_dm <-
#   metabolite_metrics_correlation_dm %>%
#   plyr::dlply(.(metrics), function(x) {
#     x$p_value_adjust <- p.adjust(x$p_value, method = "fdr")
#     x
#   }) %>%
#   do.call(rbind, .) %>%
#   as.data.frame()
#
# metabolite_metrics_correlation_ra <-
#   1:nrow(metabolite_data_ra@variable_info) %>%
#   purrr::map(function(i) {
#     cat(i, " ")
#     test_pwv <-
#       cor.test(
#         as.numeric(metabolite_data_ra@expression_data[i, ]),
#         metabolite_data_ra@sample_info$pwv,
#         method = "spearman"
#       )
#
#     test_FMD <-
#       cor.test(
#         as.numeric(metabolite_data_ra@expression_data[i, ]),
#         metabolite_data_ra@sample_info$FMD,
#         method = "spearman"
#       )
#
#     test_avgcimt <-
#       cor.test(
#         as.numeric(metabolite_data_ra@expression_data[i, ]),
#         metabolite_data_ra@sample_info$avgcimt,
#         method = "spearman"
#       )
#
#     test_maxcimt <-
#       cor.test(
#         as.numeric(metabolite_data_ra@expression_data[i, ]),
#         metabolite_data_ra@sample_info$maxcimt,
#         method = "spearman"
#       )
#
#     test_cavi_avg <-
#       cor.test(
#         as.numeric(metabolite_data_ra@expression_data[i, ]),
#         as.numeric(metabolite_data_ra@sample_info$cavi_avg),
#         method = "spearman"
#       )
#
#     rbind(
#       data.frame(
#         variable_id = metabolite_data_ra@variable_info$variable_id[i],
#         metrics = "pwv",
#         cor = test_pwv$estimate,
#         p_value = test_pwv$p.value
#       ),
#       data.frame(
#         variable_id = metabolite_data_ra@variable_info$variable_id[i],
#         metrics = "FMD",
#         cor = test_FMD$estimate,
#         p_value = test_FMD$p.value
#       ),
#       data.frame(
#         variable_id = metabolite_data_ra@variable_info$variable_id[i],
#         metrics = "avgcimt",
#         cor = test_avgcimt$estimate,
#         p_value = test_avgcimt$p.value
#       ),
#       data.frame(
#         variable_id = metabolite_data_ra@variable_info$variable_id[i],
#         metrics = "maxcimt",
#         cor = test_maxcimt$estimate,
#         p_value = test_maxcimt$p.value
#       ),
#       data.frame(
#         variable_id = metabolite_data_ra@variable_info$variable_id[i],
#         metrics = "cavi_avg",
#         cor = test_cavi_avg$estimate,
#         p_value = test_cavi_avg$p.value
#       )
#     )
#   }) %>%
#   do.call(rbind, .) %>%
#   as.data.frame()
#
# metabolite_metrics_correlation_ra$cor[is.na(metabolite_metrics_correlation_ra$cor)] <- 0
# metabolite_metrics_correlation_ra$p_value[is.na(metabolite_metrics_correlation_ra$p_value)] <- 1
#
# library(plyr)
# metabolite_metrics_correlation_ra <-
#   metabolite_metrics_correlation_ra %>%
#   plyr::dlply(.(metrics), function(x) {
#     x$p_value_adjust <- p.adjust(x$p_value, method = "fdr")
#     x
#   }) %>%
#   do.call(rbind, .) %>%
#   as.data.frame()
#
# sum(metabolite_metrics_correlation_dm$p_value_adjust < 0.05)
# sum(metabolite_metrics_correlation_ra$p_value_adjust < 0.05)
#
# max(metabolite_metrics_correlation_dm$cor)
# max(metabolite_metrics_correlation_ra$cor)
#
# save(metabolite_metrics_correlation_dm,
#      file = "metabolite_metrics_correlation_dm.rda",
#      compress = "xz")
# save(metabolite_metrics_correlation_ra,
#      file = "metabolite_metrics_correlation_ra.rda",
#      compress = "xz")

load("metabolite_metrics_correlation_dm.rda")
load("metabolite_metrics_correlation_ra.rda")

significant_dm <- which(metabolite_metrics_correlation_dm$p_value < 0.05)
significant_ra <- which(metabolite_metrics_correlation_ra$p_value < 0.05)

intersect(significant_dm, significant_ra)

all_index <-
  c(significant_dm, significant_ra) %>%
  unique() %>%
  sort()

dir.create("comparison_plot", showWarnings = FALSE)
for (i in all_index) {
  cat(i, " ")
  value1 <-
    metabolite_data_dm@expression_data[metabolite_metrics_correlation_dm[i, ]$variable_id, ] %>%
    as.numeric()
  value2 <-
    metabolite_data_dm@sample_info %>%
    dplyr::select(metabolite_metrics_correlation_dm[i, ]$metrics)
  value2 <-
    value2[, 1] %>%
    as.numeric()
  
  plot_dm <-
    data.frame(value1 = value1, value2 = value2) %>%
    ggplot(aes(x = value1, y = value2)) +
    geom_point(size = 3) +
    geom_smooth(method = "lm", color = "red") +
    labs(x = variable_info$variable_id[match(metabolite_metrics_correlation_dm[i, ]$variable_id,
                                             variable_info$variable_id)],
         y = metabolite_metrics_correlation_dm[i, ]$metrics,
         title = "DM") +
    theme_bw() +
    geom_text(
      aes(label = paste0(
        "R = ",
        round(metabolite_metrics_correlation_dm[i, ]$cor, 2),
        "\n",
        "p = ",
        formatC(
          metabolite_metrics_correlation_dm[i, ]$p_value,
          format = "e",
          digits = 2
        )
      )),
      x = Inf,
      y = Inf,
      hjust = 1.1,
      vjust = 1.1
    )
  
  value1 <-
    metabolite_data_ra@expression_data[metabolite_metrics_correlation_ra[i, ]$variable_id, ] %>%
    as.numeric()
  value2 <-
    metabolite_data_ra@sample_info %>%
    dplyr::select(metabolite_metrics_correlation_ra[i, ]$metrics)
  value2 <-
    value2[, 1] %>%
    as.numeric()
  
  plot_ra <-
    data.frame(value1 = value1, value2 = value2) %>%
    ggplot(aes(x = value1, y = value2)) +
    geom_point(size = 3) +
    geom_smooth(method = "lm", color = "red") +
    labs(x = variable_info$variable_id[match(metabolite_metrics_correlation_ra[i, ]$variable_id,
                                             variable_info$variable_id)],
         y = metabolite_metrics_correlation_ra[i, ]$metrics,
         title = "RA") +
    theme_bw() +
    geom_text(
      aes(label = paste0(
        "R = ",
        round(metabolite_metrics_correlation_ra[i, ]$cor, 2),
        "\n",
        "p = ",
        formatC(
          metabolite_metrics_correlation_ra[i, ]$p_value,
          format = "e",
          digits = 2
        )
      )),
      x = Inf,
      y = Inf,
      hjust = 1.1,
      vjust = 1.1
    )
  
  library(patchwork)
  plot <-
    plot_dm + plot_ra +
    plot_layout(ncol = 2) +
    plot_annotation(
      title = paste0(
        variable_info$full_name[match(metabolite_metrics_correlation_dm[i, ]$variable_id,
                                      variable_info$variable_id)],
        " vs. ",
        metabolite_metrics_correlation_dm[i, ]$metrics
      )
    )
  
  ggsave(
    plot,
    filename = paste0(
      "comparison_plot/",
      variable_info$variable_id[match(metabolite_metrics_correlation_dm[i, ]$variable_id,
                                      variable_info$variable_id)],
      "_",
      metabolite_metrics_correlation_dm[i, ]$metrics,
      ".pdf"
    ),
    width = 8,
    height = 4
  )
  
}
