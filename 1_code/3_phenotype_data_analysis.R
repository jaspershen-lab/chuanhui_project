library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

load("3_data_analysis/1_data_preparation_lipoprotein/lipoprotein_data.rda")

dir.create("3_data_analysis/3_phenotype_data_analysis", showWarnings = FALSE)
setwd("3_data_analysis/3_phenotype_data_analysis")

sample_info <-
  lipoprotein_data@sample_info


library(gghalves)
library(ggsignif)

plot <-
  sample_info %>%
  dplyr::mutate(disease = factor(disease, levels = c("RA", "DM"))) %>%
  ggplot(aes(x = disease, y = age)) +
  geom_half_boxplot(side = "l", aes(color = disease), show.legend = FALSE) +
  geom_half_violin(
    side = "r",
    aes(fill = disease),
    alpha = 0.5,
    show.legend = FALSE
  ) +
  geom_half_point(
    side = "l",
    aes(fill = disease),
    color = "black",
    shape = 21,
    size = 3,
    show.legend = FALSE
  ) +
  ggsignif::geom_signif(
    comparisons = list(c("RA", "DM")),
    map_signif_level = TRUE,
    textsize = 5
  ) +
  scale_color_manual(values = ra_dm_color) +
  scale_fill_manual(values = ra_dm_color) +
  theme_bw() +
  labs(x = "", y = "Age")

ggsave(plot,
       filename = "age_comparison.pdf",
       width = 4,
       height = 4)




plot <-
  sample_info %>%
  dplyr::mutate(disease = factor(disease, levels = c("RA", "DM"))) %>%
  ggplot(aes(x = disease, y = bmi)) +
  geom_half_boxplot(side = "l", aes(color = disease), show.legend = FALSE) +
  geom_half_violin(
    side = "r",
    aes(fill = disease),
    alpha = 0.5,
    show.legend = FALSE
  ) +
  geom_half_point(
    side = "l",
    aes(fill = disease),
    color = "black",
    shape = 21,
    size = 3,
    show.legend = FALSE
  ) +
  ggsignif::geom_signif(
    comparisons = list(c("RA", "DM")),
    map_signif_level = TRUE,
    textsize = 5
  ) +
  scale_color_manual(values = ra_dm_color) +
  scale_fill_manual(values = ra_dm_color) +
  theme_bw() +
  labs(x = "", y = "BMI")
plot
ggsave(plot,
       filename = "bmi_comparison.pdf",
       width = 4,
       height = 4)








plot <-
  sample_info %>%
  dplyr::mutate(disease = factor(disease, levels = c("RA", "DM"))) %>%
  ggplot(aes(x = disease, y = sbp)) +
  geom_half_boxplot(side = "l", aes(color = disease), show.legend = FALSE) +
  geom_half_violin(
    side = "r",
    aes(fill = disease),
    alpha = 0.5,
    show.legend = FALSE
  ) +
  geom_half_point(
    side = "l",
    aes(fill = disease),
    color = "black",
    shape = 21,
    size = 3,
    show.legend = FALSE
  ) +
  ggsignif::geom_signif(
    comparisons = list(c("RA", "DM")),
    map_signif_level = TRUE,
    textsize = 5
  ) +
  scale_color_manual(values = ra_dm_color) +
  scale_fill_manual(values = ra_dm_color) +
  theme_bw() +
  labs(x = "", y = "SBP")
plot
ggsave(plot,
       filename = "sbp_comparison.pdf",
       width = 4,
       height = 4)









plot <-
  sample_info %>%
  dplyr::mutate(disease = factor(disease, levels = c("RA", "DM"))) %>%
  ggplot(aes(x = disease)) +
  geom_bar(stat = "count", color = "black", aes(fill = sex)) +
  scale_fill_manual(values = sex_color) +
  theme_bw()
plot
ggsave(plot,
       filename = "sex_comparison.pdf",
       width = 4,
       height = 4)




plot <-
  sample_info %>%
  dplyr::mutate(disease = factor(disease, levels = c("RA", "DM"))) %>%
  ggplot(aes(x = disease)) +
  geom_bar(stat = "count", color = "black", aes(fill = sex)) +
  scale_fill_manual(values = sex_color) +
  theme_bw()
plot
ggsave(plot,
       filename = "sex_comparison.pdf",
       width = 4,
       height = 4)



plot <-
  sample_info %>%
  dplyr::mutate(disease = factor(disease, levels = c("RA", "DM"))) %>%
  ggplot(aes(x = disease)) +
  geom_bar(stat = "count", color = "black", aes(fill = HTN)) +
  scale_fill_manual(values = htn_color) +
  theme_bw()
plot
ggsave(plot,
       filename = "htn_comparison.pdf",
       width = 4,
       height = 4)


plot <-
  sample_info %>%
  dplyr::mutate(disease = factor(disease, levels = c("RA", "DM"))) %>%
  ggplot(aes(x = disease)) +
  geom_bar(stat = "count", color = "black", aes(fill = HLD)) +
  scale_fill_manual(values = hld_color) +
  theme_bw()
plot
ggsave(plot,
       filename = "hld_comparison.pdf",
       width = 4,
       height = 4)
