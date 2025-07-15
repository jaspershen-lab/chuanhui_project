library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

data <- readxl::read_xlsx("2_data/")
dim(data)

dir.create("3_data_analysis/1_data_preparation_lipoprotein",
           showWarnings = FALSE)
setwd("3_data_analysis/1_data_preparation_lipoprotein")

dim(data)

variable_info <-
  data[1:2, 3:114] %>%
  t() %>%
  as.data.frame() %>%
  dplyr::rename(variable_id = V1, unit = V2)


sample_info <-
  data[-1, 1:2] %>%
  as.data.frame()

colnames(sample_info) <-
  as.character(sample_info[1, ])

sample_info <-
  sample_info[-1, ]

colnames(sample_info) <-
  c("S_N", "subject_id")

expression_data <-
  data[-c(1:2), -c(1:2)] %>%
  apply(2, function(x) {
    as.numeric(x)
  }) %>%
  t() %>%
  as.data.frame()
