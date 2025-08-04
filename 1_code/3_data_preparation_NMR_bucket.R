library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

data <- readxl::read_xlsx("2_data/TTSH_NMR buckets.xlsx", sheet = 1)

load("3_data_analysis/2_data_preparation_metabolites/metabolite_data.rda")

dir.create("3_data_analysis/3_data_preparation_NMR_bucket",
           showWarnings = FALSE)
setwd("3_data_analysis/3_data_preparation_NMR_bucket")

dim(data)

head(data)

sample_info <-
  data %>%
  dplyr::select("Subject ID") %>%
  dplyr::rename(sample_id = "Subject ID")

expression_data <-
  data %>%
  dplyr::select(-"Subject ID") %>%
  t() %>%
  as.data.frame()

colnames(expression_data) <-
  sample_info$sample_id

variable_info <-
  data.frame(variable_id = row.names(expression_data),
             nmr_number = as.numeric(row.names(expression_data)))

sample_info <-
  sample_info %>%
  dplyr::left_join(metabolite_data@sample_info, by = "sample_id") %>%
  dplyr::filter(!is.na(class))


expression_data <-
  expression_data[,sample_info$sample_id]

nmr_data <-
  create_mass_dataset(
    expression_data = expression_data,
    variable_info = variable_info,
    sample_info = sample_info
  )


save(nmr_data, file = "nmr_data.rda", compress = "xz")
