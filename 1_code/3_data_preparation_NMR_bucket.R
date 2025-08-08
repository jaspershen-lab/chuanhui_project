library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

data <- readxl::read_xlsx("2_data/univariate_results_28July2025.xlsx", sheet = 4)

load("3_data_analysis/2_data_preparation_metabolites/metabolite_data.rda")

dir.create("3_data_analysis/3_data_preparation_NMR_bucket",
           showWarnings = FALSE)
setwd("3_data_analysis/3_data_preparation_NMR_bucket")

dim(data)

head(data)

sample_info <-
  data %>%
  dplyr::select("Subject ID") %>%
  dplyr::rename(sample_id = "Subject ID") %>% 
  dplyr::mutate(sample_id = stringr::str_extract(sample_id, "RDP2P[0-9]{1,4}"))

expression_data <-
  data %>%
  dplyr::select(-c(1:2)) %>%
  t() %>%
  as.data.frame()

colnames(expression_data) <-
  sample_info$sample_id

variable_info <-
  data.frame(variable_id = row.names(expression_data))

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

nmr_data@variable_info$variable_id[nmr_data@variable_info$variable_id == "b-glucose"] <- "beta-glucose"
nmr_data@variable_info$variable_id[nmr_data@variable_info$variable_id == "α-glucose"] <- "alpha-glucose"

rownames(nmr_data@expression_data) <-
  nmr_data@variable_info$variable_id

save(nmr_data, file = "nmr_data.rda", compress = "xz")

