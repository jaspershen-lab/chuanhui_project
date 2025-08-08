library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

data <- readxl::read_xlsx("2_data/TTSH_small molecular quantification.xlsx", sheet = 1)

phenotype_data <-
  readxl::read_xlsx("2_data/phenotype_data2.xlsx", sheet = 1)

dir.create("3_data_analysis/2_data_preparation_metabolites",
           showWarnings = FALSE)
setwd("3_data_analysis/2_data_preparation_metabolites")

dim(phenotype_data)

phenotype_data <-
  phenotype_data %>%
  dplyr::rename(sample_id = "record_id") %>%
  dplyr::mutate(disease = case_when(disease == "1" ~ "RA", disease == "2" ~ "DM", ))

data <-
  convert_nmr_metabolite_data(data)

expression_data <-
  data$expression_data

variable_info <-
  data$variable_info

sample_info <-
  data$sample_info

sample_info$sample_id == colnames(expression_data)

sample_info$class <- "Subject"

sample_info <-
  sample_info %>%
  dplyr::left_join(phenotype_data, by = "sample_id")

library(tidymass)

sample_info$sample_id <-
  stringr::str_replace_all(sample_info$sample_id, " \\(poor quality\\)", "")

colnames(expression_data) <-
  sample_info$sample_id

metabolite_data <-
  create_mass_dataset(
    expression_data = expression_data,
    variable_info = variable_info,
    sample_info = sample_info
  )

###remove some cases with IHD (CVD)
metabolite_data <-
  metabolite_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(
    !internal_id %in% c(
      "A05",
      "A06",
      "A11",
      "A12",
      "A31",
      "A32",
      "A33",
      "A34",
      "A45",
      "A46",
      "A51",
      "A52"
    )
  )

save(metabolite_data, file = "metabolite_data.rda", compress = "xz")

