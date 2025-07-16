library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

data <- readxl::read_xlsx("2_data/1_lipoprotein quantification.xlsx", sheet = 1)
data3 <- readxl::read_xlsx("2_data/1_lipoprotein quantification.xlsx", sheet = 3)

phenotype_data <-
  readxl::read_xlsx("2_data/phenotype_data2.xlsx", sheet = 1)

dir.create("3_data_analysis/1_data_preparation_lipoprotein",
           showWarnings = FALSE)
setwd("3_data_analysis/1_data_preparation_lipoprotein")

dim(phenotype_data)

phenotype_data <-
  phenotype_data %>%
  dplyr::rename(sample_id = "record_id") %>% 
  dplyr::mutate(disease = case_when(
    disease == "1" ~ "RA",
    disease == "2" ~ "DM",
  ))

data <-
  convert_nmr_lipoprotein_data(data)

expression_data <-
  data$expression_data

variable_info <-
  data$variable_info

variable_info <-
  variable_info %>%
  dplyr::left_join(data3[, c("Lipoprotein variables", "Full Name", "Explanation")], by = c("variable_id" = "Lipoprotein variables")) %>%
  dplyr::rename(full_name = "Full Name", explanation = Explanation)

sample_info <-
  data$sample_info

sample_info$sample_id == colnames(expression_data)

sample_info$class <- "Subject"

sample_info <-
  sample_info %>%
  dplyr::left_join(phenotype_data, by = "sample_id")

library(tidymass)

lipoprotein_data <-
  create_mass_dataset(
    expression_data = expression_data,
    variable_info = variable_info,
    sample_info = sample_info
  )

save(lipoprotein_data, 
     file = "lipoprotein_data.rda", compress = "xz")
