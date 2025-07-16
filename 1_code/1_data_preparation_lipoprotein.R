library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

data1 <- readxl::read_xlsx("2_data/1_lipoprotein quantification.xlsx", sheet = 1)
data2 <- readxl::read_xlsx("2_data/1_lipoprotein quantification.xlsx", sheet = 2)
data3 <- readxl::read_xlsx("2_data/1_lipoprotein quantification.xlsx", sheet = 3)

phenotype_data1 <-
  readxl::read_xlsx("2_data/phenotype_data1.xlsx", sheet = 1)

phenotype_data2 <-
  readxl::read_xlsx("2_data/phenotype_data2.xlsx", sheet = 1)

dir.create("3_data_analysis/1_data_preparation_lipoprotein",
           showWarnings = FALSE)
setwd("3_data_analysis/1_data_preparation_lipoprotein")

dim(phenotype_data1)
dim(phenotype_data2)

intersect(colnames(phenotype_data1), colnames(phenotype_data2))

phenotype_data1 <-
  phenotype_data1 %>%
  dplyr::rename(sample_id = "Subject")

colnames(phenotype_data2)

phenotype_data2 <-
  phenotype_data2 %>%
  dplyr::rename(sample_id = "record_id")

intersect(colnames(phenotype_data1), colnames(phenotype_data2))

phenotype_data <-
  phenotype_data1 %>% 
  dplyr::left_join(phenotype_data2, by = intersect(colnames(phenotype_data1), colnames(phenotype_data2)))


dim(data1)

data1 <-
  convert_nmr_lipoprotein_data(data1)
data2 <-
  convert_nmr_lipoprotein_data(data2)

expression_data1 <-
  data1$expression_data
expression_data2 <-
  data2$expression_data

rownames(expression_data1) == rownames(expression_data2)

expression_data <-
  cbind(expression_data1, expression_data2)

variable_info1 <-
  data1$variable_info
variable_info2 <-
  data2$variable_info


variable_info <-
  variable_info1 %>%
  dplyr::left_join(data3[, c("Lipoprotein variables", "Full Name", "Explanation")], by = c("variable_id" = "Lipoprotein variables")) %>%
  dplyr::rename(full_name = "Full Name", explanation = Explanation)

sample_info1 <-
  data1$sample_info
sample_info2 <-
  data2$sample_info

sample_info <-
  rbind(sample_info1, sample_info2)

sample_info$sample_id == colnames(expression_data)

sample_info$class <- "Subject"

library(tidymass)
lipoprotein_data <-
  create_mass_dataset(
    expression_data = expression_data,
    variable_info = variable_info,
    sample_info = sample_info
  )
save(lipoprotein_data,
     file = "lipoprotein_data.rda",
     compress = "xz")
