library(tidyverse)
library(ggplot2)

convert_nmr_lipoprotein_data <- function(x) {
  variable_info <-
    x[1:2, -c(1:2)] %>%
    t() %>%
    as.data.frame() %>%
    dplyr::rename(variable_id = V1, unit = V2)
  
  sample_info <-
    x[-1, 1:2] %>%
    as.data.frame()
  
  colnames(sample_info) <-
    as.character(sample_info[1, ])
  
  sample_info <-
    sample_info[-1, ] %>%
    dplyr::rename(internal_id = "S/N", sample_id = "Subject ID") %>%
    dplyr::select(sample_id, internal_id)
  
  expression_data <-
    x[-c(1:2), -c(1:2)] %>%
    apply(2, function(x) {
      as.numeric(x)
    }) %>%
    t() %>%
    as.data.frame()
  
  colnames(expression_data) <- sample_info$sample_id
  rownames(expression_data) <- variable_info$variable_id
  
  return(list(
    variable_info = variable_info,
    sample_info = sample_info,
    expression_data = expression_data
  ))
  
  
}