library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

load("3_data_analysis/1_data_preparation_lipoprotein/lipoprotein_data.rda")

dir.create("3_data_analysis/4_different_lipoprotein", showWarnings = FALSE)
setwd("3_data_analysis/4_different_lipoprotein")

library(tidymass)

expression_data <-
  lipoprotein_data %>%
  extract_expression_data()


sample_info <-
  lipoprotein_data %>%
  extract_sample_info()

variable_info <-
  lipoprotein_data %>%
  extract_variable_info()

lipoprotein_data_dm <-
  lipoprotein_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(disease == "DM")

lipoprotein_data_ra <-
  lipoprotein_data %>%
  activate_mass_dataset(what = "sample_info") %>%
  dplyr::filter(disease == "RA")


lipoprotein_metrics_correlation_dm <-
  1:nrow(lipoprotein_data_dm@variable_info) %>%
  purrr::map(function(i) {
    cat(i, " ")
    test_pwv <-
      cor.test(
        as.numeric(lipoprotein_data_dm@expression_data[i, ]),
        lipoprotein_data_dm@sample_info$pwv,
        method = "spearman"
      )
    
    test_FMD <-
      cor.test(
        as.numeric(lipoprotein_data_dm@expression_data[i, ]),
        lipoprotein_data_dm@sample_info$FMD,
        method = "spearman"
      )
    
    test_avgcimt <-
      cor.test(
        as.numeric(lipoprotein_data_dm@expression_data[i, ]),
        lipoprotein_data_dm@sample_info$avgcimt,
        method = "spearman"
      )
    
    test_maxcimt <-
      cor.test(
        as.numeric(lipoprotein_data_dm@expression_data[i, ]),
        lipoprotein_data_dm@sample_info$maxcimt,
        method = "spearman"
      )
    
    test_cavi_avg <-
      cor.test(
        as.numeric(lipoprotein_data_dm@expression_data[i, ]),
        as.numeric(lipoprotein_data_dm@sample_info$cavi_avg),
        method = "spearman"
      )
    
    rbind(
      data.frame(
        variable_id = lipoprotein_data_dm@variable_info$variable_id[i],
        metrics = "pmv",
        cor = test_pwv$estimate,
        p_value = test_pwv$p.value
      ),
      data.frame(
        variable_id = lipoprotein_data_dm@variable_info$variable_id[i],
        metrics = "FMD",
        cor = test_FMD$estimate,
        p_value = test_FMD$p.value
      ),
      data.frame(
        variable_id = lipoprotein_data_dm@variable_info$variable_id[i],
        metrics = "avgcimt",
        cor = test_avgcimt$estimate,
        p_value = test_avgcimt$p.value
      ),
      data.frame(
        variable_id = lipoprotein_data_dm@variable_info$variable_id[i],
        metrics = "maxcimt",
        cor = test_maxcimt$estimate,
        p_value = test_maxcimt$p.value
      ),
      data.frame(
        variable_id = lipoprotein_data_dm@variable_info$variable_id[i],
        metrics = "cavi_avg",
        cor = test_cavi_avg$estimate,
        p_value = test_cavi_avg$p.value
      )
    )
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()


library(plyr)
lipoprotein_metrics_correlation_dm <-
  lipoprotein_metrics_correlation_dm %>%
  plyr::dlply(.(metrics), function(x) {
    x$p_value_adjust <- p.adjust(x$p_value, method = "fdr")
    x
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()









lipoprotein_metrics_correlation_ra <-
  1:nrow(lipoprotein_data_ra@variable_info) %>%
  purrr::map(function(i) {
    cat(i, " ")
    test_pwv <-
      cor.test(
        as.numeric(lipoprotein_data_ra@expression_data[i, ]),
        lipoprotein_data_ra@sample_info$pwv,
        method = "spearman"
      )
    
    test_FMD <-
      cor.test(
        as.numeric(lipoprotein_data_ra@expression_data[i, ]),
        lipoprotein_data_ra@sample_info$FMD,
        method = "spearman"
      )
    
    test_avgcimt <-
      cor.test(
        as.numeric(lipoprotein_data_ra@expression_data[i, ]),
        lipoprotein_data_ra@sample_info$avgcimt,
        method = "spearman"
      )
    
    test_maxcimt <-
      cor.test(
        as.numeric(lipoprotein_data_ra@expression_data[i, ]),
        lipoprotein_data_ra@sample_info$maxcimt,
        method = "spearman"
      )
    
    test_cavi_avg <-
      cor.test(
        as.numeric(lipoprotein_data_ra@expression_data[i, ]),
        as.numeric(lipoprotein_data_ra@sample_info$cavi_avg),
        method = "spearman"
      )
    
    rbind(
      data.frame(
        variable_id = lipoprotein_data_ra@variable_info$variable_id[i],
        metrics = "pmv",
        cor = test_pwv$estimate,
        p_value = test_pwv$p.value
      ),
      data.frame(
        variable_id = lipoprotein_data_ra@variable_info$variable_id[i],
        metrics = "FMD",
        cor = test_FMD$estimate,
        p_value = test_FMD$p.value
      ),
      data.frame(
        variable_id = lipoprotein_data_ra@variable_info$variable_id[i],
        metrics = "avgcimt",
        cor = test_avgcimt$estimate,
        p_value = test_avgcimt$p.value
      ),
      data.frame(
        variable_id = lipoprotein_data_ra@variable_info$variable_id[i],
        metrics = "maxcimt",
        cor = test_maxcimt$estimate,
        p_value = test_maxcimt$p.value
      ),
      data.frame(
        variable_id = lipoprotein_data_ra@variable_info$variable_id[i],
        metrics = "cavi_avg",
        cor = test_cavi_avg$estimate,
        p_value = test_cavi_avg$p.value
      )
    )
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()


library(plyr)
lipoprotein_metrics_correlation_ra <-
  lipoprotein_metrics_correlation_ra %>%
  plyr::dlply(.(metrics), function(x) {
    x$p_value_adjust <- p.adjust(x$p_value, method = "fdr")
    x
  }) %>%
  do.call(rbind, .) %>%
  as.data.frame()

sum(lipoprotein_metrics_correlation_dm$p_value_adjust < 0.05)
sum(lipoprotein_metrics_correlation_ra$p_value_adjust < 0.05)
