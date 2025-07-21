library(r4projects)
setwd(get_project_wd())
rm(list = ls())
source('1_code/100_tools.R')

load("3_data_analysis/4_different_lipoprotein/lipoprotein_data.rda")
load("3_data_analysis/5_different_metabolites/metabolite_data.rda")

dir.create("3_data_analysis/8_correlation_network", showWarnings = FALSE)
setwd("3_data_analysis/8_correlation_network")

library(tidymass)

lipoprotein_data <-
  lipoprotein_data %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(class = "lipoprotein")

metabolite_data <-
  metabolite_data %>%
  activate_mass_dataset(what = "variable_info") %>%
  dplyr::mutate(class = "metabolite")

remove_idx <-
  apply(metabolite_data@expression_data, 1, function(x) {
    all(as.numeric(x) == 0)
  }) %>%
  which()

if (length(remove_idx) > 0) {
  metabolite_data <-
    metabolite_data[-remove_idx, ]
}

data <-
  rbind(lipoprotein_data, metabolite_data)

cor_data <-
  massstat::cor_mass_dataset(data)

correlation <-
  cor_data$correlation %>%
  as.data.frame() %>%
  tibble::rownames_to_column("from") %>%
  pivot_longer(cols = -from,
               names_to = "to",
               values_to = "correlation") %>%
  dplyr::filter(from != to) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(edge_id = paste(sort(c(from, to)), collapse = "_")) %>%
  dplyr::distinct(edge_id, .keep_all = TRUE)


p <-
  cor_data$p %>%
  as.data.frame() %>%
  tibble::rownames_to_column("from") %>%
  pivot_longer(cols = -from,
               names_to = "to",
               values_to = "p_fdr") %>%
  dplyr::filter(from != to) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(edge_id = paste(sort(c(from, to)), collapse = "_")) %>%
  dplyr::distinct(edge_id, .keep_all = TRUE)

cor_result <-
  correlation %>%
  dplyr::left_join(p, by = c("edge_id", "from", "to"))


sum(cor_result$correlation > 0.5 & cor_result$p_fdr < 0.05)


library(igraph)
library(ggraph)
library(tidygraph)

edge_data <-
  cor_result %>%
  dplyr::filter(abs(correlation) > 0.6 & p_fdr < 0.05) %>%
  dplyr::select(from, to, correlation, p_fdr)

node_data <-
  data@variable_info %>%
  dplyr::filter(variable_id %in% c(edge_data$from, edge_data$to)) %>%
  dplyr::filter(p_value < 0.05)

edge_data <-
  edge_data %>%
  dplyr::filter(from %in% node_data$variable_id &
                  to %in% node_data$variable_id)


node_data$full_name[is.na(node_data$full_name)] <-
  node_data$variable_id[is.na(node_data$full_name)]

graph_data <-
  tidygraph::tbl_graph(nodes = node_data,
                       edges = edge_data,
                       directed = FALSE) %>%
  dplyr::mutate(degree = tidygraph::centrality_degree())


plot <-
  graph_data %>%
  ggraph(layout = "stress", circular = FALSE) +
  ggraph::geom_edge_link(aes(edge_alpha = -log(p_fdr, 10), edge_color = correlation), show.legend = TRUE) +
  ggraph::geom_node_point(aes(
    size = degree,
    fill = log(fc, 2),
    shape = class
  ), show.legend = TRUE) +
  scale_shape_manual(values = c("lipoprotein" = 21, "metabolite" = 22)) +
  scale_size_continuous(range = c(2, 6), name = "Degree") +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    name = "log2FC"
  ) +
  scale_edge_color_gradient2(
    low = "#00aeef",
    mid = "white",
    high = "#c51f05",
    midpoint = 0,
    limits = c(-1, 1)
  ) +
  theme_graph() +
  geom_node_text(aes(label = full_name), repel = TRUE, size = 3)

library(extrafont)
extrafont::loadfonts()
ggsave(plot,
       filename = "correlation_network.pdf",
       width = 10,
       height = 6)


###detect module

library(igraph)
cluster_data <-
cluster_louvain(graph = graph_data, weights = abs(E(graph_data)$correlation))

module_info <-  
data.frame(
  variable_id = V(graph_data)$variable_id,
  module = cluster_data$membership
) %>%
  dplyr::mutate(module = paste0("module_", module)) %>%
  dplyr::select(variable_id, module) %>% 
  dplyr::left_join(node_data, by = "variable_id") %>% 
  dplyr::arrange(module, fc)


library(openxlsx)

openxlsx::write.xlsx(
  module_info,
  file = "module_info.xlsx",
  colNames = TRUE,
  rowNames = FALSE,
  asTable = FALSE
)

module_info <-
  readxl::read_xlsx("module_info_manual.xlsx")

#####module fc distribution
plot <-
module_info %>% 
  ggplot(aes(x = module, y = log(fc, 2))) +
  geom_hline(yintercept = 0, color = "red") +
  geom_boxplot(aes(color = module), show.legend = TRUE) +
  stat_summary(fun = mean, geom = "text", aes(label = round(..y.., 2)), 
               vjust = -0.8, color = "black", size = 5) +
  geom_jitter(aes(fill = module,
                  size = -log(p_value, 10)), 
              shape = 21,
              show.legend = TRUE) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        legend.position = "bottom") +
  labs(x = "", y = "log2FC (RA/DM)") +
  ggsci::scale_color_bmj() +
  ggsci::scale_fill_bmj() +
  scale_size_continuous(range = c(1, 5), 
                        name = "-log10(p-value)")
plot
ggsave(plot,
       filename = "module_fc_distribution.pdf",
       width = 7,
       height = 6)
