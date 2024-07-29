# Packages
library(GeneBridge)
library(readr)
library(dplyr)
library(purrr)
library(biomaRt)
library(magrittr)
library(KEGGREST)
library(UpSetR)
library(ape)
library(tidyverse)
library(data.table)
library(stringi)
library(geneplast.data)
library(AnnotationHub)

# Metastasis Genes data frame
load("results/Detach_and_disseminate_genes.RDS")
df <- result
rm(result)

#Upset Plot
pivotada <- df %>% 
  dplyr::select(external_gene_name, Signature) %>% 
  dplyr::mutate(n = 1) %>% 
  tidyr::pivot_wider(
    id_cols = external_gene_name,
    names_from = Signature,
    values_from = n,
    values_fn = list(n = length),
    values_fill = list(n = 0),
  )

UpSetR::upset(as.data.frame(pivotada), nsets = 50, nintersects = NA)

#---- GeneBridge
# NCBI Eukaryotes 
load('assets/string_eukaryotes.rda')


# Table with Orthologous Groups and their proteins
#source(https://stringdb-static.org/download/COG.mappings.v11.0.txt.gz)

cogs <- fread(
  "assets/COG.mappings.v11.0.txt",
  header           = F,
  stringsAsFactors = F,
  skip             = 1,
  sep              = "\t",
  col.names        = c("taxid.string_id","cog_id"),
  select           = c(1,4),
  quote            = ""
)

# StringDB - metastasis genes (obs: save as gene_id)
genes_hgnc <- unique(df$external_gene_name)
genes_hgnc_concatenado <- paste0(genes_hgnc, collapse = "%0d") 

# Query by STRING API
req <- RCurl::postForm(
  "https://string-db.org/api/tsv/get_string_ids",
  identifiers = genes_hgnc_concatenado,
  echo_query = "1",
  species = "9606"
)
map_ids <- read.table(text = req, sep = "\t", header = T, quote = "")

map_ids$stringId <- substring(map_ids$stringId, 6, 1000)

### Orthology Data
#----
# Spliting first column into taxid and string_id
separated_ids <- cogs %$% stri_split_fixed(taxid.string_id, pattern = ".", n = 2, simplify = T)

cogs[["taxid"]]     <- separated_ids[, 1]
cogs[["string_id"]] <- separated_ids[, 2]

# Freeing up some memory
rm(separated_ids)
gc()

# keeping only eukaryotes
cogs %<>% dplyr::select(-taxid.string_id) %>% filter(taxid %in% string_eukaryotes[["taxid"]])
gc()

cogdata2 <- cogs %>% dplyr::select(protein_id = string_id, ssp_id = taxid, og_id = cog_id)
cogdata2 <- as.data.frame(cogdata2)

# Subsetting cogs of interest - METASTASIS GENES
gene_cogs <- cogs %>%
  filter(string_id %in% map_ids[["stringId"]]) %>%
  dplyr::select(-taxid) %>%
  group_by(string_id) %>%
  summarise(n = n(), cog_id = paste(cog_id, collapse = " / "))

## Proteins with multiple COGs
gene_cogs %>% filter(n > 1)

# Resolving main proteins
gene_cogs_resolved <- tribble(
  ~string_id, ~cog_id,
  "ENSP00000239462",     "KOG1225",        
  "ENSP00000264808",     "KOG1721",          
  "ENSP00000265562",     "KOG2220",          
  "ENSP00000359085",     "KOG3512",      
  "ENSP00000418112",     "COG5599",      
  "ENSP00000422533",     "KOG3545",
  "ENSP00000494788",     "KOG0033"
)

# Removing unresolved cases and adding manual assignments
gene_cogs %<>%
  filter(n == 1) %>%
 dplyr:: select(-n) %>%
  bind_rows(gene_cogs_resolved)


# Exporting for package use
save(gene_cogs, file = 'assets/gene_cogs')

#----
##NETWORK

identifiers <- map_ids %>% pull(stringId) %>% na.omit %>% paste0(collapse="%0d")

req2 <- RCurl::postForm(
  "https://string-db.org/api/tsv/network",
  identifiers = identifiers,
  required_score = "0",
  species     = "9606"
)

string_edgelist <- read.table(text = req2, sep = "\t", header = T)
string_edgelist <- unique(string_edgelist)

## Recomputing scores
combine_scores <- function(dat, evidences = "all", confLevel = 0.4) {
  if(evidences[1] == "all"){
    edat<-dat[,-c(1,2,ncol(dat))]
  } else {
    if(!all(evidences%in%colnames(dat))){
      stop("NOTE: one or more 'evidences' not listed in 'dat' colnames!")
    }
    edat<-dat[,evidences]
  }
  if (any(edat > 1)) {
    edat <- edat/1000
  }
  edat<-1-edat
  sc<- apply(X = edat, MARGIN = 1, FUN = function(x) 1-prod(x))
  dat <- cbind(dat[,c(1,2)],combined_score = sc)
  idx <- dat$combined_score >= confLevel
  dat <-dat[idx,]
  return(dat)
}

string_edgelist <- combine_scores(string_edgelist, evidences = c("ascore", "escore", "dscore"), confLevel = 0.4)
colnames(string_edgelist) <- c("stringId_A", "stringId_B", "combined_score")

# Remove o species id
string_edgelist$stringId_A <- substring(string_edgelist$stringId_A, 6, 1000)
string_edgelist$stringId_B <- substring(string_edgelist$stringId_B, 6, 1000)

# How many edgelist proteins are absent in gene_ids? (should return 0)
setdiff(
  string_edgelist %$% c(stringId_A, stringId_B),
  map_ids %>% pull(stringId)
) 

# Exporting for package use
save(string_edgelist, file = 'results/string_edgelist')

#----
# ##Analysis
# Gene Bridge Data

ah <- AnnotationHub()
meta <- query(ah, "geneplast")
load(meta[["AH83116"]])
cogdata <- dplyr::select(cogdata, "protein_id", "ssp_id", "og_id" = "cog_id")

ogdata <- unique(rbind(cogdata, cogdata2))

## Run geneplast
ogr <- newBridge(ogdata=ogdata, phyloTree=phyloTree, ogids = gene_cogs$cog_id, refsp="9606")

ogr <- runBridge(ogr, penalty = 2, threshold = 0.5, verbose = TRUE)

ogr <- runPermutation(ogr, nPermutations=1000, verbose=FALSE)

res <- getBridge(ogr, what="results")

save(res, file = 'results/resultado_geneplast.Rdata')

## Naming the rooted clades and getting the final results table
CLADE_NAMES <- "https://raw.githubusercontent.com/dalmolingroup/neurotransmissionevolution/ctenophora_before_porifera/analysis/geneplast_clade_names.tsv"

lca_names <- read_table(CLADE_NAMES)

groot_df <- res %>%
  tibble::rownames_to_column("cog_id") %>%
  dplyr::select(cog_id, root = Root) %>%
  inner_join(lca_names) %>%
  inner_join(gene_cogs) 

View(groot_df)

save(groot_df, file = 'results/groot_df.Rdata')

#-----
library(ggplot2)
library(ggraph)
library(dplyr)
library(tidyr)
library(igraph)
library(purrr)
library(vroom)
library(paletteer)
library(easylayout)

### Ploting network

nodelist <- data.frame(node = unique(c(string_edgelist$stringId_A, string_edgelist$stringId_B)))

merged_paths <- merge(nodelist, groot_df, by.x = "node", by.y = "string_id")

save(merged_paths, file = 'results/merged_paths.Rdata')

# Get interaction network from STRINGdb
get_string_network <-
  function(ids,
           species = "9606",
           required_score = 0) {
    ids_collapsed <- paste0(ids, collapse = "%0d")
    
    jsonlite::fromJSON(
      RCurl::postForm(
        "https://string-db.org/api/json/network",
        identifiers = ids_collapsed,
        echo_query  = "1",
        required_score = as.character(required_score),
        species = species
      ),
    )
  }

net <- get_string_network(merged_paths$node, required_score = 0.4)
net <- combine_scores(net, evidences = c("ascore", "escore", "dscore"), confLevel = 0.4)

net <-  net %>%
  separate(stringId_A,
           into = c("ncbi_taxon_id", "stringId_A"),
           sep = "\\.") %>%
  separate(stringId_B,
           into = c("ncbi_taxon_id", "stringId_B"),
           sep = "\\.") 

network_filtered <- net %>%
  dplyr::select(stringId_A, stringId_B) |>
  distinct()

#---- ScatterPie Plot

pivotada <- df %>% 
  dplyr::select(external_gene_name, Signature) %>% 
  dplyr::mutate(n = 1) %>% 
  tidyr::pivot_wider(
    id_cols = external_gene_name,
    names_from = Signature,
    values_from = n,
    values_fn = list(n = length),
    values_fill = list(n = 0),
  )

source_statements <-
  colnames(pivotada)[2:length(pivotada)]

nodelist <-
  data.frame(node = unique(c(network_filtered$stringId_A, network_filtered$stringId_B))) %>%
  left_join(merged_paths, by = c("node" = "node")) %>%
  left_join(map_ids, by = c("node" = "stringId")) %>%
  left_join(pivotada, by = c("queryItem" = "external_gene_name"))

upset(dplyr::select(nodelist, 
                    "queryItem", 
                    "cell adhesion",                          
                    "extracellular matrix organization",      
                    "epithelial to mesenchymal transition",   
                    "regulation of metallopeptidase activity",
                    "cell junction organization",             
                    "cellular extravasation"),
                    nsets = 50, nintersects = NA,
      sets.bar.color = c(
        "cell adhesion"                               = "#06141FFF"
        ,"cell junction organization"                 = "#72874EFF"
        ,"epithelial to mesenchymal transition"       = "#3D4F7DFF"
        ,"extracellular matrix organization"          = "#742C14FF"
        ,"cellular extravasation"                     = "#046E8FFF"
        ,"regulation of metallopeptidase activity"    = "#E48C2AFF"
      ), mainbar.y.label = "Intersecção dos \nProcessos Biológicos",
      sets.x.label = "Tamanho do set")

# Graph Analysis ---------------------------------------------------------------
connected_nodes <- rle(sort(c(network_filtered[,1], network_filtered[,2])))
connected_nodes <- data.frame(count=connected_nodes$lengths, node=connected_nodes$values)
connected_nodes <- left_join(nodelist, connected_nodes, by = c("node" = "node"))

graph <-
  graph_from_data_frame(network_filtered, directed = FALSE, vertices = nodelist)

layout <- easylayout::vivagraph(graph)
layout <- easylayout::vivagraph(graph, layout = layout, pin_nodes = TRUE, lcc_margin_left = 10)
V(graph)$x <- layout[, 1]
V(graph)$y <- layout[, 2]

ggraph(graph, "manual", x = V(graph)$x, y = V(graph)$y) +
  geom_edge_link0(color = "#90909020") +
  geom_node_point(aes(color = -root), size = 2) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "bottom")

color_mappings <- c("#06141FFF", "#742C14FF", "#3D4F7DFF", "#E48C2AFF", "#72874EFF", "#046E8FFF")

ggraph::ggraph(graph,
               "manual",
               x = V(graph)$x,
               y = V(graph)$y) +
  ggraph:: geom_edge_link0(edge_width = 0.5, color = "#90909020") +
  scatterpie::geom_scatterpie(
    cols = source_statements,
    data = igraph::as_data_frame(graph, "vertices"),
    colour = NA,
    pie_scale = 0.40
  ) +
  geom_node_text(aes(label = nodelist$queryItem, colour = "red"), nudge_x = 0.8, nudge_y = 0.8, size = 2) +
  ggplot2::scale_fill_manual(values = color_mappings, drop = FALSE)

ggsave("results/network_plot.pdf")
ggsave("results/network_plot.png")
save(graph, nodelist, file = "results/ppi_network.rda")

#### Plot root
library(ggplot2)
library(ggraph)
library(dplyr)
library(tidyr)
library(igraph)
library(purrr)
library(tinter)

subset_graph_by_root <-
  function(geneplast_result, root_number, graph) {
    filtered <- geneplast_result %>%
      filter(root >= root_number) %>%
      pull(node)
    
    induced_subgraph(graph, which(V(graph)$name %in% filtered))
  }

adjust_color_by_root <- function(geneplast_result, root_number, graph) {
  filtered <- geneplast_result %>%
    filter(root == root_number) %>%
    pull(node)
  
  V(graph)$color <- ifelse(V(graph)$name %in% filtered, "black", "gray")
  return(graph)
}

# Configure graph collors by genes incrementation
subset_and_adjust_color_by_root <- function(geneplast_result, root_number, graph) {
  subgraph <- subset_graph_by_root(geneplast_result, root_number, graph)
  adjusted_graph <- adjust_color_by_root(geneplast_result, root_number, subgraph)
  return(adjusted_graph)
}

plot_network <- function(graph, title, nodelist, xlims, ylims, legend = "none") {
  
  # Generate color map
  source_statements <-
    colnames(nodelist)[11:length(nodelist)]
  
  color_mappings <- c("#06141FFF", "#742C14FF", "#3D4F7DFF", "#E48C2AFF", "#72874EFF", "#046E8FFF")
  
  names(color_mappings) <- source_statements
  
  vertices <- igraph::as_data_frame(graph, "vertices")
  
  ggraph:: ggraph(graph,
                  "manual",
                  x = V(graph)$x,
                  y = V(graph)$y) +
    ggraph::geom_edge_link0(edge_width = 0.2, color = "#90909020") +
    ggraph::geom_node_point(ggplot2::aes(color = I(V(graph)$color)), size = 0.5) +
    scatterpie::geom_scatterpie(
      aes(x=x, y=y, r=18),
      cols = source_statements,
      data = vertices[rownames(vertices) %in% V(graph)$name[V(graph)$color == "black"],],
      colour = NA,
      pie_scale = 1
    ) +
    geom_node_text(aes(label = V(graph)$queryItem), nudge_x = 0.8, nudge_y = 0.8, size = 0.5, colour = "red") +
    ggplot2::scale_fill_manual(values = color_mappings, drop = FALSE) +
    ggplot2::coord_fixed() +
    ggplot2::scale_x_continuous(limits = xlims) +
    ggplot2::scale_y_continuous(limits = ylims) +
    ggplot2::theme_void() +
    ggplot2::theme(
      legend.position = legend,
      legend.key.size = ggplot2::unit(0.5, 'cm'),
      legend.key.height = ggplot2::unit(0.5, 'cm'),
      legend.key.width = ggplot2::unit(0.5, 'cm'),
      legend.title = ggplot2::element_text(size=6),
      legend.text = ggplot2::element_text(size=6),
      panel.border = ggplot2::element_rect(
        colour = "#161616",
        fill = NA,
        linewidth = 1
      ),
      plot.title = ggplot2::element_text(size = 8, face = "bold")
    ) +
    ggplot2::guides(
      color = "none",
      fill = "none"
    ) +
    ggplot2::labs(fill = "Source:", title = title)
}

geneplast_roots <- merged_paths |>
  arrange(root)

buffer <- c(-50, 50)
xlims <- ceiling(range(V(graph)$x)) + buffer
ylims <- ceiling(range(V(graph)$y)) + buffer

roots <- unique(geneplast_roots$root) %>%
  set_names(unique(geneplast_roots$clade_name))

# Subset graphs by LCAs
subsets <-
  map(roots, ~ subset_and_adjust_color_by_root(geneplast_roots, .x, graph))

# Plot titles
titles <- names(roots)

plots <-
  map2(
    subsets,
    titles,
    plot_network,
    nodelist = nodelist,
    xlims = xlims,
    ylims = ylims,
    legend = "right"
  ) %>%
  discard(is.null)

design <- 
  '01234
   56789
   ABCDE
   FGHIJ
   KLMNO'

patchwork::wrap_plots(
  rev(plots),
  design = design,
  nrow = 6,
  ncol = 6
)

ggsave("results/networks_all_roots.pdf")
ggsave("results/networks_all_roots.png")

# Abundance analysiss

lca_names %<>% rename("lca" = root)

lca_spp <- ogr@spbranches %>%
  rename("taxid" = ssp_id, "species" = ssp_name, "lca" = "branch") %>%
  mutate(taxid_order = row_number()) %>%
  dplyr::select(lca, taxid, taxid_order)

clade_taxids <- lca_spp
clade_names <- lca_names

cog_annotation <- map_ids %>%
  left_join(groot_df, by = c("stringId" = "string_id")) %>%
  left_join(df, by = c("queryItem" = "external_gene_name")) %>%
  distinct(queryItem, cog_id, Signature) %>%
  dplyr::select(cog_id, Signature) %>%
  unique()


cog_abundance_by_taxid <- cogs %>%
  filter(cog_id %in% nodelist[["cog_id"]]) %>%
  count(taxid, cog_id, name = "abundance") %>%
  left_join(cog_annotation, by = "cog_id")

# Mapping species to clade info
ordered_species <- string_eukaryotes %>%
  dplyr::select(taxid, ncbi_name) %>%
  left_join(clade_taxids, by = "taxid") %>%
  left_join(clade_names, by = "lca") %>%
  na.omit() %>% unique() %>%
  mutate(
    ncbi_name  = fct_reorder(ncbi_name, -taxid_order)
    ,clade_name = fct_reorder(clade_name, -taxid_order)
  )

avg_abundance_by_function <- cog_abundance_by_taxid %>%
  group_by(taxid, Signature) %>%
  summarise(avg_abundance = mean(abundance)) %>%
  ungroup() %>%
  # Adding species and clade info
  left_join(ordered_species, by = "taxid") %>%
  unique() %>%
  na.omit()

# Plotting colors and labels
annotation_colors <- c(
  "cell adhesion"                               = "#06141FFF"
  ,"extracellular matrix organization"          = "#742C14FF"
  ,"epithelial to mesenchymal transition"       = "#3D4F7DFF"
  ,"regulation of metallopeptidase activity"    = "#E48C2AFF"
  ,"cell junction organization"                  ="#72874EFF"
  ,"cellular extravasation"                      ="#046E8FFF"
)

annotation_labels <- c(
   "cell adhesion"                              = "Adesão celular"
  ,"extracellular matrix organization"          = "Organização da matriz \nextracelular" 
  ,"epithelial to mesenchymal transition"       = "Transição \nepitélio-mesênquimal"
  ,"regulation of metallopeptidase activity"    = "Regulação da atividade \nda metalopeptidase"
  ,"cell junction organization"                 = "Organização da \njunção celular"
  ,"cellular extravasation"                     = "Extravasão celular"
)



# This vertical line indicates the first metazoan (Amphimedon queenslandica / Ctenophora)
ctenophora_line <- geom_vline(
  xintercept = "Amphimedon queenslandica"
  ,color      = "#FF0000"
  ,linetype   = "11"
  ,alpha      = 1
  ,linewidth  = 0.25
)

# Plotting
# Custom tick function
tick_function <- function(x) {
  seq(x[2], 0, length.out = 3) %>% head(-1) %>% tail(-1) %>% { ceiling(./5)*5 }
}

capped_abundance_by_function <- avg_abundance_by_function %>%
  # mutate(capped_abundance = ifelse(abundance >= 100, 100, abundance)) %>%
  group_by(Signature) %>%
  mutate(
    # max_abundance = max(abundance[lca <= 29])
    max_abundance = avg_abundance[lca <= 29] %>% { mean(.) + 3*sd(.) }
    ,abundance     = ifelse(avg_abundance >= max_abundance, pmin(max_abundance, 100), pmin(avg_abundance, 100))
  )


#Themes ---------------------------------------------------------------------
theme_main <- theme(
  panel.spacing      = unit(2.5, "pt")
  ,strip.background   = element_blank()
  ,panel.grid.major.x = element_blank()
  ,panel.grid.major.y = element_line(linewidth = 0.25, linetype = "dotted", color = "#E0E0E0")
  ,strip.text.x       = element_text(size = 9, angle = 90, hjust = 0, vjust = 0.5, color = "#757575")
  ,strip.text.y       = element_text(size = 10, angle = 0, hjust = 0, vjust = 0.5, color = "#757575")
  ,axis.title         = element_text(size = 15, color = "#424242")
  ,axis.ticks.x       = element_blank()
  ,axis.text.x        = element_blank()
  ,axis.text.y        = element_text(size = 5.5)
  ,legend.position    = "none"
)

theme_supplementary <- theme(
  panel.grid.major.x = element_line(color = "#E0E0E0", linewidth = 0.25, linetype = "dotted")
  ,panel.grid.major.y = element_blank()
  ,strip.text.y       = element_text(size = 7, angle = 0, hjust = 0, vjust = 0.5, color = "#757575")
  ,strip.text.x       = element_text(size = 7, angle = 90, hjust = 0, vjust = 0.5, color = "#757575")
  ,axis.title         = element_text(size = 12, color = "#424242")
  ,axis.ticks         = element_line(colour = "grey20")
  ,axis.text.y        = element_text(size = 6, angle = 0, hjust = 1, vjust = 0.5, color = "#757575")
  ,axis.text.x        = element_text(size = 6)
)

theme_average <- theme(
  panel.spacing      = unit(1, "pt")
  ,axis.title         = element_text(color = "#424242")
  ,axis.text          = element_text(color = "#757575")
  ,axis.text.x        = element_text(size = 7, angle = -45, vjust = 0, hjust = 0)
  ,axis.text.y        = element_text(size = 5)
  ,strip.background   = element_blank()
  ,strip.text         = element_text(color = "#757575")
  ,strip.text.y       = element_text(angle = 0, hjust = 0, vjust = 0.5)
)

theme_big <- theme(
  panel.spacing      = unit(0.5, "pt")
  ,panel.grid.major.x = element_line(linewidth = 0.1, linetype = "dashed")
  ,panel.grid.major.y = element_blank()
  ,strip.background   = element_blank()
  ,strip.text.x       = element_text(size = 8, angle = 90, hjust = 0.5, vjust = 0)
  ,strip.text.y       = element_text(size = 8, angle = 0, hjust = 0, vjust = 0.5)
  ,axis.text.x        = element_text(size = 6, angle = 90, vjust = 0, hjust = 0)
  ,axis.text.y        = element_text(size = 4.5)
  ,axis.ticks         = element_line(size = 0.1)
)

#---------------------------------------------------------------------

# Plotting by species
ggplot(avg_abundance_by_function) +
  # Geoms  ----------------
ctenophora_line +
  geom_bar(
    aes(x = ncbi_name, y = avg_abundance, fill = Signature, color = after_scale(darken(fill, 0.1)))
    ,stat = "identity"
  ) +
  # Labels  ---------------
 xlab("Espécies") +
  ylab("Abundância média de proteínas em grupos ortólogos") +
  #ylab("Average protein abundance in orthologous groups") +
  # Scales ----------------
scale_y_continuous(breaks = tick_function, minor_breaks = NULL) +
  scale_fill_manual(values = annotation_colors %>% darken(0.1)) +
  # Styling ---------------
facet_grid(
  Signature ~ clade_name
  ,scales   = "free"
  ,space    = "free"
  ,labeller = labeller(annotation = annotation_labels)
) +
  theme_classic() + 
  theme_main

# Plotting by species capped
ggplot(capped_abundance_by_function) +
  # Geoms  ----------------
chorarnoflagellata_line +
  geom_bar(
    aes(x = ncbi_name, y = abundance, fill = Signature, color = after_scale(darken(fill, 0.1)))
    ,stat = "identity"
  ) +
  # Labels  ---------------
  xlab("Espécies") +
  ylab("Abundância média de proteínas em grupos ortólogos") +
  #ylab("Average protein abundance in orthologous groups") +
  # Scales ----------------
  scale_y_continuous(breaks = tick_function, minor_breaks = NULL) +
  scale_fill_manual(values = annotation_colors %>% darken(0.1)) +
  # Styling ---------------
facet_grid(
  Signature ~ clade_name
  ,scales   = "free"
  ,space    = "free"
  ,labeller = labeller(annotation = annotation_labels)
) +
  theme_classic() + 
  theme_main

# Ploting by clade
ggplot(avg_abundance_by_function) +
  geom_bar(
    aes(x = clade_name, y = avg_abundance, fill = Signature, color = after_scale(darken(fill, 0.1)))
    ,stat = "summary"
    ,fun  = "mean"
  ) +
  scale_y_continuous(breaks = tick_function, minor_breaks = NULL) +
  scale_fill_manual(values = annotation_colors, guide = "none") +
  facet_grid(
    Signature ~ .
    ,scales   = "free"
    ,space    = "free_y"
    ,labeller = labeller(annotation = sub("\n", "", annotation_labels))
  ) +
  xlab("Clados") +
  ylab("Abundância média por clado") +
  theme_classic() + 
  theme_average

# Cumulative root -------------------------------------------------------------

# Mapping roots and proteins info
node_annotation <- nodelist %>%
  inner_join(gene_cogs, by = c("node" = "string_id", "cog_id")) %>%
  inner_join(df, by = c("queryItem" = "external_gene_name")) %>%
  distinct(queryItem, cog_id, Signature, root, clade_name)

# Função para calcular cumulativo e plotar com linha de 90%
plot_cumulative <- function(data, signature) {
  
  # Obter todas as categorias possíveis de clade_name
  all_clades <- node_annotation %>%
    arrange(desc(root)) %>%
    dplyr:: select(clade_name) %>%
    unique()
  
  # Filtrar e calcular cumulativo
  cumulative_data <- node_annotation %>%
    filter(Signature == signature) %>%
    arrange(desc(root)) %>%
    group_by(-root, clade_name) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(cumulative_sum = cumsum(count)) %>%
    rename("root" = "-root")
  
  cumulative_data <- all_clades %>%
    left_join(cumulative_data)
  
  for (i in 2:nrow(cumulative_data)) {
    if (is.na(cumulative_data$cumulative_sum[i])) {
      cumulative_data$cumulative_sum[i] <- cumulative_data$cumulative_sum[i - 1]
    }
  }
  
  # Selecionar nomes do eixo x
  clade_name <- unique(all_clades$clade_name)
  
  # Calcular 90% do total
  total_genes <- sum(na.omit(cumulative_data$count))
  threshold_90 <- 0.9 * total_genes
  root_90 <- cumulative_data$clade_name[which(cumulative_data$cumulative_sum >= threshold_90)][1]
  cumulative_90 <- cumulative_data$cumulative_sum[which(cumulative_data$cumulative_sum >= threshold_90)][1]
  
  # Criar gráfico
  p <- ggplot(cumulative_data, aes(x = factor(clade_name, levels = clade_name), y = cumulative_sum)) +
    geom_bar(stat = "identity", 
             fill = annotation_colors[signature]) +
    geom_hline(yintercept = threshold_90, linetype = "dashed", color = "red") +
    labs(x = NULL, y = NULL, title = paste(signature)) +
    theme_minimal() +
    theme(axis.text.x = element_blank())
  
  return(p)
}

plot_cumulative_2 <- function(data, signature) {
  
  # Obter todas as categorias possíveis de clade_name
  all_clades <- node_annotation %>%
    arrange(desc(root)) %>%
    dplyr:: select(clade_name) %>%
    unique()
  
  # Filtrar e calcular cumulativo
  cumulative_data <- node_annotation %>%
    filter(Signature == signature) %>%
    arrange(desc(root)) %>%
    group_by(-root, clade_name) %>%
    summarise(count = n(), .groups = "drop") %>%
    mutate(cumulative_sum = cumsum(count)) %>%
    rename("root" = "-root")
  
  cumulative_data <- all_clades %>%
    left_join(cumulative_data)
  
  for (i in 2:nrow(cumulative_data)) {
    if (is.na(cumulative_data$cumulative_sum[i])) {
      cumulative_data$cumulative_sum[i] <- cumulative_data$cumulative_sum[i - 1]
    }
  }
  
  # Selecionar nomes do eixo x
  clade_name <- unique(all_clades$clade_name)
  
  # Calcular 90% do total
  total_genes <- sum(na.omit(cumulative_data$count))
  threshold_90 <- 0.9 * total_genes
  root_90 <- cumulative_data$clade_name[which(cumulative_data$cumulative_sum >= threshold_90)][1]
  cumulative_90 <- cumulative_data$cumulative_sum[which(cumulative_data$cumulative_sum >= threshold_90)][1]
  
  # Criar gráfico
  p <- ggplot(cumulative_data, aes(x = factor(clade_name, levels = clade_name), y = cumulative_sum)) +
    geom_bar(stat = "identity", 
             fill = annotation_colors[signature]) +
    geom_hline(yintercept = threshold_90, linetype = "dashed", color = "red") +
    labs(x = NULL, y = NULL, title = paste(signature)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.margin = margin(1, 1, 0, 1))
  
  return(p)
}
# Lista de características a serem plotadas
signatures <- unique(node_annotation$Signature)

# Criar e salvar plots
plots <- lapply(signatures, function(s) plot_cumulative(node_annotation, s))
plots2 <- lapply(signatures, function(s) plot_cumulative_2(node_annotation, s))


#Criar grafico das raíses
library(hrbrthemes)

roots_seq <- node_annotation %>%
  arrange(desc(root)) %>%
  dplyr:: select(root, clade_name) %>%
  unique()

roots_seq$clade_name <- factor(roots_seq$clade_name, levels = roots_seq$clade_name)



j <- ggplot(roots_seq, aes(x = clade_name, y = root)) +
  geom_segment(aes(x = clade_name, xend = clade_name, y = 0, yend = root), 
               color = "#666666FF", 
               size = 0.7) +
  geom_point(color = "#666666FF", size = 2) +
  geom_text(aes(label = clade_name), vjust = 1, hjust = -0.1, size = 3, color = "#666666FF") +
  theme_ipsum() +
  scale_y_reverse() +
  theme_minimal() +  
  theme(
    legend.position = "none",       
    axis.title.x = element_blank(), 
    axis.title.y = element_blank(),  
    axis.text.y = element_blank(),   
    axis.ticks.y = element_blank(),  
    panel.grid = element_blank(),    
    plot.title = element_blank(),     
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  geom_hline(yintercept = 0, color = "#666666FF", size = 0.5)


plots_combined <- c(plots, list(j))
plots_combined <- c(plots[1:6], plots2[7])

combined_plot <- patchwork::wrap_plots(
  plots_combined,
  nrow = 8,
  ncol = 1
)


print(combined_plot)