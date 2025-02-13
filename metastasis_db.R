# Querying metastasis genes by GOids
library(tidyr)
library(dplyr)
library(httr)
library(vroom)

# Data frame with GOids and annotation
metastasis_ids <- read.table("assets/Detatch_and_Dissemination.csv", sep = ",", header = T)                              
metastasis_ids$Signature <- substr(metastasis_ids$Signature, 12,nchar(metastasis_ids$Signature))                            
metastasis_ids <- separate(metastasis_ids, GO_id, into = c("GO_id", "Description"), sep = 10, extra = "merge")

# Function to search genes by GOids filtering by exoerimental evidences
buscar_genes_por_GO <- function(GO_id) {
  url <- paste0("https://golr-aux.geneontology.io/solr/select?defType=edismax&qt=standard&indent=on&wt=csv&rows=100000&start=0&fl=bioentity_label&facet=true&facet.mincount=1&facet.sort=count&json.nl=arrarr&facet.limit=25&hl=true&hl.simple.pre=%3Cem%20class%3D%22hilite%22%3E&hl.snippets=1000&csv.encapsulator=&csv.separator=%09&csv.header=false&csv.mv.separator=%7C&fq=document_category:%22annotation%22&fq=isa_partof_closure:%22", GO_id,"%22&fq=taxon_subset_closure_label:%22Homo%20sapiens%22&fq=type:%22protein%22&fq=evidence_subset_closure_label:%22experimental%20evidence%22&facet.field=aspect&facet.field=taxon_subset_closure_label&facet.field=type&facet.field=evidence_subset_closure_label&facet.field=regulates_closure_label&facet.field=isa_partof_closure_label&facet.field=annotation_class_label&facet.field=qualifier&facet.field=annotation_extension_class_closure_label&facet.field=assigned_by&facet.field=panther_family_label&q=*%3A*")
  response <- tryCatch({
    vroom(url, delim = "\t", col_names = FALSE)
  }, error = function(e) {
    return(data.frame())
  })
  
  if (ncol(response) == 0) {
    return(character(0))
  }
  
  genes <- response[[1]] 
  return(unique(as.vector(genes)))
}

# Create a list with GOids and their respective genes
obter_genes_por_GO_ids <- function(go_ids) {
  genes_per_term <- list()
  for (go_id in go_ids) {
    genes_per_term[[go_id]] <- buscar_genes_por_GO(go_id)
  }
  return(genes_per_term)
}

# Get genes per GOid
genesPerTerm <- obter_genes_por_GO_ids(metastasis_ids$GO_id)

# Create gene data frame
result <- data.frame(external_gene_name = character(), Signature = character(), stringsAsFactors = FALSE)

for (go_id in names(genesPerTerm)) {
  genes <- genesPerTerm[[go_id]]
  signature <- metastasis_ids$Signature[metastasis_ids$GO_id == go_id]
  if (length(signature) > 0 && length(genes) > 0) {
    temp <- data.frame(external_gene_name = genes, Signature = rep(signature, length(genes)), stringsAsFactors = FALSE)
    result <- rbind(result, temp)
  }
}

result <- unique(result)

save(result, file = "results/Detach_and_disseminate_genes.RDS")
