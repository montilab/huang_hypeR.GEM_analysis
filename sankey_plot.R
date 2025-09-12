#' Create 3-layer Sankey Diagram
#'
#' @param hypeR_GEM_obj an hypeR-GEM object 
#' @param msets_list a list of MSETS, e.g list(refmet_superPathway, refmet_mainPathway)
#' @param gsets_list a list of gene sets, e.g. list(hallmark$genesets, kegg$genesets)
#' @param key the column indicate the metabolite label (e.g. refmet_name, HMDB, etc.)
#' @param font_size the fontsize on the plot
#' @param node_width the node width on the plot
#'
#' @returns a networkD3::sankeyNetwork object
#' @export
#' @importFrom networkD3 sankeyNetwork
sankey_plot <- function(hypeR_GEM_obj,
                        msets_list,
                        gsets_list,
                        key = "refmet_name",
                        font_size = 12,
                        node_width = 24){
  
  ## remove empty element
  msets_list <- msets_list[lengths(msets_list) > 0]
  gsets_list <- gsets_list[lengths(gsets_list) > 0]
  
  ## metabolite -> MSETs edges
  met2MSET_edge_list <- .met2MSET_links(hypeR_GEM_obj,  msets_list, key = key)
  
  ## create metabolite table
  met_table <- .create_met_table(hypeR_GEM_obj, msets_list, key = key)
  
  ## create gene table
  gene_table <- .create_gene_table(hypeR_GEM_obj, gsets_list)
  
  ## From metabolite -> MSETs table to MSETs -> metabolite table
  msets_table <- .create_msets_table(met_table)
  
  ## From gene -> GSETs table to GSETs to gene table
  gsets_table <- .create_gsets_table(gene_table)
  
  ## edge list !!!
  msets2gsets_edge_list <- .msets2gsets_links(met_table, msets_table, gene_table, gsets_table, key=key)
  
  ## sankeyNetwork
  links <- rbind(met2MSET_edge_list, msets2gsets_edge_list)
  
  nodes <- data.frame(name = union(unique(links$source), unique(links$target)))
  
  # after you build `nodes` and `links`:
  first_col_order <- c("up","down")  # desired top→bottom order for 1st layer
  nodes$name <- factor(
    nodes$name,
    levels = c(first_col_order, setdiff(nodes$name, first_col_order))
  )
  nodes <- nodes[order(as.integer(nodes$name)), , drop = FALSE]
  
  
  links$source <- match(links$source, nodes$name) - 1
  links$target <- match(links$target, nodes$name) - 1
  
  p <- networkD3::sankeyNetwork(Links = links, Nodes = nodes,
                                Source = "source", Target = "target", Value = "weight",
                                NodeID = "name", fontSize = font_size, nodeWidth = node_width,
                                sinksRight = FALSE)  # So the layers stay left-to-right
  
  return(p)
  
}



#' Create edge list between metabolite signatures and give MSETs
#'
#' @param hypeR_GEM_obj a hypeR_GEM_obj
#' @param msets_list a list of MSETS, e.g list(refmet_superPathway, refmet_mainPathway)
#' @param key the column indicate the metabolite label (e.g. refmet_name, HMDB, etc.)
#'
#' @returns a edge list between metabolite signatures (by type) and MSETs
#' @importFrom dplyr bind_rows
#' @importFrom purrr flatten
#' @keywords internal
.met2MSET_links <- function(hypeR_GEM_obj,
                            msets_list,
                            key = "refmet_name"){
  
  ## rbind with specified "type"
  ## met_df = data frame containing at least two columns of "key" and "signature_type"
  met_df <- dplyr::bind_rows(hypeR_GEM_obj$mapped_metabolite_signatures, .id = "signature_type")
  
  ## split metabolties by "group"
  met_groups <- split(met_df[[key]], met_df[["signature_type"]])
  
  ## flatten the input list
  msets <- purrr::flatten(msets_list)
  
  ## create an edge list
  edge_list <- base::expand.grid(source = base::names(met_groups),
                                 target = base::names(msets),
                                 stringsAsFactors = FALSE)
  
  # Use mapply to compute the overlap (weight)
  edge_list$weight <- mapply(function(i, j) {
    base::length(base::intersect(met_groups[[i]], msets[[j]]))
  }, edge_list$source, edge_list$target)
  
  edge_list <- edge_list %>%
    dplyr::filter(weight > 0)
  
  return(edge_list)
}


#' Find associated metabolite/gene-sets for a given "id" and a given "sets"
#'
#' @param id a character id, e.g. "APOE"
#' @param sets a given of gene/metabolite set, e.g. "hallmark$genesets"
#'
#' @returns a character indicating associated sets
#' @keywords internal
.find_association <- function(id, sets) {
  # Find associated metabolite/gene-sets for a given "id" and a given "sets"
  # output: "Amino Acid; peptides"
  
  # Ensure each element is atomic
  sets <- lapply(sets, function(x) unlist(x, use.names = FALSE))
  
  # Logical vector: which sets contain `id`
  is_member <- vapply(sets, function(x) id %in% x, logical(1))
  
  # Return the set names where id is present
  names(sets)[is_member]
  
}


#' Find associated metabolite/gene-sets for a given "key" and a list of "sets"
#'
#' @param id a character id,e.g. "APOE"
#' @param sets_list a list of gene/metabolite set, e.g. list(hallmark$genesets, kegg$genesets)
#' @param add_set_name indicate the name of each gene/metabolite set
#'
#' @returns a character vector
#' @keywords internal
.find_associated_sets <- function(id, sets_list, add_set_name = FALSE){
  
  # Find associated metabolite/gene-sets for a given "key" and a list of "sets"
  # output: "Amino Acid_M_metabolonSupPthwy;Organic acids_refmetSuperClass;"
  out <- unlist(lapply(sets_list, .find_association, id = id))
  
  return(paste(out, collapse = ";"))
}


#' Convert the gene_table in "hypR_GEM_obj" to met_table
#'
#' @param hypeR_GEM_obj a hypeR-GEM object
#' @param msets_list a list of MSETS, e.g list(refmet_superPathway, refmet_mainPathway)
#' @param key the column indicate the metabolite label (e.g. refmet_name, HMDB, etc.)
#' @param sep delimiter 
#'
#' @returns a metabolite-based table containing the associated gene and MSET pathway
#' @keywords internal
#' @importFrom dplyr bind_rows select mutate filter group_by summarise rename right_join
#' @importFrom tidyr separate_rows
.create_met_table <- function(hypeR_GEM_obj, msets_list, key='refmet_name',sep=";"){
  
  k <- sym(key)   # turn string "refmet_name" into a symbol
  
  ## Find associated MSets for each metabolite
  met_df <- bind_rows(hypeR_GEM_obj$mapped_metabolite_signatures, .id = "type") %>%
    dplyr::select(!!k, type) %>%
    dplyr::mutate(
      associated_msets = sapply(.data[[key]], .find_associated_sets, sets_list = msets_list))
  
  ## Find associated genes of each metabolite
  met_table <- dplyr::bind_rows(hypeR_GEM_obj$gene_tables, .id = "type") %>% 
    dplyr::select(symbol, associated_metabolites) %>% 
    tidyr::separate_rows(associated_metabolites, sep = sep) %>% 
    dplyr::mutate(associated_metabolites = str_trim(associated_metabolites)) %>%  # trim spaces
    dplyr::filter(associated_metabolites != "", !is.na(associated_metabolites)) %>%# drop empties
    dplyr::group_by(metabolite = associated_metabolites) %>%
    dplyr::summarise(
      # unique() preserves first-appearance order; use sort(unique(.)) for A–Z
      associated_genes = paste(unique(symbol), collapse = sep),
      .groups = "drop"
    ) %>% 
    dplyr::rename(!!key := metabolite) %>% 
    dplyr::right_join(met_df, by=key) 
  
  return(met_table)
}


#' Find associated gene-set for each mapped enzyme-coding gene
#'
#' @param hypeR_GEM_obj a hypeR-GEM object
#' @param gsets_list a list of gene sets, e.g. list(hallmark$genesets, kegg$genesets)
#'
#' @returns a gene-based table containing the associated gene-set
#' @keywords internal
#' @importFrom dplyr bind_rows select
.create_gene_table <- function(hypeR_GEM_obj, gsets_list){
  
  gene_df <- dplyr::bind_rows(hypeR_GEM_obj$gene_tables, .id = "gene_type") %>%
    dplyr::select(name, symbol)
  
  ## Find associated Gsets for each gene
  gene_df$associated_gsets <- sapply(gene_df$symbol, .find_associated_sets, sets_list = gsets_list)
  
  return(gene_df)
}

#' From metabolite-based table to MSET-based table
#'
#' @param met_table the output from ".create_met_table()"
#' @param key the column indicate the metabolite label (e.g. refmet_name, HMDB, etc.)
#'
#' @returns a MSET-based table containing its associated metabolites
#' @keywords internal
#' @importFrom dplyr filter select group_by summarise rename
#' @importFrom tidyr separate_rows
#' @importFrom stringr str_c
.create_msets_table <- function(met_table, key="refmet_name"){
  
  msets_table <- met_table %>%
    dplyr::filter(associated_msets != "") %>%
    dplyr::select(!!as.name(key), associated_msets) %>%
    tidyr::separate_rows(associated_msets, sep = ";") %>%
    dplyr::filter(nzchar(associated_msets)) %>%  # Filter out empty strings
    dplyr::group_by(associated_msets) %>%
    dplyr::summarise(
      associated_metabolite = stringr::str_c(!!as.name(key), collapse = ";"),
      .groups = "drop"
    ) %>%
    dplyr::rename(msets = associated_msets)
  
  return(msets_table)
}

#' From gene-based table to geneset-based table
#'
#' @param gene_table the output from ".create_gene_table()"
#'
#' @returns a geneset-based table containing its associated genes
#' @keywords internal
#' @importFrom dplyr filter select group_by summarise rename
#' @importFrom tidyr separate_rows
#' @importFrom stringr str_c
.create_gsets_table <- function(gene_table){
  
  gsets_table <- gene_table %>%
    dplyr::select(symbol, associated_gsets) %>%
    dplyr::filter(associated_gsets != "") %>%
    tidyr::separate_rows(associated_gsets, sep = ";") %>%
    dplyr::group_by(associated_gsets) %>%
    dplyr::summarise(
      associated_gene = stringr::str_c(symbol, collapse = ";"),
      .groups = "drop"
    ) %>%
    dplyr::rename(gsets = associated_gsets)
  
  return(gsets_table)
}

.msets2gsets_links <- function(met_table, msets_table, gene_table, gsets_table, key = "refmet_name") {
  edge_list <- data.frame(source = character(), target = character(), weight = integer(), stringsAsFactors = FALSE)
  
  for (i in seq_len(nrow(msets_table))) {
    # safe scalar extraction
    mset_i <- as.character(msets_table$msets[i])
    assoc_met_str <- as.character(msets_table$associated_metabolite[i])
    
    # skip if missing
    if (is.na(mset_i) || is.na(assoc_met_str)) next
    
    # split metabolites; handle empty → skip
    associated_metabolite <- strsplit(assoc_met_str, ";", fixed = TRUE)[[1]]
    associated_metabolite <- associated_metabolite[nzchar(associated_metabolite)]
    if (!length(associated_metabolite)) next
    
    # genes mapped from these metabolites
    met_table_subset <- met_table %>% dplyr::filter(.data[[key]] %in% associated_metabolite)
    
    # flatten associated genes safely, drop NAs/empties
    ag <- unlist(strsplit(stats::na.omit(met_table_subset$associated_genes), ";", fixed = TRUE), use.names = FALSE)
    msets_associated_genes <- unique(ag[nzchar(ag)])
    if (!length(msets_associated_genes)) next
    
    # genes → gsets
    gene_table_subset <- gene_table %>% dplyr::filter(.data$symbol %in% msets_associated_genes)
    gs_flat <- unlist(strsplit(stats::na.omit(gene_table_subset$associated_gsets), ";", fixed = TRUE), use.names = FALSE)
    associated_gsets <- unique(gs_flat[nzchar(gs_flat)])
    if (!length(associated_gsets)) next
    
    # compute weights (overlap sizes)
    weight <- vapply(associated_gsets, function(gset) {
      gs_genes_str <- gsets_table %>%
        dplyr::filter(.data$gsets == gset) %>%
        dplyr::pull(.data$associated_gene)
      gs_genes <- unique(unlist(strsplit(stats::na.omit(gs_genes_str), ";", fixed = TRUE), use.names = FALSE))
      length(intersect(msets_associated_genes, gs_genes))
    }, integer(1))
    
    # keep only positive edges
    keep <- weight > 0L
    if (any(keep)) {
      df <- data.frame(
        source = rep(mset_i, sum(keep)),
        target = associated_gsets[keep],
        weight = weight[keep],
        stringsAsFactors = FALSE
      )
      edge_list <- rbind(edge_list, df)
    }
  }
  
  dplyr::distinct(edge_list)
}
