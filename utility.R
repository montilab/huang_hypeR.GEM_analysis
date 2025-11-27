# Filter genes based on expression
filtESET <- function(eset, batchCol, threshold) {
  batchVals <- pData(eset)[, batchCol]
  batchUnique <- unique(batchVals)
  
  expr0 <- exprs(eset)
  expr0[is.na(expr0)] <- 0
  expr0 <- expr0 > 0
  
  batchMin <- do.call(cbind, lapply(batchUnique, function(b) {
    rowMeans(expr0[,  batchVals == b])
  })) %>%
    rowMin
  
  geneKeep <- batchMin >= threshold
  
  eset <- eset[geneKeep,]
  
  return(eset)
  
}

## half-min imputaiton
halfmin_impute <- function(dat) {
  halfmin <- matrix(matrixStats::rowMins(dat, na.rm = TRUE) / 2,
                    nrow = nrow(dat), ncol = ncol(dat), dimnames = dimnames(dat)
  )
  ## just checking
  stopifnot(all.equal(
    matrixStats::rowMeans2(halfmin),
    matrixStats::rowMins(dat, na.rm = TRUE) / 2
  ))
  dat[is.na(dat)] <- halfmin[is.na(dat)]
  return(dat)
}


# Outliers to NA
outlier2NA <- function(eset, trim = 0.1, nSD = 4) {
  
  exprs(eset) <- apply(exprs(eset), 1, function(x) {
    tm <- mean(x, trim = trim)
    SDs <- sd(x) * nSD
    bounds <- tm + c(-SDs, SDs)
    xt <- x
    xt[xt < bounds[1]] <- NA; xt[xt > bounds[2]] <- NA
    return(xt)
  }) %>%
    t
  return(eset)
  
}

#' Title Calculate the Jaccard similarity between two vectors
#'
#' @param x a vector
#' @param y a vector
#'
#' @return the Jaccard similarity between x and y
#' @export
jaccard_similarity <- function(x,y){
  intersect_norm <- length(intersect(x,y))
  union_norm <- length(union(x,y))

  return(intersect_norm/union_norm)
}



pair_evaluation <- function(hypeR_GEM_obj,
                            diffanal,
                            background){
  
  
  if(!all(c("symbol","gene_type","estimate") %in% colnames(diffanal))) stop("'symbol','gene_type','estimate' should be in the columns of 'diffanal' data frame")
  
  ## background,  genes in GEM X genes in proteomics/mRNA
  
  ## genes that are differentially expressed in proteomics/mRNA
  signif_gene_up <- diffanal %>% 
    dplyr::filter(gene_type == "up") %>% 
    dplyr::pull("symbol")
  
  signif_gene_dn <- diffanal %>% 
    dplyr::filter(gene_type == "dn") %>% 
    dplyr::pull("symbol")
  
  ## differentially expressed genes X backgound, set B
  signif_gene_in_bg <- diffanal %>% 
    dplyr::filter(symbol %in% background, gene_type != "ns") %>% 
    dplyr::pull("symbol")
  
  signif_gene_up_in_bg <- diffanal %>% 
    dplyr::filter(symbol %in% background, gene_type == "up") %>% 
    dplyr::pull("symbol")
  
  signif_gene_dn_in_bg <- diffanal %>% 
    dplyr::filter(symbol %in% background, gene_type == "dn") %>% 
    dplyr::pull("symbol")
  
  
  ## define the output table
  met2gene_list <- list()
  met2gene_in_bg_list <- list()
  overlap_list <- list()
  
  # ## loop over each signature
  for(i in 1:length(names(hypeR_GEM_obj$gene_tables))){
    name <- names(hypeR_GEM_obj$gene_tables)[i]
    gene_table <- hypeR_GEM_obj$gene_tables[[name]]
    
    ## statistics
    met2gene <- nrow(gene_table)
    
    ## ## GEM-mapped genes X background, set A
    met2gene_in_bg <- intersect(gene_table$symbol, background)
    
    met2gene_list[[name]] <- gene_table$symbol
    met2gene_in_bg_list[[name]] <- met2gene_in_bg
    
    ## unweighted overlap
    overlap <- intersect(met2gene_in_bg,  signif_gene_in_bg)
    overlap_list[[name]] <- overlap
  }
  
  ## grouped 
  ## grouped up & dn signatures
  ## use "unique()" because mapped up&dn genes are not mutually exclusive
  df_grouped <- data.frame()
  unique_met2gene_bg <- unique(unlist(met2gene_in_bg_list))
  unique_signif_gene_bg <- unique(signif_gene_in_bg)
  
  
  hyp_GEM <- length(unique_met2gene_bg)
  DE_gene <- length(unique_signif_gene_bg)
  overlap <- length(intersect(unique_met2gene_bg , unique_signif_gene_bg))
  
  pval <- suppressWarnings(stats::phyper(q=overlap-1,
                                         m=DE_gene,
                                         n=length(background)-DE_gene,
                                         k=hyp_GEM,
                                         lower.tail=FALSE))
  
  res <- c(hyp_GEM, DE_gene, overlap, pval, length(background))
  df_grouped <- rbind(df_grouped, res)
  colnames(df_grouped) <- c("hyp_GEM", "DE_gene", "overlap", "pvalue", "background")
  
  
  return(list(df_grouped = df_grouped,
              background = background,
              signif_gene_up = signif_gene_up,
              signif_gene_dn = signif_gene_dn,
              signif_gene_up_bg = signif_gene_up_in_bg,
              signif_gene_dn_bg = signif_gene_dn_in_bg,
              met2gene_list = met2gene_list,
              met2gene_bg_list= met2gene_in_bg_list,
              overlap_list = overlap_list))
}

pair_evaluation_enrichment <- function(hyp_met,
                                       hyp_diffanal,
                                       background,
                                       threshold=0.05){
  ## Union of significant pathway
  met_pathway <- c()
  for(i in 1:length(names(hyp_met$data))){
    name <- names(hyp_diffanal$data)[i]
    tmp <- hyp_met$data[[name]]$data %>% 
      dplyr::filter(fdr <= threshold) %>% 
      dplyr::pull("label")
    met_pathway <- c(met_pathway, tmp)
  }
  met_pathway <- unique(met_pathway)
  
  
  diffanal_pathway <- c()
  for(i in 1:length(names(hyp_diffanal$data))){
    name <- names(hyp_diffanal$data)[i]
    tmp <- hyp_diffanal$data[[name]]$data %>% 
      dplyr::filter(fdr <= threshold) %>% 
      dplyr::pull("label")
    diffanal_pathway <- c(diffanal_pathway, tmp)
  }
  diffanal_pathway <- unique(diffanal_pathway)
  
  
  ## Background
  background <- background
  
  ## overlap
  overlap <- intersect(met_pathway, diffanal_pathway)
  hit_pathway <- paste(overlap, collapse=";")
  
  
  ## hypergeometric test
  pval <- suppressWarnings(stats::phyper(q=length(overlap)-1,
                                         m=length(met_pathway),
                                         n=background-length(met_pathway),
                                         k=length(diffanal_pathway),
                                         lower.tail=FALSE))
  
  res <- c(length(met_pathway), length(diffanal_pathway), length(overlap), pval, background, hit_pathway)
  return(res)
  
}

## helper function for a given signatures to test again multiple class of MSETs/genesets
enrich_helper <- function(signatures,
                          sets_list,
                          background,
                          threshold = 0.05){
  stopifnot(is.list(signatures))
  stopifnot(is.list(sets_list))
  
  ## output
  res <- data.frame()
  
  for(j in 1:length(sets_list)){
    set_name <- names(sets_list)[j]
    sets <- sets_list[[j]]
    
    ## test
    hyp_res <- hypeR(signatures, sets,test="hypergeometric", background=background)
    
    ## output filtered data frame
    df <- data.frame()
    
    for(i in 1:length(names(hyp_res$data))){
      
      signature_name <- names(hyp_res$data)[i]
      
      tmp <- hyp_res$data[[signature_name]][["data"]] %>% 
        dplyr::mutate(signature_name = signature_name) %>% 
        dplyr::filter(fdr <= threshold)
      
      df <- rbind(df, tmp)
    }
    
    df <- df %>% 
      dplyr::mutate(set_name = set_name)
    
    res <- rbind(res, df)
  }
  
  return(res)
  
}

## Filter MSETs
mset_reduce <- function(mset,
                        background,
                        min_size){
  
  mset_flt <- lapply(mset, intersect, background)
  mset_flt <- mset_flt[lengths(mset_flt) >= min_size]
  
  return(mset_flt)
}

reduce_gsets <- function(gsets, background, min_size = 5, verbose = TRUE) {
  new_gsets <- gsets$clone()
  new_gsets$genesets <- lapply(new_gsets$genesets, intersect, background)
  new_gsets$genesets <- new_gsets$genesets[lengths(new_gsets$genesets)>=5]
  if(verbose) cat(length(gsets$genesets), "-->", length(new_gsets$genesets))
  return(new_gsets)
}


#' Read multiple sheets of an excel file
#'
#' @param path path to the excel file
#' @param excluded specified sheets to be excluded
#' @param range read documentation of "readxl::read_excel()"
#' @param col_names read documentation of "readxl::read_excel()"
#' @param col_types read documentation of "readxl::read_excel()"
#' @param na read documentation of "readxl::read_excel()"
#' @param trim_ws read documentation of "readxl::read_excel()"
#' @param skip read documentation of "readxl::read_excel()"
#' @param n_max read documentation of "readxl::read_excel()"
#' @param guess_max read documentation of "readxl::read_excel()"
#' @param progress read documentation of "readxl::read_excel()"
#' @param .name_repair read documentation of "readxl::read_excel()"

#' @importFrom readxl excel_sheets read_excel readxl_progress
#' @importFrom magrittr %>%

#' @return a list of dataframes
#' @export
read_multiple_sheets <- function(file_path,
                                 excluded = NULL,
                                 range = NULL,
                                 col_names = TRUE,
                                 col_types = NULL,
                                 na = "",
                                 trim_ws = TRUE,
                                 skip = 0,
                                 n_max = Inf,
                                 guess_max = min(1000, n_max),
                                 progress = readxl::readxl_progress(),
                                 .name_repair = "unique") {
  
  # getting info about all excel sheets
  sheets <- readxl::excel_sheets(file_path)
  if(!is.null(excluded)){sheets <- sheets[!(sheets %in% excluded)]}
  dfs <- lapply(sheets, function(x){readxl::read_excel(file_path,
                                                       sheet = x,
                                                       range = range,
                                                       col_names = col_names,
                                                       col_types = col_types,
                                                       na  = na,
                                                       trim_ws = trim_ws,
                                                       skip = skip,
                                                       n_max = n_max,
                                                       guess_max =  guess_max,
                                                       progress = progress,
                                                       .name_repair = .name_repair)}) %>%
    lapply(., as.data.frame)
  
  # assigning names to data frames
  names(dfs) <- sheets
  
  return(dfs)
}

## Downloaded and adapted from
## https://www.metabolomicsworkbench.org/databases/refmet/refmet_convert.zip

library(data.table)
library(curl)
library(stringi)

refmet_convert <- function(
    DF # data.frame w/ one metaboline name column
)
{
  # Note: DF[,1] can be any data frame column containing metabolite names
  mets <- stri_join_list(list(DF[, 1]), sep = "\n")
  h <- new_handle()
  handle_setform(h, metabolite_name = mets)
  
  # run the RefMet request on the Metabolomics Workbench server
  req <- curl_fetch_memory(
    "https://www.metabolomicsworkbench.org/databases/refmet/name_to_refmet_new_min.php",
    handle = h
  )
  # Parse the output
  x <- rawToChar(req$content)
  y <- strsplit(x, "\n")
  refmet <- data.frame(ncol = 7)
  
  for (i in 1:length(y[[1]])) {
    if (nchar(y[[1]][i]) > 1) {
      z <- strsplit(y[[1]][i], "\t")
      for (j in 1:length(z[[1]])) {
        refmet[i, j] <- z[[1]][j]
      }
    }
  }
  refmet <- refmet[rowSums(is.na(refmet)) != ncol(refmet), ]
  colnames(refmet) <- refmet[1, ]
  refmet <- refmet[-c(1), ]
  refmet[is.na(refmet)] <- ""
  
  return(
    refmet |>
      dplyr::rename_with(~ trimws(.x)) |>
      dplyr::rename_with(~ tolower(stringr::str_replace(.x,"[ ]+","_")), everything()) |>
      dplyr::rename(metabolite = "input_name", refmet_name = "standardized_name") |>
      dplyr::mutate(refmet_name = ifelse(refmet_name == "-", NA, refmet_name))
  )
}
DF <- data.frame(metabolite = c(
  "deoxycholate",
  "dummy",
  "lithocholate",
  "lithocholate sulfate (1)",
  "glycochenodeoxycholate",
  "glycodeoxycholate")
)


trim_gomf <- function(gomf){
  names(gomf$genesets) <-
    stringr::str_replace(names(gomf$genesets),"Transmembrane Transporter Activity", "Transm Transp Act")
  names(gomf$genesets) <-
    stringr::str_replace(names(gomf$genesets),"Kinase Activity", "KA")
  names(gomf$genesets) <-
    stringr::str_replace(names(gomf$genesets),"Kinase Regulator Activity", "KRA")
  names(gomf$genesets) <-
    stringr::str_replace(names(gomf$genesets),"Channel Activity", "CA")
  names(gomf$genesets) <-
    stringr::str_replace(names(gomf$genesets),"Transporter Activity", "Transp Act")
  names(gomf$genesets) <-
    stringr::str_replace(names(gomf$genesets),"Catalytic Activity Acting", "Cat Act")
  names(gomf$genesets) <-
    stringr::str_replace(names(gomf$genesets),"Transferase Activity Transferring", "Transf Act")
  names(gomf$genesets) <-
    stringr::str_replace(names(gomf$genesets),"Very Long Chain", "VLC")
  names(gomf$genesets) <-
    stringr::str_replace(names(gomf$genesets),
                         "Hydrolase Activity Acting On Carbon Nitrogen But Not Peptide Bonds",
                         "Hydro Act On Carb Na but not Pept Bonds")
  names(gomf$genesets) <-
    stringr::str_replace(names(gomf$genesets),
                         "Oxidoreductase Activity Acting On",
                         "Oxi Act On")
  names(gomf$genesets) <-
    stringr::str_replace(names(gomf$genesets),
                         "Oxi Act On Paired Donors With Incorporation Or Reduction Of Molecular Oxygen",
                         "Oxi Act On Pair Don Inc Reduct Mol Oxy")
  names(gomf$genesets) <-
    stringr::str_replace(names(gomf$genesets),"Group Of Donors", "GoD")
  
  return(gomf)
}

#####################################################
## GGPLOT VOLCANO
#####################################################
## under development
## https://erikaduan.github.io/posts/2021-01-02-volcano-plots-with-ggplot2/
ggplot_volcano <- function(
    dat, # data.frame with "genes" results
    xaxis, # column name for x-axis values
    yaxis, # column name for y-axis values
    gene_label, # column name for "gene" IDs
    sig_il_genes = NULL, # highlight significant genes
    up_il_genes = NULL, # significant up genes IDs (or logical)
    dn_il_genes = NULL, # significant dn genes IDs (or logical)
    top_n = 5, # number of top "genes" to show w/ labels
    exclude = NULL, # "gene" names to exclude (regexp or char vector)
    yintercept = 0.05, # q-value
    xintercept = c(0.5, 2), # fold-changes
    # display settings
    title = "volcano plot",
    xlab = "fold-change",
    ylab = paste0("-log10(", yaxis, ")"),
    cols = c("up" = "red", "dn" = "blue", "ns" = "grey"),
    sizes = c("up" = 2, "dn" = 2, "ns" = 1),
    alphas = c("up" = 1, "dn" = 1, "ns" = 0.5),
    ggrepel_force = 1,
    ggrepel_force_pull = 1,
    ggrepel_nudge_x = 0,
    ggrepel_nudge_y = 0.1)
{
  ## BEGIN input checks
  stopifnot(is(dat, "data.frame"))
  stopifnot(gene_label %in% colnames(dat))
  stopifnot("gene_type" %in% colnames(dat))
  stopifnot(all(dat$gene_type %in% c("up", "dn", "ns")))
  stopifnot(is.null(exclude) || is.character(exclude))
  ## END input checks

  ## establish y-axis range (handling of zero p-/q-values)
  ## creation of a y_axis column where values are thresholded
  dat <- dat |> dplyr::mutate(y_axis = .data[[yaxis]])
  y_rng <- dat |> pull(yaxis)
  if ( min(y_rng) <= .Machine$double.xmin ) {
    #cat("*** THRESHOLDING ***\n")
    min_y <- min(y_rng[y_rng > .Machine$double.xmin])
    dat <- dat |>
      dplyr::mutate(y_axis = ifelse(y_axis <= .Machine$double.xmin, min_y * 0.1, y_axis))
  }
  dat_flt <- dat

  if (!is.null(exclude)) {
    ## list of "gene" IDs
    if (length(exclude) > 1) {
      dat_flt <- dat |>
        dplyr::filter(!.data[[gene_label]] %in% exclude)
    } ## regular expression
    else {
      dat_flt <- dat |>
        dplyr::filter(!stringr::str_detect(.data[[gene_label]], pattern = exclude))
    }
  }
  ## up-regulated "genes" to display
  if (!is.null(up_il_genes)) {
    ## select top_n
    if (is.logical(up_il_genes)) {
      up_il_genes <- dat_flt |>
        dplyr::filter(gene_type == "up") |>
        dplyr::arrange(y_axis) |>
        dplyr::slice_head(n = top_n)
    } ## select "genes" in input list
    else {
      up_il_genes <- dat_flt |>
        dplyr::filter(.data[[gene_label]] %in% up_il_genes) |>
        dplyr::arrange(y_axis) |>
        dplyr::slice_head(n = top_n)
    }
  }
  ## down-regulated "genes" to display
  if (!is.null(dn_il_genes)) {
    ## select top_n
    if (is.logical(dn_il_genes)) {
      dn_il_genes <- dat_flt |>
        dplyr::filter(gene_type == "dn") |>
        dplyr::arrange(y_axis) |>
        dplyr::slice_head(n = top_n)
    } ## select "genes" in input list
    else {
      dn_il_genes <- dat_flt |>
        dplyr::filter(.data[[gene_label]] %in% dn_il_genes) |>
        dplyr::arrange(y_axis) |>
        dplyr::slice_head(n = top_n)
    }
  }
  ## up.down-regulated "genes" to display
  if (!is.null(sig_il_genes)) {
    if (is.logical(sig_il_genes) && (!is.null(up_il_genes) || !is.null(dn_il_genes))) {
      sig_il_genes <- dplyr::bind_rows(up_il_genes, dn_il_genes)
    } else {
      sig_il_genes <- dat_flt |>
        dplyr::filter(.data[[gene_label]] %in% sig_il_genes) |>
        dplyr::arrange(y_axis) |>
        dplyr::slice_head(n = top_n * 2)
    }
  }
  ## display settings
  ## establish y-axis range
  #y_rng <- dat |> pull(y_axis)
  #if (min(y_rng) <= .Machine$double.xmin) {
  #  y_rng[y_rng <= .Machine$double.xmin] <- min(y_rng[y_rng!=0]) * 0.1
  #}
  min_y <- min(dat |> dplyr::pull(y_axis))
  ## establish x-axis range (symmetric)
  dat_range <- dat |> pull(.data[[xaxis]])
  xlim <- rep(max(abs(signif(dat_range, 3))), 2) * c(-1.1, 1.1)
  ylim <- c(0, -log10(min_y) * 1.2)
  dat <- dat |>
    mutate(gene_type = fct_relevel(gene_type, "up", "dn"))

  p <- ggplot(data = dat, aes(
    x = .data[[xaxis]],
    y = -log10(y_axis)
  )) +
    geom_point(aes(colour = gene_type),
               alpha = 1,
               shape = 16,
               size = 1
    ) +
    geom_hline(
      yintercept = -log10(yintercept),
      linetype = "dashed", linewidth = 0.25
    ) +
    #    annotate("text",x=max(xlim),y=-log10(yintercept), label=yintercept, color="red") +
    geom_vline(
      xintercept = xintercept,
      linetype = "dashed", linewidth = 0.25
    ) +
    scale_colour_manual(values = cols) +
    xlim(xlim) +
    ylim(ylim) +
    labs(title = title, x = xlab, y = ylab, colour = "sig. direction") +
    theme_bw() + # Select theme with a white background
    theme(
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank()
    )
  if (!is.null(up_il_genes)) {
    #cat("up_il_genes:",paste(up_il_genes |> dplyr::pull(.data[[gene_label]]),collapse="\n"))
    p <- p + geom_point(
      data = up_il_genes,
      shape = 21,
      size = 2,
      fill = "firebrick",
      colour = "black"
    )
  }
  if (!is.null(dn_il_genes)) {
    #cat("\ndn_il_genes:",paste(dn_il_genes |> dplyr::pull(.data[[gene_label]]),collapse="\n"))
    p <- p + geom_point(
      data = dn_il_genes,
      shape = 21,
      size = 2,
      fill = "steelblue",
      colour = "black"
    )
  }
  if (!is.null(sig_il_genes)) {
    #cat("\nsig_il_genes:",paste(sig_il_genes |> dplyr::pull(.data[[gene_label]]),collapse="\n"))
    p <- p + ggrepel::geom_label_repel(
      data = sig_il_genes,
      aes(label = .data[[gene_label]]),
      force = ggrepel_force,
      force_pull = ggrepel_force_pull,
      nudge_x = ggrepel_nudge_x,
      nudge_y = ggrepel_nudge_y
    )
  }
  return(p)
}

sigmoid_transformation <- function(p,
                                   a = 1, ## smoothness
                                   b = -1,   # p = 0.1 as half-point
                                   eps = 1e-12) {
  
  # clip to avoid log10(0)
  p <- pmax(pmin(p, 1 - eps), eps)
  
  w <- 1 / (1 + exp(a * (log10(p) - b)))
  return(w)
}


