#' Contruct weighted graph using node weights
#'
#' @param edge_list_ A \code{data.frame}. First two columns will be parsed as a edge list.
#' @import limma
#' @import dplyr
construct_weighted_graph_from_edge_list <- 
  function(edge_list_, expr_matrix_, grouping_vec_) {
    lmFit_results <- lmFit(expr_matrix_, model.matrix(~ 0 + grouping_vec_))
    contrast_fit_results <- contrasts.fit(lmFit_results, c(-1, 1))
    contrast_fit_results <- eBayes(contrast_fit_results)
    node_weights <- abs(contrast_fit_results$t[,1])
    
    
    v1_tb <- data.frame(v1 = names(node_weights),
                        v1_weight = node_weights)
    
    v2_tb <- data.frame(v2 = names(node_weights),
                        v2_weight = node_weights)
    
    edge_list <- edge_list_[, 1:2]
    colnames(edge_list) <- c('v1', 'v2')
    
    # filter edge_list so that genes not defined in the node_weights are discarded
    # also calculate the edge weights
    edge_list <- edge_list %>% 
      inner_join(v1_tb, 'v1') %>% 
      inner_join(v2_tb, 'v2') %>% 
      mutate(weight = (v1_weight + v2_weight)/2) %>% 
      select(v1, v2, weight)
    
    
    v_tb <- data.frame(v = names(node_weights),
                       weight = node_weights,
                       t_stat = contrast_fit_results$t[,1])
    
    weighted_graph <- graph_from_data_frame(edge_list, vertices = v_tb)
    
    return(weighted_graph)
  }




#' Generate communities using Spin Glass algorithm
#' 
#' @import igraph
#' @import doParallel
#' @import foreach
generate_spinglass_communities <-
  function(weighted_graph_, seed_genes_, 
           n_threads_ = 1, gamma_ = 0.5, min_size_ = 1, max_size_ = 100, verbose_level_ = 1, seed_=12345) {
    cl <- makeCluster(n_threads_)
    registerDoParallel(cl)
    
    if (verbose_level_ >= 2) message('num of threads for parallel computing = ', n_threads_)
    if (verbose_level_ >= 2) message(sprintf('seed_genes_ (length: %d) = %s', length(seed_genes_), paste(seed_genes_, collapse=', ')))

    set.seed(seed_)
    modules <- foreach(i = 1:length(seed_genes_)) %dopar% {
      module <- cluster_spinglass(weighted_graph_, 
                                  weights = NULL,
                                  vertex = seed_genes_[i],
                                  spins = 25, 
                                  update.rule = c("config"),
                                  gamma = gamma_)
      module$seed <- seed_genes_[i]
      module
    }
    
    if (verbose_level_ >= 1) message('  - Generating module reports ...')
    
    modules <- modules[sapply(modules, function(x) between(length(x$community), min_size_, max_size_))]
    if (verbose_level_ >= 2) message('num of modules = ', length(modules))
    
    results <- list()
    
    genelists <- list()
    graphs <- list()
    info <- list()
    for (i in 1:length(modules)) {
      seed <- modules[[i]]$seed
      graph <- modules[[i]]$community %>% 
        induced_subgraph(weighted_graph_, .)
      vertex_indices <- modules[[i]]$community
      
      graphs[[seed]] <- graph
      
      info$seed[i] <- seed
      info$size[i] <- vertex_indices %>% length
      info$connectivity[i] <- graph %>% 
        degree_distribution() %>% 
        sum(. * 1:length(.))
      
      genelists[[seed]] <- V(weighted_graph_)$name[vertex_indices]
    }
    info <- data.frame(info)
    
    results$info <- info
    results$genelists <- genelists
    results$graphs <- graphs
    
    if (verbose_level_ >= 2) message('num of modules = ', nrow(results$info))
    if (verbose_level_ >= 2) message('module seeds = ', paste(results$info$seed, collapse = ', '))
    
    stopCluster(cl)
    
    return(results)
  }



#' Do Monte Carlo simulation
#' 
#' @param gene_weights_ A numeric vector. Should be the same one the user supplied to \code{construct_weighted_graph}.
#' @import doParallel
#' @import foreach
get_monte_carlo_simulation_p_values <-
  function(module_graphs_, 
           n_monte_carlo_ = 1000, n_threads_ = 1, verbose_level_ = 1, seed_=12345) {
    cl <- makeCluster(n_threads_)
    registerDoParallel(cl)
    
    set.seed(seed_)
    
    n_modules <- length(module_graphs_)
    
    if (verbose_level_ >= 2) message('num of modules = ', n_modules)
    if (verbose_level_ >= 2) message('num of monte carlo runs = ', n_monte_carlo_)
    if (verbose_level_ >= 2) message('num of threads for parallel computing = ', n_threads_)
    
    mean_weights <- vector(length = n_modules)
    for (i in 1:n_modules) {
      mean_weights[i] <- mean(E(module_graphs_[[i]])$weight)
    }
    
    monte_carlo_result_matrix <- matrix(nrow = n_modules, ncol = n_monte_carlo_)
    for (modules_index in 1:n_modules) {
      edgelist <- as_edgelist(module_graphs_[[modules_index]], name = F)
      vertex_weights <- V(module_graphs_[[modules_index]])$weight
      mean_list <- foreach(run_index = 1:n_monte_carlo_) %dopar% {
        vertex_weights <- sample(vertex_weights)
        mean(vertex_weights[edgelist])
      }
      monte_carlo_result_matrix[modules_index, ] <- unlist(mean_list)
    }
    
    p_values <- apply(monte_carlo_result_matrix > mean_weights, 1, mean)
    
    results <- list()
    results$seed <- names(module_graphs_)
    results$p_values <- p_values
    results <- data.frame(results)
    
    stopCluster(cl)
    
    return(results)
  }


#topGO_analysis <-
#  function(top_module_genelists_, all_genes_, eg_db_, verbose_level_ = 1) {
#    library(eg_db_, character.only = T)
#    
#    if (verbose_level_ >= 1) message('* Doing GO enrichment analysis ... ')
#    if (verbose_level_ >= 2) message('num of modules = ', length(top_module_genelists_))
#    if (verbose_level_ >= 2) message('module seeds = ', paste(names(top_module_genelists_), collapse = ', '))
#    
#    go_analysis_df <- list()
#    for (i in 1:length(top_module_genelists_)) {
#      seed <- names(top_module_genelists_)[i]
#      
#      if (verbose_level_ >= 1) message('  - Doing module with seed gene: ', seed, ' (', i, '/', length(top_module_genelists_), ') ...')
#      
#      go_analysis_df$seed[i] <- seed
#      genelist <- top_module_genelists_[[i]]
#      
#      temp_genelist <- top_module_genelists_[[i]]
#      interested_genes <- as.integer(all_genes_ %in% temp_genelist)
#      names(interested_genes) <- all_genes_
#      
#      sink("/dev/null")
#      GOdata <- invisible(
#        new("topGOdata",
#            ontology = "BP",
#            allGenes = interested_genes,
#            geneSel = function(x){x==1},
#            nodeSize = 5,
#            annot = annFUN.org,
#            mapping = eg_db_, 
#            ID = 'symbol'))
#      
#      
#      resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
#      gen_table <- GenTable(GOdata, classicFisher = resultFisher, numChar=999)
#      
#      sink()
#      
#      go_analysis_df$first_term[i] <- gen_table$Term[1]
#      go_analysis_df$first_fisher[i] <- gen_table$classicFisher[1]
#      go_analysis_df$second_term[i] <- gen_table$Term[2]
#      go_analysis_df$second_fisher[i] <- gen_table$classicFisher[2]
#    }
#    go_analysis_df <- data.frame(go_analysis_df)
#    
#    return(go_analysis_df)
#  }



#' KEGG pathway analysis
#' 
#' @param selected_genes_ A character vector of genes, by default gene symbols. The ID type can be changed
#'    using \code{gene_id_type_} parameter.
#' @param all_genes_ A character vector of genes treated as background. It should be a superset of \code{selected_genes_}.
#' @param organism_ Either "human" or "mouse".
#' @param gene_id_type_ Currently either "symbol" (default), or "entrez_id".
#' @param verbose_level_ 0 for quiet, 1 for simple messages, 2 for additional debug messages.
#' @param pvalueCutoff_,qvalueCutoff_ Cutoff threshold for either p or q values.
#' 
#' @import clusterProfiler
#' @import org.Mm.eg.db
#' @import org.Hs.eg.db
#' @import AnnotationDbi
#' @export
KEGG_analysis <-
  function(selected_genes_,
           all_genes_,
           organism_,
           gene_id_type_='symbol',
           ontology_='CC',
           verbose_level_=1,
           pvalueCutoff_=0.01,
           qvalueCutoff_=0.05
           ) {
    if (verbose_level_ >= 2) message('# [GO_analysis] is called')
    
    if (gene_id_type_ == 'symbol') {
      if (verbose_level_ >= 2) message('# genes are assumed to be symbols')
      
      if (verbose_level_ >= 2) message(sprintf('# use %s genes', organism_))
      
      organism_code <- switch(organism_,
        mouse = 'mmu',
        human = 'hsa',
        stop(sprintf('ERROR: organism_ "%s" not recognized.', organism_))
      )
      
      SYMBOL2EG <- switch(organism_,
        mouse = org.Mm.egSYMBOL2EG,
        human = org.Hs.egSYMBOL2EG,
        stop(sprintf('ERROR: organism_ "%s" not recognized.', organism_))
      )
      
      selected_genes_filtered <- intersect(selected_genes_, mappedkeys(SYMBOL2EG))
      if (verbose_level_ >= 2) message(sprintf('# %d of %d selected genes are found in SYMBOL2EG', length(selected_genes_filtered), length(selected_genes_)))
      selected_genes_ <- unique(unlist(AnnotationDbi::as.list(SYMBOL2EG[selected_genes_filtered])))
      
      all_genes_filtered <- intersect(all_genes_, mappedkeys(SYMBOL2EG))
      if (verbose_level_ >= 2) message(sprintf('# %d of %d background genes are found in SYMBOL2EG', length(all_genes_filtered), length(all_genes_)))
      all_genes_ <- unique(unlist(AnnotationDbi::as.list(SYMBOL2EG[all_genes_filtered])))
    } else if (gene_id_type_ == 'entrez_id') {
    } else {
      stop(sprintf('ERROR: gene_id_type_ "%s" is not recognized', gene_id_type_))
    }
      
    
    # ego <- enrichKEGG(gene          = selected_genes_,
    #                   organism      = organism_code,
    #                   pvalueCutoff  = pvalueCutoff_,
    #                   universe      = all_genes_,
    #                   qvalueCutoff  = qvalueCutoff_, 
    #                   readable      = TRUE)
    # 
    # ego
  }




#' Gene Ontology (GO) term analysis
#' 
#' @param selected_genes_ A character vector of genes, by default gene symbols. The ID type can be changed
#'    using \code{gene_id_type_} parameter.
#' @param all_genes_ A character vector of genes treated as background. It should be a superset of \code{selected_genes_}.
#' @param organism_ Either "human" or "mouse".
#' @param gene_id_type_ Currently either "symbol" (default), or "entrez_id".
#' @param ontology_ Either "MF" (Molecular Function), "CC" (Cellular Component, default), or "BP" (Biological Process).
#' @param verbose_level_ 0 for quiet, 1 for simple messages, 2 for additional debug messages.
#' @param pvalueCutoff_,qvalueCutoff_ Cutoff threshold for either p or q values.
#' 
#' @import clusterProfiler
#' @import org.Mm.eg.db
#' @import org.Hs.eg.db
#' @import AnnotationDbi
#' @export
GO_analysis <-
  function(selected_genes_,
           all_genes_,
           organism_,
           gene_id_type_='symbol',
           ontology_='CC',
           verbose_level_=1,
           pvalueCutoff_=0.01,
           qvalueCutoff_=0.05
           ) {
    if (verbose_level_ >= 2) message('# [GO_analysis] is called')
    
    if (gene_id_type_ == 'symbol') {
      if (verbose_level_ >= 2) message('# genes are assumed to be symbols')
      
      if (verbose_level_ >= 2) message(sprintf('# use %s genes', organism_))
      
      OrgDb <- switch(organism_,
        mouse = 'org.Mm.eg.db',
        human = 'org.Hs.eg.db',
        stop(sprintf('ERROR: organism_ "%s" not recognized.', organism_))
      )
      
      SYMBOL2EG <- switch(organism_,
        mouse = org.Mm.egSYMBOL2EG,
        human = org.Hs.egSYMBOL2EG,
        stop(sprintf('ERROR: organism_ "%s" not recognized.', organism_))
      )
      
      selected_genes_filtered <- intersect(selected_genes_, mappedkeys(SYMBOL2EG))
      if (verbose_level_ >= 2) message(sprintf('# %d of %d selected genes are found in SYMBOL2EG', length(selected_genes_filtered), length(selected_genes_)))
      selected_genes_ <- unique(unlist(AnnotationDbi::as.list(SYMBOL2EG[selected_genes_filtered])))
      
      all_genes_filtered <- intersect(all_genes_, mappedkeys(SYMBOL2EG))
      if (verbose_level_ >= 2) message(sprintf('# %d of %d background genes are found in SYMBOL2EG', length(all_genes_filtered), length(all_genes_)))
      all_genes_ <- unique(unlist(AnnotationDbi::as.list(SYMBOL2EG[all_genes_filtered])))
    } else if (gene_id_type_ == 'entrez_id') {
    } else {
      stop(sprintf('ERROR: gene_id_type_ "%s" is not recognized', gene_id_type_))
    }
      
    
    ego <- enrichGO(gene          = selected_genes_,
                    OrgDb         = OrgDb,
                    ont           = ontology_,
                    pvalueCutoff  = pvalueCutoff_,
                    pAdjustMethod = "BH",
                    universe      = all_genes_,
                    qvalueCutoff  = qvalueCutoff_, 
                    readable      = TRUE)
    
    ego
  }



#' Generate modules from seed genes using Spinglass community detection method
#' 
#' @param expr_matrix_ The expression matrix, each row is a gene and each
#'   column is a sample. The row names should be gene symbols and the column
#'   names should be sample ids. Use the correct capticalization for gene
#'   symbols. For mouse genes,only the first letter is capitalized (e.g. Tp53);
#'   for human genes, all letters are capitalized (e.g. TP53). The values will
#'   be log-transformed by default, this can be changed using the
#'   \code{log_transofrmation_} option.
#' @param grouping_vec_ A factor vector, indicating the groupping. Its length
#'   should be the same as the number of columns in \code{expr_matrix_}. Example:
#'   \code{factor(c(1,1,1,2,2,2,2,2,2))}
#' @param selected_genes_ A character vector. The genes that are used as seed
#'   genes in the FEM algorithm. Typically generated from another method such as
#'   NMF.
#' @param ppi_edges_ A \code{data.frame} containing all the protein-protein interactions (PPI).
#'    It should have two columns of gene symbols, each row indicating a pair of genes that have interaction.
#'    \code{data(hppi)} and \code{data(mppi)} are two example PPI networks, respectively for humans and mice.
#' @param organism_ Either "human" or "mouse".
#' @param log_transformation_ Whether to log-transform the expression matrix.
#'   A pseudo count \code{1} will be added to each value before log transformation
#'   to avoid infinity. That is, \code{expr_matrix_ <- log(expr_matrix_+1)}. Default: \code{True}.
#' @param verbose_level_ How much information is printed. Either 0, 1 or 2,
#'   default: 1.
#' @param n_top_modules_ Number of top modules to pick. Default: 5.
#' @param n_threads_ Number of threads to use. Default: 1.
#' @param min_size_ Module size lower cut-off. Default: 10.
#' @param max_size_ Module size upper cut-off. Default: 100.
#' @param seed_ Random seed.
#' @param enrichment_analysis_ Either "GO" or "KEGG".
#' @param ... Additional parameters for \code{KEGG_analysis} or \code{GO_analysis}.
#' 
#' @import dplyr
#' @export
#' @examples
#' res <- spinglass_procedure(fpkm, phe, leading_genes, mppi, 'mouse', n_threads=50)
spinglass_procedure <-
  function(
           expr_matrix_,
           grouping_vec_,
           selected_genes_,
           ppi_edges_,
           organism_,
           log_transformation_ = T,
           verbose_level_ = 1,
           n_top_modules_ = 5, 
           n_threads_ = 1,
           min_size_ = 10,
           max_size_ = 100,
           enrichment_analysis_ = "GO",
           gamma_ = 0.5,
           seed_ = 12345,
           ...
           ) {
    options(stringsAsFactors=F)
    
    set.seed(seed_)
    
    output <- list()
    
    if (log_transformation_) {
      if (verbose_level_ >= 1) message('* Log transformation ...')
      expr_matrix_ <- log(expr_matrix_+1)
    }
    
    if (verbose_level_ >= 1) message('* Construct weighted graph from edge list ...')
    ppi <- construct_weighted_graph_from_edge_list(ppi_edges_, expr_matrix_, grouping_vec_)
    
    if (verbose_level_ >= 1) message('* Filter genes according to graph ...')
    genes <- rownames(expr_matrix_)
    
    genes <- intersect(genes, V(ppi)$name)  # 14865 -> 8247
    ppi <- induced_subgraph(ppi, genes)
    largest_cluster <- data_frame(membership=clusters(ppi)$membership) %>% group_by(membership) %>% summarize(n=n()) %>% arrange(-n) %>% .$membership %>% .[1]
    ppi <- induced.subgraph(ppi, V(ppi)$name[clusters(ppi)$membership==largest_cluster])
    
    output$ppi <- ppi
    
    genes <- V(ppi)$name
    expr_matrix_ <- expr_matrix_[genes, ]
    
    selected_genes_ <- intersect(selected_genes_, genes)
    
    if (verbose_level_ >= 2) message(sprintf('genes = %s', paste(genes, collapse=', ')))
    if (verbose_level_ >= 2) message(sprintf('ppi = '))
    if (verbose_level_ >= 2) print(ppi)
    if (verbose_level_ >= 2) message(sprintf('selected_genes_ (length: %d) = %s', length(selected_genes_), paste(selected_genes_, collapse=', ')))

    if (verbose_level_ >= 1) message('* Perform spinglass algorithm ...')
    spinglass_results <- generate_spinglass_communities(ppi, selected_genes_, 
                                                        min_size_ = min_size_, 
                                                        max_size_ = max_size_, 
                                                        n_threads_ = n_threads_, 
                                                        verbose_level_ = verbose_level_,
                                                        gamma_ = gamma_,
                                                        seed_ = seed_)
    
    output$spinglass_results <- spinglass_results
    
    if (verbose_level_ >= 1) message('* Perform Monte Carlo simulation to estimate p-values ...')
    seed_p_values <- get_monte_carlo_simulation_p_values(spinglass_results$graphs,
                                                         n_threads_ = n_threads_,
                                                         verbose_level_ = verbose_level_,
                                                         seed_=seed_)
    
    top_modules_tb <- spinglass_results$info %>%
      inner_join(seed_p_values, 'seed') %>%
      arrange(p_values) %>%
      head(n_top_modules_)
    
    output$top_modules_tb <- top_modules_tb
    
    top_modules <- top_modules_tb$seed
    
    output$top_modules <- top_modules
    
    if (verbose_level_ >= 1) message('* Perform enrichment analysis ...')
    top_module_genelists_ <- spinglass_results$genelists[top_modules]
    enrichment_df <- list()
    for (i in seq(top_module_genelists_)) {
      seed <- names(top_module_genelists_)[i]
      if (verbose_level_ >= 1) message(sprintf('  - Module with seed gene "%s" (%d/%d) ...', seed, i, length(top_module_genelists_)))
      
      if (enrichment_analysis_ == 'GO') {
        if (verbose_level_ >= 2) message("enrichment_analysis_ == 'GO'")
        enrichment_res <- GO_analysis(top_module_genelists_[[i]], genes, organism_, verbose_level_ = verbose_level_, ...) %>% summary
      } else if (enrichment_analysis_ == 'KEGG') {
        if (verbose_level_ >= 2) message("enrichment_analysis_ == 'KEGG'")
        enrichment_res <- KEGG_analysis(top_module_genelists_[[i]], genes, organism_, verbose_level_ = verbose_level_, ...) %>% summary
      } else {
        stop(sprintf('ERROR: enrichment_analysis_ "%s" not recognized', enrichment_analysis_))
      }
      enrichment_df$seed[i] <- seed
      enrichment_df$first_term[i] <- enrichment_res$Description[1]
      enrichment_df$first_qvalue[i] <- enrichment_res$qvalue[1]
      enrichment_df$second_term[i] <- enrichment_res$Description[2]
      enrichment_df$second_qvalue[i] <- enrichment_res$qvalue[2]
    }
    enrichment_df <- dplyr::as_data_frame(enrichment_df)
    
    output$final_tb <- top_modules_tb %>%
      inner_join(enrichment_df, 'seed') %>%
      arrange(p_values)
    
    output$graph_plot <- plot_module_graph(spinglass_results$graphs[top_modules],
                                           ppi, top_modules, levels(grouping_vec_)[1],
                                           levels(grouping_vec_)[2], seed_=seed_,
                                           verbose_level_=verbose_level_)

    return(output)
  }


#' Plot module graph
#' 
#' @import ggplot2
#' @import ggrepel
#' @import dplyr
#' @export
plot_module_graph <- 
  function(
    graphs_,
    ppi_,
    top_modules_,
    from_,
    to_,
    seed_=12345,
    verbose_level_=1
  ){
    set.seed(seed_)
    
    g <- graph.union(graphs_, byname=T)
    glayout <- layout_nicely(g)
    
    t_stats <- ppi_ %>%
      igraph::as_data_frame('vertices') %>%
      tbl_df %>%
      select(name, t_stat)
    
    if (verbose_level_ >= 2) message('gv = ')
    if (verbose_level_ >= 2) print(t_stats)
    
    gv <- igraph::as_data_frame(g, 'vertices') %>%
      tbl_df %>%
      select(name) %>%
      mutate(is_seed=factor(name %in% top_modules_, levels=c(T, F))) %>%
      mutate(x=glayout[,1]) %>%
      mutate(y=glayout[,2]) %>%
      left_join(t_stats, by='name')
    
    if (verbose_level_ >= 2) message('gv = ')
    if (verbose_level_ >= 2) print(gv)
    
    ge <- igraph::as_data_frame(g, 'edges') %>%
      tbl_df %>%
      select(from, to) %>%
      left_join(gv %>% select(from=name, from_x=x, from_y=y), 'from') %>%
      left_join(gv %>% select(to=name, to_x=x, to_y=y), 'to')
    
    if (verbose_level_ >= 2) message('ge = ')
    if (verbose_level_ >= 2) print(ge)
    
    graph_plot <- ggplot() +
      geom_segment(aes(x=from_x, y=from_y, xend=to_x, yend=to_y), ge, size=.5, color='#BEBEBE') +
      geom_point(aes(x, y, color=t_stat, size=is_seed), gv) +
      geom_text_repel(aes(x, y, label=name), gv) +
      scale_size_manual(name='', values=c(6, 3), labels=c('seed', 'not seed')) +
      scale_color_gradient2(limits=c(-2.5, 2.5), name=sprintf('t-stat for\nexpression level\nchange\nfrom %s\nto %s', from_, to_),
                            low='blue', mid='#BEBEBE', high='red') +
      theme_void() + 
      theme(legend.title=element_text(size=10))
    
    graph_plot
  }

#' Detect subpopulations using NMF

#' @param expr_matrix_ The expression matrix, each row is a gene and each
#'   column is a sample. The row names should be gene symbols and the column
#'   names should be sample ids. Use the correct capticalization for gene
#'   symbols. For mouse genes,only the first letter is capitalized (e.g. Tp53);
#'   for human genes, all letters are capitalized (e.g. TP53). The values will
#'   be log-transformed by default, this can be changed using the
#'   \code{log_transofrmation_} option.
#' @param log_transformation_ Whether to log-transform the expression matrix.
#'   A pseudo count \code{1} will be added to each value before log transformation
#'   to avoid infinity. That is, \code{expr_matrix_ <- log(expr_matrix_+1)}. Default: \code{True}.
#' @param verbose_level_ How much information is printed. 0 = quiet, 1 = normal, 2 = with debug info, 3 = with extra debug info. Default: 1.
#' @param n_threads_,nrun_,method_,.options_,seed_ Parameters for \code{nmf(...)}.
#'   Default values are supplied.
#' @return A list with the following elements
#'   \describe{
#'     \item{\code{nmf_result}}{The \code{NMFfitX1} object returned by the \code{nmf} call.}
#'     \item{\code{gene_info}}{A \code{data_frame} detailing the bases, weighted bases, 
#'        and D-score for each gene. The data frame is sorted by their D-scores}
#'     \item{\code{d_score_frequency_plot}}{Kernal distribution for the D-scores, separated by the coef components they represent.}
#'     \item{\code{ordered_sample_ids}}{The pseudo-order of samples calculated by sorting the differences of coef values for each sample.} 
#'     \item{\code{coef_line_dat}}{Plotting data for \code{coef_line_plot}}
#'     \item{\code{coef_line_plot}}{A line plot showing the trend formed by the coef values, using the order in \code{ordered_sample_ids}.}
#'     \item{\code{path_dat}}{Plotting data for \code{pca_path_plot} and \code{mds_path_plot}.}
#'     \item{\code{pca_path_plot}}{Connecting the PCA plot of the expression matrix using the order in \code{ordered_sample_ids}.
#'        The thickness of the lines indicate the jump of differences of the coef values.} 
#'     \item{\code{mds_path_plot}}{Same as \code{pca_path_plot}, except for MDS plot. The distance is "1 - Pearson Correlation".}
#'     \item{\code{top_genes_dat}}{Plotting data for \code{top_genes_free_y_plot} and \code{top_genes_fixed}.}
#'     \item{\code{top_genes_free_y_plot}}{Faceted bar-plot showing the expression levels of the top 100 genes, across all samples.
#'        Samples are sorted using the order in \code{ordered_sample_ids}.}
#'     \item{\code{top_genes_fixed}}{Sample as \code{top_genes_free_y_plot}, except all y-axes are fixed across all faceted.}
#'   }
#' @import NMF
#' @import ggplot2
#' @import reshape2
#' @import tidyr
#' @import dplyr
#' @export
#' @examples
#' res <- nmf_subpopulation(fpkm)
nmf_subpopulation <-
  function(
           expr_matrix_,
           log_transformation_=T,
           verbose_level_=1,
           n_threads_=1,
           nrun_=30,
           method_="brunet",
           top_genes_facet_title_font_size_=24,
           .options_=sprintf("p%dv%d", n_threads_, verbose_level_-1),
           seed_=12345
           ) {
    options(stringsAsFactors=F)
    
    if (class(expr_matrix_) != 'matrix') stop('Error: expr_matrix_ has to be of class \'matrix\'. Stop.')
    
    library(NMF)
    
    output <- list()
    
    non_zero_rows_tf <- apply(expr_matrix_, 1, function(x)any(x>0))
    if (!all(non_zero_rows_tf)) {
      if (verbose_level_ >= 1) message(sprintf('* All-zero rows detected, getting rid of %d/%d genes ...',
                                               sum(!non_zero_rows_tf), length(non_zero_rows_tf)))
      if (verbose_level_ >= 2) message('removed genes = ')
      if (verbose_level_ >= 2) print(rownames(expr_matrix_)[!non_zero_rows_tf])
      expr_matrix_ <- expr_matrix_[non_zero_rows_tf,]
    }
    
    if (log_transformation_) {
      if (verbose_level_ >= 1) message('* Log transformation ...')
      rwl <- log(expr_matrix_+1)
    } else {
      rwl <- expr_matrix_
    }

    rwl <- rbind(rwl, artificial=mean(rwl))
    
#   -----------------------------------------------------------------------
        
    if (verbose_level_ >= 1) message('* Call nmf function ...')
    ren <- nmf(rwl, rank=2, nrun=nrun_, method=method_, .options=.options_, seed=seed_)
    output$nmf_result <- ren
    
#   -----------------------------------------------------------------------
    
    if (verbose_level_ >= 1) message('* Calculating D-scores ...')
    weighting_factors <- basis(ren)['artificial',]
    
    d_score_calc <- function(x) {
      names(x) <- paste0('basis_',1:length(x))
      weighted_x <- x/weighting_factors
      names(weighted_x) <- paste0(names(x), '_weighted')
      basis_sorted <- sort(range(weighted_x), decreasing=T)
      d_score <- basis_sorted[1] - basis_sorted[2]
      which_max <- which.max(unname(weighted_x))
      c(x, weighted_x, d_score=d_score, which_max=which_max)
    }
    
    gene_info <- basis(ren) %>% apply(1, d_score_calc) %>%
      t %>% as.data.frame %>% add_rownames('gene_name') %>%
      tbl_df %>% mutate(which_max=factor(which_max)) %>%
      arrange(-d_score) %>%
      mutate(rank=row_number())
    
    table(gene_info$which_max)
    d_score_frequency_plot <- ggplot() +
      geom_density(aes(x=d_score, color=which_max, fill=which_max), gene_info, alpha=0.3, stat='bin') +
      scale_x_log10() +
      scale_fill_manual(name='Expr. Pattern', values=c('blue', 'red', 'green')) +
      scale_color_manual(name='Expr. Pattern', values=c('blue', 'red', 'green')) +
      theme(plot.margin=rep(unit(0.2, 'inches'), 4)) +
      labs(x='D-score', y='Frequency')
    
    output$gene_info <- gene_info
    output$d_score_frequency_plot <- d_score_frequency_plot
    
    if (verbose_level_ >= 2) message('gene_info = ')
    if (verbose_level_ >= 2) print(gene_info, n=20)
    
#   -----------------------------------------------------------------------
    
    if (verbose_level_ >= 1) message('* Plot the coef lines ...')
    
    subp_order <- order(coef(ren)[2,] - coef(ren)[1,])
    ordered_sample_ids <- colnames(rwl)[subp_order]
    
    output$ordered_sample_ids <- ordered_sample_ids
    
    dat_coef_line <- melt(coef(ren)[,ordered_sample_ids], c('coef', 'sample')) %>% tbl_df
    
    coef_line_plot <- ggplot() +
      geom_line(aes(x=sample, group=factor(coef), color=factor(coef), y=value), dat_coef_line) +
      geom_point(aes(x=sample, group=factor(coef), color=factor(coef), y=value), dat_coef_line) +
      scale_color_manual(name='Expr. Pattern', values=c('blue', 'red', 'green')) +
      theme(axis.text.x=element_text(angle=90, vjust=0.5), plot.margin=rep(unit(0.2, 'inches'), 4)) +
      labs(x='samples sorted by pseudo-timeline', y='log expr. lvl.')
    
    output$coef_line_dat <- dat_coef_line
    output$coef_line_plot <- coef_line_plot
    
#   -----------------------------------------------------------------------
    
    if (verbose_level_ >= 1) message('* Plot the trend paths on MDS and PCA ...')
    
    mdsr <- cmdscale(1-cor(rwl))
    pr <- prcomp(t(rwl))
    
    dat_pca <- dat_coef_line %>%
      mutate(sample=factor(sample, ordered_sample_ids)) %>%
      spread(coef, value) %>%
      mutate(difference=`1`-`2`) %>%
      mutate(blend=`1`/(`1`+`2`)) %>%
      mutate(pc1=pr$x[as.character(sample),1]) %>%
      mutate(pc2=pr$x[as.character(sample),2]) %>%
      mutate(mds1=mdsr[as.character(sample),1]) %>%
      mutate(mds2=mdsr[as.character(sample),2])
    
    if (verbose_level_ >= 2) message('dat_pca = ')
    if (verbose_level_ >= 2) print(dat_pca, 20)
    
    dat_seg <- dat_pca %>%
      mutate(order=row_number()) %>%
      mutate(prev_blend=lag(blend)) %>%
      mutate(avg_blend=(blend + prev_blend)/2) %>%
      mutate(prev_sample=lag(sample)) %>%
      mutate(prev_difference=lag(difference)) %>%
      mutate(prev_pc1=lag(pc1)) %>%
      mutate(prev_pc2=lag(pc2)) %>%
      mutate(prev_mds1=lag(mds1)) %>%
      mutate(prev_mds2=lag(mds2)) %>%
      mutate(jump=difference - prev_difference)
    
    if (verbose_level_ >= 2) message('dat_seg = ')
    if (verbose_level_ >= 2) print(dat_seg, 20)
    
    output$path_dat <- dat_seg
    
    pca_path_plot <- ggplot() +
      geom_segment(aes(x=prev_pc1, y=prev_pc2, xend=pc1, yend=pc2, size=jump, alpha=jump, color=avg_blend), dat_seg) +
      geom_point(aes(x=pc1, y=pc2), dat_pca) +
      scale_color_gradientn(name='precentage\nof the blue\ncomponent', colors=c('red', 'blue')) +
      scale_alpha_continuous(name='Jump', range=c(0.1,0.6)) +
      scale_size_continuous(name='Jump', range=c(0.5,2)) +
      labs(x='PC1', y='PC2')
    
    output$pca_path_plot <- pca_path_plot
    
    mds_path_plot <- ggplot() +
      geom_segment(aes(x=prev_mds1, y=prev_mds2, xend=mds1, yend=mds2, size=jump, alpha=jump, color=avg_blend), dat_seg) +
      geom_point(aes(x=mds1, y=mds2), dat_pca) +
      scale_color_gradientn(name='precentage\nof the blue\ncomponent', colors=c('red', 'blue')) +
      scale_alpha_continuous(name='Jump', range=c(0.1,0.6)) +
      scale_size_continuous(name='Jump', range=c(0.5,2)) +
      labs(x='MDS1', y='MDS2')
    
    output$mds_path_plot <- mds_path_plot
    
#   -----------------------------------------------------------------------

    if (verbose_level_ >= 1) message('* Plot the expression bar-chart of the top genes ...')
    
    top_genes <- gene_info$gene_name[1:min(100, nrow(gene_info))]
    
    top_genes_dat <- expand.grid(gene_name=top_genes, sample=ordered_sample_ids) %>% tbl_df %>%
      mutate(gene_name=as.character(gene_name)) %>% rowwise %>%
      mutate(log_expr=rwl[as.character(gene_name),as.character(sample)]) %>% ungroup %>%
      left_join(gene_info %>% select(gene_name, which_max, d_score, rank), by='gene_name') %>%
      mutate(gene_name=factor(gene_name, top_genes)) %>%
      mutate(sample=factor(sample, ordered_sample_ids))
    
    output$top_genes_dat <- top_genes_dat
    
    if (verbose_level_ >= 2) message('top_genes_dat = ')
    if (verbose_level_ >= 2) print(top_genes_dat, n=20)
    
    gglabeller <- function(df){
      res <- df %>% left_join(gene_info, by='gene_name') %>% rowwise %>%
        mutate(label=sprintf('%d - %s', rank, gene_name)) %>%
        select(label)
    }
    
    top_genes_free_y_plot <- ggplot() +
      geom_bar(aes(x=sample, y=log_expr, fill=which_max), top_genes_dat, stat='identity') +
      scale_fill_manual(values=c('blue', 'red', 'green')) +
      facet_wrap(~ gene_name, labeller=gglabeller, scales='free') +
      theme(axis.text.x=element_text(angle=90, vjust=0.5)) +
      theme(strip.text=element_text(size=top_genes_facet_title_font_size_), panel.margin = unit(2, "lines"), legend.position='none')
    
    output$top_genes_free_y_plot <- top_genes_free_y_plot
    
    top_genes_fixed <- ggplot() +
      geom_bar(aes(x=sample, y=log_expr, fill=which_max), top_genes_dat, stat='identity') +
      scale_fill_manual(values=c('blue', 'red', 'green')) +
      facet_wrap(~ gene_name, labeller=gglabeller) +
      theme(axis.text.x=element_text(angle=90, vjust=0.5)) +
      theme(strip.text=element_text(size=top_genes_facet_title_font_size_), panel.margin = unit(2, "lines"), legend.position='none')
    
    output$top_genes_fixed <- top_genes_fixed
    
#   -----------------------------------------------------------------------
    
    output
  }
