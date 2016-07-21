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
    if (verbose_level_ >= 2) message('num of seed genes = ', length(seed_genes_))
    
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
#' @param enrichment_analysis_ Either "GO" or "KEGG".
#' @param gamma_ Parameter for \code{igraph}'s \code{cluster_spinglass} function.
#' @param seed_ Random seed.
#' @param ... Additional parameters for \code{KEGG_analysis} or \code{GO_analysis}.
#' 
#' @import dplyr
#' @export
#' @examples
#' res <- spinglass_procedure(fpkm, phe, leading_genes, mppi, 'mouse')
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
    ppi <- induced.subgraph(ppi, V(ppi)$name[clusters(ppi)$membership==1])  # 8247 -> 8174
    
    output$ppi <- ppi
    
    genes <- V(ppi)$name
    expr_matrix_ <- expr_matrix_[genes, ]
    
    selected_genes_ <- intersect(selected_genes_, genes)
    
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

