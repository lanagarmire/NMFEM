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
#' @param method_ Interative method for NMF. See the documentation of \code{NMF::nmf}.
#' @param seed_ Random seed.
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
#' res <- nmf_subpopulation(toy)
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
      mutate(rank=row_number()) %>%
      filter(gene_name != 'artificial')
    
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
