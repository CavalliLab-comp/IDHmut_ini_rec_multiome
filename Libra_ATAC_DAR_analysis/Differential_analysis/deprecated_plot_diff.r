#!/bioinfo/users/tberthom/miniforge3/envs/R-ArchR3/bin/R

library(ArchR)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggbeeswarm)

CODE_PATH <- '/mnt/beegfs/home/alyne/persistent/ArchR_functions_Thomas'

source(file.path(CODE_PATH, "Differential_analysis", "Tf_related.r"))
source(file.path(CODE_PATH, "Differential_analysis", "plot_diff.r"))


# from Summarized experiment object of fullmarkeranalysis
# plot top tf families of each category
plotTopFamilies <- function(se, threshold1=2, threshold2=20){
    gene_df <- extractTFperConditionDF(se, "mlog10Padj", threshold1)
    filtered_df <- gene_df %>% filter(.data[["mlog10Padj"]] >= threshold2)

    family_levels <- unique(c(filtered_df$family))
    nb.cols <- length(family_levels)
    color_palette <- viridis::turbo(nb.cols)
    p <- list()

    # Full plot
    p$full_plot <- ggplot(gene_df, aes(x = factor(Condition), y = mlog10Padj, color = family)) +
      geom_quasirandom(size = 3, alpha = 0.6) +  
      scale_color_manual(values = setNames(color_palette, family_levels)) +  # Set color scale
      theme_minimal() +
      labs(x = "Category Group", y = "mlog10Padj", title = "Swarm Plot of families, all TF") +
      theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust = 1)) 

    # Top plot
    p$top_plot <- ggplot(filtered_df, aes(x = factor(Condition), y = mlog10Padj, color = family)) +
      geom_quasirandom(size = 3, alpha = 0.6) +  
      scale_color_manual(values = setNames(color_palette, family_levels)) +  # Set color scale
      theme_minimal() +
      labs(x = "Category Group", y = "mlog10Padj", title = "Swarm Plot of families, top TF") +
      theme(axis.text.x = element_text(angle = 90, hjust = 1))
    return(p)
}



# from one Summarized Experiment list output of fullDifferential analysis,
# plot a 2 by 2 plots of all 4 objects with a title
plotSingleDifferentialAnalysis <- function(se_list,
                              annotation_control = NULL,
                              label_control = NULL,
                              threshold_FDR_regions = 0.1,
                              threshold_LogFC_regions = 0.5,
                              threshold_FDR_post = 0.05,
                              threshold_LogFC_post = 1, 
                              annotation_differential = NULL,
                              Group1 = NULL,
                              Group2 = NULL
                              ) {
  if (is.null(annotation_differential)) {
    annotation_differential <- metadata(se_list$Genes)$Params$groupBy
  }
  if (is.null(Group1)) {
    Group1 <- metadata(se_list$Genes)$Params$useGroups
  }
  if (is.null(Group2)) {
    Group2 <- metadata(se_list$Genes)$Params$bgdGroups
  }

  # Increase plot size
  options(repr.plot.width = 14, repr.plot.height = 12)
  # plotting regions
  threshold_regions = paste('FDR <=', threshold_FDR_regions, '& abs(Log2FC) >=', threshold_LogFC_regions, sep = " ")
  pv_regions <- plotMarkers(seMarker = se_list$Regions, name = colnames(se_list$Regions), cutOff = threshold_regions, plotAs = "Volcano")

  current_title <- ggplot_build(pv_regions)$plot$labels$title
  new_title <- paste("Differential peak accessibility of ", Group1, " vs ", Group2, sep = "")
  pv_regions_with_title <- pv_regions + labs(title = new_title, subtitle = current_title)
  
  # plotting gene scores 
  threshold_post_processing = paste('FDR <=', threshold_FDR_post, '& abs(Log2FC) >=', threshold_LogFC_post, sep = " ")
  if (class(se_list$Genes) %in%  c("data.frame", "tibble")) { # plot with Libra output
    pv_genes_with_title <- plotVolcanoLibra(as.data.frame(se_list$Genes), logfc_col = "avg_logFC", pval_col = "p_val_adj",
                         label_col = "gene", title = paste("Differential gene accessibility score of ", Group1, " vs ", Group2, sep = ""),
                         logfc_thresh = threshold_LogFC_post, pval_thresh = threshold_FDR_post)

    # Combine regions and genes plots 
    combined_regions_genes <- (pv_regions_with_title | pv_genes_with_title) +
    plot_layout(widths = c(1, 1)) & theme(legend.position = "bottom")
  } else { # plot with ArchR output
    pv_genes <- plotMarkers(seMarker = se_list$Genes, name = colnames(se_list$Genes), cutOff = threshold_post_processing, plotAs = "Volcano")
    current_title <- ggplot_build(pv_genes)$plot$labels$title
    new_title <- paste("Differential gene accessibility score of ", Group1, " vs ", Group2, sep = "")
    pv_genes_with_title <- pv_genes + labs(title = new_title, subtitle = current_title)

    # Combine regions and genes plots with a shared legend at the bottom
    combined_regions_genes <- (pv_regions_with_title | pv_genes_with_title) +
    plot_layout(guides = "collect", widths = c(1, 1)) & theme(legend.position = "bottom")
  }

  # plotting enriched motifs 
  # Compute global color scale limits
  all_values <- c(assay(se_list$MotifsDown)[,1], assay(se_list$MotifsUp)[,1])
  color_limits <- range(all_values)

  # Generate plots with shared color scale
  up <- plotTopDiffMotifs(se_list$MotifsUp, name = Group1, color_limits = color_limits)
  down <- plotTopDiffMotifs(se_list$MotifsDown, name = Group2, color_limits = color_limits)

  # Ensure both plots share the same y-axis scale
  y_min <- min(down$data$mlog10Padj, up$data$mlog10Padj)
  y_max <- max(down$data$mlog10Padj, up$data$mlog10Padj)

  up <- up + coord_cartesian(ylim = c(y_min, y_max))
  down <- down + coord_cartesian(ylim = c(y_min, y_max))

  # Combine plots with shared legend and put legend at bottom
  combined_motifs <- up + down + plot_layout(guides = "collect", widths = c(1, 1)) & theme(legend.position = "bottom")

  # build title 
  if (is.null(annotation_control)){
    TitlePlot <- paste("Differential analysis on", annotation_differential, Group1, "vs", Group2)
  } else {
    TitlePlot <- paste("Differential analysis on", annotation_differential, Group1, "vs", Group2, "for subgroup", annotation_control,"=", label_control)
  }

  # Final 2x2 grid
  final_plot <- (wrap_elements(grid::textGrob(TitlePlot, gp = grid::gpar(fontsize = 16, fontface = "bold"))) /
                 combined_regions_genes /
                 combined_motifs) +
                 plot_layout(ncol = 1, heights = c(0.1, 1, 1))
  return(final_plot)
}

#do previous plot on paiwise diff analysis ouput
plotPairwiseDifferentialAnalysis <- function(se_list_pairwise, 
                                             annotation_control = NULL, 
                                             label_control = NULL,
                                             threshold_FDR_regions = 0.1,
                                             threshold_LogFC_regions = 0.5,
                                             threshold_FDR_post = 0.05,
                                             threshold_LogFC_post = 1
                                            ) {
    list_out <- lapply(names(se_list_pairwise), function(name) {
        plotSingleDifferentialAnalysis(se_list = se_list_pairwise[[name]], 
                                  annotation_control = annotation_control,
                                  label_control = label_control,
                                  threshold_FDR_regions = threshold_FDR_regions,
                                  threshold_LogFC_regions = threshold_LogFC_regions,
                                  threshold_FDR_post = threshold_FDR_post,
                                  threshold_LogFC_post = threshold_LogFC_post
                                  ) 
    })
    names(list_out) <- names(se_list_pairwise)
    return(list_out)
}

# finally do previous plot on full diff analysis output
# = across subgroups 
plotFullDifferentialAnalysisSub <- function(output_diffanalysis_sub, 
                                             annotation_control, 
                                             threshold_FDR_regions = 0.1, #should match the one used for getmarkerfeatures
                                             threshold_LogFC_regions = 0.5,
                                             threshold_FDR_post = 0.05,
                                             threshold_LogFC_post = 1
                                            ) {
    list_out <- lapply(names(output_diffanalysis_sub), function(name) {
        plotPairwiseDifferentialAnalysis(se_list = output_diffanalysis_sub[[name]], 
                                         annotation_control = annotation_control,
                                          label_control = name,
                                  threshold_FDR_regions = threshold_FDR_regions,
                                  threshold_LogFC_regions = threshold_LogFC_regions,
                                  threshold_FDR_post = threshold_FDR_post,
                                  threshold_LogFC_post = threshold_LogFC_post
                                  )
    })
    names(list_out) <- names(output_diffanalysis_sub)
    return(list_out)
}