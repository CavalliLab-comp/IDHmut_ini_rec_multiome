#!/bioinfo/users/tberthom/miniforge3/envs/R-ArchR3/bin/R
library(ggplot2)
library(dplyr)
library(ArchR)
library(SummarizedExperiment)
# Author : Thomas Berthomier

CODE_PATH <- '/mnt/beegfs/home/alyne/persistent/ArchR_functions_Thomas'

source(file.path(CODE_PATH, "Miscellaneous", "SE_operations.r"))

#plot top motifs from differential enriched motifs archr output
plotTopDiffMotifs <- function(motif_se, name = NULL, n_top_labels = 20, color_limits = NULL) {
    df <- data.frame(TF = rownames(motif_se), mlog10Padj = assay(motif_se)[, 1])
    df <- df[order(df$mlog10Padj, decreasing = TRUE), ]
    df$rank <- seq_len(nrow(df))

    ## plot
    p <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) +
        geom_point(size = 1) +
        ggrepel::geom_label_repel(
            data = df[rev(seq_len(n_top_labels)), ], aes(x = rank, y = mlog10Padj, label = TF),
            size = 1.5,
            nudge_x = 2,
            color = "black",
            max.overlaps = 30
        ) +
        theme_ArchR() +
        ylab("-log10(P-adj) Motif Enrichment") +
        xlab("Rank Sorted TFs Enriched") +
        scale_color_gradientn(colors = paletteContinuous(set = "comet"), limits = color_limits) +
        theme(legend.position = "bottom") # Legend at the bottom

    if (!is.null(name)) {
        p <- p + ggtitle(paste("Top enriched motifs in", name, sep = " "))
    }

    return(p)
}


#plot proportion sig peaks per chromosome
plotSigPeaksPerChr <- function(se_regions,
                               threshold_FDR = 0.1,
                               threshold_Log2FC = 0.5) {
  reg_df <- getDFfromSESample(se_regions)

  df_summary <- reg_df %>%
    mutate(
      significance = case_when(
        FDR < threshold_FDR & Log2FC > threshold_Log2FC ~ "up",
        FDR < threshold_FDR & Log2FC < -threshold_Log2FC ~ "down",
        TRUE ~ "not significant"
      )
    ) %>%
    group_by(seqnames, significance) %>%
    summarise(count = n(), .groups = "drop") %>%
    tidyr::pivot_wider(
      names_from = significance,
      values_from = count,
      values_fill = list(count = 0)
    ) %>%
    # Ensure all expected columns exist
    mutate(
      up = if (!"up" %in% names(.)) 0 else up,
      down = if (!"down" %in% names(.)) 0 else down,
      `not significant` = if (!"not significant" %in% names(.)) 0 else `not significant`
    )

  df_long <- df_summary %>%
    mutate(total = rowSums(across(where(is.numeric)))) %>%
    pivot_longer(
      cols = c("up", "down", "not significant"),
      names_to = "Direction",
      values_to = "Count"
    ) %>%
    mutate(
      Proportion = Count / total,
      seqnames = factor(
        seqnames,
        levels = unique(seqnames[order(as.numeric(gsub("chr", "", seqnames)))])),
      Direction = factor(Direction, levels = c("up", "not significant", "down"))
    )

  ggplot(df_long, aes(x = seqnames, y = Proportion, fill = Direction)) +
    geom_bar(stat = "identity") +
    geom_text(
      aes(label = Count),
      position = position_stack(vjust = 0.5),
      angle = 90,
      hjust = 0.5,
      vjust = 0.5,
      size = 3,
      color = "black"
    ) +
    scale_y_continuous(labels = scales::percent_format()) +
    scale_fill_manual(values = c("up" = "firebrick", "down" = "steelblue", "not significant" = "gray80")) +
    labs(
      title = "Proportion of Peaks by Direction per Chromosome",
      x = "Chromosome",
      y = "Proportion of Peaks"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    )
}


#' Create Volcano Plot for scATAC-seq Differential Accessibility Results
#'
#' Generates a volcano plot to visualize differential accessibility results from
#' Libra pseudobulk analysis with customizable significance thresholds and automatic
#' labeling of top differential peaks.
#'
#' @param df Data frame containing differential accessibility results from Libra analysis
#' @param logfc_col String specifying column name containing log fold change values.
#'   Default: "avg_logFC"
#' @param pval_col String specifying column name containing adjusted p-values.
#'   Default: "p_val_adj"
#' @param label_col String specifying column name containing feature labels for
#'   annotation (genomic coordinates). Default: "gene"
#' @param title String specifying plot title. Default: "Volcano Plot"
#' @param logfc_thresh Numeric threshold for log fold change significance.
#'   Default: 0.25
#' @param pval_thresh Numeric threshold for adjusted p-value significance.
#'   Default: 0.05
#'
#' @return ggplot2 object representing the volcano plot with differential accessibility
#'   visualization
#'
#' @details The function creates a scatter plot optimized for scATAC-seq differential
#'   accessibility analysis with:
#' \itemize{
#'   \item X-axis: Log2 fold change values (accessibility difference)
#'   \item Y-axis: -log10 transformed adjusted p-values
#'   \item Color coding: More accessible (red), Less accessible (blue), 
#'     Not significant (grey)
#'   \item Threshold lines: Dashed lines indicating significance cutoffs
#'   \item Labels: Top 5 most significant differential peaks meeting both thresholds
#'   \item Professional styling with ggrepel for non-overlapping labels
#' }
#'
#' Points are classified as:
#' \itemize{
#'   \item Upregulated (More Accessible): p_val < pval_thresh AND logFC > logfc_thresh
#'   \item Downregulated (Less Accessible): p_val < pval_thresh AND logFC < -logfc_thresh
#'   \item Not significant: All other points
#' }
#'
#' @export
plotVolcanoLibra <- function(df, logfc_col = "avg_logFC", pval_col = "p_val_adj",
                         label_col = "gene", title = "Volcano Plot",
                         logfc_thresh = 0.25, pval_thresh = 0.05) {
  #think the log fold changes are the wrong direction so reverse
  df[[logfc_col]] = -df[[logfc_col]]
  
  # Add significance category
  df <- df %>%
    mutate(
      significance = case_when(
        !!sym(pval_col) < pval_thresh & !!sym(logfc_col) > logfc_thresh ~ "Upregulated",
        !!sym(pval_col) < pval_thresh & !!sym(logfc_col) < -logfc_thresh ~ "Downregulated",
        TRUE ~ "Not significant"
      )
    )
  
  # Base plot
  p <- ggplot(df, aes(x = !!sym(logfc_col), y = -log10(!!sym(pval_col)))) +
    geom_point(aes(color = significance), alpha = 0.8, size = 1) +
    scale_color_manual(values = c("Upregulated" = "#E41A1C", "Downregulated" = "#377EB8", "Not significant" = "grey70")) +
    theme_minimal(base_size = 14) +
    labs(
      title = title,
      x = "Log2 Fold Change",
      y = "-Log10 Adjusted P-Value",
      color = "Significance"
    ) +
    geom_vline(xintercept = c(-logfc_thresh, logfc_thresh), linetype = "dashed", color = "grey50") +
    geom_hline(yintercept = -log10(pval_thresh), linetype = "dashed", color = "grey50")
  
  # # Optional: Label top significant points
  # top_genes <- df %>%
  #   filter(!!sym(pval_col) < pval_thresh & abs(!!sym(logfc_col)) > logfc_thresh) %>%
  #   arrange(!!sym(pval_col)) %>%
  #   slice_head(n = 5)
  # 
  # p <- p + ggrepel::geom_text_repel(
  #   data = top_genes,
  #   aes(label = !!sym(label_col)),
  #   size = 3.5,
  #   max.overlaps = 20
  # )
  
  return(p)
}