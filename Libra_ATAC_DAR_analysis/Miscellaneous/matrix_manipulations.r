#!/bioinfo/users/tberthom/miniforge3/envs/R-ArchR4/bin/R
library(cluster)


#' Clip Outlier Values from Numeric Data
#'
#' This function clips outlier values in numeric columns of a matrix or data frame
#' using one of several methods: IQR-based clipping, quantile-based proportion clipping,
#' or custom threshold clipping.
#'
#' @param data A numeric matrix or data frame. Must have column names.
#' @param method Character string specifying the clipping method. One of:
#'   \itemize{
#'     \item "iqr" - Uses interquartile range method
#'     \item "proportion" - Uses quantile-based proportion method (default)
#'     \item "custom" - Uses user-specified min/max values
#'   }
#' @param lower_quantile Numeric value between 0 and 1. Lower quantile for IQR method (default: 0.25)
#' @param upper_quantile Numeric value between 0 and 1. Upper quantile for IQR method (default: 0.75)
#' @param iqr_multiplier Numeric multiplier for IQR method (default: 1.5)
#' @param clip_proportion Numeric value between 0 and 1. Proportion of extreme values to clip
#'   from each tail for proportion method (default: 0.01, clips 1% from each tail)
#' @param custom_min Numeric value for custom lower bound (only used with method = "custom")
#' @param custom_max Numeric value for custom upper bound (only used with method = "custom")
#'
#' @return A matrix with outliers clipped to the calculated bounds
#'
#' @details
#' The function applies the selected clipping method to each numeric column:
#' \itemize{
#'   \item IQR method: Values below Q1 - k*IQR or above Q3 + k*IQR are clipped
#'   \item Proportion method: Values below the (p/2)th percentile or above the (1-p/2)th percentile are clipped
#'   \item Custom method: Values outside the specified min/max range are clipped
#' }
#'
#' Non-numeric columns are skipped. Missing values (NA) are preserved.
#'
#' @examples
#' # Create sample data with outliers
#' set.seed(123)
#' data <- matrix(c(rnorm(100, 50, 10), c(0, 100, rnorm(98, 50, 10))), ncol = 2)
#' colnames(data) <- c("var1", "var2")
#'
#' # Clip using proportion method (default)
#' clipped1 <- clipOutliers(data, method = "proportion", clip_proportion = 0.05)
#'
#' # Clip using IQR method
#' clipped2 <- clipOutliers(data, method = "iqr", iqr_multiplier = 2)
#'
#' # Clip using custom bounds
#' clipped3 <- clipOutliers(data, method = "custom", custom_min = 20, custom_max = 80)
#' @export
clipOutliers <- function(data, method = "proportion", # 'iqr', 'custom', or 'proportion'
                          lower_quantile = 0.25, upper_quantile = 0.75,
                          iqr_multiplier = 1.5,
                          clip_proportion = 0.01,
                          custom_min = NULL, custom_max = NULL) {
  
  # Convert to data.frame if it's a matrix
  if (is.matrix(data)) {
    data <- as.data.frame(data)
  }
  
  if (is.null(colnames(data))) {
    stop("Input data must have column names.")
  }
  
  for (col in colnames(data)) {
    x <- data[[col]]
    
    if (!is.numeric(x)) next  # Skip non-numeric columns
    if (method == "iqr") {
      q1 <- quantile(x, lower_quantile, na.rm = TRUE)
      q3 <- quantile(x, upper_quantile, na.rm = TRUE)
      iqr <- q3 - q1
      lower <- q1 - iqr_multiplier * iqr
      upper <- q3 + iqr_multiplier * iqr
    } else if (method == "custom") {
      lower <- if (!is.null(custom_min)) custom_min else min(x, na.rm = TRUE)
      upper <- if (!is.null(custom_max)) custom_max else max(x, na.rm = TRUE)
    } else if (method == "proportion") {
      lower <- quantile(x, clip_proportion / 2, na.rm = TRUE)
      upper <- quantile(x, 1 - clip_proportion / 2, na.rm = TRUE)
    } else {
      stop("Method must be one of: 'iqr', 'custom', or 'proportion'.")
    }
    
    x[x < lower] <- lower
    x[x > upper] <- upper
    data[[col]] <- x
  }
  
  return(as.matrix(data))
}

#' Compare Clustering Quality Between Two Features Using Silhouette Analysis
#'
#' This function performs silhouette analysis to compare how well two different features 
#' separate samples in RNA-seq data.
#'
#' @param expr_data A numeric matrix or data frame of RNA-seq expression data.
#'                 ROWS should be genes/features and COLUMNS should be samples.
#'                 Example dimensions: 20,000 genes (rows) x 50 samples (columns).
#'                 Should already be normalized and transformed (e.g., log2(counts+1)).
#'
#' @param feature1 A named character or factor vector containing cluster assignments for the first feature.
#'                The names should match the column names of expr_data (sample IDs).
#'                Example: c("tumor"="sample1", "normal"="sample2", "tumor"="sample3", ...)
#'                Alternative format: A character vector with length equal to ncol(expr_data) where
#'                the positions correspond to columns in expr_data.
#'
#' @param feature2 A named character or factor vector containing cluster assignments for the second feature.
#'                Same format requirements as feature1.
#'
#' @param feature1_name Character string naming the first feature (default: "Feature A").
#'
#' @param feature2_name Character string naming the second feature (default: "Feature B").
#'
#' @param dist_method Character string specifying distance metric to use (default: "euclidean").
#'                   Other options: "correlation", "manhattan", etc. See ?dist for all options.
#'
#' @return A list containing silhouette objects, average silhouette widths, p-value from
#'         Wilcoxon test, and indication of which feature had better clustering.
#'
#' @examples
#' # With expression data as genes (rows) x samples (columns)
#' # And named feature vectors:
#' # feature_A <- c("A1", "A2", "A1", "A3", ...)
#' # names(feature_A) <- c("sample1", "sample2", "sample3", "sample4", ...)
#' # results <- compare_feature_clustering(rna_seq_data, feature_A, feature_B)
#'
compare_feature_clustering <- function(expr_data, feature1, feature2, 
                                      feature1_name = "Feature A",
                                      feature2_name = "Feature B",
                                      dist_method = "euclidean") {
  
  # Input validation
  if (!is.matrix(expr_data) && !is.data.frame(expr_data)) {
    stop("expr_data must be a matrix or data frame with genes as rows and samples as columns")
  }
  
  # Ensure proper orientation (genes as rows, samples as columns)
  if (ncol(expr_data) > nrow(expr_data)) {
    warning("Your expr_data appears to have more columns than rows. ",
            "Please ensure genes are rows and samples are columns.")
  }
  
  # Feature validation
  if (length(feature1) != ncol(expr_data) && 
      !all(names(feature1) %in% colnames(expr_data))) {
    stop("feature1 must have same length as ncol(expr_data) or have names matching column names")
  }
  
  if (length(feature2) != ncol(expr_data) && 
      !all(names(feature2) %in% colnames(expr_data))) {
    stop("feature2 must have same length as ncol(expr_data) or have names matching column names")
  }
  
  # If features are named, ensure they're in the right order
  if (!is.null(names(feature1))) {
    feature1 <- feature1[colnames(expr_data)]
  }
  
  if (!is.null(names(feature2))) {
    feature2 <- feature2[colnames(expr_data)]
  }
  
  # Calculate distance matrix - note we transpose to get samples as rows
  dist_mat <- dist(t(expr_data), method = dist_method)
  
  # Convert features to numeric cluster assignments
  cluster1 <- as.numeric(as.factor(feature1))
  cluster2 <- as.numeric(as.factor(feature2))
  
  # Calculate silhouette coefficients
  sil1 <- silhouette(cluster1, dist_mat)
  sil2 <- silhouette(cluster2, dist_mat)
  
  # Calculate average silhouette widths
  avg_sil1 <- mean(sil1[, "sil_width"])
  avg_sil2 <- mean(sil2[, "sil_width"])
  
  # Statistical test
  wilcox_result <- wilcox.test(sil1[, "sil_width"], sil2[, "sil_width"], paired = TRUE)
  
  # Results
  results <- list(
    silhouette1 = sil1,
    silhouette2 = sil2,
    avg_silhouette1 = avg_sil1,
    avg_silhouette2 = avg_sil2,
    p_value = wilcox_result$p.value,
    better_feature = ifelse(avg_sil1 > avg_sil2, feature1_name, feature2_name)
  )
  
  # Print summary
  cat("Average silhouette width for", feature1_name, ":", round(avg_sil1, 3), "\n")
  cat("Average silhouette width for", feature2_name, ":", round(avg_sil2, 3), "\n")
  cat("Wilcoxon test p-value:", round(wilcox_result$p.value, 4), "\n")
  
  if (wilcox_result$p.value < 0.05) {
    if (avg_sil1 > avg_sil2) {
      cat(feature1_name, "creates significantly better separated clusters\n")
    } else {
      cat(feature2_name, "creates significantly better separated clusters\n")
    }
  } else {
    cat("No significant difference in cluster separation\n")
  }
  
  # Return results object
  return(results)
}
