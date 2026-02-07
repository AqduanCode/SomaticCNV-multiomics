# =============================================================================
# Description:
#   This script performs empirical Kolmogorov-Smirnov (KS) testing to assess whether
#   β-values (methylation levels) in CNV regions differ significantly from non-CNV
#   background regions. It uses Monte Carlo simulation for observed statistics and
#   empirical null distribution generation.
#
#   Designed for multi-generation cell line methylation analysis (293T, G6, C10).
#   Supports different input RDS files and customizable simulation parameters.
# =============================================================================

library(dplyr)


# -----------------------------------------------------------------------------
# Core Functions
# -----------------------------------------------------------------------------

#' Load methylation data for multiple samples
#'
#' @param data_paths Named list: sample name → path to RDS file
#'
#' @return Named list: sample name → methylation data frame (with pair_sampleA_beta list-column)
load_methylation_data <- function(data_paths) {
  lapply(data_paths, readRDS)
}


#' Simulate observed KS statistics for a given CNV state
#'
#' @param selectFeature CNV state to test (e.g., "CNVregion")
#' @param num_simu Number of Monte Carlo simulations
#' @param df Methylation data frame with columns: state, pair_sampleA_beta (list of betas)
#' @param rate Fraction of selected fragments to use in each simulation (default 1)
#'
#' @return Matrix with num_simu rows and 2 columns: KS statistic, p-value
simu_KS <- function(selectFeature, num_simu, df, rate = 1) {
  
  # Remove rows with no beta values
  temp <- sapply(df$pair_sampleA_beta, length)
  valid_idx <- which(temp > 0)
  subdf <- df[valid_idx, ]
  
  # Select fragments of interest
  idx_selected <- which(subdf$state == selectFeature)
  num_selected <- length(idx_selected)
  
  # Number of fragments after rate adjustment
  n_select_rate <- max(1, round(num_selected * rate))
  
  simu_results <- matrix(-1, nrow = num_simu, ncol = 2)
  
  for (i in 1:num_simu) {
    # Sample rate proportion of target fragments
    if (n_select_rate < num_selected) {
      temp_idx_select <- sample(idx_selected, n_select_rate, replace = FALSE)
    } else {
      temp_idx_select <- idx_selected
    }
    
    selected_beta <- unlist(subdf$pair_sampleA_beta[temp_idx_select])
    
    # Sample equal number of background fragments
    temp_idx_bg <- sample(seq_len(nrow(subdf)), length(temp_idx_select), replace = FALSE)
    bg_beta <- unlist(subdf$pair_sampleA_beta[temp_idx_bg])
    
    # Perform KS test
    ks_res <- ks.test(selected_beta, bg_beta)
    simu_results[i, ] <- c(ks_res$statistic, ks_res$p.value)
  }
  
  simu_results
}


#' Simulate empirical null distribution (background vs background)
#'
#' @param selectFeature CNV state to test (e.g., "CNVregion")
#' @param num_simu Number of Monte Carlo simulations
#' @param df Methylation data frame
#' @param rate Fraction of fragments to use (same as in simu_KS)
#'
#' @return Matrix with num_simu rows and 2 columns: KS statistic, p-value
simu_KSsimu <- function(selectFeature, num_simu, df, rate = 1) {
  
  # Same preprocessing as simu_KS
  temp <- sapply(df$pair_sampleA_beta, length)
  valid_idx <- which(temp > 0)
  subdf <- df[valid_idx, ]
  
  idx_selected <- which(subdf$state == selectFeature)
  num_selected <- length(idx_selected)
  
  n_select_rate <- max(1, round(num_selected * rate))
  
  # Initial background sample
  temp_idx_init <- sample(seq_len(nrow(subdf)), n_select_rate, replace = FALSE)
  selected_beta <- unlist(subdf$pair_sampleA_beta[temp_idx_init])
  
  simu_results <- matrix(-1, nrow = num_simu, ncol = 2)
  
  for (i in 1:num_simu) {
    # New background sample of same size
    temp_idx_bg <- sample(seq_len(nrow(subdf)), n_select_rate, replace = FALSE)
    bg_beta <- unlist(subdf$pair_sampleA_beta[temp_idx_bg])
    
    ks_res <- ks.test(selected_beta, bg_beta)
    simu_results[i, ] <- c(ks_res$statistic, ks_res$p.value)
    
    # Chain update (use previous background as new "observed" for null)
    selected_beta <- bg_beta
  }
  
  simu_results
}


#' Run empirical KS test (observed vs empirical null)
#'
#' @param sample_name Name of the sample/cell line
#' @param df Methylation data frame for this sample
#' @param selectFeature CNV state to test (default: "CNVregion")
#' @param rate Fraction of fragments to use (default: 1)
#' @param n_obs Number of simulations for observed statistic
#' @param n_bg Number of simulations for empirical null
#'
#' @return One-row data.frame with sample, state, rate, mean_ks_obs, n_obs, n_bg, p_empirical
run_KS_empirical <- function(sample_name,
                             df,
                             selectFeature = "CNVregion",
                             rate = 1,
                             n_obs = 1000,
                             n_bg = 1000) {
  
  # Observed KS statistics
  test_obs <- simu_KS(
    selectFeature = selectFeature,
    num_simu = n_obs,
    df = df,
    rate = rate
  )
  
  # Empirical null (background vs background)
  test_bg <- simu_KSsimu(
    selectFeature = selectFeature,
    num_simu = n_bg,
    df = df,
    rate = rate
  )
  
  ks_obs <- test_obs[, 1]
  ks_bg  <- test_bg[, 1]
  
  T_obs <- mean(ks_obs)
  
  # Empirical p-value (one-sided, greater)
  p_empirical <- (sum(ks_bg >= T_obs) + 1) / (length(ks_bg) + 1)
  
  data.frame(
    sample        = sample_name,
    state         = selectFeature,
    rate          = rate,
    mean_ks_obs   = T_obs,
    n_obs         = n_obs,
    n_bg          = n_bg,
    p_empirical   = p_empirical,
    stringsAsFactors = FALSE
  )
}


# -----------------------------------------------------------------------------
# Main workflow – Run analysis across multiple samples and rates
# -----------------------------------------------------------------------------
run_ks_analysis <- function(
    data_list,
    rate_vec = c(0.1, 0.5, 0.8, 1),
    selectFeature = "CNVregion",
    n_obs = 1000,
    n_bg = 1000,
    output_path = "CNV_region_empirical_p_values.csv"
) {
  
  result_list <- lapply(names(data_list), function(samp) {
    lapply(rate_vec, function(r) {
      run_KS_empirical(
        sample_name   = samp,
        df            = data_list[[samp]],
        selectFeature = selectFeature,
        rate          = r,
        n_obs         = n_obs,
        n_bg          = n_bg
      )
    }) |> bind_rows()
  }) |> bind_rows()
  
  write.csv(
    result_list,
    file = output_path,
    row.names = FALSE
  )
  
  message("Analysis completed. Results saved to: ", output_path)
  message("Number of rows in result: ", nrow(result_list))
  
  invisible(result_list)
}


# -----------------------------------------------------------------------------
# Example execution (293T, G6, C10)
# -----------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# Configuration – Input data paths and sample names
# -----------------------------------------------------------------------------

SAMPLE_NAMES <- c("HEK293T", "G6", "C10")

DATA_PATHS <- list(
  HEK293T = "HEK293T-CNV_NonCNVregion.rds",
  G6      = "U251serialA-CNV_NonCNVregion.rds",
  C10     = "U251serialB-CNV_NonCNVregion.rds"
)

# Load data
CNVmeth_list <- lapply(DATA_PATHS, readRDS)

# Run full analysis
results <- run_ks_analysis(
  data_list     = CNVmeth_list,
  rate_vec      = c(0.1, 0.5, 0.8, 1),
  selectFeature = "CNVregion",
  n_obs         = 1000,
  n_bg          = 1000,
  output_path   = "CNV_region_empirical_p_values.csv"
)
