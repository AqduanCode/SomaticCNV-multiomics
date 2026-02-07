# =============================================================================
# Description: 
#   Creates per-chromosome read depth line plots for two samples, with sliding 
#   window smoothing and highlighted CNV regions called between the pair.
#   Designed for 293T lineage subclone comparison, but made generalizable.
# =============================================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(showtext)

# Enable better font rendering (especially useful on Windows)
showtext_auto()

# Add Arial font family (adjust paths if fonts are installed elsewhere)
font_add("Arial",
         regular    = "arial.ttf",
         bold       = "arialbd.ttf",
         italic     = "ariali.ttf",
         bolditalic = "arialbi.ttf")


#' Sample label → column name mapping table
#' 
#' This is project-specific. Users should provide their own mapping when 
#' calling main functions.
SAMPLE_MAPPING <- tibble(
  label    = c("H.1", "H.1.1", "H.1.2", "H.2", "H.2.1", "H.2.2", "H.3", "H.3.1", "H.3.2"),
  col_name = c("dfE9E10", "dfE9E10E8", "dfE9E10F5",
               "dfG6E3", "dfG6E3C3", "dfG6E3E7",
               "dfG6G7", "dfG6G7B1", "dfG6G7F5"),
  number   = 1:9
)


#' Custom sliding window mean (non-overlapping or stepped)
#'
#' @param x Numeric vector
#' @param window_size Size of the sliding window
#' @param step_size Step between window centers
#'
#' @return Numeric vector of same length as x, with NA outside valid range
sliding_mean_step <- function(x, window_size = 100, step_size = 10) {
  n <- length(x)
  out <- rep(NA_real_, n)
  
  half <- floor(window_size / 2)
  centers <- seq(
    from = ceiling(window_size / 2),
    to   = n - floor(window_size / 2),
    by   = step_size
  )
  
  for (c in centers) {
    idx <- (c - half):(c + half - 1)
    out[c] <- mean(x[idx], na.rm = TRUE)
  }
  
  out
}


#' Load and smooth read depth data for selected samples
#'
#' @param rd_file Path to read-depth matrix (.csv)
#' @param depth_col_names Character vector of column names containing raw depths
#' @param window_size Sliding window size (bins)
#' @param step_size Step size between smoothed points
#' @param chr_levels Ordered chromosome names (default chr1–chr22)
#'
#' @return Tibble with columns: Chromosome, BinLowEdge, BinUpEdge, smoothed_*
prepare_smoothed_depth <- function(rd_file,
                                   depth_col_names,
                                   window_size = 100,
                                   step_size = 10,
                                   chr_levels = paste0("chr", 1:22)) {
  
  RD <- read.csv(rd_file, stringsAsFactors = FALSE) %>%
    arrange(Chromosome, BinLowEdge) %>%
    mutate(Chromosome = factor(Chromosome, levels = chr_levels))
  
  # Apply smoothing to each selected column
  RD_smoothed <- RD %>%
    group_by(Chromosome) %>%
    mutate(
      across(
        all_of(depth_col_names),
        ~ sliding_mean_step(.x, window_size, step_size),
        .names = "smoothed_{.col}"
      )
    ) %>%
    ungroup()
  
  # Remove rows where any smoothed value is NA
  na.omit(RD_smoothed)
}


#' Prepare long-format data ready for plotting (two samples)
#'
#' @param smoothed_data Output from prepare_smoothed_depth()
#' @param sampleA_label Label of first sample
#' @param sampleB_label Label of second sample
#' @param mapping Data frame with label ↔ col_name mapping
#'
#' @return Long-format tibble with columns: Chromosome, BinLowEdge, sample, mean_depth
prepare_plot_data <- function(smoothed_data,
                              sampleA_label,
                              sampleB_label,
                              mapping = SAMPLE_MAPPING) {
  
  infoA <- mapping %>% filter(label == sampleA_label)
  infoB <- mapping %>% filter(label == sampleB_label)
  
  if (nrow(infoA) == 0 || nrow(infoB) == 0) {
    stop("One or both sample labels not found in mapping.")
  }
  
  colA <- paste0("smoothed_", infoA$col_name)
  colB <- paste0("smoothed_", infoB$col_name)
  
  smoothed_data %>%
    select(Chromosome, BinLowEdge, BinUpEdge, all_of(c(colA, colB))) %>%
    pivot_longer(
      cols      = c(all_of(colA), all_of(colB)),
      names_to  = "sample_raw",
      values_to = "mean_depth"
    ) %>%
    mutate(
      sample = case_when(
        sample_raw == colA ~ sampleA_label,
        sample_raw == colB ~ sampleB_label,
        TRUE ~ NA_character_
      )
    ) %>%
    select(-sample_raw)
}


#' Filter CNV calls for a specific sample pair
#'
#' @param cnv_file Path to CNV result table
#' @param sampleA_number Numeric ID of sample A
#' @param sampleB_number Numeric ID of sample B
#' @param chr_levels Ordered chromosome factor levels
#'
#' @return Filtered CNV tibble with factorized Chromosome
get_cnv_for_pair <- function(cnv_file,
                             sampleA_number,
                             sampleB_number,
                             chr_levels = paste0("chr", 1:22)) {
  
  read.csv(cnv_file, stringsAsFactors = FALSE) %>%
    filter(pair_sampleA == sampleA_number, pair_sampleB == sampleB_number) %>%
    mutate(Chromosome = factor(chr, levels = chr_levels))
}


#' Create read depth comparison plot for two samples with CNV regions highlighted
#'
#' @param sampleA_label Label of first sample
#' @param sampleB_label Label of second sample
#' @param rd_file Path to read depth matrix
#' @param cnv_file Path to CNV calling result table
#' @param mapping Sample label ↔ column mapping (default SAMPLE_MAPPING)
#' @param window_size Smoothing window size
#' @param step_size Smoothing step size
#' @param y_limits Numeric vector of length 2 for ylim
#' @param cnv_fill_color Fill color for CNV rectangles
#'
#' @return ggplot object
plot_depth_comparison <- function(sampleA_label,
                                  sampleB_label,
                                  rd_file,
                                  cnv_file,
                                  mapping      = SAMPLE_MAPPING,
                                  window_size  = 100,
                                  step_size    = 10,
                                  y_limits     = c(30, 250),
                                  cnv_fill_color = "#f4a7a6") {
  
  # ── Data preparation ───────────────────────────────────────────────────────
  
  depth_cols <- mapping$col_name
  
  smoothed <- prepare_smoothed_depth(
    rd_file         = rd_file,
    depth_col_names = depth_cols,
    window_size     = window_size,
    step_size       = step_size
  )
  
  plot_df <- prepare_plot_data(
    smoothed_data = smoothed,
    sampleA_label = sampleA_label,
    sampleB_label = sampleB_label,
    mapping       = mapping
  )
  
  infoA <- mapping %>% filter(label == sampleA_label)
  infoB <- mapping %>% filter(label == sampleB_label)
  
  cnv_df <- get_cnv_for_pair(
    cnv_file         = cnv_file,
    sampleA_number   = infoA$number,
    sampleB_number   = infoB$number
  )
  
  # ── Plotting ────────────────────────────────────────────────────────────────
  
  p <- ggplot(plot_df, aes(x = BinLowEdge, y = mean_depth, color = sample)) +
    
    # CNV background rectangles
    geom_rect(
      data = cnv_df,
      inherit.aes = FALSE,
      aes(xmin = start_pos, xmax = end_pos, ymin = -Inf, ymax = Inf,
          fill = "CNV region"),
      alpha = 0.4
    ) +
    scale_fill_manual(
      name   = NULL,
      values = c("CNV region" = cnv_fill_color)
    ) +
    
    # Depth lines
    geom_line(linewidth = 0.25) +
    
    scale_color_manual(
      values = c(
        setNames("#155075", sampleA_label),
        setNames("#d49053", sampleB_label)
      )
    ) +
    
    labs(
      title = NULL,
      x     = NULL,
      y     = "Read Depth"
    ) +
    
    facet_wrap(~ Chromosome, scales = "free_x", ncol = 6) +
    
    coord_cartesian(ylim = y_limits) +
    
    theme_minimal(base_family = "Arial", base_size = 8) +
    theme(
      axis.line          = element_line(color = "black", linewidth = 0.3),
      axis.ticks.y       = element_line(color = "black"),
      axis.ticks.length  = unit(0.05, "cm"),
      axis.text.x        = element_blank(),
      axis.text.y        = element_text(size = 8),
      axis.title.y       = element_text(size = 9, margin = margin(r = 10)),
      panel.grid.minor   = element_blank(),
      panel.grid.major   = element_line(color = "grey82", linewidth = 0.25),
      strip.text         = element_text(size = 8, face = "plain"),
      legend.position    = "bottom",
      legend.text        = element_text(size = 8),
      legend.title       = element_blank()
    )
  
  p
}


#' Batch generate comparison plots for multiple sample pairs
#'
#' @param pairs List of character vectors, each containing two sample labels
#' @param rd_file Read depth matrix path
#' @param cnv_file CNV result path
#' @param output_dir Directory to save PDF files
#' @param width,height Plot dimensions in inches
#' @param ... Additional arguments passed to plot_depth_comparison()
batch_plot_pairs <- function(pairs,
                             rd_file,
                             cnv_file,
                             output_dir,
                             width = 7,
                             height = 4.5,
                             ...) {
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (comb in pairs) {
    sampleA <- comb[1]
    sampleB <- comb[2]
    
    message("Generating plot: ", sampleA, " vs ", sampleB)
    
    p <- plot_depth_comparison(
      sampleA_label = sampleA,
      sampleB_label = sampleB,
      rd_file       = rd_file,
      cnv_file      = cnv_file,
      ...
    )
    
    cleanA <- gsub("\\.", "", sampleA)
    cleanB <- gsub("\\.", "", sampleB)
    filename <- sprintf("CNV_%s_vs_%s.pdf", cleanA, cleanB)
    
    ggsave(
      filename = file.path(output_dir, filename),
      plot     = p,
      width    = width,
      height   = height,
      device   = cairo_pdf
    )
  }
  
  invisible(TRUE)
}


# ──────────────────────────────────────────────────────────────────────────────
#  Example usage – 293T lineage
# ──────────────────────────────────────────────────────────────────────────────

main_293T_example <- function() {
  
  pairs <- list(
    c("H.1",   "H.1.1")
    # c("H.1",   "H.1.2"),
    # c("H.2",   "H.2.1"),
    # c("H.2",   "H.2.2"),
    # c("H.3",   "H.3.1"),
    # c("H.3",   "H.3.2")
  )
  
  batch_plot_pairs(
    pairs     = pairs,
    rd_file   = "T293metaDataDF.csv",
    cnv_file  = "CNVresult_T293.csv",
    output_dir = "figures/293T_CNV_comparisons/",
    window_size = 100,
    step_size   = 10,
    y_limits    = c(30, 250)
  )
}


main_293T_example()
