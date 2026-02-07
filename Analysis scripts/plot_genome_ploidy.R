# =============================================================================
# Description: Creates a smoothed read depth plot across chromosomes, highlighting
#              regions close to integer ploidy levels, with alternating chromosome
#              background and secondary ploidy axis.
# =============================================================================

library(dplyr)
library(ggplot2)
library(showtext)

# Optional: better font rendering (especially for Windows + Chinese path issues)
showtext_auto()

font_add("Arial", regular = "arial.ttf", bold = "arialbd.ttf",
         italic = "ariali.ttf", bolditalic = "arialbi.ttf")


#' Prepare read depth data for plotting
#'
#' @param meta_file Path to the metadata CSV file containing bin-level read depth
#' @param depth_cols Indices of columns that contain per-sample read depth values
#' @param smooth_window Number of consecutive bins to average (default: 10)
#'
#' @return A tibble with columns: chr, start, end, depth (smoothed), global_pos
prepare_depth_data <- function(meta_file, depth_cols, smooth_window = 10) {
  
  raw <- read.csv(meta_file, stringsAsFactors = FALSE)
  
  # Calculate mean read depth across selected samples
  raw$meanRD <- rowMeans(raw[, depth_cols, drop = FALSE], na.rm = TRUE)
  
  df <- raw %>%
    select(chr = 1, start = 2, end = 3, depth = meanRD) %>%
    mutate(chr = factor(chr, levels = paste0("chr", 1:22))) %>%
    arrange(chr, start)
  
  # Smooth: average every N bins within each chromosome
  df_smooth <- df %>%
    group_by(chr) %>%
    mutate(
      bin_id   = row_number(),
      group_id = (bin_id - 1) %/% smooth_window
    ) %>%
    group_by(chr, group_id) %>%
    summarise(
      depth = mean(depth, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(chr) %>%
    mutate(pos = row_number()) %>%
    ungroup()
  
  # Chromosome position offsets for genome-wide plotting
  chr_pos <- df_smooth %>%
    group_by(chr) %>%
    summarise(
      start = min(pos),
      end   = max(pos),
      .groups = "drop"
    ) %>%
    mutate(
      offset   = c(0, cumsum(end[-n()] + 10)),
      start    = start + offset,
      end      = end   + offset,
      mid      = (start + end) / 2,
      fill_bg  = rep(c("white", "grey96"), length.out = n())
    )
  
  # Add global genomic position
  df_smooth <- df_smooth %>%
    left_join(select(chr_pos, chr, offset), by = "chr") %>%
    mutate(global_pos = pos + offset)
  
  list(
    data     = df_smooth,
    chr_info = chr_pos
  )
}


#' Load ploidy estimation results and extract reference read depth values
#'
#' @param ploidy_file Path to CSV file containing ploidy estimation results
#' @param ploidy_range Integer vector of expected ploidy levels (e.g., 1:4 or 1:6)
#'
#' @return Numeric vector of mean read depths corresponding to each ploidy level
get_ploidy_reference_depths <- function(ploidy_file, ploidy_range = 1:4) {
  
  est <- read.csv(ploidy_file, stringsAsFactors = FALSE)

  values <- lapply(est$means, function(x) {
    x <- gsub("\\[|\\]", "", x)
    as.numeric(unlist(strsplit(x, ",")))
  })
  
  mat <- matrix(unlist(values), nrow = nrow(est), byrow = TRUE)
  rd_per_ploidy <- colMeans(mat, na.rm = TRUE)
  
  # Return only values corresponding to requested ploidy levels
  setNames(rd_per_ploidy[1:length(ploidy_range)], ploidy_range)
}


#' Add highlight flag for points close to integer ploidy levels
#'
#' @param df_smoothed Smoothed data frame from prepare_depth_data()
#' @param rd_ref Named numeric vector of reference read depths per ploidy
#' @param tolerance Distance threshold to consider "close" to ploidy level
#'
#' @return Updated data frame with new column 'highlight' (logical)
add_ploidy_highlight <- function(df_smoothed, rd_ref, tolerance = 0.15) {
  
  df_smoothed %>%
    mutate(
      highlight = sapply(depth, function(d) {
        any(abs(d - rd_ref) <= tolerance)
      })
    )
}


#' Create genome-wide read depth plot with ploidy annotation
#'
#' @param df_smoothed Data frame with global_pos, depth, chr, highlight
#' @param chr_info Chromosome layout information (from prepare_depth_data)
#' @param ploidy_ref Named vector of ploidy → read depth
#' @param highlight_color Color for highlighted points
#' @param point_size_normal Size of non-highlighted points
#' @param point_size_highlight Size of highlighted points
#' @param output_file Path to save the plot (PDF recommended)
#' @param width,height Plot dimensions in inches
#'
#' @return ggplot object (invisibly)
plot_genome_read_depth <- function(df_smoothed, chr_info, ploidy_ref,
                                   highlight_color = "grey70",
                                   point_size_normal   = 0.1,
                                   point_size_highlight = 0.05,
                                   output_file = NULL,
                                   width = 6.8, height = 3.5) {
  
  p <- ggplot() +
    
    # Alternating chromosome background
    geom_rect(
      data = chr_info,
      aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf, fill = fill_bg),
      inherit.aes = FALSE, alpha = 1
    ) +
    scale_fill_identity() +
    
    # Normal points
    geom_point(
      data = filter(df_smoothed, !highlight),
      aes(x = global_pos, y = depth, color = chr),
      size = point_size_normal, stroke = 0
    ) +
    
    # Highlighted points (close to integer ploidy)
    geom_point(
      data = filter(df_smoothed, highlight),
      aes(x = global_pos, y = depth),
      color = highlight_color, size = point_size_highlight, alpha = 0.8
    ) +
    
    scale_color_manual(
      values = c(
        "#40628c","#5e6c85","#6c7697","#8ca4c1",
        "#b7cdd7","#c5d1d5","#c4dbe3","#dce6ec",
        "#f3dcd1","#f6d0bb","#f5c6a9","#e2b7a8",
        "#e7a79b","#d98f84","#c97a74","#f8e8ec","#efd4db",
        "#e0b7bf","#d29da8","#ba8590","#ab737c","#9f5d62"
      ),
      guide = "none"
    ) +
    
    scale_x_continuous(
      breaks = chr_info$mid,
      labels = levels(chr_info$chr),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    
    scale_y_continuous(
      name   = "Read Depth",
      sec.axis = sec_axis(
        ~ .,
        breaks = as.numeric(ploidy_ref),
        labels = names(ploidy_ref),
        name   = "Ploidy"
      )
    ) +
    
    labs(x = NULL) +
    theme_bw(base_size = 8) +
    theme(
      panel.grid       = element_blank(),
      axis.text.x      = element_text(size = 7, family = "Arial", angle = 90, hjust = 1),
      axis.ticks.x     = element_blank(),
      axis.title.y     = element_text(size = 8, family = "Arial", margin = margin(r = 5)),
      axis.title.y.right = element_text(size = 8, family = "Arial", margin = margin(l = 5)),
      axis.text.y      = element_text(size = 7, family = "Arial"),
      axis.text.y.right = element_text(size = 7, family = "Arial"),
      legend.position  = "none"
    )
  
  if (!is.null(output_file)) {
    ggsave(
      filename = output_file,
      plot     = p,
      width    = width,
      height   = height,
      device   = cairo_pdf
    )
  }
  
  invisible(p)
}


# ──────────────────────────────────────────────────────────────────────────────
#  Example usage – U251
# ──────────────────────────────────────────────────────────────────────────────

main_U251 <- function() {
  
  data_res <- prepare_depth_data(
    meta_file   = "metaDataDF.csv",
    depth_cols  = 4:13,           # adjust according to your file
    smooth_window = 10
  )
  
  ploidy_rd <- get_ploidy_reference_depths(
    ploidy_file  = "Estimate.csv",
    ploidy_range = 1:4
  )
  
  df_plot <- add_ploidy_highlight(
    df_smoothed = data_res$data,
    rd_ref      = ploidy_rd,
    tolerance   = 0.15
  )
  
  plot_genome_read_depth(
    df_smoothed     = df_plot,
    chr_info        = data_res$chr_info,
    ploidy_ref      = ploidy_rd,
    output_file     = "Fig1-ploidy_U251.pdf",
    width           = 6.8,
    height          = 3.5
  )
}


main_U251()


# You can easily create main_293T() in the same style by changing:
#   - meta file path
#   - depth_cols = 4:12
#   - ploidy_file path
#   - ploidy_range = 1:6