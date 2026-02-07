# =============================================================================
# Description:   Merge adjacent or nearby CNV candidate regions based on read-depth
#            difference. Bridge small-to-medium genomic gaps only if the gap
#            region itself shows consistent depth anomaly (abnormal diffRD).
#            Used as a post-filtering step after initial CNV calling.
# =============================================================================

library(dplyr)
library(tibble)


# -----------------------------------------------------------------------------
# 1. Merging Function
# -----------------------------------------------------------------------------

#' Merge adjacent/nearby CNV regions with conditional gap bridging
#'
#' Processes CNV candidates on one chromosome for one sample pair.
#' Merging rules:
#'   - Small gaps (< min_gap): merge if regions are significant
#'   - Medium gaps (min_gap ≤ gap < max_gap): merge only if gap shows strong depth difference
#'   - Large gaps (≥ max_gap): never merge
#' Regions with small depth difference are kept as singletons.
#'
#' @param group_df tibble. CNV candidates for one (pair_sampleA, pair_sampleB, chr) group.
#'                 Required columns: start_pos, end_pos, length, mean_depth_A, mean_depth_B,
#'                 depth_diff, start_pos_metaData, end_pos_metaData
#' @param group_keys tibble (1 row). Grouping keys: pair_sampleA, pair_sampleB, chr
#' @param meta_data_df tibble. Bin-level depth metadata with columns:
#'                     Chromosome, BinLowEdge, BinUpEdge, and sample depth columns
#'                     
#' @param min_gap numeric. Minimum gap size requiring gap check (bp, default 1e6)
#' @param max_gap numeric. Maximum allowable gap for bridging (bp, default 3e6)
#' @param diff_threshold numeric. |depth_diff| below which region is insignificant
#' @param gap_diff_threshold numeric. |gap depth difference| required to bridge gap
#'
#' @return tibble. Merged regions with columns:
#'         start_pos, end_pos, length, start_pos_metaData, end_pos_metaData,
#'         mean_depth_A, mean_depth_B, depth_diff
#'         (Grouping columns removed)

         
merge_group <- function(group_df, group_keys, 
                        meta_data_df,
                        min_gap = 1e6, 
                        max_gap = 3e6,
                        diff_threshold = 15, 
                        gap_diff_threshold = 15) {
  
  # Return empty tibble if input group has no rows
  if (nrow(group_df) == 0) {
    return(tibble())
  }
  
  # Ensure regions are sorted by start position
  group_df <- group_df %>%
    arrange(start_pos)
  
  # Fixed column name mapping from meta_data_df (sample depth columns start from column 4)
  sample_cols <- colnames(meta_data_df[, 4:ncol(meta_data_df)])
  
  # Safely extract grouping keys from the single-row tibble
  pairA <- group_keys[["pair_sampleA"]][1]
  pairB <- group_keys[["pair_sampleB"]][1]
  chr   <- as.character(group_keys[["chr"]][1])
  
  # Map pair identifiers to corresponding column names in meta_data_df
  colA <- sample_cols[pairA]
  colB <- sample_cols[pairB]
  
  # Pre-filter meta data for the current chromosome (executed only once)
  meta_chr <- meta_data_df %>%
    filter(Chromosome == chr)
  
  # Columns to exclude from final output (grouping keys)
  group_cols <- c("pair_sampleA", "pair_sampleB", "chr")
  
  merged_rows <- list()
  i <- 1
  
  while (i <= nrow(group_df)) {
    curr <- group_df[i, ]
    
    # If the absolute depth difference is ≤ threshold, keep as a single row
    if (abs(curr$depth_diff) <= diff_threshold) {
      single_row <- curr[, !names(curr) %in% group_cols, drop = FALSE]
      merged_rows[[length(merged_rows) + 1]] <- single_row
      i <- i + 1
      next
    } else {
      # Start a new merge cluster (current region must be significant)
      cluster <- list()
      cluster[[1]] <- curr
      current_end <- curr$end_pos
      total_length <- curr$length
      sum_A <- curr$mean_depth_A * curr$length
      sum_B <- curr$mean_depth_B * curr$length
      start_meta <- curr$start_pos_metaData
      end_meta <- curr$end_pos_metaData
      
      j <- i + 1
      
      while (j <= nrow(group_df)) {
        next_r <- group_df[j, ]
        gap_genomic <- next_r$start_pos - current_end
        
        # Gap ≥ max_gap: stop merging
        if (gap_genomic >= max_gap) {
          break
        }
        
        # Calculate Gap_diffRD (only when there is a positive gap)
        can_merge_gap <- TRUE
        if (gap_genomic > 0) {
          gap_start <- current_end + 1
          gap_end <- next_r$start_pos - 1
          
          gap_bins <- meta_chr %>%
            filter(BinUpEdge > gap_start & BinLowEdge < gap_end) %>%
            mutate(
              overlap_start  = pmax(BinLowEdge, gap_start),
              overlap_end    = pmin(BinUpEdge,   gap_end),
              overlap_length = overlap_end - overlap_start
            ) %>%
            filter(overlap_length > 0)
          
          if (nrow(gap_bins) == 0 || sum(gap_bins$overlap_length) == 0) {
            can_merge_gap <- FALSE
          } else {
            sum_A_weighted <- sum(gap_bins[[colA]] * gap_bins$overlap_length)
            sum_B_weighted <- sum(gap_bins[[colB]] * gap_bins$overlap_length)
            total_gap_len  <- sum(gap_bins$overlap_length)
            mean_A_gap     <- sum_A_weighted / total_gap_len
            mean_B_gap     <- sum_B_weighted / total_gap_len
            gap_diffRD     <- mean_A_gap - mean_B_gap
            
            can_merge_gap <- abs(gap_diffRD) > gap_diff_threshold
          }
        }
        
        # For gaps >= min_gap, require the next region to be significant
        if (gap_genomic >= min_gap && abs(next_r$depth_diff) <= diff_threshold) {
          break
        }
        
        # For gaps >= min_gap, require the gap region itself to be abnormal
        if (gap_genomic >= min_gap && !can_merge_gap) {
          break
        }
        
        # All conditions satisfied: merge the next region
        cluster[[length(cluster) + 1]] <- next_r
        current_end <- next_r$end_pos
        total_length <- total_length + next_r$length
        sum_A <- sum_A + next_r$mean_depth_A * next_r$length
        sum_B <- sum_B + next_r$mean_depth_B * next_r$length
        end_meta <- next_r$end_pos_metaData
        
        # Note: Removed forced stop after medium gap merge
        # This allows continued merging of subsequent small gaps after bridging
        
        j <- j + 1
      }
      
      # Create merged row (weighted average)
      new_start  <- cluster[[1]]$start_pos
      new_end    <- current_end
      new_length <- new_end - new_start
      new_mean_A <- sum_A / total_length
      new_mean_B <- sum_B / total_length
      new_diff   <- new_mean_A - new_mean_B
      
      new_row <- tibble(
        start_pos          = new_start,
        end_pos            = new_end,
        length             = new_length,
        start_pos_metaData = start_meta,
        end_pos_metaData   = end_meta,
        mean_depth_A       = new_mean_A,
        mean_depth_B       = new_mean_B,
        depth_diff         = new_diff
      )
      
      merged_rows[[length(merged_rows) + 1]] <- new_row
      i <- j
    }
  }
  
  if (length(merged_rows) == 0) {
    return(tibble())
  }
  
  bind_rows(merged_rows)
}




# -----------------------------------------------------------------------------
# 2. Example Usage (293T)
# -----------------------------------------------------------------------------

# Load data

plotCNVregion_T293 <- read.csv("path/to/plotCNVregion_T293.csv")
T293metaDataDF     <- read.csv("path/to/T293metaDataDF.csv")


merged_T293 <- plotCNVregion_T293 %>%
  group_by(pair_sampleA, pair_sampleB, chr) %>%
  group_modify(~ merge_group(.x, .y,
                             meta_data_df = T293metaDataDF,
                             min_gap = 1e6,
                             max_gap = 3e6,
                             diff_threshold = 15,
                             gap_diff_threshold = 10)) %>%
  ungroup()



# Optional: apply final significance filter
CNVresult <- merged_T293 %>%
  filter(abs(depth_diff) >= 19)
