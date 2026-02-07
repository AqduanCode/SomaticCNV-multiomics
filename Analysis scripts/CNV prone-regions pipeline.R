# =============================================================================
# Description:
#   This script integrates scaffold-level copy number estimates into interval-based
#   CNV calls, determines a representative copy number difference (ΔCN) per region
#   (using weighted mode by scaffold length), and detects overlapping/redundant
#   CNV calls across different sample pairs using graph-based connected components.
#
#   Designed for multi-generation cell line CNV analysis (e.g. 293T lineage).
# =============================================================================

library(dplyr)
library(tidyr)
library(purrr)
library(igraph)


# -----------------------------------------------------------------------------
# 1. Configuration – Sample label to scaffold column mapping
# -----------------------------------------------------------------------------
#' Mapping from sample label (as used in pair_sampleA / pair_sampleB)
#' to the corresponding column name in the scaffoldCN data frame.
SAMPLE_COLS_293T <- c(
  "H.1"    = "E9E10",
  "H.1.1"  = "E9E10E8",
  "H.1.2"  = "E9E10F5",
  "H.2"    = "G6E3",
  "H.2.1"  = "G6E3C3",
  "H.2.2"  = "G6E3E7",
  "H.3"    = "G6G7",
  "H.3.1"  = "G6G7B1",
  "H.3.2"  = "G6G7F5"
)


# -----------------------------------------------------------------------------
# 2. Core Processing Functions
# -----------------------------------------------------------------------------

#' Attach scaffold-level copy numbers to each CNV interval
#'
#' @param cnv_df Original CNV calls data frame
#' @param scaffold_df Scaffold copy number data frame
#' @param sample_col_map Named vector: sample label → scaffold column name
#'
#' @return Long-format tibble (one row per scaffold inside each CNV region)
attach_scaffold_cn <- function(cnv_df, scaffold_df, sample_col_map) {
  
  cnv_df %>%
    rowwise() %>%
    mutate(
      scaffold_rows = list(
        scaffold_df %>%
          filter(
            scaffold_start >= start_pos_metaData,
            scaffold_end   <= end_pos_metaData
          )
      )
    ) %>%
    ungroup() %>%
    unnest(scaffold_rows) %>%
    rowwise() %>%
    mutate(
      sampleA_CN = {
        cur <- pick(everything())
        if (cur$pair_sampleA %in% names(sample_col_map)) {
          cur[[ sample_col_map[[cur$pair_sampleA]] ]]
        } else {
          NA_integer_
        }
      },
      sampleB_CN = {
        cur <- pick(everything())
        if (cur$pair_sampleB %in% names(sample_col_map)) {
          cur[[ sample_col_map[[cur$pair_sampleB]] ]]
        } else {
          NA_integer_
        }
      }
    ) %>%
    ungroup() %>%
    select(
      pair_sampleA, pair_sampleB, chr, start_pos, end_pos, length,
      start_pos_metaData, end_pos_metaData,
      mean_depth_A, mean_depth_B, depth_diff,
      scaffold_length, sampleA_CN, sampleB_CN
    )
}


#' Summarize representative CNV difference per called region
#' (prefers longest scaffold length when ambiguous)
#'
#' @param long_df Output from attach_scaffold_cn()
#'
#' @return Summarised tibble (one row per original CNV call) with CNV column
summarise_cnv_difference <- function(long_df) {
  
  long_df %>%
    group_by(pair_sampleA, pair_sampleB, chr, start_pos, end_pos) %>%
    summarise(
      CNV = {
        valid_rows <- sampleA_CN != -1 & sampleB_CN != -1
        
        if (!any(valid_rows)) {
          0L
        } else {
          diff_vals <- sampleB_CN[valid_rows] - sampleA_CN[valid_rows]
          lens      <- scaffold_length[valid_rows]
          
          nonzero         <- diff_vals != 0
          diff_nonzero    <- diff_vals[nonzero]
          len_nonzero     <- lens[nonzero]
          
          if (length(diff_nonzero) == 0) {
            0L
          } else {
            temp_df <- tibble(diff = diff_nonzero, len = len_nonzero)
            
            agg <- temp_df %>%
              group_by(diff) %>%
              summarise(
                total_len = sum(len),
                n         = n(),
                .groups   = "drop"
              )
            
            if (nrow(agg) == 2) {
              # When exactly two different diffs → prefer the one with largest total length
              agg %>%
                arrange(desc(total_len), desc(n)) %>%
                slice(1) %>%
                pull(diff) %>%
                as.integer()
            } else {
              # Otherwise → classic mode (most frequent diff)
              tab <- table(diff_nonzero)
              as.integer(names(tab)[which.max(tab)])
            }
          }
        }
      },
      .groups = "drop"
    )
}


#' Add pair identifier column (e.g. "H.1 - H.1.1")
#'
#' @param df Data frame containing pair_sampleA and pair_sampleB
#'
#' @return Same data frame with added pair_id column
add_pair_id <- function(df) {
  df %>%
    rowwise() %>%
    mutate(pair_id = paste(pair_sampleA, "-", pair_sampleB)) %>%
    ungroup()
}


#' Detect and cluster overlapping CNV calls using graph connected components
#'
#' @param df Data frame with at least chr, start_pos, end_pos, pair_id, CNV
#' @param overlap_threshold Minimum overlap ratio to consider two regions overlapping (default 0.9)
#'
#' @return List with edge_cnv, cluster_assignments, row_level, cluster_level
get_overlap_dup <- function(df, overlap_threshold = 0.9) {
  
  library(dplyr)
  library(igraph)
  library(tidyr)
  library(purrr)
  
  # 0. Prepare identifiers and lengths
  df <- df %>%
    mutate(
      length = end_pos - start_pos,
      cnv_id = paste(chr, start_pos, end_pos, pair_id, sep = "|")
    )
  
  # 1. Self-join on same chromosome
  df_joined <- df %>%
    inner_join(df, by = "chr", suffix = c("_A", "_B")) %>%
    filter(cnv_id_A != cnv_id_B)
  
  # 2. Calculate overlap
  overlap_df <- df_joined %>%
    mutate(
      overlap_len   = pmax(0, pmin(end_pos_A, end_pos_B) - pmax(start_pos_A, start_pos_B)),
      min_len       = pmin(length_A, length_B),
      overlap_ratio = overlap_len / min_len,
      overlap_start = pmax(start_pos_A, start_pos_B),
      overlap_end   = pmin(end_pos_A, end_pos_B)
    ) %>%
    filter(overlap_ratio >= overlap_threshold)
  
  # 3. Row-level summary
  row_base <- df %>%
    select(chr, start_pos, end_pos, pair_id) %>%
    mutate(
      self_region = pmap(
        list(chr, start_pos, end_pos, pair_id),
        ~ tibble(chr = ..1, start_pos = ..2, end_pos = ..3, pair_id = ..4)
      )
    )
  
  overlapping_summary <- overlap_df %>%
    group_by(chr, start_pos = start_pos_A, end_pos = end_pos_A, pair_id = pair_id_A) %>%
    summarise(
      high_overlap_pair_ids = list(unique(pair_id_B)),
      n_high_overlap_pairs  = n_distinct(pair_id_B),
      max_overlap_ratio     = max(overlap_ratio),
      overlapping_regions   = list(
        tibble(chr = chr, start_pos = start_pos_B, end_pos = end_pos_B, pair_id = pair_id_B) %>% distinct()
      ),
      .groups = "drop"
    )
  
  row_level <- row_base %>%
    left_join(overlapping_summary, by = c("chr", "start_pos", "end_pos", "pair_id")) %>%
    replace_na(list(
      high_overlap_pair_ids = list(character(0)),
      n_high_overlap_pairs  = 0L,
      max_overlap_ratio     = NA_real_,
      overlapping_regions   = list(tibble(chr=character(0), start_pos=integer(0), end_pos=integer(0), pair_id=character(0)))
    )) %>%
    mutate(
      source_CNV_regions = map2(self_region, overlapping_regions, bind_rows)
    ) %>%
    select(-self_region)
  
  # 4. Edge list for graph
  edge_cnv <- overlap_df %>%
    transmute(chr = chr, from = cnv_id_A, to = cnv_id_B) %>%
    distinct()
  
  # 5. Connected components per chromosome
  cluster_assignments <- edge_cnv %>%
    group_by(chr) %>%
    group_modify(function(.x, .y) {
      if (nrow(.x) == 0) return(tibble(cnv_id = character(), cluster_id = integer()))
      g <- graph_from_data_frame(.x %>% select(from, to), directed = FALSE)
      comp <- components(g)
      tibble(cnv_id = names(comp$membership), cluster_id = comp$membership)
    }) %>%
    ungroup()
  
  # 6. Cluster-level summaries
  cluster_source_regions <- cluster_assignments %>%
    left_join(df %>% select(chr, cnv_id, start_pos, end_pos, pair_id), by = c("chr", "cnv_id")) %>%
    group_by(chr, cluster_id) %>%
    summarise(
      source_CNV_regions = list(tibble(chr = chr, start_pos = start_pos, end_pos = end_pos, pair_id = pair_id)),
      .groups = "drop"
    )
  
  cluster_position <- cluster_assignments %>%
    left_join(df %>% select(chr, cnv_id, start_pos, end_pos, pair_id), by = c("chr", "cnv_id")) %>%
    group_by(chr, cluster_id) %>%
    summarise(
      start_pos        = min(start_pos),
      end_pos          = max(end_pos),
      consensus_length = end_pos - start_pos,
      n_pairs          = n_distinct(pair_id),
      .groups          = "drop"
    )
  
  cluster_pair_cnv <- cluster_assignments %>%
    left_join(df %>% select(chr, cnv_id, pair_id, CNV), by = c("chr", "cnv_id")) %>%
    group_by(chr, cluster_id, pair_id) %>%
    summarise(CNV = first(CNV), .groups = "drop") %>%
    group_by(chr, cluster_id) %>%
    summarise(
      pair_list = list(pair_id),
      CNV_list  = list(CNV),
      .groups   = "drop"
    )
  
  cluster_level <- cluster_position %>%
    left_join(cluster_pair_cnv, by = c("chr", "cluster_id")) %>%
    left_join(cluster_source_regions, by = c("chr", "cluster_id")) %>%
    mutate(
      n_CN_types = map_int(CNV_list, ~ {
        cn_vals <- unlist(lapply(.x, function(one) {
          unlist(strsplit(gsub(" ", "", one), ">|<|="))
        }))
        length(unique(as.integer(cn_vals)))
      })
    ) %>%
    select(
      chr, cluster_id, start_pos, end_pos, consensus_length,
      pair_list, n_pairs, CNV_list, n_CN_types, source_CNV_regions
    )
  
  list(
    edge_cnv            = edge_cnv,
    cluster_assignments = cluster_assignments,
    row_level           = row_level,
    cluster_level       = cluster_level
  )
}


# -----------------------------------------------------------------------------
# Main execution function – 293T example
# -----------------------------------------------------------------------------
run_293t_analysis <- function(
    cnv_path = "CNVresult_293T.csv",
    scaffold_path = "scaffoldCN.csv",
    sample_map = SAMPLE_COLS_293T
) {
  
  # Load data
  CNVresult  <- read.csv(cnv_path)
  scaffoldCN <- read.csv(scaffold_path)
  
  # Clean columns
  CNVresult <- CNVresult %>% select(-any_of(c("sampleA_CN", "sampleB_CN")))
  
  # 1. Attach scaffold copy numbers
  long_df <- attach_scaffold_cn(CNVresult, scaffoldCN, sample_map)
  
  # 2. Summarize representative CN difference
  CNV_with_diff <- summarise_cnv_difference(long_df)
  
  # 3. Add pair identifier
  CNVresult <- add_pair_id(CNV_with_diff)
  
  # 4. Detect overlapping CNVs
  res_overlap <- get_overlap_dup(CNVresult, overlap_threshold = 0.9)
  
  # 5. Extract results
  overlap_row    <- res_overlap$row_level
  overlap_cluster <- res_overlap$cluster_level
  
  # 6. Example filtering (customize as needed)
  overlap_cluster_filtered <- overlap_cluster %>%
    filter(
      n_pairs != 2 |
        n_CN_types != 1 |
        map_lgl(pair_list, ~ {
          if (length(.x) < 2) return(TRUE)
          prefixes <- strsplit(.x, " - ")
          prefixes[[1]][1] != prefixes[[2]][1]
        })
    )
  
  # Return all main objects
  list(
    long_df               = long_df,
    CNV_with_diff         = CNV_with_diff,
    overlap_row           = overlap_row,
    overlap_cluster       = overlap_cluster,
    overlap_cluster_filtered = overlap_cluster_filtered
  )
}


# -----------------------------------------------------------------------------
# Run the analysis
# -----------------------------------------------------------------------------
results <- run_293t_analysis()
