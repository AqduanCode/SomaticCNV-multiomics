# ============================================================
# CNV–RNAseq Dosage Analysis (293T)
# Description:
#   Integrate CNV regions with RNA-seq expression data
#   and quantify expression dosage effects.
# ============================================================

library(dplyr)
library(GenomicRanges)
library(rtracklayer)


# ------------------------------------------------------------
# 1. Load input data
# ------------------------------------------------------------
#' Load CNV results, RNA-seq logTPM matrix, and filtered gene list.
#' Rename columns to sample names and compute theoretical expression difference.
#'
#' @param cnv_path Path to CNV result CSV file
#' @param rnaseq_path Path to logTPM RNA-seq matrix CSV file
#' @param gene_filtered_path Path to filtered gene list CSV file
#' @param sample_names Character vector of sample names for column renaming
#'
#' @return A named list containing:
#'   - CNVresult: data.frame with added theoDiffexp column
#'   - rnaseq: logTPM expression matrix with renamed columns
#'   - gene_filtered: Character vector of retained gene IDs
load_input_data <- function(
    cnv_path,
    rnaseq_path,
    gene_filtered_path,
    sample_names
) {
  CNVresult <- read.csv(cnv_path)
  CNVresult <- as.data.frame(CNVresult)
  
  # Theoretical expression change from CNV
  CNVresult$theoDiffexp <- log10(
    CNVresult$sampleB_CN / CNVresult$sampleA_CN
  )
  
  rnaseq <- read.csv(rnaseq_path, row.names = 1)
  colnames(rnaseq) <- sample_names
  
  gene_filtered <- read.csv(gene_filtered_path, row.names = 1)
  colnames(gene_filtered) <- sample_names
  
  list(
    CNVresult = CNVresult,
    rnaseq = rnaseq,
    gene_filtered = gene_filtered
  )
}


# ------------------------------------------------------------
# 2. Build GRanges for CNV regions
# ------------------------------------------------------------
#' Convert CNV intervals into a single GRanges object.
#' Set sequence levels to all unique chromosomes in the data.
#'
#' @param CNVresult Data frame containing CNV calls with chr, start_pos, end_pos
#'
#' @return A GRanges object representing all CNV regions
build_cnv_granges <- function(CNVresult) {
  all_seqlevels <- unique(CNVresult$chr)
  
  gr_list <- lapply(seq_len(nrow(CNVresult)), function(i) {
    row <- CNVresult[i, ]
    gr <- GRanges(
      seqnames = row$chr,
      ranges = IRanges(
        start = as.numeric(row$start_pos),
        end = as.numeric(row$end_pos)
      )
    )
    seqlevels(gr) <- all_seqlevels
    gr
  })
  
  Reduce(c, gr_list)
}


# ------------------------------------------------------------
# 3. Extract genes fully within CNV regions
# ------------------------------------------------------------
#' Load gene annotations from GTF, find genes fully contained in CNV regions,
#' clean Ensembl IDs, and filter to genes present in the filtered list.
#'
#' @param cnv_gr GRanges object of CNV regions
#' @param CNVresult Original CNV data frame (for metadata)
#' @param gtf_path Path to GTF annotation file
#' @param gene_filtered Character vector of retained gene IDs
#'
#' @return Data frame of overlapping gene-CNV pairs (filtered)
extract_cnv_genes <- function(
    cnv_gr,
    CNVresult,
    gtf_path,
    gene_filtered
) {
  genes <- import(gtf_path, format = "gtf")
  genes <- genes[genes$type == "gene", ]
  
  hits <- findOverlaps(genes, cnv_gr, type = "within")
  
  overlap_df <- data.frame(
    GeneID = genes$gene_id[queryHits(hits)],
    CNVresult[subjectHits(hits), ]
  )
  
  overlap_df$GeneID <- gsub("\\.\\d*", "", overlap_df$GeneID)
  
  overlap_df[
    overlap_df$GeneID %in% rownames(gene_filtered),
  ]
}


# ------------------------------------------------------------
# 4. Integrate RNA-seq expression
# ------------------------------------------------------------
#' Merge CNV-overlapping genes with RNA-seq expression values.
#' Compute mean expression, remove specified samples,
#' extract sample-specific expression, and calculate differences and labels.
#'
#' @param overlap_df Data frame from extract_cnv_genes()
#' @param rnaseq RNA-seq logTPM matrix
#' @param remove_samples Character vector of sample labels to exclude
#'
#' @return Data frame with expression values, meanrnaseq, actExpdiff, group, CNV
integrate_expression <- function(
    overlap_df,
    rnaseq,
    remove_samples = c("H.3.1")
) {
  stopifnot(
    "pair_sampleA" %in% colnames(overlap_df),
    "pair_sampleB" %in% colnames(overlap_df)
  )
  
  rnaseq$GeneID <- rownames(rnaseq)
  
  CNVgene_rnaseq <- merge(
    overlap_df,
    rnaseq[, 1:ncol(rnaseq), drop = FALSE],
    by = "GeneID",
    all.x = TRUE
  )
  
  CNVgene_rnaseq <- CNVgene_rnaseq %>%
    mutate(
      meanrnaseq = rowMeans(
        across((ncol(.) - ncol(rnaseq) + 2):ncol(.)),
        na.rm = TRUE
      )
    ) %>%
    filter(
      !(pair_sampleA %in% remove_samples |
          pair_sampleB %in% remove_samples)
    ) %>%
    rowwise() %>%
    mutate(
      sampleA_exp = cur_data()[[ as.character(pair_sampleA) ]],
      sampleB_exp = cur_data()[[ as.character(pair_sampleB) ]]
    ) %>%
    ungroup()
  
  CNVgene_rnaseq$actExpdiff <-
    CNVgene_rnaseq$sampleB_exp -
    CNVgene_rnaseq$sampleA_exp
  
  CNVgene_rnaseq$group <- paste0(
    CNVgene_rnaseq$sampleA_CN, "C > ",
    CNVgene_rnaseq$sampleB_CN, "C"
  )
  
  CNVgene_rnaseq$CNV <-
    CNVgene_rnaseq$sampleB_CN -
    CNVgene_rnaseq$sampleA_CN
  
  CNVgene_rnaseq
}


# ------------------------------------------------------------
# 5. Calculate dosage effect
# ------------------------------------------------------------
#' Filter genes based on mean expression cutoff and calculate dosage effect.
#'
#' @param CNVgene_rnaseq Data frame from integrate_expression()
#' @param mean_expr_cutoff Minimum mean expression threshold (default: 1)
#'
#' @return Filtered data frame with added dosage column (10^actExpdiff)
calculate_dosage <- function(
    CNVgene_rnaseq,
    mean_expr_cutoff = 1
) {
  res <- CNVgene_rnaseq[
    CNVgene_rnaseq$meanrnaseq > mean_expr_cutoff,
  ]
  
  res$dosage <- 10 ^ res$actExpdiff
  res
}


# ============================================================
# Example: 293T analysis
# ============================================================
sample_names <- c(
  "H.1.1","H.1.2","H.1",
  "H.2.1","H.2.2","H.2",
  "H.3.2","H.3"
)

data_list <- load_input_data(
  cnv_path = "CNVresult_293T.csv",
  rnaseq_path = "T293_gene.csv",
  gene_filtered_path = "gene_filtered.csv",
  sample_names = sample_names
)

cnv_gr <- build_cnv_granges(data_list$CNVresult)

overlap_df_diff <- extract_cnv_genes(
  cnv_gr = cnv_gr,
  CNVresult = data_list$CNVresult,
  gtf_path = "gencode.v47.chr_patch_hapl_scaff.annotation.gtf",
  gene_filtered = data_list$gene_filtered
)

CNVgene_rnaseq <- integrate_expression(
  overlap_df = overlap_df_diff,
  rnaseq = data_list$rnaseq
)



CNVgene_rnaseq_del <- calculate_dosage(CNVgene_rnaseq)
