# BIC-seq2 two-sample CNV workflow

This workflow detects copy-number differences in a test sample relative to a control sample using 10-kb bins.

Official documentation: <https://www.math.pku.edu.cn/teachers/xirb/downloads/software/BICseq2/BICseq2.html>

## Requirements

The following commands must be available in `PATH`:

```text
samtools
perl
awk
NBICseq-norm.pl
NBICseq-seg.pl
```

Required inputs:

```text
coordinate-sorted test BAM
coordinate-sorted control BAM
matching hg38 reference FASTA
```

The BAM files and reference must use the same chromosome names.

## Run BIC-seq2

Copy the following code into a Bash script and edit only the variables at the top:

```bash
#!/usr/bin/env bash
set -euo pipefail

TEST_SAMPLE=TEST_SAMPLE
CONTROL_SAMPLE=CONTROL_SAMPLE
TEST_BAM=/path/to/TEST_SAMPLE.sorted.bam
CONTROL_BAM=/path/to/CONTROL_SAMPLE.sorted.bam
REFERENCE_FASTA=/path/to/hg38.fa
OUTDIR=/path/to/bicseq2_workdir

READ_LEN=100
FRAG_SIZE=300
BIN_SIZE=10000
LAMBDA=2
CHROMS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22)

mkdir -p \
  "$OUTDIR/chrom_fasta" \
  "$OUTDIR/mappability_nonN" \
  "$OUTDIR/read_positions/$TEST_SAMPLE" \
  "$OUTDIR/read_positions/$CONTROL_SAMPLE" \
  "$OUTDIR/normalized_bins_10kb/$TEST_SAMPLE" \
  "$OUTDIR/normalized_bins_10kb/$CONTROL_SAMPLE" \
  "$OUTDIR/tmp_10kb/$TEST_SAMPLE" \
  "$OUTDIR/tmp_10kb/$CONTROL_SAMPLE" \
  "$OUTDIR/tmp_10kb/seg" \
  "$OUTDIR/results" \
  "$OUTDIR/logs"

[[ -s "${REFERENCE_FASTA}.fai" ]] || samtools faidx "$REFERENCE_FASTA"
[[ -s "${TEST_BAM}.bai" || -s "${TEST_BAM%.bam}.bai" ]] || samtools index "$TEST_BAM"
[[ -s "${CONTROL_BAM}.bai" || -s "${CONTROL_BAM%.bam}.bai" ]] || samtools index "$CONTROL_BAM"

# Prepare chromosome FASTA files, non-N masks, and read positions.
for chr in "${CHROMS[@]}"; do
  samtools faidx "$REFERENCE_FASTA" "$chr" > "$OUTDIR/chrom_fasta/${chr}.fa"

  perl -ne '
    BEGIN { $pos=0; $start=0; }
    next if /^>/;
    chomp;
    for my $base (split //) {
      $pos++;
      if ($base =~ /[Nn]/) {
        if ($start) { print $start,"\t",$pos-1,"\n"; $start=0; }
      } else {
        $start ||= $pos;
      }
    }
    END { print $start,"\t",$pos,"\n" if $start; }
  ' "$OUTDIR/chrom_fasta/${chr}.fa" > "$OUTDIR/mappability_nonN/${chr}.nonN.map"

  samtools view -f 64 -F 260 "$TEST_BAM" "$chr" \
    | awk '{print $4}' > "$OUTDIR/read_positions/$TEST_SAMPLE/${chr}.seq"

  samtools view -f 64 -F 260 "$CONTROL_BAM" "$chr" \
    | awk '{print $4}' > "$OUTDIR/read_positions/$CONTROL_SAMPLE/${chr}.seq"
done

# Create one normalization configuration for each sample.
make_norm_config() {
  sample="$1"
  config="$OUTDIR/config_norm_${sample}_10kb.tsv"

  printf 'chrom\tfa\tmappability\tread_positions\tnormalized_bins\n' > "$config"
  for chr in "${CHROMS[@]}"; do
    printf '%s\t%s\t%s\t%s\t%s\n' \
      "$chr" \
      "$OUTDIR/chrom_fasta/${chr}.fa" \
      "$OUTDIR/mappability_nonN/${chr}.nonN.map" \
      "$OUTDIR/read_positions/$sample/${chr}.seq" \
      "$OUTDIR/normalized_bins_10kb/$sample/${chr}.norm.bin" \
      >> "$config"
  done
}

make_norm_config "$TEST_SAMPLE"
make_norm_config "$CONTROL_SAMPLE"

# Normalize test and control independently.
NBICseq-norm.pl \
  -l="$READ_LEN" -s="$FRAG_SIZE" -b="$BIN_SIZE" \
  --tmp="$OUTDIR/tmp_10kb/$TEST_SAMPLE" \
  "$OUTDIR/config_norm_${TEST_SAMPLE}_10kb.tsv" \
  "$OUTDIR/results/${TEST_SAMPLE}.10kb.norm_params.txt" \
  > "$OUTDIR/logs/${TEST_SAMPLE}.10kb.norm.log" 2>&1

NBICseq-norm.pl \
  -l="$READ_LEN" -s="$FRAG_SIZE" -b="$BIN_SIZE" \
  --tmp="$OUTDIR/tmp_10kb/$CONTROL_SAMPLE" \
  "$OUTDIR/config_norm_${CONTROL_SAMPLE}_10kb.tsv" \
  "$OUTDIR/results/${CONTROL_SAMPLE}.10kb.norm_params.txt" \
  > "$OUTDIR/logs/${CONTROL_SAMPLE}.10kb.norm.log" 2>&1

# Tell NBICseq-seg.pl which normalized bins are test and control.
BICSEQ2_SEG_CONFIG="$OUTDIR/config_seg_${TEST_SAMPLE}_vs_${CONTROL_SAMPLE}_10kb.tsv"
printf 'chrom\tcase_normalized_bins\tcontrol_normalized_bins\n' > "$BICSEQ2_SEG_CONFIG"

for chr in "${CHROMS[@]}"; do
  printf '%s\t%s\t%s\n' \
    "$chr" \
    "$OUTDIR/normalized_bins_10kb/$TEST_SAMPLE/${chr}.norm.bin" \
    "$OUTDIR/normalized_bins_10kb/$CONTROL_SAMPLE/${chr}.norm.bin" \
    >> "$BICSEQ2_SEG_CONFIG"
done

# Final CNV calling step.
NBICseq-seg.pl \
  --control \
  --lambda="$LAMBDA" \
  --tmp="$OUTDIR/tmp_10kb/seg" \
  "$BICSEQ2_SEG_CONFIG" \
  "$OUTDIR/results/${TEST_SAMPLE}_vs_${CONTROL_SAMPLE}.10kb.bicseq2_segments.tsv"
```

## Main output

```text
results/TEST_SAMPLE_vs_CONTROL_SAMPLE.10kb.bicseq2_segments.tsv
```

The second column of the segmentation configuration is the test/case sample; the third column is the control. Reversing them reverses the gain/loss direction.

