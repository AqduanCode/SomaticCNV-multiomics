# BIC-seq2 two-sample CNV workflow

This document describes a complete, reusable BIC-seq2 workflow for detecting copy-number differences in a test sample relative to a control sample. Official BIC-seq2 documentation: <https://www.math.pku.edu.cn/teachers/xirb/downloads/software/BICseq2/BICseq2.html>

The workflow has two required BIC-seq2 stages:

1. `NBICseq-norm.pl` normalizes the test and control samples independently.
2. `NBICseq-seg.pl --control` compares the two sets of normalized bins and calls CNV segments.


## 1. Requirements

Required software:

```text
Bash
SAMtools
Perl
awk
NBICseq-norm.pl
NBICseq-seg.pl
```

Make the BIC-seq2 executables available in `PATH`, or provide the normalization executable through `NBICSEQ_NORM`:

```bash
export PATH=/path/to/bicseq2/bin:$PATH
export NBICSEQ_NORM=/path/to/NBICseq-norm.pl

command -v samtools
command -v NBICseq-norm.pl
command -v NBICseq-seg.pl
```

Input data:

```text
coordinate-sorted test BAM and BAM index
coordinate-sorted control BAM and BAM index
matching reference FASTA and FASTA index
```

All inputs must use the same genome build and chromosome naming convention. For example, `chr1` in the BAM must correspond to `>chr1` in the FASTA.

## 2. Define the analysis

Set the repository, data, and output paths:

```bash
REPOSITORY=/path/to/cnv-calling-benchmark
SCRIPT_DIR="$REPOSITORY/scripts"

TEST_SAMPLE=TEST_SAMPLE
CONTROL_SAMPLE=CONTROL_SAMPLE

TEST_BAM=/path/to/TEST_SAMPLE.sorted.bam
CONTROL_BAM=/path/to/CONTROL_SAMPLE.sorted.bam
REFERENCE_FASTA=/path/to/hg38.fa

OUTDIR=/path/to/bicseq2_workdir
mkdir -p "$OUTDIR"
```

Select the chromosomes. The following list matches the original server workflow:

```bash
CHROMS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY"
export CHROMS
```

For an autosome-only benchmark, use the same chromosome list in every step:

```bash
CHROMS="chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22"
export CHROMS
```

Index the inputs if needed:

```bash
[[ -s "${REFERENCE_FASTA}.fai" ]] || samtools faidx "$REFERENCE_FASTA"
[[ -s "${TEST_BAM}.bai" || -s "${TEST_BAM%.bam}.bai" ]] || samtools index "$TEST_BAM"
[[ -s "${CONTROL_BAM}.bai" || -s "${CONTROL_BAM%.bam}.bai" ]] || samtools index "$CONTROL_BAM"
```

## 3. Normalize the test sample

Run the supplied preparation script:

```bash
bash "$SCRIPT_DIR/prepare_bicseq2_normalized_bins_10kb.sh" \
  "$TEST_SAMPLE" \
  "$TEST_BAM" \
  "$REFERENCE_FASTA" \
  "$OUTDIR"
```

The script performs the following operations:

1. Extracts one FASTA file per chromosome with `samtools faidx`.
2. Generates the non-N interval mask used by the original benchmark workflow.
3. Extracts mapped read-1 positions from the BAM into `chr*.seq` files.
4. Creates `config_norm_TEST_SAMPLE_10kb.tsv`.
5. Runs `NBICseq-norm.pl` with read length 100 bp, fragment size 300 bp, and bin size 10,000 bp.

Principal outputs:

```text
$OUTDIR/config_norm_TEST_SAMPLE_10kb.tsv
$OUTDIR/read_positions/TEST_SAMPLE/chr*.seq
$OUTDIR/normalized_bins_10kb/TEST_SAMPLE/chr*.norm.bin
$OUTDIR/results/TEST_SAMPLE.10kb.norm_params.txt
$OUTDIR/logs/TEST_SAMPLE.10kb.norm.log
```

The normalization configuration is tab-delimited and has five columns:

```text
chrom  fa  mappability  read_positions  normalized_bins
```

The fifth column tells `NBICseq-norm.pl` where to write each `chr*.norm.bin` file.

## 4. Normalize the control sample

Run the same script independently for the control:

```bash
bash "$SCRIPT_DIR/prepare_bicseq2_normalized_bins_10kb.sh" \
  "$CONTROL_SAMPLE" \
  "$CONTROL_BAM" \
  "$REFERENCE_FASTA" \
  "$OUTDIR"
```

Principal outputs:

```text
$OUTDIR/config_norm_CONTROL_SAMPLE_10kb.tsv
$OUTDIR/read_positions/CONTROL_SAMPLE/chr*.seq
$OUTDIR/normalized_bins_10kb/CONTROL_SAMPLE/chr*.norm.bin
$OUTDIR/results/CONTROL_SAMPLE.10kb.norm_params.txt
$OUTDIR/logs/CONTROL_SAMPLE.10kb.norm.log
```

A normalized-bin file normally begins with:

```text
start  end    obs  expected  var
10001  20000  506  551.049   561.452
```

Verify that both samples produced non-empty files:

```bash
for chr in $CHROMS; do
  test -s "$OUTDIR/normalized_bins_10kb/$TEST_SAMPLE/${chr}.norm.bin"
  test -s "$OUTDIR/normalized_bins_10kb/$CONTROL_SAMPLE/${chr}.norm.bin"
done
```

## 5. Create the test/control segmentation configuration

Generate the configuration from the actual output paths. This is safer than manually editing each line of the supplied template:

```bash
BICSEQ2_SEG_CONFIG="$OUTDIR/config_seg_${TEST_SAMPLE}_vs_${CONTROL_SAMPLE}_10kb.tsv"

{
  printf 'chrom\tcase_normalized_bins\tcontrol_normalized_bins\n'

  for chr in $CHROMS; do
    printf '%s\t%s\t%s\n' \
      "$chr" \
      "$OUTDIR/normalized_bins_10kb/$TEST_SAMPLE/${chr}.norm.bin" \
      "$OUTDIR/normalized_bins_10kb/$CONTROL_SAMPLE/${chr}.norm.bin"
  done
} > "$BICSEQ2_SEG_CONFIG"
```

The file must be tab-delimited:

```text
chrom  case_normalized_bins  control_normalized_bins
chr1   /path/to/TEST_SAMPLE/chr1.norm.bin   /path/to/CONTROL_SAMPLE/chr1.norm.bin
```

Column order is critical:

- column 2 is the test/case sample;
- column 3 is the control sample.

Reversing these columns reverses the comparison and therefore the gain/loss direction.

## 6. Call CNVs with BIC-seq2

Create the output and temporary directories:

```bash
mkdir -p "$OUTDIR/results" "$OUTDIR/tmp_10kb/seg"
```

Run the final segmentation step:

```bash
SEGMENT_OUTPUT="$OUTDIR/results/${TEST_SAMPLE}_vs_${CONTROL_SAMPLE}.10kb.bicseq2_segments.tsv"
SEGMENT_FIGURE="$OUTDIR/results/${TEST_SAMPLE}_vs_${CONTROL_SAMPLE}.10kb.bicseq2_profile.png"

NBICseq-seg.pl \
  --control \
  --lambda=2 \
  --tmp="$OUTDIR/tmp_10kb/seg" \
  --fig="$SEGMENT_FIGURE" \
  --title="$TEST_SAMPLE vs $CONTROL_SAMPLE BIC-seq2 CNV 10kb" \
  "$BICSEQ2_SEG_CONFIG" \
  "$SEGMENT_OUTPUT"
```

`--control` activates the case/control mode. It does not specify the BAM paths or sample roles; those roles are already encoded in columns 2 and 3 of `BICSEQ2_SEG_CONFIG`.

The principal result is:

```text
$OUTDIR/results/TEST_SAMPLE_vs_CONTROL_SAMPLE.10kb.bicseq2_segments.tsv
```

Typical columns include:

```text
chrom
start
end
binNum
tumor
tumor_expect
normal
normal_expect
log2.copyRatio
log2.TumorExpectRatio
```

## 7. Optional gain/loss annotation

The original benchmark classified relative copy-ratio values as:

```text
copyRatio > 1.1  => gain
copyRatio < 0.9  => loss
otherwise        => neutral
```

Add `copyRatio` and `status` columns:

```bash
ANNOTATED_OUTPUT="$OUTDIR/results/${TEST_SAMPLE}_vs_${CONTROL_SAMPLE}.10kb.bicseq2_segments.annotated.tsv"

awk 'BEGIN {
       OFS="\t"
       loss=log(0.9)/log(2)
       gain=log(1.1)/log(2)
     }
     NR==1 {
       for (i=1; i<=NF; i++) {
         if ($i=="log2.copyRatio") column=i
       }
       if (!column) {
         print "ERROR: log2.copyRatio column not found" > "/dev/stderr"
         exit 2
       }
       print $0,"copyRatio","status"
       next
     }
     {
       log2ratio=$column+0
       ratio=2^log2ratio
       if (log2ratio > gain) status="gain"
       else if (log2ratio < loss) status="loss"
       else status="neutral"
       print $0,ratio,status
     }' "$SEGMENT_OUTPUT" > "$ANNOTATED_OUTPUT"
```
