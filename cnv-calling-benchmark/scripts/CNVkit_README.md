# CNVkit two-sample CNV workflow

This workflow detects copy-number differences in a test sample relative to a control sample using CNVkit in WGS mode with 10-kb target bins.

Official documentation: <https://cnvkit.readthedocs.io/>

## Requirements

Required commands:

```text
Bash
samtools
cnvkit.py
```

Required inputs:

```text
coordinate-sorted test BAM
coordinate-sorted control BAM
matching hg38 reference FASTA
```

The BAM files and reference must use the same genome build and chromosome names.

## Run CNVkit

Copy the following code into a Bash script and edit only the variables at the top:

```bash
#!/usr/bin/env bash
set -euo pipefail

TEST_SAMPLE=TEST_SAMPLE
CONTROL_SAMPLE=CONTROL_SAMPLE
TEST_BAM=/path/to/TEST_SAMPLE.sorted.bam
CONTROL_BAM=/path/to/CONTROL_SAMPLE.sorted.bam
REFERENCE_FASTA=/path/to/hg38.fa
OUTDIR=/path/to/CNVkit_results
CNVKIT=cnvkit.py

BIN_SIZE=10000
PROCESSES=8
ACCESS="$OUTDIR/hg38.access.bed"
REFERENCE_CNN="$OUTDIR/${CONTROL_SAMPLE}_normal_reference.10kb.cnn"
LOG="$OUTDIR/${TEST_SAMPLE}_vs_${CONTROL_SAMPLE}.CNVkit.log"
TEST_BASENAME="$(basename "$TEST_BAM" .bam)"

mkdir -p "$OUTDIR"

command -v samtools >/dev/null 2>&1
command -v "$CNVKIT" >/dev/null 2>&1

[[ -s "$TEST_BAM" ]] || { echo "ERROR: missing test BAM: $TEST_BAM" >&2; exit 1; }
[[ -s "$CONTROL_BAM" ]] || { echo "ERROR: missing control BAM: $CONTROL_BAM" >&2; exit 1; }
[[ -s "$REFERENCE_FASTA" ]] || { echo "ERROR: missing reference FASTA: $REFERENCE_FASTA" >&2; exit 1; }

[[ -s "${REFERENCE_FASTA}.fai" ]] || samtools faidx "$REFERENCE_FASTA"
[[ -s "${TEST_BAM}.bai" || -s "${TEST_BAM%.bam}.bai" ]] || samtools index "$TEST_BAM"
[[ -s "${CONTROL_BAM}.bai" || -s "${CONTROL_BAM%.bam}.bai" ]] || samtools index "$CONTROL_BAM"

# Generate accessible genomic regions from the reference.
"$CNVKIT" access "$REFERENCE_FASTA" -o "$ACCESS"

# Build a control reference and call CNVs in the test sample.
"$CNVKIT" batch "$TEST_BAM" \
  --normal "$CONTROL_BAM" \
  --method wgs \
  --segment-method haar \
  --target-avg-size "$BIN_SIZE" \
  --fasta "$REFERENCE_FASTA" \
  --access "$ACCESS" \
  --output-reference "$REFERENCE_CNN" \
  --output-dir "$OUTDIR" \
  --scatter \
  --diagram \
  --processes "$PROCESSES" \
  > "$LOG" 2>&1

RESULT="$OUTDIR/${TEST_BASENAME}.call.cns"
[[ -s "$RESULT" ]] || {
  echo "ERROR: CNVkit did not produce $RESULT; inspect $LOG" >&2
  exit 1
}

cut -f1,2,3,5,8,11 "$RESULT" \
  > "$OUTDIR/${TEST_SAMPLE}_vs_${CONTROL_SAMPLE}.selected_columns.tsv"

cut -f1,2,3,8 "$RESULT" \
  > "$OUTDIR/${TEST_SAMPLE}_vs_${CONTROL_SAMPLE}.copy_number.tsv"

echo "CNVkit completed"
echo "Segment result: $RESULT"
echo "Run log: $LOG"
```

## Main outputs

```text
TEST_SAMPLE.sorted.cnr
TEST_SAMPLE.sorted.cns
TEST_SAMPLE.sorted.call.cns
TEST_SAMPLE.sorted-scatter.png
TEST_SAMPLE.sorted-diagram.pdf
CONTROL_SAMPLE_normal_reference.10kb.cnn
TEST_SAMPLE_vs_CONTROL_SAMPLE.selected_columns.tsv
TEST_SAMPLE_vs_CONTROL_SAMPLE.copy_number.tsv
```

The principal segment-level result is `TEST_SAMPLE.sorted.call.cns`. Positive `log2` values indicate relative gains in the test sample; negative values indicate relative losses. These are test-versus-control changes rather than independent absolute copy-number measurements.
