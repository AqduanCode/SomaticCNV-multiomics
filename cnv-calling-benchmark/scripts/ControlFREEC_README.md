# Control-FREEC two-sample CNV workflow

This workflow detects copy-number differences in a test sample relative to a control sample using Control-FREEC with 10-kb windows.

Official repository: <https://github.com/BoevaLab/FREEC>

## Requirements

Required commands:

```text
Bash
samtools
freec
```

Required inputs:

```text
coordinate-sorted test BAM
coordinate-sorted control BAM
matching hg38 reference FASTA
```

The BAM files and reference must use the same genome build and chromosome names.

## Run Control-FREEC

Copy the following code into a Bash script and edit only the variables at the top:

```bash
#!/usr/bin/env bash
set -euo pipefail

TEST_SAMPLE=TEST_SAMPLE
CONTROL_SAMPLE=CONTROL_SAMPLE
TEST_BAM=/path/to/TEST_SAMPLE.sorted.bam
CONTROL_BAM=/path/to/CONTROL_SAMPLE.sorted.bam
REFERENCE_FASTA=/path/to/hg38.fa
OUTDIR=/path/to/ControlFREEC_results
FREEC=freec

WINDOW=10000
THREADS=8
CHR_FASTA_DIR="$OUTDIR/chromosome_fasta"
CONFIG="$OUTDIR/controlfreec_10kb.conf"
LOG="$OUTDIR/controlfreec.log"
TEST_BASENAME="$(basename "$TEST_BAM")"

mkdir -p "$OUTDIR" "$CHR_FASTA_DIR"

command -v samtools >/dev/null 2>&1
command -v "$FREEC" >/dev/null 2>&1

[[ -s "$TEST_BAM" ]] || { echo "ERROR: missing test BAM: $TEST_BAM" >&2; exit 1; }
[[ -s "$CONTROL_BAM" ]] || { echo "ERROR: missing control BAM: $CONTROL_BAM" >&2; exit 1; }
[[ -s "$REFERENCE_FASTA" ]] || { echo "ERROR: missing reference FASTA: $REFERENCE_FASTA" >&2; exit 1; }

[[ -s "${REFERENCE_FASTA}.fai" ]] || samtools faidx "$REFERENCE_FASTA"
[[ -s "${TEST_BAM}.bai" || -s "${TEST_BAM%.bam}.bai" ]] || samtools index "$TEST_BAM"
[[ -s "${CONTROL_BAM}.bai" || -s "${CONTROL_BAM%.bam}.bai" ]] || samtools index "$CONTROL_BAM"

# Control-FREEC expects separate chromosome FASTA files.
while IFS=$'\t' read -r chromosome length remainder; do
  chromosome_fasta="$CHR_FASTA_DIR/${chromosome}.fa"
  if [[ ! -s "$chromosome_fasta" ]]; then
    samtools faidx "$REFERENCE_FASTA" "$chromosome" > "$chromosome_fasta"
  fi
done < "${REFERENCE_FASTA}.fai"

# Generate the Control-FREEC configuration.
cat > "$CONFIG" <<EOF
[general]
chrLenFile = ${REFERENCE_FASTA}.fai
ploidy = 2
breakPointThreshold = 0.8
window = $WINDOW
chrFiles = $CHR_FASTA_DIR/
outputDir = $OUTDIR/
numberOfProcesses = $THREADS
BedGraphOutput = FALSE
breakPointType = 1

[sample]
mateFile = $TEST_BAM
inputFormat = bam
mateOrientation = 0

[control]
mateFile = $CONTROL_BAM
inputFormat = bam
mateOrientation = 0
EOF

"$FREEC" -conf "$CONFIG" > "$LOG" 2>&1

RESULT="$OUTDIR/${TEST_BASENAME}_CNVs"
RATIO="$OUTDIR/${TEST_BASENAME}_ratio.txt"

[[ -s "$RESULT" ]] || {
  echo "ERROR: Control-FREEC did not produce $RESULT; inspect $LOG" >&2
  exit 1
}

echo "Control-FREEC completed"
echo "CNV segments: $RESULT"
echo "Bin-level ratios: $RATIO"
echo "Run log: $LOG"
```

## Main outputs

```text
TEST_SAMPLE.sorted.bam_CNVs
TEST_SAMPLE.sorted.bam_ratio.txt
TEST_SAMPLE.sorted.bam_info.txt
TEST_SAMPLE.sorted.bam_sample.cpn
CONTROL_SAMPLE.sorted.bam_control.cpn
controlfreec.log
```

The principal segment-level result is `TEST_SAMPLE.sorted.bam_CNVs`. The configuration places the test BAM under `[sample]` and the baseline BAM under `[control]`; reversing them reverses the comparison direction.

The supplied `controlfreec_10kb.conf` file in this directory is an editable template. The Bash code above generates the same configuration automatically from the selected paths.
