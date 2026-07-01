# ichorCNA two-sample CNV workflow

This workflow detects copy-number differences in a test sample relative to a control sample using 10-kb bins. It converts both BAM files to fixed-step WIG files and supplies the control WIG to ichorCNA through `--NORMWIG`.

Official repository: <https://github.com/broadinstitute/ichorCNA>

## Requirements

Required commands:

```text
Bash
samtools
awk
Rscript
```

The ichorCNA R package and its dependencies must already be installed. The ichorCNA directory must contain:

```text
scripts/runIchorCNA.R
inst/extdata/gc_hg38_10kb.wig
inst/extdata/map_hg38_10kb.wig
inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt
```

Required input data:

```text
coordinate-sorted test BAM
coordinate-sorted control BAM
matching BAM indexes
```

## Run ichorCNA

Copy the following code into a Bash script and edit only the variables at the top:

```bash
#!/usr/bin/env bash
set -euo pipefail

TEST_SAMPLE=TEST_SAMPLE
CONTROL_SAMPLE=CONTROL_SAMPLE
TEST_BAM=/path/to/TEST_SAMPLE.sorted.bam
CONTROL_BAM=/path/to/CONTROL_SAMPLE.sorted.bam
ICHORCNA_DIR=/path/to/ichorCNA
OUTDIR=/path/to/ichorCNA_results

BIN_SIZE=10000
RUN_ID="${TEST_SAMPLE}_vs_${CONTROL_SAMPLE}_10kb"
SIZES="$OUTDIR/hg38_chr1-22.chrom.sizes"
TEST_WIG="$OUTDIR/${TEST_SAMPLE}.reads.10kb.wig"
CONTROL_WIG="$OUTDIR/${CONTROL_SAMPLE}.reads.10kb.wig"
GC_WIG="$ICHORCNA_DIR/inst/extdata/gc_hg38_10kb.wig"
MAP_WIG="$ICHORCNA_DIR/inst/extdata/map_hg38_10kb.wig"
CENTROMERE_FILE="$ICHORCNA_DIR/inst/extdata/GRCh38.GCA_000001405.2_centromere_acen.txt"
LOG="$OUTDIR/${RUN_ID}.ichorCNA.log"

mkdir -p "$OUTDIR"

command -v samtools >/dev/null 2>&1
command -v awk >/dev/null 2>&1
command -v Rscript >/dev/null 2>&1

[[ -s "$TEST_BAM" ]] || { echo "ERROR: missing test BAM: $TEST_BAM" >&2; exit 1; }
[[ -s "$CONTROL_BAM" ]] || { echo "ERROR: missing control BAM: $CONTROL_BAM" >&2; exit 1; }
[[ -s "$ICHORCNA_DIR/scripts/runIchorCNA.R" ]] || { echo "ERROR: runIchorCNA.R not found" >&2; exit 1; }
[[ -s "$GC_WIG" ]] || { echo "ERROR: hg38 10-kb GC WIG not found" >&2; exit 1; }
[[ -s "$MAP_WIG" ]] || { echo "ERROR: hg38 10-kb mappability WIG not found" >&2; exit 1; }
[[ -s "$CENTROMERE_FILE" ]] || { echo "ERROR: hg38 centromere file not found" >&2; exit 1; }

[[ -s "${TEST_BAM}.bai" || -s "${TEST_BAM%.bam}.bai" ]] || samtools index "$TEST_BAM"
[[ -s "${CONTROL_BAM}.bai" || -s "${CONTROL_BAM%.bam}.bai" ]] || samtools index "$CONTROL_BAM"

# Use autosomes only, matching the CNV-Sim benchmark.
samtools idxstats "$CONTROL_BAM" \
  | awk '$1 ~ /^chr([1-9]|1[0-9]|2[0-2])$/ {print $1"\t"$2}' \
  > "$SIZES"

[[ -s "$SIZES" ]] || {
  echo "ERROR: no chr1-chr22 entries found; check chromosome naming" >&2
  exit 1
}

# Convert one BAM into a fixed-step read-count WIG.
make_wig() {
  local bam="$1"
  local wig="$2"
  local temporary="${wig}.tmp.$$"

  samtools view -F 3844 "$bam" $(cut -f1 "$SIZES") \
    | awk -v BS="$BIN_SIZE" '
        BEGIN {
          while ((getline line < ARGV[1]) > 0) {
            split(line, field, "\t")
            count++
            chromosome[count]=field[1]
            bins[field[1]]=int((field[2]+BS-1)/BS)
          }
          ARGV[1]=""
        }
        {
          bin=int(($4-1)/BS)+1
          reads[$3 SUBSEP bin]++
        }
        END {
          for (i=1; i<=count; i++) {
            chr=chromosome[i]
            print "fixedStep chrom=" chr " start=1 step=" BS " span=" BS
            for (bin=1; bin<=bins[chr]; bin++) {
              print reads[chr SUBSEP bin]+0
            }
          }
        }
      ' "$SIZES" - > "$temporary"

  mv "$temporary" "$wig"
}

make_wig "$TEST_BAM" "$TEST_WIG"
make_wig "$CONTROL_BAM" "$CONTROL_WIG"

Rscript "$ICHORCNA_DIR/scripts/runIchorCNA.R" \
  --id "$RUN_ID" \
  --WIG "$TEST_WIG" \
  --NORMWIG "$CONTROL_WIG" \
  --gcWig "$GC_WIG" \
  --mapWig "$MAP_WIG" \
  --normalPanel None \
  --ploidy "c(2,3)" \
  --normal "c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)" \
  --maxCN 7 \
  --includeHOMD FALSE \
  --chrs "c(1:22)" \
  --chrTrain "c(1:22)" \
  --chrNormalize "c(1:22)" \
  --estimateNormal TRUE \
  --estimatePloidy TRUE \
  --estimateScPrevalence TRUE \
  --scStates "c(1,3)" \
  --centromere "$CENTROMERE_FILE" \
  --genomeBuild hg38 \
  --genomeStyle UCSC \
  --exons.bed None \
  --txnE 0.9999999 \
  --txnStrength 10000000 \
  --plotFileType pdf \
  --plotYLim "c(-2,4)" \
  --outDir "$OUTDIR" \
  > "$LOG" 2>&1

RESULT="$OUTDIR/${RUN_ID}.seg"
[[ -s "$RESULT" ]] || {
  echo "ERROR: ichorCNA did not produce $RESULT; inspect $LOG" >&2
  exit 1
}

echo "ichorCNA completed"
echo "Segment result: $RESULT"
echo "Run log: $LOG"
```

## Main outputs

```text
TEST_SAMPLE_vs_CONTROL_SAMPLE_10kb.seg
TEST_SAMPLE_vs_CONTROL_SAMPLE_10kb.cna.seg
TEST_SAMPLE_vs_CONTROL_SAMPLE_10kb.params.txt
TEST_SAMPLE_vs_CONTROL_SAMPLE_10kb.correctedDepth.txt
TEST_SAMPLE.reads.10kb.wig
CONTROL_SAMPLE.reads.10kb.wig
```

The `.seg` file contains the main segment-level calls. Common event labels include `NEUT`, `HETD`, `HOMD`, `GAIN`, `AMP`, and `HLAMP`.

All GC, mappability, centromere, BAM, and WIG resources must use hg38 and compatible chromosome names. The result describes the test sample relative to the control WIG; the control is not automatically a biological matched normal.
