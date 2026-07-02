#!/usr/bin/env python3
import csv
from collections import defaultdict
from pathlib import Path


BASE = Path("/path/to/evaluation")
VALID = {"GAIN", "LOSS"}


def first_existing(*names):
    for name in names:
        path = BASE / name
        if path.exists():
            return path
    raise FileNotFoundError(", ".join(names))


def norm_chrom(value):
    chrom = str(value).strip().strip('"').lstrip("\ufeff")
    return chrom if chrom.lower().startswith("chr") else "chr" + chrom


def to_int(value):
    return int(float(str(value).strip().strip('"')))


def add_interval(store, chrom, start, end, state):
    if end > start and state in VALID:
        store.append((chrom, start, end, state))


def benchmark_state(Control_CN, Treat_CN):
    control = float(Control_CN)
    target = float(Treat_CN)
    if target > control:
        return "GAIN"
    if target < control:
        return "LOSS"
    return "Neutral"


def normalize_state(value):
    raw = str(value).strip().upper()
    if raw in {"GAIN", "AMP", "HLAMP", "DUP", "DUPLICATION", "AMPLIFICATION"}:
        return "GAIN"
    if raw in {"LOSS", "DEL", "HETD", "HOMD", "DELETION"}:
        return "LOSS"
    if raw in {"NEUT", "NEUTRAL", "NORMAL"}:
        return "Neutral"
    return "NA"


def load_benchmark():
    intervals = []
    with (BASE / "CNVbenchmark.csv").open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            add_interval(
                intervals,
                norm_chrom(row["chr"]),
                to_int(row["Start"]),
                to_int(row["End"]),
                benchmark_state(row["Control_CN"], row["Treat_CN"]),
            )
    return intervals


def load_bicseq2():
    intervals = []
    with first_existing("BICseq2.tsv", "bicseq2.tsv").open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        status_field = reader.fieldnames[-1]
        for row in reader:
            state = normalize_state(row[status_field])
            add_interval(intervals, norm_chrom(row["chrom"]), to_int(row["start"]), to_int(row["end"]), state)
    return intervals


def load_cnvkit():
    intervals = []
    with first_existing("CNVkit.tsv", "CNVKIT.tsv").open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if "CNV_status" in row:
                state = normalize_state(row["CNV_status"])
            else:
                cn = float(row["cn"])
                if cn > 2:
                    state = "GAIN"
                elif cn < 2:
                    state = "LOSS"
                else:
                    state = "Neutral"
            add_interval(intervals, norm_chrom(row["chromosome"]), to_int(row["start"]), to_int(row["end"]), state)
    return intervals


def load_controlfreec():
    intervals = []
    with first_existing("ControlFreec_CNVs", "Controlfreec_CNVs", "controlfreec_CNVs").open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if row:
                add_interval(intervals, norm_chrom(row[0]), to_int(row[1]), to_int(row[2]), normalize_state(row[4]))
    return intervals


def load_ichor():
    intervals = []
    with (BASE / "ichorDNA.seg").open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            add_interval(
                intervals,
                norm_chrom(row["chr"]),
                to_int(row["start"]),
                to_int(row["end"]),
                normalize_state(row["event"]),
            )
    return intervals


def load_gmm():
    intervals = []
    with (BASE / "GMM.csv").open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        last_field = reader.fieldnames[-1]
        for row in reader:
            value = float(row[last_field])
            if value > 0:
                state = "LOSS"
            elif value < 0:
                state = "GAIN"
            else:
                state = "Neutral"
            add_interval(intervals, norm_chrom(row["chr"]), to_int(row["start_pos"]), to_int(row["end_pos"]), state)
    return intervals


def merge_intervals(intervals):
    by_key = defaultdict(list)
    for chrom, start, end, state in intervals:
        by_key[(chrom, state)].append((start, end))
    merged = defaultdict(list)
    for key, ranges in by_key.items():
        ranges.sort()
        for start, end in ranges:
            if not merged[key] or start > merged[key][-1][1]:
                merged[key].append([start, end])
            elif end > merged[key][-1][1]:
                merged[key][-1][1] = end
    return merged


def total_length(merged):
    return sum(end - start for ranges in merged.values() for start, end in ranges)


def overlap_one(start, end, targets):
    overlap = 0
    for target_start, target_end in targets:
        if target_end <= start:
            continue
        if target_start >= end:
            break
        value = min(end, target_end) - max(start, target_start)
        if value > 0:
            overlap += value
    return overlap


def same_state_overlap(query_merged, target_merged):
    overlap = 0
    for key, ranges in query_merged.items():
        targets = target_merged.get(key, [])
        for start, end in ranges:
            overlap += overlap_one(start, end, targets)
    return overlap


def f1_score(recall, precision):
    return 2 * recall * precision / (recall + precision) if recall + precision else 0


def evaluate(tool, calls, benchmark_merged, benchmark_length):
    call_merged = merge_intervals(calls)
    called_length = total_length(call_merged)
    tp = same_state_overlap(call_merged, benchmark_merged)
    recall = tp / benchmark_length if benchmark_length else 0
    precision = tp / called_length if called_length else 0
    fdr = 1 - precision if called_length else 0
    return {
        "tool": tool,
        "benchmark_gain_loss_length_bp": benchmark_length,
        "called_gain_loss_length_bp": called_length,
        "intersection_with_benchmark_length_bp": tp,
        "recall_percent": f"{recall * 100:.2f}",
        "precision_percent": f"{precision * 100:.2f}",
        "fdr_percent": f"{fdr * 100:.2f}",
        "f1_score": f"{f1_score(recall, precision):.4f}",
    }


def main():
    benchmark = load_benchmark()
    benchmark_merged = merge_intervals(benchmark)
    benchmark_length = total_length(benchmark_merged)
    tools = {
        "BICseq2": load_bicseq2(),
        "CNVkit": load_cnvkit(),
        "ControlFreec": load_controlfreec(),
        "ichorDNA": load_ichor(),
        "GMM": load_gmm(),
    }
    rows = [evaluate(tool, calls, benchmark_merged, benchmark_length) for tool, calls in tools.items()]
    out = BASE / "tool_overall_recall_precision_fdr_f1.tsv"
    with out.open("w", newline="") as handle:
        fieldnames = [
            "tool",
            "benchmark_gain_loss_length_bp",
            "called_gain_loss_length_bp",
            "intersection_with_benchmark_length_bp",
            "recall_percent",
            "precision_percent",
            "fdr_percent",
            "f1_score",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    print(out)
    print(out.read_text())


if __name__ == "__main__":
    main()
