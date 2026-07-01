#!/usr/bin/env python3
import csv
from pathlib import Path


BASE = Path("/path/to/evaluation")
VALID = {"GAIN", "LOSS"}
OVERLAP_THRESHOLD = 0.5


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


def add_event(events, chrom, start, end, state):
    if end > start and state in VALID:
        events.append(
            {
                "chrom": chrom,
                "start": start,
                "end": end,
                "state": state,
                "length": end - start,
            }
        )


def benchmark_state(Control, Treat):
    control = float(Control)
    target = float(Treat)
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
    events = []
    with (BASE / "CNVbenchmark.csv").open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            add_event(
                events,
                norm_chrom(row["chr"]),
                to_int(row["BinLowEdge"]),
                to_int(row["BinUpEdge"]),
                benchmark_state(row["Control"], row["Treat"]),
            )
    return events


def load_bicseq2():
    events = []
    with first_existing("BICseq2.tsv", "bicseq2.tsv").open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        status_field = reader.fieldnames[-1]
        for row in reader:
            add_event(
                events,
                norm_chrom(row["chrom"]),
                to_int(row["start"]),
                to_int(row["end"]),
                normalize_state(row[status_field]),
            )
    return events


def load_cnvkit():
    events = []
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
            add_event(
                events,
                norm_chrom(row["chromosome"]),
                to_int(row["start"]),
                to_int(row["end"]),
                state,
            )
    return events


def load_controlfreec():
    events = []
    with first_existing("ControlFreec_CNVs", "Controlfreec_CNVs", "controlfreec_CNVs").open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if row:
                add_event(
                    events,
                    norm_chrom(row[0]),
                    to_int(row[1]),
                    to_int(row[2]),
                    normalize_state(row[4]),
                )
    return events


def load_ichor():
    events = []
    with (BASE / "ichorDNA.seg").open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            add_event(
                events,
                norm_chrom(row["chr"]),
                to_int(row["start"]),
                to_int(row["end"]),
                normalize_state(row["event"]),
            )
    return events


def load_gmm():
    events = []
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
            add_event(
                events,
                norm_chrom(row["chr"]),
                to_int(row["start_pos"]),
                to_int(row["end_pos"]),
                state,
            )
    return events


def overlap_len(a, b):
    if a["chrom"] != b["chrom"] or a["state"] != b["state"]:
        return 0
    return max(0, min(a["end"], b["end"]) - max(a["start"], b["start"]))


def is_match(call, benchmark):
    overlap = overlap_len(call, benchmark)
    if overlap <= 0:
        return False
    return (
        overlap / call["length"] > OVERLAP_THRESHOLD
        or overlap / benchmark["length"] > OVERLAP_THRESHOLD
    )


def evaluate(tool, calls, benchmarks):
    matched_benchmark = set()
    matched_calls = set()

    for call_index, call in enumerate(calls):
        for benchmark_index, benchmark in enumerate(benchmarks):
            if is_match(call, benchmark):
                matched_calls.add(call_index)
                matched_benchmark.add(benchmark_index)

    recall = len(matched_benchmark) / len(benchmarks) if benchmarks else 0
    precision = len(matched_calls) / len(calls) if calls else 0
    fdr = 1 - precision if calls else 0
    f1 = 2 * recall * precision / (recall + precision) if recall + precision else 0

    return {
        "tool": tool,
        "benchmark_events": len(benchmarks),
        "called_events": len(calls),
        "matched_benchmark_events": len(matched_benchmark),
        "matched_called_events": len(matched_calls),
        "recall_percent": f"{recall * 100:.2f}",
        "precision_percent": f"{precision * 100:.2f}",
        "fdr_percent": f"{fdr * 100:.2f}",
        "f1_score": f"{f1:.4f}",
    }


def main():
    benchmarks = load_benchmark()
    tools = {
        "BICseq2": load_bicseq2(),
        "CNVkit": load_cnvkit(),
        "ControlFreec": load_controlfreec(),
        "ichorDNA": load_ichor(),
        "GMM": load_gmm(),
    }
    rows = [evaluate(tool, calls, benchmarks) for tool, calls in tools.items()]
    out = BASE / "tool_event_recall_precision_fdr_f1.tsv"
    with out.open("w", newline="") as handle:
        fieldnames = [
            "tool",
            "benchmark_events",
            "called_events",
            "matched_benchmark_events",
            "matched_called_events",
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
