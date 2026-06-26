#!/usr/bin/env python3
"""
Per-site m6A annotation with transcript region and motif information.
"""

import csv
from collections import defaultdict
from pathlib import Path

import pysam

from .conservation.conserved_sites import (
    build_position_index,
    genomic_to_transcript_pos,
    parse_annotation,
)
from .motif import classify_motif_simple, extract_3mer_motif, extract_sequence_context

OUTPUT_COLUMNS = [
    "site_index",
    "chrom",
    "start",
    "end",
    "position",
    "bed_name",
    "ratio",
    "bed_strand",
    "assigned_strand",
    "motif",
    "motif_type",
    "sequence_context",
    "region",
    "gene_id",
    "gene_name",
    "transcript_id",
    "transcript_position",
    "overlap_count",
    "annotation_status",
]

REGION_PRIORITY = {
    "cds": 0,
    "5utr": 1,
    "3utr": 2,
    "noncoding": 3,
    "intron": 4,
    "intergenic": 5,
}


def infer_annotation_format(annotation_path, annotation_format="auto"):
    if annotation_format != "auto":
        return annotation_format
    path = str(annotation_path).lower()
    if path.endswith((".gff", ".gff3", ".gff.gz", ".gff3.gz")):
        return "gff"
    return "gtf"


def parse_ratio(value):
    if value in (None, ""):
        return ""
    try:
        return float(value)
    except ValueError:
        return ""


def iter_bed_sites(bed_path):
    with open(bed_path) as bed:
        site_index = 0
        for line in bed:
            if line.startswith("#") or line.startswith("track") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            start = int(float(fields[1]))
            end = int(float(fields[2]))
            yield {
                "site_index": site_index,
                "chrom": fields[0],
                "start": start,
                "end": end,
                "position": (start + end) // 2,
                "bed_name": fields[3] if len(fields) > 3 else "",
                "ratio": parse_ratio(fields[4]) if len(fields) > 4 else "",
                "bed_strand": fields[5] if len(fields) > 5 and fields[5] in ["+", "-"] else "",
            }
            site_index += 1


def build_binned_index(position_index, bin_size=100000):
    binned = defaultdict(lambda: defaultdict(list))
    for chrom, intervals in position_index.items():
        for interval in intervals:
            start, end = interval[0], interval[1]
            for bin_id in range(start // bin_size, end // bin_size + 1):
                binned[chrom][bin_id].append(interval)
    return binned


def candidate_intervals(site, binned_index, bin_size=100000):
    chrom = site["chrom"]
    pos = site["position"]
    intervals = binned_index.get(chrom, {}).get(pos // bin_size, [])
    return [iv for iv in intervals if iv[0] <= pos < iv[1]]


def classify_site_against_model(site, model):
    tx_pos, region = genomic_to_transcript_pos(site["position"], model)
    if tx_pos is not None:
        status = "noncoding" if region == "noncoding" else "annotated"
        return tx_pos, region, status
    return "", "intron", "intron"


def make_annotation_candidate(site, model, tx_pos, region, status, overlap_count):
    return {
        "assigned_strand": model.strand,
        "region": region,
        "gene_id": model.gene_id,
        "gene_name": model.gene_name,
        "transcript_id": model.tid,
        "transcript_position": tx_pos,
        "overlap_count": overlap_count,
        "annotation_status": status,
    }


def best_candidate(site, candidates):
    if not candidates:
        return None
    bed_strand = site["bed_strand"]

    def sort_key(candidate):
        strand_penalty = 0 if not bed_strand or candidate["assigned_strand"] == bed_strand else 1
        has_gene_name_penalty = 0 if candidate["gene_name"] else 1
        return (
            strand_penalty,
            REGION_PRIORITY.get(candidate["region"], 99),
            has_gene_name_penalty,
            candidate["gene_id"] or "",
            candidate["transcript_id"] or "",
        )

    return sorted(candidates, key=sort_key)[0]


def annotate_region(site, models, binned_index, overlap_policy="best"):
    intervals = candidate_intervals(site, binned_index)
    if not intervals:
        return [{
            "assigned_strand": site["bed_strand"] or "+",
            "region": "intergenic",
            "gene_id": "",
            "gene_name": "",
            "transcript_id": "",
            "transcript_position": "",
            "overlap_count": 0,
            "annotation_status": "intergenic",
        }]

    candidates = []
    overlap_count = len(intervals)
    for _, _, tid, _ in intervals:
        model = models[tid]
        tx_pos, region, status = classify_site_against_model(site, model)
        candidates.append(make_annotation_candidate(site, model, tx_pos, region, status, overlap_count))

    if overlap_policy == "all":
        return sorted(
            candidates,
            key=lambda item: (
                REGION_PRIORITY.get(item["region"], 99),
                item["gene_id"] or "",
                item["transcript_id"] or "",
            ),
        )
    return [best_candidate(site, candidates)]


def annotate_motif(site, genome, assigned_strand, motif_window=2):
    strand_order = []
    if assigned_strand in ["+", "-"]:
        strand_order.append(assigned_strand)
    if site["bed_strand"] in ["+", "-"] and site["bed_strand"] not in strand_order:
        strand_order.append(site["bed_strand"])
    for fallback in ["+", "-"]:
        if fallback not in strand_order:
            strand_order.append(fallback)

    for strand in strand_order:
        context = extract_sequence_context(
            genome, site["chrom"], site["position"], strand, window=motif_window
        )
        if not context:
            continue
        motif, is_valid = extract_3mer_motif(context)
        if is_valid:
            motif_type, motif_seq = classify_motif_simple(context)
            return strand, motif_seq, motif_type, context, ""
        motif_type, motif_seq = classify_motif_simple(context)
        return strand, motif_seq or "", motif_type, context, "motif_failed"

    return assigned_strand, "", "unknown", "", "motif_failed"


def write_site_table(rows, output_path):
    with open(output_path, "w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=OUTPUT_COLUMNS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def annotate_m6a_sites(
    genome_path,
    annotation_path,
    bed_path,
    output_path=None,
    annotation_format="auto",
    motif_window=2,
    overlap_policy="best",
    verbose=True,
):
    """Annotate every BED m6A site with region, motif, gene, strand, and ratio."""
    if overlap_policy not in ["best", "all"]:
        raise ValueError("overlap_policy must be 'best' or 'all'")

    genome_path = str(Path(genome_path).expanduser())
    annotation_path = str(Path(annotation_path).expanduser())
    bed_path = str(Path(bed_path).expanduser())
    fmt = infer_annotation_format(annotation_path, annotation_format)

    if verbose:
        print("=== Annotating m6A sites ===")
        print(f"Genome: {genome_path}")
        print(f"Annotation: {annotation_path} ({fmt})")
        print(f"BED: {bed_path}")

    models = parse_annotation(annotation_path, fmt=fmt)
    position_index = build_position_index(models)
    binned_index = build_binned_index(position_index)

    rows = []
    genome = pysam.FastaFile(genome_path)
    try:
        for site in iter_bed_sites(bed_path):
            annotations = annotate_region(site, models, binned_index, overlap_policy=overlap_policy)
            for annotation in annotations:
                motif_strand, motif, motif_type, context, motif_status = annotate_motif(
                    site, genome, annotation["assigned_strand"], motif_window=motif_window
                )
                status = annotation["annotation_status"]
                if motif_status and status == "annotated":
                    status = motif_status
                row = {
                    **site,
                    **annotation,
                    "assigned_strand": motif_strand or annotation["assigned_strand"],
                    "motif": motif,
                    "motif_type": motif_type,
                    "sequence_context": context,
                    "annotation_status": status,
                }
                rows.append(row)
    finally:
        genome.close()

    if output_path:
        write_site_table(rows, output_path)

    if verbose:
        print(f"Annotated rows: {len(rows):,}")
        if output_path:
            print(f"Output: {output_path}")

    return rows


__all__ = ["annotate_m6a_sites", "OUTPUT_COLUMNS"]
