#!/usr/bin/env python3
"""Extract one longest protein-coding transcript per gene from annotation and genome.

Protein and CDS sequences are extracted with gffread, then the longest protein
isoform per gene is selected. FASTA headers keep transcript IDs so OrthoFinder
pair tables can be passed directly to plantm6a conserve-pair.
"""

import argparse
import csv
import re
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path


def parse_gtf_attrs(attrs):
    parsed = {}
    for token in attrs.split(";"):
        token = token.strip()
        if not token:
            continue
        if " " in token:
            key, value = token.split(" ", 1)
            parsed[key] = value.strip().strip('"')
        elif "=" in token:
            key, value = token.split("=", 1)
            parsed[key] = value.strip().strip('"')
    return parsed


def parse_gff_attrs(attrs):
    parsed = {}
    for token in attrs.split(";"):
        token = token.strip()
        if not token or "=" not in token:
            continue
        key, value = token.split("=", 1)
        parsed[key] = value.strip().strip('"')
    return parsed


def infer_format(path, fmt):
    if fmt != "auto":
        return fmt
    lower = str(path).lower()
    return "gff" if lower.endswith((".gff", ".gff3", ".gff.gz", ".gff3.gz")) else "gtf"


def build_transcript_info(annotation_path, fmt):
    transcript_info = {}
    gff_transcript_to_gene = {}
    gff_transcript_to_name = {}

    with open(annotation_path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            feature = fields[2]
            attrs = fields[8]

            if fmt == "gtf":
                attr = parse_gtf_attrs(attrs)
                tid = attr.get("transcript_id")
                if not tid:
                    continue
                gene_id = attr.get("gene_id") or tid
                gene_name = attr.get("gene_name") or attr.get("gene") or gene_id
            else:
                attr = parse_gff_attrs(attrs)
                if feature in ("mRNA", "transcript"):
                    tid = attr.get("ID") or attr.get("transcript_id")
                    gene_id = attr.get("Parent") or attr.get("gene_id") or attr.get("gene") or tid
                    gene_name = attr.get("Name") or attr.get("gene_name") or gene_id
                    if tid:
                        gff_transcript_to_gene[tid] = gene_id
                        gff_transcript_to_name[tid] = gene_name
                    continue
                parent = attr.get("Parent", "").split(",")[0] if attr.get("Parent") else None
                tid = attr.get("transcript_id") or parent or attr.get("ID")
                if not tid:
                    continue
                gene_id = attr.get("gene_id") or attr.get("gene") or gff_transcript_to_gene.get(tid) or tid
                gene_name = attr.get("Name") or attr.get("gene_name") or gff_transcript_to_name.get(tid) or gene_id

            current = transcript_info.setdefault(
                tid,
                {
                    "transcript_id": tid,
                    "gene_id": gene_id,
                    "gene_name": gene_name,
                },
            )
            if gene_id and (not current["gene_id"] or current["gene_id"] == tid):
                current["gene_id"] = gene_id
            if gene_name and (not current["gene_name"] or current["gene_name"] == current["gene_id"]):
                current["gene_name"] = gene_name

    return transcript_info


def run_gffread(annotation_path, genome_path, cds_path, protein_path):
    cmd = [
        "gffread",
        str(annotation_path),
        "-g",
        str(genome_path),
        "-x",
        str(cds_path),
        "-y",
        str(protein_path),
    ]
    subprocess.run(cmd, check=True)


def parse_fasta(path):
    records = {}
    header = None
    seq = []

    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    rid = header[1:].split()[0]
                    records[rid] = {"id": rid, "header": header, "sequence": "".join(seq)}
                header = line
                seq = []
            else:
                seq.append(line)

    if header is not None:
        rid = header[1:].split()[0]
        records[rid] = {"id": rid, "header": header, "sequence": "".join(seq)}

    return records


def parse_header_ids(header):
    text = header.lstrip(">")
    tid = text.split()[0] if text.split() else text
    gene_id = None

    gene_match = re.search(r'gene_id[=\s]+"?([^";,\s]+)"?', text)
    if gene_match:
        gene_id = gene_match.group(1)
    if not gene_id:
        gene_match = re.search(r'gene[=:]([^\s;,]+)', text)
        if gene_match:
            gene_id = gene_match.group(1)
    if not gene_id:
        parent_match = re.search(r'Parent=([^;\s]+)', text)
        if parent_match:
            gene_id = parent_match.group(1)
    if not gene_id and tid:
        inferred = re.sub(r"\.\d+$", "", tid)
        inferred = re.sub(r"-[A-Z]+$", "", inferred)
        gene_id = inferred or tid

    return tid, gene_id


def choose_longest_per_gene(protein_records, cds_records, transcript_info, min_aa):
    by_gene = defaultdict(list)

    for tid, protein_record in protein_records.items():
        protein = protein_record["sequence"]
        if len(protein) < min_aa:
            continue
        if "*" in protein:
            continue

        parsed_tid, parsed_gene_id = parse_header_ids(protein_record["header"])
        info = transcript_info.get(tid) or transcript_info.get(parsed_tid) or {}
        gene_id = info.get("gene_id") or parsed_gene_id or tid
        gene_name = info.get("gene_name") or gene_id
        cds_seq = cds_records.get(tid, {}).get("sequence", "")

        by_gene[gene_id].append({
            "transcript_id": tid,
            "gene_id": gene_id,
            "gene_name": gene_name,
            "cds_seq": cds_seq,
            "protein": protein,
        })

    selected = []
    for records in by_gene.values():
        selected.append(sorted(records, key=lambda item: (-len(item["protein"]), item["transcript_id"]))[0])
    return sorted(selected, key=lambda item: item["transcript_id"])


def wrap(seq, width=60):
    for i in range(0, len(seq), width):
        yield seq[i:i + width]


def write_outputs(records, cds_out, protein_out, map_out, species):
    with open(cds_out, "w") as cds_handle, open(protein_out, "w") as protein_handle:
        for record in records:
            header = (
                f'>{record["transcript_id"]} gene={record["gene_id"]} '
                f'gene_name={record["gene_name"]} species={species}'
            )
            cds_handle.write(header + "\n")
            if record["cds_seq"]:
                cds_handle.write("\n".join(wrap(record["cds_seq"])) + "\n")
            protein_handle.write(header + "\n")
            protein_handle.write("\n".join(wrap(record["protein"])) + "\n")

    with open(map_out, "w", newline="") as map_handle:
        writer = csv.DictWriter(
            map_handle,
            fieldnames=["species", "gene_id", "gene_name", "transcript_id", "cds_length", "protein_length"],
            delimiter="\t",
        )
        writer.writeheader()
        for record in records:
            writer.writerow({
                "species": species,
                "gene_id": record["gene_id"],
                "gene_name": record["gene_name"],
                "transcript_id": record["transcript_id"],
                "cds_length": len(record["cds_seq"]),
                "protein_length": len(record["protein"]),
            })


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--species", required=True)
    parser.add_argument("--genome", required=True)
    parser.add_argument("--annotation", required=True)
    parser.add_argument("--format", choices=["auto", "gtf", "gff"], default="auto")
    parser.add_argument("--cds-out", required=True)
    parser.add_argument("--protein-out", required=True)
    parser.add_argument("--map-out", required=True)
    parser.add_argument("--min-aa", type=int, default=30)
    args = parser.parse_args()

    fmt = infer_format(args.annotation, args.format)
    transcript_info = build_transcript_info(args.annotation, fmt)

    Path(args.cds_out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.protein_out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.map_out).parent.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix="gffread_longest_") as tmpdir:
        all_cds = Path(tmpdir) / "all_cds.fa"
        all_proteins = Path(tmpdir) / "all_proteins.fa"
        run_gffread(args.annotation, args.genome, all_cds, all_proteins)
        cds_records = parse_fasta(all_cds)
        protein_records = parse_fasta(all_proteins)

    selected = choose_longest_per_gene(protein_records, cds_records, transcript_info, args.min_aa)
    write_outputs(selected, args.cds_out, args.protein_out, args.map_out, args.species)

    print(f"{args.species}: parsed {len(protein_records)} protein transcripts; selected {len(selected)} longest transcripts")
    print(f"CDS: {args.cds_out}")
    print(f"Protein: {args.protein_out}")
    print(f"Map: {args.map_out}")


if __name__ == "__main__":
    main()
