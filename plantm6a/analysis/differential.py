#!/usr/bin/env python3
"""Differential m6A analysis between two conditions of one species."""

from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

try:
    from intervaltree import IntervalTree
except ImportError:
    IntervalTree = None


class DifferentialGeneModel:
    __slots__ = [
        "gene_id",
        "transcript_id",
        "chrom",
        "strand",
        "tx_start",
        "tx_end",
        "exons",
        "cds_intervals",
    ]

    def __init__(self, gene_id, transcript_id, chrom, strand):
        self.gene_id = gene_id
        self.transcript_id = transcript_id
        self.chrom = chrom
        self.strand = strand
        self.tx_start = None
        self.tx_end = None
        self.exons = []
        self.cds_intervals = []

    def finalize(self):
        self.exons = sorted(set(self.exons))
        self.cds_intervals = sorted(set(self.cds_intervals))
        if self.exons:
            self.tx_start = self.exons[0][0]
            self.tx_end = self.exons[-1][1]

    def classify_position(self, pos):
        in_exon = any(start <= pos < end for start, end in self.exons)
        if not in_exon:
            return "intron"
        if not self.cds_intervals:
            return "noncoding"

        cds_start = self.cds_intervals[0][0]
        cds_end = self.cds_intervals[-1][1]
        in_cds = any(start <= pos < end for start, end in self.cds_intervals)
        if in_cds:
            return "cds"

        if self.strand == "+":
            return "5utr" if pos < cds_start else "3utr"
        return "5utr" if pos >= cds_end else "3utr"


def _parse_gtf_attrs(attrs):
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


def _parse_gff_attrs(attrs):
    parsed = {}
    for token in attrs.split(";"):
        token = token.strip()
        if not token or "=" not in token:
            continue
        key, value = token.split("=", 1)
        parsed[key] = value.strip().strip('"')
    return parsed


def parse_differential_annotation(annotation_path, verbose=True):
    """Parse GTF/GFF annotation for differential m6A site annotation."""
    annotation_path = str(annotation_path)
    is_gff3 = annotation_path.lower().endswith((".gff", ".gff3"))
    models = {}
    gene_to_transcripts = defaultdict(list)
    transcript_to_gene = {}

    if verbose:
        print(f"Parsing annotation: {annotation_path}")

    with open(annotation_path) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            feature = parts[2]
            if feature not in ("exon", "CDS", "mRNA", "transcript"):
                continue

            chrom = parts[0]
            start = int(parts[3]) - 1
            end = int(parts[4])
            strand = parts[6]
            attrs = parts[8]

            if is_gff3:
                attr = _parse_gff_attrs(attrs)
                if feature in ("mRNA", "transcript"):
                    tid = attr.get("ID") or attr.get("transcript_id")
                    gid = attr.get("Parent") or attr.get("gene_id") or tid
                    if tid:
                        transcript_to_gene[tid] = gid
                    continue
                parent = attr.get("Parent", "").split(",")[0] if attr.get("Parent") else None
                tid = attr.get("transcript_id") or parent or attr.get("ID")
                gid = attr.get("gene_id") or transcript_to_gene.get(tid) or tid
            else:
                attr = _parse_gtf_attrs(attrs)
                tid = attr.get("transcript_id")
                gid = attr.get("gene_id") or tid

            if not tid:
                continue
            if not gid:
                gid = tid.rsplit(".", 1)[0]

            if tid not in models:
                models[tid] = DifferentialGeneModel(gid, tid, chrom, strand)
                gene_to_transcripts[gid].append(tid)

            if feature == "exon":
                models[tid].exons.append((start, end))
            elif feature == "CDS":
                models[tid].cds_intervals.append((start, end))

    for model in models.values():
        model.finalize()
    models = {tid: model for tid, model in models.items() if model.exons}

    if verbose:
        print(f"  Transcripts: {len(models)}, Genes: {len(gene_to_transcripts)}")
    return models, build_differential_index(models)


def build_differential_index(models, bin_size=100000):
    if IntervalTree is not None:
        trees = defaultdict(lambda: IntervalTree())
        for tid, model in models.items():
            if model.tx_start is not None and model.tx_end is not None:
                trees[model.chrom][model.tx_start:model.tx_end] = tid
        return dict(trees)

    binned = defaultdict(lambda: defaultdict(list))
    for tid, model in models.items():
        if model.tx_start is None or model.tx_end is None:
            continue
        for bin_id in range(model.tx_start // bin_size, model.tx_end // bin_size + 1):
            binned[model.chrom][bin_id].append(tid)
    return binned


def load_bed6_sites(bed_path, verbose=True):
    """Load m6A BED6 sites using BED start as the m6A position."""
    if verbose:
        print(f"Loading BED: {bed_path}")
    sites = []
    with open(bed_path) as handle:
        for line in handle:
            if line.startswith("#") or line.startswith("track") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            try:
                ratio = float(parts[4])
            except ValueError:
                ratio = 0.0
            sites.append({
                "chrom": parts[0],
                "pos": int(parts[1]),
                "strand": parts[5],
                "ratio": ratio,
            })
    if verbose:
        print(f"  Sites loaded: {len(sites)}")
    return sites


def annotate_differential_site(chrom, pos, strand, models, binned_index, bin_size=100000):
    if IntervalTree is not None:
        if chrom not in binned_index:
            return None, None, "intergenic"
        overlaps = binned_index[chrom][pos]
        if not overlaps:
            return None, None, "intergenic"

        best_model = None
        best_len = 0
        for interval in overlaps:
            model = models[interval.data]
            if model.strand != strand:
                continue
            tx_len = model.tx_end - model.tx_start
            if tx_len > best_len:
                best_len = tx_len
                best_model = model

        if best_model is None:
            for interval in overlaps:
                model = models[interval.data]
                tx_len = model.tx_end - model.tx_start
                if tx_len > best_len:
                    best_len = tx_len
                    best_model = model
    else:
        overlaps = []
        for tid in binned_index.get(chrom, {}).get(pos // bin_size, []):
            model = models[tid]
            if model.tx_start <= pos < model.tx_end:
                overlaps.append(model)
        if not overlaps:
            return None, None, "intergenic"

        best_model = None
        best_len = 0
        for model in overlaps:
            if model.strand != strand:
                continue
            tx_len = model.tx_end - model.tx_start
            if tx_len > best_len:
                best_len = tx_len
                best_model = model

        if best_model is None:
            for model in overlaps:
                tx_len = model.tx_end - model.tx_start
                if tx_len > best_len:
                    best_len = tx_len
                    best_model = model

    if best_model is None:
        return None, None, "intergenic"
    return best_model.gene_id, best_model.transcript_id, best_model.classify_position(pos)


def annotate_differential_sites(sites, models, binned_index, verbose=True):
    """Annotate all m6A sites with gene, transcript, and region."""
    if verbose:
        print("Annotating sites...")
    for site in sites:
        gene_id, transcript_id, region = annotate_differential_site(
            site["chrom"], site["pos"], site["strand"], models, binned_index
        )
        site["gene_id"] = gene_id
        site["transcript_id"] = transcript_id
        site["region"] = region

    if verbose:
        annotated = sum(1 for site in sites if site["gene_id"] is not None)
        pct = annotated / len(sites) * 100 if sites else 0
        print(f"  Annotated: {annotated}/{len(sites)} ({pct:.1f}%)")
    return sites


def find_differential_m6a(wt_sites, mut_sites, ratio_threshold=0.1, fold_change_threshold=1.5, verbose=True):
    """Compare m6A profiles and classify each exact genomic position."""
    if verbose:
        print("Comparing m6A profiles...")

    wt_index = {(site["chrom"], site["pos"], site["strand"]): site for site in wt_sites}
    mut_index = {(site["chrom"], site["pos"], site["strand"]): site for site in mut_sites}
    all_keys = set(wt_index) | set(mut_index)
    results = []

    for chrom, pos, strand in sorted(all_keys):
        key = (chrom, pos, strand)
        in_wt = key in wt_index
        in_mut = key in mut_index
        wt_ratio = wt_index[key]["ratio"] if in_wt else 0.0
        mut_ratio = mut_index[key]["ratio"] if in_mut else 0.0

        source = wt_index[key] if in_wt else mut_index[key]
        gene_id = source.get("gene_id")
        transcript_id = source.get("transcript_id")
        region = source.get("region") or "intergenic"

        if in_wt and not in_mut:
            status = "wt_only"
        elif not in_wt and in_mut:
            status = "mut_only"
        else:
            diff = mut_ratio - wt_ratio
            fc = mut_ratio / wt_ratio if wt_ratio > 0 else (float("inf") if mut_ratio > 0 else 1.0)
            if fc >= fold_change_threshold and diff >= ratio_threshold:
                status = "hyper_in_mut"
            elif fc <= 1.0 / fold_change_threshold and diff <= -ratio_threshold:
                status = "hypo_in_mut"
            else:
                status = "shared"

        pseudo = 0.01
        results.append({
            "chrom": chrom,
            "pos": pos,
            "strand": strand,
            "gene_id": gene_id,
            "transcript_id": transcript_id,
            "region": region,
            "wt_ratio": wt_ratio,
            "mut_ratio": mut_ratio,
            "diff": mut_ratio - wt_ratio,
            "log2fc": np.log2((mut_ratio + pseudo) / (wt_ratio + pseudo)),
            "in_wt": in_wt,
            "in_mut": in_mut,
            "status": status,
        })

    df = pd.DataFrame(results, columns=DIFF_SITE_COLUMNS)
    if verbose:
        print(f"\n  Total unique positions: {len(df)}")
        if len(df):
            for status, count in df["status"].value_counts().items():
                print(f"    {status:15s}: {count:6d} ({count/len(df)*100:.1f}%)")
    return df


def gene_level_summary(df, verbose=True):
    """Summarize differential m6A classes at gene level."""
    gene_data = []
    annotated = df[df["gene_id"].notna()]
    for gene_id, sub in annotated.groupby("gene_id"):
        n_total = len(sub)
        n_wt_only = (sub["status"] == "wt_only").sum()
        n_mut_only = (sub["status"] == "mut_only").sum()
        n_shared = (sub["status"] == "shared").sum()
        n_hyper = (sub["status"] == "hyper_in_mut").sum()
        n_hypo = (sub["status"] == "hypo_in_mut").sum()
        n_in_wt = n_wt_only + n_shared + n_hyper + n_hypo
        n_in_mut = n_mut_only + n_shared + n_hyper + n_hypo
        gain_score = n_mut_only + n_hyper
        loss_score = n_wt_only + n_hypo
        n_differential = gain_score + loss_score

        if n_in_wt == 0 and n_in_mut > 0:
            gene_status = "gained"
        elif n_in_mut == 0 and n_in_wt > 0:
            gene_status = "lost"
        elif n_differential == 0:
            gene_status = "unchanged"
        elif gain_score >= 2 and loss_score >= 2 and min(gain_score, loss_score) / max(gain_score, loss_score) > 0.5:
            gene_status = "mixed"
        elif gain_score > loss_score and gain_score >= 2:
            gene_status = "hyper"
        elif loss_score > gain_score and loss_score >= 2:
            gene_status = "hypo"
        elif gain_score > loss_score:
            gene_status = "hyper"
        elif loss_score > gain_score:
            gene_status = "hypo"
        else:
            gene_status = "unchanged"

        gene_data.append({
            "gene_id": gene_id,
            "total_sites": n_total,
            "sites_in_wt": n_in_wt,
            "sites_in_mut": n_in_mut,
            "wt_only": n_wt_only,
            "mut_only": n_mut_only,
            "shared": n_shared,
            "hyper_in_mut": n_hyper,
            "hypo_in_mut": n_hypo,
            "gain_score": gain_score,
            "loss_score": loss_score,
            "net_change": gain_score - loss_score,
            "gene_status": gene_status,
            "mean_log2fc": sub["log2fc"].mean(),
        })

    gdf = pd.DataFrame(gene_data, columns=GENE_SUMMARY_COLUMNS)
    if verbose:
        print("\n  Gene-level summary:")
        for status in ["gained", "lost", "hyper", "hypo", "mixed", "unchanged"]:
            count = (gdf["gene_status"] == status).sum() if not gdf.empty else 0
            if count > 0:
                print(f"    {status:12s}: {count:6d} genes")
        print(f"    {'TOTAL':12s}: {len(gdf):6d} genes")
    return gdf


DIFF_SITE_COLUMNS = [
    "chrom",
    "pos",
    "strand",
    "gene_id",
    "transcript_id",
    "region",
    "wt_ratio",
    "mut_ratio",
    "diff",
    "log2fc",
    "in_wt",
    "in_mut",
    "status",
]

GENE_SUMMARY_COLUMNS = [
    "gene_id",
    "total_sites",
    "sites_in_wt",
    "sites_in_mut",
    "wt_only",
    "mut_only",
    "shared",
    "hyper_in_mut",
    "hypo_in_mut",
    "gain_score",
    "loss_score",
    "net_change",
    "gene_status",
    "mean_log2fc",
]


STATUS_COLORS = {
    "wt_only": "#2196f3",
    "mut_only": "#f44336",
    "shared": "#9e9e9e",
    "hyper_in_mut": "#c62828",
    "hypo_in_mut": "#1565c0",
}

STATUS_LABELS = {
    "wt_only": "WT only",
    "mut_only": "Mutant only",
    "shared": "Shared (unchanged)",
    "hyper_in_mut": "Hyper in mutant",
    "hypo_in_mut": "Hypo in mutant",
}


def _prepare_matplotlib():
    import logging

    import matplotlib
    import matplotlib.pyplot as plt

    logging.getLogger("fontTools").setLevel(logging.WARNING)
    matplotlib.rcParams["pdf.fonttype"] = 42
    matplotlib.rcParams["ps.fonttype"] = 42
    matplotlib.rcParams["font.family"] = "DejaVu Sans"
    matplotlib.rcParams["font.sans-serif"] = ["DejaVu Sans", "Liberation Sans"]
    return plt


def plot_volcano(df, wt_name, mut_name, fig_dir):
    plt = _prepare_matplotlib()
    fig, ax = plt.subplots(figsize=(7, 5.5))
    plot_df = df.copy()
    plot_df["mean_ratio"] = (plot_df["wt_ratio"] + plot_df["mut_ratio"]) / 2
    for status in ["shared", "wt_only", "mut_only", "hyper_in_mut", "hypo_in_mut"]:
        sub = plot_df[plot_df["status"] == status]
        if sub.empty:
            continue
        ax.scatter(sub["mean_ratio"], sub["log2fc"], c=STATUS_COLORS.get(status, "#999"), s=8, alpha=0.4, edgecolors="none", label=f"{STATUS_LABELS.get(status, status)} (n={len(sub)})")
    ax.axhline(0, color="grey", lw=0.8, ls=":")
    ax.axhline(np.log2(1.5), color="#c62828", lw=0.5, ls="--", alpha=0.5)
    ax.axhline(-np.log2(1.5), color="#1565c0", lw=0.5, ls="--", alpha=0.5)
    ax.set_xlabel("Mean m6A ratio", fontsize=10)
    ax.set_ylabel(f"log2FC ({mut_name} / {wt_name})", fontsize=10)
    ax.set_title(f"Differential m6A: {mut_name} vs {wt_name}", fontsize=11, fontweight="bold")
    ax.legend(fontsize=7, loc="upper left", markerscale=2)
    ax.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()
    for fmt in ["pdf", "svg", "png"]:
        fig.savefig(fig_dir / f"diff_m6a_volcano.{fmt}", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_region_distribution(df, wt_name, mut_name, fig_dir):
    plt = _prepare_matplotlib()
    regions = ["5utr", "cds", "3utr"]
    statuses = ["wt_only", "shared", "hyper_in_mut", "hypo_in_mut", "mut_only"]
    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.arange(len(regions))
    width = 0.15
    for i, status in enumerate(statuses):
        counts = [len(df[(df["region"] == region) & (df["status"] == status)]) for region in regions]
        ax.bar(x + i * width, counts, width, color=STATUS_COLORS.get(status, "#999"), label=STATUS_LABELS.get(status, status), edgecolor="white", linewidth=0.3)
    ax.set_xticks(x + width * 2)
    ax.set_xticklabels(["5' UTR", "CDS", "3' UTR"], fontsize=10)
    ax.set_ylabel("Number of m6A sites", fontsize=10)
    ax.set_title(f"Differential m6A by transcript region: {mut_name} vs {wt_name}", fontsize=11, fontweight="bold")
    ax.legend(fontsize=7)
    ax.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()
    for fmt in ["pdf", "svg", "png"]:
        fig.savefig(fig_dir / f"diff_m6a_region.{fmt}", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_ratio_comparison(df, wt_name, mut_name, fig_dir):
    plt = _prepare_matplotlib()
    shared = df[df["in_wt"] & df["in_mut"]].copy()
    if shared.empty:
        return
    fig, ax = plt.subplots(figsize=(6, 6))
    for status in ["shared", "hyper_in_mut", "hypo_in_mut"]:
        sub = shared[shared["status"] == status]
        if sub.empty:
            continue
        ax.scatter(sub["wt_ratio"], sub["mut_ratio"], c=STATUS_COLORS.get(status, "#999"), s=10, alpha=0.4, edgecolors="none", label=f"{STATUS_LABELS.get(status, status)} (n={len(sub)})")
    lim = max(ax.get_xlim()[1], ax.get_ylim()[1])
    ax.plot([0, lim], [0, lim], color="grey", lw=0.8, ls=":", alpha=0.5)
    ax.set_xlabel(f"{wt_name} m6A ratio", fontsize=10)
    ax.set_ylabel(f"{mut_name} m6A ratio", fontsize=10)
    ax.set_title("m6A ratio comparison (shared sites)", fontsize=11, fontweight="bold")
    ax.legend(fontsize=7)
    ax.set_aspect("equal")
    ax.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()
    for fmt in ["pdf", "svg", "png"]:
        fig.savefig(fig_dir / f"diff_m6a_scatter.{fmt}", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_summary_pie(df, wt_name, mut_name, fig_dir):
    plt = _prepare_matplotlib()
    status_order = ["wt_only", "hypo_in_mut", "shared", "hyper_in_mut", "mut_only"]
    counts = [len(df[df["status"] == status]) for status in status_order]
    labels = [f"{STATUS_LABELS[status]}\n({count})" for status, count in zip(status_order, counts)]
    colors = [STATUS_COLORS[status] for status in status_order]
    nonzero = [(label, count, color) for label, count, color in zip(labels, counts, colors) if count > 0]
    if not nonzero:
        return
    labels, counts, colors = zip(*nonzero)
    fig, ax = plt.subplots(figsize=(6, 6))
    _, _, autotexts = ax.pie(counts, labels=labels, colors=colors, autopct="%1.1f%%", startangle=90, pctdistance=0.75, textprops={"fontsize": 8})
    for text in autotexts:
        text.set_fontsize(7)
        text.set_color("white")
        text.set_fontweight("bold")
    ax.set_title(f"Differential m6A: {mut_name} vs {wt_name}", fontsize=11, fontweight="bold")
    plt.tight_layout()
    for fmt in ["pdf", "svg", "png"]:
        fig.savefig(fig_dir / f"diff_m6a_pie.{fmt}", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_ratio_distribution(df, wt_name, mut_name, fig_dir, verbose=True):
    plt = _prepare_matplotlib()
    fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
    wt_only = df[df["status"] == "wt_only"]["wt_ratio"]
    mut_only = df[df["status"] == "mut_only"]["mut_ratio"]
    bins = np.arange(0, 1.05, 0.05)

    axes[0].hist(wt_only, bins=bins, color=STATUS_COLORS["wt_only"], alpha=0.8, edgecolor="white", linewidth=0.3)
    axes[0].set_xlabel(f"{wt_name} m6A ratio", fontsize=10)
    axes[0].set_ylabel("Number of sites", fontsize=10)
    axes[0].set_title(f"{wt_name}-only sites (n={len(wt_only)})", fontsize=11, fontweight="bold", color=STATUS_COLORS["wt_only"])
    if len(wt_only):
        axes[0].axvline(wt_only.median(), color="black", ls="--", lw=1, label=f"Median={wt_only.median():.3f}")
    axes[0].legend(fontsize=8)
    axes[0].spines[["top", "right"]].set_visible(False)

    axes[1].hist(mut_only, bins=bins, color=STATUS_COLORS["mut_only"], alpha=0.8, edgecolor="white", linewidth=0.3)
    axes[1].set_xlabel(f"{mut_name} m6A ratio", fontsize=10)
    axes[1].set_title(f"{mut_name}-only sites (n={len(mut_only)})", fontsize=11, fontweight="bold", color=STATUS_COLORS["mut_only"])
    if len(mut_only):
        axes[1].axvline(mut_only.median(), color="black", ls="--", lw=1, label=f"Median={mut_only.median():.3f}")
    axes[1].legend(fontsize=8)
    axes[1].spines[["top", "right"]].set_visible(False)
    plt.suptitle("Ratio distribution of condition-specific m6A sites", fontsize=12, fontweight="bold")
    plt.tight_layout()
    for fmt in ["pdf", "svg", "png"]:
        fig.savefig(fig_dir / f"diff_m6a_ratio_hist.{fmt}", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 5.5))
    plot_data = []
    plot_labels = []
    plot_colors = []
    status_plot_order = [
        ("wt_only", "wt_ratio", f"{wt_name}-only"),
        ("hypo_in_mut", "wt_ratio", f"Hypo ({wt_name})"),
        ("hypo_in_mut", "mut_ratio", f"Hypo ({mut_name})"),
        ("shared", "wt_ratio", f"Shared ({wt_name})"),
        ("shared", "mut_ratio", f"Shared ({mut_name})"),
        ("hyper_in_mut", "wt_ratio", f"Hyper ({wt_name})"),
        ("hyper_in_mut", "mut_ratio", f"Hyper ({mut_name})"),
        ("mut_only", "mut_ratio", f"{mut_name}-only"),
    ]
    color_map = {
        f"{wt_name}-only": STATUS_COLORS["wt_only"],
        f"Hypo ({wt_name})": "#64b5f6",
        f"Hypo ({mut_name})": "#1565c0",
        f"Shared ({wt_name})": "#bdbdbd",
        f"Shared ({mut_name})": "#757575",
        f"Hyper ({wt_name})": "#ef9a9a",
        f"Hyper ({mut_name})": "#c62828",
        f"{mut_name}-only": STATUS_COLORS["mut_only"],
    }
    for status, ratio_col, label in status_plot_order:
        vals = df[df["status"] == status][ratio_col].values
        if len(vals) < 5:
            continue
        plot_data.append(vals)
        plot_labels.append(f"{label}\n(n={len(vals)})")
        plot_colors.append(color_map.get(label, "#999"))
    if plot_data:
        violin = ax.violinplot(plot_data, positions=range(len(plot_data)), showmedians=False, showextrema=False)
        for i, body in enumerate(violin["bodies"]):
            body.set_facecolor(plot_colors[i])
            body.set_alpha(0.5)
        box = ax.boxplot(plot_data, positions=range(len(plot_data)), widths=0.2, showfliers=False, patch_artist=True, zorder=3)
        for i, patch in enumerate(box["boxes"]):
            patch.set_facecolor(plot_colors[i])
            patch.set_alpha(0.8)
        for element in ["whiskers", "caps", "medians"]:
            for line in box[element]:
                line.set_color("black")
                line.set_linewidth(0.8)
        ax.set_xticks(range(len(plot_data)))
        ax.set_xticklabels(plot_labels, fontsize=7, rotation=30, ha="right")
        ax.set_ylabel("m6A ratio", fontsize=10)
        ax.set_title("m6A ratio distribution by site category", fontsize=11, fontweight="bold")
        ax.spines[["top", "right"]].set_visible(False)
        plt.tight_layout()
        for fmt in ["pdf", "svg", "png"]:
            fig.savefig(fig_dir / f"diff_m6a_ratio_violin.{fmt}", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 5))
    for status, ratio_col, label, color in [
        ("wt_only", "wt_ratio", f"{wt_name}-only", STATUS_COLORS["wt_only"]),
        ("mut_only", "mut_ratio", f"{mut_name}-only", STATUS_COLORS["mut_only"]),
        ("shared", "wt_ratio", f"Shared ({wt_name})", "#9e9e9e"),
        ("shared", "mut_ratio", f"Shared ({mut_name})", "#616161"),
    ]:
        vals = df[df["status"] == status][ratio_col].values
        if len(vals) < 10:
            continue
        sorted_vals = np.sort(vals)
        cdf = np.arange(1, len(sorted_vals) + 1) / len(sorted_vals)
        ax.plot(sorted_vals, cdf, color=color, lw=1.5, alpha=0.8, label=f"{label} (n={len(vals)}, med={np.median(vals):.3f})")
    ax.set_xlabel("m6A ratio", fontsize=10)
    ax.set_ylabel("Cumulative fraction", fontsize=10)
    ax.set_title("Cumulative distribution of m6A ratios", fontsize=11, fontweight="bold")
    handles, labels = ax.get_legend_handles_labels()
    if handles:
        ax.legend(handles, labels, fontsize=7, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
    ax.set_xlim(0, 1)
    plt.tight_layout()
    for fmt in ["pdf", "svg", "png"]:
        fig.savefig(fig_dir / f"diff_m6a_ratio_cdf.{fmt}", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)

    if verbose:
        print("\n  Ratio distribution summary:")
        print(f"    {'Category':20s} {'n':>8} {'median':>8} {'mean':>8} {'<0.2':>8} {'≥0.5':>8}")
        print(f"    {'-'*56}")
        for status, ratio_col, label in [
            ("wt_only", "wt_ratio", f"{wt_name}-only"),
            ("mut_only", "mut_ratio", f"{mut_name}-only"),
            ("shared", "wt_ratio", f"Shared ({wt_name})"),
            ("shared", "mut_ratio", f"Shared ({mut_name})"),
            ("hyper_in_mut", "mut_ratio", f"Hyper ({mut_name})"),
            ("hypo_in_mut", "wt_ratio", f"Hypo ({wt_name})"),
        ]:
            vals = df[df["status"] == status][ratio_col]
            if len(vals) == 0:
                continue
            print(f"    {label:20s} {len(vals):8d} {vals.median():8.3f} {vals.mean():8.3f} {(vals < 0.2).sum():8d} {(vals >= 0.5).sum():8d}")


def plot_differential_m6a(df, wt_name, mut_name, fig_dir, verbose=True):
    fig_dir = Path(fig_dir)
    fig_dir.mkdir(parents=True, exist_ok=True)
    plot_volcano(df, wt_name, mut_name, fig_dir)
    plot_region_distribution(df, wt_name, mut_name, fig_dir)
    plot_ratio_comparison(df, wt_name, mut_name, fig_dir)
    plot_summary_pie(df, wt_name, mut_name, fig_dir)
    plot_ratio_distribution(df, wt_name, mut_name, fig_dir, verbose=verbose)


def write_go_inputs(df, gdf, output_dir, verbose=True):
    go_dir = Path(output_dir) / "go_input"
    go_dir.mkdir(parents=True, exist_ok=True)
    outputs = {}
    for status in ["hyper", "hypo", "gained", "lost", "mixed"]:
        genes = gdf[gdf["gene_status"] == status]["gene_id"].tolist() if not gdf.empty else []
        path = go_dir / f"genes_{status}.txt"
        path.write_text("".join(f"{gene}\n" for gene in sorted(set(genes))))
        outputs[status] = str(path)
        if verbose:
            print(f"  GO input ({status}): {len(genes)} genes")

    hyper_all = gdf[gdf["gene_status"].isin(["hyper", "gained"])]["gene_id"].tolist() if not gdf.empty else []
    path = go_dir / "genes_hyper_all.txt"
    path.write_text("".join(f"{gene}\n" for gene in sorted(set(hyper_all))))
    outputs["hyper_all"] = str(path)
    if verbose:
        print(f"  GO input (hyper+gained): {len(set(hyper_all))} genes")

    hypo_all = gdf[gdf["gene_status"].isin(["hypo", "lost"])]["gene_id"].tolist() if not gdf.empty else []
    path = go_dir / "genes_hypo_all.txt"
    path.write_text("".join(f"{gene}\n" for gene in sorted(set(hypo_all))))
    outputs["hypo_all"] = str(path)
    if verbose:
        print(f"  GO input (hypo+lost): {len(set(hypo_all))} genes")

    bg_genes = set(df[df["gene_id"].notna()]["gene_id"].tolist())
    path = go_dir / "background_all_m6a.txt"
    path.write_text("".join(f"{gene}\n" for gene in sorted(bg_genes)))
    outputs["background"] = str(path)
    if verbose:
        print(f"  GO background: {len(bg_genes)} genes")
    return outputs


def run_differential_m6a(
    wt_bed,
    mut_bed,
    annotation,
    output_dir="results/diff_m6a",
    wt_name="WT",
    mut_name="Mutant",
    ratio_threshold=0.1,
    fold_change_threshold=1.5,
    make_plots=True,
    verbose=True,
):
    """Run the complete differential m6A workflow and write legacy-compatible outputs."""
    output_dir = Path(output_dir)
    fig_dir = output_dir / "figures"
    output_dir.mkdir(parents=True, exist_ok=True)
    if make_plots:
        fig_dir.mkdir(parents=True, exist_ok=True)

    if verbose:
        print("=" * 60)
        print(f"Differential m6A: {mut_name} vs {wt_name}")
        print("=" * 60)

    models, binned_index = parse_differential_annotation(annotation, verbose=verbose)
    wt_sites = annotate_differential_sites(load_bed6_sites(wt_bed, verbose=verbose), models, binned_index, verbose=verbose)
    mut_sites = annotate_differential_sites(load_bed6_sites(mut_bed, verbose=verbose), models, binned_index, verbose=verbose)
    df = find_differential_m6a(wt_sites, mut_sites, ratio_threshold=ratio_threshold, fold_change_threshold=fold_change_threshold, verbose=verbose)

    out_all = output_dir / "diff_m6a_all.tsv"
    df.to_csv(out_all, sep="\t", index=False, float_format="%.4f")
    if verbose:
        print(f"\n  All sites: {out_all}")

    region_counts = df.groupby(["region", "status"]).size().unstack(fill_value=0) if not df.empty else pd.DataFrame()
    out_region = output_dir / "diff_m6a_by_region.tsv"
    region_counts.to_csv(out_region, sep="\t")
    if verbose:
        print(f"  By region: {out_region}")

    gdf = gene_level_summary(df, verbose=verbose)
    out_gene = output_dir / "diff_m6a_by_gene.tsv"
    gdf.to_csv(out_gene, sep="\t", index=False, float_format="%.4f")
    if verbose:
        print(f"  By gene: {out_gene}")

    summary = {
        "wt_name": wt_name,
        "mut_name": mut_name,
        "total_wt_sites": len(wt_sites),
        "total_mut_sites": len(mut_sites),
        "total_unique_positions": len(df),
    }
    for status in ["wt_only", "mut_only", "shared", "hyper_in_mut", "hypo_in_mut"]:
        summary[status] = int((df["status"] == status).sum()) if not df.empty else 0
    summary["total_differential"] = summary["wt_only"] + summary["mut_only"] + summary["hyper_in_mut"] + summary["hypo_in_mut"]
    summary["genes_with_diff_m6a"] = int((gdf["gene_status"] != "unchanged").sum()) if not gdf.empty else 0
    out_summary = output_dir / "diff_m6a_summary.tsv"
    pd.DataFrame([summary]).to_csv(out_summary, sep="\t", index=False)
    if verbose:
        print(f"  Summary: {out_summary}")

    go_outputs = write_go_inputs(df, gdf, output_dir, verbose=verbose)

    plot_outputs = []
    if make_plots:
        if verbose:
            print("\nGenerating figures...")
        plot_differential_m6a(df, wt_name, mut_name, fig_dir, verbose=verbose)
        plot_outputs = [str(path) for path in sorted(fig_dir.glob("diff_m6a_*.*"))]

    if verbose:
        print(f"\n{'='*60}")
        print("Done!")
        print(f"  Results: {output_dir}/")
        if make_plots:
            print(f"  Figures: {fig_dir}/")
        print(f"  GO input: {output_dir / 'go_input'}/")
        print(f"{'='*60}")

    return {
        "all_sites": str(out_all),
        "by_region": str(out_region),
        "by_gene": str(out_gene),
        "summary": str(out_summary),
        "go_inputs": go_outputs,
        "figures": plot_outputs,
        "site_counts": summary,
    }


__all__ = [
    "DifferentialGeneModel",
    "annotate_differential_site",
    "annotate_differential_sites",
    "find_differential_m6a",
    "gene_level_summary",
    "load_bed6_sites",
    "parse_differential_annotation",
    "plot_differential_m6a",
    "run_differential_m6a",
    "write_go_inputs",
]
