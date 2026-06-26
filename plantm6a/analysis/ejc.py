#!/usr/bin/env python3
"""Exon-junction-exon triplet m6A stacking analysis."""

import bisect
import math
import sys
from collections import defaultdict
from pathlib import Path


DEFAULT_EJC_SPECIES_MAP = {
    "human": ("Homo_sapiens.GRCh38.115.bed12", "human_ratioover0.1.bed6"),
    "arabidopsis": ("TAIR10.release55.bed12", "Col_seedling7days_mR_coverageover20_ratioover0.1.bed6"),
    "rice": ("IRGSP-1.0.release57.bed12", "Nip_7daysseddlings_cam_mR_coverageover20_ratioover0.1.bed6"),
    "maize": ("Zm-B73-REFERENCE-NAM-5.0.release55.bed12", "B73_cam_mR_variety_coverageover20_ratioover0.1.bed6"),
    "soybean": ("Glycine_max_v2.1.release55.bed12", "W82_cam_mR_variety_coverageover20_ratioover0.1.bed6"),
    "tomato": ("SL3.0.release57.bed12", "tomato_cam_mR_variety_coverageover20_ratioover0.1_fixed.bed6"),
    "tobacco": ("NbT2T.annot.bed12", "Tobacco_cam_mR_variety_coverageover20_ratioover0.1.bed6"),
    "potato": ("A157.bed12", "Cam_Potato_32L_A157_leaf_coverageover20_ratioover0.1.bed6"),
    "cucumber": ("ChineseLong_v3.bed12", "cucumber_mR_coverageover20_ratioover0.1_fixed.bed6"),
    "marchantia": ("MpTak2_v7.1.fa.pseudo_label.bed12", "Marchantia_mR_coverageover20_ratioover0.1.bed6"),
    "chlamydomonas": ("GCF_000002595.2_Chlamydomonas_reinhardtii_v5.5_genomic.bed12", "AlageWT_mR_coverageover20_ratioover0.1_fixed.bed6"),
    "fly": ("BDGP6.32.release110.bed12", "Fly_merge_all_sort_coverageover20_ratioover0.1.bed6"),
}


class ExonJunctionTripletAnalyzer:
    """Analyze m6A enrichment around exon-junction-exon triplets."""

    def __init__(self, window_size=10, upstream_range=500, downstream_range=500, max_range=2000, smooth_zero=False):
        self.window_size = window_size
        self.upstream_range = upstream_range
        self.downstream_range = downstream_range
        self.max_range = max_range
        self.smooth_zero = smooth_zero
        self.internal_m6a = defaultdict(int)
        self.internal_coverage = defaultdict(int)
        self.last_m6a = defaultdict(int)
        self.last_coverage = defaultdict(int)
        self.stats = {
            "total_genes": 0,
            "internal_junctions": 0,
            "last_junctions": 0,
            "skipped_single_exon": 0,
        }

    def parse_bed12_line(self, line):
        fields = line.strip().split("\t")
        if len(fields) < 12:
            return None

        chrom = fields[0]
        tx_start = int(fields[1])
        strand = fields[5]
        exon_count = int(fields[9])
        exon_sizes = [int(x) for x in fields[10].rstrip(",").split(",") if x]
        exon_starts = [int(x) for x in fields[11].rstrip(",").split(",") if x]
        gene_id = fields[3]

        exons = []
        for i in range(exon_count):
            exon_start = tx_start + exon_starts[i]
            exon_end = exon_start + exon_sizes[i]
            exons.append((exon_start, exon_end))

        return {
            "chrom": chrom,
            "strand": strand,
            "exons": exons,
            "gene_id": gene_id,
            "exon_count": exon_count,
        }

    def extract_junction_triplets(self, gene_info):
        exons = gene_info["exons"]
        strand = gene_info["strand"]
        if len(exons) < 2:
            return []

        triplets = []
        for i in range(len(exons) - 1):
            upstream_exon = exons[i]
            downstream_exon = exons[i + 1]
            junction_pos = [upstream_exon[1], downstream_exon[0]]
            if strand == "+":
                is_last = i == len(exons) - 2
            else:
                is_last = i == 0
            triplets.append({
                "upstream_exon": upstream_exon,
                "downstream_exon": downstream_exon,
                "junction_pos": junction_pos,
                "is_last": is_last,
                "strand": strand,
                "junction_index": i,
            })
        return triplets

    def map_position_to_bin(self, genomic_pos, junction_pos, strand):
        if genomic_pos >= junction_pos[1]:
            relative_pos = genomic_pos - junction_pos[1]
        elif genomic_pos <= junction_pos[0]:
            relative_pos = genomic_pos - junction_pos[0]
        else:
            relative_pos = 0
        if strand == "-":
            relative_pos = -relative_pos
        return math.floor(relative_pos / self.window_size)

    def process_triplet(self, triplet, m6a_sites, is_internal):
        upstream_exon = triplet["upstream_exon"]
        downstream_exon = triplet["downstream_exon"]
        junction_pos = triplet["junction_pos"]
        strand = triplet["strand"]

        if is_internal:
            m6a_dict = self.internal_m6a
            coverage_dict = self.internal_coverage
        else:
            m6a_dict = self.last_m6a
            coverage_dict = self.last_coverage

        for m6a_pos in m6a_sites:
            if (
                abs(upstream_exon[0] - m6a_pos) + abs(m6a_pos - upstream_exon[1]) == abs(upstream_exon[0] - upstream_exon[1])
                or abs(downstream_exon[0] - m6a_pos) + abs(m6a_pos - downstream_exon[1]) == abs(downstream_exon[0] - downstream_exon[1])
            ):
                bin_idx = self.map_position_to_bin(m6a_pos, junction_pos, strand)
                if bin_idx is not None:
                    m6a_dict[bin_idx] += 1

        for pos in range(upstream_exon[0], upstream_exon[1]):
            bin_idx = self.map_position_to_bin(pos, junction_pos, strand)
            if bin_idx is not None:
                coverage_dict[bin_idx] += 1 / self.window_size

        for pos in range(downstream_exon[0], downstream_exon[1]):
            bin_idx = self.map_position_to_bin(pos, junction_pos, strand)
            if bin_idx is not None:
                coverage_dict[bin_idx] += 1 / self.window_size

    def load_m6a_sites(self, bed6_file):
        m6a_dict = defaultdict(list)
        print(f"Loading m6A sites from {bed6_file}...", file=sys.stderr)
        with open(bed6_file) as f:
            for line in f:
                fields = line.strip().split("\t")
                if len(fields) < 6:
                    continue
                chrom = fields[0]
                start = int(fields[1])
                m6a_dict[chrom].append(start)
        for chrom in m6a_dict:
            m6a_dict[chrom].sort()
        total_sites = sum(len(sites) for sites in m6a_dict.values())
        print(f"Loaded {total_sites} m6A sites from {len(m6a_dict)} chromosomes", file=sys.stderr)
        return m6a_dict

    def analyze(self, bed12_file, bed6_file):
        print(f"Analyzing junctions from {bed12_file}...", file=sys.stderr)
        print(f"Window size: {self.window_size} bp", file=sys.stderr)
        print(f"Smooth zero: {self.smooth_zero}", file=sys.stderr)
        m6a_dict = self.load_m6a_sites(bed6_file)

        with open(bed12_file) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                gene_info = self.parse_bed12_line(line)
                if not gene_info:
                    continue
                self.stats["total_genes"] += 1
                if gene_info["exon_count"] < 2:
                    self.stats["skipped_single_exon"] += 1
                    continue
                if self.stats["total_genes"] % 5000 == 0:
                    print(f"  Processed {self.stats['total_genes']} genes...", file=sys.stderr)

                triplets = self.extract_junction_triplets(gene_info)
                chrom = gene_info["chrom"]
                m6a_sites = m6a_dict.get(chrom, [])
                for triplet in triplets:
                    if triplet["is_last"]:
                        self.process_triplet(triplet, m6a_sites, is_internal=False)
                        self.stats["last_junctions"] += 1
                    else:
                        self.process_triplet(triplet, m6a_sites, is_internal=True)
                        self.stats["internal_junctions"] += 1

        print("\n" + "=" * 60, file=sys.stderr)
        print("Analysis Statistics:", file=sys.stderr)
        print("=" * 60, file=sys.stderr)
        print(f"Total genes processed: {self.stats['total_genes']}", file=sys.stderr)
        print(f"Single-exon genes skipped: {self.stats['skipped_single_exon']}", file=sys.stderr)
        print(f"Internal junctions: {self.stats['internal_junctions']}", file=sys.stderr)
        print(f"Last junctions: {self.stats['last_junctions']}", file=sys.stderr)
        print("=" * 60 + "\n", file=sys.stderr)

    def smooth_bin_zero(self, m6a_dict, coverage_dict, label=""):
        if -1 in coverage_dict and 1 in coverage_dict:
            m6a_m1 = m6a_dict.get(-1, 0)
            m6a_p1 = m6a_dict.get(1, 0)
            cov_m1 = coverage_dict.get(-1, 0)
            cov_p1 = coverage_dict.get(1, 0)
            orig_m6a = m6a_dict.get(0, 0)
            orig_cov = coverage_dict.get(0, 0)
            m6a_dict[0] = round((m6a_m1 + m6a_p1) / 2)
            coverage_dict[0] = (cov_m1 + cov_p1) / 2
            orig_like = orig_m6a / orig_cov if orig_cov > 0 else 0
            new_like = m6a_dict[0] / coverage_dict[0] if coverage_dict[0] > 0 else 0
            print(f"  Smoothed bin 0 ({label}):", file=sys.stderr)
            print(f"    m6a_level:      {orig_m6a:>10} → {m6a_dict[0]:>10}", file=sys.stderr)
            print(f"    coverage:       {orig_cov:>10.1f} → {coverage_dict[0]:>10.1f}", file=sys.stderr)
            print(f"    likelihood:     {orig_like:>10.4f} → {new_like:>10.4f}", file=sys.stderr)

    def calculate_metrics(self, m6a_dict, coverage_dict, label=""):
        if self.smooth_zero:
            self.smooth_bin_zero(m6a_dict, coverage_dict, label)
        results = []
        all_bins = sorted(set(list(m6a_dict.keys()) + list(coverage_dict.keys())))
        for bin_idx in all_bins:
            m6a_level = m6a_dict.get(bin_idx, 0)
            coverage = coverage_dict.get(bin_idx, 0)
            likelihood = m6a_level / coverage if coverage > 0 else 0
            density = m6a_level / self.window_size if self.window_size > 0 else 0
            distance = bin_idx * self.window_size
            results.append({
                "bin": bin_idx,
                "distance": distance,
                "m6a_level": m6a_level,
                "coverage": coverage,
                "m6a_density": density,
                "m6a_likelihood": likelihood,
            })
        return results

    def write_output(self, output_prefix):
        output_path = Path(output_prefix)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        internal_results = self.calculate_metrics(self.internal_m6a, self.internal_coverage, "internal")
        internal_file = f"{output_prefix}_internal_junctions.txt"
        self._write_metrics(internal_file, internal_results)
        last_results = self.calculate_metrics(self.last_m6a, self.last_coverage, "last")
        last_file = f"{output_prefix}_last_junctions.txt"
        self._write_metrics(last_file, last_results)
        print("\n=== Bin symmetry check ===", file=sys.stderr)
        for label, results in [("Internal", internal_results), ("Last", last_results)]:
            bin_0 = next((r for r in results if r["bin"] == 0), None)
            bin_m1 = next((r for r in results if r["bin"] == -1), None)
            bin_1 = next((r for r in results if r["bin"] == 1), None)
            if bin_0 and bin_m1 and bin_1:
                print(
                    f"  {label}: cov[-1]={bin_m1['coverage']:.1f}  cov[0]={bin_0['coverage']:.1f}  cov[+1]={bin_1['coverage']:.1f}",
                    file=sys.stderr,
                )
                print(
                    f"  {label}: like[-1]={bin_m1['m6a_likelihood']:.4f}  like[0]={bin_0['m6a_likelihood']:.4f}  like[+1]={bin_1['m6a_likelihood']:.4f}",
                    file=sys.stderr,
                )
        print("\nResults written to:", file=sys.stderr)
        print(f"  - {internal_file} ({len(internal_results)} bins)", file=sys.stderr)
        print(f"  - {last_file} ({len(last_results)} bins)", file=sys.stderr)
        return internal_file, last_file

    @staticmethod
    def _write_metrics(output_file, results):
        with open(output_file, "w") as f:
            f.write("bin\tdistance_from_junction\tm6a_level\tcoverage\tm6a_density\tm6a_likelihood\n")
            for item in results:
                f.write(
                    f"{item['bin']}\t{item['distance']}\t{item['m6a_level']}\t"
                    f"{item['coverage']}\t{item['m6a_density']:.6f}\t{item['m6a_likelihood']:.6f}\n"
                )


def run_ejc_triplet(annotation_bed12, m6a_bed6, output_prefix, window_size=10, upstream_range=500, downstream_range=500, smooth_zero=False):
    analyzer = ExonJunctionTripletAnalyzer(
        window_size=window_size,
        upstream_range=upstream_range,
        downstream_range=downstream_range,
        smooth_zero=smooth_zero,
    )
    analyzer.analyze(annotation_bed12, m6a_bed6)
    internal_file, last_file = analyzer.write_output(output_prefix)
    return {
        "internal_file": internal_file,
        "last_file": last_file,
        "stats": analyzer.stats,
    }


def parse_species_map_file(path):
    species_map = {}
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                raise ValueError("Species map rows must have: species<TAB>bed12_file<TAB>bed6_file")
            species_map[parts[0]] = (parts[1], parts[2])
    return species_map


def _prepare_ejc_plot_data(path, x_limit=1000):
    pd = __import__("pandas")
    np = __import__("numpy")
    data = pd.read_csv(path, sep="\t")
    data = data[
        (data["distance_from_junction"] >= -x_limit)
        & (data["distance_from_junction"] <= x_limit)
    ].copy()
    data["se"] = np.sqrt(
        data["m6a_likelihood"] * (1 - data["m6a_likelihood"]) / data["coverage"].clip(lower=1)
    )
    data["ci_lower"] = (data["m6a_likelihood"] - 1.96 * data["se"]).clip(lower=0)
    data["ci_upper"] = data["m6a_likelihood"] + 1.96 * data["se"]
    return data


def _style_ejc_axis(ax, title, ylabel, x_limit):
    ax.axvline(x=0, linestyle="--", color="black", linewidth=0.5)
    ax.set_xlim(-x_limit, x_limit)
    ax.set_xticks(list(range(-x_limit, x_limit + 1, 500)))
    ax.set_title(title, fontsize=11, fontstyle="italic")
    ax.set_xlabel("Distance to exon junction", fontsize=10)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.grid(True, which="major", color="0.90", linewidth=0.3)
    ax.grid(False, which="minor")
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.5)
        spine.set_color("black")


def _plot_level_panel(ax, data, species_name, junction_label, color, x_limit):
    ax.plot(data["distance_from_junction"], data["m6a_level"], color=color, linewidth=0.8)
    _style_ejc_axis(ax, f"{species_name} - {junction_label} Junctions", "m6A level", x_limit)


def _plot_likelihood_panel(ax, data, species_name, junction_label, color, fill, x_limit):
    x = data["distance_from_junction"].to_numpy(dtype=float)
    lower = data["ci_lower"].to_numpy(dtype=float)
    upper = data["ci_upper"].to_numpy(dtype=float)
    y = data["m6a_likelihood"].to_numpy(dtype=float)
    ax.fill_between(x, lower, upper, color=fill, alpha=0.3, linewidth=0)
    ax.plot(x, y, color=color, linewidth=0.8)
    _style_ejc_axis(ax, f"{species_name} - {junction_label} Junctions", "m6A likelihood", x_limit)


def plot_ejc_triplet(internal_file, last_file, species_name, output_dir=".", x_limit=1000):
    """Draw EJC triplet SVG figures matching the legacy visualize_triplet.R outputs."""
    import matplotlib.pyplot as plt
    import pandas as pd

    internal_data = _prepare_ejc_plot_data(internal_file, x_limit=x_limit)
    last_data = _prepare_ejc_plot_data(last_file, x_limit=x_limit)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    outputs = []

    plt.style.use("seaborn-v0_8-whitegrid")

    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    _plot_level_panel(axes[0], internal_data, species_name, "Internal", "#2C3E50", x_limit)
    _plot_level_panel(axes[1], last_data, species_name, "Last", "#E74C3C", x_limit)
    fig.suptitle(f"{species_name} - m6A Level at Junctions", fontsize=14, fontweight="bold")
    fig.tight_layout()
    path = output_dir / f"{species_name}_junction_m6a_level.svg"
    fig.savefig(path, format="svg", bbox_inches="tight")
    outputs.append(str(path))
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    _plot_likelihood_panel(axes[0], internal_data, species_name, "Internal", "#2C3E50", "#3498DB", x_limit)
    _plot_likelihood_panel(axes[1], last_data, species_name, "Last", "#E74C3C", "#E67E22", x_limit)
    fig.suptitle(f"{species_name} - m6A Likelihood at Junctions", fontsize=14, fontweight="bold")
    fig.tight_layout()
    path = output_dir / f"{species_name}_junction_m6a_likelihood.svg"
    fig.savefig(path, format="svg", bbox_inches="tight")
    outputs.append(str(path))
    plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    _plot_level_panel(axes[0, 0], internal_data, species_name, "Internal", "#2C3E50", x_limit)
    _plot_level_panel(axes[0, 1], last_data, species_name, "Last", "#E74C3C", x_limit)
    _plot_likelihood_panel(axes[1, 0], internal_data, species_name, "Internal", "#2C3E50", "#3498DB", x_limit)
    _plot_likelihood_panel(axes[1, 1], last_data, species_name, "Last", "#E74C3C", "#E67E22", x_limit)
    fig.suptitle(f"{species_name} - Junction m6A Analysis", fontsize=16, fontweight="bold")
    fig.tight_layout()
    path = output_dir / f"{species_name}_junction_combined.svg"
    fig.savefig(path, format="svg", bbox_inches="tight")
    outputs.append(str(path))
    plt.close(fig)

    combined_level = pd.concat([
        internal_data[["distance_from_junction", "m6a_level"]].assign(type="Internal"),
        last_data[["distance_from_junction", "m6a_level"]].assign(type="Last"),
    ])
    fig, ax = plt.subplots(figsize=(8, 6))
    for label, color in [("Internal", "#3498DB"), ("Last", "#E74C3C")]:
        subset = combined_level[combined_level["type"] == label]
        ax.plot(subset["distance_from_junction"], subset["m6a_level"], color=color, linewidth=1, label=label)
    _style_ejc_axis(ax, f"{species_name} - m6A Level Comparison", "m6A level", x_limit)
    ax.legend(title="Junction Type", loc="upper center", bbox_to_anchor=(0.5, 1.12), ncol=2)
    fig.tight_layout()
    path = output_dir / f"{species_name}_level_comparison.svg"
    fig.savefig(path, format="svg", bbox_inches="tight")
    outputs.append(str(path))
    plt.close(fig)

    combined_likelihood = pd.concat([
        internal_data[["distance_from_junction", "m6a_likelihood", "ci_lower", "ci_upper"]].assign(type="Internal"),
        last_data[["distance_from_junction", "m6a_likelihood", "ci_lower", "ci_upper"]].assign(type="Last"),
    ])
    fig, ax = plt.subplots(figsize=(8, 6))
    for label, color in [("Internal", "#3498DB"), ("Last", "#E74C3C")]:
        subset = combined_likelihood[combined_likelihood["type"] == label]
        x = subset["distance_from_junction"].to_numpy(dtype=float)
        lower = subset["ci_lower"].to_numpy(dtype=float)
        upper = subset["ci_upper"].to_numpy(dtype=float)
        y = subset["m6a_likelihood"].to_numpy(dtype=float)
        ax.fill_between(x, lower, upper, color=color, alpha=0.2, linewidth=0)
        ax.plot(x, y, color=color, linewidth=1, label=label)
    _style_ejc_axis(ax, f"{species_name} - m6A Likelihood Comparison", "m6A likelihood", x_limit)
    ax.legend(title="Junction Type", loc="upper center", bbox_to_anchor=(0.5, 1.12), ncol=2)
    fig.tight_layout()
    path = output_dir / f"{species_name}_likelihood_comparison.svg"
    fig.savefig(path, format="svg", bbox_inches="tight")
    outputs.append(str(path))
    plt.close(fig)

    return outputs


def run_ejc_batch(
    bed12_dir,
    bed6_dir,
    output_dir,
    species_map=None,
    species=None,
    window_size=50,
    upstream_range=500,
    downstream_range=500,
    smooth_zero=True,
    strict=False,
):
    species_map = species_map or DEFAULT_EJC_SPECIES_MAP
    selected = set(species) if species else set(species_map)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    results = {}

    for sp in sorted(species_map):
        if sp not in selected:
            continue
        bed12_file, bed6_file = species_map[sp]
        bed12_path = Path(bed12_dir) / bed12_file
        bed6_path = Path(bed6_dir) / bed6_file
        output_prefix = output_dir / sp
        print(f"Processing: {sp}", file=sys.stderr)
        if not bed12_path.is_file():
            message = f"Annotation not found: {bed12_path}"
            if strict:
                raise FileNotFoundError(message)
            print(f"  ERROR: {message}", file=sys.stderr)
            results[sp] = {"status": "missing_annotation", "annotation": str(bed12_path), "m6a": str(bed6_path)}
            continue
        if not bed6_path.is_file():
            message = f"m6A file not found: {bed6_path}"
            if strict:
                raise FileNotFoundError(message)
            print(f"  ERROR: {message}", file=sys.stderr)
            results[sp] = {"status": "missing_m6a", "annotation": str(bed12_path), "m6a": str(bed6_path)}
            continue
        result = run_ejc_triplet(
            str(bed12_path),
            str(bed6_path),
            str(output_prefix),
            window_size=window_size,
            upstream_range=upstream_range,
            downstream_range=downstream_range,
            smooth_zero=smooth_zero,
        )
        result.update({"status": "completed", "annotation": str(bed12_path), "m6a": str(bed6_path)})
        results[sp] = result
        print(f"  ✓ Completed: {sp}", file=sys.stderr)
    return results


__all__ = [
    "DEFAULT_EJC_SPECIES_MAP",
    "ExonJunctionTripletAnalyzer",
    "parse_species_map_file",
    "plot_ejc_triplet",
    "run_ejc_batch",
    "run_ejc_triplet",
]
