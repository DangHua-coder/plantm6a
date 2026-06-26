#!/usr/bin/env python3
"""
Command line interface for plantm6a.
"""

import argparse
import sys
from pathlib import Path



def _require_import(import_callback, package_hint):
    try:
        return import_callback()
    except ImportError as exc:
        raise SystemExit(
            f"Missing optional dependency for this command: {exc}. "
            f"Install with: pip install {package_hint}"
        ) from exc


def _write_key_value_tsv(result, output_file):
    if not output_file:
        return
    with open(output_file, "w") as out:
        out.write("metric\tvalue\n")
        for key, value in result.items():
            if not isinstance(value, (list, dict, tuple)):
                out.write(f"{key}\t{value}\n")


def run_stats(args):
    from plantm6a.analysis.statistics import analyze

    result = analyze(args.genome, args.gtf, args.bed, verbose=True)
    _write_key_value_tsv(result, args.output)


def run_batch(args):
    batch_analyze = _require_import(
        lambda: __import__(
            "plantm6a.analysis.batch", fromlist=["batch_analyze"]
        ).batch_analyze,
        "plantm6a",
    )
    batch_analyze(args.config, args.output)


def run_motif(args):
    analyze_motifs = _require_import(
        lambda: __import__(
            "plantm6a.analysis.motif", fromlist=["analyze_motifs"]
        ).analyze_motifs,
        "plantm6a[motif]",
    )
    result = analyze_motifs(
        genome_path=args.genome,
        annotation_path=args.annotation,
        bed_path=args.bed,
        mode=args.mode,
        verbose=True,
    )
    _write_key_value_tsv(result, args.output)


def run_annotate_sites(args):
    annotate_m6a_sites = _require_import(
        lambda: __import__(
            "plantm6a.analysis.site_annotation", fromlist=["annotate_m6a_sites"]
        ).annotate_m6a_sites,
        "plantm6a[motif]",
    )
    annotate_m6a_sites(
        genome_path=args.genome,
        annotation_path=args.annotation,
        bed_path=args.bed,
        output_path=args.output,
        annotation_format=args.fmt,
        motif_window=args.motif_window,
        overlap_policy=args.overlap_policy,
        verbose=True,
    )


def run_translate_cds(args):
    translate_cds = _require_import(
        lambda: __import__(
            "plantm6a.analysis.conservation", fromlist=["translate_cds"]
        ).translate_cds,
        "plantm6a[conservation]",
    )
    translate_cds(args.input, args.output, verbose=True)


def run_extract_orthologs(args):
    extract_pairwise_orthologs = _require_import(
        lambda: __import__(
            "plantm6a.analysis.conservation", fromlist=["extract_pairwise_orthologs"]
        ).extract_pairwise_orthologs,
        "plantm6a[conservation]",
    )
    extract_pairwise_orthologs(
        args.orthogroups,
        args.orthologues_dir,
        args.output_dir,
        verbose=True,
    )


def run_conserve_pair(args):
    analyze_species_pair = _require_import(
        lambda: __import__(
            "plantm6a.analysis.conservation", fromlist=["analyze_species_pair"]
        ).analyze_species_pair,
        "plantm6a[conservation]",
    )
    analyze_species_pair(
        spA=args.sp_a,
        spB=args.sp_b,
        genome_A=args.genome_a,
        genome_B=args.genome_b,
        annot_A=args.annot_a,
        annot_B=args.annot_b,
        fmt_A=args.fmt_a,
        fmt_B=args.fmt_b,
        m6a_A_file=args.m6a_a,
        m6a_B_file=args.m6a_b,
        ortholog_file=args.ortholog_file,
        output_file=args.output,
        orthologs_only=not args.include_inparalogs,
        verbose=True,
    )


def run_summarize_conservation(args):
    summarize_conserved_m6a = _require_import(
        lambda: __import__(
            "plantm6a.analysis.conservation", fromlist=["summarize_conserved_m6a"]
        ).summarize_conserved_m6a,
        "plantm6a[conservation]",
    )
    summarize_conserved_m6a(args.conserved_dir, args.output_dir, verbose=True)


def run_plot_chord(args):
    plot_chord = _require_import(
        lambda: __import__(
            "plantm6a.analysis.conservation", fromlist=["plot_chord"]
        ).plot_chord,
        "plantm6a[visualization]",
    )
    pd = _require_import(lambda: __import__("pandas"), "plantm6a[visualization]")
    df = pd.read_csv(args.input, sep="\t")
    plot_chord(df, args.output, min_sites=args.min_sites)


def run_metagene_bin(args):
    run_region2bin = _require_import(
        lambda: __import__(
            "plantm6a.analysis.metagene", fromlist=["run_region2bin"]
        ).run_region2bin,
        "plantm6a[visualization]",
    )
    run_region2bin(
        input_bed=args.input,
        annotation_bed12=args.annotation,
        output_file=args.output,
        keep_tmp=args.keep_tmp,
        strand=args.strand,
        len_scale=args.len_scale,
        loc=args.loc,
        pct=args.pct,
        rpm=args.rpm,
        bin_sum=args.bin,
        bin_numbers=args.bin_numbers,
        bin_output_file=args.bin_output,
    )


def run_metagene_plot(args):
    plot_metagene = _require_import(
        lambda: __import__(
            "plantm6a.analysis.metagene", fromlist=["plot_metagene"]
        ).plot_metagene,
        "plantm6a[visualization]",
    )
    plot_metagene(args.input, args.output_dir, args.labels, adjust=args.adjust)


def run_ejc_triplet(args):
    run_ejc_triplet_analysis = _require_import(
        lambda: __import__(
            "plantm6a.analysis.ejc", fromlist=["run_ejc_triplet"]
        ).run_ejc_triplet,
        "plantm6a",
    )
    run_ejc_triplet_analysis(
        annotation_bed12=args.annotation,
        m6a_bed6=args.m6a,
        output_prefix=args.output,
        window_size=args.window,
        upstream_range=args.upstream,
        downstream_range=args.downstream,
        smooth_zero=args.smooth_zero,
    )


def run_ejc_batch(args):
    ejc_batch = _require_import(
        lambda: __import__(
            "plantm6a.analysis.ejc", fromlist=["run_ejc_batch", "parse_species_map_file", "DEFAULT_EJC_SPECIES_MAP"]
        ),
        "plantm6a",
    )
    species_map = ejc_batch.DEFAULT_EJC_SPECIES_MAP
    if args.map_file:
        species_map = ejc_batch.parse_species_map_file(args.map_file)
    result = ejc_batch.run_ejc_batch(
        bed12_dir=args.bed12_dir,
        bed6_dir=args.bed6_dir,
        output_dir=args.output_dir,
        species_map=species_map,
        species=args.species,
        window_size=args.window,
        upstream_range=args.upstream,
        downstream_range=args.downstream,
        smooth_zero=args.smooth_zero,
        strict=args.strict,
    )
    if args.summary:
        import csv

        Path(args.summary).parent.mkdir(parents=True, exist_ok=True)
        with open(args.summary, "w", newline="") as out:
            writer = csv.writer(out, delimiter="\t")
            writer.writerow(["species", "status", "internal_file", "last_file", "total_genes", "internal_junctions", "last_junctions", "skipped_single_exon"])
            for sp, info in sorted(result.items()):
                stats = info.get("stats", {})
                writer.writerow([
                    sp,
                    info.get("status", "unknown"),
                    info.get("internal_file", ""),
                    info.get("last_file", ""),
                    stats.get("total_genes", ""),
                    stats.get("internal_junctions", ""),
                    stats.get("last_junctions", ""),
                    stats.get("skipped_single_exon", ""),
                ])


def run_ejc_plot(args):
    plot_ejc_triplet = _require_import(
        lambda: __import__(
            "plantm6a.analysis.ejc", fromlist=["plot_ejc_triplet"]
        ).plot_ejc_triplet,
        "plantm6a[visualization]",
    )
    plot_ejc_triplet(
        internal_file=args.internal,
        last_file=args.last,
        species_name=args.species_name,
        output_dir=args.output_dir,
        x_limit=args.x_limit,
    )


def run_diff_m6a(args):
    run_differential_m6a = _require_import(
        lambda: __import__(
            "plantm6a.analysis.differential", fromlist=["run_differential_m6a"]
        ).run_differential_m6a,
        "plantm6a[visualization]",
    )
    run_differential_m6a(
        wt_bed=args.wt,
        mut_bed=args.mut,
        annotation=args.annotation,
        output_dir=args.output,
        wt_name=args.wt_name,
        mut_name=args.mut_name,
        ratio_threshold=args.ratio_threshold,
        fold_change_threshold=args.fc_threshold,
        make_plots=not args.no_plots,
        verbose=True,
    )


def build_parser():
    parser = argparse.ArgumentParser(
        description="PlantM6A: a toolkit for plant m6A RNA modification analysis"
    )
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    stats_parser = subparsers.add_parser(
        "stats", help="Statistical analysis for a single species or sample"
    )
    stats_parser.add_argument("--genome", required=True, help="Genome FASTA file")
    stats_parser.add_argument("--gtf", required=True, help="Annotation GTF/GFF file")
    stats_parser.add_argument("--bed", required=True, help="m6A sites BED file")
    stats_parser.add_argument("--output", help="Optional key-value TSV output")
    stats_parser.set_defaults(func=run_stats)

    batch_parser = subparsers.add_parser(
        "batch", help="Batch statistical analysis from a YAML config"
    )
    batch_parser.add_argument("--config", required=True, help="Configuration YAML file")
    batch_parser.add_argument("--output", help="Output summary TSV file")
    batch_parser.set_defaults(func=run_batch)

    motif_parser = subparsers.add_parser(
        "motif", help="Analyze sequence motifs around m6A sites"
    )
    motif_parser.add_argument("--genome", required=True, help="Genome FASTA file")
    motif_parser.add_argument("--annotation", required=True, help="Annotation GTF/GFF file")
    motif_parser.add_argument("--bed", required=True, help="m6A sites BED file")
    motif_parser.add_argument(
        "--mode", choices=["simple", "complete"], default="simple",
        help="Motif mode: simple RAC/GAT/Others or complete 16 centered-A 3-mers",
    )
    motif_parser.add_argument("--output", help="Optional key-value TSV output")
    motif_parser.set_defaults(func=run_motif)

    annotate_parser = subparsers.add_parser(
        "annotate-sites", help="Annotate each BED m6A site with region, motif, gene, strand, and ratio"
    )
    annotate_parser.add_argument("--genome", required=True, help="Genome FASTA file")
    annotate_parser.add_argument("--annotation", required=True, help="Annotation GTF/GFF file")
    annotate_parser.add_argument("--bed", required=True, help="m6A sites BED file")
    annotate_parser.add_argument("--output", required=True, help="Annotated per-site TSV output")
    annotate_parser.add_argument(
        "--fmt", choices=["auto", "gtf", "gff"], default="auto",
        help="Annotation format; default auto-detects from file extension",
    )
    annotate_parser.add_argument(
        "--overlap-policy", choices=["best", "all"], default="best",
        help="Use one best transcript per BED site or emit all overlapping transcript annotations",
    )
    annotate_parser.add_argument(
        "--motif-window", type=int, default=2,
        help="Bases on each side of the site for motif context; 2 gives 5 bp context",
    )
    annotate_parser.set_defaults(func=run_annotate_sites)

    translate_parser = subparsers.add_parser(
        "translate-cds", help="Translate CDS FASTA records to protein FASTA"
    )
    translate_parser.add_argument("--input", required=True, help="Input CDS FASTA")
    translate_parser.add_argument("--output", required=True, help="Output protein FASTA")
    translate_parser.set_defaults(func=run_translate_cds)

    ortholog_parser = subparsers.add_parser(
        "extract-orthologs", help="Extract pairwise gene tables from OrthoFinder output"
    )
    ortholog_parser.add_argument("--orthogroups", required=True, help="Orthogroups.tsv")
    ortholog_parser.add_argument(
        "--orthologues-dir", required=True, help="OrthoFinder Orthologues directory"
    )
    ortholog_parser.add_argument("--output-dir", required=True, help="Output directory")
    ortholog_parser.set_defaults(func=run_extract_orthologs)

    conserve_parser = subparsers.add_parser(
        "conserve-pair", help="Identify conserved m6A sites for one species pair"
    )
    conserve_parser.add_argument("--sp-a", required=True, help="Species A label")
    conserve_parser.add_argument("--sp-b", required=True, help="Species B label")
    conserve_parser.add_argument("--genome-a", required=True, help="Species A genome FASTA")
    conserve_parser.add_argument("--genome-b", required=True, help="Species B genome FASTA")
    conserve_parser.add_argument("--annot-a", required=True, help="Species A annotation")
    conserve_parser.add_argument("--annot-b", required=True, help="Species B annotation")
    conserve_parser.add_argument("--fmt-a", default="gtf", choices=["gtf", "gff"], help="Species A annotation format")
    conserve_parser.add_argument("--fmt-b", default="gtf", choices=["gtf", "gff"], help="Species B annotation format")
    conserve_parser.add_argument("--m6a-a", required=True, help="Species A m6A BED file")
    conserve_parser.add_argument("--m6a-b", required=True, help="Species B m6A BED file")
    conserve_parser.add_argument("--ortholog-file", required=True, help="Pairwise ortholog TSV")
    conserve_parser.add_argument("--output", required=True, help="Output conserved m6A TSV")
    conserve_parser.add_argument(
        "--include-inparalogs",
        action="store_true",
        help="Include inparalog_og rows in addition to direct ortholog rows",
    )
    conserve_parser.set_defaults(func=run_conserve_pair)

    summarize_parser = subparsers.add_parser(
        "summarize-conservation", help="Summarize conserved m6A result files"
    )
    summarize_parser.add_argument("--conserved-dir", required=True, help="Directory with *_conserved_m6A.tsv files")
    summarize_parser.add_argument("--output-dir", required=True, help="Summary output directory")
    summarize_parser.set_defaults(func=run_summarize_conservation)

    plot_parser = subparsers.add_parser(
        "plot-chord", help="Plot a chord diagram from pair_region_breakdown.tsv"
    )
    plot_parser.add_argument("--input", required=True, help="pair_region_breakdown.tsv")
    plot_parser.add_argument("--output", required=True, help="Output prefix for PDF/SVG")
    plot_parser.add_argument("--min-sites", type=int, default=100, help="Minimum sites per species pair")
    plot_parser.set_defaults(func=run_plot_chord)

    metagene_bin_parser = subparsers.add_parser(
        "metagene-bin", help="Map BED6 peaks to BED12 transcript regions and write metagene bins"
    )
    metagene_bin_parser.add_argument("--input", required=True, help="Input BED6 peak file")
    metagene_bin_parser.add_argument("--annotation", required=True, help="BED12 transcript annotation")
    metagene_bin_parser.add_argument("--output", required=True, help="Output region-bin table")
    metagene_bin_parser.add_argument("--keep-tmp", action="store_true", help="Keep intermediate intersectBed output")
    metagene_bin_parser.add_argument("--strand", action="store_true", help="Use strand-specific bedtools intersection")
    metagene_bin_parser.add_argument("--len-scale", action="store_true", help="Scale region widths by average region length")
    metagene_bin_parser.add_argument("--loc", action="store_true", help="Count peak locations instead of peak values")
    metagene_bin_parser.add_argument("--pct", action="store_true", help="Normalize by total peak count")
    metagene_bin_parser.add_argument("--rpm", action="store_true", help="RPM normalization; mutually exclusive with --loc/--pct")
    metagene_bin_parser.add_argument("--bin", type=int, default=100, help="Number of bins per region")
    metagene_bin_parser.add_argument("--bin-numbers", help="Comma-separated bin indexes to write extra details for")
    metagene_bin_parser.add_argument("--bin-output", help="Output file for --bin-numbers details")
    metagene_bin_parser.set_defaults(func=run_metagene_bin)

    metagene_plot_parser = subparsers.add_parser(
        "metagene-plot", help="Draw metagene profile from a region-bin table"
    )
    metagene_plot_parser.add_argument("--input", required=True, help="region-bin table from metagene-bin")
    metagene_plot_parser.add_argument("--output-dir", required=True, help="Output directory for PDF/PNG/SVG")
    metagene_plot_parser.add_argument("--labels", required=True, help="Sample labels, comma-separated for multiple columns")
    metagene_plot_parser.add_argument("--adjust", type=float, default=0.5, help="KDE bandwidth adjustment")
    metagene_plot_parser.set_defaults(func=run_metagene_plot)

    ejc_parser = subparsers.add_parser(
        "ejc-triplet", help="Analyze m6A stacking around exon-junction-exon triplets"
    )
    ejc_parser.add_argument("--annotation", "-a", required=True, help="Gene annotation in BED12 format")
    ejc_parser.add_argument("--m6a", "-m", required=True, help="m6A sites in BED6 format")
    ejc_parser.add_argument("--output", "-o", required=True, help="Output prefix")
    ejc_parser.add_argument("--window", "-w", type=int, default=10, help="Window size in bp")
    ejc_parser.add_argument("--upstream", "-u", type=int, default=500, help="Upstream range")
    ejc_parser.add_argument("--downstream", "-d", type=int, default=500, help="Downstream range")
    ejc_parser.add_argument("--smooth-zero", action="store_true", help="Replace bin 0 with mean of bins -1 and +1")
    ejc_parser.set_defaults(func=run_ejc_triplet)

    ejc_batch_parser = subparsers.add_parser(
        "ejc-batch", help="Run exon-junction triplet analysis across multiple species"
    )
    ejc_batch_parser.add_argument("--bed12-dir", required=True, help="Directory containing species BED12 annotations")
    ejc_batch_parser.add_argument("--bed6-dir", required=True, help="Directory containing species BED6 m6A files")
    ejc_batch_parser.add_argument("--output-dir", required=True, help="Directory for per-species output prefixes")
    ejc_batch_parser.add_argument("--map-file", help="Optional TSV: species<TAB>bed12_file<TAB>bed6_file")
    ejc_batch_parser.add_argument("--species", nargs="+", help="Optional subset of species names to run")
    ejc_batch_parser.add_argument("--window", type=int, default=50, help="Window size in bp")
    ejc_batch_parser.add_argument("--upstream", type=int, default=500, help="Upstream range")
    ejc_batch_parser.add_argument("--downstream", type=int, default=500, help="Downstream range")
    ejc_batch_parser.add_argument("--smooth-zero", action="store_true", help="Replace bin 0 with mean of bins -1 and +1")
    ejc_batch_parser.add_argument("--strict", action="store_true", help="Fail on missing annotation or m6A file instead of skipping")
    ejc_batch_parser.add_argument("--summary", help="Optional TSV summary output")
    ejc_batch_parser.set_defaults(func=run_ejc_batch)

    ejc_plot_parser = subparsers.add_parser(
        "ejc-plot", help="Draw SVG visualizations from EJC internal/last junction outputs"
    )
    ejc_plot_parser.add_argument("--internal", required=True, help="*_internal_junctions.txt file")
    ejc_plot_parser.add_argument("--last", required=True, help="*_last_junctions.txt file")
    ejc_plot_parser.add_argument("--species-name", required=True, help="Species label used in titles and output filenames")
    ejc_plot_parser.add_argument("--output-dir", required=True, help="Directory for SVG outputs")
    ejc_plot_parser.add_argument("--x-limit", type=int, default=1000, help="Plot distance range on each side of junction")
    ejc_plot_parser.set_defaults(func=run_ejc_plot)

    diff_parser = subparsers.add_parser(
        "diff-m6a", help="Compare differential m6A sites between WT and mutant/condition BED6 files"
    )
    diff_parser.add_argument("--wt", required=True, help="WT/control m6A BED6 file")
    diff_parser.add_argument("--mut", required=True, help="Mutant/treatment m6A BED6 file")
    diff_parser.add_argument("--annotation", "--gtf", required=True, help="GTF/GFF annotation file")
    diff_parser.add_argument("--output", default="results/diff_m6a", help="Output directory")
    diff_parser.add_argument("--wt-name", default="WT", help="WT/control sample name")
    diff_parser.add_argument("--mut-name", default="Mutant", help="Mutant/treatment sample name")
    diff_parser.add_argument("--ratio-threshold", type=float, default=0.1, help="Minimum ratio difference for hyper/hypo calls")
    diff_parser.add_argument("--fc-threshold", type=float, default=1.5, help="Minimum fold-change for hyper/hypo calls")
    diff_parser.add_argument("--no-plots", action="store_true", help="Skip figure generation")
    diff_parser.set_defaults(func=run_diff_m6a)

    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if not hasattr(args, "func"):
        parser.print_help()
        return 1

    try:
        args.func(args)
    except FileNotFoundError as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 2
    except KeyError as exc:
        print(f"Error: missing expected field {exc}", file=sys.stderr)
        return 2

    return 0


if __name__ == "__main__":
    sys.exit(main())
