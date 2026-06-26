import pytest

from plantm6a.cli import main


def write_stats_fixture(tmp_path):
    genome = tmp_path / "genome.fa"
    gtf = tmp_path / "annotation.gtf"
    bed = tmp_path / "sites.bed"

    genome.write_text(">chr1\nAAAACCCCGGGG\n")
    gtf.write_text(
        "chr1\ttest\texon\t1\t4\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";\n"
    )
    bed.write_text("chr1\t1\t2\tsite1\nchr1\t9\t10\tsite2\n")
    return genome, gtf, bed


def write_motif_fixture(tmp_path):
    genome = tmp_path / "motif.fa"
    gtf = tmp_path / "motif.gtf"
    bed = tmp_path / "motif.bed"

    genome.write_text(">chr1\nTTGAACCTAGATTT\n")
    gtf.write_text(
        "chr1\ttest\texon\t1\t14\t.\t+\t.\tgene_id \"g1\"; gene_name \"G1\"; transcript_id \"t1\";\n"
    )
    bed.write_text("chr1\t4\t5\tsite1\t0.1\t+\nchr1\t10\t11\tsite2\t0.2\t+\n")
    return genome, gtf, bed


def test_cli_help(capsys):
    with pytest.raises(SystemExit) as excinfo:
        main(["--help"])
    assert excinfo.value.code == 0
    captured = capsys.readouterr()
    assert "PlantM6A" in captured.out
    assert "stats" in captured.out
    assert "motif" in captured.out


def test_stats_help(capsys):
    with pytest.raises(SystemExit) as excinfo:
        main(["stats", "--help"])
    assert excinfo.value.code == 0
    captured = capsys.readouterr()
    assert "--genome" in captured.out
    assert "--gtf" in captured.out
    assert "--bed" in captured.out


def test_annotate_sites_help(capsys):
    with pytest.raises(SystemExit) as excinfo:
        main(["annotate-sites", "--help"])
    assert excinfo.value.code == 0
    captured = capsys.readouterr()
    assert "--genome" in captured.out
    assert "--annotation" in captured.out
    assert "--overlap-policy" in captured.out


def test_ejc_triplet_help(capsys):
    with pytest.raises(SystemExit) as excinfo:
        main(["ejc-triplet", "--help"])
    assert excinfo.value.code == 0
    captured = capsys.readouterr()
    assert "--annotation" in captured.out
    assert "--smooth-zero" in captured.out

    with pytest.raises(SystemExit) as excinfo:
        main(["diff-m6a", "--help"])
    assert excinfo.value.code == 0
    captured = capsys.readouterr()
    assert "--wt" in captured.out
    assert "--fc-threshold" in captured.out

    with pytest.raises(SystemExit) as excinfo:
        main(["ejc-batch", "--help"])
    assert excinfo.value.code == 0
    captured = capsys.readouterr()
    assert "--bed12-dir" in captured.out
    assert "--map-file" in captured.out

    with pytest.raises(SystemExit) as excinfo:
        main(["ejc-plot", "--help"])
    assert excinfo.value.code == 0
    captured = capsys.readouterr()
    assert "--internal" in captured.out
    assert "--species-name" in captured.out


def test_metagene_commands_help(capsys):
    with pytest.raises(SystemExit) as excinfo:
        main(["metagene-bin", "--help"])
    assert excinfo.value.code == 0
    captured = capsys.readouterr()
    assert "--annotation" in captured.out
    assert "--len-scale" in captured.out

    with pytest.raises(SystemExit) as excinfo:
        main(["metagene-plot", "--help"])
    assert excinfo.value.code == 0
    captured = capsys.readouterr()
    assert "--labels" in captured.out
    assert "--output-dir" in captured.out


def test_stats_command_writes_key_value_output(tmp_path):
    genome, gtf, bed = write_stats_fixture(tmp_path)
    output = tmp_path / "stats.tsv"

    exit_code = main([
        "stats",
        "--genome", str(genome),
        "--gtf", str(gtf),
        "--bed", str(bed),
        "--output", str(output),
    ])

    assert exit_code == 0
    lines = output.read_text().splitlines()
    assert "metric\tvalue" in lines
    assert "exon_A_count\t4" in lines
    assert "exon_m6A_sites\t1" in lines


def test_batch_command_writes_summary(tmp_path):
    pytest.importorskip("pandas")
    pytest.importorskip("yaml")
    genome, gtf, bed = write_stats_fixture(tmp_path)
    config = tmp_path / "config.yaml"
    output = tmp_path / "batch.tsv"
    config.write_text(
        f"species:\n"
        f"  - name: tiny\n"
        f"    genome: {genome}\n"
        f"    annotation: {gtf}\n"
        f"    bed: {bed}\n"
    )

    exit_code = main(["batch", "--config", str(config), "--output", str(output)])

    assert exit_code == 0
    lines = output.read_text().splitlines()
    assert lines[0].startswith("Species\tExon_Bases")
    assert lines[1].startswith("tiny\t4\t4")


def test_motif_command_writes_key_value_output(tmp_path):
    pytest.importorskip("pysam")
    genome, gtf, bed = write_motif_fixture(tmp_path)
    output = tmp_path / "motif.tsv"

    exit_code = main([
        "motif",
        "--genome", str(genome),
        "--annotation", str(gtf),
        "--bed", str(bed),
        "--mode", "simple",
        "--output", str(output),
    ])

    assert exit_code == 0
    lines = output.read_text().splitlines()
    assert "RAC_count\t1" in lines
    assert "GAT_count\t1" in lines


def test_annotate_sites_command_writes_table(tmp_path):
    pytest.importorskip("pysam")
    genome, gtf, bed = write_motif_fixture(tmp_path)
    output = tmp_path / "annotated.tsv"

    exit_code = main([
        "annotate-sites",
        "--genome", str(genome),
        "--annotation", str(gtf),
        "--bed", str(bed),
        "--output", str(output),
    ])

    assert exit_code == 0
    lines = output.read_text().splitlines()
    assert lines[0].startswith("site_index\tchrom\tstart\tend")
    assert "site1" in lines[1]
    assert "RAC" in lines[1]


def test_translate_cds_command_writes_protein_fasta(tmp_path):
    cds = tmp_path / "cds.fa"
    protein = tmp_path / "protein.fa"
    cds.write_text(">gene1\nATGGCTTAA\n")

    exit_code = main(["translate-cds", "--input", str(cds), "--output", str(protein)])

    assert exit_code == 0
    assert protein.read_text() == ">gene1\nMA\n"


def test_extract_orthologs_command_writes_pair_tables(tmp_path):
    orthogroups = tmp_path / "Orthogroups.tsv"
    orthologues_dir = tmp_path / "Orthologues"
    pair_dir = orthologues_dir / "Orthologues_Ara"
    output_dir = tmp_path / "pairs"
    pair_dir.mkdir(parents=True)
    orthogroups.write_text("Orthogroup\tAra\tRice\nOG0001\tAra1\tRice1\n")
    (pair_dir / "Ara__v__Rice.tsv").write_text("Orthogroup\tAra\tRice\nOG0001\tAra1\tRice1\n")

    exit_code = main([
        "extract-orthologs",
        "--orthogroups", str(orthogroups),
        "--orthologues-dir", str(orthologues_dir),
        "--output-dir", str(output_dir),
    ])

    assert exit_code == 0
    assert (output_dir / "Ara_vs_Rice.tsv").exists()
    assert (output_dir / "summary_statistics.tsv").exists()


def test_summarize_conservation_command_writes_outputs(tmp_path):
    conserved_dir = tmp_path / "conserved"
    output_dir = tmp_path / "summary"
    conserved_dir.mkdir()
    (conserved_dir / "Ara_vs_Rice_conserved_m6A.tsv").write_text(
        "OG_ID\tAra_transcript\tRice_transcript\tAra_tx_pos\tRice_tx_pos\t"
        "Ara_region\tRice_region\tAra_m6a_score\tRice_m6a_score\t"
        "Ara_aln_col\tRice_aln_col\tcol_diff\taln_score\taln_len\tAra_context\tRice_context\n"
        "OG0001\ttA\ttR\t5\t6\tcds\tcds\t0.7\t0.8\t10\t11\t1\t12.0\t20\tAAAA\tAAAA\n"
    )

    exit_code = main([
        "summarize-conservation",
        "--conserved-dir", str(conserved_dir),
        "--output-dir", str(output_dir),
    ])

    assert exit_code == 0
    assert (output_dir / "pair_summary.tsv").exists()
    assert (output_dir / "region_distribution.tsv").exists()


def test_metagene_bin_command_passes_expected_arguments(monkeypatch, tmp_path):
    import plantm6a.analysis.metagene as metagene

    calls = []

    def fake_run_region2bin(**kwargs):
        calls.append(kwargs)

    monkeypatch.setattr(metagene, "run_region2bin", fake_run_region2bin)
    output = tmp_path / "bins.tsv"

    exit_code = main([
        "metagene-bin",
        "--input", "sites.bed6",
        "--annotation", "tx.bed12",
        "--output", str(output),
        "--strand",
        "--len-scale",
        "--loc",
        "--bin", "50",
        "--bin-numbers", "1,2",
        "--bin-output", "details.tsv",
    ])

    assert exit_code == 0
    assert calls == [{
        "input_bed": "sites.bed6",
        "annotation_bed12": "tx.bed12",
        "output_file": str(output),
        "keep_tmp": False,
        "strand": True,
        "len_scale": True,
        "loc": True,
        "pct": False,
        "rpm": False,
        "bin_sum": 50,
        "bin_numbers": "1,2",
        "bin_output_file": "details.tsv",
    }]


def test_metagene_plot_command_passes_expected_arguments(monkeypatch, tmp_path):
    import plantm6a.analysis.metagene as metagene

    calls = []

    def fake_plot_metagene(*args, **kwargs):
        calls.append((args, kwargs))

    monkeypatch.setattr(metagene, "plot_metagene", fake_plot_metagene)
    output_dir = tmp_path / "figures"

    exit_code = main([
        "metagene-plot",
        "--input", "bins.tsv",
        "--output-dir", str(output_dir),
        "--labels", "sampleA,sampleB",
        "--adjust", "0.7",
    ])

    assert exit_code == 0
    assert calls == [(('bins.tsv', str(output_dir), 'sampleA,sampleB'), {"adjust": 0.7})]


def test_diff_m6a_command_passes_expected_arguments(monkeypatch, tmp_path):
    import plantm6a.analysis.differential as differential

    calls = []

    def fake_run_differential_m6a(**kwargs):
        calls.append(kwargs)

    monkeypatch.setattr(differential, "run_differential_m6a", fake_run_differential_m6a)
    output_dir = tmp_path / "diff_out"

    exit_code = main([
        "diff-m6a",
        "--wt", "wt.bed6",
        "--mut", "mut.bed6",
        "--annotation", "annotation.gtf",
        "--output", str(output_dir),
        "--wt-name", "Col",
        "--mut-name", "C11",
        "--ratio-threshold", "0.2",
        "--fc-threshold", "2.0",
        "--no-plots",
    ])

    assert exit_code == 0
    assert calls == [{
        "wt_bed": "wt.bed6",
        "mut_bed": "mut.bed6",
        "annotation": "annotation.gtf",
        "output_dir": str(output_dir),
        "wt_name": "Col",
        "mut_name": "C11",
        "ratio_threshold": 0.2,
        "fold_change_threshold": 2.0,
        "make_plots": False,
        "verbose": True,
    }]


def test_ejc_plot_command_passes_expected_arguments(monkeypatch, tmp_path):
    import plantm6a.analysis.ejc as ejc

    calls = []

    def fake_plot_ejc_triplet(**kwargs):
        calls.append(kwargs)

    monkeypatch.setattr(ejc, "plot_ejc_triplet", fake_plot_ejc_triplet)
    output_dir = tmp_path / "figures"

    exit_code = main([
        "ejc-plot",
        "--internal", "internal.txt",
        "--last", "last.txt",
        "--species-name", "Arabidopsis",
        "--output-dir", str(output_dir),
        "--x-limit", "500",
    ])

    assert exit_code == 0
    assert calls == [{
        "internal_file": "internal.txt",
        "last_file": "last.txt",
        "species_name": "Arabidopsis",
        "output_dir": str(output_dir),
        "x_limit": 500,
    }]


def test_ejc_batch_command_passes_expected_arguments(monkeypatch, tmp_path):
    import plantm6a.analysis.ejc as ejc

    calls = []

    def fake_run_ejc_batch(**kwargs):
        calls.append(kwargs)
        return {"ok": {"status": "completed"}}

    monkeypatch.setattr(ejc, "run_ejc_batch", fake_run_ejc_batch)
    output_dir = tmp_path / "ejc_out"

    exit_code = main([
        "ejc-batch",
        "--bed12-dir", "bed12s",
        "--bed6-dir", "bed6s",
        "--output-dir", str(output_dir),
        "--species", "human", "arabidopsis",
        "--window", "50",
        "--upstream", "600",
        "--downstream", "700",
        "--strict",
        "--summary", str(tmp_path / "summary.tsv"),
    ])

    assert exit_code == 0
    assert calls == [{
        "bed12_dir": "bed12s",
        "bed6_dir": "bed6s",
        "output_dir": str(output_dir),
        "species_map": ejc.DEFAULT_EJC_SPECIES_MAP,
        "species": ["human", "arabidopsis"],
        "window_size": 50,
        "upstream_range": 600,
        "downstream_range": 700,
        "smooth_zero": False,
        "strict": True,
    }]


def test_ejc_triplet_command_passes_expected_arguments(monkeypatch, tmp_path):
    import plantm6a.analysis.ejc as ejc

    calls = []

    def fake_run_ejc_triplet(**kwargs):
        calls.append(kwargs)

    monkeypatch.setattr(ejc, "run_ejc_triplet", fake_run_ejc_triplet)
    output = tmp_path / "ejc_out"

    exit_code = main([
        "ejc-triplet",
        "--annotation", "annotation.bed12",
        "--m6a", "sites.bed6",
        "--output", str(output),
        "--window", "50",
        "--upstream", "600",
        "--downstream", "700",
        "--smooth-zero",
    ])

    assert exit_code == 0
    assert calls == [{
        "annotation_bed12": "annotation.bed12",
        "m6a_bed6": "sites.bed6",
        "output_prefix": str(output),
        "window_size": 50,
        "upstream_range": 600,
        "downstream_range": 700,
        "smooth_zero": True,
    }]


def test_conserve_pair_command_passes_expected_arguments(monkeypatch, tmp_path):
    import plantm6a.analysis.conservation as conservation

    calls = []

    def fake_analyze_species_pair(**kwargs):
        calls.append(kwargs)

    monkeypatch.setattr(conservation, "analyze_species_pair", fake_analyze_species_pair)
    output = tmp_path / "conserved.tsv"

    exit_code = main([
        "conserve-pair",
        "--sp-a", "Ara",
        "--sp-b", "Rice",
        "--genome-a", "ara.fa",
        "--genome-b", "rice.fa",
        "--annot-a", "ara.gtf",
        "--annot-b", "rice.gtf",
        "--fmt-a", "gtf",
        "--fmt-b", "gff",
        "--m6a-a", "ara.bed",
        "--m6a-b", "rice.bed",
        "--ortholog-file", "pairs.tsv",
        "--output", str(output),
        "--include-inparalogs",
    ])

    assert exit_code == 0
    assert calls == [{
        "spA": "Ara",
        "spB": "Rice",
        "genome_A": "ara.fa",
        "genome_B": "rice.fa",
        "annot_A": "ara.gtf",
        "annot_B": "rice.gtf",
        "fmt_A": "gtf",
        "fmt_B": "gff",
        "m6a_A_file": "ara.bed",
        "m6a_B_file": "rice.bed",
        "ortholog_file": "pairs.tsv",
        "output_file": str(output),
        "orthologs_only": False,
        "verbose": True,
    }]
