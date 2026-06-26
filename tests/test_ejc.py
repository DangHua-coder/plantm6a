from plantm6a.analysis.ejc import ExonJunctionTripletAnalyzer, plot_ejc_triplet, run_ejc_batch, run_ejc_triplet


def write_ejc_fixture(tmp_path):
    bed12 = tmp_path / "annotation.bed12"
    bed6 = tmp_path / "sites.bed6"
    bed12.write_text(
        "chr1\t0\t250\ttx1\t0\t+\t0\t250\t0\t3\t50,50,50,\t0,100,200,\n"
        "chr1\t300\t460\ttx2\t0\t-\t300\t460\t0\t2\t60,60,\t0,100,\n"
        "chr1\t500\t550\tsingle\t0\t+\t500\t550\t0\t1\t50,\t0,\n"
    )
    bed6.write_text(
        "chr1\t45\t46\tsite1\t1\t+\n"
        "chr1\t105\t106\tsite2\t1\t+\n"
        "chr1\t205\t206\tsite3\t1\t+\n"
        "chr1\t355\t356\tsite4\t1\t-\n"
    )
    return bed12, bed6


def test_ejc_triplet_bin_mapping_uses_floor():
    analyzer = ExonJunctionTripletAnalyzer(window_size=50)

    assert analyzer.map_position_to_bin(99, [100, 150], "+") == -1
    assert analyzer.map_position_to_bin(150, [100, 150], "+") == 0
    assert analyzer.map_position_to_bin(99, [100, 150], "-") == 0


def test_ejc_batch_skips_missing_species(tmp_path):
    bed12_dir = tmp_path / "bed12"
    bed6_dir = tmp_path / "bed6"
    output_dir = tmp_path / "out"
    bed12_dir.mkdir()
    bed6_dir.mkdir()
    (bed12_dir / "a.bed12").write_text("chr1\t0\t250\ttx1\t0\t+\t0\t250\t0\t2\t50,50,\t0,100,\n")
    (bed6_dir / "a.bed6").write_text("chr1\t45\t46\tsite1\t1\t+\n")

    result = run_ejc_batch(
        str(bed12_dir),
        str(bed6_dir),
        str(output_dir),
        species_map={"ok": ("a.bed12", "a.bed6"), "missing": ("missing.bed12", "missing.bed6")},
        species=["ok", "missing"],
        window_size=50,
        smooth_zero=True,
    )

    assert result["ok"]["status"] == "completed"
    assert result["missing"]["status"] == "missing_annotation"
    assert (output_dir / "ok_internal_junctions.txt").exists()


def test_plot_ejc_triplet_writes_svg_outputs(tmp_path):
    header = "bin\tdistance_from_junction\tm6a_level\tcoverage\tm6a_density\tm6a_likelihood\n"
    rows = "".join(
        f"{i}\t{i * 50}\t{10 + i}\t100\t0.2\t0.1\n"
        for i in range(-3, 4)
    )
    internal = tmp_path / "internal.txt"
    last = tmp_path / "last.txt"
    internal.write_text(header + rows)
    last.write_text(header + rows)
    output_dir = tmp_path / "figures"

    outputs = plot_ejc_triplet(str(internal), str(last), "Ara", str(output_dir), x_limit=150)

    assert len(outputs) == 5
    assert (output_dir / "Ara_junction_m6a_level.svg").exists()
    assert (output_dir / "Ara_junction_m6a_likelihood.svg").exists()
    assert (output_dir / "Ara_junction_combined.svg").exists()
    assert (output_dir / "Ara_level_comparison.svg").exists()
    assert (output_dir / "Ara_likelihood_comparison.svg").exists()


def test_ejc_triplet_analysis_writes_outputs(tmp_path):
    bed12, bed6 = write_ejc_fixture(tmp_path)
    output_prefix = tmp_path / "ejc"

    result = run_ejc_triplet(str(bed12), str(bed6), str(output_prefix), window_size=50)

    internal = tmp_path / "ejc_internal_junctions.txt"
    last = tmp_path / "ejc_last_junctions.txt"
    assert result["internal_file"] == str(internal)
    assert result["last_file"] == str(last)
    assert result["stats"] == {
        "total_genes": 3,
        "internal_junctions": 1,
        "last_junctions": 2,
        "skipped_single_exon": 1,
    }
    assert internal.exists()
    assert last.exists()
    lines = internal.read_text().splitlines()
    assert lines[0] == "bin\tdistance_from_junction\tm6a_level\tcoverage\tm6a_density\tm6a_likelihood"
    minus_one = lines[1].split("\t")
    assert minus_one[:3] == ["-1", "-50", "1"]
    assert float(minus_one[3]) == 1.0000000000000004
    assert minus_one[4:] == ["0.020000", "1.000000"]
