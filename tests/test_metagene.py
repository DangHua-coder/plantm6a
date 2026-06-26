from plantm6a.analysis.metagene import (
    RegionBinCalculator,
    default_weights,
    plot_metagene,
    process_cds_region,
    process_intersection,
    process_ncrna_region,
    process_utr3_region,
    process_utr5_region,
    write_output,
)


def test_region_processors_assign_expected_bins():
    exon_lengths = [100, 100]
    exon_starts = [0, 200]

    assert process_utr5_region(25, 0, "+", 50, 2, exon_lengths, exon_starts, 100) == (25, 50, "UTR5", False)
    assert process_cds_region(75, 0, "+", 50, 250, 2, exon_lengths, exon_starts, 100) == (25, 25, "CDS", False)
    assert process_utr3_region(275, 0, "+", 250, 2, exon_lengths, exon_starts, 100) == (26, 50, "UTR3", False)
    assert process_ncrna_region(25, 0, "+", 2, exon_lengths, exon_starts, 100) == (25, 13, "exonFirst", False)


def test_process_intersection_and_write_output(tmp_path):
    tmp_intersect = tmp_path / "intersect.tmp"
    tmp_intersect.write_text(
        "chr1\t24\t25\tsite1\t2\t+\tchr1\t0\t300\ttx1\t0\t+\t50\t250\t0\t2\t100,100,\t0,200,\n"
        "chr1\t74\t75\tsite2\t3\t+\tchr1\t0\t300\ttx1\t0\t+\t50\t250\t0\t2\t100,100,\t0,200,\n"
    )
    calculator = RegionBinCalculator(bin_sum=100, pos_num={"chr125+": 1, "chr175+": 1})

    process_intersection(str(tmp_intersect), calculator)

    assert calculator.region_bin["UTR5"][50] == 2
    assert calculator.region_bin["CDS"][25] == 3

    output = tmp_path / "bins.tsv"
    write_output(str(output), calculator, default_weights(), 100)
    lines = output.read_text().splitlines()
    assert len(lines) == 800
    assert lines[49] == "UTR5\t50\t2.0"
    assert lines[124] == "CDS\t125\t3.0"


def test_plot_metagene_writes_three_formats(tmp_path):
    bin_file = tmp_path / "bins.tsv"
    lines = []
    for region in ["UTR5", "CDS", "UTR3", "exonFirst", "exonInternal", "exonLast", "mRNAIntron", "ncRNAIntron"]:
        for i in range(1, 101):
            value = 5 if i in (20, 50, 80) else 0
            lines.append(f"{region}\t{i}\t{value}")
    bin_file.write_text("\n".join(lines) + "\n")
    output_dir = tmp_path / "figures"

    outputs = plot_metagene(str(bin_file), str(output_dir), "sample1")

    assert len(outputs) == 3
    assert (output_dir / "sample1.pdf").exists()
    assert (output_dir / "sample1.png").exists()
    assert (output_dir / "sample1.svg").exists()
