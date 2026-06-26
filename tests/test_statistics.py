from plantm6a.analysis.statistics import analyze, merge_overlapping_intervals


def test_merge_overlapping_intervals():
    assert merge_overlapping_intervals([(5, 8), (1, 3), (2, 6), (10, 12)]) == [
        (1, 8),
        (10, 12),
    ]


def test_analyze_tiny_fixture(tmp_path):
    genome = tmp_path / "genome.fa"
    gtf = tmp_path / "annotation.gtf"
    bed = tmp_path / "sites.bed"

    genome.write_text(">chr1\nAAAACCCCGGGG\n")
    gtf.write_text(
        "chr1\ttest\texon\t1\t4\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";\n"
        "chr1\ttest\texon\t3\t8\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";\n"
    )
    bed.write_text("chr1\t1\t2\tsite1\nchr1\t9\t10\tsite2\n")

    result = analyze(str(genome), str(gtf), str(bed), verbose=False)

    assert result["exon_bases"] == 8
    assert result["exon_A_count"] == 4
    assert result["total_m6A_sites"] == 2
    assert result["exon_m6A_sites"] == 1
    assert result["non_exon_sites"] == 1
    assert result["m6A_per_1000A"] == 250
