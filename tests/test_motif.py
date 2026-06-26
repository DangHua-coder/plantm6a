import pytest


def test_generate_3mers_with_center_a():
    pysam = pytest.importorskip("pysam")
    assert pysam is not None

    from plantm6a.analysis import generate_all_3mers_with_A

    motifs = generate_all_3mers_with_A()

    assert len(motifs) == 16
    assert motifs == sorted(motifs)
    assert "AAC" in motifs
    assert "GAT" in motifs
    assert all(len(motif) == 3 and motif[1] == "A" for motif in motifs)


def test_motif_simple_mode_on_tiny_fixture(tmp_path):
    pytest.importorskip("pysam")

    from plantm6a.analysis import analyze_motifs

    genome = tmp_path / "genome.fa"
    gtf = tmp_path / "annotation.gtf"
    bed = tmp_path / "sites.bed"

    genome.write_text(">chr1\nTTGAACCTAGATTT\n")
    gtf.write_text(
        "chr1\ttest\texon\t1\t14\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";\n"
    )
    bed.write_text("chr1\t4\t5\tsite1\t1\t+\nchr1\t10\t11\tsite2\t1\t+\n")

    result = analyze_motifs(
        genome_path=str(genome),
        annotation_path=str(gtf),
        bed_path=str(bed),
        mode="simple",
        verbose=False,
    )

    assert result["total_sites"] == 2
    assert result["unique_sites"] == 2
    assert result["valid_motifs"] == 2
    assert result["RAC_count"] == 1
    assert result["GAT_count"] == 1
    assert result["others_count"] == 0
