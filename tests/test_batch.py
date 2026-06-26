import pytest

pd = pytest.importorskip("pandas")
pytest.importorskip("yaml")

from plantm6a.analysis.batch import batch_analyze


def test_batch_analyze_writes_summary(tmp_path):
    genome = tmp_path / "genome.fa"
    gtf = tmp_path / "annotation.gtf"
    bed = tmp_path / "sites.bed"
    config = tmp_path / "config.yaml"
    output = tmp_path / "summary.tsv"

    genome.write_text(">chr1\nAAAACCCCGGGG\n")
    gtf.write_text(
        "chr1\ttest\texon\t1\t4\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";\n"
    )
    bed.write_text("chr1\t1\t2\tsite1\n")
    config.write_text(
        f"species:\n"
        f"  - name: tiny\n"
        f"    genome: {genome}\n"
        f"    annotation: {gtf}\n"
        f"    bed: {bed}\n"
    )

    batch_analyze(str(config), str(output))

    df = pd.read_csv(output, sep="\t")
    assert list(df["Species"]) == ["tiny"]
    assert df.loc[0, "Exon_A_Count"] == 4
    assert df.loc[0, "Exon_m6A_Sites"] == 1
    assert "m6A_per_1000A" in df.columns
