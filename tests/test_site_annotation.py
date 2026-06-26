import csv

import pytest


def write_fixture(tmp_path):
    genome = tmp_path / "genome.fa"
    gtf = tmp_path / "annotation.gtf"
    bed = tmp_path / "sites.bed"

    seq = list("C" * 120)
    for pos, base in {
        4: "A", 10: "A", 25: "A", 74: "A",
        44: "T", 49: "T", 54: "T", 84: "T",
    }.items():
        seq[pos] = base
    genome.write_text(">chr1\n" + "".join(seq) + "\n")

    gtf.write_text(
        "chr1\ttest\texon\t1\t15\t.\t+\t.\tgene_id \"g_plus\"; gene_name \"PLUS\"; transcript_id \"t_plus\";\n"
        "chr1\ttest\texon\t21\t30\t.\t+\t.\tgene_id \"g_plus\"; gene_name \"PLUS\"; transcript_id \"t_plus\";\n"
        "chr1\ttest\tCDS\t6\t15\t.\t+\t0\tgene_id \"g_plus\"; gene_name \"PLUS\"; transcript_id \"t_plus\";\n"
        "chr1\ttest\tCDS\t21\t25\t.\t+\t0\tgene_id \"g_plus\"; gene_name \"PLUS\"; transcript_id \"t_plus\";\n"
        "chr1\ttest\texon\t41\t60\t.\t-\t.\tgene_id \"g_minus\"; gene_name \"MINUS\"; transcript_id \"t_minus\";\n"
        "chr1\ttest\tCDS\t46\t55\t.\t-\t0\tgene_id \"g_minus\"; gene_name \"MINUS\"; transcript_id \"t_minus\";\n"
        "chr1\ttest\texon\t71\t80\t.\t+\t.\tgene_id \"g_nc\"; gene_name \"NONCODING\"; transcript_id \"t_nc\";\n"
    )

    bed.write_text(
        "chr1\t4\t5\tplus_5utr\t0.10\t+\n"
        "chr1\t10\t11\tplus_cds\t0.20\t+\n"
        "chr1\t25\t26\tplus_3utr\t0.30\t+\n"
        "chr1\t17\t18\tplus_intron\t0.40\t+\n"
        "chr1\t44\t45\tminus_3utr\t0.50\t-\n"
        "chr1\t54\t55\tminus_cds\t0.60\t-\n"
        "chr1\t84\t85\tintergenic\t0.70\t-\n"
        "chr1\t74\t75\tnoncoding\t0.80\t+\n"
    )
    return genome, gtf, bed


def test_annotate_m6a_sites_regions_and_fields(tmp_path):
    pytest.importorskip("pysam")

    from plantm6a.analysis.site_annotation import annotate_m6a_sites

    genome, gtf, bed = write_fixture(tmp_path)
    rows = annotate_m6a_sites(str(genome), str(gtf), str(bed), verbose=False)
    by_name = {row["bed_name"]: row for row in rows}

    assert len(rows) == 8
    assert by_name["plus_5utr"]["region"] == "5utr"
    assert by_name["plus_5utr"]["gene_name"] == "PLUS"
    assert by_name["plus_5utr"]["ratio"] == 0.10
    assert by_name["plus_5utr"]["motif_type"] == "others"
    assert by_name["plus_cds"]["region"] == "cds"
    assert by_name["plus_3utr"]["region"] == "3utr"
    assert by_name["plus_intron"]["region"] == "intron"
    assert by_name["minus_3utr"]["region"] == "3utr"
    assert by_name["minus_cds"]["region"] == "cds"
    assert by_name["minus_cds"]["assigned_strand"] == "-"
    assert by_name["intergenic"]["region"] == "intergenic"
    assert by_name["intergenic"]["transcript_id"] == ""
    assert by_name["noncoding"]["region"] == "noncoding"
    assert by_name["noncoding"]["annotation_status"] == "noncoding"


def test_annotate_m6a_sites_writes_tsv(tmp_path):
    pytest.importorskip("pysam")

    from plantm6a.analysis.site_annotation import OUTPUT_COLUMNS, annotate_m6a_sites

    genome, gtf, bed = write_fixture(tmp_path)
    output = tmp_path / "annotated.tsv"

    rows = annotate_m6a_sites(str(genome), str(gtf), str(bed), output_path=str(output), verbose=False)

    with output.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        written = list(reader)

    assert reader.fieldnames == OUTPUT_COLUMNS
    assert len(written) == len(rows)
    assert written[0]["bed_name"] == "plus_5utr"
    assert written[0]["region"] == "5utr"
