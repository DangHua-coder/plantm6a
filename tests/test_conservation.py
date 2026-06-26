from plantm6a.analysis.conservation import conserved_sites
from plantm6a.analysis.conservation.conserved_sites import (
    TranscriptModel,
    build_position_index,
    find_conserved_sites,
    genomic_to_transcript_pos,
    load_m6a_sites,
    map_m6a_to_transcripts,
    parse_annotation,
)
from plantm6a.analysis.conservation.extract_orthologs import extract_pairwise_orthologs
from plantm6a.analysis.conservation.summarize import summarize_conserved_m6a
from plantm6a.analysis.conservation.translate import translate, translate_cds


def test_translate_stops_at_stop_codon():
    assert translate("ATGGCTTAAATG") == "MA"


def test_translate_cds_writes_protein_fasta(tmp_path):
    cds = tmp_path / "cds.fa"
    protein = tmp_path / "protein.fa"

    cds.write_text(">gene1\nATGGCTTAA\n>empty\nTAA\n")

    written, skipped = translate_cds(str(cds), str(protein), verbose=False)

    assert written == 1
    assert skipped == 1
    assert protein.read_text() == ">gene1\nMA\n"


def test_extract_pairwise_orthologs_writes_deduplicated_tables(tmp_path):
    orthogroups = tmp_path / "Orthogroups.tsv"
    orthologues_dir = tmp_path / "Orthologues"
    pair_dir = orthologues_dir / "Orthologues_Ara"
    output_dir = tmp_path / "pairs"
    pair_dir.mkdir(parents=True)

    orthogroups.write_text(
        "Orthogroup\tAra\tRice\tMaize\n"
        "OG0001\tAra1\tRice1\tMaize1\n"
        "OG0002\tAra2, Ara3\tRice2\t\n"
    )
    (pair_dir / "Ara__v__Rice.tsv").write_text(
        "Orthogroup\tAra\tRice\n"
        "OG0001\tAra1\tRice1\n"
        "OG0003\tAraX\tRiceX, RiceY\n"
    )

    stats = extract_pairwise_orthologs(
        str(orthogroups), str(orthologues_dir), str(output_dir), verbose=False
    )

    ara_rice = output_dir / "Ara_vs_Rice.tsv"
    assert ara_rice.exists()
    assert ara_rice.read_text().splitlines() == [
        "OG_ID\tAra_gene\tRice_gene\ttype",
        "OG0001\tAra1\tRice1\tortholog",
        "OG0003\tAraX\tRiceX\tortholog",
        "OG0003\tAraX\tRiceY\tortholog",
        "OG0002\tAra2\tRice2\tinparalog_og",
        "OG0002\tAra3\tRice2\tinparalog_og",
    ]
    assert stats[("Ara", "Rice")] == {"orthologs": 3, "inparalogs": 2, "total": 5}
    assert (output_dir / "summary_statistics.tsv").exists()


def test_parse_annotation_maps_plus_minus_and_noncoding_regions(tmp_path):
    annotation = tmp_path / "annotation.gtf"
    annotation.write_text(
        "chr1\ttest\texon\t1\t10\t.\t+\t.\tgene_id \"g_plus\"; gene_name \"PLUS\"; transcript_id \"t_plus\";\n"
        "chr1\ttest\texon\t21\t30\t.\t+\t.\tgene_id \"g_plus\"; gene_name \"PLUS\"; transcript_id \"t_plus\";\n"
        "chr1\ttest\tCDS\t4\t10\t.\t+\t0\tgene_id \"g_plus\"; gene_name \"PLUS\"; transcript_id \"t_plus\";\n"
        "chr1\ttest\tCDS\t21\t25\t.\t+\t0\tgene_id \"g_plus\"; gene_name \"PLUS\"; transcript_id \"t_plus\";\n"
        "chr1\ttest\texon\t41\t50\t.\t-\t.\tgene_id \"g_minus\"; gene_name \"MINUS\"; transcript_id \"t_minus\";\n"
        "chr1\ttest\tCDS\t44\t47\t.\t-\t0\tgene_id \"g_minus\"; gene_name \"MINUS\"; transcript_id \"t_minus\";\n"
        "chr1\ttest\texon\t61\t70\t.\t+\t.\tgene_id \"g_nc\"; gene_name \"NC\"; transcript_id \"t_nc\";\n"
    )

    models = parse_annotation(str(annotation), fmt="gtf")

    assert models["t_plus"].gene_name == "PLUS"
    assert genomic_to_transcript_pos(2, models["t_plus"]) == (2, "5utr")
    assert genomic_to_transcript_pos(5, models["t_plus"]) == (5, "cds")
    assert genomic_to_transcript_pos(25, models["t_plus"]) == (15, "3utr")
    assert genomic_to_transcript_pos(49, models["t_minus"]) == (0, "5utr")
    assert genomic_to_transcript_pos(45, models["t_minus"]) == (4, "cds")
    assert genomic_to_transcript_pos(40, models["t_minus"]) == (9, "3utr")
    assert genomic_to_transcript_pos(60, models["t_nc"]) == (0, "noncoding")


def test_map_m6a_sites_to_transcripts(tmp_path):
    bed = tmp_path / "sites.bed"
    bed.write_text("chr1\t5\t6\tsite1\t0.5\t+\nchr1\t25\t26\tsite2\t0.8\t+\nchr2\t1\t2\toff\t1\t+\n")
    model = TranscriptModel("tx1", "chr1", "+", "gene1", "GENE1")
    model.exons = [(0, 10), (20, 30)]
    model.cds_intervals = [(3, 10), (20, 25)]
    model.finalize()

    raw_sites = load_m6a_sites(str(bed))
    mapped = map_m6a_to_transcripts(raw_sites, build_position_index({"tx1": model}), {"tx1": model})

    assert raw_sites["chr1"] == [(5, "+", 0.5), (25, "+", 0.8)]
    assert mapped["tx1"] == [(5, "cds", 0.5), (15, "3utr", 0.8)]


def test_find_conserved_sites_uses_alignment_columns(monkeypatch):
    def fake_sw_align(seq_a, seq_b):
        assert seq_a == "AAAAAA"
        assert seq_b == "AAAAAA"
        return "AAAAAA-", "AAA-AAA", 12, 0, 0

    monkeypatch.setattr(conserved_sites, "sw_align", fake_sw_align)
    conserved = find_conserved_sites(
        "txA",
        "txB",
        "AAAAAA",
        "AAAAAA",
        [(3, "cds", 0.7)],
        [(3, "cds", 0.9)],
        tolerance=1,
    )

    assert len(conserved) == 1
    assert conserved[0]["tidA"] == "txA"
    assert conserved[0]["tidB"] == "txB"
    assert conserved[0]["col_diff"] == 1
    assert conserved[0]["scoreB"] == 0.9


def test_summarize_conserved_m6a_writes_expected_outputs(tmp_path):
    conserved_dir = tmp_path / "conserved"
    output_dir = tmp_path / "summary"
    conserved_dir.mkdir()
    (conserved_dir / "Ara_vs_Rice_conserved_m6A.tsv").write_text(
        "OG_ID\tAra_transcript\tRice_transcript\tAra_tx_pos\tRice_tx_pos\t"
        "Ara_region\tRice_region\tAra_m6a_score\tRice_m6a_score\t"
        "Ara_aln_col\tRice_aln_col\tcol_diff\taln_score\taln_len\tAra_context\tRice_context\n"
        "OG0001\ttA\ttR\t5\t6\tcds\tcds\t0.7\t0.8\t10\t11\t1\t12.0\t20\tAAAA\tAAAA\n"
        "OG0002\ttB\ttS\t7\t8\t3utr\t3utr\t0.4\t0.5\t14\t14\t0\t10.0\t18\tCCCC\tCCCC\n"
    )
    (conserved_dir / "Ara_vs_Maize_conserved_m6A.tsv").write_text(
        "OG_ID\tAra_transcript\tMaize_transcript\tAra_tx_pos\tMaize_tx_pos\t"
        "Ara_region\tMaize_region\tAra_m6a_score\tMaize_m6a_score\t"
        "Ara_aln_col\tMaize_aln_col\tcol_diff\taln_score\taln_len\tAra_context\tMaize_context\n"
        "OG0001\ttA\ttM\t5\t5\tcds\tcds\t0.7\t0.9\t10\t10\t0\t13.0\t20\tAAAA\tAAAA\n"
    )

    result = summarize_conserved_m6a(str(conserved_dir), str(output_dir), verbose=False)

    assert result["summary"] == [("Ara", "Maize", 1), ("Ara", "Rice", 2)]
    assert result["region_distribution"] == {"cds-cds": 2, "3utr-3utr": 1}
    assert set(result["multi_conserved_ogs"]) == {"OG0001"}
    assert (output_dir / "pair_summary.tsv").read_text().splitlines() == [
        "SpeciesA\tSpeciesB\tConserved_m6A_sites",
        "Ara\tMaize\t1",
        "Ara\tRice\t2",
    ]
    assert (output_dir / "region_distribution.tsv").exists()
    assert (output_dir / "pair_region_breakdown.tsv").exists()
    assert (output_dir / "multi_species_conserved_OGs.tsv").exists()
