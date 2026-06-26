import pandas as pd

from plantm6a.analysis.differential import run_differential_m6a


def write_differential_fixture(tmp_path):
    annotation = tmp_path / "annotation.gtf"
    wt = tmp_path / "wt.bed6"
    mut = tmp_path / "mut.bed6"

    annotation.write_text(
        "chr1\ttest\texon\t1\t100\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";\n"
        "chr1\ttest\tCDS\t21\t80\t.\t+\t0\tgene_id \"g1\"; transcript_id \"t1\";\n"
    )
    wt.write_text(
        "chr1\t10\t11\twt_5utr\t0.20\t+\n"
        "chr1\t30\t31\tshared\t0.20\t+\n"
        "chr1\t40\t41\thyper\t0.20\t+\n"
        "chr1\t50\t51\thypo\t0.60\t+\n"
        "chr1\t90\t91\twt_3utr\t0.30\t+\n"
    )
    mut.write_text(
        "chr1\t30\t31\tshared\t0.22\t+\n"
        "chr1\t40\t41\thyper\t0.50\t+\n"
        "chr1\t50\t51\thypo\t0.20\t+\n"
        "chr1\t95\t96\tmut_3utr\t0.40\t+\n"
    )
    return annotation, wt, mut


def test_run_differential_m6a_writes_expected_outputs(tmp_path):
    annotation, wt, mut = write_differential_fixture(tmp_path)
    output_dir = tmp_path / "diff"

    result = run_differential_m6a(
        str(wt),
        str(mut),
        str(annotation),
        output_dir=str(output_dir),
        wt_name="Col",
        mut_name="C11",
        make_plots=False,
        verbose=False,
    )

    all_sites = pd.read_csv(result["all_sites"], sep="\t")
    assert set(all_sites["status"]) == {"wt_only", "mut_only", "shared", "hyper_in_mut", "hypo_in_mut"}
    assert result["site_counts"]["total_unique_positions"] == 6
    assert result["site_counts"]["wt_only"] == 2
    assert result["site_counts"]["mut_only"] == 1
    assert result["site_counts"]["shared"] == 1
    assert result["site_counts"]["hyper_in_mut"] == 1
    assert result["site_counts"]["hypo_in_mut"] == 1

    by_region = pd.read_csv(result["by_region"], sep="\t", index_col=0)
    assert by_region.loc["5utr", "wt_only"] == 1
    assert by_region.loc["cds", "shared"] == 1
    assert by_region.loc["cds", "hyper_in_mut"] == 1
    assert by_region.loc["cds", "hypo_in_mut"] == 1
    assert by_region.loc["3utr", "mut_only"] == 1

    by_gene = pd.read_csv(result["by_gene"], sep="\t")
    assert by_gene.loc[0, "gene_id"] == "g1"
    assert by_gene.loc[0, "gene_status"] == "mixed"
    assert (output_dir / "go_input" / "genes_mixed.txt").read_text() == "g1\n"
    assert (output_dir / "go_input" / "background_all_m6a.txt").read_text() == "g1\n"
    assert result["figures"] == []


def test_run_differential_m6a_can_generate_figures(tmp_path):
    annotation, wt, mut = write_differential_fixture(tmp_path)
    output_dir = tmp_path / "diff_figures"

    result = run_differential_m6a(
        str(wt),
        str(mut),
        str(annotation),
        output_dir=str(output_dir),
        wt_name="Col",
        mut_name="C11",
        make_plots=True,
        verbose=False,
    )

    assert len(result["figures"]) >= 15
    assert (output_dir / "figures" / "diff_m6a_volcano.svg").exists()
    assert (output_dir / "figures" / "diff_m6a_pie.png").exists()
