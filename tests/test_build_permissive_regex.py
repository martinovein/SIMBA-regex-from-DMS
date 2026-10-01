import hashlib
import re
import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from scripts.build_permissive_regex import AA_ORDER, AnalysisError, class_for_allowed, decide_substitution, main


def decision(effects, flags=None, wt=None, cutoff=0.5, candidate="A", wt_residue="G", **kwargs):
    flags = flags or {rep: True for rep in ("rep1", "rep2", "rep3")}
    wt = wt or {rep: 1.0 for rep in ("rep1", "rep2", "rep3")}
    return decide_substitution(
        candidate_residue=candidate,
        wt_residue=wt_residue,
        replicate_effects=effects,
        valid_flags=flags,
        wt_values=wt,
        loss_cutoff=cutoff,
        **kwargs,
    )


def test_two_valid_negative_replicates_loss_above_cutoff_excluded():
    d = decision(
        {"rep1": -0.7, "rep2": -0.8, "rep3": 0.3},
        {"rep1": True, "rep2": True, "rep3": False},
        cutoff=0.7,
    )
    assert not d.allowed
    assert d.valid_count == 2


def test_three_valid_negative_replicates_loss_below_cutoff_retained():
    d = decision({"rep1": -0.2, "rep2": -0.3, "rep3": -0.4}, cutoff=0.5)
    assert d.allowed
    assert d.all_negative
    assert not d.loss_cutoff_passed


def test_mixed_negative_and_positive_replicate_effects_retained():
    d = decision({"rep1": -1.0, "rep2": 0.1, "rep3": -1.0}, cutoff=0.5)
    assert d.allowed
    assert not d.all_negative
    assert d.reason == "retained: one or more valid paired replicate effects are zero or positive"


def test_only_one_valid_paired_replicate_retained():
    d = decision(
        {"rep1": -2.0, "rep2": -2.0, "rep3": -2.0},
        {"rep1": True, "rep2": False, "rep3": False},
        cutoff=0.5,
    )
    assert d.allowed
    assert d.valid_count == 1
    assert d.unassessable


def test_pooled_only_data_retained():
    d = decision(
        {"rep1": None, "rep2": None, "rep3": None},
        {"rep1": False, "rep2": False, "rep3": False},
        pooled_only=True,
    )
    assert d.allowed
    assert "pooled-only" in d.reason


def test_missing_mutant_retained():
    d = decision(
        {"rep1": None, "rep2": None, "rep3": None},
        {"rep1": False, "rep2": False, "rep3": False},
        mutation_missing=True,
    )
    assert d.allowed
    assert "missing" in d.reason


def test_missing_wt_denominator_retained_and_flagged():
    d = decision({"rep1": -1.0, "rep2": -1.0, "rep3": -1.0}, wt={"rep1": 1.0, "rep2": None, "rep3": 1.0})
    assert d.allowed
    assert d.unassessable
    assert "denominator missing" in d.reason


def test_nonpositive_wt_denominator_retained_and_flagged():
    d = decision({"rep1": -1.0, "rep2": -1.0, "rep3": -1.0}, wt={"rep1": 0.0, "rep2": 0.0, "rep3": 0.0})
    assert d.allowed
    assert d.unassessable
    assert "nonpositive" in d.reason


def test_value_exactly_at_unrounded_cutoff_excluded():
    cutoff = 0.851234567
    d = decision({"rep1": -cutoff, "rep2": -cutoff, "rep3": -cutoff}, cutoff=cutoff)
    assert not d.allowed
    assert d.loss_cutoff_passed


def test_wt_residue_always_retained():
    d = decision({"rep1": -10.0, "rep2": -10.0, "rep3": -10.0}, candidate="G", wt_residue="G", cutoff=0.1)
    assert d.allowed
    assert "WT residue" in d.reason


def test_regex_matches_parental_wt_motif():
    allowed_by_pos = ["A", AA_ORDER, "ILM", "S", "ST", "V", AA_ORDER, "Y"]
    regex = "^" + "".join(class_for_allowed(x) for x in allowed_by_pos) + "$"
    assert re.fullmatch(regex, "AILSTVCY")


def regex_position_count(regex):
    body = regex.strip("^").strip("$")
    count = 0
    i = 0
    while i < len(body):
        if body[i] == "[":
            j = body.index("]", i)
            i = j + 1
        else:
            i += 1
        count += 1
    return count


def run_analysis(tmp_path, *extra):
    out = tmp_path / ("out_" + "_".join(str(x).strip("-") for x in extra if str(x).startswith("--")) or "full")
    args = [
        "simba_t0_t20.xlsx",
        "--dms-only-workbook", "simba_t0_t20_DMS_only.xlsx",
        "--independent-benchmark", "S2-IE1-independant-data.xlsx",
        "--output-dir", str(out),
        *extra,
    ]
    assert main(args) == 0
    return {
        "out": out,
        "regex": pd.read_csv(out / "regex_summary.csv"),
        "decisions": pd.read_csv(out / "residue_decisions.csv"),
        "positions": pd.read_csv(out / "position_summary.csv"),
        "independent": pd.read_csv(out / "S2_IE1_independent_validation.csv"),
        "report": (out / "analysis_report.md").read_text(),
    }


def test_default_execution_analyzes_positions_1_to_18(tmp_path):
    result = run_analysis(tmp_path)
    regex = result["regex"]
    assert set(regex["analysis_positions"]) == {"1-18"}
    assert set(regex["analysis_length"]) == {18}
    assert result["decisions"]["construct_position"].min() == 1
    assert result["decisions"]["construct_position"].max() == 18


def test_full_regexes_are_18_positions_and_match_complete_wt(tmp_path):
    result = run_analysis(tmp_path)
    for row in result["regex"].itertuples(index=False):
        assert regex_position_count(row.exact_regex) == 18
        assert len(row.WT_motif) == 18
        assert re.fullmatch(row.exact_regex, row.WT_motif)


def test_full_residue_decision_row_counts_and_all_amino_acids(tmp_path):
    result = run_analysis(tmp_path)
    decisions = result["decisions"]
    assert len(decisions) == 2160
    counts = decisions.groupby(["domain", "parent"]).size()
    assert set(counts) == {360}
    aa_counts = decisions.groupby(["domain", "parent", "construct_position"])["candidate_residue"].nunique()
    assert len(aa_counts) == 6 * 18
    assert set(aa_counts) == {20}


def test_no_coverage_positions_retain_all_20_and_are_flagged(tmp_path):
    result = run_analysis(tmp_path)
    positions = result["positions"]
    no_cov = positions[positions["DMS_coverage_status"] == "no_dms_coverage"]
    assert not no_cov.empty
    assert set(no_cov["allowed_residues"]) == {AA_ORDER}
    assert set(no_cov["number_allowed"]) == {20}
    assert set(no_cov["number_excluded_among_19_non_WT"]) == {0}
    assert set(no_cov["has_supported_constraint"]) == {False}


def test_no_unassessable_or_wt_residue_is_excluded(tmp_path):
    result = run_analysis(tmp_path)
    decisions = result["decisions"]
    assert decisions[(decisions["unassessable"] == True) & (decisions["decision"] == "excluded")].empty
    assert decisions[(decisions["candidate_residue"] == decisions["WT_residue"]) & (decisions["decision"] != "allowed")].empty


def test_explicit_start_end_can_generate_smaller_window(tmp_path):
    result = run_analysis(tmp_path, "--motif-start", "6", "--motif-end", "13")
    assert set(result["regex"]["analysis_positions"]) == {"6-13"}
    assert set(result["regex"]["analysis_length"]) == {8}
    assert len(result["decisions"]) == 960
    for row in result["regex"].itertuples(index=False):
        assert regex_position_count(row.exact_regex) == 8


def test_positions_6_to_13_decisions_match_smaller_window(tmp_path):
    full = run_analysis(tmp_path / "full")
    small = run_analysis(tmp_path / "small", "--motif-start", "6", "--motif-end", "13")
    key = ["domain", "parent", "construct_position", "candidate_residue"]
    cols = key + [
        "replicate_DMS_effect_rep1", "replicate_DMS_effect_rep2", "replicate_DMS_effect_rep3",
        "replicate_DMS_pair_valid_rep1", "replicate_DMS_pair_valid_rep2", "replicate_DMS_pair_valid_rep3",
        "valid_paired_replicate_count", "mean_paired_DMS_effect", "mean_paired_WT_denominator",
        "L_a_j", "all_valid_effects_negative", "sufficient_replicates", "loss_cutoff_passed",
        "decision", "decision_reason", "unassessable",
    ]
    f = full["decisions"][full["decisions"]["construct_position"].between(6, 13)][cols].sort_values(key).reset_index(drop=True)
    s = small["decisions"][cols].sort_values(key).reset_index(drop=True)
    pd.testing.assert_frame_equal(f, s, check_dtype=False, check_exact=False, atol=1e-12, rtol=1e-12)


def test_dms_only_cross_check_and_independent_position_5_inside_full_window(tmp_path):
    result = run_analysis(tmp_path)
    assert "DMS-only cross-check matched" in result["report"]
    pos5 = result["independent"][result["independent"]["construct_position"] == 5]
    assert not pos5.empty
    assert set(pos5["inside_regex_window"]) == {True}


def test_s2_ie1_position_8_spot_checks_preserved(tmp_path):
    result = run_analysis(tmp_path)
    dec = result["decisions"]
    rows = dec[(dec["domain"] == "S2") & (dec["parent"] == "IE1") & (dec["construct_position"] == 8)]
    assert rows[rows["candidate_residue"] == "L"].iloc[0]["decision"] == "allowed"
    assert rows[rows["candidate_residue"] == "I"].iloc[0]["decision"] == "allowed"
    assert rows[rows["candidate_residue"] == "M"].iloc[0]["decision"] == "allowed"


def test_full_sequence_rejects_explicit_range(tmp_path):
    with pytest.raises(AnalysisError, match="full-sequence cannot be combined"):
        main([
            "simba_t0_t20.xlsx",
            "--output-dir", str(tmp_path / "bad"),
            "--full-sequence",
            "--motif-start", "6",
            "--motif-end", "13",
        ])


def test_independent_activities_do_not_influence_regex_generation(tmp_path):
    original = pd.read_excel("S2-IE1-independant-data.xlsx")
    modified = original.copy()
    modified["score"] = list(reversed(modified["score"].tolist()))
    modified_path = tmp_path / "modified_independent.xlsx"
    modified.to_excel(modified_path, index=False)

    files = [
        "regex_summary.csv",
        "residue_decisions.csv",
        "position_summary.csv",
        "permissive_regexes.txt",
        "permissive_regexes.json",
    ]
    out_a = tmp_path / "orig_ind"
    out_b = tmp_path / "mod_ind"
    common = [
        "simba_t0_t20.xlsx",
        "--dms-only-workbook", "simba_t0_t20_DMS_only.xlsx",
        "--output-dir",
    ]
    assert main(common + [str(out_a), "--independent-benchmark", "S2-IE1-independant-data.xlsx"]) == 0
    assert main(common + [str(out_b), "--independent-benchmark", str(modified_path)]) == 0

    for name in files:
        assert hashlib.sha256((out_a / name).read_bytes()).hexdigest() == hashlib.sha256((out_b / name).read_bytes()).hexdigest()
