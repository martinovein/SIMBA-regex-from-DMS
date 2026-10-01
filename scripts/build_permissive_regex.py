#!/usr/bin/env python3
"""Build permissive amino-acid regexes from CLN2-only SIMBA DMS workbooks."""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd


AA_ORDER = "ACDEFGHIKLMNPQRSTVWY"
REQUESTED_PAIRS = (
    ("S2", "V"),
    ("S2", "IE1"),
    ("S3", "IL20RA"),
    ("S3", "IE1"),
    ("S4", "IE1"),
    ("S4", "V"),
)
REPLICATES = ("rep1", "rep2", "rep3")
CALIBRATION_DOMAIN = "S2"
CALIBRATION_PARENT = "IE1"
INDEPENDENT_WT_WINDOW = "TIPLSSVIV"
INDEPENDENT_WINDOW_START = 5


class AnalysisError(RuntimeError):
    """Raised when workbook contents fail required validation."""


@dataclass(frozen=True)
class Decision:
    allowed: bool
    reason: str
    valid_ids: tuple[str, ...]
    valid_count: int
    mean_dms: float | None
    mean_wt: float | None
    loss: float | None
    all_negative: bool
    sufficient_replicates: bool
    loss_cutoff_passed: bool
    unassessable: bool


def is_missing(value: Any) -> bool:
    return value is None or (isinstance(value, float) and math.isnan(value)) or pd.isna(value)


def as_bool(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if is_missing(value):
        return False
    if isinstance(value, str):
        return value.strip().lower() in {"true", "1", "yes", "y"}
    return bool(value)


def finite_or_none(value: Any) -> float | None:
    if is_missing(value):
        return None
    try:
        f = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(f):
        return None
    return f


def decide_substitution(
    *,
    candidate_residue: str,
    wt_residue: str,
    replicate_effects: dict[str, float | None],
    valid_flags: dict[str, bool],
    wt_values: dict[str, float | None],
    loss_cutoff: float,
    mutation_missing: bool = False,
    pooled_only: bool = False,
) -> Decision:
    """Apply the established permissive exclusion rule to one candidate substitution."""
    if candidate_residue == wt_residue:
        return Decision(True, "WT residue always retained", tuple(), 0, 0.0, None, 0.0, False, False, False, False)

    valid_ids = tuple(
        rep
        for rep in REPLICATES
        if valid_flags.get(rep, False)
        and finite_or_none(replicate_effects.get(rep)) is not None
    )
    valid_count = len(valid_ids)
    effects = [float(replicate_effects[rep]) for rep in valid_ids]
    wt_subset = [finite_or_none(wt_values.get(rep)) for rep in valid_ids]

    if mutation_missing:
        return Decision(True, "retained: mutant observation missing", valid_ids, valid_count, None, None, None, False, False, False, True)
    if pooled_only:
        return Decision(True, "retained: pooled-only observation is unassessable", valid_ids, valid_count, None, None, None, False, False, False, True)
    if valid_count < 2:
        mean_dms = sum(effects) / valid_count if valid_count else None
        return Decision(True, "retained: fewer than two valid paired replicates", valid_ids, valid_count, mean_dms, None, None, all(x < 0 for x in effects) if effects else False, False, False, True)
    if any(v is None for v in wt_subset):
        return Decision(True, "retained: paired WT denominator missing", valid_ids, valid_count, sum(effects) / valid_count, None, None, all(x < 0 for x in effects), True, False, True)

    wt_numeric = [float(v) for v in wt_subset if v is not None]
    mean_wt = sum(wt_numeric) / len(wt_numeric)
    mean_dms = sum(effects) / len(effects)
    if mean_wt <= 0:
        return Decision(True, "retained: paired WT denominator nonpositive", valid_ids, valid_count, mean_dms, mean_wt, None, all(x < 0 for x in effects), True, False, True)

    loss = -mean_dms / mean_wt
    all_negative = all(x < 0 for x in effects)
    loss_pass = loss >= loss_cutoff
    if not all_negative:
        reason = "retained: one or more valid paired replicate effects are zero or positive"
    elif not loss_pass:
        reason = "retained: parent-relative loss below calibrated cutoff"
    else:
        reason = "excluded: sufficient paired replicates, all negative, and loss reaches cutoff"
    return Decision(not (all_negative and loss_pass), reason, valid_ids, valid_count, mean_dms, mean_wt, loss, all_negative, True, loss_pass, False)


def require_columns(df: pd.DataFrame, sheet: str, columns: list[str]) -> None:
    missing = [c for c in columns if c not in df.columns]
    if missing:
        raise AnalysisError(f"{sheet} is missing required columns: {missing}")


def read_required_sheets(path: Path) -> dict[str, pd.DataFrame]:
    required = ["Run Parameters", "DMS Variant Summary", "Replicate SIMBA Summary", "DMS Matrix Long"]
    try:
        xl = pd.ExcelFile(path)
    except Exception as exc:
        raise AnalysisError(f"Could not open {path}: {exc}") from exc
    missing = [s for s in required if s not in xl.sheet_names]
    if missing:
        raise AnalysisError(f"{path} is missing required worksheets: {missing}")
    return {sheet: pd.read_excel(path, sheet_name=sheet) for sheet in required}


def validate_run_parameters(params: pd.DataFrame) -> dict[str, Any]:
    require_columns(params, "Run Parameters", ["parameter", "value"])
    values = dict(zip(params["parameter"].astype(str), params["value"]))
    slope_ref = str(values.get("slope_reference_strains", "")).strip()
    if slope_ref != "CLN2":
        raise AnalysisError(f"Expected slope_reference_strains exactly CLN2, found {slope_ref!r}")
    command_line = str(values.get("command_line", ""))
    if "--slope-reference-strains CLN2" not in command_line:
        raise AnalysisError("Run Parameters command_line does not show original CLN2-only normalization")
    if any(term in slope_ref for term in ("IKK", "TBK")) or "--slope-reference-strains CLN2,IKK,TBK" in command_line:
        raise AnalysisError("Workbook appears to use later multi-reference normalization, not CLN2-only")
    return values


def get_parent_wt_sequences(rs: pd.DataFrame) -> dict[str, str]:
    wt_rows = rs[rs["variant_id"].astype(str).str.endswith("_WT", na=False)]
    seqs: dict[str, str] = {}
    for parent, group in wt_rows.groupby("parent"):
        unique = sorted(set(group["protein_sequence"].dropna().astype(str)))
        if len(unique) != 1:
            raise AnalysisError(f"Expected one WT sequence for parent {parent}, found {unique}")
        seqs[str(parent)] = unique[0]
    return seqs


def calculate_calibration(rs: pd.DataFrame, tau: float) -> tuple[float, float, list[float]]:
    row = rs[
        (rs["parent"] == CALIBRATION_PARENT)
        & (rs["tested_strain"] == CALIBRATION_DOMAIN)
        & (rs["variant_id"].astype(str) == f"{CALIBRATION_PARENT}_WT")
    ]
    if len(row) != 1:
        raise AnalysisError("Could not identify unique S2-IE1 WT row in Replicate SIMBA Summary")
    vals = [finite_or_none(row.iloc[0][f"replicate_SIMBA_slope_delta_{rep}"]) for rep in REPLICATES]
    if any(v is None for v in vals):
        raise AnalysisError(f"S2-IE1 WT calibration replicate values are incomplete: {vals}")
    wt_mean = sum(float(v) for v in vals if v is not None) / len(vals)
    if wt_mean <= 0:
        raise AnalysisError(f"S2-IE1 WT calibration denominator must be positive, found {wt_mean}")
    return -tau / wt_mean, wt_mean, [float(v) for v in vals if v is not None]


def build_decision_rows(
    dms: pd.DataFrame,
    matrix_long: pd.DataFrame,
    wt_sequences: dict[str, str],
    loss_cutoff: float,
    analysis_ranges: dict[str, tuple[int, int]],
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    indexed = {
        (str(r.parent), str(r.tested_strain), int(r.position), str(r.mutant_residue)): r
        for r in dms.itertuples(index=False)
        if str(getattr(r, "mutant_residue")) in AA_ORDER
    }
    designed_counts = (
        dms[dms["mutant_residue"].astype(str).isin(list(AA_ORDER))]
        .groupby(["parent", "tested_strain", "position"])["mutant_residue"]
        .nunique()
        .to_dict()
    )
    matrix_best = matrix_long[matrix_long["analysis_type"] == "best_available"] if "analysis_type" in matrix_long.columns else matrix_long.iloc[0:0]
    matrix_index = {
        (str(r.parent), str(r.tested_strain), int(r.position), str(r.amino_acid)): r
        for r in matrix_best.itertuples(index=False)
        if str(getattr(r, "amino_acid", "")) in AA_ORDER and not is_missing(getattr(r, "position", None))
    }

    for domain, parent in REQUESTED_PAIRS:
        if parent not in wt_sequences:
            raise AnalysisError(f"Missing WT sequence for requested parent {parent}")
        wt_seq = wt_sequences[parent]
        start, end = analysis_ranges[parent]
        if len(wt_seq) < end:
            raise AnalysisError(f"WT sequence for parent {parent} is shorter than analysis end {end}")
        sub = dms[(dms["parent"] == parent) & (dms["tested_strain"] == domain)]
        if sub.empty:
            raise AnalysisError(f"Missing requested parent-domain combination {domain}-{parent}")

        for pos in range(start, end + 1):
            wt_residue = wt_seq[pos - 1]
            if wt_residue not in AA_ORDER:
                raise AnalysisError(f"{parent} WT residue at position {pos} is not standard amino acid: {wt_residue}")
            annotated = set(sub.loc[sub["position"] == pos, "wt_residue"].dropna().astype(str))
            if annotated and annotated != {wt_residue}:
                raise AnalysisError(f"{domain}-{parent} position {pos} WT mismatch: sequence {wt_residue}, DMS annotations {sorted(annotated)}")
            position_designed_count = int(designed_counts.get((parent, domain, pos), 0))
            no_position_coverage = position_designed_count == 0
            for aa in AA_ORDER:
                source = "DMS Variant Summary"
                mutation_missing = False
                pooled_only = False
                row = indexed.get((parent, domain, pos, aa))
                matrix_row = matrix_index.get((parent, domain, pos, aa))
                matrix_missing_reason = None if matrix_row is None else getattr(matrix_row, "missing_data_reason", None)
                matrix_best_source = None if matrix_row is None else getattr(matrix_row, "best_available_source", None)
                matrix_failure_reason = None if matrix_row is None else getattr(matrix_row, "best_available_failure_reason", None)
                variant_design_status = "WT residue" if aa == wt_residue else "designed standard-AA mutant"
                if aa == wt_residue:
                    effects = {rep: 0.0 for rep in REPLICATES}
                    flags = {rep: False for rep in REPLICATES}
                    wt_values = {rep: None for rep in REPLICATES}
                    mutation = f"{wt_residue}{pos}{aa}"
                elif row is None:
                    effects = {rep: None for rep in REPLICATES}
                    flags = {rep: False for rep in REPLICATES}
                    wt_values = {rep: None for rep in REPLICATES}
                    mutation = f"{wt_residue}{pos}{aa}"
                    mutation_missing = True
                    if no_position_coverage:
                        variant_design_status = "not present in designed standard-AA library"
                        source = "constructed no-DMS-coverage candidate"
                    else:
                        variant_design_status = "not present in DMS Variant Summary or score unavailable"
                        source = "constructed missing/unavailable candidate"
                else:
                    effects = {rep: finite_or_none(getattr(row, f"replicate_DMS_effect_{rep}")) for rep in REPLICATES}
                    flags = {rep: as_bool(getattr(row, f"replicate_DMS_pair_valid_{rep}")) for rep in REPLICATES}
                    wt_values = {rep: finite_or_none(getattr(row, f"{rep}_wt_SIMBA_slope_delta")) for rep in REPLICATES}
                    mutation = str(getattr(row, "mutation"))
                    best_source = str(getattr(row, "best_available_source", ""))
                    pooled_only = best_source == "pooled_fallback" and sum(flags.values()) < 2
                decision = decide_substitution(
                    candidate_residue=aa,
                    wt_residue=wt_residue,
                    replicate_effects=effects,
                    valid_flags=flags,
                    wt_values=wt_values,
                    loss_cutoff=loss_cutoff,
                    mutation_missing=mutation_missing,
                    pooled_only=pooled_only,
                )
                reason = decision.reason
                rows.append(
                    {
                        "domain": domain,
                        "parent": parent,
                        "construct_position": pos,
                        "motif_relative_position": pos - start + 1,
                        "WT_residue": wt_residue,
                        "candidate_residue": aa,
                        "mutation_identifier": mutation,
                        "replicate_DMS_effect_rep1": effects["rep1"],
                        "replicate_DMS_effect_rep2": effects["rep2"],
                        "replicate_DMS_effect_rep3": effects["rep3"],
                        "replicate_DMS_pair_valid_rep1": flags["rep1"],
                        "replicate_DMS_pair_valid_rep2": flags["rep2"],
                        "replicate_DMS_pair_valid_rep3": flags["rep3"],
                        "valid_paired_replicate_ids": ";".join(decision.valid_ids),
                        "valid_paired_replicate_count": decision.valid_count,
                        "mean_paired_DMS_effect": decision.mean_dms,
                        "paired_WT_value_rep1": wt_values["rep1"],
                        "paired_WT_value_rep2": wt_values["rep2"],
                        "paired_WT_value_rep3": wt_values["rep3"],
                        "mean_paired_WT_denominator": decision.mean_wt,
                        "L_a_j": decision.loss,
                        "all_valid_effects_negative": decision.all_negative,
                        "sufficient_replicates": decision.sufficient_replicates,
                        "loss_cutoff_passed": decision.loss_cutoff_passed,
                        "decision": "allowed" if decision.allowed else "excluded",
                        "decision_reason": reason,
                        "unassessable": decision.unassessable,
                        "variant_design_status": variant_design_status,
                        "position_has_dms_coverage": not no_position_coverage,
                        "matrix_missing_data_reason": matrix_missing_reason,
                        "matrix_best_available_source": matrix_best_source,
                        "matrix_best_available_failure_reason": matrix_failure_reason,
                        "data_provenance_source": source,
                    }
                )
    decisions = pd.DataFrame(rows)
    bad_exclusions = decisions[(decisions["decision"] == "excluded") & (decisions["unassessable"])]
    if not bad_exclusions.empty:
        raise AnalysisError("Internal error: unassessable residues were excluded")
    return decisions


def class_for_allowed(allowed: str) -> str:
    return allowed if len(allowed) == 1 else f"[{allowed}]"


def compact_piece(allowed: str) -> str:
    return "." if allowed == AA_ORDER else class_for_allowed(allowed)


def summarize_regexes(decisions: pd.DataFrame, wt_sequences: dict[str, str], tau: float, loss_cutoff: float, analysis_ranges: dict[str, tuple[int, int]]) -> tuple[pd.DataFrame, pd.DataFrame, list[dict[str, Any]]]:
    regex_rows: list[dict[str, Any]] = []
    pos_rows: list[dict[str, Any]] = []
    json_rows: list[dict[str, Any]] = []
    for domain, parent in REQUESTED_PAIRS:
        pieces: list[str] = []
        compact: list[str] = []
        start, end = analysis_ranges[parent]
        wt_motif = wt_sequences[parent][start - 1 : end]
        pair = decisions[(decisions["domain"] == domain) & (decisions["parent"] == parent)]
        for pos in range(start, end + 1):
            sub = pair[pair["construct_position"] == pos]
            wt = wt_sequences[parent][pos - 1]
            allowed = "".join(aa for aa in AA_ORDER if ((sub["candidate_residue"] == aa) & (sub["decision"] == "allowed")).any())
            excluded = "".join(aa for aa in AA_ORDER if aa != wt and ((sub["candidate_residue"] == aa) & (sub["decision"] == "excluded")).any())
            if not allowed:
                raise AnalysisError(f"{domain}-{parent} position {pos} retains no residues")
            if wt not in allowed:
                raise AnalysisError(f"{domain}-{parent} position {pos} lost WT residue {wt}")
            pieces.append(class_for_allowed(allowed))
            compact.append(compact_piece(allowed))
            non_wt = sub[sub["candidate_residue"] != wt]
            non_wt_excluded = len(excluded)
            designed_mutants = int(non_wt[non_wt["variant_design_status"] == "designed standard-AA mutant"].shape[0])
            usable_dms = int(non_wt[non_wt["sufficient_replicates"] & non_wt["mean_paired_WT_denominator"].notna() & (non_wt["mean_paired_WT_denominator"] > 0)].shape[0])
            unassessable_count = int(non_wt[non_wt["unassessable"]].shape[0])
            if designed_mutants == 0:
                coverage_status = "no_dms_coverage"
            elif usable_dms == 19:
                coverage_status = "complete_usable_dms_coverage"
            else:
                coverage_status = "incomplete_usable_dms_coverage"
            pos_rows.append(
                {
                    "domain": domain,
                    "parent": parent,
                    "construct_position": pos,
                    "motif_relative_position": pos - start + 1,
                    "WT_residue": wt,
                    "allowed_residues": allowed,
                    "excluded_residues": excluded,
                    "number_allowed": len(allowed),
                    "number_excluded_among_19_non_WT": non_wt_excluded,
                    "fraction_non_WT_substitutions_excluded": non_wt_excluded / 19,
                    "number_designed_mutants": designed_mutants,
                    "number_with_usable_DMS_data": usable_dms,
                    "number_unassessable": unassessable_count,
                    "DMS_coverage_status": coverage_status,
                    "has_supported_constraint": non_wt_excluded > 0,
                    "position_constrained": len(allowed) < 20,
                }
            )
        regex = "^" + "".join(pieces) + "$"
        compact_regex = "^" + "".join(compact) + "$"
        if re.fullmatch(regex, wt_motif) is None:
            raise AnalysisError(f"{domain}-{parent} regex does not match its WT sequence/window")
        if len(wt_motif) != end - start + 1:
            raise AnalysisError(f"{domain}-{parent} WT sequence/window length mismatch")
        total_excluded = int(pair[pair["decision"] == "excluded"].shape[0])
        unassessable = int(pair[(pair["candidate_residue"] != pair["WT_residue"]) & (pair["unassessable"])].shape[0])
        constrained_positions = sum(1 for p in pieces if p != f"[{AA_ORDER}]")
        row = {
            "domain": domain,
            "parent": parent,
            "WT_motif": wt_motif,
            "WT_sequence_or_window": wt_motif,
            "exact_regex": regex,
            "compact_regex": compact_regex,
            "motif_positions": f"{start}-{end}",
            "analysis_positions": f"{start}-{end}",
            "analysis_length": end - start + 1,
            "calibration_tau": tau,
            "exact_L_cut": loss_cutoff,
            "rounded_L_cut": round(loss_cutoff, 3),
            "number_constrained_positions": constrained_positions,
            "total_excluded_substitutions": total_excluded,
            "missing_or_unassessable_substitutions": unassessable,
        }
        regex_rows.append(row)
        json_rows.append(row | {"positions": [r for r in pos_rows if r["domain"] == domain and r["parent"] == parent]})
    return pd.DataFrame(regex_rows), pd.DataFrame(pos_rows), json_rows



def s2_ie1_matrix_median(dms: pd.DataFrame, motif_start: int, motif_end: int) -> float:
    sub = dms[
        (dms["parent"] == CALIBRATION_PARENT)
        & (dms["tested_strain"] == CALIBRATION_DOMAIN)
        & (dms["position"].between(motif_start, motif_end))
        & (~dms["contains_stop"])
    ]
    vals = pd.to_numeric(sub["mean_replicate_DMS_effect"], errors="coerce").dropna()
    if vals.empty:
        raise AnalysisError("No S2-IE1 mean replicate DMS values available for calibration QC")
    return float(vals.median())


def compare_dms_only(full_matrix: pd.DataFrame, dms_only_path: Path, motif_start: int, motif_end: int) -> dict[str, Any]:
    xls = pd.ExcelFile(dms_only_path)
    if "DMS Long Mean Rep" not in xls.sheet_names:
        raise AnalysisError(f"{dms_only_path} lacks DMS Long Mean Rep")
    other = pd.read_excel(dms_only_path, sheet_name="DMS Long Mean Rep")
    full = full_matrix[
        (full_matrix["parent"] == CALIBRATION_PARENT)
        & (full_matrix["tested_strain"] == CALIBRATION_DOMAIN)
        & (full_matrix["analysis_type"] == "mean_replicate")
        & (full_matrix["position"].between(motif_start, motif_end))
    ]
    other = other[
        (other["parent"] == CALIBRATION_PARENT)
        & (other["tested_strain"] == CALIBRATION_DOMAIN)
        & (other["analysis_type"] == "mean_replicate")
        & (other["position"].between(motif_start, motif_end))
    ]
    merged = full.merge(other, on=["parent", "tested_strain", "analysis_type", "position", "amino_acid"], suffixes=("_full", "_dms_only"))
    if merged.empty:
        raise AnalysisError("Could not match S2-IE1 full and DMS-only matrix values")
    diffs = (pd.to_numeric(merged["DMS_effect_full"], errors="coerce") - pd.to_numeric(merged["DMS_effect_dms_only"], errors="coerce")).abs()
    max_diff = float(diffs.dropna().max()) if not diffs.dropna().empty else 0.0
    if max_diff > 1e-12:
        raise AnalysisError(f"S2-IE1 full and DMS-only mean replicate matrices differ; max abs diff {max_diff}")
    return {"matched_cells": int(len(merged)), "max_abs_diff": max_diff}


def independent_validation(independent_path: Path, decisions: pd.DataFrame, analysis_ranges: dict[str, tuple[int, int]]) -> pd.DataFrame:
    ind = pd.read_excel(independent_path)
    require_columns(ind, independent_path.name, ["peptide", "score"])
    out: list[dict[str, Any]] = []
    s2_ie1 = decisions[(decisions["domain"] == "S2") & (decisions["parent"] == "IE1")]
    start, end = analysis_ranges["IE1"]
    for _, row in ind.iterrows():
        raw = str(row["peptide"])
        activity = finite_or_none(row["score"])
        peptide = raw.removeprefix("WT_")
        null_like = activity is not None and activity <= 0.25
        base = {
            "independent_peptide_sequence": raw,
            "independent_retained_activity": activity,
            "null_like_activity_le_0_25": null_like,
        }
        if peptide != INDEPENDENT_WT_WINDOW or raw.startswith("WT_"):
            if len(peptide) != len(INDEPENDENT_WT_WINDOW) or peptide == "c-_STOP":
                out.append(base | {"construct_position": None, "inside_regex_window": False, "decision_reason": "outside regex classification: non-single-mutant or control"})
                continue
            diffs = [i for i, (a, b) in enumerate(zip(peptide, INDEPENDENT_WT_WINDOW)) if a != b]
            if len(diffs) != 1:
                out.append(base | {"construct_position": None, "inside_regex_window": False, "decision_reason": "outside regex classification: not a single mutant"})
                continue
            idx = diffs[0]
            wt = INDEPENDENT_WT_WINDOW[idx]
            mut = peptide[idx]
            pos = INDEPENDENT_WINDOW_START + idx
        else:
            pos = None
            wt = None
            mut = None
        inside = pos is not None and start <= pos <= end
        record = base | {
            "construct_position": pos,
            "WT_residue": wt,
            "mutant_residue": mut,
            "inside_regex_window": inside,
        }
        if not inside:
            out.append(record | {"decision_reason": f"outside analyzed positions {start}-{end}"})
            continue
        match = s2_ie1[(s2_ie1["construct_position"] == pos) & (s2_ie1["candidate_residue"] == mut)]
        if len(match) != 1:
            raise AnalysisError(f"Could not map independent mutation {raw} to S2-IE1 decision")
        d = match.iloc[0].to_dict()
        out.append(
            record
            | {
                "replicate_DMS_effect_rep1": d["replicate_DMS_effect_rep1"],
                "replicate_DMS_effect_rep2": d["replicate_DMS_effect_rep2"],
                "replicate_DMS_effect_rep3": d["replicate_DMS_effect_rep3"],
                "paired_replicate_count": d["valid_paired_replicate_count"],
                "mean_DMS_effect": d["mean_paired_DMS_effect"],
                "paired_WT_denominator": d["mean_paired_WT_denominator"],
                "L_a_j": d["L_a_j"],
                "every_valid_effect_negative": d["all_valid_effects_negative"],
                "fixed_regex_decision": d["decision"],
                "decision_reason": d["decision_reason"],
            }
        )
    return pd.DataFrame(out)


def write_report(
    outdir: Path,
    args: argparse.Namespace,
    run_values: dict[str, Any],
    regex_summary: pd.DataFrame,
    position_summary: pd.DataFrame,
    independent: pd.DataFrame | None,
    calibration: dict[str, Any],
    warnings: list[str],
    dms_only_qc: dict[str, Any] | None,
) -> None:
    lines = [
        "# Permissive SIMBA DMS Regex Analysis",
        "",
        "## Formulas",
        "",
        f"- Calibration uses tau_IE1 = {args.tau_ie1}.",
        "- D_WT,cal is the mean of S2-IE1 WT replicate SIMBA slope-delta values.",
        "- L_cut = -tau_IE1 / D_WT,cal, using the unrounded value for comparisons.",
        "- For each mutant, mean DMS effect and paired WT denominator are computed over the same valid paired replicate subset.",
        "- L_a_j = -mean_paired_DMS_effect / mean_paired_WT_denominator.",
        "",
        "## Calibration",
        "",
        f"- S2-IE1 WT replicate values: {calibration['s2_ie1_wt_replicates']}",
        f"- D_WT,cal: {calibration['s2_ie1_wt_mean']:.12g}",
        f"- Exact L_cut: {calibration['loss_cutoff']:.12g}",
        f"- Rounded L_cut: {round(calibration['loss_cutoff'], 3):.3f}",
        f"- S2-IE1 positions 6-13 matrix median: {calibration['s2_ie1_matrix_median']:.12g}",
        "",
        "## Input Validation",
        "",
        f"- slope_reference_strains: {run_values.get('slope_reference_strains')}",
        f"- command_line: `{run_values.get('command_line')}`",
        "- Confirmed original CLN2-only normalization; no CLN2/IKK/TBK normalization was accepted.",
        "",
        "## Missing-Data Policy",
        "",
        "- Stop codons are excluded from regex construction.",
        "- Missing, pooled-only, insufficient-replicate, or invalid-denominator substitutions are retained and flagged.",
        "- WT residues are always retained.",
        "",
        "## Replicate Rules",
        "",
        "- A non-WT residue is excluded only with at least two valid paired replicates, all valid DMS effects strictly negative, and L_a_j >= exact L_cut.",
        "- Pooled-only observations are never used to exclude residues.",
        "",
        "## Regex Construction",
        "",
    ]
    for r in regex_summary.itertuples(index=False):
        lines.append(f"- {r.domain}-{r.parent}: WT {r.WT_motif}, regex `{r.exact_regex}`, compact `{r.compact_regex}`")
    lines += ["", "## Numerical QC", ""]
    if dms_only_qc:
        lines.append(f"- DMS-only cross-check matched {dms_only_qc['matched_cells']} S2-IE1 cells; max absolute difference {dms_only_qc['max_abs_diff']}.")
    else:
        lines.append("- DMS-only cross-check was not run because no DMS-only workbook was supplied.")
    lines += ["", "## Independent Benchmark Results", ""]
    if independent is None:
        lines.append("- Independent benchmark was not supplied.")
    else:
        inside = independent[independent["inside_regex_window"] == True]
        lines.append(f"- Independent rows: {len(independent)}; inside analyzed window: {len(inside)}.")
        if "fixed_regex_decision" in inside:
            counts = inside["fixed_regex_decision"].value_counts(dropna=False).to_dict()
            lines.append(f"- Fixed-rule decisions inside window: {counts}.")
    lines += ["", "## Warnings", ""]
    lines.extend([f"- {w}" for w in warnings] or ["- None."])
    lines += ["", "## Constrained-Position Summaries", ""]
    for r in position_summary[position_summary["has_supported_constraint"]].itertuples(index=False):
        lines.append(f"- {r.domain}-{r.parent} position {r.construct_position} ({r.WT_residue}): allowed {r.allowed_residues}; excluded {r.excluded_residues}; fraction excluded {r.fraction_non_WT_substitutions_excluded:.3f}")
    lines += [
        "",
        "## Reproduction",
        "",
        "```bash",
        f"PYTHONPATH=.pythondeps python3 scripts/build_permissive_regex.py {args.simba_workbook} --dms-only-workbook {args.dms_only_workbook or ''} --independent-benchmark {args.independent_benchmark or ''} --output-dir {args.output_dir} --full-sequence --tau-ie1 {args.tau_ie1}",
        "```",
    ]
    (outdir / "analysis_report.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("simba_workbook", type=Path)
    p.add_argument("--dms-only-workbook", type=Path)
    p.add_argument("--independent-benchmark", type=Path)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--motif-start", type=int, default=None, help="Optional 1-based start position for a user-selected window")
    p.add_argument("--motif-end", type=int, default=None, help="Optional 1-based end position for a user-selected window")
    p.add_argument("--full-sequence", action="store_true", help="Analyze each parental peptide's complete inferred WT sequence; this is the default when no range is supplied")
    p.add_argument("--tau-ie1", type=float, default=-0.04684)
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    explicit_range = args.motif_start is not None or args.motif_end is not None
    if args.full_sequence and explicit_range:
        raise AnalysisError("--full-sequence cannot be combined with --motif-start/--motif-end")
    if explicit_range and (args.motif_start is None or args.motif_end is None):
        raise AnalysisError("Provide both --motif-start and --motif-end, or neither for full-sequence analysis")
    if explicit_range and args.motif_end < args.motif_start:
        raise AnalysisError("motif-end must be >= motif-start")
    outdir = args.output_dir
    outdir.mkdir(parents=True, exist_ok=True)

    sheets = read_required_sheets(args.simba_workbook)
    run_values = validate_run_parameters(sheets["Run Parameters"])
    dms = sheets["DMS Variant Summary"]
    rs = sheets["Replicate SIMBA Summary"]
    matrix_long = sheets["DMS Matrix Long"]
    require_columns(dms, "DMS Variant Summary", [
        "parent", "tested_strain", "position", "wt_residue", "mutant_residue", "contains_stop",
        "mutation", "mean_replicate_DMS_effect",
        "replicate_DMS_effect_rep1", "replicate_DMS_pair_valid_rep1", "rep1_wt_SIMBA_slope_delta",
        "replicate_DMS_effect_rep2", "replicate_DMS_pair_valid_rep2", "rep2_wt_SIMBA_slope_delta",
        "replicate_DMS_effect_rep3", "replicate_DMS_pair_valid_rep3", "rep3_wt_SIMBA_slope_delta",
    ])
    require_columns(rs, "Replicate SIMBA Summary", ["protein_sequence", "parent", "variant_id", "tested_strain"] + [f"replicate_SIMBA_slope_delta_{rep}" for rep in REPLICATES])

    wt_sequences = get_parent_wt_sequences(rs)
    if explicit_range:
        analysis_ranges = {parent: (args.motif_start, args.motif_end) for _, parent in REQUESTED_PAIRS}
    else:
        analysis_ranges = {parent: (1, len(wt_sequences[parent])) for _, parent in REQUESTED_PAIRS}
        args.full_sequence = True
        args.motif_start = min(start for start, _ in analysis_ranges.values())
        args.motif_end = max(end for _, end in analysis_ranges.values())
    loss_cutoff, wt_mean, wt_reps = calculate_calibration(rs, args.tau_ie1)
    median = s2_ie1_matrix_median(dms, 6, 13)
    warnings: list[str] = []
    if abs(median - args.tau_ie1) > 5e-4:
        warnings.append(f"S2-IE1 matrix median {median:.8g} differs from tau {args.tau_ie1:.8g} by more than 5e-4")

    decisions = build_decision_rows(dms, matrix_long, wt_sequences, loss_cutoff, analysis_ranges)
    regex_summary, position_summary, json_rows = summarize_regexes(decisions, wt_sequences, args.tau_ie1, loss_cutoff, analysis_ranges)
    if decisions["candidate_residue"].eq("*").any():
        raise AnalysisError("Stop codon appeared in residue decisions")
    expected_rows = sum((analysis_ranges[parent][1] - analysis_ranges[parent][0] + 1) * len(AA_ORDER) for _, parent in REQUESTED_PAIRS)
    if len(decisions) != expected_rows:
        raise AnalysisError(f"Residue decision output has {len(decisions)} rows; expected {expected_rows}")

    dms_only_qc = None
    if args.dms_only_workbook:
        dms_only_qc = compare_dms_only(matrix_long, args.dms_only_workbook, args.motif_start, args.motif_end)

    independent = None
    if args.independent_benchmark:
        independent = independent_validation(args.independent_benchmark, decisions, analysis_ranges)
        independent.to_csv(outdir / "S2_IE1_independent_validation.csv", index=False)

    regex_summary.to_csv(outdir / "regex_summary.csv", index=False)
    decisions.to_csv(outdir / "residue_decisions.csv", index=False)
    position_summary.to_csv(outdir / "position_summary.csv", index=False)
    with (outdir / "permissive_regexes.json").open("w", encoding="utf-8") as fh:
        json.dump({"regexes": json_rows, "calibration": {"tau_ie1": args.tau_ie1, "loss_cutoff": loss_cutoff, "s2_ie1_wt_mean": wt_mean}}, fh, indent=2)
    with (outdir / "permissive_regexes.txt").open("w", encoding="utf-8") as fh:
        for r in regex_summary.itertuples(index=False):
            fh.write(f"{r.domain}-{r.parent}\tWT={r.WT_motif}\tregex={r.exact_regex}\tcompact={r.compact_regex}\n")
    write_report(
        outdir,
        args,
        run_values,
        regex_summary,
        position_summary,
        independent,
        {"loss_cutoff": loss_cutoff, "s2_ie1_wt_mean": wt_mean, "s2_ie1_wt_replicates": wt_reps, "s2_ie1_matrix_median": median},
        warnings,
        dms_only_qc,
    )
    print(f"L_cut exact={loss_cutoff:.12g} rounded={round(loss_cutoff, 3):.3f}")
    print(f"S2-IE1 matrix median positions 6-13={median:.12g}")
    for r in regex_summary.itertuples(index=False):
        print(f"{r.domain}-{r.parent}: WT {r.WT_motif} {r.exact_regex}")
    if warnings:
        print("Warnings:")
        for w in warnings:
            print(f"- {w}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except AnalysisError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2)
