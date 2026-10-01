# SIMBA DMS Permissive Regexes

This repository builds permissive amino-acid regular expressions from the original CLN2-only SIMBA DMS analysis. By default, the script now analyzes each requested parental peptide's complete inferred WT sequence. In this dataset, all requested parents are 18 amino acids long, so the principal regexes cover construct positions 1-18.

The main implementation is [`scripts/build_permissive_regex.py`](scripts/build_permissive_regex.py). Focused rule and end-to-end regression tests are in [`tests/test_build_permissive_regex.py`](tests/test_build_permissive_regex.py). Generated outputs from the current run are under [`outputs/permissive_regex/`](outputs/permissive_regex/).

## Inputs

The script expects the original SIMBA workbooks as explicit command-line paths:

- `simba_t0_t20.xlsx`: authoritative source for paired replicate-level mutant and WT scores. The script reads `Run Parameters`, `DMS Variant Summary`, `Replicate SIMBA Summary`, and `DMS Matrix Long`.
- `simba_t0_t20_DMS_only.xlsx`: optional numerical cross-check of exported DMS matrices.
- `S2-IE1-independant-data.xlsx`: optional independent S2-IE1 benchmark validation. This is used only after regex construction and never changes thresholds, residue decisions, or regexes.

The input workbooks are read-only inputs. The script writes all results to a separate output directory.

## Required Normalization Check

Before computing regexes, the script checks `Run Parameters` and aborts unless:

- `slope_reference_strains` is exactly `CLN2`;
- the recorded command line is consistent with `--slope-reference-strains CLN2`;
- the workbook does not indicate later CLN2/IKK/TBK multi-reference normalization.

This keeps the analysis tied to the original CLN2-only SIMBA normalization.

## Parent-Domain Regexes

Six regexes are generated independently:

- S2-V
- S2-IE1
- S3-IL20RA
- S3-IE1
- S4-IE1
- S4-V

Matrices from different parental peptides are never merged. Each parent-domain pair gets its own residue-level audit table and position summary.

## Analysis Window

Default execution uses full-sequence mode. The script infers the WT sequence independently for each parent from `Replicate SIMBA Summary` and analyzes positions `1-len(WT sequence)`. For the current dataset this is positions 1-18, producing 18-residue regexes.

A smaller user-selected window can still be requested with `--motif-start` and `--motif-end`. The explicit `--full-sequence` flag documents full-sequence intent, but it cannot be combined with an explicit start/end range.

The fixed S2-IE1 calibration QC still reports the matrix median for positions 6-13 because that was the original calibration motif check. That QC window does not define a core region and does not limit regex generation.

Residues are ordered with the canonical amino-acid alphabet:

```text
ACDEFGHIKLMNPQRSTVWY
```

Stop codons are excluded from regex construction.

## Calibration

The fixed S2-IE1 calibration uses:

```text
tau_IE1 = -0.04684
```

The S2-IE1 WT SIMBA denominator is the mean of the three S2-IE1 WT replicate SIMBA signals:

```text
D_WT,cal = mean(D_WT,S2-IE1,rep1, D_WT,S2-IE1,rep2, D_WT,S2-IE1,rep3)
```

The loss cutoff is then:

```text
L_cut = -tau_IE1 / D_WT,cal
```

For the current workbook, the exact value is:

```text
L_cut = 0.8507418657305194
```

The rounded reporting value is `0.851`, but comparisons use the exact unrounded value. The S2-IE1 DMS matrix median for positions 6-13 is reported as a QC check only; it is not used to tune the cutoff.

## Per-Mutant Calculation

For a candidate amino acid at a parent-domain-position, the script identifies valid paired replicates using:

- `replicate_DMS_pair_valid_rep1`
- `replicate_DMS_pair_valid_rep2`
- `replicate_DMS_pair_valid_rep3`

For the same valid replicate subset, it computes:

```text
mean_paired_DMS_effect = mean(replicate_DMS_effect_r for valid paired replicates r)
mean_paired_WT_denominator = mean(D_WT,r for the same valid paired replicates r)
L_a_j = -mean_paired_DMS_effect / mean_paired_WT_denominator
```

The WT denominator is intentionally replicate-matched to the mutant evidence. A global WT denominator is not used when only a subset of paired replicates is valid.

If the paired WT denominator is missing, zero, or negative, the substitution is marked unassessable and retained.

## Exclusion Rule

A non-WT residue is excluded only when all three conditions are true:

1. At least two valid paired replicates are available.
2. Every valid paired replicate DMS effect is strictly negative.
3. The original mean parent-relative loss reaches the exact calibrated cutoff: `L_a_j >= L_cut`.

This is a strict conjunction. A residue remains allowed if any condition fails, including insufficient paired replicates, mixed signs, zero or positive effects, loss below cutoff, missing mutant data, pooled-only data, or invalid WT denominator. Positive replicate effects are not clipped, adjusted, shrunk, rewarded, or reinterpreted inside the regex-generation script.

The WT residue is always restored and retained.

## Positions Without DMS Coverage

Positions lacking DMS coverage retain all 20 standard amino acids and are represented explicitly as `[ACDEFGHIKLMNPQRSTVWY]`. These positions are flagged in `position_summary.csv` with `DMS_coverage_status = no_dms_coverage`.

Where the workbook permits it, the residue audit distinguishes variants that were not present in the designed standard-AA library from variants present but unavailable or otherwise unassessable. No unassessable residue can be excluded.

## Regex Construction

For each analyzed position:

1. Start with all 20 standard amino acids allowed.
2. Remove only non-WT residues satisfying the full exclusion rule.
3. Restore the WT residue.
4. Sort allowed residues in canonical amino-acid order.
5. Write a literal residue if only one amino acid is allowed, otherwise write an explicit bracket class.

Principal regexes use explicit amino-acid classes, including unconstrained positions as `[ACDEFGHIKLMNPQRSTVWY]`. They are anchored with `^` and `$`. In the current default run, each regex represents exactly 18 residues and matches the complete parental WT sequence.

## Current Full-Length Regexes

```text
S2-V  WT PQTSRNVNLDSIKLYTSD  ^[ACDEFGHIKLMNPQRSTVWY][ACDEFGHIKLMNPQRSTVWY][ACDEFGHIKLMNPQRSTVWY][ACDEFGHIKLMNPQRSTVWY][ADEFGHKLQRSTVWY][DEGLNPQ][HINPRSVY][ADGHLNPTVY][HLMPW][ADEGHKLNQRSY][CDHISTY][HILPTV][ACFHIKLMNQRSTVY][CDEFHIKLMNTWY][ACDEFGHILMNQRSTWY][ACDEFGHIKLNPRSTVWY][ACDEFGHIKLMNPQRSTVWY][ACDEFGHIKLMNPQRSTVWY]$
S2-IE1  WT PVPATIPLSSVIVAENSD  ^[ACDEFGHIKLMNPQRSTVWY][ACDEFGHIKLMNPQRSTVWY][ACDEFGHIKLMNPQRSTVWY][ACDEFGHIKLMNPQRSTVWY][ACDEFGHIKLMNPQRSTVWY][ADHIKPV][ADEFKLMNPQRSV][FHILMY][ACDEFGHILMPQRSTVWY][ACDKLMNQSTWY][DEFILNPVY][ACEFHIKLMNQRTVWY][CFLMVY][ACDFGHIKLMNPQRSTVWY][ACDEFGHIKLMNPQRSTVWY][ACDEFGHIKLMNPQRSTVWY][ACDEFGHIKLMNPQRSTVWY][ACDEFGHIKLMNPQRSTVWY]$
S3-IL20RA  WT TYLMQFMEEWGLYVQMEN  ^[ACDEFGHIKLMNPQRSTVWY][ACDEFGHIKLMNPQRSTVWY][ACDEFGHIKLMNPQRSTVWY][ACDEFGHIKLMNPQRSTVWY][ADEHIKMNPQRSTWY][CFILMQV][ACDEFILMNPQRSTVWY][ACDEGHKLMNRSTVWY][ADEGHKLMNQRSTW][AGHIMNQSVW][AEGKNPR][ACLMNQR][ACDFGHIKLMNQRSTVWY][HINPQVW][ACDEGHIKMNQRTW][ACDEFGHKLMNPQSW][ACDEFGHIKLMNPQRSTVWY][ACDEFGHIKLMNPQRSTVWY]$
S3-IE1  WT PVPATIPLSSVIVAENSD  ^[ACDEFGHIKLMNPQRSTVWY][ACDEFGHIKLMNPQRSTVWY][ACDEFGHIKLMNPQRSTVWY][ACFGIKMNPRSTVY][ACEGHIKLNPQRSTVWY][DGHIMSVW][AFGKLMPRST][CFHLNSY][ACDEGHIMNPQRSVWY][ADFIKLMNSTWY][DEFHILMNRV][CHILMNPSTVW][CEFHLPQSTV][ACDFGHIKLMNPQRTVWY][CDEGHILMQRSVWY][ACDEFGHIKLMNPQRSTVWY][ACDEFGHIKLMNPQRSTVWY][ACDEFGHIKLMNPQRSTVWY]$
S4-IE1  WT PVPATIPLSSVIVAENSD  ^[ACDEFGHIKLMNPQRSTVWY][ACDEFGHIKLMNPQRSTVWY][ACDEFGHIKLMNPQRSTVWY][ACEFGHILMNPRSTVWY][ACDEFGHIKLMNPQRSTVWY][ADHIKPVW][ACEFLMPR][CDFHLSVY][ACEFGHILMNPQRSTVWY][AFLMNQRSTWY][EFILNTV][ACHIKLMNPQRTVWY][ACEFHILMQRSTVW][ACDFGHIKLMNPQRSVWY][ACDEFGHIKLMNPQRSVWY][ACDEFGHIKLMNPQRSTVWY][ACDEFGHIKLMNPQRSTVWY][ACDEFGHIKLMNPQRSTVWY]$
S4-V  WT PQTSRNVNLDSIKLYTSD  ^[ACDEFGHIKLMNPQRSTVWY][ACDEFGHIKLMNPQRSTVWY][ACDEFGHIKLMNPQRSTVWY][ACDEFGHIKLMNPQRSTVWY][ADEFGHKLNQRSVWY][ADHKLMNPQRV][CFHIKNPQRSVWY][ACHNPQRVY][CFIKLMWY][ADEFGHKLNQRSVW][CDEFHIKMRSTY][DIMRTV][EFGHIKLMNQRTV][CDEFGILMNRTWY][ACDFGHILMRSTWY][ACDEFGHIKLMPQRSTVWY][ACDEFGHIKLMNPQRSTVWY][ACDEFGHIKLMNPQRSTVWY]$
```

## Outputs

The script writes:

- `regex_summary.csv`: one row per parent-domain regex, including WT sequence/window, exact regex, compact regex, analysis positions, analysis length, cutoff, constrained-position count, excluded-substitution count, and missing/unassessable count.
- `residue_decisions.csv`: one row per parent-domain-position-amino-acid decision, including replicate effects, valid flags, paired replicate IDs, denominator, original mean loss, rule flags, decision, reason, coverage/provenance fields, and source.
- `position_summary.csv`: allowed and excluded residues per position, continuous exclusion fraction, designed/usable/unassessable counts, coverage status, and whether the position has a supported constraint.
- `permissive_regexes.txt`: concise human-readable regex list.
- `permissive_regexes.json`: machine-readable regex summary with position-level counts.
- `S2_IE1_independent_validation.csv`: independent S2-IE1 mutation mapping and fixed-rule classification when an independent workbook is supplied.
- `analysis_report.md`: formulas, calibration, validation checks, QC, independent benchmark summary, warnings, constrained-position summary, and reproduction command.

For the current full-sequence run, `residue_decisions.csv` contains `2160` rows: 6 parent-domain pairs x 18 positions x 20 amino acids.

## Independent S2-IE1 Validation

The independent benchmark window is `TIPLSSVIV`, corresponding to construct positions 5-13. Because default output now analyzes positions 1-18, this independent window is inside the analyzed range.

Independent single mutants are mapped to construct positions and annotated with their fixed-rule DMS decision. The independent data are not used to alter residue decisions, set thresholds, choose a core region, or add/remove residues from regexes. A practical independence audit confirmed that `regex_summary.csv`, `residue_decisions.csv`, `position_summary.csv`, `permissive_regexes.txt`, and `permissive_regexes.json` are identical with and without the independent workbook.

## Quality Checks

The script checks that:

- the input uses CLN2-only normalization;
- all six parent-domain combinations exist;
- WT sequences can be inferred for each requested parent;
- WT residues agree between parental WT sequence and DMS annotations wherever DMS rows exist;
- every output parent-domain-position contains all 20 standard amino acids;
- every regex position retains at least one residue;
- every regex position retains its WT residue;
- each exact regex matches its WT sequence/window;
- each default regex represents exactly 18 residues in this dataset;
- no stop codon appears in regex decisions;
- no unassessable residue is excluded;
- the full default decision table has 2160 rows in this dataset;
- S2-IE1 calibration and positions 6-13 matrix-median QC are reported;
- S2-IE1 values in the full and DMS-only workbooks agree when the DMS-only workbook is supplied.

The tests additionally confirm that an explicit positions 6-13 run remains numerically identical to the corresponding positions in the full-sequence default output and that independent activity values do not influence regex-generating outputs.

## Reproducing the Analysis

A local dependency directory `.pythondeps` was used for this run because the system pandas/NumPy install was incompatible. To reproduce the current full-sequence outputs from this workspace:

```bash
PYTHONPATH=.pythondeps PYTHONNOUSERSITE=1 python3 scripts/build_permissive_regex.py \
  simba_t0_t20.xlsx \
  --dms-only-workbook simba_t0_t20_DMS_only.xlsx \
  --independent-benchmark S2-IE1-independant-data.xlsx \
  --output-dir outputs/permissive_regex \
  --full-sequence \
  --tau-ie1 -0.04684
```

The `--full-sequence` flag is optional when no start/end range is supplied, because full-sequence analysis is the default.

Run tests with:

```bash
PYTHONPATH=.pythondeps PYTHONNOUSERSITE=1 python3 -m pytest tests/test_build_permissive_regex.py -q
```

## Architectural Note for Future Upstream Work

Reproducible-benefit processing belongs upstream of this regex generator. Raw fields such as `replicate_DMS_effect_rep1`, `replicate_DMS_effect_rep2`, and `replicate_DMS_effect_rep3` should remain unchanged and auditable. Future processing should preferably create separately named adjusted fields and adjusted DMS/PSSM matrices. `build_permissive_regex.py` should continue using the original paired DMS effects unless its biological method is explicitly changed later.
