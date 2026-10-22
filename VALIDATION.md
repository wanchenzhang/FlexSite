# FlexSite AdK validation

The reliability update is described in [CHANGELOG.md](CHANGELOG.md). The
archived numerical results below were reproduced with the repaired code.
New exports include `possible_pairs` and `pair_coverage`; the original root CSV
files retain their historical schema.

FlexSite was tested with chain A from `1AKE` (closed reference), `4AKE`
(open), and `7APU` (ADP-bound). The binding site was defined as reference
residues within 6 A of AP5 in `1AKE`.

## Alignment and coverage

- `adk_open`: 214 aligned residues, 100% identity, 100% reference coverage.
- `adk_adp`: 214 aligned residues, 100% identity, 100% reference coverage.
- 55 binding-site residues were scored and all 55 passed the default 80%
  structure-coverage threshold.

The structural superposition uses a stable core fit. A separate unrefined
sequence alignment supplies residue correspondence, so mobile LID and NMP
residues are not lost as structural outliers.

## Domain-level result

Using the mean pairwise backbone RMSD of each scored residue:

| Region | Residues | Mean of residue scores | Highest residue |
|---|---:|---:|---|
| NMP (30-67) | 16 | 5.244 A | GLY56, 10.574 A |
| LID (118-167) | 15 | 5.904 A | VAL132, 9.140 A |
| CORE | 24 | 0.958 A | ALA11, 1.674 A |

The expected moving NMP and LID domains are therefore substantially more
variable than the CORE in this test ensemble.

## Independent residue check

For GLY56 in `adk_open` versus `adk_adp`, the direct aligned-atom distances
were:

- N: 16.353532 A
- CA: 16.026193 A
- C: 15.264448 A
- O: 15.498405 A

Their root mean square is 15.791460 A. This exactly matches the maximum
pairwise backbone RMSD stored for GLY56, and the CSV correctly identifies
`adk_open | adk_adp` as the producing pair.

## Reproducible artifacts

- `adk_validation.csv`: full numerical output, coverage, and maximum-pair IDs.
- `adk_validation.png`: fixed 0-12 A blue-white-red scale with a labeled color
  bar.
- `adk_sidechain_validation.csv`: side-chain heavy-atom result for the same
  three-structure ensemble.
- `adk_sidechain_validation.png`: side-chain line representation with a fixed
  0-16 A scale and labeled color bar.
- `examples/validate_adk.py`: portable backbone validation script.
- `examples/validate_adk_sidechain.py`: portable side-chain validation script.

Run these scripts with a Python interpreter that provides PyMOL. They load the
bundled mmCIF files from `structures/` and write regenerated artifacts to
`results/`, preserving the original validation files.

## Side-chain validation

The side-chain analysis scored 46 binding-site residues; glycine residues were
omitted because they contain no side-chain heavy atom. LYS57 was highest, with
a mean pairwise RMSD of 12.557313 A. Its maximum was produced by `adk_open` and
`adk_adp`.

An independent calculation for that pair used CB, CG, CD, CE, and NZ distances
of 14.165584, 15.609756, 18.391821, 20.169420, and 22.580064 A. Their RMS is
18.435033 A, exactly matching the tool's stored maximum. The validation also
reconstructed the view using `flexsite_show`, confirming that cached side-chain
results retain the correct line representation and color scale.

This validation exercises the prototype's calculation, correspondence mapping,
filtering, reporting, and visualization on this test set. The manual checker
shares atom-matching code with the analysis and is a diagnostic, not an
independent implementation. Numerical fixtures and archived AdK values are also
checked by the regression suite. This does not establish a general biological interpretation
for every protein ensemble; chain selection, ligand/site definition, alignment
protocol, crystal packing, mutations, and missing atoms still need to be
controlled for each study.
