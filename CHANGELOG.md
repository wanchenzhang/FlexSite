# Changelog

## Unreleased — reliability fixes (2026-10-22)

### Fixed

- Structure coverage now counts only structures contributing to a valid pair,
  rather than every occurrence of a residue key.
- Side-chain extraction skips mapped residues with a different amino-acid type.
- Coordinate extraction no longer mixes mutually exclusive alternate conformers.
  Zero-occupancy and non-finite atoms are excluded. C-alpha fitting uses the
  same conformer policy and one atom per residue.
- `flexsite_show` restores the last output after a script reload and remembers
  the representation and color scale of each output, including custom names.
- Constant-valued heatmaps can be displayed again without a range error.
- Disabling a color bar removes stale bars even with `hide_others=0`.
- Output names cannot delete input structures or unrelated existing objects.
  Object ownership is checked against a marker on the current molecular object,
  so stale cache entries do not authorize overwriting a replacement object.
  Color-bar and label name collisions are also rejected.
- Duplicate mobile names and invalid numeric settings fail before alignment.
- Manual residue checks reject ambiguous multichain selections.

### Added

- `possible_pairs` and `pair_coverage` columns in new CSV exports; console output
  now reports valid/possible pairs.
- Portable examples and source mmCIF files previously referenced by the README
  but absent from the uploaded repository.
- Numerical tests and real PyMOL regression tests, including AdK missing-atom
  and mutated-residue cases.
- A GitHub Actions workflow for standard-library tests on Windows and Linux.

### Verification scope

The complete suite was run locally with PyMOL 3.1.8. The three-structure AdK
ensemble retains 55 backbone scores and 46 side-chain scores. GLY56's maximum
backbone RMSD remains 15.791460 A and LYS57's maximum side-chain RMSD remains
18.435033 A. Both headless examples generated CSV and PNG outputs. Results are
specific to this ensemble and do not establish accuracy for every protein.

The GitHub Actions configuration was added locally; a hosted workflow run has
not been observed. No software license was selected as part of these repairs.

### Remaining limitations

- CGO color-bar ownership is inferred from its type and its marked label
  companion because PyMOL CGOs do not support object properties. Give manually
  substituted CGOs a different name rather than reusing a generated bar name.
- Residue keys contain chain, residue identifier, and residue name; repeated
  identifiers across segments remain unsupported. Analyze distinct chain copies.
