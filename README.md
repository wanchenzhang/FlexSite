# FlexSite

FlexSite compares the same binding-site residues across an ensemble of protein
structures. It calculates every valid structure pair, averages the residue-level
pairwise difference, and colors the result on a copy of a selected reference
structure.

This is an MVP intended for structures of the same protein. It supports C-alpha
distances, backbone RMSD (`N`, `CA`, `C`, `O`), and side-chain heavy-atom RMSD.
Side-chain mode handles alternate conformations, missing atoms, and common
chemically equivalent atom-name swaps.

## Requirements

- Python 3.7 or newer.
- A working PyMOL installation with Python scripting support for analysis and
  visualization. The Python interpreter used for example scripts must be able
  to `import pymol`.
- The numerical unit tests require only the Python standard library.
- An internet connection is needed for `fetch`; bundled mmCIF files allow the
  AdK examples to run offline.

FlexSite measures conformational variability across static structures. Its
side-chain scores include domain and backbone displacement; they are not a
measurement of isolated rotamer motion or dynamic transition rates.

## Install and load

Download or clone this repository. In PyMOL, change to the downloaded project
directory, then load the script:

```pml
cd /path/to/FlexSite
run flexsite.py
```

Replace `/path/to/FlexSite` with the actual directory on your computer. You can
also use PyMOL's File > Run menu to select `flexsite.py`.

Seven commands are registered:

- `flexsite_define_site`: define binding-site residues around a ligand selection;
- `flexsite_analyze`: align structures, calculate all pairs, and draw the heatmap;
- `flexsite_export`: export the most recent result as CSV;
- `flexsite_last`: print the ten most flexible residues;
- `flexsite_show`: isolate, recolor, and focus an existing heatmap object;
- `flexsite_check_residue`: independently recalculate one aligned residue pair;
- `flexsite_adk_report`: summarize an AdK result by NMP, LID, and CORE regions.

## Minimal example

The following example uses three adenylate kinase (AdK) structures that were
tested with this version of FlexSite. `1AKE` supplies the closed reference and
the AP5 ligand, while `4AKE` and `7APU` provide alternative conformations:

```pml
fetch 1ake, async=0
fetch 4ake, async=0
fetch 7apu, async=0
create adk_closed, 1ake and chain A
create adk_open, 4ake and chain A
create adk_adp, 7apu and chain A
```

Inspect the reference and select its intended ligand. Do not assume that every
`organic` object is the biologically relevant ligand:

```pml
select adk_ligand, 1ake and chain A and resn AP5
flexsite_define_site adk_closed, adk_ligand, 6.0, adk_site
```

Analyze chain A in all three structures:

```pml
flexsite_analyze adk_closed, adk_open+adk_adp, adk_site, A, A+A, backbone
flexsite_last
flexsite_export flexsite_results.csv
```

By default, only residues for which at least 80% of the input structures
contribute to a valid atom-matched pair are colored. Residues with too few usable
atoms do not count as contributing structures. Lower-coverage residues remain
gray but are retained in the CSV. The
CSV and `flexsite_last` also identify the two structures that produce each
residue's maximum pairwise RMSD.

The reference is included automatically. Do not include it again in `mobiles`.
The output object is named `adk_closed_flex_backbone` by default, and its scored
site selection is named `adk_closed_flex_backbone_site`.

Choose an output name containing only letters, digits, and underscores. The
tool refuses to replace an input structure or an unrelated existing object.
Reusing a FlexSite output name replaces that result. Alignment transforms the
mobile objects in place, so use working copies if you need the original poses.
When upgrading from an earlier version, start a fresh PyMOL session and rerun
the analysis so that new output ownership markers and coverage values are used.

After analysis, the input ensemble is disabled so it cannot cover the heatmap.
The output is focused on the binding site. Backbone RMSD is shown by default as
thin, atom-colored `N-CA-C-O` lines on a gray cartoon; side-chain atoms are not
used as the colored representation. Blue is the minimum pairwise value, white is
the middle of the current range, and red is the maximum; the numerical range is
printed in the PyMOL console.

A labeled blue-white-red color bar is created automatically. Its object names
are `<output>_colorbar` and `<output>_colorbar_labels`.

## Fixed scale and coverage filtering

The automatic scale is useful within one run, but separate analyses must use the
same fixed range for visual comparison. Use the Python block in PyMOL:

```pml
python
flexsite_show(
    output_object="adk_closed_flex_backbone",
    color_min="0",
    color_max="12",
    colorbar=1,
)
python end
```

`flexsite_show` remembers each output's representation and color range within
the current PyMOL process, including after reloading `flexsite.py`. Omitting the
object name restores the last analyzed output. These cached settings are not
guaranteed to survive closing PyMOL or a new process; run the analysis again in
that case.

To require complete structural coverage during a new analysis:

```pml
python
flexsite_analyze(
    reference="adk_closed",
    mobiles="adk_open+adk_adp",
    site="adk_site",
    min_structure_fraction=1.0,
    color_min="0",
    color_max="12",
)
python end
```

## Validation and AdK domain report

After the AdK analysis:

```pml
flexsite_adk_report
flexsite_check_residue adk_open, adk_adp, 56, A, 56, A
```

The first command groups binding-site residues into the NMP domain (30-67), LID
domain (118-167), and CORE. The second prints the direct `N`, `CA`, `C`, and `O`
distances and their RMSD. It assumes the structures are already aligned and the
two structures use the residue numbers supplied by the user.

If an older run is still covered by the original structures, reload the updated
script and repair the view without recomputing:

```pml
flexsite_show adk_closed_flex_backbone
```

Results calculated by a version from before full sequence mapping was added
must instead be recomputed. Recoloring cannot recover residues that the older
alignment mapping rejected as structural outliers.

To add the optional transparent whole-residue surface overlay:

```pml
flexsite_show adk_closed_flex_backbone, adk_closed_flex_backbone_site, 1
```

For C-alpha pairwise distances:

```pml
flexsite_analyze adk_closed, adk_open+adk_adp, adk_site, A, A+A, ca
```

## Side-chain flexibility

Use `sidechain` as the metric:

```pml
flexsite_analyze adk_closed, adk_open+adk_adp, adk_site, A, A+A, sidechain
flexsite_last
flexsite_export adk_sidechain.csv
```

The score uses only common non-hydrogen side-chain atoms; backbone atoms do not
contribute. The drawing includes CA as an anchor so the CA-CB connection and an
alanine CB remain visible, but CA is not included in the side-chain RMSD. Glycine
has no side-chain heavy atom and is omitted. By default one common side-chain
atom is sufficient; use the Python API's `min_atoms` argument to require more.

Equivalent terminal or branch names are symmetry-corrected for ASP, GLU, ARG,
LEU, VAL, PHE, and TYR. For example, swapping ASP `OD1` and `OD2` labels alone
does not create an artificial difference.

Only pairs with identical amino-acid types are compared in side-chain mode.
Mutated positions in a mobile structure are excluded from that residue's
side-chain observations, while backbone mode may still compare their backbones.
Symmetry swaps require the complete equivalent atom group in both structures;
incomplete groups use conservative same-name matching.

To independently check one aligned side chain:

```pml
flexsite_check_residue adk_open, adk_adp, 57, A, 57, A, sidechain
```

An existing side-chain result can be restored with:

```pml
flexsite_show adk_closed_flex_sidechain
```

Arguments can also be supplied from PyMOL's Python console, which is preferable
when paths or selections contain spaces:

```python
flexsite_analyze(
    reference="adk_closed",
    mobiles="adk_open+adk_adp",
    site="adk_site",
    ref_chain="A",
    mobile_chains="A+A",
    metric="backbone",
    csv_file="flexsite_results.csv",
)
```

## What the score means

For residue `i` and structures `p` and `q`:

```text
RMSD_i(p,q) = sqrt(mean_a(||x_i,a,p - x_i,a,q||^2))
```

The reported `mean_pairwise_rmsd` is the arithmetic mean of this value over all
valid `N(N-1)/2` structure pairs. For the `ca` metric, there is only one atom, so
the value is a C-alpha distance rather than an atom-set RMSD.

The CSV also reports maximum, population standard deviation, number of valid
pairs, number of contributing structures, and the minimum/maximum atom count
used by any valid pair. A backbone pair requires at least three of the four
requested backbone atoms by default. `possible_pairs` is the total number of
ensemble pairs, while `pair_coverage = n_pairs / possible_pairs` exposes how many
comparisons actually contributed. Structure coverage and pair coverage differ:
every structure can participate in at least one pair without all pairs being
valid. The default display threshold applies to structure coverage only; inspect
pair coverage and atom counts before comparing residues with missing data.

## Scientific and technical limitations

- All structures are independently aligned to the reference using PyMOL
  `align` on C-alpha atoms. The structural fit uses iterative outlier rejection
  to establish a stable core, while a separate sequence alignment without
  structural rejection preserves residue correspondence in moving domains. The
  reference controls the common coordinate frame, while every structure,
  including the reference, contributes to the pairwise score.
- PyMOL's default iterative alignment rejects structural outliers. This is useful
  for finding a common core, but the alignment protocol should be agreed with the
  supervisor before interpreting the result biologically.
- Only state 1 is analyzed.
- Corresponding chains must be supplied for multichain objects.
- By default, each mobile must have at least 70% sequence identity and cover at
  least 50% of the selected reference chain. This prevents unrelated PDB entries
  from producing a visually plausible but meaningless heatmap.
- Residues are mapped from the C-alpha alignment. Unmapped or missing residues
  do not contribute to a pair.
- Shared blank-alt atoms are combined with one consistent residue conformer.
  Alternate `A` is preferred; if it is absent, the highest mean occupancy wins
  (alphabetical tie breaking). Missing atoms are not filled from another
  conformer. The same choice is used for coordinate extraction and C-alpha
  alignment; zero-occupancy and non-finite atoms are excluded.
- Side-chain RMSD is calculated after the same global structural alignment as
  backbone RMSD; no per-residue local fit is applied. It therefore captures
  both domain/backbone displacement and side-chain conformational change.
- A pair contributes only its common side-chain atom names. The CSV's
  `min_atoms_per_pair` and `max_atoms_per_pair` fields should be inspected when
  missing density differs between structures.
- Raw crystallographic B-factors are not used. The output is a copy of the
  reference, and its B-factor field is repurposed only for visualization.
- The colors are scaled to the minimum and maximum value in the current result;
  two separate runs therefore need a shared fixed scale before their colors can
  be compared visually.

The latest score table is cached in PyMOL's persistent `stored` namespace, so
re-running `flexsite.py` does not discard `flexsite_last` or `flexsite_export`.
This persistence begins with the first analysis performed by this version.

## Run the numerical tests outside PyMOL

```powershell
python -m unittest discover -v
```

With a standard Python interpreter, numerical tests run and PyMOL integration
tests are explicitly skipped. To run the complete suite, use a Python interpreter
that can import PyMOL. GitHub Actions runs numerical tests on Windows and Linux
with Python 3.9 and 3.13; it does not install PyMOL.

## Reproduce the AdK example

From the repository directory, using PyMOL's Python interpreter:

```sh
python examples/validate_adk.py
python examples/validate_adk_sidechain.py
```

These scripts load the bundled files from `structures/` and generate CSV and PNG
files in `results/`. The original example outputs in the repository root are
retained as a reference. For manual offline use in PyMOL, replace the three
`fetch` commands above with:

```pml
load structures/1ake.cif, 1ake
load structures/4ake.cif, 4ake
load structures/7apu.cif, 7apu
```

## Structure sources and references

The bundled coordinate files come from RCSB PDB and retain their original
metadata. Chain A is used in the example analyses.

- [1AKE](https://www.rcsb.org/structure/1AKE), Ap5A-bound reference:
  Müller CW, Schulz GE (1992). Journal of Molecular Biology 224, 159–177.
  [DOI](https://doi.org/10.1016/0022-2836(92)90582-5).
- [4AKE](https://www.rcsb.org/structure/4AKE), unligated open structure:
  Müller CW, Schlauderer GJ, Reinstein J, Schulz GE (1996). Structure 4, 147–156.
  [DOI](https://doi.org/10.1016/S0969-2126(96)00018-4).
- [7APU](https://www.rcsb.org/structure/7APU), two-ADP-bound structure:
  Ojeda-May P et al. (2021). Biochemistry 60, 2246–2258.
  [DOI](https://doi.org/10.1021/acs.biochem.1c00221).

See [VALIDATION.md](VALIDATION.md) for this project's numerical results. These
results are generated by FlexSite and are distinct from the cited publications.

## Repository contents

- `flexsite.py`: PyMOL commands and numerical core.
- `test_flexsite.py`: numerical tests.
- `test_pymol_integration.py`: real PyMOL and AdK regression tests.
- `examples/`: portable AdK validation scripts.
- `structures/`: source mmCIF files for offline reproduction.
- `adk*.csv` and `adk*.png`: original validation data and figures.
- `VALIDATION.md`: validation methods, results, and interpretation.
- `CHANGELOG.md`: repair notes and verification scope.

No software license has been selected yet. Before reusing or distributing the
code, obtain permission from the rights holder. PDB data have their own terms;
see the [wwPDB usage policy](https://www.wwpdb.org/about/usage-policies).
