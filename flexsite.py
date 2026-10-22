"""FlexSite: residue-level pairwise flexibility analysis inside PyMOL.

Load this file with ``run flexsite.py`` from PyMOL.  The numerical core is
kept independent of PyMOL so that it can be tested from a normal Python
interpreter.
"""

from __future__ import annotations

import csv
import math
import os
import re
from dataclasses import dataclass
from itertools import combinations
from uuid import uuid4
from typing import Dict, Iterable, List, Mapping, MutableMapping, Optional, Sequence, Tuple

try:
    from pymol import cmd, stored
except ImportError:  # Allows unit tests outside PyMOL.
    cmd = None
    stored = None


ResidueKey = Tuple[str, str, str]
Coord = Tuple[float, float, float]
AtomCoords = Dict[str, Coord]
StructureCoords = Dict[ResidueKey, AtomCoords]

BACKBONE_ATOMS = ("N", "CA", "C", "O")
CA_ATOMS = ("CA",)
SIDECHAIN_ATOMS = (
    "CB", "CG", "CG1", "CG2", "CD", "CD1", "CD2", "CE", "CE1", "CE2",
    "CE3", "CZ", "CZ2", "CZ3", "CH2", "ND1", "ND2", "NE", "NE1", "NE2",
    "NH1", "NH2", "NZ", "OD1", "OD2", "OE1", "OE2", "OG", "OG1", "OH",
    "SD", "SG", "SE",
)

# These atom-name pairs describe chemically equivalent positions. Aromatic
# swaps must happen together so that the ring topology remains consistent.
SIDECHAIN_SYMMETRY_SWAPS = {
    "ASP": (("OD1", "OD2"),),
    "GLU": (("OE1", "OE2"),),
    "ARG": (("NH1", "NH2"),),
    "LEU": (("CD1", "CD2"),),
    "VAL": (("CG1", "CG2"),),
    "PHE": (("CD1", "CD2"), ("CE1", "CE2")),
    "TYR": (("CD1", "CD2"), ("CE1", "CE2")),
}


@dataclass(frozen=True)
class ResidueScore:
    chain: str
    resi: str
    resn: str
    metric: str
    mean_pairwise: float
    max_pairwise: float
    std_pairwise: float
    n_pairs: int
    n_structures: int
    total_structures: int
    min_atoms_per_pair: int
    max_atoms_per_pair: int
    max_pair_structure_1: str
    max_pair_structure_2: str

    @property
    def structure_coverage(self) -> float:
        return self.n_structures / self.total_structures if self.total_structures else 0.0

    @property
    def possible_pairs(self) -> int:
        return self.total_structures * (self.total_structures - 1) // 2

    @property
    def pair_coverage(self) -> float:
        return self.n_pairs / self.possible_pairs if self.possible_pairs else 0.0


_LAST_SCORES: List[ResidueScore] = []
_LAST_OUTPUT_OBJECT: Optional[str] = None
_LAST_MIN_STRUCTURE_FRACTION: float = 0.0
_VIEW_SETTINGS: Dict[str, dict] = {}
_OWNER_PROPERTY = "_flexsite_owner"


def _view_settings(output: str) -> dict:
    if stored is not None:
        return dict(getattr(stored, "flexsite_views", {}).get(output, {}))
    return dict(_VIEW_SETTINGS.get(output, {}))


def _remember_view(output: str, site: str, minimum: float, maximum: float,
                   representation: str) -> None:
    settings = dict(site=site, minimum=minimum, maximum=maximum,
                    representation=representation,
                    owner=cmd.get_property(_OWNER_PROPERTY, output))
    _VIEW_SETTINGS[output] = settings
    if stored is not None:
        views = dict(getattr(stored, "flexsite_views", {}))
        views[output] = settings
        stored.flexsite_views = views


def _owns_output(output: str) -> bool:
    owner = _view_settings(output).get("owner")
    return bool(owner and output in cmd.get_names("objects")
                and cmd.get_property(_OWNER_PROPERTY, output) == owner)


def _validate_colorbar_names(output: str) -> None:
    names = cmd.get_names("all")
    bar = f"{output}_colorbar"
    labels = f"{bar}_labels"
    owner = cmd.get_property(_OWNER_PROPERTY, output)
    owned_labels = (
        labels in names and owner
        and cmd.get_type(labels) == "object:molecule"
        and cmd.get_property(_OWNER_PROPERTY, labels) == owner
    )
    if labels in names and not owned_labels:
        raise ValueError(f"Color-bar name belongs to another object: {labels}")
    # CGO objects do not support object properties; verify their type and the
    # marked molecular label companion before deleting or regenerating them.
    if bar in names and (cmd.get_type(bar) != "object:cgo" or not owned_labels):
        raise ValueError(f"Color-bar name belongs to another object: {bar}")


def _validate_output_name(output: str, inputs: Sequence[str]) -> None:
    if not re.fullmatch(r"[A-Za-z0-9_]+", output) or output.lower() in ("all", "none"):
        raise ValueError("output_object must contain only letters, digits, or underscores")
    generated_names = (
        output, f"{output}_site", f"{output}_site_backbone",
        f"{output}_site_sidechain", f"{output}_colorbar", f"{output}_colorbar_labels",
    )
    if any(name in inputs for name in generated_names):
        raise ValueError("output_object must not replace an input structure")
    names = cmd.get_names("all")
    if not _owns_output(output):
        for name in generated_names:
            if name in names:
                raise ValueError(f"Output name already exists: {name}; choose another output_object")
    else:
        for name in generated_names[1:4]:
            if name in names and (
                cmd.get_type(name) != "selection"
                or set(cmd.get_object_list(name)) != {output}
            ):
                raise ValueError(f"Output selection name already exists: {name}")
        _validate_colorbar_names(output)


def _cached_scores() -> List[ResidueScore]:
    if _LAST_SCORES:
        return _LAST_SCORES
    if stored is not None:
        return list(getattr(stored, "flexsite_scores", []))
    return []


def _cached_output_object() -> Optional[str]:
    if _LAST_OUTPUT_OBJECT:
        return _LAST_OUTPUT_OBJECT
    if stored is not None:
        return getattr(stored, "flexsite_output_object", None)
    return None


def _cached_min_structure_fraction() -> float:
    if _LAST_MIN_STRUCTURE_FRACTION:
        return _LAST_MIN_STRUCTURE_FRACTION
    if stored is not None:
        return float(getattr(stored, "flexsite_min_structure_fraction", 0.0))
    return 0.0


def _cache_result(
    scores: Sequence[ResidueScore],
    output_object: str,
    min_structure_fraction: float,
) -> None:
    global _LAST_MIN_STRUCTURE_FRACTION, _LAST_OUTPUT_OBJECT, _LAST_SCORES
    _LAST_SCORES = list(scores)
    _LAST_OUTPUT_OBJECT = output_object
    _LAST_MIN_STRUCTURE_FRACTION = float(min_structure_fraction)
    if stored is not None:
        stored.flexsite_scores = list(scores)
        stored.flexsite_output_object = output_object
        stored.flexsite_min_structure_fraction = float(min_structure_fraction)


def _distance(a: Coord, b: Coord) -> float:
    return math.sqrt(sum((x - y) ** 2 for x, y in zip(a, b)))


def _valid_coord(coord: Sequence[float]) -> bool:
    return len(coord) == 3 and all(math.isfinite(value) for value in coord)


def _blue_white_red(value: float, minimum: float, maximum: float) -> Tuple[float, float, float]:
    """Return an RGB value using a linear blue-white-red palette."""
    if maximum <= minimum:
        return (1.0, 1.0, 1.0)
    fraction = max(0.0, min(1.0, (value - minimum) / (maximum - minimum)))
    if fraction <= 0.5:
        channel = 2.0 * fraction
        return (channel, channel, 1.0)
    channel = 2.0 * (1.0 - fraction)
    return (1.0, channel, channel)


def _matched_atom_pairs(
    first: Mapping[str, Coord],
    second: Mapping[str, Coord],
    atom_names: Sequence[str],
    min_atoms: int,
    residue_name: str = "",
    symmetry_aware: bool = False,
) -> Optional[List[Tuple[str, str, float]]]:
    common = [
        name for name in atom_names
        if name in first and name in second
        and _valid_coord(first[name]) and _valid_coord(second[name])
    ]
    if len(common) < min_atoms:
        return None

    mappings: List[Dict[str, str]] = [{}]
    swaps = SIDECHAIN_SYMMETRY_SWAPS.get(residue_name.upper(), ())
    swap_atoms = {name for pair in swaps for name in pair}
    if symmetry_aware and swaps and swap_atoms.issubset(common):
        swapped: Dict[str, str] = {}
        for left, right in swaps:
            swapped[left] = right
            swapped[right] = left
        mappings.append(swapped)

    candidates: List[List[Tuple[str, str, float]]] = []
    for mapping in mappings:
        pairs = [
            (name, mapping.get(name, name), _distance(first[name], second[mapping.get(name, name)]))
            for name in common
        ]
        candidates.append(pairs)
    return min(candidates, key=lambda pairs: sum(distance * distance for _, _, distance in pairs))


def _pair_rmsd(
    first: Mapping[str, Coord],
    second: Mapping[str, Coord],
    atom_names: Sequence[str],
    min_atoms: int,
    residue_name: str = "",
    symmetry_aware: bool = False,
) -> Optional[Tuple[float, int]]:
    pairs = _matched_atom_pairs(
        first, second, atom_names, min_atoms, residue_name, symmetry_aware
    )
    if pairs is None:
        return None
    squared = [distance * distance for _, _, distance in pairs]
    return math.sqrt(sum(squared) / len(squared)), len(pairs)


def compute_pairwise_scores(
    structures: Mapping[str, StructureCoords],
    atom_names: Sequence[str] = BACKBONE_ATOMS,
    min_atoms: int = 3,
    metric_name: str = "backbone",
    symmetry_aware: bool = False,
) -> List[ResidueScore]:
    """Calculate mean residue RMSD over every valid structure pair.

    Coordinates must already be in a common reference frame and residue keys
    must already be mapped to the reference structure.
    """
    if len(structures) < 2:
        raise ValueError("At least two structures are required")
    if min_atoms < 1 or min_atoms > len(atom_names):
        raise ValueError("min_atoms must be between 1 and len(atom_names)")

    residue_keys = sorted(
        {key for residue_map in structures.values() for key in residue_map},
        key=lambda key: (key[0], _resi_sort_key(key[1]), key[2]),
    )
    scores: List[ResidueScore] = []

    for key in residue_keys:
        observations: List[Tuple[str, AtomCoords]] = [
            (name, residue_map[key])
            for name, residue_map in structures.items()
            if key in residue_map
        ]
        pair_records: List[Tuple[float, int, str, str]] = []
        for (first_name, first), (second_name, second) in combinations(observations, 2):
            result = _pair_rmsd(
                first,
                second,
                atom_names,
                min_atoms,
                residue_name=key[2],
                symmetry_aware=symmetry_aware,
            )
            if result is None:
                continue
            value, atom_count = result
            pair_records.append((value, atom_count, first_name, second_name))

        if not pair_records:
            continue
        pair_values = [record[0] for record in pair_records]
        atom_counts = [record[1] for record in pair_records]
        max_record = max(pair_records, key=lambda record: record[0])
        contributing_structures = {
            name for record in pair_records for name in record[2:4]
        }
        mean = sum(pair_values) / len(pair_values)
        variance = sum((value - mean) ** 2 for value in pair_values) / len(pair_values)
        scores.append(
            ResidueScore(
                chain=key[0],
                resi=key[1],
                resn=key[2],
                metric=metric_name,
                mean_pairwise=mean,
                max_pairwise=max(pair_values),
                std_pairwise=math.sqrt(variance),
                n_pairs=len(pair_values),
                n_structures=len(contributing_structures),
                total_structures=len(structures),
                min_atoms_per_pair=min(atom_counts),
                max_atoms_per_pair=max(atom_counts),
                max_pair_structure_1=max_record[2],
                max_pair_structure_2=max_record[3],
            )
        )
    return scores


def write_scores_csv(scores: Iterable[ResidueScore], filename: str) -> str:
    path = os.path.abspath(os.path.expanduser(filename))
    parent = os.path.dirname(path)
    if parent:
        os.makedirs(parent, exist_ok=True)
    fields = [
        "chain",
        "resi",
        "resn",
        "metric",
        "mean_pairwise_rmsd",
        "max_pairwise_rmsd",
        "std_pairwise_rmsd",
        "n_pairs",
        "n_structures",
        "total_structures",
        "structure_coverage",
        "possible_pairs",
        "pair_coverage",
        "min_atoms_per_pair",
        "max_atoms_per_pair",
        "max_pair_structure_1",
        "max_pair_structure_2",
    ]
    with open(path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for score in scores:
            writer.writerow(
                {
                    "chain": score.chain,
                    "resi": score.resi,
                    "resn": score.resn,
                    "metric": score.metric,
                    "mean_pairwise_rmsd": f"{score.mean_pairwise:.6f}",
                    "max_pairwise_rmsd": f"{score.max_pairwise:.6f}",
                    "std_pairwise_rmsd": f"{score.std_pairwise:.6f}",
                    "n_pairs": score.n_pairs,
                    "n_structures": score.n_structures,
                    "total_structures": score.total_structures,
                    "structure_coverage": f"{score.structure_coverage:.6f}",
                    "possible_pairs": score.possible_pairs,
                    "pair_coverage": f"{score.pair_coverage:.6f}",
                    "min_atoms_per_pair": score.min_atoms_per_pair,
                    "max_atoms_per_pair": score.max_atoms_per_pair,
                    "max_pair_structure_1": score.max_pair_structure_1,
                    "max_pair_structure_2": score.max_pair_structure_2,
                }
            )
    return path


def _require_pymol() -> None:
    if cmd is None:
        raise RuntimeError("This command must be run inside PyMOL")


def _split_values(value: str) -> List[str]:
    return [part for part in re.split(r"[,+\s]+", value.strip()) if part]


def _optional_float(value: object) -> Optional[float]:
    if value is None or str(value).strip() == "":
        return None
    result = float(value)
    if not math.isfinite(result):
        raise ValueError("Numeric limits must be finite")
    return result


def _resi_sort_key(resi: str) -> Tuple[int, str]:
    match = re.match(r"^(-?\d+)(.*)$", str(resi))
    if not match:
        return (10**9, str(resi))
    return (int(match.group(1)), match.group(2))


def _chain_clause(chain: str) -> str:
    return f" and chain {chain}" if chain else ""


def _atom_metadata(selection: str) -> Dict[Tuple[str, int], ResidueKey]:
    metadata: Dict[Tuple[str, int], ResidueKey] = {}
    cmd.iterate(
        selection,
        "metadata[(model, index)] = (chain, resi, resn)",
        space={"metadata": metadata},
    )
    return metadata


def _alignment_residue_map(
    reference: str,
    mobile: str,
    reference_selection: str,
    mobile_selection: str,
    cycles: int,
    cutoff: float,
) -> Dict[ResidueKey, ResidueKey]:
    reference_ca = _alignment_ca_selection(reference_selection)
    mobile_ca = _alignment_ca_selection(mobile_selection)
    fit_name = cmd.get_unused_name(f"_flexsite_fit_{mobile}")
    alignment_name = cmd.get_unused_name(f"_flexsite_map_{mobile}")

    # First transform the mobile onto a stable structural core. Iterative
    # rejection is useful here because genuine domain motions should not drag
    # the common reference frame.
    cmd.align(
        mobile_ca,
        reference_ca,
        cutoff=float(cutoff),
        cycles=int(cycles),
        object=fit_name,
        transform=1,
        mobile_state=1,
        target_state=1,
        quiet=1,
    )
    cmd.delete(fit_name)

    # Build a second, non-transforming sequence alignment without structural
    # outlier rejection. This preserves correspondence for the mobile regions
    # whose movement we actually want to measure.
    cmd.align(
        mobile_ca,
        reference_ca,
        cycles=0,
        object=alignment_name,
        transform=0,
        mobile_state=1,
        target_state=1,
        quiet=1,
    )
    metadata = _atom_metadata(
        f"(({reference}) or ({mobile})) and polymer.protein and name CA"
    )
    residue_map: Dict[ResidueKey, ResidueKey] = {}
    for column in cmd.get_raw_alignment(alignment_name):
        reference_key = None
        mobile_key = None
        for atom_id in column:
            model, _ = atom_id
            if model == reference:
                reference_key = metadata.get(tuple(atom_id))
            elif model == mobile:
                mobile_key = metadata.get(tuple(atom_id))
        if reference_key is not None and mobile_key is not None:
            residue_map[mobile_key] = reference_key
    cmd.delete(alignment_name)
    if not residue_map:
        raise ValueError(f"No aligned CA residue pairs found for {mobile} -> {reference}")
    return residue_map


def _selected_reference_residues(reference_selection: str, site: str) -> List[ResidueKey]:
    residues: MutableMapping[ResidueKey, bool] = {}
    cmd.iterate(
        f"byres (({reference_selection}) and ({site}) and polymer.protein)",
        "residues[(chain, resi, resn)] = True",
        space={"residues": residues},
    )
    return sorted(residues, key=lambda key: (key[0], _resi_sort_key(key[1]), key[2]))


def _extract_structure(
    object_name: str,
    site_keys: Sequence[ResidueKey],
    residue_map: Optional[Mapping[ResidueKey, ResidueKey]] = None,
    require_same_resn: bool = False,
) -> StructureCoords:
    allowed = set(site_keys)
    output: StructureCoords = {}
    model = cmd.get_model(f"({object_name}) and polymer.protein", state=1)
    grouped: Dict[ResidueKey, list] = {}
    for atom in model.atom:
        original_key = (atom.chain, atom.resi, atom.resn)
        reference_key = residue_map.get(original_key) if residue_map is not None else original_key
        if reference_key not in allowed:
            continue
        if require_same_resn and original_key[2] != reference_key[2]:
            continue
        grouped.setdefault(reference_key, []).append(atom)
    for key, atoms in grouped.items():
        coordinates = _conformer_coordinates(atoms)
        if coordinates:
            output[key] = coordinates
    return output


def _conformer_atoms(atoms: Iterable[object]) -> dict:
    """Use shared blank-alt atoms plus one coherent residue conformer.

    A is preferred when available; otherwise use the conformer with the highest
    mean occupancy, with alphabetical tie breaking. Never fill missing atoms
    from a different conformer.
    """
    valid = [
        atom for atom in atoms
        if math.isfinite(float(atom.q)) and float(atom.q) > 0.0
        and _valid_coord(atom.coord)
    ]
    occupancies: Dict[str, List[float]] = {}
    for atom in valid:
        alt = (atom.alt or "").strip()
        if alt:
            occupancies.setdefault(alt, []).append(float(atom.q))
    selected_alt = "A" if "A" in occupancies else min(
        occupancies,
        key=lambda alt: (-sum(occupancies[alt]) / len(occupancies[alt]), alt),
        default="",
    )
    selected = {}
    priorities: Dict[str, Tuple[int, float]] = {}
    for atom in valid:
        alt = (atom.alt or "").strip()
        if alt and alt != selected_alt:
            continue
        name = atom.name.strip().upper()
        priority = (int(alt == ""), float(atom.q))
        if priority > priorities.get(name, (-1, -1.0)):
            selected[name] = atom
            priorities[name] = priority
    return selected


def _conformer_coordinates(atoms: Iterable[object]) -> AtomCoords:
    return {
        name: tuple(float(v) for v in atom.coord)
        for name, atom in _conformer_atoms(atoms).items()
    }


def _alignment_ca_selection(selection: str) -> str:
    grouped: Dict[ResidueKey, list] = {}
    for atom in cmd.get_model(selection, state=1).atom:
        grouped.setdefault((atom.chain, atom.resi, atom.resn), []).append(atom)
    indices = []
    for atoms in grouped.values():
        ca = _conformer_atoms(atoms).get("CA")
        if ca is not None:
            indices.append(str(ca.index))
    if not indices:
        raise ValueError(f"No valid protein CA atoms found in {selection}")
    return f"({selection}) and index {'+'.join(indices)}"


def flexsite_define_site(
    reference: str,
    ligand: str,
    cutoff: float = 5.0,
    name: str = "flexsite_site",
) -> str:
    """Create a reference binding-site selection around a ligand selection."""
    _require_pymol()
    expression = (
        f"byres ((({reference}) and polymer.protein) within {float(cutoff)} of ({ligand}))"
    )
    count = cmd.select(name, expression)
    if count == 0:
        raise ValueError("The binding-site selection is empty; check reference and ligand")
    print(f"FlexSite: selected {count} atoms as {name}")
    return name


def flexsite_analyze(
    reference: str,
    mobiles: str,
    site: str = "flexsite_site",
    ref_chain: str = "",
    mobile_chains: str = "",
    metric: str = "backbone",
    min_atoms: int = 0,
    cycles: int = 5,
    cutoff: float = 2.0,
    min_identity: float = 0.70,
    min_coverage: float = 0.50,
    output_object: str = "",
    csv_file: str = "",
    min_structure_fraction: float = 0.80,
    color_min: str = "",
    color_max: str = "",
    colorbar: int = 1,
) -> List[ResidueScore]:
    """Align an ensemble and color mean pairwise residue differences.

    ``mobiles`` and ``mobile_chains`` accept ``+``, comma, or whitespace as
    separators. All objects are expected to use state 1.
    """
    _require_pymol()
    mobile_objects = _split_values(mobiles)
    if not mobile_objects:
        raise ValueError("Provide at least one mobile object")
    if reference in mobile_objects:
        raise ValueError("Do not repeat the reference in mobiles; it is included automatically")
    if len(set(mobile_objects)) != len(mobile_objects):
        raise ValueError("Each mobile object must appear only once")

    objects = cmd.get_object_list("all")
    for object_name in [reference] + mobile_objects:
        if object_name not in objects:
            raise ValueError(f"Expected a loaded molecular object: {object_name}")

    chain_values = _split_values(mobile_chains)
    if chain_values and len(chain_values) != len(mobile_objects):
        raise ValueError("mobile_chains must contain one chain per mobile object")
    if not chain_values:
        chain_values = [""] * len(mobile_objects)

    metric = metric.lower().strip()
    if metric == "ca":
        atom_names = CA_ATOMS
        effective_min_atoms = 1
        metric_name = "ca_distance"
    elif metric == "backbone":
        atom_names = BACKBONE_ATOMS
        effective_min_atoms = int(min_atoms) if int(min_atoms) else 3
        metric_name = "backbone_rmsd"
        symmetry_aware = False
    elif metric in ("sidechain", "side-chain", "side_chain"):
        metric = "sidechain"
        atom_names = SIDECHAIN_ATOMS
        # Ala has only CB. Gly has no side-chain heavy atom and is therefore
        # correctly omitted rather than assigned an artificial zero.
        effective_min_atoms = int(min_atoms) if int(min_atoms) else 1
        metric_name = "sidechain_rmsd"
        symmetry_aware = True
    else:
        raise ValueError("metric must be 'ca', 'backbone', or 'sidechain'")

    if metric in ("ca", "backbone"):
        symmetry_aware = False

    output = output_object.strip() or f"{reference}_flex_{metric}"
    _validate_output_name(output, [reference] + mobile_objects)
    if effective_min_atoms < 1 or effective_min_atoms > len(atom_names):
        raise ValueError("min_atoms is outside the valid range for this metric")
    min_fraction = float(min_structure_fraction)
    if not 0.0 <= min_fraction <= 1.0:
        raise ValueError("min_structure_fraction must be between 0 and 1")
    for name, value in (("min_identity", min_identity), ("min_coverage", min_coverage)):
        if not 0.0 <= float(value) <= 1.0:
            raise ValueError(f"{name} must be between 0 and 1")
    if int(cycles) < 0 or not math.isfinite(float(cutoff)) or float(cutoff) <= 0.0:
        raise ValueError("cycles must be non-negative and cutoff must be finite and positive")
    fixed_min = _optional_float(color_min)
    fixed_max = _optional_float(color_max)
    if (fixed_min is None) != (fixed_max is None):
        raise ValueError("Provide both color_min and color_max, or leave both empty")
    if fixed_min is not None and fixed_max <= fixed_min:
        raise ValueError("color_max must be greater than color_min")
    show_colorbar = bool(int(colorbar))

    reference_selection = f"({reference}) and polymer.protein{_chain_clause(ref_chain)}"
    site_keys = _selected_reference_residues(reference_selection, site)
    if not site_keys:
        raise ValueError("No protein residues from the reference were found in site")

    structures: Dict[str, StructureCoords] = {
        reference: _extract_structure(reference, site_keys)
    }
    reference_ca_count = cmd.count_atoms(_alignment_ca_selection(reference_selection))
    for mobile, chain in zip(mobile_objects, chain_values):
        mobile_selection = f"({mobile}) and polymer.protein{_chain_clause(chain)}"
        residue_map = _alignment_residue_map(
            reference,
            mobile,
            reference_selection,
            mobile_selection,
            int(cycles),
            float(cutoff),
        )
        identical = sum(
            1 for mobile_key, reference_key in residue_map.items()
            if mobile_key[2] == reference_key[2]
        )
        identity = identical / len(residue_map)
        coverage = len(residue_map) / reference_ca_count if reference_ca_count else 0.0
        print(
            f"FlexSite alignment {mobile} -> {reference}: "
            f"aligned={len(residue_map)}, identity={identity:.1%}, "
            f"reference coverage={coverage:.1%}"
        )
        if identity < float(min_identity) or coverage < float(min_coverage):
            raise ValueError(
                f"{mobile} does not appear to be the same protein/chain as {reference} "
                f"(identity={identity:.1%}, reference coverage={coverage:.1%}). "
                "Check the PDB IDs and chain selections."
            )
        structures[mobile] = _extract_structure(
            mobile, site_keys, residue_map, require_same_resn=(metric == "sidechain")
        )

    scores = compute_pairwise_scores(
        structures,
        atom_names=atom_names,
        min_atoms=effective_min_atoms,
        metric_name=metric_name,
        symmetry_aware=symmetry_aware,
    )
    if not scores:
        raise ValueError("No residue had enough matched atoms for a pairwise calculation")

    display_scores = [score for score in scores if score.structure_coverage >= min_fraction]
    if not display_scores:
        raise ValueError(
            "No residue passed min_structure_fraction; lower the threshold or inspect missing data"
        )

    _color_scores_on_reference(
        reference,
        output,
        display_scores,
        # Source PDB objects may still be enabled when chain-specific working
        # copies are analyzed. Hide the complete scene so they cannot cover the
        # heatmap; the output object is re-enabled by _show_heatmap.
        hide_objects=cmd.get_object_list("all"),
        minimum=fixed_min,
        maximum=fixed_max,
        show_colorbar=show_colorbar,
        representation=metric,
    )
    _cache_result(scores, output, min_fraction)

    if csv_file:
        saved = write_scores_csv(scores, csv_file)
        print(f"FlexSite: wrote {saved}")
    expected_pairs = len(structures) * (len(structures) - 1) // 2
    print(
        f"FlexSite: {len(structures)} structures, {expected_pairs} possible pairs, "
        f"{len(scores)} scored residues, {len(display_scores)} displayed at "
        f">={min_fraction:.0%} structure coverage; output={output}"
    )
    return scores


def _score_selection(object_name: str, score: ResidueScore) -> str:
    chain = f" and chain {score.chain}" if score.chain else ' and chain ""'
    return f"({object_name}{chain} and resi {score.resi})"


def _color_scores_on_reference(
    reference: str,
    output_object: str,
    scores: Sequence[ResidueScore],
    hide_objects: Sequence[str] = (),
    minimum: Optional[float] = None,
    maximum: Optional[float] = None,
    show_colorbar: bool = True,
    representation: str = "backbone",
) -> None:
    _validate_output_name(output_object, [reference])
    owner = _view_settings(output_object).get("owner") if _owns_output(output_object) else uuid4().hex
    # These names were validated above. Remove old representation selections
    # before recreating the object, otherwise unused selections become empty
    # and cannot be distinguished from unrelated names on the next analysis.
    for suffix in ("_site", "_site_backbone", "_site_sidechain"):
        cmd.delete(output_object + suffix)
    cmd.delete(output_object)
    cmd.create(output_object, reference, source_state=1, target_state=1)
    cmd.set_property(_OWNER_PROPERTY, owner, output_object)
    cmd.alter(output_object, "b=-1.0")
    residue_selections = []
    for score in scores:
        selection = _score_selection(output_object, score)
        cmd.alter(selection, f"b={score.mean_pairwise:.10f}")
        residue_selections.append(selection)
    site_name = f"{output_object}_site"
    cmd.select(site_name, " or ".join(residue_selections))
    values = [score.mean_pairwise for score in scores]
    minimum = min(values) if minimum is None else float(minimum)
    maximum = max(values) if maximum is None else float(maximum)
    if math.isclose(minimum, maximum):
        maximum = minimum + 1e-6
    _show_heatmap(
        output_object,
        site_name,
        minimum,
        maximum,
        hide_objects=hide_objects,
        show_colorbar=show_colorbar,
        representation=representation,
    )
    _remember_view(output_object, site_name, minimum, maximum,
                   "sidechain" if representation == "sidechain" else "backbone")


def _create_colorbar(
    output_object: str,
    site_name: str,
    minimum: float,
    maximum: float,
    segments: int = 24,
) -> Optional[str]:
    """Create a vertical 3D CGO color bar beside the analyzed site."""
    try:
        from pymol.cgo import CYLINDER

        (minimum_xyz, maximum_xyz) = cmd.get_extent(site_name)
        x = float(maximum_xyz[0]) + 4.0
        z = float(maximum_xyz[2])
        site_height = max(1.0, float(maximum_xyz[1]) - float(minimum_xyz[1]))
        height = max(6.0, min(14.0, site_height * 0.75))
        y0 = (float(minimum_xyz[1]) + float(maximum_xyz[1]) - height) / 2.0
        radius = max(0.10, min(0.22, height * 0.018))
        cgo: List[float] = []
        for index in range(segments):
            start_fraction = index / segments
            end_fraction = (index + 1) / segments
            start_value = minimum + start_fraction * (maximum - minimum)
            end_value = minimum + end_fraction * (maximum - minimum)
            start_rgb = _blue_white_red(start_value, minimum, maximum)
            end_rgb = _blue_white_red(end_value, minimum, maximum)
            cgo.extend(
                [
                    CYLINDER,
                    x,
                    y0 + start_fraction * height,
                    z,
                    x,
                    y0 + end_fraction * height,
                    z,
                    radius,
                    *start_rgb,
                    *end_rgb,
                ]
            )

        bar_object = f"{output_object}_colorbar"
        label_object = f"{bar_object}_labels"
        cmd.delete(bar_object)
        cmd.delete(label_object)
        cmd.load_cgo(cgo, bar_object)
        for index, (fraction, value) in enumerate(
            ((0.0, minimum), (0.5, (minimum + maximum) / 2.0), (1.0, maximum))
        ):
            cmd.pseudoatom(
                label_object,
                name=f"L{index}",
                pos=(x + 0.7, y0 + fraction * height, z),
                label=f"{value:.2f} A",
            )
        cmd.hide("nonbonded", label_object)
        cmd.set("label_color", "white", label_object)
        cmd.set("label_size", 14, label_object)
        cmd.set_property(_OWNER_PROPERTY, cmd.get_property(_OWNER_PROPERTY, output_object),
                         label_object)
        return label_object
    except Exception as error:
        print(f"FlexSite warning: could not create color bar: {error}")
        return None


def _show_heatmap(
    output_object: str,
    site_name: str,
    minimum: float,
    maximum: float,
    hide_objects: Sequence[str] = (),
    show_surface: bool = False,
    focus: bool = True,
    show_colorbar: bool = True,
    representation: str = "backbone",
) -> None:
    """Apply an uncluttered atom-colored heatmap view."""
    _validate_colorbar_names(output_object)
    cmd.delete(f"{output_object}_colorbar")
    cmd.delete(f"{output_object}_colorbar_labels")
    if not cmd.get_property(_OWNER_PROPERTY, output_object):
        cmd.set_property(_OWNER_PROPERTY, uuid4().hex, output_object)
    for object_name in hide_objects:
        if object_name != output_object:
            cmd.disable(object_name)
    cmd.enable(output_object)
    cmd.hide("everything", output_object)
    cmd.color("gray70", output_object)

    representation = representation.lower().strip()
    if representation == "sidechain":
        display_name = f"{site_name}_sidechain"
        # CA is a drawing anchor for the CA-CB bond but is not part of the
        # side-chain RMSD calculation.
        display_expression = (
            f"({site_name}) and polymer.protein and not hydro "
            "and not name N+C+O+OXT"
        )
        scored_atom_expression = "not hydro and not name N+CA+C+O+OXT"
        if cmd.count_atoms(f"({site_name}) and ({scored_atom_expression})") == 0:
            raise ValueError(f"No binding-site side-chain heavy atoms found in {site_name}")
    else:
        representation = "backbone"
        display_name = f"{site_name}_backbone"
        display_expression = f"({site_name}) and name N+CA+C+O"
    cmd.select(display_name, display_expression)
    if cmd.count_atoms(display_name) == 0:
        raise ValueError(f"No atoms available for {representation} display in {site_name}")

    # Keep the schematic whole-protein cartoon gray. Only the actual N-CA-C-O
    # lines (and an optional surface) should carry the pairwise heatmap.
    cmd.set("cartoon_color", "gray70", output_object)
    cmd.set("stick_color", -1, output_object)
    cmd.set("surface_color", -1, output_object)

    cmd.show("cartoon", output_object)
    # The N-CA-C bonds occupy the same space as the cartoon and would otherwise
    # be hidden inside it. Remove only the scored cartoon segments.
    cmd.hide("cartoon", site_name)
    cmd.show("lines", display_name)
    cmd.unset_bond("line_color", display_name, display_name)
    cmd.unset_bond("line_width", display_name, display_name)
    cmd.set_bond("line_width", 2.5, display_name, display_name)
    cmd.show("sticks", f"({output_object}) and organic")
    if show_surface:
        cmd.show("surface", site_name)
        cmd.set("transparency", 0.20, site_name)
    # Apply spectrum last. Some PyMOL builds assign CPK element colors when a
    # line representation is first shown, which would otherwise overwrite the
    # flexibility colors.
    color_selection = site_name if show_surface else display_name
    cmd.spectrum(
        "b",
        "blue_white_red",
        color_selection,
        minimum=float(minimum),
        maximum=float(maximum),
    )

    # `line_color` is a bond-level setting. Ordinary `set line_color` on an
    # atom selection is accepted by PyMOL but has no visible effect in several
    # builds. Assign a named color to every scored residue and its bonds.
    residue_values: Dict[Tuple[str, str], float] = {}
    cmd.iterate(
        display_name,
        "residue_values[(chain, resi)] = float(b)",
        space={"residue_values": residue_values},
    )
    color_prefix = re.sub(r"[^A-Za-z0-9_]", "_", output_object)
    for index, ((chain, resi), value) in enumerate(
        sorted(residue_values.items(), key=lambda item: (item[0][0], _resi_sort_key(item[0][1])))
    ):
        color_name = f"_flexsite_{color_prefix}_{index:03d}"
        cmd.set_color(color_name, list(_blue_white_red(value, minimum, maximum)))
        chain_clause = f" and chain {chain}" if chain else ' and chain ""'
        if representation == "sidechain":
            residue_selection = (
                f"({output_object}{chain_clause} and resi {resi} and polymer.protein "
                "and not hydro and not name N+C+O+OXT)"
            )
        else:
            residue_selection = (
                f"({output_object}{chain_clause} and resi {resi} and name N+CA+C+O)"
            )
        cmd.color(color_name, residue_selection)
        cmd.set_bond(
            "line_color",
            color_name,
            residue_selection,
            display_name,
        )
    label_object = None
    if show_colorbar:
        label_object = _create_colorbar(output_object, site_name, minimum, maximum)
    if focus:
        focus_selection = (
            f"({site_name}) or ({label_object})" if label_object else site_name
        )
        cmd.zoom(focus_selection, buffer=5.0)
    cmd.rebuild(output_object)
    print(
        f"FlexSite heatmap: blue={minimum:.3f} A, "
        f"white={(minimum + maximum) / 2.0:.3f} A, red={maximum:.3f} A"
    )


def flexsite_show(
    output_object: str = "",
    site_selection: str = "",
    show_surface: int = 0,
    focus: int = 1,
    hide_others: int = 1,
    color_min: str = "",
    color_max: str = "",
    colorbar: int = 1,
    representation: str = "",
) -> None:
    """Redisplay an existing FlexSite object as an uncluttered heatmap.

    This can repair the view after reloading the script without recomputing the
    scores. By default, the site selection is ``<output_object>_site``.
    """
    _require_pymol()
    output = output_object.strip() or (_cached_output_object() or "")
    if not output:
        raise ValueError("Provide the FlexSite output object name")
    settings = _view_settings(output)
    if settings and not _owns_output(output):
        raise ValueError("Cached FlexSite output was deleted or replaced; run a new analysis")
    site = site_selection.strip() or settings.get("site", f"{output}_site")
    if cmd.count_atoms(site) == 0:
        raise ValueError(f"Site selection is empty or missing: {site}")
    values: List[float] = []
    cmd.iterate(site, "values.append(float(b))", space={"values": values})
    values = [value for value in values if value >= 0.0 and math.isfinite(value)]
    if not values:
        raise ValueError(f"No non-negative FlexSite values found in {site}")
    fixed_min = _optional_float(color_min)
    fixed_max = _optional_float(color_max)
    if (fixed_min is None) != (fixed_max is None):
        raise ValueError("Provide both color_min and color_max, or leave both empty")
    if fixed_min is not None and fixed_max <= fixed_min:
        raise ValueError("color_max must be greater than color_min")
    minimum = settings.get("minimum", min(values)) if fixed_min is None else fixed_min
    maximum = settings.get("maximum", max(values)) if fixed_max is None else fixed_max
    if math.isclose(minimum, maximum):
        maximum = minimum + 1e-6
    objects_to_hide = cmd.get_object_list("all") if int(hide_others) else []
    display_mode = representation.lower().strip() or settings.get("representation", "")
    if not display_mode:
        cached = _cached_scores()
        if output == _cached_output_object() and cached:
            display_mode = "sidechain" if cached[0].metric == "sidechain_rmsd" else "backbone"
        else:
            display_mode = "sidechain" if output.lower().endswith("_flex_sidechain") else "backbone"
    if display_mode not in ("backbone", "sidechain"):
        raise ValueError("representation must be 'backbone' or 'sidechain'")
    _show_heatmap(
        output,
        site,
        minimum,
        maximum,
        hide_objects=objects_to_hide,
        show_surface=bool(int(show_surface)),
        focus=bool(int(focus)),
        show_colorbar=bool(int(colorbar)),
        representation=display_mode,
    )
    _remember_view(output, site, minimum, maximum, display_mode)


def flexsite_export(filename: str = "flexsite_results.csv") -> str:
    """Export the most recent analysis to CSV."""
    scores = _cached_scores()
    if not scores:
        raise RuntimeError("No FlexSite result is available; run flexsite_analyze first")
    path = write_scores_csv(scores, filename)
    print(f"FlexSite: wrote {path}")
    return path


def _residue_atom_coordinates(
    object_name: str,
    chain: str,
    resi: str,
    atom_names: Sequence[str],
) -> AtomCoords:
    chain_clause = f" and chain {chain}" if chain else ""
    model = cmd.get_model(
        f"({object_name}){chain_clause} and resi {resi} and polymer.protein",
        state=1,
    )
    residues = {(atom.segi, atom.chain, atom.resi, atom.resn) for atom in model.atom}
    if len(residues) > 1:
        raise ValueError("Residue selection is ambiguous; specify one chain/residue")
    return {
        name: coord for name, coord in _conformer_coordinates(model.atom).items()
        if name in atom_names
    }


def _residue_name(object_name: str, chain: str, resi: str) -> str:
    chain_clause = f" and chain {chain}" if chain else ""
    model = cmd.get_model(
        f"({object_name}){chain_clause} and resi {resi} and polymer.protein",
        state=1,
    )
    return model.atom[0].resn.upper() if model.atom else ""


def flexsite_check_residue(
    structure_1: str,
    structure_2: str,
    resi_1: str,
    chain_1: str = "",
    resi_2: str = "",
    chain_2: str = "",
    atoms: str = "N+CA+C+O",
) -> float:
    """Independently print atom distances and RMSD for one aligned residue pair."""
    _require_pymol()
    sidechain_mode = atoms.lower().strip() in ("sidechain", "side-chain", "side_chain")
    atom_names = SIDECHAIN_ATOMS if sidechain_mode else tuple(
        name.upper() for name in _split_values(atoms)
    )
    second_resi = resi_2.strip() or str(resi_1)
    first = _residue_atom_coordinates(structure_1, chain_1, str(resi_1), atom_names)
    second = _residue_atom_coordinates(structure_2, chain_2, second_resi, atom_names)
    first_resn = _residue_name(structure_1, chain_1, str(resi_1))
    second_resn = _residue_name(structure_2, chain_2, second_resi)
    if sidechain_mode and first_resn != second_resn:
        raise ValueError("Side-chain comparisons require the same amino-acid type")
    residue_name = first_resn if first_resn == second_resn else ""
    pairs = _matched_atom_pairs(
        first,
        second,
        atom_names,
        1,
        residue_name=residue_name,
        symmetry_aware=sidechain_mode,
    )
    if not pairs:
        raise ValueError("No common requested atoms were found for the residue pair")
    squared: List[float] = []
    print(
        f"FlexSite residue check: {structure_1}/{chain_1}/{resi_1} vs "
        f"{structure_2}/{chain_2}/{second_resi}"
    )
    for first_atom, second_atom, distance in pairs:
        squared.append(distance * distance)
        atom_label = first_atom if first_atom == second_atom else f"{first_atom}->{second_atom}"
        print(f"  {atom_label}: {distance:.6f} A")
    rmsd = math.sqrt(sum(squared) / len(squared))
    metric_label = "side-chain" if sidechain_mode else "requested-atom"
    print(f"  direct {metric_label} RMSD ({len(pairs)} atoms): {rmsd:.6f} A")
    return rmsd


def flexsite_adk_report(min_coverage: float = 0.80) -> None:
    """Summarize cached E. coli adenylate-kinase scores by known domains."""
    scores = _cached_scores()
    if not scores:
        raise RuntimeError("No FlexSite result is available; run flexsite_analyze first")
    threshold = float(min_coverage)
    regions: Dict[str, List[ResidueScore]] = {"NMP (30-67)": [], "LID (118-167)": [], "CORE": []}
    skipped = 0
    for score in scores:
        if score.structure_coverage < threshold:
            skipped += 1
            continue
        match = re.match(r"^-?\d+", score.resi)
        if not match:
            continue
        number = int(match.group())
        if 30 <= number <= 67:
            regions["NMP (30-67)"].append(score)
        elif 118 <= number <= 167:
            regions["LID (118-167)"].append(score)
        else:
            regions["CORE"].append(score)

    print(f"FlexSite AdK domain report (coverage >= {threshold:.0%})")
    for region, members in regions.items():
        if not members:
            print(f"  {region}: no binding-site residues")
            continue
        region_mean = sum(score.mean_pairwise for score in members) / len(members)
        highest = max(members, key=lambda score: score.mean_pairwise)
        print(
            f"  {region}: n={len(members)}, mean of residue means={region_mean:.3f} A, "
            f"highest={highest.resn}{highest.resi} ({highest.mean_pairwise:.3f} A)"
        )
    if skipped:
        print(f"  excluded for low coverage: {skipped} residues")


def flexsite_last() -> None:
    """Print a compact summary of the most recent result."""
    scores = _cached_scores()
    if not scores:
        print("FlexSite: no result available")
        return
    print(f"FlexSite output object: {_cached_output_object()}")
    threshold = _cached_min_structure_fraction()
    print(f"Display coverage threshold: {threshold:.0%}")
    for score in sorted(scores, key=lambda value: value.mean_pairwise, reverse=True)[:10]:
        coverage_flag = "" if score.structure_coverage >= threshold else "  [not displayed]"
        print(
            f"{score.chain}:{score.resn}{score.resi}  "
            f"mean={score.mean_pairwise:.3f} A  max={score.max_pairwise:.3f} A  "
            f"pairs={score.n_pairs}/{score.possible_pairs}  "
            f"coverage={score.n_structures}/{score.total_structures}  "
            f"max_pair={score.max_pair_structure_1}|{score.max_pair_structure_2}"
            f"{coverage_flag}"
        )


if cmd is not None:
    cmd.extend("flexsite_define_site", flexsite_define_site)
    cmd.extend("flexsite_analyze", flexsite_analyze)
    cmd.extend("flexsite_export", flexsite_export)
    cmd.extend("flexsite_last", flexsite_last)
    cmd.extend("flexsite_show", flexsite_show)
    cmd.extend("flexsite_check_residue", flexsite_check_residue)
    cmd.extend("flexsite_adk_report", flexsite_adk_report)
    print(
        "FlexSite loaded. Commands: flexsite_define_site, flexsite_analyze, "
        "flexsite_export, flexsite_last, flexsite_show, flexsite_check_residue, "
        "flexsite_adk_report"
    )
