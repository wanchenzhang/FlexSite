import csv
import math
import os
import tempfile
import unittest
from types import SimpleNamespace
from unittest.mock import patch

import flexsite

from flexsite import (
    SIDECHAIN_ATOMS,
    _blue_white_red,
    compute_pairwise_scores,
    write_scores_csv,
)


def residue(offset, missing_o=False):
    atoms = {
        "N": (0.0, offset, 0.0),
        "CA": (1.0, offset, 0.0),
        "C": (2.0, offset, 0.0),
        "O": (3.0, offset, 0.0),
    }
    if missing_o:
        del atoms["O"]
    return atoms


class PairwiseScoreTests(unittest.TestCase):
    def test_blue_white_red_palette(self):
        self.assertEqual(_blue_white_red(0.0, 0.0, 2.0), (0.0, 0.0, 1.0))
        self.assertEqual(_blue_white_red(1.0, 0.0, 2.0), (1.0, 1.0, 1.0))
        self.assertEqual(_blue_white_red(2.0, 0.0, 2.0), (1.0, 0.0, 0.0))

    def test_three_structure_mean_pairwise_backbone_rmsd(self):
        key = ("A", "10", "GLY")
        structures = {
            "s1": {key: residue(0.0)},
            "s2": {key: residue(1.0)},
            "s3": {key: residue(2.0)},
        }
        score = compute_pairwise_scores(structures)[0]
        self.assertAlmostEqual(score.mean_pairwise, 4.0 / 3.0)
        self.assertAlmostEqual(score.max_pairwise, 2.0)
        self.assertEqual(score.n_pairs, 3)
        self.assertEqual(score.n_structures, 3)
        self.assertEqual(score.total_structures, 3)
        self.assertEqual(score.structure_coverage, 1.0)
        self.assertEqual(score.min_atoms_per_pair, 4)
        self.assertEqual(
            (score.max_pair_structure_1, score.max_pair_structure_2),
            ("s1", "s3"),
        )

    def test_ca_metric_is_distance(self):
        key = ("A", "10", "GLY")
        structures = {
            "s1": {key: residue(0.0)},
            "s2": {key: residue(3.0)},
        }
        score = compute_pairwise_scores(
            structures, atom_names=("CA",), min_atoms=1, metric_name="ca_distance"
        )[0]
        self.assertAlmostEqual(score.mean_pairwise, 3.0)
        self.assertEqual(score.metric, "ca_distance")

    def test_missing_oxygen_uses_three_common_backbone_atoms(self):
        key = ("A", "10", "GLY")
        structures = {
            "s1": {key: residue(0.0)},
            "s2": {key: residue(1.0, missing_o=True)},
        }
        score = compute_pairwise_scores(structures, min_atoms=3)[0]
        self.assertAlmostEqual(score.mean_pairwise, 1.0)
        self.assertEqual(score.min_atoms_per_pair, 3)

    def test_residue_with_only_one_observation_is_skipped(self):
        shared = ("A", "10", "GLY")
        unique = ("A", "11", "ALA")
        structures = {
            "s1": {shared: residue(0.0), unique: residue(0.0)},
            "s2": {shared: residue(1.0)},
        }
        scores = compute_pairwise_scores(structures)
        self.assertEqual([(score.resi, score.n_pairs) for score in scores], [("10", 1)])

    def test_sidechain_alanine_uses_cb_and_glycine_is_skipped(self):
        ala = ("A", "10", "ALA")
        gly = ("A", "11", "GLY")
        structures = {
            "s1": {ala: {"CB": (0.0, 0.0, 0.0)}, gly: residue(0.0)},
            "s2": {ala: {"CB": (0.0, 2.0, 0.0)}, gly: residue(1.0)},
        }
        scores = compute_pairwise_scores(
            structures,
            atom_names=SIDECHAIN_ATOMS,
            min_atoms=1,
            metric_name="sidechain_rmsd",
            symmetry_aware=True,
        )
        self.assertEqual(len(scores), 1)
        self.assertEqual(scores[0].resn, "ALA")
        self.assertAlmostEqual(scores[0].mean_pairwise, 2.0)
        self.assertEqual(scores[0].min_atoms_per_pair, 1)

    def test_sidechain_symmetry_avoids_false_aspartate_motion(self):
        key = ("A", "20", "ASP")
        first = {
            "CB": (0.0, 0.0, 0.0),
            "CG": (1.0, 0.0, 0.0),
            "OD1": (2.0, 1.0, 0.0),
            "OD2": (2.0, -1.0, 0.0),
        }
        second = dict(first)
        second["OD1"], second["OD2"] = first["OD2"], first["OD1"]
        structures = {"s1": {key: first}, "s2": {key: second}}
        raw = compute_pairwise_scores(
            structures,
            atom_names=SIDECHAIN_ATOMS,
            min_atoms=1,
            metric_name="sidechain_rmsd",
        )[0]
        corrected = compute_pairwise_scores(
            structures,
            atom_names=SIDECHAIN_ATOMS,
            min_atoms=1,
            metric_name="sidechain_rmsd",
            symmetry_aware=True,
        )[0]
        self.assertGreater(raw.mean_pairwise, 1.0)
        self.assertAlmostEqual(corrected.mean_pairwise, 0.0)

    def test_csv_export(self):
        key = ("A", "10", "GLY")
        scores = compute_pairwise_scores(
            {"s1": {key: residue(0.0)}, "s2": {key: residue(1.0)}}
        )
        with tempfile.TemporaryDirectory() as directory:
            path = write_scores_csv(scores, os.path.join(directory, "scores.csv"))
            with open(path, newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle))
        self.assertEqual(rows[0]["chain"], "A")
        self.assertEqual(rows[0]["mean_pairwise_rmsd"], "1.000000")
        self.assertEqual(rows[0]["structure_coverage"], "1.000000")
        self.assertEqual(rows[0]["max_pair_structure_1"], "s1")
        self.assertEqual(rows[0]["max_pair_structure_2"], "s2")

    def test_coverage_counts_only_structures_in_valid_pairs(self):
        key = ("A", "10", "GLY")
        score = compute_pairwise_scores({
            "complete1": {key: residue(0.0)},
            "complete2": {key: residue(1.0)},
            "incomplete": {key: {"CA": (0.0, 0.0, 0.0)}},
        })[0]
        self.assertEqual(score.n_structures, 2)
        self.assertAlmostEqual(score.structure_coverage, 2.0 / 3.0)

    def test_pair_coverage_reports_missing_comparisons(self):
        key = ("A", "10", "GLY")
        score = compute_pairwise_scores({
            "s1": {key: residue(0.0)},
            "s2": {key: residue(1.0)},
            "s3": {},
        })[0]
        self.assertEqual(score.possible_pairs, 3)
        self.assertAlmostEqual(score.pair_coverage, 1.0 / 3.0)

    def test_nonfinite_coordinates_do_not_produce_nan_scores(self):
        key = ("A", "10", "GLY")
        damaged = residue(1.0)
        damaged["O"] = (math.nan, 0.0, 0.0)
        score = compute_pairwise_scores({
            "s1": {key: residue(0.0)}, "s2": {key: damaged},
        })[0]
        self.assertAlmostEqual(score.mean_pairwise, 1.0)
        self.assertEqual(score.min_atoms_per_pair, 3)

    def test_invalid_color_limits_are_rejected(self):
        for value in ("nan", "inf", "-inf"):
            with self.subTest(value=value):
                with self.assertRaises(ValueError):
                    flexsite._optional_float(value)


class ExtractionTests(unittest.TestCase):
    def atom(self, name, alt, coord, occupancy=1.0, resn="SER"):
        return SimpleNamespace(
            name=name, alt=alt, coord=coord, q=occupancy,
            chain="A", resi="10", resn=resn,
        )

    def extract(self, atoms, **kwargs):
        model = SimpleNamespace(atom=atoms)
        with patch.object(flexsite, "cmd", SimpleNamespace(get_model=lambda *a, **k: model)):
            return flexsite._extract_structure("mobile", [("A", "10", "SER")], **kwargs)

    def test_alternate_conformers_are_not_mixed_atom_by_atom(self):
        coordinates = self.extract([
            self.atom("CA", "", (0.0, 0.0, 0.0)),
            self.atom("CB", "A", (1.0, 0.0, 0.0)),
            self.atom("CB", "B", (2.0, 0.0, 0.0)),
            self.atom("OG", "B", (3.0, 0.0, 0.0)),
        ])[('A', '10', 'SER')]
        self.assertEqual(coordinates["CB"], (1.0, 0.0, 0.0))
        self.assertNotIn("OG", coordinates)

    def test_zero_occupancy_atoms_are_not_observations(self):
        self.assertEqual(self.extract([
            self.atom("CB", "", (1.0, 0.0, 0.0), occupancy=0.0),
        ]), {})

    def test_non_a_alternate_uses_highest_residue_occupancy(self):
        coordinates = self.extract([
            self.atom("CB", "B", (1.0, 0.0, 0.0), occupancy=0.8),
            self.atom("CB", "C", (2.0, 0.0, 0.0), occupancy=0.2),
        ])[('A', '10', 'SER')]
        self.assertEqual(coordinates["CB"], (1.0, 0.0, 0.0))

    def test_sidechain_mapping_skips_mutated_residue(self):
        mapping = {("A", "10", "THR"): ("A", "10", "SER")}
        self.assertEqual(self.extract([
            self.atom("CB", "", (1.0, 0.0, 0.0), resn="THR"),
        ], residue_map=mapping, require_same_resn=True), {})


if __name__ == "__main__":
    unittest.main()
