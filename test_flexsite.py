import csv
import math
import os
import tempfile
import unittest

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


if __name__ == "__main__":
    unittest.main()
