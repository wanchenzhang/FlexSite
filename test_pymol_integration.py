"""Real PyMOL regression tests; skipped when the interpreter lacks PyMOL."""

import importlib
from pathlib import Path
import unittest

import flexsite


@unittest.skipIf(flexsite.cmd is None, "PyMOL is not installed in this interpreter")
class PyMOLIntegrationTests(unittest.TestCase):
    def setUp(self):
        self.cmd = flexsite.cmd
        self.cmd.reinitialize()
        flexsite._LAST_SCORES = []
        flexsite._LAST_OUTPUT_OBJECT = None
        flexsite.stored.flexsite_views = {}
        self.cmd.fragment("lys", "reference")
        self.cmd.create("mobile", "reference")
        atom = self.cmd.get_model("reference and name CA").atom[0]
        self.key = (atom.chain, atom.resi, atom.resn)
        self.cmd.select("site", "reference and polymer.protein")
        self.scores = flexsite.compute_pairwise_scores({
            "reference": {self.key: {"CB": (0.0, 0.0, 0.0)}},
            "mobile": {self.key: {"CB": (1.0, 0.0, 0.0)}},
        }, atom_names=flexsite.SIDECHAIN_ATOMS, min_atoms=1, metric_name="sidechain_rmsd")

    def draw(self):
        flexsite._color_scores_on_reference(
            "reference", "custom_result", self.scores,
            minimum=0.0, maximum=5.0, representation="sidechain",
        )
        flexsite._cache_result(self.scores, "custom_result", 0.8)

    def tearDown(self):
        self.cmd.reinitialize()

    def test_reload_restores_custom_sidechain_name_without_arguments(self):
        self.draw()
        importlib.reload(flexsite)
        flexsite.flexsite_show(focus=0)
        self.assertGreater(self.cmd.count_atoms("custom_result_site_sidechain"), 0)

    def test_restore_preserves_fixed_color_range(self):
        self.draw()
        flexsite.flexsite_show("custom_result", focus=0)
        labels = []
        self.cmd.iterate("custom_result_colorbar_labels", "labels.append(label)",
                         space={"labels": labels})
        self.assertIn("5.00 A", labels)

    def test_automatic_constant_scores_can_be_redisplayed(self):
        flexsite._color_scores_on_reference("reference", "constant_result", self.scores)
        flexsite.flexsite_show("constant_result", focus=0)

    def test_colorbar_can_be_disabled_without_hiding_other_objects(self):
        self.draw()
        flexsite.flexsite_show("custom_result", colorbar=0, hide_others=0, focus=0,
                              color_min="0", color_max="5")
        self.assertNotIn("custom_result_colorbar", self.cmd.get_names("objects", enabled_only=1))

    def test_output_name_cannot_replace_reference(self):
        before = self.cmd.count_atoms("reference")
        with self.assertRaises(ValueError):
            flexsite.flexsite_analyze("reference", "mobile", "site", output_object="reference")
        self.assertEqual(self.cmd.count_atoms("reference"), before)

    def test_duplicate_mobile_names_are_rejected(self):
        with self.assertRaises(ValueError):
            flexsite.flexsite_analyze("reference", "mobile+mobile", "site")

    def test_unrelated_existing_object_is_not_overwritten(self):
        self.cmd.fragment("ala", "keep_me")
        before = self.cmd.count_atoms("keep_me")
        with self.assertRaises(ValueError):
            flexsite.flexsite_analyze("reference", "mobile", "site", output_object="keep_me")
        self.assertEqual(self.cmd.count_atoms("keep_me"), before)

    def test_colorbar_name_collision_cannot_delete_input(self):
        self.cmd.set_name("mobile", "result_colorbar")
        before = self.cmd.count_atoms("result_colorbar")
        with self.assertRaises(ValueError):
            flexsite.flexsite_analyze("reference", "result_colorbar", "site",
                                     output_object="result")
        self.assertEqual(self.cmd.count_atoms("result_colorbar"), before)

    def test_alignment_temporary_names_do_not_delete_existing_objects(self):
        self.cmd.fragment("ala", "_flexsite_fit_mobile")
        self.cmd.fragment("ala", "_flexsite_map_mobile")
        flexsite.flexsite_analyze("reference", "mobile", "site", colorbar=0)
        self.assertGreater(self.cmd.count_atoms("_flexsite_fit_mobile"), 0)
        self.assertGreater(self.cmd.count_atoms("_flexsite_map_mobile"), 0)

    def test_stale_cache_does_not_authorize_replacing_unrelated_object(self):
        self.draw()
        self.cmd.delete("custom_result")
        self.cmd.fragment("ala", "custom_result")
        before = self.cmd.count_atoms("custom_result")
        with self.assertRaises(ValueError):
            flexsite.flexsite_analyze("reference", "mobile", "site",
                                     output_object="custom_result")
        self.assertEqual(self.cmd.count_atoms("custom_result"), before)

    def test_owned_output_can_be_recomputed(self):
        self.draw()
        scores = flexsite.flexsite_analyze("reference", "mobile", "site",
                                          output_object="custom_result", colorbar=0)
        self.assertEqual(len(scores), 1)
        self.assertGreater(self.cmd.count_atoms("custom_result"), 0)

    def test_switching_metric_allows_repeated_output_reuse(self):
        self.draw()
        for metric in ("backbone", "backbone", "sidechain", "sidechain"):
            scores = flexsite.flexsite_analyze(
                "reference", "mobile", "site", output_object="custom_result",
                metric=metric, colorbar=0,
            )
            self.assertEqual(len(scores), 1)

    def test_replaced_colorbar_is_protected_when_hiding_bar(self):
        self.draw()
        self.cmd.delete("custom_result_colorbar")
        self.cmd.fragment("ala", "custom_result_colorbar")
        before = self.cmd.count_atoms("custom_result_colorbar")
        with self.assertRaises(ValueError):
            flexsite.flexsite_show("custom_result", colorbar=0)
        self.assertEqual(self.cmd.count_atoms("custom_result_colorbar"), before)

    def test_replaced_labels_are_protected_when_reanalyzing(self):
        self.draw()
        self.cmd.delete("custom_result_colorbar_labels")
        self.cmd.fragment("ala", "custom_result_colorbar_labels")
        before = self.cmd.count_atoms("custom_result_colorbar_labels")
        with self.assertRaises(ValueError):
            flexsite.flexsite_analyze("reference", "mobile", "site",
                                     output_object="custom_result")
        self.assertEqual(self.cmd.count_atoms("custom_result_colorbar_labels"), before)

    def test_invalid_threshold_does_not_transform_mobile(self):
        self.cmd.translate([10.0, 0.0, 0.0], "mobile")
        before = self.cmd.get_atom_coords("mobile and name CA")
        with self.assertRaises(ValueError):
            flexsite.flexsite_analyze("reference", "mobile", "site",
                                     min_structure_fraction=1.5)
        self.assertEqual(self.cmd.get_atom_coords("mobile and name CA"), before)

    def test_actual_pymol_run_preserves_last_result(self):
        self.draw()
        self.cmd.run(str(Path(flexsite.__file__).resolve()))
        self.cmd.keyword["flexsite_show"][0](focus=0)
        self.assertGreater(self.cmd.count_atoms("custom_result_site_sidechain"), 0)

    def test_each_output_remembers_its_own_representation(self):
        self.draw()
        flexsite._color_scores_on_reference("reference", "backbone_result", self.scores)
        flexsite._cache_result(self.scores, "backbone_result", 0.8)
        importlib.reload(flexsite)
        flexsite.flexsite_show("custom_result", focus=0)
        self.assertNotIn("custom_result_site_backbone", self.cmd.get_names("selections"))

    def test_manual_check_rejects_ambiguous_multichain_residue(self):
        self.cmd.create("other_chain", "reference")
        self.cmd.alter("other_chain", "chain='B'")
        self.cmd.create("multichain", "reference or other_chain")
        with self.assertRaises(ValueError):
            flexsite.flexsite_check_residue("multichain", "mobile", self.key[1])

    def test_alignment_uses_one_ca_per_alternate_residue(self):
        self.cmd.remove("reference and hydro")
        self.cmd.alter("reference and name CA", "alt='A'")
        model = self.cmd.get_model("reference")
        import copy
        alternate_ca = copy.deepcopy(next(atom for atom in model.atom if atom.name == "CA"))
        alternate_ca.alt = "B"
        alternate_ca.coord[0] += 1.0
        model.add_atom(alternate_ca)
        self.cmd.delete("reference")
        self.cmd.load_model(model, "reference")
        selection = flexsite._alignment_ca_selection("reference and polymer.protein")
        self.assertEqual(self.cmd.count_atoms(selection), 1)
        self.assertEqual(self.cmd.get_model(selection).atom[0].alt, "A")


@unittest.skipIf(flexsite.cmd is None, "PyMOL is not installed in this interpreter")
class AdKReferenceTests(unittest.TestCase):
    """Pin the real alignment/extraction pipeline to the archived AdK results."""

    def setUp(self):
        self.cmd = flexsite.cmd
        self.cmd.reinitialize()
        structures = Path(__file__).resolve().parent / "structures"
        for name in ("1ake", "4ake", "7apu"):
            self.cmd.load(str(structures / (name + ".cif")), name)
        for output, original in (("closed", "1ake"), ("opened", "4ake"), ("adp", "7apu")):
            self.cmd.create(output, original + " and chain A")
        flexsite.flexsite_define_site("closed", "1ake and chain A and resn AP5", 6.0, "site")

    def tearDown(self):
        self.cmd.reinitialize()

    def analyze(self, metric="backbone"):
        return flexsite.flexsite_analyze("closed", "opened+adp", "site", metric=metric, colorbar=0)

    def test_backbone_pipeline_matches_archived_gly56(self):
        scores = self.analyze()
        self.assertEqual(len(scores), 55)
        score = next(score for score in scores if score.resi == "56")
        self.assertAlmostEqual(score.max_pairwise, 15.791460, places=4)
        self.assertEqual(score.n_pairs, 3)
        self.assertEqual(score.structure_coverage, 1.0)

    def test_sidechain_pipeline_matches_archived_lys57(self):
        scores = self.analyze("sidechain")
        self.assertEqual(len(scores), 46)
        top = max(scores, key=lambda score: score.mean_pairwise)
        self.assertEqual((top.resn, top.resi), ("LYS", "57"))
        self.assertAlmostEqual(top.mean_pairwise, 12.557313, places=4)
        self.assertAlmostEqual(top.max_pairwise, 18.435033, places=4)

    def test_incomplete_backbone_is_retained_but_not_colored(self):
        self.cmd.remove("opened and resi 56 and name N+O")
        scores = self.analyze()
        score = next(score for score in scores if score.resi == "56")
        self.assertEqual(score.n_structures, 2)
        self.assertEqual(score.n_pairs, 1)
        b_values = []
        self.cmd.iterate("closed_flex_backbone and chain A and resi 56",
                         "b_values.append(b)", space={"b_values": b_values})
        self.assertTrue(b_values)
        self.assertTrue(all(value == -1.0 for value in b_values))

    def test_mutated_sidechain_is_excluded_from_coverage(self):
        self.cmd.alter("opened and resi 57", "resn='ARG'")
        scores = self.analyze("sidechain")
        score = next(score for score in scores if score.resi == "57")
        self.assertEqual(score.n_structures, 2)
        self.assertEqual(score.n_pairs, 1)


if __name__ == "__main__":
    unittest.main()
