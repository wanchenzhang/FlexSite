import importlib.util
from pathlib import Path
import sys

import pymol
from pymol import cmd


pymol.finish_launching(["pymol", "-cqk"])

root = Path(__file__).resolve().parents[1]
tool_path = root / "flexsite.py"
output_dir = root / "results"
output_dir.mkdir(parents=True, exist_ok=True)
spec = importlib.util.spec_from_file_location("flexsite_validation", tool_path)
flexsite = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = flexsite
spec.loader.exec_module(flexsite)

cmd.reinitialize()
for pdb_id in ("1ake", "4ake", "7apu"):
    cmd.load(str(root / "structures" / f"{pdb_id}.cif"), pdb_id)

cmd.create("adk_closed", "1ake and chain A")
cmd.create("adk_open", "4ake and chain A")
cmd.create("adk_adp", "7apu and chain A")
cmd.select("adk_ligand", "1ake and chain A and resn AP5")
flexsite.flexsite_define_site("adk_closed", "adk_ligand", 6.0, "adk_site")
scores = flexsite.flexsite_analyze(
    "adk_closed",
    "adk_open+adk_adp",
    "adk_site",
    min_structure_fraction=0.8,
    csv_file=str(output_dir / "adk_validation.csv"),
)
flexsite.flexsite_last()
flexsite.flexsite_adk_report()
flexsite.flexsite_check_residue(
    "adk_open", "adk_adp", "56", "A", "56", "A", "N+CA+C+O"
)
flexsite.flexsite_show(
    "adk_closed_flex_backbone", color_min="0", color_max="12", colorbar=1
)
cmd.png(str(output_dir / "adk_validation.png"), 1200, 900, ray=0)
print(f"VALIDATION_RESIDUES={len(scores)}")
cmd.quit()
