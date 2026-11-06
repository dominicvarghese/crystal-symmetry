"""Optional independent reference check: install spglib, then run this script.

spglib is a development check only; the browser app has no Python dependency.
"""
import json
from pathlib import Path
import subprocess
import numpy as np
import spglib

root = Path(__file__).resolve().parent.parent
extract = r"""
const fs=require('node:fs'),vm=require('node:vm');
const html=fs.readFileSync('index.html','utf8');
const source=[...html.matchAll(/<script\b[^>]*>([\s\S]*?)<\/script>/g)].map(m=>m[1]).find(s=>s.includes('function analyzeSingleRep'));
const ctx=vm.createContext({});vm.runInContext(source,ctx);
const data=Object.fromEntries(Object.entries(ctx.PRESET_CRYSTALS).map(([key,d])=>[key,{
  ...d,site_point_groups:d.atoms.map(a=>ctx.Crystallography.getAtomSiteSymmetry(d,a).pointGroup.hm)
}]));
process.stdout.write(JSON.stringify(data));
"""
presets = json.loads(subprocess.check_output(["node", "-e", extract], cwd=root, text=True))
print(f"Independent reference: spglib {spglib.__version__}")
for key, data in presets.items():
    lattice = np.array(data["cell_vectors"])
    cell = (lattice, [a["frac"] for a in data["atoms"]], [a["atomic_number"] for a in data["atoms"]])
    ref = spglib.get_symmetry_dataset(cell, symprec=data["symprec"])
    assert ref is not None, key
    assert ref.number == data["spacegroup_number"], (key, ref.number)
    assert len(ref.rotations) == len(data["operations"]), key
    for op in data["operations"]:
        matches = []
        for rotation, translation in zip(ref.rotations, ref.translations):
            if np.array_equal(rotation, op["W"]):
                delta = np.array(op["w_raw"]) - translation
                delta -= np.round(delta)
                matches.append(np.linalg.norm(delta @ lattice) <= data["symprec"])
        assert any(matches), (key, op["xyz"])
    for actual, expected in zip(data["site_point_groups"], ref.site_symmetry_symbols):
        # spglib includes axis-position dots; mm2 may be expressed in another orientation.
        expected = expected.replace(".", "")
        expected = {"m2m": "mm2", "2mm": "mm2"}.get(expected, expected)
        assert actual == expected, (key, actual, expected)
    print(f"{key}: {ref.international} #{ref.number}; all operations and site point groups agree")
