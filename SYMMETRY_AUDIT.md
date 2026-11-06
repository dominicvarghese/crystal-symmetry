# Crystal symmetry implementation review

Reviewed upstream commit `b138f7bde77f7022b4df1ee7eb9480a64dc2022a` on 2026-11-06.

## ZnS question

The observed corner-to-edge-midpoint motion is valid for the displayed offset screw axis. A half-turn about `(1/2, 1/4, z)` maps `(0,0,0)` to `(1,1/2,0)`; translating by `(0,0,1/2)` gives `(1,1/2,1/2)`. The transverse displacement reverses about the screw axis, not necessarily the center of a cell face. About the face-center axis `(1/2,1/2,z)`, the same corner would rotate to `(1,1,0)`. The intermediate position need not be an occupied atomic site: only the complete screw operation must map the crystal onto itself.

The inspected unstrained operations already had correct animation endpoints. The fix adds an explicit A → rotation/reflection milestone → B coordinate path and explains the axis anchor.

## Changes

| Issue | Previous behavior | Correction |
| --- | --- | --- |
| Quartz oxygen preset | Incorrect z values in the second oxygen triplet removed all twofold operations; detected point group `3` incorrectly implied a polar crystal | Corrected triplet using `(y,x,-z)`; six operations, PG `32`, Si site PG `2`, nonpolar |
| Screw handedness labels | Absolute pitch and an asymmetric sign correction could interchange `3_1/3_2` or `4_1/4_3` | Signed axial translation modulo the rotation order, with the axis oriented consistently with the rotation |
| Decorative screw helix | Arbitrary 1.5 turns regardless of screw translation and cell size | Angular advance derives from the actual axial rise and operation angle; arrow label says step translation |
| Tolerance and hydrogen | Low symmetry could silently trigger 0.12 Å tolerance or use a hydrogen-free skeleton | Every species must match at the requested atomic tolerance |
| Atom correspondence | Independent nearest-atom matches could assign multiple source atoms to the same target | Bipartite matching requires a bijective species-preserving permutation |
| CIF operation bypass | Extra explicit operations could be appended without a distance or metric check | Every explicit operation must preserve the cell and map all atoms bijectively within tolerance |
| Near-site merging | Fixed 0.06 Å deduplication could erase distinct sites even at a smaller selected tolerance | Deduplication capped at the smaller of the selected tolerance and 0.001 Å |
| Cell metric | A dimensionless 0.01 allowance admitted incompatible rotations under small strains | Strict 1e-6 metric threshold, also used by strain classification |
| Broken-operation animation | Nonorthogonal `A W A^-1` under strain disagreed with the rigid Rodrigues animation | Use the rigid parent operation in the lab frame for metric-breaking strain; disclose coordinate basis |
| Pair matching and mismatch | Fractional 0.02/0.26 cutoffs varied in physical size with the lattice; affine endpoints could hide strain mismatch | Compare actual animation endpoints in Cartesian Å, using selected tolerance; broken pairs retain parent correspondence |
| Wyckoff labels | Orbit letters assigned by input order were presented as crystallographic Wyckoff labels | UI calls them local atomic-orbit IDs; no claim of International Tables assignments |
| CIF numeric and loop parsing | Fabricated default cell, quoted numbers mishandled, one physical line assumed per row, malformed xyz terms partially accepted | Require valid cell parameters, parse quoted values, regroup loop tokens, handle inline comments, reject invalid xyz and incomplete rows |
| Missing scripts | Three nonexistent local scripts requested as an advertised fallback | Removed missing references and documented the actual CDN dependency |

## Validation

`node tests/symmetry.test.cjs`: 19 passing tests. These include all seven preset counts and point groups, operation closure modulo lattice translations, bijective atom mappings, endpoint equivalence to the affine operator for both animation styles and all packed atoms, the ZnS half-turn geometry, quartz site symmetry, signed screws, input rejection and strain rigidity. Known independent-component counts are checked for polarization, dielectric, piezoelectric and elastic tensors in ZnS, diamond and quartz. Every inline script is also compiled.

Expected operation counts: ZnS 96; diamond 192; rutile 16; chalcopyrite 16; quartz 6; baddeleyite 4; rubrene 16.

A Chromium/Playwright smoke test also loaded every preset, selected a screw operation, inspected its coordinate path, scrubbed the timeline, aligned the camera, selected an atomic site and exercised the strain tab without JavaScript errors. The test served cached copies of the pinned Three.js/OrbitControls libraries; CDN availability was not tested.

## Remaining limitations and recommended next work

This is a targeted correctness review, not certification of the full crystallography workbench.

1. Replace the finite `{-1,0,1}` integer-matrix search and componentwise minimum-image rounding with a crystallographic engine such as spglib in WASM. General unreduced/skew cells can need larger matrix coefficients and a metric-aware closest lattice vector. Benchmark equivalent settings against spglib and published CIFs.
2. A CIF space-group name/number alone is not enough for the current parser to expand a general asymmetric unit. Only explicit operations or a small hardcoded fallback table are supported. Require explicit operations for those uploads, or supply a complete setting-aware database. Metadata labels are supplied labels, not independently verified space-group identification.
3. CIF support remains a subset: multiple data blocks, multiline scalar values, disorder and occupancy-aware chemistry need dedicated handling. Do not treat missing/partial occupancies as fully occupied without a clear model.
4. Genuine Wyckoff assignments need a setting/origin-aware space-group database. Keep local orbit IDs until that is implemented.
5. The named order-parameter modes are illustrative displacement/strain patterns, not general ISODISTORT or Bilbao irrep calculations. Parent-dependent irrep labels, full subgroup identification and domain enumeration need independent reference validation. Magnetic/time-reversal tensor restrictions are not included.
6. Atomic tolerance can produce approximate operation sets that are not closed for arbitrary noisy uploads. The tests prove closure for presets; a future importer should validate closure and report/refine inconsistent approximate groups. The strict lattice threshold now favors exact rigid geometry and may report lower symmetry for noisy cell parameters.
7. Vendor pinned Three.js assets if fully offline use is required. The single HTML file currently needs jsDelivr at runtime.
