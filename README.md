# 🔮 3D Crystal Symmetry, Site-Symmetry, ISOTROPY & Jahn Tensor Lab

👉 **Live Interactive Web App:** [https://dominicvarghese.github.io/crystal-symmetry/](https://dominicvarghese.github.io/crystal-symmetry/)

An interactive, browser-based 3D crystallography and solid-state physics workbench for visualizing **space-group symmetry operations**, **local atomic site symmetries**, **ISOTROPY / Bilbao order-parameter & strain-induced phase transitions**, and **Jahn symbol physical property tensors** directly from `.cif` or VASP `POSCAR`/`CONTCAR` files.

---

## ✨ Key Features

### 1. 🔮 Interactive 2-Stage 3D Symmetry Element Visualizer (Tab 1)
Standard crystal viewers (such as VESTA) display crystal structures but do not geometrically render or animate 3D compound space-group operations inside the unit cell. This tool decomposes every space-group operation $(W, \mathbf{w})$ following *International Tables for Crystallography, Volume A* (Wondratschek & Neubüser) into its **exact 3D geometric element** and **two-stage physical motion ($0\% \rightarrow 50\% \rightarrow 100\%$):**

* **Rotoinversions ($\bar{3}, \bar{4}, \bar{6}$):**
  * **Stage 1 ($0\% \rightarrow 50\%$):** Rotates the atom (and its 4-colored chiral reference tetrahedron) by $120^\circ, 90^\circ,$ or $60^\circ$ around the rotoinversion axis to an intermediate wireframe milestone (`[1]`).
  * **Stage 2 ($50\% \rightarrow 100\%$):** Inverts radially through the glowing rotoinversion point $\mathbf{r}_0$ on the axis, flipping chirality/handedness (`+1 → -1`) and landing on the target crystal site (`[B]`).
* **Glide Planes ($a, b, c, n, d$):**
  * Renders the exact crystallographic plane clipped to the unit cell walls with in-plane glide vectors.
  * **Stage 1 ($0\% \rightarrow 50\%$):** Reflects perpendicularly across the plane (flipping handedness).
  * **Stage 2 ($50\% \rightarrow 100\%$):** Slides parallel to the plane by the fractional glide vector $\mathbf{w}_\parallel$ (e.g., $\frac{1}{2}\mathbf{c}$ for $c$-glide, $\frac{1}{2}\mathbf{a}+\frac{1}{2}\mathbf{b}$ for diagonal $n$-glide, or $\frac{1}{4}\mathbf{a}\pm\frac{1}{4}\mathbf{b}\pm\frac{1}{4}\mathbf{c}$ for diamond $d$-glide).
* **Screw Axes ($2_1, 3_1, 3_2, 4_1, 4_2, 4_3, 6_1, \dots$):**
  * Renders the axis rod, polygon symmetry badges, a 3D helical spiral ribbon showing screw handedness, and the pitch vector arrow $\mathbf{w}_\parallel$.
  * Animates both **2-Stage Step-by-Step** (*Rotate $\rightarrow$ Slide along axis*) and **Simultaneous Helical Motion**.
* **Proper Rotations ($2, 3, 4, 6$), Mirror Planes ($m$), and Inversion Centers ($\bar{1}$).**
* **3D Chiral Tetrahedron Motif:** Attaches an asymmetric 4-colored chiral tetrahedron to animated atoms so you can visually distinguish proper operations ($\det W = +1$, handedness preserved) from improper operations ($\det W = -1$, enantiomorphic handedness flip).

---

### 2. 🔬 Atomic Site-Symmetry & Coordination Shell Inspector
* **Atomic Orbit & Site Point Group Table:** Identifies the exact local **Site-Symmetry Point Group** (in both Hermann-Mauguin and Schoenflies notation, e.g., `-43m (Td)`, `mmm (D2h)`, `2/m (C2h)`, `2 (C2)`, `-1 (Ci)`, `1 (C1)`) for every computed atomic orbit in the crystal.
* **Single-Atom Site Explorer (`🔬 Site Sym (1 Atom)`):** Click any atom in 3D (or any row in the Wyckoff table) to isolate the stabilizer subgroup of symmetry operations passing directly through that atom, highlight its complete first-shell coordination polyhedron, and overlay all intersecting site-symmetry axes and mirror planes on that atom.
* **Atom-to-Atom SymOp Finder (`🎯 Pair A→B (2 Atoms)`):** Click any two equivalent atoms $A$ and $B$ in the unit cell to filter the operation list to only the symmetry operations that map Atom $A$ directly onto Atom $B$.

---

### 3. ⚡ ISOTROPY / Bilbao Order Parameter & Strain Lab (Tab 2)
Inspired by the **ISOTROPY Software Suite (`ISODISTORT`)** and the **Bilbao Crystallographic Server (`AMPLIMODES` / `SUBGROUPGRAPH`)**, Tab 2 allows experimentalists and theorists to apply macroscopic strain tensors $\varepsilon_{ij}$ and microscopic optical-phonon order parameters $\eta$ in real time:

* **Built-in Irrep & Order Parameter Direction (OPD) Modes:**
  1. **$\Gamma_1^+ (A_{1g})$** — Hydrostatic Volume Strain ($e_a$: preserves full parent symmetry $G_0$)
  2. **$\Gamma_3^+ (E_g)$ OPD $P1(a,0)$** — Uniaxial Tetragonal Bain Strain ($\varepsilon_{zz} \neq \varepsilon_{xx}=\varepsilon_{yy}$, e.g., $F\bar{4}3m \rightarrow I\bar{4}m2$, $Fd\bar{3}m \rightarrow I4_1/amd$)
  3. **$\Gamma_3^+ / B_{1g}$ OPD $C1(a,b)$** — Orthorhombic Strain ($\varepsilon_{xx} = -\varepsilon_{yy}$, e.g., Rutile $P4_2/mnm \rightarrow \text{CaCl}_2$-type $Pnnm$)
  4. **$\Gamma_5^+ (T_{2g} / B_{2g})$ OPD $P1(0,0,a)$** — In-Plane Shear Strain ($\varepsilon_{xy} \neq 0$, $\gamma \neq 90^\circ$)
  5. **$\Gamma_5^+ (T_{2g})$ OPD $P3(a,a,a)$** — Rhombohedral $[111]$ Trigonal Shear ($\varepsilon_{xy}=\varepsilon_{yz}=\varepsilon_{xz}$, e.g., $F\bar{4}3m \rightarrow R3m$, $Fd\bar{3}m \rightarrow R\bar{3}m$)
  6. **$\Gamma_5^+ (E_g)$ OPD $C(0,a,0)$** — Monoclinic Shear Strain ($\varepsilon_{xz} \neq 0$, $\beta \neq 90^\circ$)
  7. **$\Gamma_4^- (T_{1u} / A_{2u})$ OPD $P1(0,0,a)$** — Zone-Center Polar Ferroelectric Mode ($\mathbf{P} \parallel [001]$)
  8. **$\Gamma_4^- (T_{1u})$ OPD $P3(a,a,a)$** — Rhombohedral Polar Ferroelectric Mode ($\mathbf{P} \parallel [111]$)
  9. **$\Gamma_3^+ / B_{1g}$** — Antiferrodistortive Octahedral / Sublattice Rotation ($\phi_z$ around $[001]$)
  * *Plus manual sliders for all 6 Voigt strain components ($\varepsilon_{xx}, \varepsilon_{yy}, \varepsilon_{zz}, \varepsilon_{xy}, \varepsilon_{xz}, \varepsilon_{yz}$), polar vectors ($P_x, P_y, P_z$), and internal rotations ($\phi_z$).*
* **Live Group $\rightarrow$ Subgroup & Domain Report:**
  * Displays the **Parent Point Group ($G_0$) $\rightarrow$ Distorted Isotropy Subgroup ($G_{\text{sub}}$)**, crystal system change, **Subgroup Index $[G_0 : G_{\text{sub}}]$** (= number of ferroelastic/ferroelectric twin domains), and physical selection rules (**Centrosymmetric**, **Polar / Ferroelectric $\mathbf{P}_s$**, **Piezoelectric / SHG**, **Chiral / Enantiomorphic**).
* **3D Visualization of 🔴 Broken Symmetry Operations:**
  * Filters parent operations into **🟢 Kept ($g \in G_{\text{sub}}$)** vs. **🔴 Broken ($g \in G_0 \setminus G_{\text{sub}}$)**. Clicking any broken symmetry operation and pressing **Play** shows the rigid operation missing the strained/displaced crystal atom in 3D with a red mismatch vector **`Δr` (in Å)**.

---

### 4. 📐 Jahn Symbol & Physical Property Tensor Calculator (Neumann's Principle)
Following *International Tables for Crystallography, Volume D* and Bilbao's **`MTENSOR`**, the visualizer constructs the Reynolds group-averaging projection operator and Reduced Row Echelon Form (RREF) directly from the active point-group matrices $\{R_g\}$ to generate the **exact symbolic non-zero/zero tensor matrix**:

* **11 Built-in Jahn Symbols:**
  * **`V`** ($1\times 3$): Pyroelectricity / Spontaneous Polarization ($P_i$)
  * **`eV`** ($1\times 3$): Axial Vector / Ferromagnetic Moment ($M_i$)
  * **`[V²]`** ($3\times 3$): Dielectric Permittivity ($\varepsilon_{ij}$), Thermal Expansion ($\alpha_{ij}$), Conductivity ($\sigma_{ij}$)
  * **`e[V²]`** ($3\times 3$): Optical Activity / Gyration Tensor ($g_{ij}$)
  * **`V²`** ($3\times 3$): General 2nd-Rank Polar / Thermoelectric Tensor ($T_{ij}$)
  * **`eV²`** ($3\times 3$): Magnetoelectric / Axial 2nd-Rank Tensor ($G_{ij}$)
  * **`V[V²]`** ($3\times 6$ Voigt): Piezoelectricity ($d_{i\alpha}$) & Second-Harmonic Generation SHG ($\chi^{(2)}_{i\alpha}$)
  * **`[V³]`** ($3\times 6$ Voigt): Kleinman-Symmetric Nonlinear Optical SHG ($\beta_{i\alpha}$)
  * **`eV[V²]`** ($3\times 6$ Voigt): Electrogyration / Piezomagnetic Tensor ($\gamma_{i\alpha}$)
  * **`[[V²]²]`** ($6\times 6$ Voigt): Elastic Stiffness ($C_{\alpha\beta}$) & Compliance ($S_{\alpha\beta}$)
  * **`[V²]²`** ($6\times 6$ Voigt): Elasto-Optic / Photoelastic ($p_{\alpha\beta}$), Electrostriction, Piezoresistance
* **Color-Coded Entries:**
  * **Cyan:** Independent free tensor parameters
  * **Green:** Symmetry-constrained equal/opposite or linear combinations (e.g., $C_{66} = \frac{1}{2}C_{11} - \frac{1}{2}C_{12}$)
  * **★ Gold Border:** Tensor components that were strictly `0` in the unstrained parent crystal $G_0$ and became **unlocked (non-zero)** by the applied strain or order parameter!
* **3 Tensor Scopes:** Evaluate tensors for the **Active Distorted Subgroup ($G_{\text{sub}}$)**, **Unstrained Parent Group ($G_0$)**, or **Selected Atomic Site Symmetry ($G_{\text{site}}$)**.

---

## 📂 Supported File Formats & Built-in Presets

* **Drag-and-Drop File Support:**
  * **`.cif` (Crystallographic Information File):** Supports both Inorganic and Organic/Molecular crystals (with symmetry expansion from explicit operations and validation against all atomic species).
  * **`POSCAR` / `CONTCAR` (VASP 4 & 5):** Detects candidate space-group and site-symmetry operations directly from the raw lattice vectors and atomic coordinates using an adjustable tolerance (`tol: 0.005 Å – 0.10 Å`).
* **7 Curated Built-in Presets:**
  1. **Sphalerite $\text{ZnS}$ ($F\bar{4}3m$, #216)** — $\bar{4}$ rotoinversion, $d$-glides, $2_1$ screws
  2. **Diamond $\text{Si}$ ($Fd\bar{3}m$, #227)** — $4_1/4_3$ screws, $d$-glides, $\bar{3}$ rotoinversions, $\bar{1}$ centers
  3. **Rutile $\text{TiO}_2$ ($P4_2/mnm$, #136)** — $4_2$ screw, $n$-glides, $\bar{4}$, octahedral modes
  4. **Chalcopyrite $\text{CuFeS}_2$ ($I\bar{4}2d$, #122)** — $\bar{4}$ rotoinversions, $d$-glides
  5. **$\alpha$-Quartz $\text{SiO}_2$ ($P3_121$, #152)** — Chiral $3_1$ screw axis, $2$-fold site symmetry
  6. **Baddeleyite $\text{ZrO}_2$ ($P2_1/c$, #14)** — Monoclinic $2_1$ screw, $c$-glide, $\bar{1}$
  7. **Rubrene $\text{C}_{42}\text{H}_{28}$ ($Cmca$, #64)** — Organic semiconductor with molecular $(100)$ mirror plane, $[100]$ $2$-fold axis, inversion center, and herringbone $b/c/a/n$ glide planes

---

## 🚀 Running Locally

The application code and presets are embedded in `index.html`; Three.js and OrbitControls load from jsDelivr and require Internet access. No build step or backend server is required:
1. Clone or download `index.html` from this repository.
2. Open `index.html` directly in any modern browser (Chrome, Firefox, Safari, Edge), or serve locally via:
   ```bash
   python3 -m http.server 8000

   ```


## Validation and interpretation

Run the regression suite with Node.js 20 or newer:

```bash
node tests/symmetry.test.cjs
```

The suite checks all seven presets, symmetry group closure, species-preserving atom permutations, animation endpoints, the offset ZnS screw axis, quartz site symmetry, CIF validation, strain behavior and known tensor dimensions.

Orbit letters in the UI are local identifiers, **not International Tables Wyckoff letters**. Symmetry detection uses the selected atomic tolerance for all species, including H, and a strict dimensionless cell-metric tolerance of `1e-6`. It does not refine an approximately symmetric experimental cell. Raise the atomic tolerance explicitly when appropriate; the detector no longer silently changes it.

A 180° rotation about an offset screw axis can take a corner to an edge midpoint. For ZnS, the axis `(1/2, 1/4, z)` maps `(0,0,0)` to `(1,1/2,0)` before the `c/2` translation. The inspector shows the actual fractional path for the selected atom pair. Symmetry operations differing by a lattice translation can have different axis locations; the inspector displays the chosen representative.

For metric-breaking strain, the animation applies the original rigid parent operation in the laboratory frame. The operation triplet and anchor then refer to the parent cell; path coordinates refer to the current cell. A displayed mismatch is a Cartesian distance in Å.

See [SYMMETRY_AUDIT.md](SYMMETRY_AUDIT.md) for the fixes and remaining scientific limitations.


### What a site rotation means

An axis through an atom does not mean the atom physically spins. A site-symmetry operation fixes that atom **and maps the entire periodic crystal onto equivalent same-element sites**. A stationary central sphere alone is insufficient evidence.

In quartz, each Si site has point group `2`: a 180° rotation leaves the Si fixed and exchanges the four nearby O atoms in two pairs. A 90° rotation also fixes the central point, but fails to reproduce the oxygen environment. Quartz O sites have point group `1`, so they have no nontrivial site rotation. In Site Symmetry mode, the animation and inspector use the same actual neighbor pair. The colored reference motif is an artificial orientation marker, not a physical atom or atomic spin.

For an optional independent reference check, install Python's `spglib` package and run:

```bash
python -m pip install spglib
python tests/reference-spglib.py
```

This compares every preset space-group number, every operation and every atomic site point group against spglib. It is a development check, not a browser runtime dependency.
