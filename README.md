# BathyModel v1.0

**© 2026 Katfish Heavy Industries**
Design: Katfish

A tool for converting bathymetric survey data (XYZ point clouds) into 3D-printable
STL models, with an optional recessed text label on the base.

Developed for generating physical harbor and waterway models from survey data.

---

## Scope and intended use

BathyModel is optimised for **facility-level surveys** — individual piers, berths,
slips, and similar bounded areas typically covered by surveys up to a few thousand
points at 3–15ft grid spacing.

It is **not intended for portfolio or city-scale datasets**. Large area surveys
(full waterfronts, harbor-wide coverage, 1m UTM grids spanning hundreds of acres)
will produce point counts and triangle counts that exceed what the pipeline handles
practically. For those datasets, consider splitting into facility-level tiles before
processing.

A practical upper limit is roughly 10,000–50,000 input points. Beyond that,
Delaunay triangulation time becomes the bottleneck and runtime grows significantly.

---

## License

BathyModel is released under the **GNU General Public License v3.0 (GPL-3.0)**.

You are free to use, modify, and distribute this software under the following conditions:

- Any distributed version of this software, or software that incorporates it, must also be released under GPL-3.0
- The source code must be made available to anyone who receives the software
- The original copyright and license notice must be preserved

This software is provided without warranty of any kind.

Full license text: https://www.gnu.org/licenses/gpl-3.0.en.html

A `LICENSE` file containing the full GPL-3.0 text is included in the GitHub repository.

---

## Attribution

BathyModel is built on the following open-source libraries:

- **NumPy** — fundamental array and numerical operations
  https://numpy.org — BSD License

- **SciPy** — Delaunay triangulation (via Qhull) and KD-tree nearest-neighbor search
  https://scipy.org — BSD License

- **Pillow (PIL Fork)** — font rendering for the recessed label
  https://python-pillow.org — HPND License

- **Python** — language runtime
  https://python.org — PSF License

- **DejaVu Sans Bold** — default label font (bundled with most Linux/Windows systems)
  https://dejavu-fonts.github.io — DejaVu Fonts License

The Delaunay triangulation and binary STL format handling follow well-established
computational geometry methods. The grid-based point-in-polygon and boundary-walking
algorithms are original implementations written for this project.

---

## Installation

### Executable (recommended)

Download `BathyModel.exe` and `README.md` from the GitHub releases page. Place them
in a folder of your choice and double-click the exe — no Python or dependencies required.

XYZ Files and STL Files folders are created automatically alongside the exe on first run.

### From source

Requires Python 3.8 or later.

```
pip install numpy scipy pillow
python BathyModel.py
```

---

## Input data format

BathyModel accepts XYZ files with no header row. Delimiters are auto-detected — tab,
comma, and space/whitespace are all supported.

```
X_coordinate    Y_coordinate    depth
6011026.98      2123259.99      33.4
6010441.97      2123394.97      23.5
...
```

- Coordinates may be in any projected system (State Plane, UTM, etc.)
- Depth values in feet, positive values = depth below surface
- Negative values (elevation convention) are also handled automatically
- **Optimised for 15ft shoal bias grid spacing** — the recommended input format
- 3ft grid files are supported but practical only for small individual berths or slips
- Recommended point count: 2,000–50,000 points for comfortable runtimes

---

## Using the GUI

1. Double-click `BathyModel.exe`
2. Click **Browse** next to XYZ file(s) and select your survey file
   - The file browser opens to `XYZ Files\` by default
   - The output path auto-fills to `STL Files\<filename>.stl`
   - Multiple files can be selected for batch processing
3. Enter a **Label** to be recessed into the base (e.g. `Pier 27 SF`)
   - Leave blank for no label
   - Auto-populated from filename if a pier name is detected
4. Enter the **Nominal Grid Spacing** in feet (e.g. `15`)
   - Optimised for 15x15 ft Shoal Bias files
   - Leave blank to auto-estimate from point spacing
5. Adjust **Parameters** if needed (see Parameters section below)
6. Click **Generate STL ▶**
7. Watch the log panel for progress — large files may take several minutes
8. The finished STL appears in `STL Files\`
9. Click **Open Output File** to open the STL in your default viewer

---

## Using the command line

```
python BathyModel.py input.xyz "Label Text" output.stl [nominal_spacing_ft]
```

Launching with no arguments opens the GUI:

```
python BathyModel.py
```

Examples:

```
python BathyModel.py "XYZ Files\Pier27.xyz" "Pier 27 SF" "STL Files\Pier27.stl" 15
python BathyModel.py "XYZ Files\Waterfront.xyz" "SF Waterfront" "STL Files\Waterfront.stl" 15
```

The output path argument is optional — if omitted, the file is saved to
`STL Files\<input_name>.stl` automatically.

---

## Parameters

These are set at the top of `BathyModel.py` and also editable in the GUI:

| Parameter | Default | Description |
|---|---|---|
| `BASE_MM` | 5.0 | Thickness of the flat base under the deepest point (mm) |
| `MAX_DIM_MM` | 200.0 | Maximum footprint dimension — model is scaled to fit (mm) |
| `Z_EXAG` | 3.0 | Vertical exaggeration — increase for flat/featureless surveys |
| `EDGE_FACTOR` | 3.0 | Edge-length filter multiplier (controls boundary triangle removal) |
| `RECESS_MM` | 1.5 | Depth of the recessed label carve-in (mm) |
| `CELL_MM` | 1.0 | Grid cell size for bottom face — increase for large files (mm) |
| `FONT_PT` | 14 | Label font size in points |

### Tips for large files

- Increase `CELL_MM` to `2.0` or `3.0` to reduce bottom face processing time
- Increase `Z_EXAG` to `8.0`–`10.0` for harbor surveys with low vertical relief
- The Delaunay triangulation step scales with point count and cannot be sped up
  further — a 1M+ point file will take several minutes regardless

---

## Output

BathyModel produces a binary STL file ready for slicing and 3D printing.

**Recommended print settings:**
- No supports required
- Print with seafloor side up (the label on the base faces down)
- 0.2mm layer height works well for typical pier-scale models
- The recessed label reads best with a contrasting filament or paint fill

**Slicer compatibility:**
The bottom face grid may produce minor T-junction gaps at the survey boundary.
These are cosmetic and can be repaired automatically with:
- PrusaSlicer → **Fix** (on import)
- Meshmixer → **Make Solid**
- Netfabb → **Repair**

---

## Troubleshooting

**"Boundary loops: 1" but model looks incomplete**
The edge filter may be bridging gaps between survey regions. Try reducing
`EDGE_FACTOR` from 3.0 to 2.0.

**Very flat model with no visible relief**
Increase `Z_EXAG`. Harbor surveys often need 8×–15× to show meaningful relief
at 200mm scale.

**Label is missing or partial**
The label placement algorithm finds the widest band inside the survey footprint.
On very irregular shapes the best available band may still be narrow. Try a
shorter label string.

**Application appears to hang at "Triangulating..."**
This is normal for large files. A 1M+ point survey will take several minutes
at this step. The process is alive — check CPU usage to confirm.

---

## File naming convention

Survey files from NOAA/port authorities typically encode useful metadata:

```
CDIM_20260204_POSF_Pier35_15x15_Depth_ShoalBias.xyz
      ↑ date   ↑ project  ↑ pier  ↑ grid   ↑ depth type
```

BathyModel uses the filename stem as the default output name, so keeping
descriptive filenames in `XYZ Files\` will produce descriptive STL names
in `STL Files\` automatically.
