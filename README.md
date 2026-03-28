# BathyModel v1.1

**© Katfish Heavy Industries**

Facility-scale bathymetric survey to 3D-printable STL converter.

Optimised for individual piers, berths, slips, and harbour approaches.
Not intended for portfolio or city-scale datasets.

---

## Contents

- [Overview](#overview)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Input File Formats](#input-file-formats)
- [Using the GUI](#using-the-gui)
- [Parameters](#parameters)
- [Multi-File Merge](#multi-file-merge)
- [Log Output Reference](#log-output-reference)
- [Troubleshooting](#troubleshooting)
- [File Naming Convention](#file-naming-convention)
- [Building from Source](#building-from-source)
- [License](#license)
- [Attributions](#attributions)

---

## Overview

BathyModel reads XYZ point-cloud sounding data, triangulates the seafloor
surface using a Delaunay mesh, extrudes a solid base, optionally recesses a
text label into the underside, and writes a binary STL ready for slicing and
printing.

**What it produces:**

```
┌─────────────────────────────────┐  ← seafloor surface (top)
│  shallows ↑     deeps ↓        │
│                                 │
├─────────────────────────────────┤  ← base (flat underside)
│       PIER 45  · SF BAY        │  ← recessed label
└─────────────────────────────────┘
```

Deeper areas print lower; shallower areas print higher. The 3× default
vertical exaggeration makes relief visible at print scale.

---

## Installation

### Distributed executable (recommended)

1. Copy `BathyModelV1.1.rar` within a folder of your choice, e.g. `F:\BathyModel\`
2. Uncompress the file, and run the .exe
3. Example Files are included from various sources in the \XYZ Files Folder.

No Python installation required.

### Running from source

**Requirements:** Python 3.9+

```
pip install numpy scipy pillow pypdf h5py
```

**Run:**

```
python BathyModel.py
```

---

## Quick Start

1. Launch `BathyModel.exe`
2. Click **Browse…** next to *.XYZ File(s) Location* and select your survey file
   — `.xyz`, `.txt`, `.pdf`, and `.bag` are all accepted
3. The Label and STL output path auto-fill from the filename
4. Adjust parameters if needed (defaults work well for 15 ft grid files)
5. Click **Generate STL ▶**
6. When complete, click **Open Output File** to preview in your STL viewer

---

## Input File Formats

BathyModel accepts four file types. All format detection is automatic — just
browse to the file and BathyModel handles the rest.

### Standard XYZ  (`.xyz`)

Tab- or space-delimited, three columns, no header:

```
X           Y           Z
1234567.89  456789.01   -12.5
1234568.12  456789.45   -13.1
```

Coordinates must be in a projected system (State Plane, UTM, etc.) where
X and Y are in feet or metres. Decimal-degree lat/lon will produce a
distorted or unusably small model — use NOAA format for those files.

Z values may be negative (depth below datum) or positive (depth below surface).
BathyModel detects the sign convention automatically.

### NOAA GeoImage TXT  (`.txt`)

Comma-delimited, three columns, no header:

```
latitude,longitude,depth
37.864517,-122.399565,14.78
37.854329,-122.399506,20.12
```

Depths are in metres, positive downward. BathyModel automatically detects the
NOAA format from the coordinate ranges, swaps columns to the expected order,
and projects decimal degrees to metres using an equirectangular approximation
centred on the survey (accurate to within 0.1% for facility-scale areas).

The log confirms detection and conversion:

```
⚑  NOAA GeoImage format detected (lat, lon, depth — comma-delimited)
✓  Converted: columns remapped (lon, lat, depth) and projected to metres
```

The Nominal Grid Spacing field is cleared automatically when a `.txt` file is
selected, since the spacing is in metres and auto-estimation is more reliable
than the foot-based default.

### NOAA GeoImage PDF  (`.pdf`)

NOAA distributes hydrographic surveys as GeoImage PDFs with the sounding
data embedded as an attachment (`H#####_Depths.txt`). BathyModel extracts
the attachment automatically — you do not need to open the PDF or save the
file separately.

```
⚑  PDF detected — extracting embedded sounding file ...
✓  Extracted: H13541_Depths.txt
✓  Converted: columns remapped (lon, lat, depth) and projected to metres
```

**Note:** PDF support requires the `pypdf` library, included in the distributed
exe. If running from source: `pip install pypdf`

**Note:** The sounding data embedded in a GeoImage PDF is a coarse reduced
dataset (typically 200m+ point spacing). For full-resolution data, download
the BAG file from NOAA's National Centers for Environmental Information:
https://www.ncei.noaa.gov/maps/bathymetry/

### BAG — Bathymetric Attributed Grid  (`.bag`)

NOAA's standard format for gridded multibeam survey products. BAG files
contain a regular elevation grid at the survey's native resolution and are
significantly higher quality than the sounding data embedded in GeoImage PDFs.

BathyModel reads the supergrid elevation layer, extracts the bounding box
from the embedded ISO 19115 metadata, and projects geographic coordinates
to metres where needed.

```
⚑  BAG format detected — Variable Resolution (VR) — supergrid only
    Full VR refinement support coming in v1.2
✓  Supergrid: 512 × 512 cells   approx. resolution: 29 m
   Valid soundings: 64,525
```

**Variable Resolution (VR) BAGs** — modern NOAA surveys use the VR variant
which stores data at multiple resolutions. BathyModel currently reads the
coarse supergrid layer only. Full native-resolution VR support is planned
for v1.2 and will dramatically improve model detail.

**Where to get BAG files** — search for your survey number at:
https://www.ncei.noaa.gov/maps/bathymetry/

**Note:** BAG support requires the `h5py` library, included in the distributed
exe. If running from source: `pip install h5py`

The Nominal Grid Spacing field is cleared automatically when a `.bag` file
is selected.

---

## Using the GUI

### Input section

**.XYZ File(s) Location** — browse to one or more survey files. Accepts
`.xyz`, `.txt`, `.pdf`, and `.bag`. Multi-select is supported for batch or
merge processing.

**Nominal Grid Spacing (ft)** — the survey's point spacing. Leave at 15 for
standard 15×15 ft Shoal Bias files. Leave blank to auto-estimate from the
data. Cleared automatically for NOAA `.txt`, `.pdf`, and `.bag` files since
those work in metres. Used to set the maximum triangle edge length filter;
too small causes holes at survey boundaries, too large allows phantom
triangles spanning gaps.

**Multi-file mode** — appears automatically when more than one file is
selected. See [Multi-File Merge](#multi-file-merge) below.

### Parameters section

| Parameter | Default | Description |
|---|---|---|
| Z-Height Exaggeration | 3.0 | Vertical scale multiplier. Increase for very flat surveys. |
| Base Padding Thickness (mm) | 5.0 | Thickness of the flat base below the deepest point. |

### Label section

| Parameter | Default | Description |
|---|---|---|
| Label Text | (from filename) | Text recessed into the bottom face. Leave blank for a flat base. |
| Font Size (pt) | 14 | Rendered using DejaVu Sans Bold. |
| Text Recess Depth (mm) | 1.5 | Depth of the recessed label into the base. |

### Output section

**.STL File(s) Output Location** — auto-fills with the survey filename in the
`STL Files\` folder. In batch mode this field is locked; files are named
automatically from each input filename.

### Buttons

**Generate STL ▶** — starts processing. The button is disabled during the run.
Multiple files run sequentially in a background thread; the GUI stays responsive.

**Open Output File** — opens the most recently completed STL in the system's
default viewer (e.g. Windows 3D Viewer, PrusaSlicer, Bambu Studio). Works
correctly in single, batch, and merge modes. Enabled after a successful run.

---

## Parameters

### Vertical exaggeration

The default 3× exaggeration suits most harbour and pier surveys at 200mm print
size. As a guide:

- Flat anchorage/approach surveys: try 6×–10×
- Active scour zones (e.g. under ferry ramps): 2×–3×
- Structural surveys with large relief: 1×–2×

### Edge factor

Hardcoded at 3.0× nominal spacing. Controls the triangle edge-length filter
that removes long boundary triangles at the survey perimeter. If the model has
large flat fins or skirts extending beyond the survey area, reduce `EDGE_FACTOR`
in the source to 2.0.

### Auto-thinning

If the survey point spacing is finer than the 1mm bottom-face cell grid,
BathyModel automatically thins the point cloud using a shoal-bias grid filter
(keeping the shallowest point per cell). This is logged when it occurs:

```
Auto-thinning: spacing 0.312mm < 1.0mm — thinning to 1.0mm grid (shoal bias)
12,450 → 3,891 points after thinning
```

---

## Multi-File Merge

When two or more files are selected, a **Multi-file mode** toggle appears in
the Input section with two options:

**Process as separate files** — the original batch behaviour. Each file is
processed independently and produces its own STL in `STL Files\`.

**Merge into single model** — all selected files are concatenated into a
single point cloud and processed as one job, producing one STL. The Label and
output path fields re-enable so you can set them for the merged result. The
output filename defaults to `<first_file>_merged.stl`.

### When to use merge

Merge is intended for complementary surveys of the same geographic area —
for example, a port authority survey of one side of a pier combined with a
NOAA survey of the channel alongside it. All file formats (`.xyz`, `.txt`,
`.pdf`, `.bag`) can be mixed in a merge — NOAA formats are projected to metres
before concatenation.

### Coordinate compatibility check

Before merging, BathyModel checks that all selected files have bounding boxes
that overlap or abut within a 5,000-unit tolerance. If any file appears to be
from a different geographic area or coordinate system, a warning dialog appears:

```
Coordinate compatibility check failed —
File 2 coordinates appear incompatible — X gap: 42,310, Y gap: 0.
Files may be in different coordinate systems or areas.
Proceed anyway?
```

You can proceed or cancel. Proceeding with incompatible coordinates will
produce a distorted or empty model.

### Merge log output

```
Merging 2 files:
  [1] AEW_20240628_POSF_FishermansWharf_3x3_depth.xyz  (122,292 pts)
  [2] CDIM_20250219_POSF_Pier45_15x15_ShoalBias.xyz  (2,879 pts)
  Coordinate check: compatible (X overlap: 1168 ft, Y overlap: 988 ft)
  Combined: 125,171 points

Loading <2-file merge> ...
  125,171 points   Z ∈ [0.01, 59.60] ft
  ...
```

---

## Log Output Reference

A typical successful run produces output like this:

```
Loading CDIM_20260204_POSF_Pier35_15x15_Depth_ShoalBias.xyz ...
  2,879 points   Z ∈ [0.01, 29.30] ft
  Footprint : 121.0 × 200.0 mm
  Scale     : 1 mm = 47.2 ft  (1 : 566)
  Z range   : 5.0 – 8.7 mm  (×3.0 exaggeration)
Triangulating 2,879 points ...
  7,241 → 6,903 triangles after filter
  Finding boundary ...
  Boundary loops: 2
Building bottom face ...
  Label 'Pier 35 · SF Bay' — 128×14 cells at 14pt
  Placed 1,247/1,251 pixels inside footprint
  6,280 bottom cells
Writing STL Files\CDIM_20260204_POSF_Pier35_15x15_Depth_ShoalBias.stl ...

✓  CDIM_20260204_POSF_Pier35_15x15_Depth_ShoalBias.stl  (28,902 triangles)
```

**Scale line** — reports the real-world distance per mm of model and the
dimensionless ratio. Unit is `ft` for standard XYZ files and `m` for
NOAA-projected and BAG data.

**Boundary loops** — the number of closed boundary edges found. A single
convex survey typically has 1 loop. Pier surveys with structural gaps (e.g.
a pier cutting through the survey) may have 2 or more.

**Track-line survey detected** — shown when the point spacing analysis finds
a bimodal edge-length distribution, indicating a multibeam track-line survey.
BathyModel uses the cross-track spacing for the edge filter in this case:

```
  Track-line survey detected — using cross-track spacing: 4.21mm
  (along-track: 0.38mm, cross-track: 4.21mm)
```

**BAG log example:**

```
Loading H13541_MB_VR_MLLW_1of1.bag ...
⚑  BAG format detected — Variable Resolution (VR) — supergrid only
    Full VR refinement support coming in v1.2
✓  Supergrid: 512 × 512 cells   approx. resolution: 29 m
   Valid soundings: 64,525
  64,525 points   Z ∈ [-56.42, -1.10] m
  Footprint : 145.0 × 200.0 mm
  Scale     : 1 mm = 75.4 m  (1 : 75,400)
  Z range   : 5.0 – 13.2 mm  (×3.0 exaggeration)
```

---

## Troubleshooting

**Holes or missing areas in the model**

The edge-length filter is trimming too aggressively. Try increasing the
Nominal Grid Spacing field, or leave it blank to let BathyModel auto-estimate.
If the survey is a track-line file, the auto-detection should handle this, but
a very irregular track pattern may need manual tuning.

**Large triangular fins extending beyond the survey boundary**

The edge filter is too permissive. Reduce Nominal Grid Spacing, or if running
from source, reduce `EDGE_FACTOR` from 3.0 to 2.0.

**Very flat model with no visible relief**

Increase Z-Height Exaggeration. Harbour surveys often need 6×–15× to show
meaningful relief at 200mm scale. BAG supergrid data at 29m resolution will
also appear smoother than a high-resolution pier survey — this is a data
resolution limitation, not a BathyModel issue.

**Label is missing or only partially visible**

The label placement algorithm finds the widest horizontal band inside the
survey footprint. On very irregular shapes the best available band may be
narrow. Try a shorter label string or reduce Font Size.

**Application appears to hang at "Triangulating..."**

This is normal for large files. A survey with 100,000+ points will take
several minutes at this step. The application is alive — check CPU usage to
confirm. The GUI remains responsive throughout.

**Open Output File does nothing after a batch run**

Fixed in v1.1. The button now tracks the last successfully written file
across all modes (single, batch, and merge).

**NOAA PDF: "No .txt attachment found"**

The PDF may be a chart rather than a GeoImage survey PDF. GeoImage PDFs
contain an embedded `H#####_Depths.txt` file. Verify by opening the PDF in
Adobe Acrobat and checking the Attachments panel (View → Navigation Panes →
Attachments).

**BAG: "Could not parse bounding box from BAG metadata XML"**

The metadata XML in the BAG file uses an unusual structure. Please report the
survey number so the parser can be updated.

**BAG: h5py ImportError**

Run `pip install h5py` then rebuild the exe with PyInstaller.

**Coordinate compatibility warning on merge**

The selected files may be from different areas or in different coordinate
systems. Verify that all files cover the same geographic region. NOAA `.txt`,
`.pdf`, and `.bag` files are projected to metres before the compatibility
check, so they can be safely merged with other metre-based files.

---

## File Naming Convention

Survey files from NOAA and port authorities typically encode useful metadata:

```
CDIM_20260204_POSF_Pier35_15x15_Depth_ShoalBias.xyz
      ↑ date   ↑ project  ↑ pier  ↑ grid   ↑ depth type
```

BathyModel uses the filename stem as the default label and output name, so
keeping descriptive filenames in `XYZ Files\` will produce descriptive STL
names in `STL Files\` automatically.

---

## Building from Source

**Dependencies:**

```
pip install numpy scipy pillow pypdf h5py
```

**PyInstaller (single-file Windows exe):**

Run from your `Scripts\` folder in Command Prompt:

```
pyinstaller --onefile --windowed --icon=BathyModel.ico --name=BathyModel BathyModel.py
```

The exe will be at `dist\BathyModel.exe`. Copy it to your release folder
alongside the README.

---

## License

BathyModel — Copyright (C) 2026 Katfish Heavy Industries

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, version 3.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

https://www.gnu.org/licenses/gpl-3.0.en.html

---

## Attributions

BathyModel is built on the following open-source libraries:

**NumPy** — Harris et al., *Nature* 585, 357–362 (2020)
BSD 3-Clause License — https://numpy.org

**SciPy** — Virtanen et al., *Nature Methods* 17, 261–272 (2020)
BSD 3-Clause License — https://scipy.org

**Pillow** — Jeff Widman et al.
HPND License — https://python-pillow.org

**pypdf** — Mathieu Fenniak et al.
BSD 3-Clause License — https://pypdf.readthedocs.io

**h5py** — Andrew Collette et al.
BSD 3-Clause License — https://www.h5py.org

**DejaVu Fonts** — Bitstream Vera Fonts / Arev Fonts
DejaVu Fonts License — https://dejavu-fonts.github.io
