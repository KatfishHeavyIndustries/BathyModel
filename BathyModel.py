#!/usr/bin/env python3
"""
BathyModel v1.1
Facility-scale bathymetric survey to 3D-printable STL converter.

Copyright (C) 2026 Katfish Heavy Industries

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <https://www.gnu.org/licenses/>.

Design: Katfish

Optimised for facility-level surveys (individual piers, berths, slips).
Not intended for portfolio or city-scale datasets.

Usage (CLI):
    python BathyModel.py input.xyz "Label Text" output.stl [nominal_spacing_ft]

Usage (GUI):
    python BathyModel.py

Dependencies:
    pip install numpy scipy pillow pypdf h5py
"""

__version__ = "1.1"
__author__  = "Katfish Heavy Industries"
__license__ = "GPL-3.0"

# BathyModel — © Katfish Heavy Industries
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3.
# https://www.gnu.org/licenses/gpl-3.0.en.html

import sys
import os
import re
import struct
import threading
import numpy as np
from scipy.spatial import Delaunay, KDTree
from collections import defaultdict
from PIL import Image, ImageDraw, ImageFont

# ── embedded icon (base64) ───────────────────────────────────────────────────
_ICON_B64 = "AAABAAEAEBAAAAAAIAB5AwAAFgAAAIlQTkcNChoKAAAADUlIRFIAAAAQAAAAEAgGAAAAH/P/YQAAA0BJREFUeJxdk1FMW2Ucxc/3fe0tt6XdljFHTBSLsAzHOgtOYE0fZLIsPPC23gfiWEhcyMgeTZxmybKYTHxtXGcXEvukK5psRs3MIspkJbJ0DwIPQtc61pYsA1JvKXCv/b6/DwLBnafzcH4nOQ+HYUvJZFJEIhEJABcvXDj9akNDj67rTUIImKaZyReL96LR6N0XszswAIyMjISTX32depRO0/LzZbIti0zTpJWVFXo4PU3fjn2Tunzpcng3s2M+u3bt7M0v4iqfz9PS0pIqFArVhYUFSUQyl8tV19bW1NzsHI0lk+qTq1fPbrPcMAw5ODjY9e6pUwlJCjMzMzIej7NMJiNisRi3bZsnEgmRyWTYne/uyKbmZrzT3Z0YPn++yzAMyYnI8fbx49FgWxv27N2risWiKJVKCIVCUEpBCAFd1+F0OiGEEC6XS50IhfBme3uUiBx8eGjoZCAQaAegXJom0uk06uvrIYSAZVmQUkJKCSICANi2LQCoo4FA+/DQ0Enu9zf2+v1+AkCapmFychKdnZ0AACEENE1DbW3tTsGW6DW/nxr8jb281udtdHs8DAArl8soFAro6OgAAKyursIwDIyOjsLr9UJKCcYYADCPx8N8Pm8j313rdDrBOYdpmgAAXdfR39+PYDAI27bB+f/iAABumuVspVL5b0KNjjeOHMGvExOQBGi6G319fTgWDGLTsgDOUVUEAFSpVMg0y1n+JJf9cTGXZQAYXzfRHQ7hj0dpCAa4/lkHlA21vgYHZ3BULexzEACwxVyWPf0r9wP//MaNn6dn/0yPl8C/L3vkK2+FYcKJsfllpA/34npe4LFyw5KElQPNiC175C+rxH+fm09HY7FxBsbQ8P6VLtHzXqrEdXpZk2qzsi5K5IAkYH/dfpSKeRw6uA9P/96QVV8drzGfM2s8eaJ48+Mpx5lbt8RYJDIV8PkGeKD7S1V3WLj1TXJXLaWqVeZ2ADUHXyLT6eLuA3sEe/aExMOfBnLxj6bObP9h2xw9dynccn0i1Xp7kYL3N6ntQZWCv9kUnNig1tuL1BK7n2o590F4N8N2XVJg66KHBj48Ta+39jDd18QByA0zwx7P3ptPfHr3xey/lGST0byGBlMAAAAASUVORK5CYII="
# ─────────────────────────────────────────────────────────────────────────────

# ── parameters ────────────────────────────────────────────────────────────────
BASE_MM     = 5.0    # flat base thickness under deepest survey point (mm)
MAX_DIM_MM  = 200.0  # maximum XY footprint dimension (mm)
Z_EXAG      = 3.0    # vertical exaggeration
EDGE_FACTOR = 3.0    # edge-length filter multiplier × nominal spacing
RECESS_MM   = 1.5    # label recess depth into the base (mm)
CELL_MM     = 1.0    # bottom-face grid cell size (mm)
FONT_PT     = 14     # label font size (pt) — DejaVu Sans Bold
# ─────────────────────────────────────────────────────────────────────────────


# ══════════════════════════════════════════════════════════════════════════════
#  CORE PIPELINE
# ══════════════════════════════════════════════════════════════════════════════

def write_stl(tri_array, path):
    """Write binary STL from (N,3,3) float32 array. Header exactly 80 bytes."""
    assert tri_array.ndim == 3 and tri_array.shape[1:] == (3, 3)
    n_tris = len(tri_array)
    p0, p1, p2 = tri_array[:, 0], tri_array[:, 1], tri_array[:, 2]
    norms = np.cross(p1 - p0, p2 - p0).astype(np.float32)
    mag   = np.linalg.norm(norms, axis=1, keepdims=True)
    mag[mag == 0] = 1.0
    norms /= mag
    record = np.zeros(n_tris, dtype=[
        ('normal', np.float32, (3,)),
        ('v0',     np.float32, (3,)),
        ('v1',     np.float32, (3,)),
        ('v2',     np.float32, (3,)),
        ('attr',   np.uint16),
    ])
    record['normal'] = norms
    record['v0']     = p0
    record['v1']     = p1
    record['v2']     = p2
    hdr = b"BathyModel output".ljust(80, b'\x00')[:80]
    with open(path, 'wb') as f:
        f.write(hdr)
        f.write(struct.pack('<I', n_tris))
        f.write(record.tobytes())


def order_boundary(b_edges_arr):
    """
    Chain boundary edges into ordered loops.
    Handles fragmented or branching edges (e.g. from incompatible merged
    surveys) by capping walk length and discarding degenerate fragments.
    """
    adj = defaultdict(list)
    for a, b in b_edges_arr:
        adj[int(a)].append(int(b))
        adj[int(b)].append(int(a))

    max_steps = len(b_edges_arr) + 1   # safety cap — can never need more
    remaining = set(map(int, b_edges_arr.ravel().tolist()))
    loops = []
    while remaining:
        start = next(iter(remaining))
        loop, prev, curr = [start], None, start
        visited = {start}
        steps = 0
        while steps < max_steps:
            # prefer unvisited neighbours; fall back to start to close the loop
            candidates = [nb for nb in adj[curr] if nb != prev and nb not in visited]
            if not candidates:
                # try to close back to start
                if start in adj[curr] and curr != start:
                    loop.append(start)
                break
            nxt = candidates[0]
            if nxt == start:
                break
            loop.append(nxt)
            visited.add(nxt)
            prev, curr = curr, nxt
            steps += 1
        for v in loop:
            remaining.discard(v)
        if len(loop) >= 3:
            loops.append(loop)
    return loops


def find_boundary_edges(simplices):
    e01 = np.sort(simplices[:, [0, 1]], axis=1)
    e12 = np.sort(simplices[:, [1, 2]], axis=1)
    e20 = np.sort(simplices[:, [2, 0]], axis=1)
    all_edges = np.vstack([e01, e12, e20])
    n_verts   = simplices.max() + 1
    encoded   = all_edges[:, 0].astype(np.int64) * n_verts + all_edges[:, 1]
    uniq, cnt = np.unique(encoded, return_counts=True)
    boundary  = uniq[cnt == 1]
    b_v0 = (boundary // n_verts).astype(np.int32)
    b_v1 = (boundary  % n_verts).astype(np.int32)
    return np.column_stack([b_v0, b_v1])


def batch_pip(px, py, poly_segments):
    inside = np.zeros(len(px), dtype=bool)
    for xi, yi, xj, yj in poly_segments:
        cross_y  = (yi > py) != (yj > py)
        x_inters = (xj - xi) * (py - yi) / (yj - yi + 1e-15) + xi
        inside  ^= cross_y & (px < x_inters)
    return inside


def build_poly_segments(b_loops, xn, yn):
    segs = []
    for loop in b_loops:
        n = len(loop)
        for k in range(n):
            i, j = loop[k], loop[(k + 1) % n]
            segs.append((xn[i], yn[i], xn[j], yn[j]))
    return segs


def render_text_bitmap(text, font_pt, cell_mm):
    font = ImageFont.load_default()
    for fp in [
        "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
        "/usr/share/fonts/truetype/liberation/LiberationSans-Bold.ttf",
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
    ]:
        try:
            font = ImageFont.truetype(fp, font_pt)
            break
        except OSError:
            pass
    probe = Image.new("L", (4000, 400))
    bbox  = ImageDraw.Draw(probe).textbbox((0, 0), text, font=font)
    tw, th = bbox[2] - bbox[0], bbox[3] - bbox[1]
    img = Image.new("L", (tw + 4, th + 4), 0)
    ImageDraw.Draw(img).text((2, 2), text, font=font, fill=255)
    pix = (np.array(img) >= 128).astype(bool)
    return pix, pix.shape[1], pix.shape[0]


def widest_run(sorted_cols):
    if not sorted_cols:
        return (0, 0)
    best_start, best_len = sorted_cols[0], 1
    cur_start,  cur_len  = sorted_cols[0], 1
    for i in range(1, len(sorted_cols)):
        if sorted_cols[i] == sorted_cols[i - 1] + 1:
            cur_len += 1
        else:
            if cur_len > best_len:
                best_start, best_len = cur_start, cur_len
            cur_start, cur_len = sorted_cols[i], 1
    if cur_len > best_len:
        best_start, best_len = cur_start, cur_len
    return (best_start, best_len)


def place_label(pix, w_cells, h_cells, interior_cols_by_row, n_rows):
    row_run = {row: widest_run(cols) for row, cols in interior_cols_by_row.items()}
    best_score, best_top_row, best_anchor = -1, 0, (0, 0)
    for top_row in range(n_rows - h_cells + 1):
        band = range(top_row, top_row + h_cells)
        if not all(r in row_run for r in band):
            continue
        min_len = min(row_run[r][1] for r in band)
        if min_len > best_score:
            best_score   = min_len
            best_top_row = top_row
            narrowest    = min(band, key=lambda r: row_run[r][1])
            best_anchor  = row_run[narrowest]
    if best_score < 0:
        if not row_run:
            return set()
        best_top_row = max(row_run, key=lambda r: row_run[r][1])
        best_anchor  = row_run[best_top_row]
    run_start, run_width = best_anchor
    col_offset = run_start + max(0, (run_width - w_cells) // 2)
    row_offset = best_top_row
    cells = set()
    for r, c in zip(*np.where(pix)):
        cells.add((col_offset + int(c), row_offset + int(r)))
    return cells


def build_bottom_face(interior_arr, text_cells, cell_mm, recess_mm):
    rd   = recess_mm
    cols = interior_arr[:, 0].astype(float)
    rows = interior_arr[:, 1].astype(float)
    x0   = cols * cell_mm;  x1 = x0 + cell_mm
    y0   = rows * cell_mm;  y1 = y0 + cell_mm
    is_recess = np.array([(int(c), int(r)) in text_cells
                          for c, r in zip(cols, rows)], dtype=bool)
    flat_mask = ~is_recess
    tris = []
    if flat_mask.any():
        fx0 = x0[flat_mask]; fx1 = x1[flat_mask]
        fy0 = y0[flat_mask]; fy1 = y1[flat_mask]
        z   = np.zeros(flat_mask.sum(), dtype=np.float32)
        t1  = np.stack([np.column_stack([fx0, fy0, z]),
                        np.column_stack([fx0, fy1, z]),
                        np.column_stack([fx1, fy0, z])], axis=1)
        t2  = np.stack([np.column_stack([fx1, fy0, z]),
                        np.column_stack([fx0, fy1, z]),
                        np.column_stack([fx1, fy1, z])], axis=1)
        tris.append(np.vstack([t1, t2]))
    recess_tris = []
    for idx in np.where(is_recess)[0]:
        cx, cy = int(cols[idx]), int(rows[idx])
        rx0, rx1 = x0[idx], x1[idx]
        ry0, ry1 = y0[idx], y1[idx]
        recess_tris += [
            [[rx0, ry0, rd], [rx0, ry1, rd], [rx1, ry0, rd]],
            [[rx1, ry0, rd], [rx0, ry1, rd], [rx1, ry1, rd]],
        ]
        if (cx, cy - 1) not in text_cells:
            recess_tris += [
                [[rx0, ry0, 0.], [rx1, ry0, 0.], [rx1, ry0, rd]],
                [[rx0, ry0, 0.], [rx1, ry0, rd], [rx0, ry0, rd]],
            ]
        if (cx, cy + 1) not in text_cells:
            recess_tris += [
                [[rx0, ry1, 0.], [rx1, ry1, rd], [rx1, ry1, 0.]],
                [[rx0, ry1, 0.], [rx0, ry1, rd], [rx1, ry1, rd]],
            ]
        if (cx - 1, cy) not in text_cells:
            recess_tris += [
                [[rx0, ry0, 0.], [rx0, ry0, rd], [rx0, ry1, rd]],
                [[rx0, ry0, 0.], [rx0, ry1, rd], [rx0, ry1, 0.]],
            ]
        if (cx + 1, cy) not in text_cells:
            recess_tris += [
                [[rx1, ry0, 0.], [rx1, ry1, rd], [rx1, ry0, rd]],
                [[rx1, ry0, 0.], [rx1, ry1, 0.], [rx1, ry1, rd]],
            ]
    if recess_tris:
        tris.append(np.array(recess_tris, dtype=np.float32))
    return np.vstack(tris).astype(np.float32) if tris else np.zeros((0, 3, 3), dtype=np.float32)


def _extract_noaa_from_pdf(pdf_path):
    """
    Extract the NOAA sounding data file from a GeoImage PDF attachment.
    Looks for an embedded file whose name contains 'depth' (case-insensitive)
    and ends with '.txt'.  Falls back to the first .txt attachment found.
    Returns the extracted content as a bytes object, or raises ValueError
    if no suitable attachment is found.
    """
    try:
        from pypdf import PdfReader
    except ImportError:
        raise ImportError(
            "pypdf is required to open NOAA GeoImage PDFs.\n"
            "Install it with:  pip install pypdf"
        )
    reader = PdfReader(pdf_path)
    attachments = reader.attachments          # dict: name -> [bytes, ...]
    if not attachments:
        raise ValueError("No embedded file attachments found in PDF.")

    # Prefer a .txt file with 'depth' in the name
    candidate = None
    for name, contents in attachments.items():
        if name.lower().endswith('.txt'):
            if 'depth' in name.lower():
                candidate = (name, contents[0])
                break
            if candidate is None:
                candidate = (name, contents[0])   # fallback: first .txt

    if candidate is None:
        raise ValueError(
            "No .txt attachment found in PDF. "
            "Expected a file like 'H12345_Depths.txt'."
        )
    return candidate   # (filename, bytes)


def _detect_noaa_format(filepath):
    """
    Returns True if the file appears to be a NOAA GeoImage TXT sounding file:
      - comma-delimited
      - 3 columns
      - column 0 in plausible CONUS latitude range  (24–50 °N)
      - column 1 in plausible CONUS longitude range (65–125 °W, stored negative)
    """
    try:
        with open(filepath) as f:
            first = f.readline().strip()
        parts = first.split(',')
        if len(parts) != 3:
            return False
        col0, col1 = float(parts[0]), float(parts[1])
        return 24.0 <= col0 <= 50.0 and -125.0 <= col1 <= -65.0
    except Exception:
        return False


def _project_lonlat_to_meters(lon, lat):
    """
    Equirectangular projection to metres, centred on the survey.
    Accurate enough for facility-scale surveys (<50 km extent).
    Returns (x_m, y_m) as numpy arrays.
    """
    R    = 6_371_000.0                      # Earth radius in metres
    lat0 = np.radians(np.mean(lat))
    x    = np.radians(lon - np.mean(lon)) * R * np.cos(lat0)
    y    = np.radians(lat - np.mean(lat)) * R
    return x, y


def _read_bag_supergrid(filepath, log=lambda _: None):
    """
    Read a BAG (Bathymetric Attributed Grid) file and return (x, y, z, coord_unit).
    Reads the supergrid elevation layer only — valid for both standard and VR BAGs.
    Full VR refinement layer support is planned for v1.2.

    Requires h5py:  pip install h5py
    """
    try:
        import h5py
    except ImportError:
        raise ImportError(
            "h5py is required to open BAG files.\n"
            "Install it with:  pip install h5py"
        )
    import xml.etree.ElementTree as ET

    with h5py.File(filepath, 'r') as f:
        elev   = f['/BAG_root/elevation'][:]          # (rows, cols) float32
        meta_bytes = f['/BAG_root/metadata'][()]
        is_vr  = 'varres_refinements' in f['/BAG_root']

    # metadata can come back as bytes, ndarray of bytes, or ndarray of objects
    if hasattr(meta_bytes, 'decode'):
        meta_xml = meta_bytes.decode('utf-8', errors='replace')
    else:
        # numpy array — flatten and join
        flat = meta_bytes.ravel()
        if flat.dtype.kind in ('S', 'O'):
            meta_xml = b''.join(
                v if isinstance(v, bytes) else v.encode('utf-8')
                for v in flat
            ).decode('utf-8', errors='replace')
        else:
            meta_xml = bytes(flat.tolist()).decode('utf-8', errors='replace')

    # ── Parse bounding box from ISO 19115 XML ───────────────────────────────
    import re

    # BAG metadata XML often has an XML declaration with encoding attribute
    # and heavy namespace usage — strip both carefully before parsing
    # Remove XML declaration entirely (ET chokes on encoding= in strings)
    meta_xml = re.sub(r'<[?]xml[^?]*[?]>', '', meta_xml).strip()
    # Remove namespace declarations (xmlns:xx="..." and xmlns="...")
    meta_xml = re.sub(r'\s+xmlns(?::\w+)?="[^"]*"', '', meta_xml)
    # Remove namespace prefixes from tags: <gmd:foo> → <foo>, </gmd:foo> → </foo>
    meta_xml = re.sub(r'<(\w+):', '<', meta_xml)
    meta_xml = re.sub(r'</(\w+):', '</', meta_xml)

    try:
        root = ET.fromstring(meta_xml)
    except ET.ParseError as e:
        # Last resort: extract bounding box with regex directly from raw XML
        log(f"    XML parse warning: {e} — falling back to regex extraction")
        root = None

    def find_text(tag):
        if root is not None:
            el = root.find('.//' + tag)
            if el is not None:
                # value may be direct text or in a child <Decimal> element
                dec = el.find('Decimal')
                txt = dec.text if dec is not None else el.text
                if txt:
                    return float(txt.strip())
        # Regex fallback — handles both direct value and <Decimal> child
        # Direct:  <westBoundLongitude>-122.5</westBoundLongitude>
        # Wrapped: <westBoundLongitude><Decimal>-122.5</Decimal></westBoundLongitude>
        m = re.search(
            rf'<{tag}[^>]*>[ \t\n\r]*(?:<Decimal>)?[ \t\n\r]*([-0-9.]+)[ \t\n\r]*(?:</Decimal>)?[ \t\n\r]*</{tag}>',
            meta_xml)
        if m:
            return float(m.group(1))
        return None

    west  = find_text('westBoundLongitude')
    east  = find_text('eastBoundLongitude')
    south = find_text('southBoundLatitude')
    north = find_text('northBoundLatitude')

    if any(v is None for v in [west, east, south, north]):
        # Log a snippet of the cleaned XML to help diagnose tag names
        snippet = meta_xml[:2000].replace('\n', ' ').replace('\r', '')
        log(f"    XML snippet: {snippet}")
        raise ValueError(
            "Could not parse bounding box from BAG metadata XML. "
            f"Found: W={west} E={east} S={south} N={north}"
        )

    rows, cols = elev.shape
    nodata     = 1_000_000.0

    # Build X/Y grid from bounding box
    lons = np.linspace(west, east, cols)
    lats = np.linspace(south, north, rows)
    lon_grid, lat_grid = np.meshgrid(lons, lats)

    # Mask nodata
    mask = elev != nodata
    if not np.any(mask):
        raise ValueError("BAG elevation grid contains no valid soundings.")

    lon_pts = lon_grid[mask]
    lat_pts = lat_grid[mask]
    z_pts   = elev[mask].astype(float)

    # Determine if geographic (lat/lon) or already projected
    if abs(west) <= 180 and abs(south) <= 90:
        x, y = _project_lonlat_to_meters(lon_pts, lat_pts)
        coord_unit = "m"
    else:
        x, y       = lon_pts, lat_pts
        coord_unit = "m"   # projected BAGs are always metric

    variant = "Variable Resolution (VR) — supergrid only" if is_vr else "Standard"
    log(f"⚑  BAG format detected — {variant}")
    if is_vr:
        log("    Full VR refinement support coming in v1.2")
    res_x = (east - west) / cols
    res_y = (north - south) / rows
    # Approximate resolution in metres at survey centre
    R = 6_371_000.0
    res_m = np.radians(max(res_x, res_y)) * R
    log(f"✓  Supergrid: {rows} × {cols} cells   approx. resolution: {res_m:.0f} m")
    log(f"   Valid soundings: {int(np.sum(mask)):,}")

    return x, y, z_pts, coord_unit


def _load_xyz_array(filepath, log=lambda _: None):
    """
    Load an XYZ survey file into numpy arrays, handling NOAA GeoImage TXT
    format and NOAA GeoImage PDF (with embedded .txt attachment) transparently.
    Returns (x, y, z, coord_unit).
    """
    import io, tempfile, os

    # ── BAG: supergrid read ───────────────────────────────────────────────────
    if filepath.lower().endswith('.bag'):
        return _read_bag_supergrid(filepath, log=log)

    # ── PDF: extract embedded Depths.txt first ────────────────────────────────
    txt_data = None
    if filepath.lower().endswith('.pdf'):
        log("⚑  PDF detected — extracting embedded sounding file ...")
        att_name, att_bytes = _extract_noaa_from_pdf(filepath)
        log(f"✓  Extracted: {att_name}")
        txt_data = att_bytes

    if txt_data is not None:
        # Parse from in-memory bytes
        text = txt_data.decode('utf-8', errors='replace')
        pts  = np.loadtxt(io.StringIO(text), delimiter=',')
        lat, lon, depth = pts[:, 0], pts[:, 1], pts[:, 2]
        x, y = _project_lonlat_to_meters(lon, lat)
        z    = depth
        log("✓  Converted: columns remapped (lon, lat, depth) and projected to metres")
        return x, y, z, "m"

    # ── Plain .txt or .xyz ────────────────────────────────────────────────────
    noaa_fmt = _detect_noaa_format(filepath)
    if noaa_fmt:
        log("⚑  NOAA GeoImage format detected (lat, lon, depth — comma-delimited)")
        pts = np.loadtxt(filepath, delimiter=',')
        lat, lon, depth = pts[:, 0], pts[:, 1], pts[:, 2]
        x, y = _project_lonlat_to_meters(lon, lat)
        z    = depth
        log("✓  Converted: columns remapped (lon, lat, depth) and projected to metres")
        coord_unit = "m"
    else:
        with open(filepath) as f:
            first_line = f.readline().strip()
        delim = '\t' if '\t' in first_line else (',' if ',' in first_line else None)
        pts   = np.loadtxt(filepath, delimiter=delim)
        x, y, z    = pts[:, 0], pts[:, 1], pts[:, 2]
        coord_unit = "ft"
    return x, y, z, coord_unit


def _check_coordinate_compatibility(file_list, tolerance=5000):
    """
    Load bounding boxes for each file and verify they overlap or abut within
    *tolerance* coordinate units.  Returns (compatible: bool, message: str).
    Operates on post-projection coordinates so NOAA .txt files are safe to mix.
    """
    bounds = []
    for f in file_list:
        x, y, z, _ = _load_xyz_array(f)
        bounds.append((x.min(), x.max(), y.min(), y.max()))

    ref = bounds[0]
    for i, b in enumerate(bounds[1:], 1):
        x_gap = max(0.0, max(ref[0], b[0]) - min(ref[1], b[1]))
        y_gap = max(0.0, max(ref[2], b[2]) - min(ref[3], b[3]))
        if x_gap > tolerance or y_gap > tolerance:
            return False, (
                f"File {i+1} coordinates appear incompatible — "
                f"X gap: {x_gap:.0f}, Y gap: {y_gap:.0f}. "
                f"Files may be in different coordinate systems or areas."
            )
    return True, "compatible"


def run_pipeline(xyz_file, label, out_file, nom_ft=None,
                 base_mm=None, z_exag=None, recess_mm=None,
                 cell_mm=None, font_pt=None,
                 xyz_points=None, coord_unit=None,
                 log=print):
    """
    Core conversion pipeline. log= accepts any callable (print, GUI logger, etc.)
    All parameter overrides default to module-level constants if not supplied.

    xyz_points : optional (x, y, z) tuple of pre-loaded numpy arrays.
                 When supplied, xyz_file is used only for display in the log
                 and file loading is skipped entirely.
    coord_unit : "ft" or "m" — required when xyz_points is provided.
    """
    _base    = base_mm   if base_mm   is not None else BASE_MM
    _zexag   = z_exag    if z_exag    is not None else Z_EXAG
    _recess  = recess_mm if recess_mm is not None else RECESS_MM
    _cell    = cell_mm   if cell_mm   is not None else CELL_MM
    _font    = font_pt   if font_pt   is not None else FONT_PT

    # ── Load ──────────────────────────────────────────────────────────────────────────────
    log(f"Loading {xyz_file} ...")
    if xyz_points is not None:
        # Pre-loaded data supplied by merge path — skip file I/O
        x, y, z = xyz_points
        # coord_unit already set by caller
    else:
        x, y, z, coord_unit = _load_xyz_array(xyz_file, log=log)
    log(f"  {len(x):,} points   Z ∈ [{z.min():.2f}, {z.max():.2f}] {coord_unit}")

    xy_span = max(x.max() - x.min(), y.max() - y.min())
    scale   = MAX_DIM_MM / xy_span
    xn = (x - x.min()) * scale
    yn = (y - y.min()) * scale

    if np.median(z) >= 0:
        z_height = (z.max() - z) * scale * _zexag
    else:
        z_height = (z - z.min()) * scale * _zexag
    zn = z_height + _base

    log(f"  Footprint : {xn.max():.1f} × {yn.max():.1f} mm")
    _real_per_mm = 1.0 / scale
    _ratio       = _real_per_mm * (1000.0 if coord_unit == "m" else 304.8)
    log(f"  Scale     : 1 mm = {_real_per_mm:.1f} {coord_unit}  (1 : {int(_ratio):,})")
    log(f"  Z range   : {_base:.1f} – {zn.max():.1f} mm  (×{_zexag} exaggeration)")

    # ── Thinning ──────────────────────────────────────────────────────────────
    xy_tmp   = np.column_stack([xn, yn])
    tree_tmp = KDTree(xy_tmp[:min(5000, len(xy_tmp))])
    d_tmp, _ = tree_tmp.query(xy_tmp[:min(5000, len(xy_tmp))], k=2)
    est_spacing = float(np.median(d_tmp[:, 1]))

    if est_spacing < _cell:
        thin_cell = _cell
        log(f"  Auto-thinning: spacing {est_spacing:.3f}mm < {_cell}mm — "
            f"thinning to {thin_cell:.1f}mm grid (shoal bias)")
    else:
        thin_cell = max(est_spacing / 4.0, 1e-6)

    gx  = np.floor(xn / thin_cell).astype(np.int64)
    gy  = np.floor(yn / thin_cell).astype(np.int64)
    key = gx * 10_000_000 + gy
    order      = np.argsort(key * 1_000_000 - zn)
    key_sorted = key[order]
    _, first   = np.unique(key_sorted, return_index=True)
    keep       = order[first]
    if len(keep) < len(xn):
        log(f"  {len(xn):,} → {len(keep):,} points after thinning")
        xn, yn, zn = xn[keep], yn[keep], zn[keep]

    # ── Delaunay ──────────────────────────────────────────────────────────────
    xy = np.column_stack([xn, yn])
    log(f"Triangulating {len(xn):,} points ...")
    tri_obj = Delaunay(xy)

    s   = tri_obj.simplices
    tp0, tp1, tp2 = xy[s[:, 0]], xy[s[:, 1]], xy[s[:, 2]]
    e01 = np.linalg.norm(tp1 - tp0, axis=1)
    e12 = np.linalg.norm(tp2 - tp1, axis=1)
    e20 = np.linalg.norm(tp0 - tp2, axis=1)

    # Bimodal edge-length detection (track-line surveys have two distinct spacings)
    all_edges = np.concatenate([e01, e12, e20])
    med       = np.median(all_edges)
    short_    = all_edges[all_edges <= 2 * med]
    long_     = all_edges[all_edges >  2 * med]
    if len(long_) > 0.1 * len(all_edges):
        nom_mm = float(np.median(long_))
        log(f"  Track-line survey detected — using cross-track spacing: {nom_mm:.2f}mm")
        log(f"  (along-track: {np.median(short_):.2f}mm, cross-track: {nom_mm:.2f}mm)")
    else:
        nom_mm = float(med)
    if nom_ft is not None:
        nom_mm = max(nom_mm, nom_ft * scale)
    max_edge = EDGE_FACTOR * nom_mm

    good = s[np.maximum(np.maximum(e01, e12), e20) <= max_edge]
    log(f"  {len(s):,} → {len(good):,} triangles after filter")

    # ── Boundary ──────────────────────────────────────────────────────────────
    log("  Finding boundary ...")
    b_edges   = find_boundary_edges(good)
    b_loops   = order_boundary(b_edges)
    log(f"  Boundary loops: {len(b_loops)}")
    poly_segs = build_poly_segments(b_loops, xn, yn)

    # ── Top surface ───────────────────────────────────────────────────────────
    verts = np.column_stack([xn, yn, zn]).astype(np.float32)
    gp0 = verts[good[:, 0]]
    gp1 = verts[good[:, 1]]
    gp2 = verts[good[:, 2]]
    nz  = np.cross(gp1 - gp0, gp2 - gp0)[:, 2]
    gp1c = np.where(nz[:, None] >= 0, gp1, gp2)
    gp2c = np.where(nz[:, None] >= 0, gp2, gp1)
    top_tris = np.stack([gp0, gp1c, gp2c], axis=1)

    # ── Side walls ────────────────────────────────────────────────────────────
    wall_list = []
    for b_loop in b_loops:
        idx      = np.array(b_loop, dtype=np.int32)
        idx_next = np.roll(idx, -1)
        T0 = verts[idx]
        T1 = verts[idx_next]
        B0 = np.column_stack([xn[idx],      yn[idx],      np.zeros(len(idx))]).astype(np.float32)
        B1 = np.column_stack([xn[idx_next], yn[idx_next], np.zeros(len(idx))]).astype(np.float32)
        wall_list.append(np.vstack([np.stack([T0, T1, B1], axis=1),
                                    np.stack([T0, B1, B0], axis=1)]))
    wall_tris = np.vstack(wall_list) if wall_list else np.zeros((0, 3, 3), dtype=np.float32)

    # ── Bottom face ───────────────────────────────────────────────────────────
    log("Building bottom face ...")
    n_cols = int(np.ceil(xn.max() / _cell)) + 1
    n_rows = int(np.ceil(yn.max() / _cell)) + 1
    col_g, row_g = np.meshgrid(np.arange(n_cols), np.arange(n_rows))
    cx_all = (col_g.ravel() * _cell + _cell / 2).astype(float)
    cy_all = (row_g.ravel() * _cell + _cell / 2).astype(float)
    inside = batch_pip(cx_all, cy_all, poly_segs)
    inside_idx   = np.where(inside)[0]
    interior_arr = np.column_stack([col_g.ravel()[inside_idx],
                                    row_g.ravel()[inside_idx]]).astype(np.int32)

    interior_cols_by_row = defaultdict(list)
    for col, row in interior_arr:
        interior_cols_by_row[int(row)].append(int(col))
    for row in interior_cols_by_row:
        interior_cols_by_row[row].sort()

    text_cells: set = set()
    if label.strip():
        pix, w_cells, h_cells = render_text_bitmap(label, _font, _cell)
        log(f"  Label '{label}' — {w_cells}×{h_cells} cells at {_font}pt")
        text_cells = place_label(pix, w_cells, h_cells, interior_cols_by_row, n_rows)
        placed = sum(1 for c in text_cells
                     if (c[0], c[1]) in {(int(r[0]), int(r[1])) for r in interior_arr})
        log(f"  Placed {placed}/{int(np.count_nonzero(pix))} pixels inside footprint")
    else:
        log("  No label — flat bottom face")

    bottom_tris = build_bottom_face(interior_arr, text_cells, _cell, _recess)
    log(f"  {len(interior_arr):,} bottom cells")

    # ── Write ─────────────────────────────────────────────────────────────────
    log(f"Writing {out_file} ...")
    all_tris = np.vstack([top_tris, wall_tris, bottom_tris]).astype(np.float32)
    write_stl(all_tris, out_file)
    log(f"\n✓  {os.path.basename(out_file)}  ({len(all_tris):,} triangles)")


# ══════════════════════════════════════════════════════════════════════════════
#  GUI
# ══════════════════════════════════════════════════════════════════════════════

def label_from_stem(stem):
    m = re.search(r'(Pier\s*\d+)', stem, re.IGNORECASE)
    return m.group(1) if m else ''


def get_base_dir():
    """Return the directory containing the exe or script, handling PyInstaller."""
    if getattr(sys, 'frozen', False):
        return os.path.dirname(sys.executable)
    return os.path.dirname(os.path.abspath(__file__))


def get_data_dir(subfolder):
    """
    Return the path for XYZ Files / STL Files folders.
    When running as exe: same folder as the exe.
    When running as script: one level up from Scripts (alongside Scripts).
    """
    if getattr(sys, 'frozen', False):
        base = os.path.dirname(sys.executable)
        return os.path.normpath(os.path.join(base, subfolder))
    else:
        base = os.path.dirname(os.path.abspath(__file__))
        return os.path.normpath(os.path.join(base, '..', subfolder))


def launch_gui():
    import tkinter as tk
    from tkinter import ttk, filedialog, scrolledtext

    class App(tk.Tk):
        def __init__(self):
            super().__init__()
            self.title("BathyModel v1.1")
            self.resizable(True, True)
            self.columnconfigure(0, weight=1)
            self._queue       = []
            self._merge_var   = tk.IntVar(value=0)  # 0=batch, 1=merge
            self._last_output = None

            # ── Icon (embedded) ───────────────────────────────────────────────
            try:
                import base64, tempfile
                ico_data  = base64.b64decode(_ICON_B64)
                tmp       = tempfile.NamedTemporaryFile(delete=False, suffix='.ico')
                tmp.write(ico_data)
                tmp.close()
                self.iconbitmap(tmp.name)
            except Exception:
                pass

            pad      = dict(padx=10, pady=4)
            LABEL_W  = 28
            ENTRY_W  = 8
            BROWSE_W = 10
            self.param_vars = {}

            # Consistent 4-column layout across all sections:
            # col 0: label (LABEL_W)  col 1: primary entry (expands)
            # col 2: secondary label  col 3: secondary entry or Browse button (fixed right)
            COL2_W = 22   # second column label width — keeps col 3 aligned

            def cfg(frame):
                frame.columnconfigure(1, weight=1)
                frame.columnconfigure(2, minsize=COL2_W * 7)  # approx px per char
                frame.columnconfigure(3, minsize=BROWSE_W * 8)

            # ── Input ─────────────────────────────────────────────────────────
            fr = ttk.LabelFrame(self, text="Input")
            fr.grid(row=0, column=0, sticky="ew", **pad)
            cfg(fr)

            ttk.Label(fr, text=".XYZ File(s) Location:", width=LABEL_W, anchor="w").grid(
                row=0, column=0, sticky="w", **pad)
            self.xyz_var = tk.StringVar()
            ttk.Entry(fr, textvariable=self.xyz_var).grid(
                row=0, column=1, columnspan=2, sticky="ew", **pad)
            ttk.Button(fr, text="Browse…", width=BROWSE_W,
                       command=self._browse_xyz).grid(row=0, column=3, sticky="e", **pad)

            ttk.Label(fr, text="Nominal Grid Spacing (ft):", width=LABEL_W, anchor="w").grid(
                row=1, column=0, sticky="w", **pad)
            self.spacing_var = tk.StringVar(value="15")
            ttk.Entry(fr, textvariable=self.spacing_var, width=ENTRY_W, justify="center").grid(
                row=1, column=1, sticky="w", **pad)
            ttk.Label(fr, text="(Optimised for 15x15 ft Shoal Bias Files — leave blank to auto-estimate)").grid(
                row=1, column=2, columnspan=2, sticky="w", **pad)

            # Merge / batch toggle — hidden until >1 file selected
            self._merge_frame = ttk.Frame(fr)
            self._merge_frame.grid(row=2, column=0, columnspan=4, sticky="w", **pad)
            self._merge_frame.grid_remove()
            ttk.Label(self._merge_frame, text="Multi-file mode:", width=LABEL_W, anchor="w").grid(
                row=0, column=0, sticky="w")
            ttk.Radiobutton(self._merge_frame, text="Process as separate files",
                            variable=self._merge_var, value=0,
                            command=self._on_mode_change).grid(row=0, column=1, padx=(0, 20))
            ttk.Radiobutton(self._merge_frame, text="Merge into single model",
                            variable=self._merge_var, value=1,
                            command=self._on_mode_change).grid(row=0, column=2)

            # ── Parameters ────────────────────────────────────────────────────
            fr3 = ttk.LabelFrame(self, text="Parameters")
            fr3.grid(row=1, column=0, sticky="ew", **pad)
            cfg(fr3)

            ttk.Label(fr3, text="Z-Height Exaggeration:", width=LABEL_W, anchor="w").grid(
                row=0, column=0, sticky="w", **pad)
            var = tk.StringVar(value="3.0")
            ttk.Entry(fr3, textvariable=var, width=ENTRY_W, justify="center").grid(
                row=0, column=1, sticky="w", **pad)
            self.param_vars["z_exag"] = var

            ttk.Label(fr3, text="Base Padding Thickness (mm):", width=COL2_W, anchor="w").grid(
                row=0, column=2, sticky="w", **pad)
            var = tk.StringVar(value="5.0")
            ttk.Entry(fr3, textvariable=var, width=ENTRY_W, justify="center").grid(
                row=0, column=3, sticky="e", **pad)
            self.param_vars["base_mm"] = var

            # ── Label ─────────────────────────────────────────────────────────
            fr2 = ttk.LabelFrame(self, text="Label")
            fr2.grid(row=2, column=0, sticky="ew", **pad)
            cfg(fr2)

            ttk.Label(fr2, text="Label Text:", width=LABEL_W, anchor="w").grid(
                row=0, column=0, sticky="w", **pad)
            self.label_var = tk.StringVar()
            self.label_entry = ttk.Entry(fr2, textvariable=self.label_var)
            self.label_entry.grid(row=0, column=1, columnspan=3, sticky="ew", **pad)

            ttk.Label(fr2, text="Font Size (pt):", width=LABEL_W, anchor="w").grid(
                row=1, column=0, sticky="w", **pad)
            var = tk.StringVar(value="14")
            ttk.Entry(fr2, textvariable=var, width=ENTRY_W, justify="center").grid(
                row=1, column=1, sticky="w", **pad)
            self.param_vars["font_pt"] = var

            ttk.Label(fr2, text="Text Recess Depth (mm):", width=COL2_W, anchor="w").grid(
                row=1, column=2, sticky="w", **pad)
            var = tk.StringVar(value="1.5")
            ttk.Entry(fr2, textvariable=var, width=ENTRY_W, justify="center").grid(
                row=1, column=3, sticky="e", **pad)
            self.param_vars["recess_mm"] = var

            # ── Output ────────────────────────────────────────────────────────
            fr4 = ttk.LabelFrame(self, text="Output")
            fr4.grid(row=3, column=0, sticky="ew", **pad)
            cfg(fr4)

            ttk.Label(fr4, text=".STL File(s) Output Location:", width=LABEL_W, anchor="w").grid(
                row=0, column=0, sticky="w", **pad)
            self.stl_var = tk.StringVar()
            self.stl_entry = ttk.Entry(fr4, textvariable=self.stl_var)
            self.stl_entry.grid(row=0, column=1, columnspan=2, sticky="ew", **pad)
            self.stl_browse_btn = ttk.Button(fr4, text="Browse…", width=BROWSE_W,
                                             command=self._browse_stl)
            self.stl_browse_btn.grid(row=0, column=3, sticky="e", **pad)

            # ── Run ───────────────────────────────────────────────────────────
            btn_frame = ttk.Frame(self)
            btn_frame.grid(row=4, column=0, pady=8)

            self.run_btn = ttk.Button(btn_frame, text="Generate STL  ▶", command=self._run)
            self.run_btn.grid(row=0, column=0, padx=10)

            self.open_btn = ttk.Button(btn_frame, text="Open Output File",
                                       command=self._open_output, state="disabled")
            self.open_btn.grid(row=0, column=1, padx=10)

            # ── Log ───────────────────────────────────────────────────────────
            fr5 = ttk.LabelFrame(self, text="Log")
            fr5.grid(row=5, column=0, sticky="nsew", **pad)
            fr5.columnconfigure(0, weight=1)
            fr5.rowconfigure(0, weight=1)
            self.rowconfigure(5, weight=1)
            self.log = scrolledtext.ScrolledText(fr5, height=14, state="disabled",
                                                 font=("Courier", 9))
            self.log.grid(row=0, column=0, sticky="nsew", padx=4, pady=4)
            self.minsize(580, 500)

        def _open_output(self):
            path = self._last_output or self.stl_var.get().strip()
            if path and os.path.exists(path):
                os.startfile(path)

        def _browse_xyz(self):
            xyz_dir = get_data_dir('XYZ Files')
            os.makedirs(xyz_dir, exist_ok=True)
            paths = filedialog.askopenfilenames(
                title="Select XYZ survey file(s)",
                initialdir=xyz_dir,
                filetypes=[("Survey files", "*.xyz *.txt *.pdf *.bag"), ("XYZ files", "*.xyz"),
                           ("NOAA TXT files", "*.txt"), ("NOAA GeoImage PDF", "*.pdf"),
                           ("BAG files", "*.bag"), ("All files", "*.*")]
            )
            if not paths:
                return
            self._queue = list(paths)
            stl_dir = get_data_dir('STL Files')
            os.makedirs(stl_dir, exist_ok=True)
            if len(paths) == 1:
                path = paths[0]
                self.xyz_var.set(path)
                stem = os.path.splitext(os.path.basename(path))[0]
                self.stl_var.set(os.path.normpath(os.path.join(stl_dir, stem + '.stl')))
                self.label_var.set(label_from_stem(stem))
                self._merge_var.set(0)
                self._merge_frame.grid_remove()
                self._set_batch_mode(False)
                # Clear spacing for NOAA formats — units are metres after
                # projection so a ft-based default would shred the mesh
                if path.lower().endswith(('.txt', '.pdf', '.bag')):
                    self.spacing_var.set("")
            else:
                stems = [os.path.splitext(os.path.basename(p))[0] for p in paths]
                self.xyz_var.set(f"{len(paths)} files selected")
                self._merge_frame.grid()          # show merge toggle
                self._on_mode_change()            # apply current mode to UI state

        def _set_batch_mode(self, enabled):
            """Lock label/STL fields when running separate-file batch."""
            state = "disabled" if enabled else "normal"
            self.label_entry.configure(state=state)
            self.stl_entry.configure(state=state)
            self.stl_browse_btn.configure(state=state)

        def _on_mode_change(self):
            """Called when merge/batch radio changes or when multi-file is first selected."""
            stl_dir = get_data_dir('STL Files')
            if self._merge_var.get() == 1:
                # Merge mode — single output, editable label + STL path
                stems = [os.path.splitext(os.path.basename(p))[0] for p in self._queue]
                merged_stem = stems[0] + "_merged"
                self.stl_var.set(os.path.normpath(os.path.join(stl_dir, merged_stem + '.stl')))
                self.label_var.set(label_from_stem(stems[0]))
                self._set_batch_mode(False)
            else:
                # Batch mode — separate outputs, lock label/STL fields
                self.stl_var.set(f"Auto → STL Files\\  ({len(self._queue)} files)")
                self.label_var.set("(auto-generated from filename)")
                self._set_batch_mode(True)

        def _browse_stl(self):
            script_dir = get_base_dir()
            stl_dir    = get_data_dir('STL Files')
            os.makedirs(stl_dir, exist_ok=True)
            path = filedialog.asksaveasfilename(
                title="Save STL as",
                initialdir=stl_dir,
                defaultextension=".stl",
                filetypes=[("STL files", "*.stl"), ("All files", "*.*")]
            )
            if path:
                self.stl_var.set(path)

        def _log(self, text):
            self.log.configure(state="normal")
            self.log.insert("end", text + "\n")
            self.log.see("end")
            self.log.configure(state="disabled")

        def _run(self):
            if not self._queue:
                self._log("ERROR: No XYZ file selected.")
                return
            stl_path = self.stl_var.get().strip()
            if len(self._queue) == 1 and not stl_path:
                self._log("ERROR: No output STL path set.")
                return
            if len(self._queue) > 1 and self._merge_var.get() == 1 and not stl_path:
                self._log("ERROR: No output STL path set.")
                return
            self.run_btn.configure(state="disabled")
            self.log.configure(state="normal")
            self.log.delete("1.0", "end")
            self.log.configure(state="disabled")
            spc = self.spacing_var.get().strip()
            params = dict(
                nom_ft    = float(spc) if spc else None,
                base_mm   = float(self.param_vars["base_mm"].get()),
                z_exag    = float(self.param_vars["z_exag"].get()),
                recess_mm = float(self.param_vars["recess_mm"].get()),
                font_pt   = int(self.param_vars["font_pt"].get()),
            )
            if len(self._queue) > 1 and self._merge_var.get() == 1:
                # Merge mode — single run with concatenated point cloud
                threading.Thread(
                    target=self._run_merge_worker,
                    args=(self._queue[:], self.label_var.get().strip(), stl_path, params),
                    daemon=True).start()
            else:
                # Single file or batch mode
                if len(self._queue) == 1:
                    jobs = [(self._queue[0], self.label_var.get().strip(), stl_path)]
                else:
                    stl_dir = get_data_dir('STL Files')
                    jobs = []
                    for path in self._queue:
                        stem = os.path.splitext(os.path.basename(path))[0]
                        jobs.append((path, label_from_stem(stem),
                                     os.path.normpath(os.path.join(stl_dir, stem + '.stl'))))
                threading.Thread(target=self._run_worker,
                                 args=(jobs, params), daemon=True).start()

        def _run_worker(self, jobs, params):
            def gui_log(text):
                self.after(0, self._log, text)

            try:
                for idx, (xyz, lbl, stl) in enumerate(jobs):
                    if len(jobs) > 1:
                        gui_log(f"\n── File {idx+1}/{len(jobs)}: "
                                f"{os.path.basename(xyz)} ──")
                    try:
                        run_pipeline(xyz_file=xyz, label=lbl, out_file=stl,
                                     log=gui_log, **params)
                        self._last_output = stl
                        gui_log(f"✓ {os.path.basename(stl)}")
                    except Exception as e:
                        import traceback
                        gui_log(f"ERROR on {os.path.basename(xyz)}: {e}\n"
                                f"{traceback.format_exc()}")
                total = len(jobs)
                gui_log(f"\n✓ Done — {total} file{'s' if total > 1 else ''} processed.")
                self.after(0, lambda: self.open_btn.configure(state="normal"))
            finally:
                self.after(0, lambda: self.run_btn.configure(state="normal"))
                self.after(0, lambda: self._set_batch_mode(False))

        def _run_merge_worker(self, file_list, label, out_stl, params):
            import traceback
            def gui_log(text):
                self.after(0, self._log, text)

            try:
                # ── Merge log header ──────────────────────────────────────────
                gui_log(f"Merging {len(file_list)} files:")
                arrays_x, arrays_y, arrays_z = [], [], []
                coord_unit = None
                for i, fpath in enumerate(file_list):
                    x, y, z, cu = _load_xyz_array(fpath, log=gui_log)
                    arrays_x.append(x)
                    arrays_y.append(y)
                    arrays_z.append(z)
                    if coord_unit is None:
                        coord_unit = cu
                    gui_log(f"  [{i+1}] {os.path.basename(fpath)}  ({len(x):,} pts)")

                # ── Coordinate sanity check ───────────────────────────────────
                from tkinter import messagebox
                ok, msg = _check_coordinate_compatibility(file_list)
                if ok:
                    # compute overlap for log
                    bounds = []
                    for fpath in file_list:
                        xi, yi, zi, _ = _load_xyz_array(fpath)
                        bounds.append((xi.min(), xi.max(), yi.min(), yi.max()))
                    x_ov = min(b[1] for b in bounds) - max(b[0] for b in bounds)
                    y_ov = min(b[3] for b in bounds) - max(b[2] for b in bounds)
                    x_ov = max(x_ov, 0); y_ov = max(y_ov, 0)
                    gui_log(f"  Coordinate check: compatible "
                            f"(X overlap: {x_ov:.0f} {coord_unit}, "
                            f"Y overlap: {y_ov:.0f} {coord_unit})")
                else:
                    proceed = messagebox.askyesno(
                        "Coordinate Warning",
                        f"Coordinate compatibility check failed:\n{msg}\n\n"
                        "Files may be from different areas or coordinate systems.\n"
                        "Proceed anyway?",
                        icon="warning"
                    )
                    if not proceed:
                        gui_log("Merge cancelled by user.")
                        return

                xs = np.concatenate(arrays_x)
                ys = np.concatenate(arrays_y)
                zs = np.concatenate(arrays_z)
                gui_log(f"  Combined: {len(xs):,} points\n")

                run_pipeline(
                    xyz_file   = f"<{len(file_list)}-file merge>",
                    label      = label,
                    out_file   = out_stl,
                    xyz_points = (xs, ys, zs),
                    coord_unit = coord_unit,
                    log        = gui_log,
                    **params,
                )
                gui_log(f"\n✓ Done — merge complete.")
                self._last_output = out_stl
                self.after(0, lambda: self.open_btn.configure(state="normal"))
            except Exception as e:
                gui_log(f"ERROR: {e}\n{traceback.format_exc()}")
            finally:
                self.after(0, lambda: self.run_btn.configure(state="normal"))

    app = App()
    app.mainloop()


# ══════════════════════════════════════════════════════════════════════════════
#  ENTRY POINT
# ══════════════════════════════════════════════════════════════════════════════

def cli_main():
    xyz_file = sys.argv[1]
    label    = sys.argv[2]
    nom_ft   = float(sys.argv[4]) if len(sys.argv) > 4 else None

    if len(sys.argv) > 3:
        out_file = sys.argv[3]
    else:
        script_dir = get_base_dir()
        stl_dir    = get_data_dir('STL Files')
        os.makedirs(stl_dir, exist_ok=True)
        stem     = os.path.splitext(os.path.basename(xyz_file))[0]
        out_file = os.path.join(stl_dir, stem + '.stl')
        print(f"  Output: {os.path.normpath(out_file)}")

    run_pipeline(xyz_file=xyz_file, label=label, out_file=out_file, nom_ft=nom_ft)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        cli_main()
    else:
        launch_gui()
