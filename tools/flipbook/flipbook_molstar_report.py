#!/usr/bin/env python3
"""Build a Mol* Flipbook report for RMSX PDB slice outputs."""

import argparse
import csv
import html
import json
from pathlib import Path

from flipbook_report_common import (
    build_residue_payload,
    read_rmsx_table,
    read_slices,
    summarize_slices,
)


MOLSTAR_VERSION = "5.4.2"
MOLSTAR_JS_URL = f"https://cdn.jsdelivr.net/npm/molstar@{MOLSTAR_VERSION}/build/viewer/molstar.js"
MOLSTAR_CSS_URL = f"https://cdn.jsdelivr.net/npm/molstar@{MOLSTAR_VERSION}/build/viewer/molstar.css"
MOLSTAR_STATE_TRANSFORMS_URL = f"https://cdn.jsdelivr.net/npm/molstar@{MOLSTAR_VERSION}/lib/mol-plugin-state/transforms.js"
MOLSTAR_STATE_TRANSFORMS_BUNDLE_URL = f"https://esm.sh/molstar@{MOLSTAR_VERSION}/lib/mol-plugin-state/transforms.js?bundle"
MANIFEST_SCHEMA_VERSION = "flipbook-molstar-viewer/v1"
MASK_OPACITY = 0.30
MIN_TILE_SPACING_FACTOR = 0.0
MAX_TILE_SPACING_FACTOR = 2.5
DEFAULT_TILE_SPACING_FACTOR = 1.0

FLIPBOOK_PALETTES = {
    "magma": [
        "#000004", "#120D32", "#331068", "#5A167E", "#7D2482",
        "#A3307E", "#C83E73", "#E95562", "#F97C5D", "#FEA873",
        "#FED395", "#FCFDBF",
    ],
    "inferno": [
        "#000004", "#140B35", "#3A0963", "#60136E", "#85216B",
        "#A92E5E", "#CB4149", "#E65D2F", "#F78311", "#FCAD12",
        "#F5DB4B", "#FCFFA4",
    ],
    "plasma": [
        "#0D0887", "#3E049C", "#6300A7", "#8707A6", "#A62098",
        "#C03A83", "#D5546E", "#E76F5A", "#F58C46", "#FDAD32",
        "#FCD225", "#F0F921",
    ],
    "viridis": [
        "#440154", "#482173", "#433E85", "#38598C", "#2D708E",
        "#25858E", "#1E9B8A", "#2BB07F", "#51C56A", "#85D54A",
        "#C2DF23", "#FDE725",
    ],
    "cividis": [
        "#00204D", "#00306F", "#2A406C", "#48526B", "#5E626E",
        "#727374", "#878479", "#9E9677", "#B6A971", "#D0BE67",
        "#EAD357", "#FFEA46",
    ],
    "rocket": [
        "#03051A", "#221331", "#451C47", "#6A1F56", "#921C5B",
        "#B91657", "#D92847", "#ED513E", "#F47C56", "#F6A47B",
        "#F7C9AA", "#FAEBDD",
    ],
    "mako": [
        "#0B0405", "#231526", "#35264C", "#403A75", "#3D526D",
        "#366DA0", "#3487A6", "#35A1AB", "#43BBAD", "#6CD3AD",
        "#ADE3C0", "#DEF5E5",
    ],
    "turbo": [
        "#30123B", "#4454C4", "#4490FE", "#1FC8DE", "#29EFA2",
        "#7DFF56", "#C1F334", "#F1CA3A", "#FE922A", "#EA4F0D",
        "#BE2102", "#7A0403",
    ],
}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdb-dir", required=True, help="Directory containing slice_*_first_frame.pdb files.")
    parser.add_argument("--rmsx-table", required=True, help="RMSX CSV table.")
    parser.add_argument("--mask-table", required=True, help="Mask metadata CSV table.")
    parser.add_argument("--output", help="Optional standalone HTML report output path for development/debugging.")
    parser.add_argument("--manifest-output", help="Optional native Galaxy visualization manifest JSON output path.")
    parser.add_argument("--palette", default="viridis", choices=sorted(FLIPBOOK_PALETTES), help="Flipbook color palette.")
    parser.add_argument("--title", default="RMSX Flipbook viewer", help="Report title.")
    return parser.parse_args()


def parse_mask_bool(value):
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


def read_mask_summary(path):
    masked = []
    total = 0
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            total += 1
            if not parse_mask_bool(row.get("Masked", "")):
                continue
            residue_id = (row.get("ResidueID") or "").strip()
            if not residue_id:
                continue
            chain_id = (row.get("ChainID") or "").strip()
            key = f"{chain_id}:{residue_id}" if chain_id else residue_id
            masked.append(
                {
                    "id": residue_id,
                    "chain": chain_id,
                    "key": key,
                    "label": f"{residue_id} / chain {chain_id}" if chain_id else residue_id,
                }
            )
    return {
        "maskedResidues": len(masked),
        "totalResidues": total,
        "masked": masked,
        "maskedKeys": [residue["key"] for residue in masked],
    }


def build_flipbook_reference(slices, domain, palette_name, palette_colors):
    interval = 0
    if len(palette_colors) > 1:
        interval = (domain["max"] - domain["min"]) / (len(palette_colors) - 1)
    color_stops = [
        {
            "bfactor": round(domain["min"] + (index * interval), 2),
            "color": color,
        }
        for index, color in enumerate(palette_colors)
    ]
    color_mapping = ":".join(f"{stop['bfactor']},{stop['color']}" for stop in color_stops)
    num_models = len(slices)
    return {
        "viewer": "chimerax",
        "defaultColumns": num_models,
        "minimumSpacingFactor": MIN_TILE_SPACING_FACTOR,
        "maximumSpacingFactor": MAX_TILE_SPACING_FACTOR,
        "defaultSpacingFactor": DEFAULT_TILE_SPACING_FACTOR,
        "tilePaddingFactor": 1.55,
        "palette": palette_name,
        "colorStops": color_stops,
        "commands": [
            "view",
            "define axis",
            "turn x 90",
            "color byattribute bfactor",
            "worm bfactor",
            "lighting soft",
            "graphics silhouettes true",
            "set bgColor white",
            f"color byattribute a:bfactor #1-{num_models} target absc palette {color_mapping}",
            f"tile all columns {num_models} spacingFactor {DEFAULT_TILE_SPACING_FACTOR:g}",
        ],
    }


def build_viewer_payload(title, slices, summaries, domain, mask_summary, residues, palette_name):
    palette_colors = FLIPBOOK_PALETTES[palette_name]
    domain_span = max(0.000001, domain["max"] - domain["min"])
    color_domain_step = round(max(0.1, domain_span / 50), 3)
    pdb_byte_sizes = [len(slice_entry["pdb"].encode("utf-8")) for slice_entry in slices]
    payload = {
        "schemaVersion": MANIFEST_SCHEMA_VERSION,
        "title": title,
        "citation": {
            "text": "RMSX/Flipbook paper, Scientific Reports (2026)",
            "doi": "10.1038/s41598-026-39869-7",
            "url": "https://doi.org/10.1038/s41598-026-39869-7",
        },
        "slices": slices,
        "summaries": summaries,
        "domain": domain,
        "maskSummary": mask_summary,
        "residues": residues,
        "palette": {
            "name": palette_name,
            "colors": palette_colors,
        },
        "availablePalettes": {
            name: FLIPBOOK_PALETTES[name]
            for name in sorted(FLIPBOOK_PALETTES)
        },
        "maskOpacity": MASK_OPACITY,
        "presentation": {
            "defaultLayout": "tiled",
            "availableLayouts": ["tiled", "overlay", "flip"],
            "layoutUrlParam": "layout",
        },
        "playback": {
            "defaultDelayMs": 900,
            "minDelayMs": 150,
            "maxDelayMs": 5000,
            "delayStepMs": 50,
            "delayUrlParam": "delayMs",
        },
        "rotationModel": {
            "mode": "per-slice local coordinate transform",
            "pivot": "geometric center of each full slice before mask splitting",
            "layoutOrder": "rotate around local center, then place that center on a shared Flipbook slot anchor plus tile offset",
            "dragSemantics": "convert Molstar screen-axis drag deltas into coordinate-space rotation matrices, then apply the same delta independently to every visible slice",
            "verticalDragSign": "screen dy is applied directly around the current screen-right axis",
            "defaultRotation": {"x": 90, "y": 0, "z": 0},
            "defaultRotationSource": "ChimeraX Flipbook command: turn x 90",
        },
        "molstarRenderStyle": {
            "preset": "clean-interactive",
            "backgroundColor": "#ffffff",
            "outline": True,
            "ambientOcclusion": False,
            "illumination": False,
            "softLightingUrlParam": "render=soft",
            "outlineUrlParam": "outline=1",
            "outlineDisableUrlParam": "outline=0",
            "multiSample": False,
        },
        "visualMapping": {
            "mode": "normalize RMSX into the PDB B-factor column, then map 0..1 onto explicit Molstar putty worm radii",
            "bfactorDomain": [0, 1],
            "defaultRadiusMin": 0.63,
            "defaultRadiusMax": 3.18,
            "defaultThicknessScale": 1.0,
            "radiusStep": 0.05,
            "defaultColorMin": domain["min"],
            "defaultColorMax": domain["max"],
            "colorDomainStep": color_domain_step,
            "colorTheme": "uncertainty",
            "molstarUncertaintyReversesColorList": True,
            "paletteOrder": "Flipbook low-to-high palette is reversed before passing to Molstar because Molstar uncertainty coloring reverses its color list internally",
            "radiusUrlParams": ["radiusMin", "radiusMax"],
            "thicknessUrlParam": "thickness",
            "colorDomainUrlParams": ["colorMin", "colorMax"],
            "paletteUrlParam": "palette",
        },
        "selectedResidueMarker": {
            "enabledDefault": False,
            "type": "spacefill",
            "color": "#111827",
            "alpha": 0.86,
            "sizeFactor": 0.36,
            "toggleUrlParam": "marker",
            "residueUrlParam": "residue",
        },
        "keyboardShortcuts": {
            "enabled": True,
            "source": "VMD Flipbook hotkey parity subset",
            "rotationStepDegrees": 5,
            "spacingStep": 0.05,
            "thicknessStep": 0.05,
            "bindings": [
                {"keys": ["u", "i"], "action": "rotate all visible slices around display X by +/-5 degrees"},
                {"keys": ["n", "m"], "action": "rotate all visible slices around display Y by +/-5 degrees"},
                {"keys": ["j", "k"], "action": "rotate all visible slices around display Z by +/-5 degrees"},
                {"keys": ["=", "+", "-"], "action": "adjust tiled spacing factor"},
                {"keys": ["[", "]"], "action": "adjust RMSX worm thickness scale"},
                {"keys": [",", "."], "action": "adjust RMSX color/radius domain"},
                {"keys": ["ArrowLeft", "ArrowRight"], "action": "step the active Flipbook slice"},
                {"keys": ["t", "o", "f"], "action": "switch tiled, overlay, and Flip layouts"},
            ],
        },
        "reportPayload": {
            "embeddedPdbCount": len(slices),
            "embeddedPdbBytes": sum(pdb_byte_sizes),
            "largestEmbeddedPdbBytes": max(pdb_byte_sizes) if pdb_byte_sizes else 0,
            "residueCount": len(residues),
            "largeReportWarningBytes": 10_000_000,
            "strategy": "self-contained Galaxy HTML with embedded PDB slices",
        },
        "flipbookReference": build_flipbook_reference(slices, domain, palette_name, palette_colors),
        "molstar": {
            "version": MOLSTAR_VERSION,
            "jsUrl": MOLSTAR_JS_URL,
            "cssUrl": MOLSTAR_CSS_URL,
            "stateTransformsUrl": MOLSTAR_STATE_TRANSFORMS_URL,
            "stateTransformsUrls": [
                MOLSTAR_STATE_TRANSFORMS_BUNDLE_URL,
                MOLSTAR_STATE_TRANSFORMS_URL,
            ],
        },
    }
    payload["reportPayload"]["estimatedJsonBytes"] = len(json.dumps(payload).encode("utf-8"))
    return payload


def html_report(payload):
    palette_colors = payload["palette"]["colors"]
    payload_json = json.dumps(payload).replace("</", "<\\/")
    escaped_title = html.escape(payload["title"])
    palette_gradient = ", ".join(palette_colors)
    template = """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>__TITLE__</title>
  <link rel="stylesheet" type="text/css" href="__MOLSTAR_CSS_URL__">
  <script src="__MOLSTAR_JS_URL__"></script>
  <style>
    :root {
      color-scheme: light;
      --bg: #f7f8fa;
      --panel: #ffffff;
      --scene: #fbfcfd;
      --ink: #1d2630;
      --muted: #5f6b7a;
      --line: #d7dce2;
      --accent: #1f6feb;
      --accent-ink: #ffffff;
      --soft: #eef4ff;
      --warn: #b42318;
      --green: #13795b;
    }

    * {
      box-sizing: border-box;
    }

    body {
      margin: 0;
      min-height: 100vh;
      background: var(--bg);
      color: var(--ink);
      font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      letter-spacing: 0;
    }

    .app {
      display: grid;
      grid-template-rows: auto 1fr;
      min-height: 100vh;
    }

    header {
      display: grid;
      grid-template-columns: minmax(220px, 1fr) auto;
      gap: 16px;
      align-items: center;
      padding: 12px 16px;
      border-bottom: 1px solid var(--line);
      background: var(--panel);
    }

    h1 {
      margin: 0;
      font-size: 16px;
      font-weight: 650;
    }

    .controls,
    .button-row,
    .layout-tabs {
      display: flex;
      flex-wrap: wrap;
      gap: 8px;
      align-items: center;
    }

    .controls {
      justify-content: flex-end;
    }

    button,
    select,
    input[type="range"] {
      height: 32px;
    }

    button,
    select {
      border: 1px solid var(--line);
      border-radius: 6px;
      background: #fff;
      color: var(--ink);
      font-size: 13px;
    }

    button {
      padding: 0 10px;
      cursor: pointer;
    }

    button.primary {
      border-color: var(--accent);
      background: var(--accent);
      color: var(--accent-ink);
    }

    button.active {
      border-color: var(--accent);
      background: var(--soft);
      color: var(--accent);
      font-weight: 650;
    }

    select {
      min-width: 120px;
      padding: 0 8px;
    }

    input[type="number"] {
      width: 76px;
      height: 32px;
      border: 1px solid var(--line);
      border-radius: 6px;
      padding: 0 8px;
      color: var(--ink);
      font-size: 13px;
      font-variant-numeric: tabular-nums;
    }

    input[type="range"] {
      width: 100%;
      accent-color: var(--accent);
    }

    .input-pair {
      display: grid;
      grid-template-columns: minmax(0, 1fr) auto;
      gap: 8px;
      align-items: center;
    }

    main {
      display: grid;
      grid-template-columns: minmax(280px, 340px) minmax(0, 1fr);
      height: calc(100vh - 57px);
      min-height: 0;
      overflow: hidden;
    }

    aside {
      padding: 14px 16px;
      overflow: auto;
      border-right: 1px solid var(--line);
      background: var(--panel);
    }

    .control-block {
      display: grid;
      gap: 9px;
      padding: 0 0 16px;
      margin: 0 0 14px;
      border-bottom: 1px solid var(--line);
    }

    .control-block:last-child {
      margin-bottom: 0;
      border-bottom: 0;
    }

    .section-title {
      margin: 0;
      color: var(--muted);
      font-size: 12px;
      font-weight: 650;
      text-transform: uppercase;
    }

    .field {
      display: grid;
      gap: 5px;
    }

    label {
      color: var(--muted);
      font-size: 13px;
    }

    .metric {
      display: grid;
      grid-template-columns: 1fr auto;
      gap: 10px;
      align-items: baseline;
      padding: 9px 0;
      border-bottom: 1px solid var(--line);
      font-size: 13px;
    }

    .metric span {
      color: var(--muted);
    }

    .metric strong {
      font-variant-numeric: tabular-nums;
    }

    .legend {
      display: grid;
      gap: 7px;
      color: var(--muted);
      font-size: 12px;
    }

    .bar {
      height: 10px;
      border: 1px solid var(--line);
      border-radius: 999px;
      background: linear-gradient(90deg, __PALETTE_GRADIENT__);
    }

    .legend-values,
    .radius-legend {
      display: grid;
      grid-template-columns: repeat(3, minmax(0, 1fr));
      gap: 8px;
      align-items: center;
    }

    .legend-stop,
    .radius-stop {
      display: inline-flex;
      align-items: center;
      min-width: 0;
      gap: 6px;
      font-variant-numeric: tabular-nums;
    }

    .legend-stop:nth-child(2),
    .radius-stop:nth-child(2) {
      justify-content: center;
    }

    .legend-stop:nth-child(3),
    .radius-stop:nth-child(3) {
      justify-content: flex-end;
    }

    .legend-swatch {
      width: 12px;
      height: 12px;
      flex: 0 0 auto;
      border-radius: 999px;
      border: 1px solid rgba(29, 38, 48, 0.25);
      background: var(--muted);
    }

    .radius-dot {
      width: 8px;
      height: 8px;
      flex: 0 0 auto;
      border-radius: 999px;
      border: 1px solid rgba(29, 38, 48, 0.25);
      background: var(--muted);
    }

    .viewer-shell {
      display: grid;
      grid-template-rows: auto minmax(0, 1fr);
      min-width: 0;
      min-height: 0;
      overflow: hidden;
      background: var(--scene);
    }

    .viewer-toolbar {
      display: grid;
      grid-template-columns: minmax(0, 1fr) auto;
      gap: 12px;
      align-items: center;
      padding: 10px 12px;
      border-bottom: 1px solid var(--line);
      background: #fff;
    }

    .slice-chip {
      display: inline-flex;
      align-items: center;
      min-height: 26px;
      padding: 0 9px;
      border: 1px solid var(--line);
      border-radius: 999px;
      background: #fff;
      color: var(--ink);
      font-size: 12px;
      font-weight: 650;
    }

    .scene-labels {
      display: flex;
      flex-wrap: wrap;
      gap: 6px;
      min-width: 0;
    }

    .scene-label {
      display: inline-flex;
      align-items: center;
      height: 26px;
      min-height: 26px;
      gap: 6px;
      padding: 0 9px;
      border: 1px solid var(--line);
      border-radius: 999px;
      background: #fff;
      color: var(--muted);
      font-size: 12px;
      font-weight: 650;
      cursor: pointer;
    }

    .scene-label.active {
      border-color: var(--accent);
      background: var(--soft);
      color: var(--accent);
    }

    .scene-label.hidden {
      opacity: 0.48;
      background: #f3f5f7;
    }

    .swatch {
      width: 9px;
      height: 9px;
      border-radius: 999px;
      border: 1px solid rgba(29, 38, 48, 0.18);
    }

    .color-value {
      display: inline-flex;
      align-items: center;
      gap: 7px;
    }

    .color-swatch {
      width: 14px;
      height: 14px;
      box-shadow: inset 0 0 0 1px rgba(255, 255, 255, 0.55);
    }

    .status {
      justify-self: end;
      max-width: 430px;
      padding: 6px 9px;
      border: 1px solid var(--line);
      border-radius: 6px;
      background: rgba(255, 255, 255, 0.94);
      color: var(--muted);
      font-size: 12px;
      line-height: 1.35;
    }

    .status.ok {
      color: var(--green);
    }

    .status.error {
      color: var(--warn);
      border-color: #f2c7c7;
      background: #fff5f5;
    }

    #molstarViewport {
      position: relative;
      min-width: 0;
      min-height: 0;
    }

    #molstarViewport.local-rotate-enabled {
      cursor: grab;
      touch-action: none;
    }

    #molstarViewport.local-rotate-dragging {
      cursor: grabbing;
    }

    .rotation-buttons {
      grid-template-columns: repeat(3, minmax(0, 1fr));
    }

    @media (max-width: 920px) {
      header {
        grid-template-columns: 1fr;
      }

      .controls {
        justify-content: flex-start;
      }

      main {
        grid-template-columns: 1fr;
        grid-template-rows: auto minmax(560px, 1fr);
        height: auto;
        overflow: visible;
      }

      aside {
        border-right: 0;
        border-bottom: 1px solid var(--line);
      }

      .viewer-toolbar {
        grid-template-columns: 1fr;
      }

      .status {
        justify-self: stretch;
        max-width: none;
      }

      #molstarViewport {
        min-height: 560px;
      }
    }
  </style>
</head>
<body>
  <div class="app">
    <header>
      <h1>__TITLE__</h1>
      <div class="controls">
        <button id="playButton" class="primary" type="button" data-testid="molstar-play">Play</button>
        <div class="layout-tabs" role="group" aria-label="Layout">
          <button id="tiledButton" class="active" type="button" data-testid="molstar-layout-tiled">Tiled</button>
          <button id="overlayButton" type="button" data-testid="molstar-layout-overlay">Overlay</button>
          <button id="flipButton" type="button" data-testid="molstar-layout-flip">Flip</button>
        </div>
        <button id="resetButton" type="button" data-testid="molstar-reset">Reset View</button>
      </div>
    </header>
    <main>
      <aside>
        <div class="control-block">
          <p class="section-title">Flipbook</p>
          <div class="field">
            <label for="sliceSelect">Slice</label>
            <select id="sliceSelect" data-testid="molstar-slice-select"></select>
          </div>
          <div class="field">
            <label for="sliceRange">Frame</label>
            <input id="sliceRange" type="range" min="0" max="0" value="0" step="1" data-testid="molstar-slice-range">
          </div>
          <div class="field">
            <label for="playDelayRange">Delay</label>
            <div class="input-pair">
              <input id="playDelayRange" type="range" min="__MIN_PLAY_DELAY_MS__" max="__MAX_PLAY_DELAY_MS__" value="__DEFAULT_PLAY_DELAY_MS__" step="__PLAY_DELAY_STEP_MS__" data-testid="molstar-play-delay-range">
              <input id="playDelayNumber" type="number" min="__MIN_PLAY_DELAY_MS__" max="__MAX_PLAY_DELAY_MS__" value="__DEFAULT_PLAY_DELAY_MS__" step="__PLAY_DELAY_STEP_MS__" aria-label="Playback delay milliseconds" data-testid="molstar-play-delay-number">
            </div>
          </div>
          <div class="field">
            <label for="thicknessRange">Thickness</label>
            <div class="input-pair">
              <input id="thicknessRange" type="range" min="0.250" max="2.500" value="1.000" step="0.050" data-testid="molstar-thickness-range">
              <input id="thicknessNumber" type="number" min="0.250" max="2.500" value="1.000" step="0.050" aria-label="Thickness value" data-testid="molstar-thickness-number">
            </div>
          </div>
          <div class="field">
            <label for="spacingRange">Spacing</label>
            <div class="input-pair">
              <input id="spacingRange" type="range" min="__MIN_TILE_SPACING_FACTOR__" max="__MAX_TILE_SPACING_FACTOR__" value="__DEFAULT_TILE_SPACING_FACTOR__" step="0.050" data-testid="molstar-spacing-range">
              <input id="spacingNumber" type="number" min="__MIN_TILE_SPACING_FACTOR__" max="__MAX_TILE_SPACING_FACTOR__" value="__DEFAULT_TILE_SPACING_FACTOR__" step="0.050" aria-label="Spacing value" data-testid="molstar-spacing-number">
            </div>
          </div>
          <div class="field">
            <label for="columnsNumber">Columns</label>
            <input id="columnsNumber" type="number" min="1" max="1" value="1" step="1" aria-label="Tile columns" data-testid="molstar-columns-number">
          </div>
          <div class="field">
            <label for="rotationXRange">Rotation X</label>
            <div class="input-pair">
              <input id="rotationXRange" type="range" min="-180" max="180" value="0" step="1" data-testid="molstar-rotation-x-range">
              <input id="rotationXNumber" type="number" min="-180" max="180" value="0" step="1" aria-label="Rotation X degrees" data-testid="molstar-rotation-x-number">
            </div>
          </div>
          <div class="field">
            <label for="rotationYRange">Rotation Y</label>
            <div class="input-pair">
              <input id="rotationYRange" type="range" min="-180" max="180" value="0" step="1" data-testid="molstar-rotation-y-range">
              <input id="rotationYNumber" type="number" min="-180" max="180" value="0" step="1" aria-label="Rotation Y degrees" data-testid="molstar-rotation-y-number">
            </div>
          </div>
          <div class="field">
            <label for="rotationZRange">Rotation Z</label>
            <div class="input-pair">
              <input id="rotationZRange" type="range" min="-180" max="180" value="0" step="1" data-testid="molstar-rotation-z-range">
              <input id="rotationZNumber" type="number" min="-180" max="180" value="0" step="1" aria-label="Rotation Z degrees" data-testid="molstar-rotation-z-number">
            </div>
          </div>
          <div class="button-row rotation-buttons" role="group" aria-label="Rotation steps">
            <button id="rotateXButton" type="button" data-testid="molstar-rotate-x">X +15</button>
            <button id="rotateYButton" type="button" data-testid="molstar-rotate-y">Y +15</button>
            <button id="rotateZButton" type="button" data-testid="molstar-rotate-z">Z +15</button>
          </div>
          <div class="button-row" role="group" aria-label="Rotation mode">
            <button id="localRotateButton" class="active" type="button" data-testid="molstar-local-rotate">Local Drag</button>
            <button id="resetRotationButton" type="button" data-testid="molstar-reset-rotation">Reset Rotation</button>
          </div>
          <div class="button-row" role="group" aria-label="Step controls">
            <button id="previousButton" type="button" data-testid="molstar-previous">Previous</button>
            <button id="nextButton" type="button" data-testid="molstar-next">Next</button>
          </div>
        </div>
        <div class="control-block">
          <p class="section-title">Metrics</p>
          <div class="metric"><span>Current slice</span><strong id="currentSlice">-</strong></div>
          <div class="metric"><span>Mean RMSX</span><strong id="meanRmsx">-</strong></div>
          <div class="metric"><span>Peak RMSX</span><strong id="maxRmsx">-</strong></div>
          <div class="metric"><span>Peak residue</span><strong id="maxResidue">-</strong></div>
          <div class="metric"><span>Residues</span><strong id="residueCount">-</strong></div>
          <div class="metric"><span>Masked</span><strong id="maskedResidues">-</strong></div>
        </div>
        <div class="control-block">
          <p class="section-title">Residue</p>
          <div class="field">
            <label for="residueSelect">Residue</label>
            <select id="residueSelect" data-testid="molstar-residue-select"></select>
          </div>
          <div class="button-row" role="group" aria-label="Residue marker">
            <button id="markerToggleButton" class="active" type="button" aria-pressed="true" data-testid="molstar-residue-marker-toggle">Marker</button>
          </div>
          <div class="metric"><span>Selected RMSX</span><strong id="selectedRmsx">-</strong></div>
          <div class="metric"><span>Visual radius</span><strong id="selectedWeight">-</strong></div>
          <div class="metric"><span>Visual color</span><strong class="color-value"><span id="selectedColorSwatch" class="swatch color-swatch" aria-hidden="true"></span><span id="selectedColor">-</span></strong></div>
          <div class="metric"><span>Mol* style</span><strong id="molstarStyle">-</strong></div>
        </div>
        <div class="control-block">
          <p class="section-title">Scale</p>
          <div class="field">
            <label for="paletteSelect">Palette</label>
            <select id="paletteSelect" data-testid="molstar-palette-select"></select>
          </div>
          <div class="legend" data-testid="molstar-rmsx-legend">
            <div id="legendColorBar" class="bar" aria-hidden="true"></div>
            <div class="legend-values" aria-label="RMSX color domain">
              <span class="legend-stop"><span id="legendLowSwatch" class="legend-swatch" aria-hidden="true"></span><span id="domainMin">-</span></span>
              <span class="legend-stop"><span id="legendMidSwatch" class="legend-swatch" aria-hidden="true"></span><span id="domainMid">-</span></span>
              <span class="legend-stop"><span id="legendHighSwatch" class="legend-swatch" aria-hidden="true"></span><span id="domainMax">-</span></span>
            </div>
            <div class="radius-legend" aria-label="RMSX radius domain" data-testid="molstar-radius-legend">
              <span class="radius-stop"><span id="legendLowRadius" class="radius-dot" aria-hidden="true"></span><span id="legendLowRadiusLabel">-</span></span>
              <span class="radius-stop"><span id="legendMidRadius" class="radius-dot" aria-hidden="true"></span><span id="legendMidRadiusLabel">-</span></span>
              <span class="radius-stop"><span id="legendHighRadius" class="radius-dot" aria-hidden="true"></span><span id="legendHighRadiusLabel">-</span></span>
            </div>
          </div>
          <div class="field">
            <label for="colorMinNumber">Color low</label>
            <input id="colorMinNumber" type="number" value="0.000" step="0.100" aria-label="Color low RMSX" data-testid="molstar-color-min-number">
          </div>
          <div class="field">
            <label for="colorMaxNumber">Color high</label>
            <input id="colorMaxNumber" type="number" value="1.000" step="0.100" aria-label="Color high RMSX" data-testid="molstar-color-max-number">
          </div>
          <div class="field">
            <label for="radiusMinNumber">Radius low</label>
            <input id="radiusMinNumber" type="number" min="0.050" max="5.000" value="0.630" step="0.050" aria-label="Low RMSX worm radius" data-testid="molstar-radius-min-number">
          </div>
          <div class="field">
            <label for="radiusMaxNumber">Radius high</label>
            <input id="radiusMaxNumber" type="number" min="0.100" max="8.000" value="3.180" step="0.050" aria-label="High RMSX worm radius" data-testid="molstar-radius-max-number">
          </div>
          <div class="button-row" role="group" aria-label="Scale controls">
            <button id="resetScaleButton" type="button" data-testid="molstar-reset-scale">Reset Scale</button>
          </div>
        </div>
      </aside>
      <section class="viewer-shell" aria-label="RMSX Flipbook viewer" data-testid="molstar-report">
        <div class="viewer-toolbar">
          <div id="sceneLabels" class="scene-labels" aria-label="Loaded slices"></div>
          <div id="status" class="status">Loading Molstar...</div>
        </div>
        <div id="molstarViewport" data-testid="molstar-viewport"></div>
        <pre id="molstarDiagnostics" data-testid="molstar-diagnostics" hidden>{}</pre>
      </section>
    </main>
  </div>
  <script>
    const REPORT = __PAYLOAD__;
    const SLICE_COLORS = ["#2c7bb6", "#13795b", "#b42318", "#6f42c1", "#d97706", "#0f766e", "#be123c", "#4b5563", "#7c3aed"];
    const VISUAL_MIN = Number(REPORT.visualMapping?.bfactorDomain?.[0] ?? 0);
    const VISUAL_MAX = Number(REPORT.visualMapping?.bfactorDomain?.[1] ?? 1);
    const DEFAULT_RADIUS_MIN = Number(REPORT.visualMapping?.defaultRadiusMin ?? 0.63);
    const DEFAULT_RADIUS_MAX = Number(REPORT.visualMapping?.defaultRadiusMax ?? 3.18);
    const DEFAULT_THICKNESS_SCALE = Number(REPORT.visualMapping?.defaultThicknessScale ?? 1);
    const DEFAULT_COLOR_MIN = Number(REPORT.visualMapping?.defaultColorMin ?? REPORT.domain.min);
    const DEFAULT_COLOR_MAX = Number(REPORT.visualMapping?.defaultColorMax ?? REPORT.domain.max);
    const DEFAULT_MARKER_ENABLED = REPORT.selectedResidueMarker?.enabledDefault !== false;
    const MARKER_COLOR_HEX = String(REPORT.selectedResidueMarker?.color || "#111827").toUpperCase();
    const MASKED_RESIDUE_KEYS = new Set(REPORT.maskSummary.maskedKeys || []);
    const URL_PARAMS = new URLSearchParams(window.location.search);
    const DEFAULT_PALETTE_NAME = defaultPaletteName();
    const KEYBOARD_ROTATION_STEP = Number(REPORT.keyboardShortcuts?.rotationStepDegrees ?? 5);
    const KEYBOARD_SPACING_STEP = Number(REPORT.keyboardShortcuts?.spacingStep ?? 0.05);
    const KEYBOARD_THICKNESS_STEP = Number(REPORT.keyboardShortcuts?.thicknessStep ?? 0.05);
    const KEYBOARD_COLOR_STEP = Number(REPORT.visualMapping?.colorDomainStep ?? 0.5);
    const RADIUS_STEP = Number(REPORT.visualMapping?.radiusStep ?? 0.05);
    const MIN_PLAY_DELAY_MS = Number(REPORT.playback?.minDelayMs ?? 150);
    const MAX_PLAY_DELAY_MS = Number(REPORT.playback?.maxDelayMs ?? 5000);
    const DEFAULT_PLAY_DELAY_MS = Number(REPORT.playback?.defaultDelayMs ?? 900);
    const PLAY_DELAY_STEP_MS = Number(REPORT.playback?.delayStepMs ?? 50);
    const MIN_TILE_SPACING_FACTOR = Number(REPORT.flipbookReference?.minimumSpacingFactor ?? 0.1);
    const MAX_TILE_SPACING_FACTOR = Number(REPORT.flipbookReference?.maximumSpacingFactor ?? 2.5);
    const DEFAULT_TILE_SPACING_FACTOR = Number(REPORT.flipbookReference?.defaultSpacingFactor ?? 1);

    function defaultRotation() {
      const rotation = REPORT.rotationModel?.defaultRotation || {};
      return {
        x: Number.isFinite(Number(rotation.x)) ? Number(rotation.x) : 0,
        y: Number.isFinite(Number(rotation.y)) ? Number(rotation.y) : 0,
        z: Number.isFinite(Number(rotation.z)) ? Number(rotation.z) : 0
      };
    }

    const DEFAULT_ROTATION = defaultRotation();

    function numericParam(name, min, max, fallback) {
      const raw = URL_PARAMS.get(name);
      if (raw === null || raw.trim() === "") {
        return fallback;
      }
      const value = Number(raw);
      if (Number.isFinite(value)) {
        return clamp(value, min, max);
      }
      return fallback;
    }

    function integerParam(name, min, max, fallback) {
      const raw = URL_PARAMS.get(name);
      if (raw === null || raw.trim() === "") {
        return fallback;
      }
      const value = Math.round(Number(raw));
      if (Number.isFinite(value)) {
        return clamp(value, min, max);
      }
      return fallback;
    }

    function booleanParam(name, fallback = false) {
      const raw = URL_PARAMS.get(name);
      if (raw === null || raw.trim() === "") {
        return fallback;
      }
      return ["1", "true", "yes", "on"].includes(raw.toLowerCase());
    }

    const LAYOUTS = new Set(["tiled", "overlay", "flip"]);
    const MANAGED_URL_PARAMS = [
      "layout", "slice", "slices", "sliceA", "sliceB", "sliceC",
      "thickness", "radiusMin", "radiusMax", "colorMin", "colorMax", "palette",
      "spacing", "columns", "rotX", "rotY", "rotZ", "residue", "marker", "delayMs"
    ];

    function compactNumber(value, precision = 3) {
      const numeric = Number(value);
      if (!Number.isFinite(numeric)) {
        return "";
      }
      const rounded = Number(numeric.toFixed(precision));
      return Object.is(rounded, -0) ? "0" : String(rounded);
    }

    function numbersClose(left, right, tolerance = 0.0005) {
      return Math.abs(Number(left) - Number(right)) <= tolerance;
    }

    function defaultLayoutName() {
      return LAYOUTS.has(REPORT.presentation?.defaultLayout) ? REPORT.presentation.defaultLayout : "tiled";
    }

    function availablePalettes() {
      const palettes = REPORT.availablePalettes || {};
      const names = Object.keys(palettes);
      if (names.length) {
        return palettes;
      }
      return { [REPORT.palette.name || "viridis"]: REPORT.palette.colors || [] };
    }

    function paletteNames() {
      return Object.keys(availablePalettes()).sort((left, right) => left.localeCompare(right));
    }

    function defaultPaletteName() {
      const requested = REPORT.palette?.name || "viridis";
      return availablePalettes()[requested] ? requested : paletteNames()[0] || requested;
    }

    function initialPaletteName() {
      const requested = (URL_PARAMS.get("palette") || "").toLowerCase();
      return availablePalettes()[requested] ? requested : DEFAULT_PALETTE_NAME;
    }

    function paletteLabel(name) {
      return String(name || "")
        .replace(/[-_]+/g, " ")
        .replace(/\\b\\w/g, (letter) => letter.toUpperCase());
    }

    function currentPaletteColors() {
      const palettes = availablePalettes();
      const colors = palettes[stageState?.paletteName] || palettes[DEFAULT_PALETTE_NAME] || REPORT.palette.colors || [];
      return colors.map((hex) => String(hex).toUpperCase());
    }

    function currentMolstarUncertaintyColors() {
      return [...currentPaletteColors()].reverse().map((hex) => Number.parseInt(hex.slice(1), 16));
    }

    function currentPaletteGradient() {
      return currentPaletteColors().join(", ");
    }

    function defaultTileColumns() {
      return clamp(
        Math.round(Number(REPORT.flipbookReference?.defaultColumns ?? REPORT.slices.length)),
        1,
        Math.max(1, REPORT.slices.length)
      );
    }

    function initialLayout() {
      const requested = (URL_PARAMS.get("layout") || "").toLowerCase();
      const fallback = defaultLayoutName();
      return LAYOUTS.has(requested) ? requested : fallback;
    }

    function allSliceIndexes() {
      return REPORT.slices.map((_, index) => index);
    }

    function sliceIndexesFromListParam(name) {
      const raw = URL_PARAMS.get(name);
      if (!raw) {
        return [];
      }
      const indexes = [];
      const seen = new Set();
      raw.split(/[,\\s]+/).forEach((entry) => {
        const requested = Number(entry);
        if (!Number.isInteger(requested) || requested < 1 || requested > REPORT.slices.length) {
          return;
        }
        const index = requested - 1;
        if (!seen.has(index)) {
          seen.add(index);
          indexes.push(index);
        }
      });
      return indexes;
    }

    function initialVisibleSliceIndexes() {
      const listed = sliceIndexesFromListParam("slices");
      if (listed.length) {
        return listed;
      }
      const explicit = ["sliceA", "sliceB", "sliceC"]
        .map(sliceIndexFromParam)
        .filter((index) => index !== null);
      if (explicit.length) {
        return [...new Set(explicit)];
      }
      return allSliceIndexes();
    }

    function expectedManagedUrlParams() {
      const params = {};
      const allVisible = stageState.layout === "flip" || visibleSliceIndexList().length === REPORT.slices.length;
      params.layout = stageState.layout === defaultLayoutName() ? null : stageState.layout;
      params.slice = stageState.currentIndex === 0 ? null : String(stageState.currentIndex + 1);
      params.slices = allVisible ? null : visibleSliceIndexList().map((index) => index + 1).join(",");
      params.sliceA = null;
      params.sliceB = null;
      params.sliceC = null;
      params.thickness = numbersClose(stageState.thicknessScale, DEFAULT_THICKNESS_SCALE)
        ? null
        : compactNumber(stageState.thicknessScale);
      params.radiusMin = numbersClose(stageState.radiusMin, DEFAULT_RADIUS_MIN)
        ? null
        : compactNumber(stageState.radiusMin);
      params.radiusMax = numbersClose(stageState.radiusMax, DEFAULT_RADIUS_MAX)
        ? null
        : compactNumber(stageState.radiusMax);
      params.colorMin = numbersClose(colorDomainMin(), DEFAULT_COLOR_MIN)
        ? null
        : compactNumber(colorDomainMin());
      params.colorMax = numbersClose(colorDomainMax(), DEFAULT_COLOR_MAX)
        ? null
        : compactNumber(colorDomainMax());
      params.palette = stageState.paletteName === DEFAULT_PALETTE_NAME ? null : stageState.paletteName;
      params.spacing = numbersClose(stageState.spacingFactor, DEFAULT_TILE_SPACING_FACTOR)
        ? null
        : compactNumber(stageState.spacingFactor);
      params.columns = stageState.tileColumns === defaultTileColumns() ? null : String(stageState.tileColumns);
      params.rotX = numbersClose(stageState.rotation.x, DEFAULT_ROTATION.x) ? null : compactNumber(stageState.rotation.x, 0);
      params.rotY = numbersClose(stageState.rotation.y, DEFAULT_ROTATION.y) ? null : compactNumber(stageState.rotation.y, 0);
      params.rotZ = numbersClose(stageState.rotation.z, DEFAULT_ROTATION.z) ? null : compactNumber(stageState.rotation.z, 0);
      params.residue = selectedResidue()?.key === defaultResidueKey() ? null : selectedResidue()?.key;
      params.marker = stageState.selectionMarkerEnabled === DEFAULT_MARKER_ENABLED
        ? null
        : stageState.selectionMarkerEnabled ? "1" : "0";
      params.delayMs = stageState.playDelayMs === DEFAULT_PLAY_DELAY_MS ? null : String(stageState.playDelayMs);
      return params;
    }

    function syncUrlState() {
      if (!window.history?.replaceState) {
        return;
      }
      const url = new URL(window.location.href);
      const expected = expectedManagedUrlParams();
      Object.entries(expected).forEach(([name, value]) => {
        if (value === null || value === "") {
          url.searchParams.delete(name);
        } else {
          url.searchParams.set(name, value);
        }
      });
      const current = `${window.location.pathname}${window.location.search}${window.location.hash}`;
      const next = `${url.pathname}${url.search}${url.hash}`;
      if (next !== current) {
        window.history.replaceState({ flipbookMolstarState: true }, "", next);
      }
    }

    function urlStateSummary() {
      const params = new URLSearchParams(window.location.search);
      const expected = expectedManagedUrlParams();
      const current = {};
      const matches = {};
      MANAGED_URL_PARAMS.forEach((name) => {
        current[name] = params.get(name);
        matches[name] = expected[name] === null ? !params.has(name) : params.get(name) === expected[name];
      });
      return {
        href: window.location.href,
        managedParams: current,
        expectedManagedParams: expected,
        synced: Object.values(matches).every(Boolean),
        matches
      };
    }

    const stageState = {
      viewer: null,
      currentIndex: 0,
      timer: null,
      reloadTimer: null,
      interactiveFrame: null,
      pendingInteractiveAutoView: false,
      loadVersion: 0,
      layout: initialLayout(),
      autoplay: booleanParam("autoplay") || booleanParam("play"),
      thicknessScale: numericParam("thickness", 0.25, 2.5, DEFAULT_THICKNESS_SCALE),
      radiusMin: numericParam("radiusMin", 0.05, 5, DEFAULT_RADIUS_MIN),
      radiusMax: numericParam("radiusMax", 0.1, 8, DEFAULT_RADIUS_MAX),
      colorDomainMin: numericParam("colorMin", REPORT.domain.min, REPORT.domain.max, DEFAULT_COLOR_MIN),
      colorDomainMax: numericParam("colorMax", REPORT.domain.min, REPORT.domain.max, DEFAULT_COLOR_MAX),
      paletteName: initialPaletteName(),
      spacingFactor: numericParam("spacing", MIN_TILE_SPACING_FACTOR, MAX_TILE_SPACING_FACTOR, DEFAULT_TILE_SPACING_FACTOR),
      tileColumns: integerParam("columns", 1, Math.max(1, REPORT.slices.length), REPORT.flipbookReference?.defaultColumns ?? REPORT.slices.length),
      playDelayMs: integerParam("delayMs", MIN_PLAY_DELAY_MS, MAX_PLAY_DELAY_MS, DEFAULT_PLAY_DELAY_MS),
      rotation: {
        x: numericParam("rotX", -180, 180, DEFAULT_ROTATION.x),
        y: numericParam("rotY", -180, 180, DEFAULT_ROTATION.y),
        z: numericParam("rotZ", -180, 180, DEFAULT_ROTATION.z)
      },
      rotationMatrix: null,
      localRotateMode: true,
      selectionMarkerEnabled: booleanParam("marker", DEFAULT_MARKER_ENABLED),
      rotationSensitivity: numericParam("rotateSensitivity", 0.1, 3, 1.05),
      rotationDrag: null,
      lastRotationGesture: null,
      keyboardShortcutCount: 0,
      lastKeyboardAction: null,
      lastInteractiveReloadAt: 0,
      interactiveFlushCount: 0,
      lastTransformFlush: null,
      lastVisibilityUpdate: null,
      visibleSliceIndexes: new Set(initialVisibleSliceIndexes()),
      preferStructureTransforms: URL_PARAMS.get("geometry") === "state",
      preferRepresentationTransforms: URL_PARAMS.get("geometry") !== "pdb",
      forcePdbCoordinateTransforms: false,
      stateTransformsModuleStatus: "not-requested",
      stateTransformsModuleUrl: "",
      geometryMode: "pdb-coordinate-rewrite",
      capabilities: {},
      records: [],
      representationMode: "-"
    };

    const playButton = document.getElementById("playButton");
    const resetButton = document.getElementById("resetButton");
    const tiledButton = document.getElementById("tiledButton");
    const overlayButton = document.getElementById("overlayButton");
    const flipButton = document.getElementById("flipButton");
    const previousButton = document.getElementById("previousButton");
    const nextButton = document.getElementById("nextButton");
    const sliceSelect = document.getElementById("sliceSelect");
    const sliceRange = document.getElementById("sliceRange");
    const playDelayRange = document.getElementById("playDelayRange");
    const playDelayNumber = document.getElementById("playDelayNumber");
    const thicknessRange = document.getElementById("thicknessRange");
    const thicknessNumber = document.getElementById("thicknessNumber");
    const spacingRange = document.getElementById("spacingRange");
    const spacingNumber = document.getElementById("spacingNumber");
    const columnsNumber = document.getElementById("columnsNumber");
    const rotationXRange = document.getElementById("rotationXRange");
    const rotationXNumber = document.getElementById("rotationXNumber");
    const rotationYRange = document.getElementById("rotationYRange");
    const rotationYNumber = document.getElementById("rotationYNumber");
    const rotationZRange = document.getElementById("rotationZRange");
    const rotationZNumber = document.getElementById("rotationZNumber");
    const rotateXButton = document.getElementById("rotateXButton");
    const rotateYButton = document.getElementById("rotateYButton");
    const rotateZButton = document.getElementById("rotateZButton");
    const localRotateButton = document.getElementById("localRotateButton");
    const resetRotationButton = document.getElementById("resetRotationButton");
    const molstarViewport = document.getElementById("molstarViewport");
    const residueSelect = document.getElementById("residueSelect");
    const markerToggleButton = document.getElementById("markerToggleButton");
    const paletteSelect = document.getElementById("paletteSelect");
    const colorMinNumber = document.getElementById("colorMinNumber");
    const colorMaxNumber = document.getElementById("colorMaxNumber");
    const radiusMinNumber = document.getElementById("radiusMinNumber");
    const radiusMaxNumber = document.getElementById("radiusMaxNumber");
    const resetScaleButton = document.getElementById("resetScaleButton");
    const statusElement = document.getElementById("status");
    const sceneLabels = document.getElementById("sceneLabels");

    function formatNumber(value) {
      if (typeof value !== "number" || !Number.isFinite(value)) {
        return "-";
      }
      return value.toFixed(3);
    }

    function clamp(value, min, max) {
      return Math.max(min, Math.min(max, value));
    }

    function wrapAngle(value) {
      if (!Number.isFinite(Number(value))) {
        return 0;
      }
      let angle = Number(value);
      while (angle > 180) {
        angle -= 360;
      }
      while (angle < -180) {
        angle += 360;
      }
      return angle;
    }

    function updateRotationControls() {
      const fields = [
        ["x", rotationXRange, rotationXNumber],
        ["y", rotationYRange, rotationYNumber],
        ["z", rotationZRange, rotationZNumber]
      ];
      fields.forEach(([axis, range, number]) => {
        const value = stageState.rotation[axis].toFixed(0);
        range.value = value;
        number.value = value;
      });
      localRotateButton.classList.toggle("active", stageState.localRotateMode);
      molstarViewport.classList.toggle("local-rotate-enabled", stageState.localRotateMode);
    }

    function updateResidueMarkerControls() {
      markerToggleButton.classList.toggle("active", stageState.selectionMarkerEnabled);
      markerToggleButton.setAttribute("aria-pressed", stageState.selectionMarkerEnabled ? "true" : "false");
    }

    function updatePlaybackControls() {
      playButton.textContent = stageState.timer ? "Pause" : "Play";
      playDelayRange.value = String(stageState.playDelayMs);
      playDelayNumber.value = String(stageState.playDelayMs);
    }

    function setStatus(text, isError = false) {
      statusElement.textContent = text;
      statusElement.classList.toggle("error", isError);
      statusElement.classList.toggle("ok", !isError);
      updateDiagnostics();
    }

    function sliceIndexFromParam(name) {
      const params = new URLSearchParams(window.location.search);
      const requested = Number(params.get(name));
      if (Number.isInteger(requested) && requested >= 1 && requested <= REPORT.slices.length) {
        return requested - 1;
      }
      return null;
    }

    function initialSliceIndex() {
      return sliceIndexFromParam("slice") ?? sliceIndexFromParam("sliceA") ?? 0;
    }

    function defaultResidueKey() {
      return REPORT.residues[0]?.key || "";
    }

    function initialResidueKey() {
      const requested = URL_PARAMS.get("residue");
      if (!requested) {
        return defaultResidueKey();
      }
      const match = REPORT.residues.find((residue) => {
        return residue.key === requested || residue.id === requested || residue.label === requested;
      });
      return match?.key || defaultResidueKey();
    }

    function selectedResidue() {
      return REPORT.residues.find((residue) => residue.key === residueSelect.value) || REPORT.residues[0];
    }

    function visibleSliceIndexList() {
      if (stageState.layout === "flip") {
        return [stageState.currentIndex];
      }
      return allSliceIndexes().filter((index) => stageState.visibleSliceIndexes.has(index));
    }

    function hiddenSliceIndexList() {
      return allSliceIndexes().filter((index) => !visibleSliceIndexList().includes(index));
    }

    function isSliceVisible(index) {
      return stageState.layout === "flip"
        ? index === stageState.currentIndex
        : stageState.visibleSliceIndexes.has(index);
    }

    function firstVisibleSliceIndex() {
      return visibleSliceIndexList()[0] ?? 0;
    }

    function colorDomainMin() {
      return Math.min(stageState.colorDomainMin, stageState.colorDomainMax - 0.000001);
    }

    function colorDomainMax() {
      return Math.max(stageState.colorDomainMax, stageState.colorDomainMin + 0.000001);
    }

    function normalizedRmsx(value) {
      if (typeof value !== "number" || !Number.isFinite(value)) {
        return 0;
      }
      const low = colorDomainMin();
      const high = colorDomainMax();
      const span = Math.max(0.000001, high - low);
      return clamp((value - low) / span, 0, 1);
    }

    function visualWeightForResidue(residue, slice) {
      const value = residue?.values?.[slice.rmsxColumn];
      if (typeof value !== "number" || !Number.isFinite(value)) {
        return VISUAL_MIN;
      }
      return VISUAL_MIN + normalizedRmsx(value) * (VISUAL_MAX - VISUAL_MIN);
    }

    function wormRadiusMin() {
      const low = Math.min(stageState.radiusMin, stageState.radiusMax - 0.01);
      return Math.max(0.01, low * stageState.thicknessScale);
    }

    function wormRadiusMax() {
      const high = Math.max(stageState.radiusMax, stageState.radiusMin + 0.01);
      return Math.max(wormRadiusMin() + 0.01, high * stageState.thicknessScale);
    }

    function wormRadiusSpan() {
      return Math.max(0.01, wormRadiusMax() - wormRadiusMin());
    }

    function visualRadiusForResidue(residue, slice) {
      return wormRadiusMin() + normalizedRmsx(residue?.values?.[slice.rmsxColumn]) * wormRadiusSpan();
    }

    function hexToRgb(hex) {
      const normalized = String(hex || "").replace("#", "");
      if (!/^[0-9a-fA-F]{6}$/.test(normalized)) {
        return { r: 0, g: 0, b: 0 };
      }
      return {
        r: Number.parseInt(normalized.slice(0, 2), 16),
        g: Number.parseInt(normalized.slice(2, 4), 16),
        b: Number.parseInt(normalized.slice(4, 6), 16)
      };
    }

    function rgbToHex(rgb) {
      return `#${[rgb.r, rgb.g, rgb.b].map((value) => {
        return clamp(Math.round(value), 0, 255).toString(16).padStart(2, "0").toUpperCase();
      }).join("")}`;
    }

    function interpolateHexColor(leftHex, rightHex, fraction) {
      const left = hexToRgb(leftHex);
      const right = hexToRgb(rightHex);
      const t = clamp(fraction, 0, 1);
      return rgbToHex({
        r: left.r + ((right.r - left.r) * t),
        g: left.g + ((right.g - left.g) * t),
        b: left.b + ((right.b - left.b) * t)
      });
    }

    function expectedColorForNormalizedRmsx(normalized) {
      const colors = currentPaletteColors();
      if (!colors.length) {
        return "#000000";
      }
      if (colors.length === 1) {
        return colors[0];
      }
      const scaled = clamp(normalized, 0, 1) * (colors.length - 1);
      const lower = Math.floor(scaled);
      const upper = Math.ceil(scaled);
      if (lower === upper) {
        return colors[lower];
      }
      return interpolateHexColor(colors[lower], colors[upper], scaled - lower);
    }

    function expectedColorForResidue(residue, slice) {
      return expectedColorForNormalizedRmsx(normalizedRmsx(residue?.values?.[slice.rmsxColumn]));
    }

    function hasMaskedResidues() {
      return MASKED_RESIDUE_KEYS.size > 0;
    }

    function isMaskedResidueKey(chainId, residueId) {
      const chainedKey = chainId ? `${chainId}:${residueId}` : residueId;
      return MASKED_RESIDUE_KEYS.has(chainedKey) || MASKED_RESIDUE_KEYS.has(residueId);
    }

    function atomResidueKeyFromPaddedLine(padded) {
      const residueId = padded.slice(22, 26).trim();
      const chainId = padded.slice(21, 22).trim();
      return {
        residueId,
        chainId,
        key: chainId ? `${chainId}:${residueId}` : residueId
      };
    }

    function atomLineMatchesResidue(padded, residue, options = {}) {
      if (!residue) {
        return false;
      }
      const atomResidue = atomResidueKeyFromPaddedLine(padded);
      if (atomResidue.key === residue.key
        || (residue.chain && atomResidue.chainId === residue.chain && atomResidue.residueId === residue.id)) {
        return true;
      }
      if (options.strictChain && residue.chain) {
        return false;
      }
      return atomResidue.residueId === residue.id
        || atomResidue.residueId === residue.key;
    }

    function uncertaintyColorParams() {
      return {
        domain: [VISUAL_MIN, VISUAL_MAX],
        list: {
          kind: "interpolate",
          colors: currentMolstarUncertaintyColors()
        }
      };
    }

    function mappingLegendStops() {
      const min = colorDomainMin();
      const max = colorDomainMax();
      const mid = min + ((max - min) / 2);
      return [
        { key: "Low", rmsx: min, normalized: 0, radius: wormRadiusMin() },
        { key: "Mid", rmsx: mid, normalized: 0.5, radius: wormRadiusMin() + (wormRadiusSpan() / 2) },
        { key: "High", rmsx: max, normalized: 1, radius: wormRadiusMax() }
      ].map((stop) => ({
        ...stop,
        color: expectedColorForNormalizedRmsx(stop.normalized)
      }));
    }

    function radiusDotSize(radius) {
      return clamp(radius * 5, 7, 24);
    }

    function updateMappingLegend() {
      const colorBar = document.getElementById("legendColorBar");
      if (colorBar) {
        colorBar.style.background = `linear-gradient(90deg, ${currentPaletteGradient()})`;
      }
      mappingLegendStops().forEach((stop) => {
        const swatch = document.getElementById(`legend${stop.key}Swatch`);
        const value = document.getElementById(stop.key === "Low" ? "domainMin" : stop.key === "Mid" ? "domainMid" : "domainMax");
        const dot = document.getElementById(`legend${stop.key}Radius`);
        const label = document.getElementById(`legend${stop.key}RadiusLabel`);
        if (swatch) {
          swatch.style.background = stop.color;
        }
        if (value) {
          value.textContent = formatNumber(stop.rmsx);
        }
        if (dot) {
          const size = radiusDotSize(stop.radius);
          dot.style.width = `${size.toFixed(1)}px`;
          dot.style.height = `${size.toFixed(1)}px`;
          dot.style.background = stop.color;
        }
        if (label) {
          label.textContent = stop.radius.toFixed(2);
        }
      });
    }

    function legendSummary() {
      return {
        activeRmsxColorDomain: {
          min: Number(colorDomainMin().toFixed(4)),
          max: Number(colorDomainMax().toFixed(4))
        },
        radiusRange: {
          min: Number(wormRadiusMin().toFixed(4)),
          max: Number(wormRadiusMax().toFixed(4))
        },
        stops: mappingLegendStops().map((stop) => ({
          key: stop.key.toLowerCase(),
          rmsx: Number(stop.rmsx.toFixed(4)),
          normalized: stop.normalized,
          radius: Number(stop.radius.toFixed(4)),
          color: stop.color
        })),
        elements: {
          legend: Boolean(document.querySelector('[data-testid="molstar-rmsx-legend"]')),
          colorBar: Boolean(document.getElementById("legendColorBar")),
          radiusLegend: Boolean(document.querySelector('[data-testid="molstar-radius-legend"]')),
          lowValue: document.getElementById("domainMin")?.textContent || "",
          midValue: document.getElementById("domainMid")?.textContent || "",
          highValue: document.getElementById("domainMax")?.textContent || "",
          lowRadius: document.getElementById("legendLowRadiusLabel")?.textContent || "",
          midRadius: document.getElementById("legendMidRadiusLabel")?.textContent || "",
          highRadius: document.getElementById("legendHighRadiusLabel")?.textContent || ""
        }
      };
    }

    function updateMetrics() {
      const slice = REPORT.slices[stageState.currentIndex];
      const summary = REPORT.summaries[slice.rmsxColumn] || {};
      const residue = selectedResidue();
      const rmsx = residue?.values?.[slice.rmsxColumn];
      const radius = visualRadiusForResidue(residue, slice);
      const color = expectedColorForResidue(residue, slice);
      document.getElementById("currentSlice").textContent = slice.label;
      document.getElementById("meanRmsx").textContent = formatNumber(summary.mean);
      document.getElementById("maxRmsx").textContent = formatNumber(summary.max);
      document.getElementById("maxResidue").textContent = summary.maxResidue || "-";
      document.getElementById("residueCount").textContent = summary.residueCount || "-";
      document.getElementById("selectedRmsx").textContent = formatNumber(rmsx);
      document.getElementById("selectedWeight").textContent = Number.isFinite(radius) ? radius.toFixed(2) : "-";
      document.getElementById("selectedColor").textContent = color;
      document.getElementById("selectedColorSwatch").style.background = color;
      document.getElementById("molstarStyle").textContent = stageState.representationMode;
      updateMappingLegend();
      updateResidueMarkerControls();
      updatePlaybackControls();
      sliceSelect.value = String(stageState.currentIndex);
      sliceRange.value = String(stageState.currentIndex);
      spacingRange.value = stageState.spacingFactor.toFixed(3);
      spacingNumber.value = stageState.spacingFactor.toFixed(3);
      columnsNumber.value = String(stageState.tileColumns);
      colorMinNumber.value = colorDomainMin().toFixed(3);
      colorMaxNumber.value = colorDomainMax().toFixed(3);
      paletteSelect.value = stageState.paletteName;
      radiusMinNumber.value = stageState.radiusMin.toFixed(3);
      radiusMaxNumber.value = stageState.radiusMax.toFixed(3);
      updateRotationControls();
      updateSceneLabels();
    }

    function populateControls() {
      REPORT.slices.forEach((slice, index) => {
        const option = document.createElement("option");
        option.value = String(index);
        option.textContent = slice.label;
        sliceSelect.appendChild(option);
      });
      REPORT.residues.forEach((residue) => {
        const option = document.createElement("option");
        option.value = residue.key;
        option.textContent = `Residue ${residue.label}`;
        residueSelect.appendChild(option);
      });
      paletteNames().forEach((name) => {
        const option = document.createElement("option");
        option.value = name;
        option.textContent = paletteLabel(name);
        paletteSelect.appendChild(option);
      });
      residueSelect.value = initialResidueKey();
      paletteSelect.value = stageState.paletteName;
      sliceRange.max = String(Math.max(0, REPORT.slices.length - 1));
      updateMappingLegend();
      document.getElementById("maskedResidues").textContent = `${REPORT.maskSummary.maskedResidues} / ${REPORT.maskSummary.totalResidues}`;
      playDelayRange.min = String(MIN_PLAY_DELAY_MS);
      playDelayRange.max = String(MAX_PLAY_DELAY_MS);
      playDelayRange.step = String(PLAY_DELAY_STEP_MS);
      playDelayNumber.min = String(MIN_PLAY_DELAY_MS);
      playDelayNumber.max = String(MAX_PLAY_DELAY_MS);
      playDelayNumber.step = String(PLAY_DELAY_STEP_MS);
      updatePlaybackControls();
      thicknessRange.value = stageState.thicknessScale.toFixed(3);
      thicknessNumber.value = stageState.thicknessScale.toFixed(3);
      spacingRange.value = stageState.spacingFactor.toFixed(3);
      spacingNumber.value = stageState.spacingFactor.toFixed(3);
      colorMinNumber.min = String(REPORT.domain.min);
      colorMinNumber.max = String(REPORT.domain.max);
      colorMinNumber.step = String(KEYBOARD_COLOR_STEP);
      colorMinNumber.value = colorDomainMin().toFixed(3);
      colorMaxNumber.min = String(REPORT.domain.min);
      colorMaxNumber.max = String(REPORT.domain.max);
      colorMaxNumber.step = String(KEYBOARD_COLOR_STEP);
      colorMaxNumber.value = colorDomainMax().toFixed(3);
      radiusMinNumber.step = String(RADIUS_STEP);
      radiusMinNumber.value = stageState.radiusMin.toFixed(3);
      radiusMaxNumber.step = String(RADIUS_STEP);
      radiusMaxNumber.value = stageState.radiusMax.toFixed(3);
      columnsNumber.max = String(Math.max(1, REPORT.slices.length));
      columnsNumber.value = String(stageState.tileColumns);
      updateRotationControls();
      updateResidueMarkerControls();
    }

    function disposeViewer() {
      if (stageState.viewer?.dispose) {
        stageState.viewer.dispose();
      } else if (stageState.viewer?.plugin?.dispose) {
        stageState.viewer.plugin.dispose();
      }
      stageState.viewer = null;
      stageState.records = [];
      document.getElementById("molstarViewport").replaceChildren();
    }

    function hexColorToMolstarNumber(hex, fallback) {
      if (typeof hex !== "string") {
        return fallback;
      }
      const normalized = hex.trim().replace(/^#/, "");
      if (!/^[0-9a-fA-F]{6}$/.test(normalized)) {
        return fallback;
      }
      return Number.parseInt(normalized, 16);
    }

    function activeRenderPreset() {
      const style = REPORT.molstarRenderStyle || {};
      return (URL_PARAMS.get("render") || style.preset || "clean-interactive").toLowerCase();
    }

    function flipbookCanvasProps(options = {}) {
      const style = REPORT.molstarRenderStyle || {};
      const renderPreset = activeRenderPreset();
      const softRender = ["soft", "studio", "cinematic"].includes(renderPreset);
      const backgroundColor = hexColorToMolstarNumber(style.backgroundColor, 0xffffff);
      const outlineParam = (URL_PARAMS.get("outline") || "").toLowerCase();
      const outlineRequested = outlineParam === "1" || outlineParam === "true" || outlineParam === "on";
      const outlineDisabled = outlineParam === "0" || outlineParam === "false" || outlineParam === "off";
      const outlineEnabled = !outlineDisabled && (outlineRequested || style.outline === true || ["outlined", "soft", "studio", "cinematic"].includes(renderPreset));
      const occlusionEnabled = !options.interactive && (softRender ? style.ambientOcclusion !== "never" : style.ambientOcclusion === true);
      const illuminationEnabled = softRender ? style.illumination !== "never" : style.illumination === true;
      return {
        transparentBackground: false,
        dpoitIterations: 1,
        userInteractionReleaseMs: 0,
        multiSample: {
          mode: "off",
          sampleLevel: 0,
          reduceFlicker: false,
          reuseOcclusion: false
        },
        cameraFog: { name: "off", params: {} },
        renderer: {
          backgroundColor,
          ambientIntensity: 0.78
        },
        postprocessing: {
          enabled: outlineEnabled || occlusionEnabled,
          outline: outlineEnabled
            ? { name: "on", params: { scale: 0.55, threshold: 0.22, color: 0x1f2937, includeTransparent: true } }
            : { name: "off", params: {} },
          occlusion: occlusionEnabled
            ? {
                name: "on",
                params: {
                  samples: 8,
                  multiScale: { name: "off", params: {} },
                  radius: 3.2,
                  bias: 0.85,
                  blurKernelSize: 11,
                  blurDepthBias: 0.5,
                  resolutionScale: 0.5,
                  color: 0x000000,
                  transparentThreshold: 0.4
                }
              }
            : { name: "off", params: {} },
          antialiasing: {
            name: "smaa",
            params: {}
          },
          shadow: { name: "off", params: {} },
          dof: { name: "off", params: {} },
          sharpening: { name: "off", params: {} },
          bloom: { name: "off", params: {} },
          background: { variant: { name: "off", params: {} } }
        },
        marking: {
          enabled: false,
          highlightEdgeColor: 0x000000,
          selectEdgeColor: 0x000000,
          ghostEdgeStrength: 0,
          innerEdgeFactor: 1
        },
        illumination: {
          enabled: !options.interactive && illuminationEnabled,
          maxIterations: 4,
          denoise: true
        }
      };
    }

    function applyMolstarRenderStyle(plugin, options = {}) {
      if (!plugin?.canvas3d?.setProps) {
        return false;
      }
      try {
        plugin.canvas3d.setProps(flipbookCanvasProps(options));
        return true;
      } catch (error) {
        console.warn("Molstar Flipbook render style could not be applied.", error);
        return false;
      }
    }

    async function createViewer() {
      disposeViewer();
      stageState.viewer = await molstar.Viewer.create("molstarViewport", {
        layoutIsExpanded: false,
        layoutShowControls: false,
        layoutShowRemoteState: false,
        layoutShowSequence: false,
        layoutShowLog: false,
        layoutShowLeftPanel: false,
        viewportShowExpand: true,
        viewportShowSelectionMode: false,
        viewportShowAnimation: false,
        pdbProvider: "rcsb",
        emdbProvider: "rcsb"
      });
      applyMolstarRenderStyle(stageState.viewer.plugin);
      updateDiagnostics();
      return stageState.viewer;
    }

    function structureStats(pdb) {
      const stats = {
        minX: Infinity,
        maxX: -Infinity,
        minY: Infinity,
        maxY: -Infinity,
        minZ: Infinity,
        maxZ: -Infinity,
        sumX: 0,
        sumY: 0,
        sumZ: 0,
        count: 0
      };
      pdb.split(/\\r?\\n/).forEach((line) => {
        if (!line.startsWith("ATOM") && !line.startsWith("HETATM")) {
          return;
        }
        const x = Number(line.slice(30, 38));
        const y = Number(line.slice(38, 46));
        const z = Number(line.slice(46, 54));
        if (!Number.isFinite(x) || !Number.isFinite(y) || !Number.isFinite(z)) {
          return;
        }
        stats.minX = Math.min(stats.minX, x);
        stats.maxX = Math.max(stats.maxX, x);
        stats.minY = Math.min(stats.minY, y);
        stats.maxY = Math.max(stats.maxY, y);
        stats.minZ = Math.min(stats.minZ, z);
        stats.maxZ = Math.max(stats.maxZ, z);
        stats.sumX += x;
        stats.sumY += y;
        stats.sumZ += z;
        stats.count += 1;
      });
      if (!stats.count) {
        return {
          minX: 0,
          maxX: 28,
          minY: 0,
          maxY: 22,
          minZ: 0,
          maxZ: 22,
          width: 28,
          height: 22,
          depth: 22,
          center: { x: 14, y: 11, z: 11 }
        };
      }
      return {
        ...stats,
        width: Math.max(1, stats.maxX - stats.minX),
        height: Math.max(1, stats.maxY - stats.minY),
        depth: Math.max(1, stats.maxZ - stats.minZ),
        center: {
          x: stats.sumX / stats.count,
          y: stats.sumY / stats.count,
          z: stats.sumZ / stats.count
        }
      };
    }

    function coordinateBounds(pdb) {
      const stats = structureStats(pdb);
      return {
        minX: stats.minX,
        maxX: stats.maxX,
        minY: stats.minY,
        maxY: stats.maxY,
        width: stats.width,
        height: stats.height
      };
    }

    function estimatedPuttyRadius() {
      return Math.max(0.5, wormRadiusMax());
    }

    function visualEnvelopePadding() {
      return Math.max(16, (estimatedPuttyRadius() * 4) + 8);
    }

    function visualEnvelopeSummary() {
      const extents = REPORT.slices.map((slice) => {
        const stats = structureStats(slice.pdb);
        return Math.max(stats.width, stats.height, stats.depth);
      });
      const maxExtent = Math.max(30, ...extents);
      const padding = visualEnvelopePadding();
      return {
        maxSliceExtent: Number(maxExtent.toFixed(3)),
        estimatedPuttyRadius: Number(estimatedPuttyRadius().toFixed(3)),
        padding: Number(padding.toFixed(3)),
        slotSize: Number((maxExtent + padding).toFixed(3))
      };
    }

    function visualMappingSummary() {
      const sampleSlice = REPORT.slices[stageState.currentIndex] || REPORT.slices[0];
      const sampleResidues = REPORT.residues.filter((_, index) => {
        return index === 0 || index === Math.floor((REPORT.residues.length - 1) / 2) || index === REPORT.residues.length - 1;
      });
      return {
        mode: REPORT.visualMapping?.mode || "normalized RMSX B-factor to putty radius",
        sourceRmsxDomain: {
          min: REPORT.domain.min,
          max: REPORT.domain.max
        },
        activeRmsxColorDomain: {
          min: Number(colorDomainMin().toFixed(4)),
          max: Number(colorDomainMax().toFixed(4)),
          defaultMin: DEFAULT_COLOR_MIN,
          defaultMax: DEFAULT_COLOR_MAX
        },
        pdbBfactorDomain: [VISUAL_MIN, VISUAL_MAX],
        radiusRange: {
          configuredMin: stageState.radiusMin,
          configuredMax: stageState.radiusMax,
          thicknessScale: stageState.thicknessScale,
          effectiveMin: Number(wormRadiusMin().toFixed(4)),
          effectiveMax: Number(wormRadiusMax().toFixed(4)),
          effectiveSpan: Number(wormRadiusSpan().toFixed(4))
        },
        colorPalette: stageState.paletteName,
        defaultColorPalette: DEFAULT_PALETTE_NAME,
        availableColorPalettes: paletteNames(),
        colorTheme: {
          requestedTheme: REPORT.visualMapping?.colorTheme || "uncertainty",
          molstarUncertaintyReversesColorList: REPORT.visualMapping?.molstarUncertaintyReversesColorList === true,
          flipbookLowToHighColors: currentPaletteColors(),
          sentToMolstarUncertaintyColors: [...currentPaletteColors()].reverse().map((hex) => hex.toUpperCase()),
          effectiveOrder: REPORT.visualMapping?.paletteOrder || "low-to-high"
        },
        colorStops: currentPaletteColors().map((color, index, colors) => ({
          bfactor: Number((VISUAL_MIN + ((VISUAL_MAX - VISUAL_MIN) * (index / Math.max(1, colors.length - 1)))).toFixed(4)),
          color
        })),
        legend: legendSummary(),
        sampleResidues: sampleResidues.map((residue) => {
          const rmsx = residue.values?.[sampleSlice?.rmsxColumn];
          return {
            residue: residue.label,
            slice: sampleSlice?.label || "",
            rmsx,
            normalizedBfactor: Number(visualWeightForResidue(residue, sampleSlice).toFixed(4)),
            radius: Number(visualRadiusForResidue(residue, sampleSlice).toFixed(4)),
            expectedColor: expectedColorForResidue(residue, sampleSlice)
          };
        })
      };
    }

    function degreesToRadians(degrees) {
      return degrees * Math.PI / 180;
    }

    function radiansToDegrees(radians) {
      return radians * 180 / Math.PI;
    }

    function multiplyMatrices(left, right) {
      return [
        [
          left[0][0] * right[0][0] + left[0][1] * right[1][0] + left[0][2] * right[2][0],
          left[0][0] * right[0][1] + left[0][1] * right[1][1] + left[0][2] * right[2][1],
          left[0][0] * right[0][2] + left[0][1] * right[1][2] + left[0][2] * right[2][2]
        ],
        [
          left[1][0] * right[0][0] + left[1][1] * right[1][0] + left[1][2] * right[2][0],
          left[1][0] * right[0][1] + left[1][1] * right[1][1] + left[1][2] * right[2][1],
          left[1][0] * right[0][2] + left[1][1] * right[1][2] + left[1][2] * right[2][2]
        ],
        [
          left[2][0] * right[0][0] + left[2][1] * right[1][0] + left[2][2] * right[2][0],
          left[2][0] * right[0][1] + left[2][1] * right[1][1] + left[2][2] * right[2][1],
          left[2][0] * right[0][2] + left[2][1] * right[1][2] + left[2][2] * right[2][2]
        ]
      ];
    }

    function rotationMatrixFromEuler(rotation) {
      const x = degreesToRadians(rotation.x);
      const y = degreesToRadians(rotation.y);
      const z = degreesToRadians(rotation.z);
      const sx = Math.sin(x);
      const cx = Math.cos(x);
      const sy = Math.sin(y);
      const cy = Math.cos(y);
      const sz = Math.sin(z);
      const cz = Math.cos(z);
      const rx = [[1, 0, 0], [0, cx, -sx], [0, sx, cx]];
      const ry = [[cy, 0, sy], [0, 1, 0], [-sy, 0, cy]];
      const rz = [[cz, -sz, 0], [sz, cz, 0], [0, 0, 1]];
      return multiplyMatrices(rz, multiplyMatrices(ry, rx));
    }

    function identityRotationMatrix() {
      return [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
    }

    function cloneRotationMatrix(matrix) {
      return (matrix || identityRotationMatrix()).map((row) => row.slice());
    }

    function eulerFromRotationMatrix(matrix) {
      const m = matrix || identityRotationMatrix();
      const y = Math.asin(clamp(-m[2][0], -1, 1));
      const cy = Math.cos(y);
      let x = 0;
      let z = 0;
      if (Math.abs(cy) > 0.000001) {
        x = Math.atan2(m[2][1], m[2][2]);
        z = Math.atan2(m[1][0], m[0][0]);
      } else {
        z = Math.atan2(-m[0][1], m[1][1]);
      }
      return {
        x: radiansToDegrees(x),
        y: radiansToDegrees(y),
        z: radiansToDegrees(z)
      };
    }

    function setRotationMatrix(matrix, options = {}) {
      stageState.rotationMatrix = cloneRotationMatrix(matrix);
      if (options.updateEuler !== false) {
        const euler = eulerFromRotationMatrix(stageState.rotationMatrix);
        stageState.rotation = {
          x: wrapAngle(euler.x),
          y: wrapAngle(euler.y),
          z: wrapAngle(euler.z)
        };
      }
    }

    function syncRotationMatrixFromEuler() {
      stageState.rotationMatrix = rotationMatrixFromEuler(stageState.rotation);
    }

    function currentRotationMatrix() {
      if (!stageState.rotationMatrix) {
        syncRotationMatrixFromEuler();
      }
      return stageState.rotationMatrix;
    }

    function vectorFromArray(value, fallback) {
      const vector = {
        x: Number(value?.[0]),
        y: Number(value?.[1]),
        z: Number(value?.[2])
      };
      if (Number.isFinite(vector.x) && Number.isFinite(vector.y) && Number.isFinite(vector.z)) {
        return vector;
      }
      return { ...fallback };
    }

    function subtractVectors(left, right) {
      return {
        x: left.x - right.x,
        y: left.y - right.y,
        z: left.z - right.z
      };
    }

    function crossVectors(left, right) {
      return {
        x: (left.y * right.z) - (left.z * right.y),
        y: (left.z * right.x) - (left.x * right.z),
        z: (left.x * right.y) - (left.y * right.x)
      };
    }

    function normalizeVector(vector, fallback) {
      const length = Math.sqrt((vector.x * vector.x) + (vector.y * vector.y) + (vector.z * vector.z));
      if (!Number.isFinite(length) || length < 0.000001) {
        return { ...fallback };
      }
      return {
        x: vector.x / length,
        y: vector.y / length,
        z: vector.z / length
      };
    }

    function cameraSnapshot() {
      const camera = stageState.viewer?.plugin?.canvas3d?.camera;
      try {
        return camera?.getSnapshot?.() || camera?.state || null;
      } catch (error) {
        console.debug("Molstar camera snapshot unavailable for local rotation.", error);
        return null;
      }
    }

    function currentScreenRotationAxes() {
      const snapshot = cameraSnapshot();
      const position = vectorFromArray(snapshot?.position, { x: 0, y: 0, z: 100 });
      const target = vectorFromArray(snapshot?.target, { x: 0, y: 0, z: 0 });
      const view = normalizeVector(subtractVectors(target, position), { x: 0, y: 0, z: -1 });
      let up = normalizeVector(vectorFromArray(snapshot?.up, { x: 0, y: 1, z: 0 }), { x: 0, y: 1, z: 0 });
      const right = normalizeVector(crossVectors(view, up), { x: 1, y: 0, z: 0 });
      up = normalizeVector(crossVectors(right, view), up);
      return { right, up, view };
    }

    function axisAngleRotationMatrix(axis, degrees) {
      const unit = normalizeVector(axis, { x: 0, y: 1, z: 0 });
      const angle = degreesToRadians(degrees);
      const c = Math.cos(angle);
      const s = Math.sin(angle);
      const t = 1 - c;
      const x = unit.x;
      const y = unit.y;
      const z = unit.z;
      return [
        [(t * x * x) + c, (t * x * y) - (s * z), (t * x * z) + (s * y)],
        [(t * y * x) + (s * z), (t * y * y) + c, (t * y * z) - (s * x)],
        [(t * z * x) - (s * y), (t * z * y) + (s * x), (t * z * z) + c]
      ];
    }

    function rotationDeltaMatrixForScreenDrag(dx, dy, axes) {
      const horizontal = axisAngleRotationMatrix(axes.up, dx * stageState.rotationSensitivity);
      const vertical = axisAngleRotationMatrix(axes.right, dy * stageState.rotationSensitivity);
      return multiplyMatrices(vertical, horizontal);
    }

    function roundedVector(vector) {
      return {
        x: Number(vector.x.toFixed(4)),
        y: Number(vector.y.toFixed(4)),
        z: Number(vector.z.toFixed(4))
      };
    }

    function applyScreenRotationDrag(dx, dy, axes) {
      const delta = rotationDeltaMatrixForScreenDrag(dx, dy, axes);
      setRotationMatrix(multiplyMatrices(delta, currentRotationMatrix()));
      stageState.lastRotationGesture = {
        mode: "screen-axis delta, per-slice local pivot",
        dx: Number(dx.toFixed(2)),
        dy: Number(dy.toFixed(2)),
        axes: {
          right: roundedVector(axes.right),
          up: roundedVector(axes.up),
          view: roundedVector(axes.view)
        }
      };
    }

    function applyLocalRotation(matrix, center, x, y, z) {
      const localX = x - center.x;
      const localY = y - center.y;
      const localZ = z - center.z;
      return {
        x: center.x + (matrix[0][0] * localX) + (matrix[0][1] * localY) + (matrix[0][2] * localZ),
        y: center.y + (matrix[1][0] * localX) + (matrix[1][1] * localY) + (matrix[1][2] * localZ),
        z: center.z + (matrix[2][0] * localX) + (matrix[2][1] * localY) + (matrix[2][2] * localZ)
      };
    }

    function layoutAnchorCenter() {
      return structureStats(REPORT.slices[0]?.pdb || "").center;
    }

    function layoutTargetCenter(index) {
      const anchor = layoutAnchorCenter();
      const offset = tilePosition(index);
      return {
        x: anchor.x + offset.x,
        y: anchor.y + offset.y,
        z: anchor.z
      };
    }

    function sceneTransformForSlice(slice, index) {
      const stats = structureStats(slice.pdb);
      const rotationMatrix = currentRotationMatrix();
      const center = stats.center;
      const targetCenter = layoutTargetCenter(index);
      const rotatedCenter = {
        x: (rotationMatrix[0][0] * center.x) + (rotationMatrix[0][1] * center.y) + (rotationMatrix[0][2] * center.z),
        y: (rotationMatrix[1][0] * center.x) + (rotationMatrix[1][1] * center.y) + (rotationMatrix[1][2] * center.z),
        z: (rotationMatrix[2][0] * center.x) + (rotationMatrix[2][1] * center.y) + (rotationMatrix[2][2] * center.z)
      };
      const translation = {
        x: targetCenter.x - rotatedCenter.x,
        y: targetCenter.y - rotatedCenter.y,
        z: targetCenter.z - rotatedCenter.z
      };
      return { rotationMatrix, center, targetCenter, translation };
    }

    function applySceneTransform(transform, x, y, z) {
      const r = transform.rotationMatrix;
      const t = transform.translation;
      return {
        x: (r[0][0] * x) + (r[0][1] * y) + (r[0][2] * z) + t.x,
        y: (r[1][0] * x) + (r[1][1] * y) + (r[1][2] * z) + t.y,
        z: (r[2][0] * x) + (r[2][1] * y) + (r[2][2] * z) + t.z
      };
    }

    function sceneBoundsForEntry(entry) {
      const stats = structureStats(entry.slice.pdb);
      const transform = sceneTransformForSlice(entry.slice, entry.sceneIndex);
      const visualRadius = estimatedPuttyRadius();
      const corners = [
        [stats.minX, stats.minY, stats.minZ],
        [stats.minX, stats.minY, stats.maxZ],
        [stats.minX, stats.maxY, stats.minZ],
        [stats.minX, stats.maxY, stats.maxZ],
        [stats.maxX, stats.minY, stats.minZ],
        [stats.maxX, stats.minY, stats.maxZ],
        [stats.maxX, stats.maxY, stats.minZ],
        [stats.maxX, stats.maxY, stats.maxZ]
      ].map(([x, y, z]) => applySceneTransform(transform, x, y, z));
      const bounds = corners.reduce((acc, point) => ({
        minX: Math.min(acc.minX, point.x),
        maxX: Math.max(acc.maxX, point.x),
        minY: Math.min(acc.minY, point.y),
        maxY: Math.max(acc.maxY, point.y),
        minZ: Math.min(acc.minZ, point.z),
        maxZ: Math.max(acc.maxZ, point.z)
      }), {
        minX: Infinity,
        maxX: -Infinity,
        minY: Infinity,
        maxY: -Infinity,
        minZ: Infinity,
        maxZ: -Infinity
      });
      return {
        minX: bounds.minX - visualRadius,
        maxX: bounds.maxX + visualRadius,
        minY: bounds.minY - visualRadius,
        maxY: bounds.maxY + visualRadius,
        minZ: bounds.minZ - visualRadius,
        maxZ: bounds.maxZ + visualRadius
      };
    }

    function sceneFocusSphere(entries = focusSliceEntries()) {
      if (!entries.length) {
        return null;
      }
      const bounds = entries.map(sceneBoundsForEntry).reduce((acc, bounds) => ({
        minX: Math.min(acc.minX, bounds.minX),
        maxX: Math.max(acc.maxX, bounds.maxX),
        minY: Math.min(acc.minY, bounds.minY),
        maxY: Math.max(acc.maxY, bounds.maxY),
        minZ: Math.min(acc.minZ, bounds.minZ),
        maxZ: Math.max(acc.maxZ, bounds.maxZ)
      }), {
        minX: Infinity,
        maxX: -Infinity,
        minY: Infinity,
        maxY: -Infinity,
        minZ: Infinity,
        maxZ: -Infinity
      });
      if (!Number.isFinite(bounds.minX) || !Number.isFinite(bounds.maxX)) {
        return null;
      }
      const center = [
        (bounds.minX + bounds.maxX) / 2,
        (bounds.minY + bounds.maxY) / 2,
        (bounds.minZ + bounds.maxZ) / 2
      ];
      const dx = bounds.maxX - bounds.minX;
      const dy = bounds.maxY - bounds.minY;
      const dz = bounds.maxZ - bounds.minZ;
      return {
        center,
        radius: Math.max(1, Math.sqrt((dx * dx) + (dy * dy) + (dz * dz)) / 2)
      };
    }

    function roundedCoordinate(value) {
      return Number(value.toFixed(4));
    }

    function roundedPoint(point) {
      return {
        x: roundedCoordinate(point.x),
        y: roundedCoordinate(point.y),
        z: roundedCoordinate(point.z)
      };
    }

    function centerAlignmentSummary() {
      const anchor = layoutAnchorCenter();
      const entries = visibleSliceEntries();
      return {
        mode: "rotate around each slice center, then place center on shared Flipbook slot",
        anchor: roundedPoint(anchor),
        entries: entries.map((entry) => {
          const stats = structureStats(entry.slice.pdb);
          const transform = sceneTransformForSlice(entry.slice, entry.sceneIndex);
          const transformedCenter = applySceneTransform(transform, stats.center.x, stats.center.y, stats.center.z);
          return {
            sourceIndex: entry.sourceIndex,
            label: entry.slice.label,
            sourceCenter: roundedPoint(stats.center),
            targetCenter: roundedPoint(transform.targetCenter),
            transformedCenter: roundedPoint(transformedCenter),
            delta: roundedPoint({
              x: transformedCenter.x - transform.targetCenter.x,
              y: transformedCenter.y - transform.targetCenter.y,
              z: transformedCenter.z - transform.targetCenter.z
            })
          };
        })
      };
    }

    function molstarMat4FromTransform(transform) {
      const r = transform.rotationMatrix;
      const t = transform.translation;
      return [
        r[0][0], r[1][0], r[2][0], 0,
        r[0][1], r[1][1], r[2][1], 0,
        r[0][2], r[1][2], r[2][2], 0,
        t.x, t.y, t.z, 1
      ];
    }

    function matrixTransformParamsForSlice(slice, index) {
      return {
        transform: {
          name: "matrix",
          params: {
            data: molstarMat4FromTransform(sceneTransformForSlice(slice, index)),
            transpose: false
          }
        }
      };
    }

    function sceneSpacing() {
      const envelope = visualEnvelopeSummary();
      const slotSize = envelope.slotSize;
      return {
        x: slotSize * stageState.spacingFactor,
        y: Math.max(30, slotSize * 0.82) * stageState.spacingFactor
      };
    }

    function tilePosition(index) {
      if (stageState.layout === "overlay" || stageState.layout === "flip") {
        return { x: 0, y: 0 };
      }
      const columns = clamp(Math.round(stageState.tileColumns), 1, Math.max(1, REPORT.slices.length));
      const rows = Math.ceil(REPORT.slices.length / columns);
      const row = Math.floor(index / columns);
      const column = index % columns;
      const rowLength = row === rows - 1 ? REPORT.slices.length - (row * columns) : columns;
      const spacing = sceneSpacing();
      return {
        x: (column - ((rowLength - 1) / 2)) * spacing.x,
        y: (((rows - 1) / 2) - row) * spacing.y
      };
    }

    function visibleSliceEntries() {
      return REPORT.slices.map((slice, index) => ({
        slice,
        sourceIndex: index,
        sceneIndex: index
      }));
    }

    function focusSliceEntries() {
      const entries = visibleSliceEntries();
      if (stageState.layout === "flip") {
        return entries.filter((entry) => entry.sourceIndex === stageState.currentIndex);
      }
      return entries.filter((entry) => isSliceVisible(entry.sourceIndex));
    }

    function transformedPdbForSlice(slice, index, maskMode = "all", options = {}) {
      const shouldApplySceneTransform = options.applySceneTransform ?? true;
      const sceneTransform = sceneTransformForSlice(slice, index);
      const weights = {};
      REPORT.residues.forEach((residue) => {
        const weight = clamp(visualWeightForResidue(residue, slice), VISUAL_MIN, VISUAL_MAX);
        weights[residue.id] = weight;
        weights[residue.key] = weight;
      });
      let atomCount = 0;
      const lines = slice.pdb.split(/\\r?\\n/).map((line) => {
        if (!line.startsWith("ATOM") && !line.startsWith("HETATM")) {
          return line;
        }
        const padded = line.padEnd(80, " ");
        const x = Number(padded.slice(30, 38));
        const y = Number(padded.slice(38, 46));
        const z = Number(padded.slice(46, 54));
        const residueId = padded.slice(22, 26).trim();
        const chainId = padded.slice(21, 22).trim();
        const chainedKey = chainId ? `${chainId}:${residueId}` : residueId;
        const isMasked = isMaskedResidueKey(chainId, residueId);
        if (maskMode === "unmasked" && isMasked) {
          return null;
        }
        if (maskMode === "masked" && !isMasked) {
          return null;
        }
        atomCount += 1;
        const weight = weights[chainedKey] ?? weights[residueId];
        const transformed = shouldApplySceneTransform && Number.isFinite(x) && Number.isFinite(y) && Number.isFinite(z)
          ? applySceneTransform(sceneTransform, x, y, z)
          : null;
        const xText = transformed ? transformed.x.toFixed(3).padStart(8) : padded.slice(30, 38);
        const yText = transformed ? transformed.y.toFixed(3).padStart(8) : padded.slice(38, 46);
        const zText = transformed ? transformed.z.toFixed(3).padStart(8) : padded.slice(46, 54);
        const bfactorText = Number.isFinite(weight)
          ? weight.toFixed(2).padStart(6)
          : padded.slice(60, 66);
        return `${padded.slice(0, 30)}${xText}${yText}${zText}${padded.slice(54, 60)}${bfactorText}${padded.slice(66)}`.trimEnd();
      }).filter((line) => line !== null);
      return { pdb: lines.join("\\n"), atomCount };
    }

    function selectedResidueAtomCountForSlice(slice, residue = selectedResidue(), options = {}) {
      if (!residue) {
        return 0;
      }
      return slice.pdb.split(/\\r?\\n/).filter((line) => {
        if (!line.startsWith("ATOM") && !line.startsWith("HETATM")) {
          return false;
        }
        return atomLineMatchesResidue(line.padEnd(80, " "), residue, options);
      }).length;
    }

    function selectedResidueMarkerMatchMode(slice, residue = selectedResidue()) {
      return selectedResidueAtomCountForSlice(slice, residue, { strictChain: true }) > 0
        ? "chain-and-residue"
        : "residue-id-fallback";
    }

    function selectedResidueMarkerPdbForSlice(slice, index, options = {}) {
      if (!stageState.selectionMarkerEnabled) {
        return { pdb: "", atomCount: 0, residueKey: "" };
      }
      const residue = selectedResidue();
      if (!residue) {
        return { pdb: "", atomCount: 0, residueKey: "" };
      }
      const shouldApplySceneTransform = options.applySceneTransform ?? true;
      const sceneTransform = sceneTransformForSlice(slice, index);
      const weight = clamp(visualWeightForResidue(residue, slice), VISUAL_MIN, VISUAL_MAX);
      const matchMode = selectedResidueMarkerMatchMode(slice, residue);
      const strictChain = matchMode === "chain-and-residue";
      let atomCount = 0;
      const lines = slice.pdb.split(/\\r?\\n/).map((line) => {
        if (!line.startsWith("ATOM") && !line.startsWith("HETATM")) {
          return null;
        }
        const padded = line.padEnd(80, " ");
        if (!atomLineMatchesResidue(padded, residue, { strictChain })) {
          return null;
        }
        const x = Number(padded.slice(30, 38));
        const y = Number(padded.slice(38, 46));
        const z = Number(padded.slice(46, 54));
        atomCount += 1;
        const transformed = shouldApplySceneTransform && Number.isFinite(x) && Number.isFinite(y) && Number.isFinite(z)
          ? applySceneTransform(sceneTransform, x, y, z)
          : null;
        const xText = transformed ? transformed.x.toFixed(3).padStart(8) : padded.slice(30, 38);
        const yText = transformed ? transformed.y.toFixed(3).padStart(8) : padded.slice(38, 46);
        const zText = transformed ? transformed.z.toFixed(3).padStart(8) : padded.slice(46, 54);
        const bfactorText = weight.toFixed(2).padStart(6);
        return `${padded.slice(0, 30)}${xText}${yText}${zText}${padded.slice(54, 60)}${bfactorText}${padded.slice(66)}`.trimEnd();
      }).filter((line) => line !== null);
      if (atomCount) {
        lines.push("END");
      }
      return { pdb: lines.join("\\n"), atomCount, residueKey: residue.key, matchMode };
    }

    function puttyTypeParams(alpha) {
      const params = {
        sizeFactor: 1,
        quality: "high"
      };
      if (alpha < 1) {
        params.alpha = alpha;
      }
      return params;
    }

    function cartoonTypeParams(alpha) {
      const params = {
        aspectRatio: 1.2,
        sizeFactor: Math.max(0.22, wormRadiusMax() / 2.7),
        quality: "high"
      };
      if (alpha < 1) {
        params.alpha = alpha;
      }
      return params;
    }

    function markerTypeParams() {
      const sizeFactor = Number(REPORT.selectedResidueMarker?.sizeFactor ?? 0.36);
      const alpha = Number(REPORT.selectedResidueMarker?.alpha ?? 0.86);
      return {
        sizeFactor: Number.isFinite(sizeFactor) ? sizeFactor : 0.36,
        quality: "high",
        alpha: Number.isFinite(alpha) ? clamp(alpha, 0.05, 1) : 0.86
      };
    }

    function markerColorParams() {
      return {
        value: hexColorToMolstarNumber(MARKER_COLOR_HEX, 0x111827)
      };
    }

    async function addRmsxRepresentation(plugin, structure, index, options = {}) {
      const alpha = options.alpha ?? 1;
      const tagSuffix = options.tagSuffix || "all";
      try {
        const representation = await plugin.builders.structure.representation.addRepresentation(
          structure,
          {
            type: "putty",
            typeParams: puttyTypeParams(alpha),
            color: "uncertainty",
            colorParams: uncertaintyColorParams(),
            size: "uncertainty",
            sizeParams: {
              bfactorFactor: wormRadiusSpan(),
              rmsfFactor: 0,
              baseSize: wormRadiusMin()
            }
          },
          { tag: `flipbook-putty-${index}-${tagSuffix}` }
        );
        stageState.representationMode = "putty";
        return representation;
      } catch (error) {
        console.warn("Molstar putty representation failed, falling back to cartoon.", error);
        stageState.representationMode = "cartoon";
        return plugin.builders.structure.representation.addRepresentation(
          structure,
          {
            type: "cartoon",
            typeParams: cartoonTypeParams(alpha),
            color: "uncertainty",
            colorParams: uncertaintyColorParams(),
            size: "uncertainty",
            sizeParams: {
              bfactorFactor: wormRadiusSpan(),
              rmsfFactor: 0,
              baseSize: wormRadiusMin()
            }
          },
          { tag: `flipbook-cartoon-${index}-${tagSuffix}` }
        );
      }
    }

    async function addSelectedResidueMarkerRepresentation(plugin, structure, entry) {
      const markerType = REPORT.selectedResidueMarker?.type || "spacefill";
      const markerParams = {
        type: markerType,
        typeParams: markerTypeParams(),
        color: "uniform",
        colorParams: markerColorParams()
      };
      try {
        return await plugin.builders.structure.representation.addRepresentation(
          structure,
          markerParams,
          { tag: `flipbook-selected-residue-marker-${entry.sourceIndex}` }
        );
      } catch (error) {
        console.warn("Molstar selected-residue marker failed; falling back to ball-and-stick.", error);
        return plugin.builders.structure.representation.addRepresentation(
          structure,
          {
            ...markerParams,
            type: "ball-and-stick"
          },
          { tag: `flipbook-selected-residue-marker-fallback-${entry.sourceIndex}` }
        );
      }
    }

    function transformCandidate(label, resolver) {
      try {
        const value = resolver();
        return value ? { label, value } : null;
      } catch (error) {
        console.debug(`Molstar transform candidate ${label} unavailable.`, error);
        return null;
      }
    }

    function transformerFromRegistry(registry) {
      if (!registry) {
        return null;
      }
      if (typeof registry.get === "function") {
        const directMatch = registry.get("transform-structure-conformation")
          || registry.get("TransformStructureConformation")
          || registry.get("model.transform-structure-conformation")
          || null;
        if (directMatch) {
          return directMatch;
        }
      }
      const propertyMatch = registry["transform-structure-conformation"]
        || registry.TransformStructureConformation
        || registry["model.transform-structure-conformation"]
        || null;
      if (propertyMatch) {
        return propertyMatch;
      }
      return findTransformerInRegistry(registry);
    }

    function transformerDescriptor(candidate) {
      if (!candidate) {
        return null;
      }
      let text = "";
      try {
        text = String(candidate).slice(0, 120);
      } catch (error) {
        text = `<string unavailable: ${error.message}>`;
      }
      return {
        text,
        id: candidate.id || null,
        name: candidate.name || null,
        definitionName: candidate.definition?.name || null,
        displayName: candidate.definition?.display?.name || candidate.display?.name || null,
        keys: objectKeys(candidate).slice(0, 12)
      };
    }

    function transformerSearchText(candidate) {
      const descriptor = transformerDescriptor(candidate);
      if (!descriptor) {
        return "";
      }
      return [
        descriptor.text,
        descriptor.id,
        descriptor.name,
        descriptor.definitionName,
        descriptor.displayName,
        descriptor.keys.join(" ")
      ].filter(Boolean).join(" ").toLowerCase();
    }

    function isStructureConformationTransformer(candidate) {
      const haystack = transformerSearchText(candidate);
      return haystack.includes("transform-structure-conformation")
        || haystack.includes("transformstructureconformation")
        || haystack.includes("transform conformation");
    }

    function findTransformerInRegistry(registry) {
      let found = null;
      const testCandidate = (candidate) => {
        if (!found && isStructureConformationTransformer(candidate)) {
          found = candidate;
        }
      };
      try {
        if (typeof registry.forEach === "function") {
          registry.forEach((value, key) => {
            testCandidate(value);
            testCandidate(key);
          });
        }
      } catch (error) {
        console.debug("Molstar registry forEach probe failed.", error);
      }
      try {
        if (!found && typeof registry.valueSeq === "function") {
          registry.valueSeq().toArray().forEach(testCandidate);
        }
      } catch (error) {
        console.debug("Molstar registry valueSeq probe failed.", error);
      }
      try {
        if (!found && typeof registry.keySeq === "function") {
          registry.keySeq().toArray().forEach(testCandidate);
        }
      } catch (error) {
        console.debug("Molstar registry keySeq probe failed.", error);
      }
      return found;
    }

    function structureTransformCandidate(plugin) {
      const candidates = [
        transformCandidate("window.FLIPBOOK_MOLSTAR_STATE_TRANSFORMS.Model.TransformStructureConformation", () => window.FLIPBOOK_MOLSTAR_STATE_TRANSFORMS?.Model?.TransformStructureConformation),
        transformCandidate("window.molstar.lib.StateTransforms.Model.TransformStructureConformation", () => window.molstar?.lib?.StateTransforms?.Model?.TransformStructureConformation),
        transformCandidate("window.molstar.StateTransforms.Model.TransformStructureConformation", () => window.molstar?.StateTransforms?.Model?.TransformStructureConformation),
        transformCandidate("plugin.state.data.transforms", () => transformerFromRegistry(plugin?.state?.data?.transforms)),
        transformCandidate("plugin.state.data.registry.transforms", () => transformerFromRegistry(plugin?.state?.data?.registry?.transforms)),
        transformCandidate("plugin.state.transforms", () => transformerFromRegistry(plugin?.state?.transforms)),
        transformCandidate("plugin.transforms", () => transformerFromRegistry(plugin?.transforms))
      ].filter(Boolean);
      return candidates.find((candidate) => candidate.value)?.value || null;
    }

    function structureTransformApi(plugin) {
      if (!plugin?.build) {
        return null;
      }
      return structureTransformCandidate(plugin);
    }

    async function loadMolstarStateTransforms() {
      if (window.FLIPBOOK_MOLSTAR_STATE_TRANSFORMS?.Model?.TransformStructureConformation) {
        stageState.stateTransformsModuleStatus = "loaded";
        updateDiagnostics();
        return true;
      }
      const moduleUrls = REPORT.molstar?.stateTransformsUrls || [REPORT.molstar?.stateTransformsUrl].filter(Boolean);
      if (!moduleUrls.length) {
        stageState.stateTransformsModuleStatus = "missing-url";
        stageState.stateTransformsModuleUrl = "";
        updateDiagnostics();
        return false;
      }
      const failures = [];
      for (const moduleUrl of moduleUrls) {
        stageState.stateTransformsModuleStatus = "loading";
        stageState.stateTransformsModuleUrl = moduleUrl;
        updateDiagnostics();
        try {
          const module = await import(moduleUrl);
          window.FLIPBOOK_MOLSTAR_STATE_TRANSFORMS = module.StateTransforms;
          stageState.stateTransformsModuleStatus = "loaded";
          updateDiagnostics();
          return Boolean(window.FLIPBOOK_MOLSTAR_STATE_TRANSFORMS?.Model?.TransformStructureConformation);
        } catch (error) {
          failures.push(`${moduleUrl}: ${error.message}`);
          console.warn(`Molstar StateTransforms module could not be imported from ${moduleUrl}.`, error);
        }
      }
      stageState.stateTransformsModuleStatus = `unavailable: ${failures.join(" | ")}`;
      updateDiagnostics();
      return false;
    }

    function canUseStructureTransforms(plugin) {
      return stageState.preferStructureTransforms
        && !stageState.forcePdbCoordinateTransforms
        && Boolean(structureTransformApi(plugin));
    }

    function canUseRepresentationTransforms() {
      return stageState.preferRepresentationTransforms && !stageState.forcePdbCoordinateTransforms;
    }

    async function insertSliceTransform(plugin, structure, entry, tagSuffix) {
      const transform = structureTransformApi(plugin);
      if (!transform) {
        return null;
      }
      try {
        return await plugin.build()
          .to(structure)
          .insert(transform, matrixTransformParamsForSlice(entry.slice, entry.sceneIndex), {
            tags: [`flipbook-scene-transform-${entry.sourceIndex}-${tagSuffix}`]
          })
          .commit();
      } catch (error) {
        const wrapped = new Error(`FLIPBOOK_MOLSTAR_STRUCTURE_TRANSFORM_FAILED: ${error.message}`);
        wrapped.cause = error;
        throw wrapped;
      }
    }

    async function updateRecordTransform(plugin, record) {
      if (!record.transformCell) {
        return;
      }
      await plugin.build()
        .to(record.transformCell)
        .update(matrixTransformParamsForSlice(record.slice, record.sceneIndex))
        .commit();
    }

    function hasStructureTransformRecords() {
      return stageState.records.some((record) => record?.transformCell);
    }

    function hasRepresentationTransformRecords() {
      return stageState.records.some((record) => record?.representationTransform);
    }

    function hasLiveTransformRecords() {
      return hasStructureTransformRecords() || hasRepresentationTransformRecords();
    }

    function representationObject(record) {
      return record?.representation?.cell?.obj?.data?.repr
        || record?.representation?.obj?.data?.repr
        || record?.representation?.data?.repr
        || null;
    }

    function flushMolstarTransformDraw(plugin, representations, options = {}) {
      const canvas = plugin?.canvas3d;
      const summary = {
        mode: options.fast ? "interactive" : "settled",
        visibleOnly: options.visibleOnly === true,
        representationUpdates: 0,
        canvasUpdates: 0,
        commit: false,
        requestDraw: false
      };
      if (!canvas) {
        stageState.lastTransformFlush = { ...summary, error: "missing-canvas3d" };
        return summary;
      }
      representations.forEach((repr) => {
        if (!repr) {
          return;
        }
        summary.representationUpdates += 1;
        if (typeof canvas.update === "function") {
          try {
            canvas.update(repr, options.keepBoundingSphere === true);
            summary.canvasUpdates += 1;
          } catch (error) {
            console.debug("Molstar canvas update after representation transform failed.", error);
          }
        }
      });
      if (typeof canvas.commit === "function") {
        try {
          canvas.commit(options.synchronous === true);
          summary.commit = true;
        } catch (error) {
          console.debug("Molstar canvas commit after representation transform failed.", error);
        }
      }
      if (typeof canvas.requestDraw === "function") {
        try {
          canvas.requestDraw();
          summary.requestDraw = true;
        } catch (error) {
          console.debug("Molstar canvas draw request after representation transform failed.", error);
        }
      }
      stageState.lastTransformFlush = summary;
      return summary;
    }

    function applyRecordRepresentationTransform(record) {
      const repr = representationObject(record);
      if (!repr?.setState) {
        return false;
      }
      repr.setState({ transform: molstarMat4FromTransform(sceneTransformForSlice(record.slice, record.sceneIndex)) });
      return true;
    }

    function updateRepresentationTransformRecords(plugin, options = {}) {
      let updated = 0;
      const representations = [];
      stageState.records.forEach((record) => {
        if (options.visibleOnly && record && !isSliceVisible(record.sourceIndex)) {
          return;
        }
        if (record?.representationTransform && applyRecordRepresentationTransform(record)) {
          updated += 1;
          representations.push(representationObject(record));
        }
      });
      if (updated) {
        flushMolstarTransformDraw(plugin, representations, {
          fast: options.fast === true,
          synchronous: options.fast === true,
          keepBoundingSphere: options.fast === true,
          visibleOnly: options.visibleOnly === true
        });
      }
      return updated;
    }

    function setRecordVisibility(record, visible) {
      const repr = representationObject(record);
      if (!repr?.setState) {
        return null;
      }
      repr.setState({ visible, pickable: visible });
      return repr;
    }

    function applySliceVisibilityToRecords(options = {}) {
      const plugin = stageState.viewer?.plugin;
      const representations = [];
      stageState.records.forEach((record) => {
        if (!record) {
          return;
        }
        const repr = setRecordVisibility(record, isSliceVisible(record.sourceIndex));
        if (repr) {
          representations.push(repr);
        }
      });
      const flush = representations.length
        ? flushMolstarTransformDraw(plugin, representations, {
            fast: options.fast !== false,
            synchronous: options.fast !== false,
            keepBoundingSphere: true
          })
        : null;
      stageState.lastVisibilityUpdate = {
        visibleSliceIndexes: visibleSliceIndexList(),
        hiddenSliceIndexes: hiddenSliceIndexList(),
        representationUpdates: representations.length,
        flush
      };
      return representations.length;
    }

    function loadedSliceIndexList() {
      return [...new Set(stageState.records
        .filter((record) => record && Number.isInteger(record.sourceIndex))
        .map((record) => record.sourceIndex))]
        .sort((left, right) => left - right);
    }

    function selectedResidueMarkerRecords() {
      return stageState.records.filter((record) => record?.kind === "selected-residue-marker");
    }

    function selectedResidueMarkerSummary() {
      const residue = selectedResidue();
      const markerRecords = selectedResidueMarkerRecords();
      const expectedSliceIndexes = stageState.selectionMarkerEnabled
        ? REPORT.slices
            .map((slice, index) => selectedResidueAtomCountForSlice(slice, residue) > 0 ? index : null)
            .filter((index) => index !== null)
        : [];
      return {
        enabled: stageState.selectionMarkerEnabled,
        residueKey: residue?.key || "",
        residueLabel: residue?.label || "",
        color: MARKER_COLOR_HEX,
        type: REPORT.selectedResidueMarker?.type || "spacefill",
        atomCountsBySlice: REPORT.slices.map((slice, index) => ({
          sourceIndex: index,
          atomCount: selectedResidueAtomCountForSlice(slice, residue),
          matchMode: selectedResidueMarkerMatchMode(slice, residue)
        })),
        expectedSliceIndexes,
        markerSliceIndexes: [...new Set(markerRecords.map((record) => record.sourceIndex))].sort((left, right) => left - right),
        visibleMarkerSliceIndexes: [...new Set(markerRecords
          .filter((record) => isSliceVisible(record.sourceIndex))
          .map((record) => record.sourceIndex))]
          .sort((left, right) => left - right),
        markerRecords: markerRecords.length,
        visibleMarkerRecords: markerRecords.filter((record) => isSliceVisible(record.sourceIndex)).length,
        markerAtomCounts: markerRecords.map((record) => ({
          sourceIndex: record.sourceIndex,
          atomCount: record.atomCount || 0,
          matchMode: record.matchMode || ""
        }))
      };
    }

    function objectKeys(value) {
      try {
        return value && typeof value === "object" ? Object.keys(value).sort() : [];
      } catch (error) {
        return [`<keys unavailable: ${error.message}>`];
      }
    }

    function registrySummary(registry) {
      if (!registry) {
        return { present: false, keys: [] };
      }
      const sampleValues = [];
      const matches = [];
      const addSample = (candidate) => {
        const descriptor = transformerDescriptor(candidate);
        if (!descriptor) {
          return;
        }
        if (sampleValues.length < 20) {
          sampleValues.push(descriptor);
        }
        if (isStructureConformationTransformer(candidate) && matches.length < 10) {
          matches.push(descriptor);
        }
      };
      try {
        if (typeof registry.forEach === "function") {
          registry.forEach((value, key) => {
            addSample(value);
            addSample(key);
          });
        }
      } catch (error) {
        sampleValues.push({ text: `<forEach unavailable: ${error.message}>` });
      }
      try {
        if (!sampleValues.length && typeof registry.valueSeq === "function") {
          registry.valueSeq().toArray().slice(0, 20).forEach(addSample);
        }
      } catch (error) {
        sampleValues.push({ text: `<valueSeq unavailable: ${error.message}>` });
      }
      return {
        present: true,
        keys: objectKeys(registry).slice(0, 30),
        protoKeys: objectKeys(Object.getPrototypeOf(registry)).slice(0, 30),
        size: registry.size ?? null,
        hasGet: typeof registry.get === "function",
        hasStructureConformation: Boolean(transformerFromRegistry(registry)),
        sampleValues,
        matches
      };
    }

    function invariantCheck(name, passed, details = {}) {
      return {
        name,
        passed: Boolean(passed),
        details
      };
    }

    function allExpectedSlicesLoaded(loadedIndexes) {
      return REPORT.slices.every((_, index) => loadedIndexes.includes(index));
    }

    function centerAlignmentMaxDelta(alignment) {
      const deltas = (alignment?.entries || []).flatMap((entry) => {
        return [Math.abs(entry.delta?.x || 0), Math.abs(entry.delta?.y || 0), Math.abs(entry.delta?.z || 0)];
      });
      return deltas.length ? Math.max(...deltas) : 0;
    }

    function sameNumberArray(left, right) {
      const leftValues = left || [];
      const rightValues = right || [];
      return leftValues.length === rightValues.length
        && leftValues.every((value, index) => value === rightValues[index]);
    }

    function cleanRenderInvariant(renderStyle) {
      const outlineParam = (URL_PARAMS.get("outline") || "").toLowerCase();
      const outlineDisabled = outlineParam === "0" || outlineParam === "false" || outlineParam === "off";
      const expectedOutline = outlineDisabled ? "off" : "on";
      return renderStyle.multiSampleMode === "off"
        && renderStyle.occlusion === "off"
        && renderStyle.illumination === false
        && renderStyle.cameraFog === "off"
        && Math.abs(Number(renderStyle.ambientIntensity) - 0.78) < 0.0001
        && renderStyle.outline === expectedOutline;
    }

    function flipbookInvariantSummary(summary) {
      const loadedIndexes = summary.visibility.loadedSliceIndexes || [];
      const visibleIndexes = summary.visibility.visibleSliceIndexes || [];
      const maxCenterDelta = centerAlignmentMaxDelta(summary.centerAlignment);
      const sampleResidues = summary.visualMapping.sampleResidues || [];
      const markerSummary = summary.selectedResidueMarker || {};
      const expectedMarkerIndexes = markerSummary.expectedSliceIndexes || [];
      const expectedVisibleMarkerIndexes = expectedMarkerIndexes.filter((index) => visibleIndexes.includes(index));
      const playDelayRangeElement = document.getElementById("playDelayRange");
      const playDelayMin = Number(playDelayRangeElement?.min);
      const playDelayMax = Number(playDelayRangeElement?.max);
      const playDelayStep = Number(playDelayRangeElement?.step);
      const playDelayDefault = Number(playDelayRangeElement?.defaultValue);
      const spacingRangeElement = document.getElementById("spacingRange");
      const spacingMin = Number(spacingRangeElement?.min);
      const spacingDefault = Number(spacingRangeElement?.defaultValue);
      const checks = [
        invariantCheck("molstar-viewer-loaded", summary.hasViewer && Boolean(stageState.viewer?.plugin), {
          hasViewer: summary.hasViewer,
          hasPlugin: Boolean(stageState.viewer?.plugin)
        }),
        invariantCheck("default-presentation-is-tiled", REPORT.presentation?.defaultLayout === "tiled", {
          defaultLayout: REPORT.presentation?.defaultLayout,
          currentLayout: stageState.layout,
          availableLayouts: REPORT.presentation?.availableLayouts || []
        }),
        invariantCheck("all-slices-loaded-in-one-scene", loadedIndexes.length === REPORT.slices.length && allExpectedSlicesLoaded(loadedIndexes), {
          loadedSliceIndexes: loadedIndexes,
          expectedSliceCount: REPORT.slices.length,
          totalRecords: summary.totalRecords
        }),
        invariantCheck("flip-layout-shows-one-active-slice", stageState.layout !== "flip" || (visibleIndexes.length === 1 && visibleIndexes[0] === stageState.currentIndex), {
          layout: stageState.layout,
          currentIndex: stageState.currentIndex,
          visibleSliceIndexes: visibleIndexes,
          hiddenSliceIndexes: summary.visibility.hiddenSliceIndexes
        }),
        invariantCheck("per-slice-centers-stay-on-layout-targets", maxCenterDelta <= 0.01, {
          maxCenterDelta: Number(maxCenterDelta.toFixed(6)),
          mode: summary.centerAlignment.mode
        }),
        invariantCheck("render-style-is-clean-and-stable", cleanRenderInvariant(summary.canvasRenderStyle), {
          outline: summary.canvasRenderStyle.outline,
          ambientIntensity: summary.canvasRenderStyle.ambientIntensity,
          multiSampleMode: summary.canvasRenderStyle.multiSampleMode,
          occlusion: summary.canvasRenderStyle.occlusion,
          illumination: summary.canvasRenderStyle.illumination,
          cameraFog: summary.canvasRenderStyle.cameraFog
        }),
        invariantCheck("rmsx-radius-and-color-mapping-is-explicit", summary.visualMapping.pdbBfactorDomain[0] === 0
          && summary.visualMapping.pdbBfactorDomain[1] === 1
          && summary.visualMapping.radiusRange.effectiveMax > summary.visualMapping.radiusRange.effectiveMin
          && Boolean(document.getElementById("radiusMinNumber"))
          && Boolean(document.getElementById("radiusMaxNumber"))
          && summary.visualMapping.colorTheme.molstarUncertaintyReversesColorList === true
          && sampleResidues.every((residue) => Number.isFinite(residue.radius) && /^#[0-9A-F]{6}$/.test(residue.expectedColor)), {
            pdbBfactorDomain: summary.visualMapping.pdbBfactorDomain,
            radiusRange: summary.visualMapping.radiusRange,
            radiusMinInput: document.getElementById("radiusMinNumber")?.value,
            radiusMaxInput: document.getElementById("radiusMaxNumber")?.value,
            colorTheme: summary.visualMapping.colorTheme.requestedTheme,
            sampleResidues
          }),
        invariantCheck("rmsx-radius-range-is-calibratable", Number.isFinite(summary.visualMapping.radiusRange.configuredMin)
          && Number.isFinite(summary.visualMapping.radiusRange.configuredMax)
          && summary.visualMapping.radiusRange.configuredMax > summary.visualMapping.radiusRange.configuredMin
          && Boolean(document.getElementById("radiusMinNumber"))
          && Boolean(document.getElementById("radiusMaxNumber")), {
          radiusRange: summary.visualMapping.radiusRange,
          radiusMinInput: document.getElementById("radiusMinNumber")?.value,
          radiusMaxInput: document.getElementById("radiusMaxNumber")?.value
        }),
        invariantCheck("rmsx-color-domain-is-calibratable", Number.isFinite(summary.visualMapping.activeRmsxColorDomain.min)
          && Number.isFinite(summary.visualMapping.activeRmsxColorDomain.max)
          && summary.visualMapping.activeRmsxColorDomain.max > summary.visualMapping.activeRmsxColorDomain.min
          && Boolean(document.getElementById("colorMinNumber"))
          && Boolean(document.getElementById("colorMaxNumber")), {
          activeRmsxColorDomain: summary.visualMapping.activeRmsxColorDomain,
          colorMinInput: document.getElementById("colorMinNumber")?.value,
          colorMaxInput: document.getElementById("colorMaxNumber")?.value
        }),
        invariantCheck("rmsx-palette-is-switchable", Boolean(document.getElementById("paletteSelect"))
          && summary.visualMapping.availableColorPalettes.includes("mako")
          && summary.visualMapping.availableColorPalettes.includes("turbo")
          && summary.visualMapping.availableColorPalettes.includes(stageState.paletteName)
          && summary.visualMapping.colorPalette === stageState.paletteName
          && currentPaletteColors().length > 1, {
          paletteSelect: document.getElementById("paletteSelect")?.value,
          activePalette: summary.visualMapping.colorPalette,
          defaultPalette: summary.visualMapping.defaultColorPalette,
          availablePalettes: summary.visualMapping.availableColorPalettes
        }),
        invariantCheck("rmsx-legend-matches-active-visual-domain", summary.visualMapping.legend?.elements?.legend === true
          && summary.visualMapping.legend?.elements?.colorBar === true
          && summary.visualMapping.legend?.elements?.radiusLegend === true
          && numbersClose(summary.visualMapping.legend.activeRmsxColorDomain.min, summary.visualMapping.activeRmsxColorDomain.min)
          && numbersClose(summary.visualMapping.legend.activeRmsxColorDomain.max, summary.visualMapping.activeRmsxColorDomain.max)
          && numbersClose(summary.visualMapping.legend.radiusRange.min, summary.visualMapping.radiusRange.effectiveMin)
          && numbersClose(summary.visualMapping.legend.radiusRange.max, summary.visualMapping.radiusRange.effectiveMax)
          && summary.visualMapping.legend.stops.length === 3
          && summary.visualMapping.legend.stops.every((stop) => /^#[0-9A-F]{6}$/.test(stop.color) && Number.isFinite(stop.radius)), {
          legend: summary.visualMapping.legend,
          activeRmsxColorDomain: summary.visualMapping.activeRmsxColorDomain,
          radiusRange: summary.visualMapping.radiusRange
        }),
        invariantCheck("selected-residue-marker-tracks-visible-slices", markerSummary.enabled !== true
          || (expectedMarkerIndexes.length > 0
            && sameNumberArray(markerSummary.markerSliceIndexes, expectedMarkerIndexes)
            && sameNumberArray(markerSummary.visibleMarkerSliceIndexes, expectedVisibleMarkerIndexes)), {
          selectedResidueMarker: markerSummary,
          visibleSliceIndexes: visibleIndexes,
          expectedVisibleMarkerIndexes
        }),
        invariantCheck("playback-delay-control-is-calibratable", Number.isFinite(playDelayMin)
          && Number.isFinite(playDelayMax)
          && Number.isFinite(playDelayStep)
          && numbersClose(playDelayMin, MIN_PLAY_DELAY_MS)
          && numbersClose(playDelayMax, MAX_PLAY_DELAY_MS)
          && numbersClose(playDelayDefault, DEFAULT_PLAY_DELAY_MS)
          && numbersClose(playDelayStep, PLAY_DELAY_STEP_MS)
          && stageState.playDelayMs >= MIN_PLAY_DELAY_MS
          && stageState.playDelayMs <= MAX_PLAY_DELAY_MS, {
          playDelayRangeMin: playDelayMin,
          playDelayRangeMax: playDelayMax,
          playDelayRangeDefault: playDelayDefault,
          playDelayStep,
          playback: summary.playback
        }),
        invariantCheck("visual-mapping-reset-is-available", Boolean(document.getElementById("resetScaleButton"))
          && typeof resetVisualMapping === "function"
          && Number.isFinite(DEFAULT_THICKNESS_SCALE)
          && Number.isFinite(DEFAULT_RADIUS_MIN)
          && Number.isFinite(DEFAULT_RADIUS_MAX)
          && Number.isFinite(DEFAULT_COLOR_MIN)
          && Number.isFinite(DEFAULT_COLOR_MAX)
          && DEFAULT_RADIUS_MAX > DEFAULT_RADIUS_MIN
          && DEFAULT_COLOR_MAX > DEFAULT_COLOR_MIN, {
          resetScaleButton: Boolean(document.getElementById("resetScaleButton")),
          defaults: {
            thicknessScale: DEFAULT_THICKNESS_SCALE,
            radiusMin: DEFAULT_RADIUS_MIN,
            radiusMax: DEFAULT_RADIUS_MAX,
            colorMin: DEFAULT_COLOR_MIN,
            colorMax: DEFAULT_COLOR_MAX
          }
        }),
        invariantCheck("viewer-url-reflects-current-state", summary.urlState?.synced === true, {
          managedParams: summary.urlState?.managedParams,
          expectedManagedParams: summary.urlState?.expectedManagedParams,
          matches: summary.urlState?.matches
        }),
        invariantCheck("local-rotation-uses-direct-screen-dy", REPORT.rotationModel?.verticalDragSign === "screen dy is applied directly around the current screen-right axis", {
          verticalDragSign: REPORT.rotationModel?.verticalDragSign,
          lastRotationGesture: summary.lastRotationGesture
        }),
        invariantCheck("tight-spacing-control-is-available", Number.isFinite(spacingMin)
          && numbersClose(spacingMin, MIN_TILE_SPACING_FACTOR)
          && numbersClose(spacingDefault, DEFAULT_TILE_SPACING_FACTOR)
          && DEFAULT_TILE_SPACING_FACTOR > spacingMin
          && stageState.spacingFactor >= spacingMin, {
          spacingRangeMin: spacingMin,
          spacingRangeDefault: spacingDefault,
          defaultSpacingFactor: DEFAULT_TILE_SPACING_FACTOR,
          spacingFactor: stageState.spacingFactor,
          sceneSpacing: summary.sceneSpacing
        }),
        invariantCheck("report-payload-is-diagnosed", Number.isFinite(summary.reportPayload?.embeddedPdbBytes)
          && summary.reportPayload.embeddedPdbBytes > 0
          && summary.reportPayload.embeddedPdbCount === REPORT.slices.length
          && Number.isFinite(summary.reportPayload?.estimatedJsonBytes)
          && Number.isFinite(summary.reportPayload?.largeReportWarningBytes), {
          reportPayload: summary.reportPayload
        }),
        invariantCheck("vmd-hotkey-parity-subset-is-available", summary.keyboardShortcuts?.enabled === true
          && summary.keyboardShortcuts.rotationStepDegrees === 5
          && summary.keyboardShortcuts.bindings.includes("u/i")
          && summary.keyboardShortcuts.bindings.includes("[/]")
          && summary.keyboardShortcuts.bindings.includes(",/."), {
          source: summary.keyboardShortcuts?.source,
          rotationStepDegrees: summary.keyboardShortcuts?.rotationStepDegrees,
          spacingStep: summary.keyboardShortcuts?.spacingStep,
          thicknessStep: summary.keyboardShortcuts?.thicknessStep,
          bindings: summary.keyboardShortcuts?.bindings
        })
      ];
      return {
        passed: checks.every((check) => check.passed),
        failed: checks.filter((check) => !check.passed).map((check) => check.name),
        checks
      };
    }

    function molstarCapabilitySummary() {
      const plugin = stageState.viewer?.plugin;
      const dataState = plugin?.state?.data;
      const summary = {
        layout: stageState.layout,
        presentation: REPORT.presentation,
        autoplay: stageState.autoplay,
        playback: {
          active: Boolean(stageState.timer),
          delayMs: stageState.playDelayMs,
          minDelayMs: MIN_PLAY_DELAY_MS,
          maxDelayMs: MAX_PLAY_DELAY_MS,
          defaultDelayMs: DEFAULT_PLAY_DELAY_MS,
          delayStepMs: PLAY_DELAY_STEP_MS,
          delayUrlParam: REPORT.playback?.delayUrlParam || "delayMs"
        },
        geometryMode: stageState.geometryMode,
        forcePdbCoordinateTransforms: stageState.forcePdbCoordinateTransforms,
        stateTransformsModuleStatus: stageState.stateTransformsModuleStatus,
        stateTransformsModuleUrl: stageState.stateTransformsModuleUrl,
        hasViewer: Boolean(window.molstar?.Viewer),
        stateTransformsGlobalAvailable: Boolean(window.FLIPBOOK_MOLSTAR_STATE_TRANSFORMS?.Model?.TransformStructureConformation),
        reportPayload: REPORT.reportPayload,
        molstarRootKeys: objectKeys(window.molstar),
        pluginKeys: objectKeys(plugin).slice(0, 30),
        stateKeys: objectKeys(plugin?.state).slice(0, 30),
        dataStateKeys: objectKeys(dataState).slice(0, 30),
        transformApiAvailable: Boolean(structureTransformApi(plugin)),
        transformRecords: stageState.records.filter((record) => record?.transformCell).length,
        representationTransformRecords: stageState.records.filter((record) => record?.representationTransform).length,
        totalRecords: stageState.records.length,
        rotationModel: REPORT.rotationModel,
        molstarRenderStyle: REPORT.molstarRenderStyle,
        transformFlush: stageState.lastTransformFlush,
        visibility: {
          loadedSliceIndexes: loadedSliceIndexList(),
          flipKeepsAllSlicesLoaded: stageState.layout === "flip" && loadedSliceIndexList().length === REPORT.slices.length,
          visibleSliceIndexes: visibleSliceIndexList(),
          hiddenSliceIndexes: hiddenSliceIndexList(),
          visibleRecords: stageState.records.filter((record) => record && isSliceVisible(record.sourceIndex)).length,
          hiddenRecords: stageState.records.filter((record) => record && !isSliceVisible(record.sourceIndex)).length,
          lastVisibilityUpdate: stageState.lastVisibilityUpdate
        },
        canvasRenderStyle: {
          renderPreset: activeRenderPreset(),
          backgroundColor: plugin?.canvas3d?.props?.renderer?.backgroundColor ?? null,
          ambientIntensity: plugin?.canvas3d?.props?.renderer?.ambientIntensity ?? null,
          outline: plugin?.canvas3d?.props?.postprocessing?.outline?.name ?? null,
          occlusion: plugin?.canvas3d?.props?.postprocessing?.occlusion?.name ?? null,
          antialiasing: plugin?.canvas3d?.props?.postprocessing?.antialiasing?.name ?? null,
          illumination: plugin?.canvas3d?.props?.illumination?.enabled ?? null,
          multiSampleMode: plugin?.canvas3d?.props?.multiSample?.mode ?? null,
          multiSampleLevel: plugin?.canvas3d?.props?.multiSample?.sampleLevel ?? null,
          reduceFlicker: plugin?.canvas3d?.props?.multiSample?.reduceFlicker ?? null,
          reuseOcclusion: plugin?.canvas3d?.props?.multiSample?.reuseOcclusion ?? null,
          dpoitIterations: plugin?.canvas3d?.props?.dpoitIterations ?? null,
          postprocessingEnabled: plugin?.canvas3d?.props?.postprocessing?.enabled ?? null,
          shadow: plugin?.canvas3d?.props?.postprocessing?.shadow?.name ?? null,
          dof: plugin?.canvas3d?.props?.postprocessing?.dof?.name ?? null,
          sharpening: plugin?.canvas3d?.props?.postprocessing?.sharpening?.name ?? null,
          bloom: plugin?.canvas3d?.props?.postprocessing?.bloom?.name ?? null,
          cameraFog: plugin?.canvas3d?.props?.cameraFog?.name ?? null,
          canUpdateRepresentations: typeof plugin?.canvas3d?.update === "function",
          canCommit: typeof plugin?.canvas3d?.commit === "function",
          canRequestDraw: typeof plugin?.canvas3d?.requestDraw === "function"
        },
        interaction: {
          localRotateMode: stageState.localRotateMode,
          rotationSensitivity: stageState.rotationSensitivity,
          interactiveFrameQueued: stageState.interactiveFrame !== null,
          interactiveFlushCount: stageState.interactiveFlushCount,
          lastInteractiveFlushAt: stageState.lastInteractiveReloadAt
        },
        keyboardShortcuts: {
          enabled: REPORT.keyboardShortcuts?.enabled === true,
          source: REPORT.keyboardShortcuts?.source || "",
          rotationStepDegrees: KEYBOARD_ROTATION_STEP,
          spacingStep: KEYBOARD_SPACING_STEP,
          thicknessStep: KEYBOARD_THICKNESS_STEP,
          bindings: (REPORT.keyboardShortcuts?.bindings || []).map((binding) => (binding.keys || []).join("/")),
          shortcutCount: stageState.keyboardShortcutCount,
          lastAction: stageState.lastKeyboardAction
        },
        urlState: urlStateSummary(),
        rotationMatrix: currentRotationMatrix(),
        screenRotationAxes: currentScreenRotationAxes(),
        lastRotationGesture: stageState.lastRotationGesture,
        visualMapping: visualMappingSummary(),
        selectedResidueMarker: selectedResidueMarkerSummary(),
        sceneSpacing: sceneSpacing(),
        visualEnvelope: visualEnvelopeSummary(),
        sceneFocusSphere: sceneFocusSphere(),
        centerAlignment: centerAlignmentSummary(),
        registries: {
          dataTransforms: registrySummary(dataState?.transforms),
          dataRegistryTransforms: registrySummary(dataState?.registry?.transforms),
          stateTransforms: registrySummary(plugin?.state?.transforms),
          pluginTransforms: registrySummary(plugin?.transforms)
        }
      };
      summary.flipbookInvariants = flipbookInvariantSummary(summary);
      stageState.capabilities = summary;
      return summary;
    }

    function publishReportApi() {
      const root = typeof globalThis !== "undefined" ? globalThis : window;
      root.FLIPBOOK_MOLSTAR_REPORT = {
        REPORT,
        stageState,
        loadAllSlices,
        setActiveSlice,
        setLayout,
        updateSpacing,
        updateTileColumns,
        setRotationAxis,
        rotateBy,
        resetRotation,
        resetVisualMapping,
        updateRadiusRange,
        updatePlayDelay,
        setSelectionMarkerEnabled,
        syncUrlState,
        updateMappingLegend,
        legendSummary,
        selectedResidueMarkerSummary,
        setLocalRotateMode,
        setSliceVisibility,
        toggleSliceVisibility,
        stepActiveSlice,
        currentScreenRotationAxes,
        applyScreenRotationDrag,
        updateSceneTransforms,
        queueGeometryUpdate,
        visualMappingSummary,
        flipbookInvariantSummary,
        molstarCapabilitySummary,
        stopPlayback
      };
    }

    function updateDiagnostics() {
      const diagnostics = document.getElementById("molstarDiagnostics");
      if (!diagnostics) {
        return;
      }
      diagnostics.textContent = JSON.stringify(molstarCapabilitySummary(), null, 2);
    }

    async function loadPdbIntoScene(plugin, entry, maskMode, alpha, tagSuffix) {
      const useStructureTransforms = canUseStructureTransforms(plugin);
      const useRepresentationTransforms = !useStructureTransforms && canUseRepresentationTransforms();
      const transformed = transformedPdbForSlice(entry.slice, entry.sceneIndex, maskMode, {
        applySceneTransform: !(useStructureTransforms || useRepresentationTransforms)
      });
      if (!transformed.atomCount) {
        return null;
      }
      const data = await plugin.builders.data.rawData({
        data: transformed.pdb,
        label: tagSuffix === "masked" ? `${entry.slice.filename} masked residues` : entry.slice.filename
      });
      const trajectory = await plugin.builders.structure.parseTrajectory(data, "pdb");
      const model = await plugin.builders.structure.createModel(trajectory);
      const structure = await plugin.builders.structure.createStructure(model);
      const transformCell = useStructureTransforms ? await insertSliceTransform(plugin, structure, entry, tagSuffix) : null;
      const displayStructure = transformCell || structure;
      const representation = await addRmsxRepresentation(plugin, displayStructure, entry.sourceIndex, { alpha, tagSuffix });
      const record = {
        ...entry,
        kind: "rmsx-worm",
        maskMode,
        alpha,
        data,
        trajectory,
        model,
        structure,
        transformCell,
        representationTransform: useRepresentationTransforms,
        displayStructure,
        representation
      };
      if (record.representationTransform) {
        applyRecordRepresentationTransform(record);
      }
      return record;
    }

    async function loadSelectedResidueMarkerIntoScene(plugin, entry) {
      if (!stageState.selectionMarkerEnabled) {
        return null;
      }
      const useStructureTransforms = canUseStructureTransforms(plugin);
      const useRepresentationTransforms = !useStructureTransforms && canUseRepresentationTransforms();
      const marker = selectedResidueMarkerPdbForSlice(entry.slice, entry.sceneIndex, {
        applySceneTransform: !(useStructureTransforms || useRepresentationTransforms)
      });
      if (!marker.atomCount) {
        return null;
      }
      const data = await plugin.builders.data.rawData({
        data: marker.pdb,
        label: `${entry.slice.filename} selected residue ${marker.residueKey}`
      });
      const trajectory = await plugin.builders.structure.parseTrajectory(data, "pdb");
      const model = await plugin.builders.structure.createModel(trajectory);
      const structure = await plugin.builders.structure.createStructure(model);
      const transformCell = useStructureTransforms ? await insertSliceTransform(plugin, structure, entry, "selected-residue-marker") : null;
      const displayStructure = transformCell || structure;
      const representation = await addSelectedResidueMarkerRepresentation(plugin, displayStructure, entry);
      const record = {
        ...entry,
        kind: "selected-residue-marker",
        maskMode: "selected-residue",
        alpha: markerTypeParams().alpha,
        residueKey: marker.residueKey,
        atomCount: marker.atomCount,
        matchMode: marker.matchMode,
        data,
        trajectory,
        model,
        structure,
        transformCell,
        representationTransform: useRepresentationTransforms,
        displayStructure,
        representation
      };
      if (record.representationTransform) {
        applyRecordRepresentationTransform(record);
      }
      return record;
    }

    async function loadSliceIntoScene(plugin, entry) {
      const records = [];
      if (!hasMaskedResidues()) {
        const fullRecord = await loadPdbIntoScene(plugin, entry, "all", 1, "all");
        if (fullRecord) {
          records.push(fullRecord);
        }
      } else {
        const unmaskedRecord = await loadPdbIntoScene(plugin, entry, "unmasked", 1, "unmasked");
        if (unmaskedRecord) {
          records.push(unmaskedRecord);
        }
        const maskedRecord = await loadPdbIntoScene(plugin, entry, "masked", REPORT.maskOpacity, "masked");
        if (maskedRecord) {
          records.push(maskedRecord);
        }
      }
      const markerRecord = await loadSelectedResidueMarkerIntoScene(plugin, entry);
      if (markerRecord) {
        records.push(markerRecord);
      }
      return records;
    }

    async function loadAllSlices(autoView = true) {
      const version = ++stageState.loadVersion;
      const entries = visibleSliceEntries();
      updateMetrics();
      if (stageState.layout === "flip") {
        setStatus(`Loading ${entries.length} RMSX slices in one Molstar scene for Flipbook playback...`);
      } else {
        setStatus(`Loading ${entries.length} RMSX slice${entries.length === 1 ? "" : "s"} in one Molstar scene...`);
      }
      const viewer = await createViewer();
      const plugin = viewer.plugin;
      const records = [];
      try {
        for (const entry of entries) {
          records.push(...await loadSliceIntoScene(plugin, entry));
        }
      } catch (error) {
        if (!stageState.forcePdbCoordinateTransforms && String(error.message || "").includes("FLIPBOOK_MOLSTAR_STRUCTURE_TRANSFORM_FAILED")) {
          console.warn("Molstar state transforms failed; retrying with browser-side PDB coordinate transforms.", error);
          stageState.forcePdbCoordinateTransforms = true;
          stageState.geometryMode = "pdb-coordinate-rewrite";
          updateDiagnostics();
          return loadAllSlices(autoView);
        }
        throw error;
      }
      if (version !== stageState.loadVersion) {
        return;
      }
      stageState.records = records.filter(Boolean);
      stageState.geometryMode = hasStructureTransformRecords()
        ? "molstar-state-transform"
        : hasRepresentationTransformRecords()
          ? "molstar-representation-transform"
          : "pdb-coordinate-rewrite";
      if (hasRepresentationTransformRecords()) {
        updateRepresentationTransformRecords(plugin, { fast: false });
      }
      applySliceVisibilityToRecords({ fast: false });
      if (autoView) {
        resetView();
      }
      updateMetrics();
      setLoadedSceneStatus();
      publishReportApi();
      updateDiagnostics();
    }

    function resetView() {
      const plugin = stageState.viewer?.plugin;
      const sphere = sceneFocusSphere();
      if (sphere && plugin?.managers?.camera?.focusSphere) {
        plugin.managers.camera.focusSphere(sphere, { durationMs: 0, extraRadius: Math.max(4, sphere.radius * 0.08) });
      } else if (plugin?.managers?.camera?.reset) {
        plugin.managers.camera.reset();
      } else if (plugin?.canvas3d?.requestCameraReset) {
        plugin.canvas3d.requestCameraReset();
      }
    }

    function maskStatusText() {
      if (!hasMaskedResidues()) {
        return "";
      }
      const percentOpaque = Math.round((REPORT.maskOpacity ?? 0.3) * 100);
      return `; ${REPORT.maskSummary.maskedResidues} masked residue${REPORT.maskSummary.maskedResidues === 1 ? "" : "s"} rendered at ${percentOpaque}% opacity`;
    }

    function rotationStatusText() {
      const rotation = stageState.rotation;
      if (Math.abs(rotation.x) < 0.5 && Math.abs(rotation.y) < 0.5 && Math.abs(rotation.z) < 0.5) {
        return "";
      }
      return `; screen-axis local rotation x/y/z ${rotation.x.toFixed(0)}/${rotation.y.toFixed(0)}/${rotation.z.toFixed(0)} deg`;
    }

    function geometryModeText() {
      if (stageState.geometryMode === "molstar-state-transform") {
        return "; slice geometry updated through Molstar state transforms";
      }
      if (stageState.geometryMode === "molstar-representation-transform") {
        return "; slice geometry updated through Molstar representation transforms";
      }
      return "; slice geometry baked into browser-side PDB copies";
    }

    function visibilityStatusText() {
      if (stageState.layout === "flip") {
        return "";
      }
      const visibleCount = visibleSliceIndexList().length;
      if (visibleCount === REPORT.slices.length) {
        return "";
      }
      return `; ${visibleCount}/${REPORT.slices.length} slices visible`;
    }

    function setLoadedSceneStatus() {
      if (stageState.layout === "flip") {
        setStatus(`${REPORT.slices.length} slices loaded in one Molstar scene; ${REPORT.slices[stageState.currentIndex].label} is displayed as the current Flipbook frame using RMSX-normalized B-factors with ${stageState.representationMode} uncertainty styling and the ${stageState.paletteName} palette${maskStatusText()}${rotationStatusText()}${geometryModeText()}.`);
      } else {
        setStatus(`${REPORT.slices.length} slice${REPORT.slices.length === 1 ? "" : "s"} loaded together in Molstar using RMSX-normalized B-factors with ${stageState.representationMode} uncertainty styling and the ${stageState.paletteName} palette${maskStatusText()}${rotationStatusText()}${visibilityStatusText()}${geometryModeText()}.`);
      }
    }

    function updateSceneLabels() {
      if (!sceneLabels) {
        return;
      }
      sceneLabels.replaceChildren();
      REPORT.slices.forEach((slice, index) => {
        const label = document.createElement("button");
        const swatch = document.createElement("span");
        const visible = isSliceVisible(index);
        swatch.className = "swatch";
        swatch.style.background = SLICE_COLORS[index % SLICE_COLORS.length];
        label.type = "button";
        label.className = "scene-label";
        label.classList.toggle("active", index === stageState.currentIndex);
        label.classList.toggle("hidden", !visible);
        label.setAttribute("aria-pressed", visible ? "true" : "false");
        label.setAttribute("data-testid", `molstar-slice-visible-${index + 1}`);
        label.addEventListener("click", () => {
          if (stageState.layout === "flip") {
            setActiveSlice(index);
            return;
          }
          toggleSliceVisibility(index, { focus: true });
        });
        label.appendChild(swatch);
        label.append(slice.label);
        sceneLabels.appendChild(label);
      });
    }

    function setActiveSlice(index, options = {}) {
      stageState.currentIndex = (index + REPORT.slices.length) % REPORT.slices.length;
      if (stageState.layout !== "flip" && options.ensureVisible !== false && !isSliceVisible(stageState.currentIndex)) {
        stageState.visibleSliceIndexes.add(stageState.currentIndex);
        applySliceVisibilityToRecords({ fast: true });
      }
      if (stageState.layout === "flip") {
        if (stageState.records.length) {
          applySliceVisibilityToRecords({ fast: true });
        } else if (options.reload !== false) {
          syncUrlState();
          loadAllSlices(false).catch((error) => setStatus(error.message, true));
          return;
        }
      }
      syncUrlState();
      updateMetrics();
      if (options.status !== false) {
        setLoadedSceneStatus();
      }
      updateDiagnostics();
    }

    function setSliceVisibility(index, visible, options = {}) {
      if (stageState.layout === "flip") {
        setActiveSlice(index);
        return;
      }
      const normalizedIndex = (index + REPORT.slices.length) % REPORT.slices.length;
      const nextVisible = new Set(stageState.visibleSliceIndexes);
      if (visible) {
        nextVisible.add(normalizedIndex);
      } else if (nextVisible.size > 1) {
        nextVisible.delete(normalizedIndex);
      } else {
        return;
      }
      stageState.visibleSliceIndexes = nextVisible;
      if (!stageState.visibleSliceIndexes.has(stageState.currentIndex)) {
        stageState.currentIndex = firstVisibleSliceIndex();
      }
      applySliceVisibilityToRecords({ fast: true });
      if (options.focus) {
        resetView();
      }
      syncUrlState();
      updateMetrics();
      setLoadedSceneStatus();
      updateDiagnostics();
    }

    function toggleSliceVisibility(index, options = {}) {
      setSliceVisibility(index, !stageState.visibleSliceIndexes.has(index), options);
    }

    function stepActiveSlice(delta) {
      const candidates = stageState.layout === "flip" ? allSliceIndexes() : visibleSliceIndexList();
      if (!candidates.length) {
        return;
      }
      const currentPosition = candidates.indexOf(stageState.currentIndex);
      const start = currentPosition === -1 ? 0 : currentPosition;
      const nextPosition = (start + delta + candidates.length) % candidates.length;
      setActiveSlice(candidates[nextPosition], { ensureVisible: false });
    }

    function stopPlayback() {
      if (stageState.timer) {
        window.clearInterval(stageState.timer);
        stageState.timer = null;
      }
      updatePlaybackControls();
      updateDiagnostics();
    }

    function startPlayback() {
      stopPlayback();
      stageState.timer = window.setInterval(() => {
        stepActiveSlice(1);
      }, stageState.playDelayMs);
      updatePlaybackControls();
      updateDiagnostics();
    }

    function togglePlayback() {
      if (stageState.timer) {
        stopPlayback();
        return;
      }
      startPlayback();
    }

    function updateLayoutControls() {
      tiledButton.classList.toggle("active", stageState.layout === "tiled");
      overlayButton.classList.toggle("active", stageState.layout === "overlay");
      flipButton.classList.toggle("active", stageState.layout === "flip");
    }

    function setLayout(layout) {
      if (stageState.layout === layout) {
        return;
      }
      stageState.layout = layout;
      if (layout !== "flip" && !isSliceVisible(stageState.currentIndex)) {
        stageState.currentIndex = firstVisibleSliceIndex();
      }
      updateLayoutControls();
      syncUrlState();
      loadAllSlices(true).catch((error) => setStatus(error.message, true));
    }

    function updateThickness(value) {
      const scale = clamp(Number(value), 0.25, 2.5);
      if (!Number.isFinite(scale)) {
        return;
      }
      stageState.thicknessScale = scale;
      thicknessRange.value = scale.toFixed(3);
      thicknessNumber.value = scale.toFixed(3);
      updateMetrics();
      syncUrlState();
      queueSceneReload(false);
    }

    function updateRadiusRange(bound, value) {
      const numeric = clamp(Number(value), 0.05, 8);
      if (!Number.isFinite(numeric)) {
        return;
      }
      if (bound === "min") {
        stageState.radiusMin = Math.min(numeric, stageState.radiusMax - 0.01);
      } else if (bound === "max") {
        stageState.radiusMax = Math.max(numeric, stageState.radiusMin + 0.01);
      } else {
        return;
      }
      radiusMinNumber.value = stageState.radiusMin.toFixed(3);
      radiusMaxNumber.value = stageState.radiusMax.toFixed(3);
      updateMetrics();
      syncUrlState();
      queueSceneReload(false);
    }

    function updateColorDomain(bound, value) {
      const numeric = clamp(Number(value), REPORT.domain.min, REPORT.domain.max);
      if (!Number.isFinite(numeric)) {
        return;
      }
      if (bound === "min") {
        stageState.colorDomainMin = Math.min(numeric, colorDomainMax() - 0.000001);
      } else if (bound === "max") {
        stageState.colorDomainMax = Math.max(numeric, colorDomainMin() + 0.000001);
      } else {
        return;
      }
      colorMinNumber.value = colorDomainMin().toFixed(3);
      colorMaxNumber.value = colorDomainMax().toFixed(3);
      updateMetrics();
      syncUrlState();
      queueSceneReload(false);
    }

    function handleSelectedResidueChange() {
      updateMetrics();
      syncUrlState();
      if (stageState.selectionMarkerEnabled) {
        queueSceneReload(false);
      } else {
        updateDiagnostics();
      }
    }

    function setSelectionMarkerEnabled(enabled) {
      stageState.selectionMarkerEnabled = Boolean(enabled);
      updateResidueMarkerControls();
      syncUrlState();
      updateDiagnostics();
      queueSceneReload(false);
    }

    function updatePalette(value) {
      const name = String(value || "").toLowerCase();
      if (!availablePalettes()[name]) {
        return;
      }
      stageState.paletteName = name;
      paletteSelect.value = stageState.paletteName;
      updateMetrics();
      syncUrlState();
      queueSceneReload(false);
    }

    function resetVisualMapping(options = {}) {
      stageState.thicknessScale = DEFAULT_THICKNESS_SCALE;
      stageState.radiusMin = DEFAULT_RADIUS_MIN;
      stageState.radiusMax = DEFAULT_RADIUS_MAX;
      stageState.colorDomainMin = DEFAULT_COLOR_MIN;
      stageState.colorDomainMax = DEFAULT_COLOR_MAX;
      thicknessRange.value = stageState.thicknessScale.toFixed(3);
      thicknessNumber.value = stageState.thicknessScale.toFixed(3);
      radiusMinNumber.value = stageState.radiusMin.toFixed(3);
      radiusMaxNumber.value = stageState.radiusMax.toFixed(3);
      colorMinNumber.value = colorDomainMin().toFixed(3);
      colorMaxNumber.value = colorDomainMax().toFixed(3);
      updateMetrics();
      syncUrlState();
      if (options.reload !== false) {
        queueSceneReload(false);
      }
    }

    function updateSpacing(value) {
      const spacing = clamp(Number(value), MIN_TILE_SPACING_FACTOR, MAX_TILE_SPACING_FACTOR);
      if (!Number.isFinite(spacing)) {
        return;
      }
      stageState.spacingFactor = spacing;
      spacingRange.value = spacing.toFixed(3);
      spacingNumber.value = spacing.toFixed(3);
      syncUrlState();
      if (stageState.layout === "tiled") {
        queueGeometryUpdate(false);
      }
    }

    function updateTileColumns(value) {
      const columns = clamp(Math.round(Number(value)), 1, Math.max(1, REPORT.slices.length));
      if (!Number.isFinite(columns)) {
        return;
      }
      stageState.tileColumns = columns;
      columnsNumber.value = String(columns);
      syncUrlState();
      if (stageState.layout === "tiled") {
        queueGeometryUpdate(true);
      }
    }

    function updatePlayDelay(value) {
      const delay = clamp(Math.round(Number(value)), MIN_PLAY_DELAY_MS, MAX_PLAY_DELAY_MS);
      if (!Number.isFinite(delay)) {
        return;
      }
      const wasPlaying = Boolean(stageState.timer);
      stageState.playDelayMs = delay;
      updatePlaybackControls();
      syncUrlState();
      if (wasPlaying) {
        startPlayback();
      } else {
        updateDiagnostics();
      }
    }

    function isKeyboardEditableTarget(target) {
      if (!target) {
        return false;
      }
      const tagName = String(target.tagName || "").toLowerCase();
      return target.isContentEditable || ["input", "select", "textarea"].includes(tagName);
    }

    function recordKeyboardAction(action) {
      stageState.keyboardShortcutCount += 1;
      stageState.lastKeyboardAction = action;
      updateDiagnostics();
    }

    function handleFlipbookKeyboard(event) {
      if (event.defaultPrevented || event.ctrlKey || event.metaKey || event.altKey || isKeyboardEditableTarget(event.target)) {
        return;
      }

      const key = event.key.length === 1 ? event.key.toLowerCase() : event.key;
      let action = "";
      switch (key) {
        case "u":
          rotateBy("x", KEYBOARD_ROTATION_STEP);
          action = "rotate-x-positive";
          break;
        case "i":
          rotateBy("x", -KEYBOARD_ROTATION_STEP);
          action = "rotate-x-negative";
          break;
        case "n":
          rotateBy("y", KEYBOARD_ROTATION_STEP);
          action = "rotate-y-positive";
          break;
        case "m":
          rotateBy("y", -KEYBOARD_ROTATION_STEP);
          action = "rotate-y-negative";
          break;
        case "j":
          rotateBy("z", KEYBOARD_ROTATION_STEP);
          action = "rotate-z-positive";
          break;
        case "k":
          rotateBy("z", -KEYBOARD_ROTATION_STEP);
          action = "rotate-z-negative";
          break;
        case "=":
        case "+":
          updateSpacing(stageState.spacingFactor + KEYBOARD_SPACING_STEP);
          action = "spacing-increase";
          break;
        case "-":
          updateSpacing(stageState.spacingFactor - KEYBOARD_SPACING_STEP);
          action = "spacing-decrease";
          break;
        case "[":
          updateThickness(stageState.thicknessScale + KEYBOARD_THICKNESS_STEP);
          action = "thickness-increase";
          break;
        case "]":
          updateThickness(stageState.thicknessScale - KEYBOARD_THICKNESS_STEP);
          action = "thickness-decrease";
          break;
        case ",":
          updateColorDomain("min", stageState.colorDomainMin + KEYBOARD_COLOR_STEP);
          action = "color-domain-low-increase";
          break;
        case ".":
          updateColorDomain("max", stageState.colorDomainMax - KEYBOARD_COLOR_STEP);
          action = "color-domain-high-decrease";
          break;
        case "ArrowLeft":
          stepActiveSlice(-1);
          action = "previous-slice";
          break;
        case "ArrowRight":
          stepActiveSlice(1);
          action = "next-slice";
          break;
        case "t":
          setLayout("tiled");
          action = "layout-tiled";
          break;
        case "o":
          setLayout("overlay");
          action = "layout-overlay";
          break;
        case "f":
          setLayout("flip");
          action = "layout-flip";
          break;
        default:
          return;
      }
      event.preventDefault();
      event.stopPropagation();
      recordKeyboardAction(action);
    }

    function setRotationAxis(axis, value, options = {}) {
      if (!["x", "y", "z"].includes(axis)) {
        return;
      }
      stageState.rotation[axis] = wrapAngle(Number(value));
      syncRotationMatrixFromEuler();
      updateRotationControls();
      syncUrlState();
      queueGeometryUpdate(options.autoView ?? false);
    }

    function rotateBy(axis, delta, options = {}) {
      if (!["x", "y", "z"].includes(axis)) {
        return;
      }
      setRotationAxis(axis, stageState.rotation[axis] + Number(delta), options);
    }

    function resetRotation(options = {}) {
      stageState.rotation = { ...DEFAULT_ROTATION };
      syncRotationMatrixFromEuler();
      stageState.lastRotationGesture = null;
      updateRotationControls();
      syncUrlState();
      queueGeometryUpdate(options.autoView ?? false);
    }

    function setLocalRotateMode(enabled) {
      stageState.localRotateMode = Boolean(enabled);
      updateRotationControls();
    }

    function queueSceneReload(autoView, delay = 120) {
      window.clearTimeout(stageState.reloadTimer);
      stageState.reloadTimer = window.setTimeout(() => {
        loadAllSlices(autoView).catch((error) => setStatus(error.message, true));
      }, delay);
    }

    async function updateSceneTransforms(autoView = false, options = {}) {
      const plugin = stageState.viewer?.plugin;
      if (!plugin || !hasLiveTransformRecords()) {
        await loadAllSlices(autoView);
        return;
      }
      if (hasStructureTransformRecords()) {
        const records = stageState.records.filter((record) => record?.transformCell && (!options.visibleOnly || isSliceVisible(record.sourceIndex)));
        for (const record of records) {
          await updateRecordTransform(plugin, record);
        }
      } else if (hasRepresentationTransformRecords()) {
        updateRepresentationTransformRecords(plugin, options);
      }
      if (autoView) {
        resetView();
      }
      if (options.fast) {
        return;
      }
      updateMetrics();
      setLoadedSceneStatus();
    }

    function queueGeometryUpdate(autoView, delay = 60) {
      if (!hasLiveTransformRecords()) {
        queueSceneReload(autoView, delay);
        return;
      }
      window.clearTimeout(stageState.reloadTimer);
      stageState.reloadTimer = window.setTimeout(() => {
        updateSceneTransforms(autoView).catch((error) => {
          console.warn("Molstar state transform update failed; reloading scene.", error);
          stageState.forcePdbCoordinateTransforms = true;
          loadAllSlices(autoView).catch((reloadError) => setStatus(reloadError.message, true));
        });
      }, delay);
    }

    function queueInteractiveGeometryUpdate(autoView) {
      if (!hasLiveTransformRecords()) {
        queueInteractiveSceneReload(autoView);
        return;
      }
      stageState.pendingInteractiveAutoView = stageState.pendingInteractiveAutoView || Boolean(autoView);
      if (stageState.interactiveFrame !== null) {
        return;
      }
      const scheduleFrame = window.requestAnimationFrame || ((callback) => window.setTimeout(callback, 16));
      stageState.interactiveFrame = scheduleFrame(() => {
        const pendingAutoView = stageState.pendingInteractiveAutoView;
        stageState.pendingInteractiveAutoView = false;
        stageState.interactiveFrame = null;
        stageState.lastInteractiveReloadAt = window.performance.now();
        stageState.interactiveFlushCount += 1;
        updateSceneTransforms(pendingAutoView, { fast: true, visibleOnly: true }).catch((error) => {
          console.warn("Molstar interactive transform update failed; reloading scene.", error);
          stageState.forcePdbCoordinateTransforms = true;
          loadAllSlices(pendingAutoView).catch((reloadError) => setStatus(reloadError.message, true));
        });
      });
    }

    function queueInteractiveSceneReload(autoView) {
      const now = window.performance.now();
      const delay = Math.max(0, 180 - (now - stageState.lastInteractiveReloadAt));
      window.clearTimeout(stageState.reloadTimer);
      stageState.reloadTimer = window.setTimeout(() => {
        stageState.lastInteractiveReloadAt = window.performance.now();
        loadAllSlices(autoView).catch((error) => setStatus(error.message, true));
      }, delay);
    }

    function handleRotationPointerDown(event) {
      if (!stageState.localRotateMode || event.button !== 0) {
        return;
      }
      event.preventDefault();
      event.stopPropagation();
      stageState.rotationDrag = {
        pointerId: event.pointerId,
        x: event.clientX,
        y: event.clientY,
        axes: currentScreenRotationAxes()
      };
      applyMolstarRenderStyle(stageState.viewer?.plugin, { interactive: true });
      molstarViewport.classList.add("local-rotate-dragging");
      if (molstarViewport.setPointerCapture) {
        try {
          molstarViewport.setPointerCapture(event.pointerId);
        } catch (error) {
          console.debug("Molstar local rotation pointer capture skipped.", error);
        }
      }
    }

    function handleRotationPointerMove(event) {
      const drag = stageState.rotationDrag;
      if (!stageState.localRotateMode || !drag || drag.pointerId !== event.pointerId) {
        return;
      }
      event.preventDefault();
      event.stopPropagation();
      const coalesced = typeof event.getCoalescedEvents === "function" ? event.getCoalescedEvents() : [];
      const samples = coalesced.length ? coalesced : [event];
      const latest = samples[samples.length - 1] || event;
      const dx = latest.clientX - drag.x;
      const dy = latest.clientY - drag.y;
      const axes = drag.axes || currentScreenRotationAxes();
      stageState.rotationDrag = { pointerId: event.pointerId, x: latest.clientX, y: latest.clientY, axes };
      if (dx === 0 && dy === 0) {
        return;
      }
      applyScreenRotationDrag(dx, dy, axes);
      updateRotationControls();
      queueInteractiveGeometryUpdate(false);
    }

    function handleRotationPointerUp(event) {
      const drag = stageState.rotationDrag;
      if (!drag || drag.pointerId !== event.pointerId) {
        return;
      }
      event.preventDefault();
      event.stopPropagation();
      stageState.rotationDrag = null;
      molstarViewport.classList.remove("local-rotate-dragging");
      applyMolstarRenderStyle(stageState.viewer?.plugin, { interactive: false });
      if (molstarViewport.releasePointerCapture) {
        try {
          molstarViewport.releasePointerCapture(event.pointerId);
        } catch (error) {
          console.debug("Molstar local rotation pointer release skipped.", error);
        }
      }
      syncUrlState();
      queueGeometryUpdate(false, 20);
    }

    function wireEvents() {
      playButton.addEventListener("click", togglePlayback);
      resetButton.addEventListener("click", resetView);
      tiledButton.addEventListener("click", () => setLayout("tiled"));
      overlayButton.addEventListener("click", () => setLayout("overlay"));
      flipButton.addEventListener("click", () => setLayout("flip"));
      previousButton.addEventListener("click", () => stepActiveSlice(-1));
      nextButton.addEventListener("click", () => stepActiveSlice(1));
      sliceSelect.addEventListener("change", (event) => setActiveSlice(Number(event.target.value)));
      sliceRange.addEventListener("input", (event) => setActiveSlice(Number(event.target.value)));
      playDelayRange.addEventListener("input", (event) => updatePlayDelay(event.target.value));
      playDelayRange.addEventListener("change", (event) => updatePlayDelay(event.target.value));
      playDelayNumber.addEventListener("input", (event) => updatePlayDelay(event.target.value));
      playDelayNumber.addEventListener("change", (event) => updatePlayDelay(event.target.value));
      thicknessRange.addEventListener("input", (event) => updateThickness(event.target.value));
      thicknessNumber.addEventListener("input", (event) => updateThickness(event.target.value));
      spacingRange.addEventListener("input", (event) => updateSpacing(event.target.value));
      spacingRange.addEventListener("change", (event) => updateSpacing(event.target.value));
      spacingNumber.addEventListener("input", (event) => updateSpacing(event.target.value));
      spacingNumber.addEventListener("change", (event) => updateSpacing(event.target.value));
      colorMinNumber.addEventListener("input", (event) => updateColorDomain("min", event.target.value));
      colorMinNumber.addEventListener("change", (event) => updateColorDomain("min", event.target.value));
      colorMaxNumber.addEventListener("input", (event) => updateColorDomain("max", event.target.value));
      colorMaxNumber.addEventListener("change", (event) => updateColorDomain("max", event.target.value));
      paletteSelect.addEventListener("change", (event) => updatePalette(event.target.value));
      radiusMinNumber.addEventListener("input", (event) => updateRadiusRange("min", event.target.value));
      radiusMinNumber.addEventListener("change", (event) => updateRadiusRange("min", event.target.value));
      radiusMaxNumber.addEventListener("input", (event) => updateRadiusRange("max", event.target.value));
      radiusMaxNumber.addEventListener("change", (event) => updateRadiusRange("max", event.target.value));
      resetScaleButton.addEventListener("click", () => resetVisualMapping());
      columnsNumber.addEventListener("input", (event) => updateTileColumns(event.target.value));
      columnsNumber.addEventListener("change", (event) => updateTileColumns(event.target.value));
      rotationXRange.addEventListener("input", (event) => setRotationAxis("x", event.target.value));
      rotationXNumber.addEventListener("input", (event) => setRotationAxis("x", event.target.value));
      rotationYRange.addEventListener("input", (event) => setRotationAxis("y", event.target.value));
      rotationYNumber.addEventListener("input", (event) => setRotationAxis("y", event.target.value));
      rotationZRange.addEventListener("input", (event) => setRotationAxis("z", event.target.value));
      rotationZNumber.addEventListener("input", (event) => setRotationAxis("z", event.target.value));
      rotateXButton.addEventListener("click", () => rotateBy("x", 15));
      rotateYButton.addEventListener("click", () => rotateBy("y", 15));
      rotateZButton.addEventListener("click", () => rotateBy("z", 15));
      localRotateButton.addEventListener("click", () => setLocalRotateMode(!stageState.localRotateMode));
      resetRotationButton.addEventListener("click", () => resetRotation());
      molstarViewport.addEventListener("pointerdown", handleRotationPointerDown, true);
      molstarViewport.addEventListener("pointermove", handleRotationPointerMove, true);
      molstarViewport.addEventListener("pointerup", handleRotationPointerUp, true);
      molstarViewport.addEventListener("pointercancel", handleRotationPointerUp, true);
      markerToggleButton.addEventListener("click", () => setSelectionMarkerEnabled(!stageState.selectionMarkerEnabled));
      residueSelect.addEventListener("change", handleSelectedResidueChange);
      document.addEventListener("keydown", handleFlipbookKeyboard);
    }

    async function init() {
      populateControls();
      wireEvents();
      stageState.currentIndex = initialSliceIndex();
      if (stageState.layout !== "flip" && !isSliceVisible(stageState.currentIndex)) {
        stageState.currentIndex = firstVisibleSliceIndex();
      }
      updateLayoutControls();
      syncUrlState();
      updateMetrics();
      if (!window.molstar?.Viewer) {
        setStatus("Molstar did not load. Check network access to the pinned Molstar browser bundle.", true);
        return;
      }
      if (stageState.preferStructureTransforms) {
        await loadMolstarStateTransforms();
      }
      await loadAllSlices(true);
      if (stageState.autoplay) {
        togglePlayback();
      }
    }

    document.addEventListener("DOMContentLoaded", () => {
      init().catch((error) => setStatus(error.message, true));
    });
  </script>
</body>
</html>
"""
    return (
        template.replace("__TITLE__", escaped_title)
        .replace("__MOLSTAR_JS_URL__", MOLSTAR_JS_URL)
        .replace("__MOLSTAR_CSS_URL__", MOLSTAR_CSS_URL)
        .replace("__PALETTE_GRADIENT__", palette_gradient)
        .replace("__PAYLOAD__", payload_json)
        .replace("__MIN_PLAY_DELAY_MS__", str(payload["playback"]["minDelayMs"]))
        .replace("__MAX_PLAY_DELAY_MS__", str(payload["playback"]["maxDelayMs"]))
        .replace("__DEFAULT_PLAY_DELAY_MS__", str(payload["playback"]["defaultDelayMs"]))
        .replace("__PLAY_DELAY_STEP_MS__", str(payload["playback"]["delayStepMs"]))
        .replace("__MIN_TILE_SPACING_FACTOR__", f"{MIN_TILE_SPACING_FACTOR:.3f}")
        .replace("__MAX_TILE_SPACING_FACTOR__", f"{MAX_TILE_SPACING_FACTOR:.3f}")
        .replace("__DEFAULT_TILE_SPACING_FACTOR__", f"{DEFAULT_TILE_SPACING_FACTOR:.3f}")
    )


def main():
    args = parse_args()
    slices = read_slices(args.pdb_dir)
    rmsx_rows, slice_columns = read_rmsx_table(args.rmsx_table)
    summaries, domain = summarize_slices(rmsx_rows, slice_columns)
    residues = build_residue_payload(rmsx_rows, slice_columns)
    mask_summary = read_mask_summary(args.mask_table)
    payload = build_viewer_payload(args.title, slices, summaries, domain, mask_summary, residues, args.palette)
    if args.output:
        output = Path(args.output)
        output.write_text(
            html_report(payload),
            encoding="utf-8",
        )
    if args.manifest_output:
        manifest = Path(args.manifest_output)
        manifest.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
