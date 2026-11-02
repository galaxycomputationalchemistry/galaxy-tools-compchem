# Changelog

## 0.1.0+galaxy1

- First public-server candidate includes multi-chain analysis and Analysis metrics.
- Use the real two-chain protease fixture with 180 source frames and retained timestamps.
- Record actual analyzed frame/time bounds instead of interpolated slice bounds.
- Cover selected-chain frame windows and reject non-chronological RMSD metrics.

- Embed optional Analysis metrics without changing the v1 viewer schema,
  downsampling RMSD to at most 2,048 chronological extrema points while
  retaining every RMSF residue.
- Record per-slice time bounds and logical chain atom ranges for aligned
  chain-only and full-assembly Molstar lanes.
- Add focused generator tests for metric ordering, extrema preservation, time
  annotations, combined PDB serials, and multi-character logical chain IDs.

## 0.1.0+galaxy0

- Pin the RMSX runtime scaffold to upstream RMSX `v0.2.3`.
- Document that the `v0.2.3` source still installs package metadata as
  `rmsx==0.1.0`.
- Emit the Flipbook manifest as typed Galaxy `rmsx.json` for the native Molstar
  visualization proposed in `galaxyproject/galaxy#23009`.
- Make the complete 316-frame `1UBQ.pdb` and `mon_sys.xtc` example the default
  input path and cover it with a Galaxy tool test.
- Add the Flipbook logo, static plots, preflight validation, explicit test-data
  provenance, and Tool Shed publication notes.
- Link the wrapper to the `rmsx_and_flipbook` bio.tools registry entry.
- Analyze all valid chains by default, retain an explicit selected-chain mode,
  and return per-chain tables and plots with one combined Molstar manifest.
- Add a compact two-chain regression fixture that verifies duplicate residue
  numbers remain distinct through the RMSX, PDB, and viewer-manifest paths.
- Normalize combined PDB slices to one structure boundary so older RMSX
  runtimes cannot expose separate chain records as separate Molstar entries.
- Use the pinned public container as the sole runtime dependency until an RMSX
  Bioconda recipe is available.
