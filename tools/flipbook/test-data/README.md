# Flipbook Galaxy Test Data

This directory contains the PDB/XTC fixtures used by the Galaxy wrapper tests.

- `1UBQ.pdb`: ubiquitin structure fixture used with the RMSX example path.
- `mon_sys.xtc`: compressed trajectory fixture used with `1UBQ.pdb`; current
  size is 1,002,408 bytes.
- `protease_backbone.pdb` and `protease_compact.xtc`: the real two-chain
  protease assembly, retaining 790 backbone atoms and 99 C-alpha residues per
  chain, with overlapping residue IDs 1–99 in chains A and B. The trajectory
  contains 180 uniformly sampled source frames, giving 20 frames per slice in
  nine-slice acceptance tests. No chain is duplicated or translated.

## Source And Provenance

The fixture represents the ubiquitin case study described in the RMSX/Flipbook
Scientific Reports manuscript. The paper identifies the system as an SMD
trajectory of ubiquitin using PDB `1UBQ`, sourced from the NAMD case-study
materials for ubiquitin. It reports that the simulation used NAMD 2.14 with PME
electrostatics and periodic boundary conditions, fixed Lys48, steered Met1 at
0.05 A/ps with a 208.4 pN/A spring constant, recorded frames every 10 ps, and
stripped waters for analysis.

The paper citation for the case-study source is:

> Cruz-Chu, E. & Gumbart, J. C. Case study: Ubiquitin.
> https://www.ks.uiuc.edu/Training/CaseStudies/, 2016.

The exact source archive is the Ubiquitin "Required case study files" archive
from the Theoretical and Computational Biophysics Group case-studies page:

- Source page: `https://www.ks.uiuc.edu/Training/CaseStudies/`
- Archive URL: `https://www.ks.uiuc.edu/Training/CaseStudies/files/ubq-files.tgz`
- Archive name: `ubq-files.tgz`
- Archive size from server and local download: 42,158,765 bytes
- Server last-modified header: Thu, 24 Jan 2013 23:43:31 GMT
- Server ETag header: `"2834aad-4d41160d9678b"`
- Archive SHA256: `0c52e4727fc824c52fac42ad308a698132757835b5741996d0bac8543b004cc1`
- Archive SHA1: `86627940a57ba7bb547b7f9bd038fde6eae4350c`

The archive contains `1UBQ.pdb` and `mon_sys.dcd` at top level. Checksums for
the exact extracted source files are:

- `1UBQ.pdb`: 94,837 bytes; SHA256
  `3a9efa1922fb8c0b4967bcc29ce574780eb885f3c3b8c5b977f88fa1b06e7e25`
- `mon_sys.dcd`: 4,693,508 bytes; SHA256
  `040099c1725214333e899f742b93c507304c946ee0fc7e95808058f719ae4154`

The checked-in `1UBQ.pdb` is byte-identical to the archive copy. The checked-in
`mon_sys.xtc` is a reduced-size conversion derived from the archive's
`mon_sys.dcd`; its SHA256 is
`367c424cd9ff7506c671f5c8ee8f25e00c27d4a2f9aee91707dc4194bdc0676b`.

Redistribution note: the case-studies page links to the TCBG copyright
statement, which says the materials are copyrighted and may be reproduced and
distributed for educational use with credit. This appears compatible with a
small educational test fixture, but it is not a standard OSI/open-data license.
Before a community Tool Shed submission, confirm with the repository
maintainers whether this attribution-based educational-use statement is
acceptable for bundled test data, or replace the trajectory fixture with one
carrying a clearer open-data license.

## Community Repository Readiness

The current XTC fixture is useful for cofest reproducibility and is below the
community repository's 1 MB file-size check. It was generated from the original
full `mon_sys.dcd` by preserving all 316 frames and all atoms, then writing XTC
with precision 2. The command shape is:

```bash
docker run --rm \
  -v "$PWD:/work" \
  -w /work \
  ghcr.io/antuneslab/flipbook-galaxy:0.2.3-galaxy0 \
  python scripts/create_reduced_rmsx_fixture.py \
    --topology tools/flipbook/test-data/1UBQ.pdb \
    --trajectory /path/to/original/mon_sys.dcd \
    --output tools/flipbook/test-data/mon_sys.xtc \
    --frames 316 \
    --xtc-precision 2
```

The fixture still exercises:

- PDB/XTC loading through MDAnalysis/RMSX.
- Three trajectory slices.
- Static heatmap and triple-plot generation.
- PDB slice collection output.
- Molstar manifest generation.
- Two-chain discovery, per-chain outputs, combined PDB slices, and distinct
  `A:residue`/`B:residue` viewer keys.

## Real protease multi-chain fixture

The inputs are from the verified RMSX 0.1.5 source distribution (SHA256
`48a90b8b412d9cff95e108e7ff704ec16350a17025f7baacee51a81fb81cc3b4`),
under `rmsx/test_files/`. This is the same protease system used in the existing
viewer regression. The upstream RMSX distribution carries the MIT license;
its bundled protease data have no separate license notice. Retain upstream
attribution and request reviewer confirmation of data redistribution terms.
The educational-use notice above applies to the ubiquitin single-chain example.

Source and compact-fixture SHA256 checksums:

- `source protease_backbone.pdb` (64116 bytes): `45f98d39a5507cdf8b86a19b05147d053c054cb88f0f04b054c57aa322feddfe`
- `source short_protease_backbone.dcd` (47809916 bytes): `bd6207f770e8362725f2bfcf41039af5a7864871c4f6e6bd48b8eb1f2d7a4bc3`
- `protease_backbone.pdb` (64116 bytes): `742f100899a91e3b7e2bef946a796868639506c30cffc4f2de9e1243a2e8c100`
- `protease_compact.xtc` (559596 bytes): `19acb5966a6eeabd7ab88d856d75a230329f31b9aa2db1a3eea22438ddc863d0`

Regenerate from the extracted RMSX distribution in the pinned runtime:

```bash
python tools/flipbook/test-data/create_multichain_rmsx_fixture.py \
  --topology /path/to/rmsx/test_files/protease_backbone.pdb \
  --trajectory /path/to/rmsx/test_files/short_protease_backbone.dcd \
  --output-topology tools/flipbook/test-data/protease_backbone.pdb \
  --output-trajectory tools/flipbook/test-data/protease_compact.xtc --frames 180
```

The generator retains the source chain identities, coordinates, dimensions, and
chronological timestamps; XTC uses precision 3. Both output files are below 1 MB.
