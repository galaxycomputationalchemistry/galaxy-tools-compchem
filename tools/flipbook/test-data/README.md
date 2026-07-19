# Flipbook Galaxy Test Data

This directory contains the PDB/XTC fixture used by the Galaxy wrapper tests.

- `1UBQ.pdb`: ubiquitin structure fixture used with the RMSX example path.
- `mon_sys.xtc`: compressed trajectory fixture used with `1UBQ.pdb`; current
  size is 1,002,408 bytes.

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
