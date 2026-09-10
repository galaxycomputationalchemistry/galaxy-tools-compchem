#!/usr/bin/env bash
# Use the exact upstream registry tested by this workflow, without patching it.
set -euo pipefail

galaxy_fork=$1
galaxy_sha=$2
registry=$3
[[ "$galaxy_fork" =~ ^[A-Za-z0-9_.-]+$ ]]
[[ "$galaxy_sha" =~ ^[0-9a-f]{40}$ ]]

curl --fail --show-error --silent --location --retry 3 \
    "https://raw.githubusercontent.com/$galaxy_fork/galaxy/$galaxy_sha/lib/galaxy/config/sample/datatypes_conf.xml.sample" \
    --output "$registry"
python - "$registry" <<'PY'
import sys
import xml.etree.ElementTree as ET

root = ET.parse(sys.argv[1]).getroot()
datatype = root.find("./registration/datatype[@extension='rmsx.json']")
if datatype is None or datatype.find("./visualization[@plugin='rmsxflipbook']") is None:
    raise SystemExit("The selected upstream Galaxy registry lacks RMSX Flipbook registration")
PY
sha256sum "$registry"
printf 'Galaxy commit: %s\n' "$galaxy_sha"
