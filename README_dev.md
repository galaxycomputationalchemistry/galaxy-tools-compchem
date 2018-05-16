
# Notes on dev

## Planemo - found using in a .venv to be best for me
Make sure to be out out of conda environments - `source deactivate`
```
$ virtualenv .venv; . .venv/bin/activate
$ pip install "pip>=7" # Upgrade pip if needed.
$ pip install planemo
```
## Using conda with planemo for tooldeps
```
$ planemo conda_init # NOT NEEDED IF conda is already installed on machine, if unsure you will see error message about it
planemo conda_install --conda_channels=omnia --conda_prefix /scratch/cbarnett/conda_for_planemo .
```
 
## Lint/test

`planemo lint *.xml`
`planemo test --conda_dependency_resolution --conda_prefix /scratch/cbarnett/conda_for_planemo/ galaxy_packmol.xml`

## Serve
`planemo serve --port 9091 --galaxy_branch release_18.01 --conda_dependency_resolution --conda_prefix /scratch/cbarnett/conda_for_planemo/`
