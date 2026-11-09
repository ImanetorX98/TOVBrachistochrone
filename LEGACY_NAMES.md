# Legacy names (v1 → v2)

Version 1 (tag `v1.0-submission`, Zenodo concept DOI 10.5281/zenodo.21797700) contained four files at the root of
the repository. In v2 they are kept unchanged in function, moved to `legacy_v1/` and given English names. The tag
`v1.0-submission` still contains the original files.

| v1 path | v2 path | Change |
|---|---|---|
| `tov_brachistochrone.py` | `legacy_v1/tov_brachistochrone.py` | moved only |
| `genera_grafici_tov.py` | `legacy_v1/make_figures.py` | renamed; four comment lines translated from Italian; provenance string updated |
| `brach_anisotropic.py` | `legacy_v1/anisotropic_comparison.py` | renamed; one three-line comment translated from Italian; module docstring updated |
| `requirements.txt` | `legacy_v1/requirements.txt` | moved only |

The executable code is unchanged: the abstract syntax trees of the old and new files coincide once string constants
are blanked, and the only changed strings are the provenance line of `make_figures.py` and the module docstring of
`anisotropic_comparison.py`. `make_figures.py` imports `tov_brachistochrone.py` from its own folder, which Python puts on the module path when
the script is run by its path (for example `python3 legacy_v1/make_figures.py` from the repository root).
