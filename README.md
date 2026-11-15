# TOVBrachistochrone

Coordinate-time and proper-time brachistochrones inside Tolman–Oppenheimer–Volkoff stellar interiors.

A test particle released from rest at the stellar surface moves on a frictionless rail through the interior. Minimising
the coordinate time t or the proper time τ reduces to two Fermat problems on the static slice, with different optical
indices. Version 2 of this repository is the reproducibility package of a mathematical study of these problems. For
the uniform-density (Schwarzschild interior) star the curves, the aperture and the travel times are given in closed
form: abelian integrals on a curve of genus two (for regular turning radii with a square-free sextic), and complete
elliptic integrals for the expansion near the centre. The suite also checks the limiting model at the Buchdahl bound,
a numerical candidate for a proper-time coexistence separation, and the anisotropic Florides interior beyond the
Buchdahl bound.

> **Status.** Release 2.0.0 accompanies the paper *Brachistochrones inside static stars: existence, global minimality
> and exact results for two uniform-density interiors*. Version 1 (numerical code of the 2026 submission) is frozen at
> the tag `v1.0-submission`; its scripts are kept in `legacy_v1/` (see `LEGACY_NAMES.md`).

## Citation

Archived on Zenodo, concept DOI [10.5281/zenodo.21797699](https://doi.org/10.5281/zenodo.21797699) (all versions);
each release has its own version DOI. See [`CITATION.cff`](CITATION.cff).

## Where to start

- [`RESULTS_MAP.md`](RESULTS_MAP.md): the results implemented so far, the script that checks each of them, the kind
  of evidence (exact, numerical check, diagnostic, negative control) and the analytic steps that the suite does not
  execute.
- `python3 tools/run_verification.py --profile quick`: reruns the fast checks in a few seconds.

## Layout

| Path | Content |
|---|---|
| `verification/symbolic/` | Wolfram Language scripts: exact identities and decisions, and closed forms checked against independent quadratures or ODEs |
| `verification/diagnostics/` | sampled statements that the paper labels as diagnostics, not proofs |
| `verification/common.wl` | shared result recording (versioned JSON, written atomically) |
| `verification/third_kind_lemma.wl` | the complete third-kind lemma on a genus-two sextic, shared by the aperture and the travel times |
| `verification/manifests/` | `jobs.json` (reviewed contract: expected check identifiers, evidence kinds and the claim each supports), `profiles.json` (the jobs of each profile) |
| `tools/run_verification.py` | the single entry point; validates every JSON result against the manifest |
| `tools/test_run_verification.py` | negative tests of the validator |
| `legacy_v1/` | the numerical code of version 1 (general equations of state, figures) |
| `results/` | run outputs (`results/runs/` is not tracked); `results/release/v2.0.0/` holds the report of the release run |

## Running the verification

Requirements: Wolfram Language with `wolframscript` on the `PATH` (tested with 13.3.1) and Python 3.10 or later; the
v2 runner and its tests use the Python standard library only (the legacy v1 code has its own requirements, below). The
release profile needs a Git clone with a clean working tree: from the Zenodo archive, which has no `.git`, run the
`full` profile, or clone the tag from GitHub to run `release`. The release report records the versions of the
environment; no installable environment is shipped.

```bash
python3 tools/run_verification.py --profile quick      # exact identities, centre slopes and threshold (seconds)
python3 tools/run_verification.py --profile full       # every job, including the diagnostics (a few minutes)
python3 tools/run_verification.py --profile release    # full, on a clean tree, with SHA-256 digests recorded
python3 tools/test_run_verification.py                 # tests of the result validator
```

Options: `--output-dir DIR`, `--jobs N` (parallel jobs, default 2), `--progress auto|always|never`,
`--heartbeat SECONDS`. Every run reserves a new or empty output directory atomically and stamps a unique run
identifier into every result; use a directory that only you control. Progress goes to stderr. The run writes one JSON
result, one stdout log and one stderr log per job under `jobs/`, then an English `report.txt`, and last `summary.json`:
a `passed` summary attests that every selected job ran to completion with a valid, complete result and that the
manifests used are the ones whose SHA-256 is recorded (a change during the run fails it). The exit status is 0 only if every job passed; 1 for a failed check or an invalid
result; 2 for a usage or manifest error; 3 if the profile could not be completed (missing engine, timeout, unclean
or unknown Git state for `release`).

The release profile prepares evidence; it does not tag, publish or deposit anything. The report of the release run of
version 2.0.0, made on a fresh clone of the release-candidate commit, is archived in `results/release/v2.0.0/`
(`summary.json` with the commit, the engine version and the SHA-256 of every verification file, and `report.txt`).

## Conventions

Geometrized units G = c = 1 with R = 1 in the scripts; physical times are obtained by multiplying by R/c.
`mu = 2GM/(Rc^2)`, `y = sqrt(1 - mu r^2/R^2)`, `a = sqrt(1 - mu)`; the isotropic uniform-density interior requires
`mu < 8/9` (Buchdahl), the limiting model `mu = 8/9` is used only for limits, and the Florides interior allows `mu < 1`.
Complete elliptic integrals use the parameter convention m = k^2, and Π(n|m) has the denominator 1 - n sin^2θ
(DLMF §19.2).

## Version 1 code

```bash
python3 -m pip install -r legacy_v1/requirements.txt
python3 legacy_v1/make_figures.py --outdir figures
```

The tabulated equations of state (CompOSE tables) used by some v1 figures are not part of this repository; those
modes need the tables as explicit inputs.

## Licence

MIT, see [`LICENSE`](LICENSE).
