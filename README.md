# TOVBrachistochrone

Python scripts for computing and plotting brachistochrone trajectories inside TOV stellar interiors.

## Contents

- `tov_brachistochrone.py`: core TOV integration, EOS handling, effective indices, brachistochrone curves and travel times.
- `genera_grafici_tov.py`: APS-style figure generator for the TOV brachistochrone figures, with PDF/SVG/PNG/CSV outputs.
- `brach_anisotropic.py`: comparison of isotropic and Bowers-Liang anisotropic models at fixed compactness and opening angle.

## Requirements

Install the Python dependencies with:

```bash
pip install -r requirements.txt
```

## Example usage

Generate the default polytropic figures:

```bash
python genera_grafici_tov.py --outdir figures
```

Generate only the TOV profile figure:

```bash
python genera_grafici_tov.py --figure profiles --outdir figures
```
