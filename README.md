# Schwarzschild Black Hole Ray Tracer

A Python-based, illustrative Schwarzschild black-hole thin-disk ray tracer. The simulator traces null photon trajectories, locates the first visible disk crossing, applies relativistic frequency-shift factors, and renders a supersampled image as well as numerical arrays.

**Status:** Version 1.0 release candidate. A successful fixed-grid regression does not prove image convergence or astrophysical realism.

## Example output

After rendering, open `figures/disk_standalone_300px_4x4_h0.008.png` and inspect the image. Add a verified representative render here before publishing the repository.

## Requirements

- Python 3.10 or later
- NumPy
- Matplotlib
- No SciPy dependency

Install dependencies:

```bash
python -m pip install -r requirements.txt
```

## Run

Small demo:

```bash
python render_v010s.py --size 48 --samples 2
```

Higher-resolution render:

```bash
python render_v010s.py --size 300 --samples 4
```

Produces `.png`, `.npz`, and `.json` output in `figures/`. The `.npz` file contains `linear_rgb`, `display_rgb`, and `hit_fraction` arrays.

## Tests

```bash
python release_selftest.py
```

For the local developer release gate, use `verify_release_v010v.py` on the packaged release ZIP. That script is a **prepublication tool** and is not required to run the standalone renderer.

## Physical model

- Nonrotating, uncharged Schwarzschild spacetime, with geometrized units \(G=c=M=1\).
- A static finite-distance observer at \(r=100M\), with a pinhole camera.
- Thin circular disk with radial extent \(6M\) to \(20M\), inclined at 70°.
- Photon orbit equation \(d^2u/d\phi^2=3Mu^2-u\), with \(u=1/r\), advanced with fourth-order Runge–Kutta.
- Circular-geodesic disk matter, relativistic Doppler and gravitational frequency shifts.
- Illustrative disk emissivity and RGB colors; bolometric \(g^4\) multiplier and nonlinear display tone mapping.

## What this simulation does not claim

The rendered RGB values are **not** a calibrated thermal spectrum or telescope observation. The model omits spin, disk thickness, absorption, scattering, and full relativistic radiative transfer. The central dark region is not automatically the precise vacuum shadow boundary. Very narrow higher-order features can be undersampled even when integration-step and fixed-grid regression tests pass.

## Reproducibility

The 300×300, 4×4 fixed-grid reference produced a camera-plane model brightness integral of `51.1497137695` in the validated environment. The 48×48, 2×2 reference produced `52.0340669609`. The integral is a model-specific image statistic, **not a physical luminosity**. The fixed-grid tests validate preservation of the implementation, not the continuum image.

See `VALIDATION.json` and `CLEAN_INSTALL_VALIDATION.json` for the recorded results. The release ZIP's original hashes were verified during the clean-install test.

## License

A license has **not** been selected. Choose an appropriate license before public distribution; do not imply that a license exists until you add one.
