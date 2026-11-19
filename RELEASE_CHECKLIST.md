# Final release checklist

- [ ] Run `python verify_release_v010v.py` on `releases/SchwarzschildRayTracer_v1.0rc1.zip` and verify `PASS_CLEAN_INSTALL`.
- [ ] Run `python prepare_github_v010v.py` to generate a separate GitHub folder.
- [ ] Open the final PNG output and visually inspect it for cropping, colour and disk/ring artefacts.
- [ ] Inspect the GitHub folder: `bh_engine_v010s.py`, `render_v010s.py`, `release_selftest.py`, `README.md`, `requirements.txt`, validation JSON files.
- [ ] Run `python release_selftest.py` from a new clean environment with only declared dependencies installed.
- [ ] Review scientific assumptions and limitations in the README; do not claim complete convergence or a calibrated astrophysical prediction.
- [ ] Add a representative image with a clear caption and appropriate permissions; update README with the image path.
- [ ] Choose a license deliberately (or document that no license is granted); do not invent a license attribution.
- [ ] Confirm no secrets, personal documents, giant simulation arrays, virtual environments or historical development archives are committed.
- [ ] Initialize the local Git repository, review changes, then publish to GitHub manually.
- [ ] Tag the first published version only when the release package and documentation are reviewed.

Note: The GitHub folder is not a new physics implementation. Its numerical baseline is the validated Version 0.10S standalone engine.
