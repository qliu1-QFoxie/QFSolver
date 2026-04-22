# Internal Layers

The package is four computational layers.

- `background/`: solves the homogeneous inflationary background and derived
  background observables.
- `perturbations/`: solves scalar and tensor perturbation modes on a fixed
  background.
- `spectra/`: constructs dense endpoint primordial spectra, diagnostics, and
  explicit physical-`k` calibration.
- `workflow/`: provides the public request/result orchestration layer.

The CLASS export layer consumes an already solved calibrated spectrum and writes
a table suitable for CLASS `external_Pk`; it does not re-solve modes or alter
the physical calibration.