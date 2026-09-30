# FIACH for MATLAB/SPM — R-compatible package

This repository contains a MATLAB/SPM implementation of the FIACH fMRI
noise-control pipeline. The default `GMMMethod='R'` follows the upstream R
implementation, including its DCT high-pass filter, robust tSNR calculation,
two-component GMM, spatial noise-mask expansion, PCA regressors, outlier
detection, short-gap interpolation, and median replacement.

The upstream reference is [tierneytim/FIACH](https://github.com/tierneytim/FIACH).

## Requirements

- MATLAB R2022b or newer (R2022b is the validated version).
- SPM12 on the MATLAB path or supplied to `setup_fiach`.
- Statistics and Machine Learning Toolbox only when using the optional
  `GMMMethod='ModifiedICH'`; it is not required for the R-compatible default.
- Git LFS for cloning or pushing the packaged NIfTI test data.

SPM12 and MATLAB are not bundled.

## Setup and GUI

From the repository root:

```matlab
setup_fiach('SPMRoot', 'C:\path\to\spm12')
fiach_gui
```

For an R-compatible command-line run:

```matlab
t = fiach_bold_contrast(3, 26);
result = fiach(functionalFiles, t, 1.25, ...
    'RP', motionFile, ...
    'MaxGap', 1, ...
    'Freq', 128, ...
    'NMads', 1.96, ...
    'GMMMethod', 'R');
```

`fiach_bold_contrast` also offers `Typical`, `Conservative` (the default and
upstream R settings), and `AlternativeTypical` physiology profiles. Their
assumptions, literature basis, and the limits of the packaged test data are
documented in [BOLD contrast parameters](docs/BOLD_CONTRAST_PARAMETERS.md).

By default, existing derivatives are not replaced. Add `'Overwrite', true`
only when replacement is intended.

## Full overt_0209 validation

The private test dataset contains 206 realigned input volumes, motion
parameters, and 206 upstream R-reference filtered volumes. Run:

```matlab
setup_fiach('SPMRoot', 'C:\path\to\spm12')
report = run_overt_0209_validation;
```

The test runs the packaged R-compatible MATLAB pipeline in the ignored `work/`
directory, compares every output volume with the R reference, writes the full
CSV and figures, and requires:

- 206 matched frames;
- mean in-mask correlation greater than 0.9999; and
- pooled in-mask mean absolute error below 2 intensity units.

During development the mean in-mask correlation was approximately 0.99999984.
The R reference writes zero outside its processing mask, while MATLAB retains
the original realigned values there, so whole-grid differences should not be
used as the primary compatibility measure. The reference NIfTI affine also
differs from the matching MATLAB array grid; comparisons are therefore made by
voxel-array index and report the header discrepancy.

## Visualisation tools

- `fiach_compare_outputs`: synchronized side-by-side movie, linked crosshair,
  voxel timecourses, and volume-level error metrics.
- `fiach_compare_gmm_options`: R-compatible versus ModifiedICH diagnostic and
  output comparison.
- `fiach_compare_pca_regressors`: sign/order-aligned PCA comparison.
- `fiach_processing_steps_viewer`: original, FIACH-filtered/interpolated, and
  motion-plus-FIACH nuisance-GLM stages with three linked voxel traces.

## Repository layout

```text
src/                       Core FIACH implementation and GUI
tools/                     Comparison and interactive viewer functions
tests/                     Automated sample and overt_0209 validation
examples/                  Reproducible batch scripts
legacy_reference/          Older ModifiedICH MATLAB files for provenance
test_data/overt_0209/      Inputs, R references, and development diagnostics
work/                      Generated outputs; ignored by Git
```

## GitHub and large files

Install Git LFS before the first commit:

```text
git lfs install
git add .gitattributes
git add .
git commit -m "Add R-compatible MATLAB FIACH implementation and validation data"
```

All NIfTI and MATLAB data files are already assigned to Git LFS by
`.gitattributes`. Keep the repository private unless the data-governance and
licensing checks in `test_data/README.md` and `LICENSING.md` have been resolved.
