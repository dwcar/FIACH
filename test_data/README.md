# Packaged overt_0209 validation data

`overt_0209/input` contains 206 realigned three-dimensional NIfTI volumes and
the associated six-column SPM realignment file. `overt_0209/reference_R`
contains the 206 filtered volumes produced by the upstream R FIACH workflow.
`reference_matlab_diagnostics` records the mask and diagnostic maps used during
development; the validation test generates fresh MATLAB diagnostics.

The NIfTI repetition time is 1.25 seconds. The packaged validation uses B0 = 3
T, TE = 26 ms, a 128-second high-pass period, maximum gap 1, and 1.96 MADs.

These data are included at the request of the data custodian for a private
repository. Before adding collaborators, forking, publishing, or making the
repository public, confirm that the data are de-identified and that the
original consent and institutional data-governance terms permit that use.
