# FIACH BOLD contrast parameter choices

## What this parameter controls

`fiach_bold_contrast(B0, TEms, ...)` ports upstream `FIACH::boldContrast`. FIACH passes its result to the deterministic temporal outlier threshold in `fiach`. It is a biophysical ceiling used to distinguish plausible signal changes from large artefacts; it is not an estimate of the task response in this particular run.

The implementation combines a baseline and activated effective transverse relaxation rate with a mono-exponential signal model. The contrast depends on field strength and echo time, vessel orientation, haematocrit, baseline and activated oxygen extraction, flow, and the assumed flow-to-volume exponent. Thus, a single contrast value is not portable across acquisitions.

## What the packaged data can and cannot tell us

The `overt_0209` validation package documents TR = 1.25 s and uses B0 = 3 T and TE = 26 ms. It includes realigned EPI volumes, motion estimates, and R-FIACH filtered reference volumes. It does not include a BIDS JSON acquisition sidecar, field-map or calibration acquisition, ASL CBF, blood haematocrit, measured oxygen extraction, or task-state physiological measurements. The B0 and TE here are the validation script's declared acquisition inputs, not estimates derived from the image intensities. The task contrast and the R output volumes cannot identify the physiological parameters in this forward model.

The available run therefore lets us calculate and compare model-implied thresholds at the specified B0 and TE, but it does not justify fitting an individualized threshold. Keep the R-compatible settings when exact reproduction of earlier FIACH analyses is the priority; choose a named profile when the assumed physiology is part of the analysis question.

## Named profiles

Values not listed in the table are held at the upstream FIACH defaults: random vessel orientation, haematocrit 0.40, susceptibility difference `1.8e-7` cgs, baseline OEF 0.40, and proton gyromagnetic ratio `2.675e8 rad/s/T`. `CBFBase` is in ml/100 g/min. The upstream function's documentation says ml per mg per minute, which is a unit typo; its numerical default of 55 corresponds to the conventional ml/100 g/min scale.

| Profile | alpha | Baseline CBF | Activated OEF | Activated CBF | Interpretation |
|---|---:|---:|---:|---:|---|
| `Typical` | 0.38 | 55 | 0.20 | 110 (2.0x) | FIACH's suggested less conservative OEF assumption; leaves its large, twofold flow increase in place. This is a practical middle setting within the FIACH model, not a validated population mean. |
| `Conservative` | 0.38 | 55 | 0.10 | 110 (2.0x) | The upstream FIACH defaults. The lower activated OEF and doubled flow maximize the modeled deoxyhaemoglobin reduction and give the more permissive theoretical signal-change ceiling described in FIACH. |
| `AlternativeTypical` | 0.29 | 50 | 0.33 | 68.15 (1.363x) | Alternative human-activation scenario: the CBF and OEF values are based on a published calibrated-fMRI activation example; alpha 0.29 is a later human average, rather than Grubb's original 0.38. This is an illustrative literature-based scenario, not a single internally fitted or universally representative population model. |

For the validation acquisition (3 T, TE 26 ms), the MATLAB port returns approximately 3.95, 7.87, and 0.77, respectively, on its upstream FIACH threshold scale for `Typical`, `Conservative`, and `AlternativeTypical`. These are model outputs, not values measured from the overt-speaking dataset. They are strongly profile-dependent, particularly because activated OEF and flow jointly determine the assumed deoxyhaemoglobin reduction.

### A units/normalization point to keep visible

The source-compatible routine computes `100 * (S_active - S_base)` when both signals have first been normalized to 100 at TE = 0. Conventional fractional signal change at the acquisition echo time is `100 * (S_active - S_base) / S_base`. At 3 T and TE 26 ms, the latter is approximately 6.66%, 13.27%, and 1.23% for the three profiles in the table. The returned values and these conventional percentages are not interchangeable. This MATLAB port retains the upstream returned-value convention to preserve FIACH compatibility; before interpreting the `fiach` threshold as a literal conventional percent signal change, compare this normalization with the exact upstream R implementation and intended threshold definition. This distinction should be resolved explicitly before using the threshold for a new scientific claim.

## Rationale for each input

- **B0 and TE:** pass the acquisition's actual field strength and echo time. The validation configuration is 3 T / 26 ms. BOLD amplitude and calibrated-BOLD parameters depend on TE, field strength, and vascular weighting; the calibration literature cautions against comparing M values across different acquisition conditions.
- **Baseline flow and OEF:** FIACH uses 55 and 0.40. Typical awake-human modelling values are around 50 ml/100 g/min and OEF 0.40. We retain 55 for the first two profiles to match FIACH; the alternative uses 50.
- **Flow-to-volume exponent (`Alpha`):** 0.38 is the classic Grubb et al. empirical relation. Later human measurements suggest a lower mean near 0.29; values can be lower still when specifically modelling deoxygenated venous blood. Alpha is therefore a population and vascular-compartment assumption, not a constant measured in this dataset.
- **Activated flow and OEF:** FIACH's defaults assume a twofold flow increase and OEF falling from 0.40 to 0.10; its R help explicitly calls 0.10 conservative and notes 0.20 as more realistic. In a calibrated human visual-activation experiment, Hoge et al. reported CBF +36.3% and an estimated OEF of 33% from a 40% baseline. That motivates `AlternativeTypical`; task, region, stimulus, and subject differences make these values illustrative rather than universal.
- **Haematocrit and susceptibility:** 0.40 and `1.8e-7` cgs are upstream model assumptions. No person-specific haematocrit or blood-susceptibility measurement is supplied with the test data, so these remain fixed across profiles.
- **Vessel geometry:** random orientation is the upstream default. The alternative parallel-cylinder assumption is useful as a sensitivity analysis, but the bundled data do not identify vessel geometry.

## Use

```matlab
tTypical = fiach_bold_contrast(3, 26, 'Profile', 'Typical');
tCeiling = fiach_bold_contrast(3, 26, 'Profile', 'Conservative');
tAlternative = fiach_bold_contrast(3, 26, 'Profile', 'AlternativeTypical');

result = fiach(functionalFiles, tTypical, 1.25, ...
    'RP', motionFile, 'GMMMethod', 'R');
```

Explicit `Alpha`, `CBFBase`, `EAct`, or `CBFAct` values override the corresponding selected profile values. The `Random`, `Hct`, `Chi`, `EBase`, and `W0` options continue to set the shared model assumptions.

## References

1. Tierney TM, Weiss-Croft LJ, Centeno M, et al. FIACH: A biophysical model for automatic retrospective noise control in fMRI. *NeuroImage*. 2016;124(Pt A):1009-1020. doi:[10.1016/j.neuroimage.2015.09.034](https://doi.org/10.1016/j.neuroimage.2015.09.034). [Open-access UCL record and paper](https://discovery.ucl.ac.uk/id/eprint/1475208/).
2. Tierney TM et al. `boldContrast` documentation in the FIACH R package. [Upstream function documentation](https://rdrr.io/github/tierneytim/FIACH/man/boldContrast.html). Documents defaults alpha = 0.38, hct = 0.4, baseline flow = 55, susceptibility difference = `1.8e-7`, baseline OEF = 0.4, activated OEF = 0.1, activated flow = twice baseline, and the note that 0.2 may give a more realistic threshold.
3. Grubb RL, Raichle ME, Eichling JO, Ter-Pogossian MM. The effects of changes in PaCO2 on cerebral blood volume, blood flow, and vascular mean transit time. *Stroke*. 1974;5(5):630-639. doi:[10.1161/01.STR.5.5.630](https://doi.org/10.1161/01.STR.5.5.630). Classic CBV-CBF exponent (0.38), established in rhesus monkeys during CO2 manipulation.
4. Ito H, Kanno I, Ibaraki M, Hatazawa J, Miura S. Changes in human cerebral blood flow and cerebral blood volume during hypercapnia and hypocapnia measured by positron emission tomography. *Journal of Cerebral Blood Flow & Metabolism*. 2003;23(6):665-670. doi:[10.1097/01.WCB.0000067721.64998.F5](https://doi.org/10.1097/01.WCB.0000067721.64998.F5). Reports CBV = 1.09 CBF^0.29 in healthy humans.
5. Hoge RD, Atkinson J, Gill B, Crelier GR, Marrett S, Pike GB. Coupling of cerebral blood flow and oxygen consumption during physiological activation and deactivation measured with fMRI. *NeuroImage*. 2004;23(1):148-155. doi:[10.1016/j.neuroimage.2004.05.013](https://doi.org/10.1016/j.neuroimage.2004.05.013). Reports average activation changes of 4.4% BOLD and 36.3% CBF, with estimated OEF decreasing from an assumed 40% baseline to 33%.
6. Buxton RB. The physics of functional magnetic resonance imaging (fMRI). *Reports on Progress in Physics*. 2013;76:096601. [Open-access review](https://pmc.ncbi.nlm.nih.gov/articles/PMC4376284/). Reviews the dependence of BOLD on flow, volume, oxygen metabolism, field strength, and echo time and gives common illustrative activation values.
7. Gauthier CJ, Hoge RD. A general review of calibrated fMRI models and interpretation. *NeuroImage*. 2013. [Calibrated fMRI review](https://www.sciencedirect.com/science/article/abs/pii/S1053811912001991). Notes that calibration M is acquisition-specific and summarizes reported 3 T M values and their variation.
