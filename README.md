# Supplementary code for "Alternative EEG pre-processing pipelines can lead to conflicting conclusions regarding cortical excitation/inhibition ratio" by Racz et al. (2025).

This repository will contain the complete analysis pipeline for the aformentioned project, including scripts reproducing the main analysis outcomes, scripts performing all statistical analyses reported in the manuscript, and scripts creating figures published in the manuscript and its supplementary material.

The following toolboxes were utilized during the analysis pipelines:
- Irregular Resampling Auto-Spectral Analysis (IRASA) by Wen & Liu (Wen, H. and Liu, Z., 2016. Separating fractal and oscillatory components in the power spectrum of neurophysiological signal. Brain topography, 29, pp.13-26.). **NOTE**: this toolbox had been extended by an additional function (_func_IRASA_plawfit_multi.m_) to allow for estimating the spectral slope from multiple frequency ranges simultaneously.
- Fitting oscillatons & one over f (FOOOF) or _spectral parametrization_ toolbox by Donoghue et al. (Donoghue, T., Haller, M., Peterson, E.J., Varma, P., Sebastian, P., Gao, R., Noto, T., Lara, A.H., Wallis, J.D., Knight, R.T. and Shestyuk, A., 2020. Parameterizing neural power spectra into periodic and aperiodic components. Nature neuroscience, 23(12), pp.1655-1665.)
- MultiModal Spectral Parametrization Tool (MMSPM), a method for estimating the spectral slope in a model-free, multimodal fashion using constrained piece-wise regression. This method is currently being developed by the authors and has not been published yet.

The repository contains data from two individual datasets:
- Cohort #1 (ErrCap dataset): EEG data collected with a custom electrode montage with high electrode density (10-5 standard) over the midline and parasagittal lines, such as illustrated on Fig. 4 in (Kumar, S., Liu, D.H., Racz, F.S., Retana, M., Sharma, S., Iwane, F., Murphy, B.P., O'Keeffe, R., Atashzar, S.F., Alambeigi, F. and Millán, J.D.R., 2023, May. Cognidavinci: Towards estimating mental workload modulated by visual delays during telerobotic surgery-an eeg-based analysis. In 2023 IEEE International Conference on Robotics and Automation (ICRA) (pp. 6789-6794). IEEE.)
- Cohort #2 (StimCap dataset): EEG data collected from 32 standard 10-10 electrode locations, such as in (Liu, D.H., Kumar, S., Alawieh, H., Racz, F.S. and del R Millán, J., 2025. Personalized µ-transcranial alternating current stimulation improves online brain–computer interface control. Journal of neural engineering, 22(1), p.016037.) 

Currently, the repository is only for **review and reproducibility purposes** and therefore contains the following:

- Outputs of IRASA analysis for each participant, both datasets
- Outputs of FOOOF analysis for each participant, both datasets
- Outputs of MMSPM analysis for each participant, both datasets
- Scripts performing statistical analysis and generating of figures for all three analysis techniques, both datasets
- Corresponding custom functions
- IRASA toolbox
- FOOOF toolbox

**When using this resource, please cite** Racz, F.S., Mac-Auliffe, D., Mukli, P., Milton, J., Cabrera, J.L. and Millán, J.D.R., 2024. Alternative EEG pre-processing pipelines can lead to conflicting conclusions regarding cortical excitation/inhibition ratio. bioRxiv, pp.2024-11. DOI:https://doi.org/10.1101/2024.11.08.622698

Frigyes Samuel Racz
