# RVHRcorr Regressors for Denoising BOLD Data

The script **`RVHR_corr_regressors.m`** processes physiological data and scan parameters to generate convolution-based regressors using the respective physiological response functions.

## Dependencies
The script calls the following auxiliary functions:  
- **`RVTestimate.m`** – Estimates respiratory variance time series.  
- **`HRinterp.m`** – Interpolates heart rate data.  
- **`zhangHRF.m`** – Models the heart rate response function.  
- **`birnRRF.m`** – Models the respiratory response function.  

## Acknowledgments
Some of the functions were  written by **Harrison Fisher** from the **Napadow Laboratory** Spaulding Rehabilitation Hospital in 2020.

## RETROICOR Citation
This repository incorporates physiological noise correction using **RETROICOR**, method for removing cardiac and respiratory artifacts from BOLD data. RETROICOR was first introduced in:

**Glover, G. H., Li, T. Q., & Ress, D. (2000).** *Image-based method for retrospective correction of physiological motion effects in fMRI: RETROICOR.*  
*Magnetic Resonance in Medicine, 44*(1), 162–167.  
[https://doi.org/10.1002/1522-2594(200007)44:1<162::AID-MRM23>3.0.CO;2-E](https://doi.org/10.1002/1522-2594(200007)44:1<162::AID-MRM23>3.0.CO;2-E)
