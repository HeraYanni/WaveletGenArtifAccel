# Wavelet-based stochastic model for the generation of fully non-stationary seismic accelerograms

Tool for generating target spectrum-compatible fully-nonstationary artificial seismic ground motions using the Continuous Wavelet Transform and a seed record. 

Further details are provided in the following document:

Yanni H., Fragiadakis M., and Mitseas I.P. 
"Wavelet-based stochastic model for the generation of fully non-stationary bidirectional seismic accelerograms". 
Earthquake Engineering and Structural Dynamics. 
https://doi.org/10.1002/eqe.4315

The model operates in the time-frequency domain and combines spectral representation techniques with signal processing tools. The basis of the methodology involves the generation of spectrum-compatible stationary artificial accelerogram signals whose non-stationarity is then modeled with a time-frequency modulating function that is based on a seed ground motion record. 
At the core of the proposed methodology lies the use of the Continuous Wavelet Transform (CWT). Specifically, the CWT method is used to perform time-frequency analysis and to define the non-stationary component. The proposed methodology provides any required number of seismic accelerograms whose temporal and spectral modulation is consistent with the characteristics of the site of interest.

Version 1.0 created by Hera Yanni, first release: 11th of November, 2024.

## How to run

* Download the files and run the main MAIN_gen_single_accel.m in MATLAB

* Run the MAIN_gen_single_accel.m in MATLAB online (no license is needed for MATLAB online) by clicking on this link https://matlab.mathworks.com/open/github/v1?repo=HeraYanni/WaveletGenArtifAccel,

or pressing this button Open in MATLAB Online 
[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=HeraYanni/WaveletGenArtifAccel)

Needs to login to a MATLAB account, the working interface is this:

![image](https://github.com/user-attachments/assets/b8fb67b3-2cec-4b2b-9df1-deb084f519ce)

## Copyright Notice

The present software files, created by Hera Yanni, can be freely employed, downloaded, and shared with others as long as proper credit is attributed to the creator:

Hera Yanni (c) 2024

Structural Engineer NTUA, MSc in ADERS

Ph.D. Candidate, Laboratory for Earthquake Engineering NTUA

email: heragian@mail.ntua.gr, heragian@gmail.com

Please note that there is no warranty or liability of the creator associated with this free software.
