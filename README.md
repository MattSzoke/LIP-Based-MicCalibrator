# LIP-Based-MicCalibrator
This repository demonstrates how to process sound generated using laser-induced plasma (LIP) for microphone calibration purposes. 


We demonstrate how to pre-process LIP acoustic data using Matlab. 
The following processing steps are demonstrated:
- arrange data into blocks,
- determine if LIP was formulated or not, 
- find arrival time, 
- positive and negative peak pressure determination, 
- time span determination, 
- gating and zero padding, 
- Fourier transform gated and padded signal
- scaling the signal,
- distance correction determination,
- microphone calibration, 
- microphone response function determination. 

The process is split up to two portions. First, run Script1_ReduceLIPData.m and then run Script2_CalibrationGenerator.m.
As their names suggest, Script1 digests data and Script2 determines microphone responses. 
Script1 relies on a sample dataset, which is available using this link:
https://www.dropbox.com/s/yx3h1ezuule0xx5/SampleData_v2.mat?dl=0

For any questions or comments, please do not hesitate to reach out!

Matt Szoke
m.szoke@vt.edu
January 2023.
