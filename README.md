# Stability-dependant uniform shear model

Numerical implementation of the uniform shear model including stability by Chougule et al. (2018)

[![View Uniform shear model including atmospheric stability on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://se.mathworks.com/matlabcentral/fileexchange/74480-uniform-shear-model-including-atmospheric-stability)
[![DOI](https://zenodo.org/badge/249147073.svg)](https://zenodo.org/badge/latestdoi/249147073)

The uniform shear model for non-neutral, atmospheric-boundary-layer turbulence by  Chougule et al. (2018) [1] is here implemented in Matlab. The present implementation differs likely from the one used by ref. [1]. In particular, it relies on a more intuitive approach. Thus it should be referred to as [2]

The present submission contains:

- The function ChouguleTurb that implements the spectral tensor of non-neutral, atmospheric-boundary-layer turbulence.
- A function fitChougule.m that fits the stability-corrected uniform shear model to the target velocity spectra.
- An interactive example file Example.mlx that reproduces Fig. 1 in Chougule et al. [1]
- An interactive example file Example_fitting.mlx that illustrates how the function fitChougule.mcan be used


This is the second version of the submission. Typos and bugs may still be present. Any comment, suggestion or question is welcomed.


References:

[1] Chougule, A., Mann, J., Kelly, M., & Larsen, G. C. (2018). Simplification and validation of a spectral-tensor model for turbulence including atmospheric stability. Boundary-Layer Meteorology, 167(3), 371-397.

[2]  E. Cheynet (2022). Uniform shear model including atmospheric stability. Zenodo.  Retrieved yyyy-mm-dd.  doi:10.5281/zenodo.3774066
