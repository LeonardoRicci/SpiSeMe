__MATLAB__ implementation of the SpiSeMe package.

The package provides four surrogate-generating functions:
```matlab
spiSeMe_surrogate_jodi
spiSeMe_surrogate_iaaft
spiSeMe_surrogate_sa
spiSeMe_surrogate_dither
```
In addition, the package provides auxiliary functions to analyze event sequences:
```matlab
spiSeMe_event_autocorrelation
spiSeMe_event_cross_correlation
spiSeMe_distribution_test
```
Function documentation can be found in `/docs/manual.pdf`.

### Setup

No setup is required. The functions source files (`.m`) have to be stored in a directory that is part of the Matlab _path_, or that is added to the Matlab path for the current session, for example by using
```matlab
codePath = '../../code/MATLAB';
addpath(codePath);
```

### Examples

The `/examples/MATLAB/` directory stores two Matlab scripts containing examples.
