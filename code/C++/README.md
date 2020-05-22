__C++__ implementation of the SpiSeMe package.

The package provides four surrogate-generating functions:
```cpp
spiSeMe_surrogate_jodi();
spiSeMe_surrogate_iaaft();
spiSeMe_surrogate_sa();
spiSeMe_surrogate_dither();
```
In addition, the package provides auxiliary functions to analyze event sequences:
```cpp
spiSeMe_event_autocorrelation();
spiSeMe_event_cross_correlation();
spiSeMe_distribution_test();
```
Function documentation can be found in `/docs/manual.pdf`.

### Setup

Functions have to be compiled under the C++11 standard.

The package requires the [GNU Scientific Libraries](https://www.gnu.org/software/gsl/) to run the C++ implementation of the functions.

Each function is provided within a dedicated source file (`.cpp`). Function declarations are reported in the header file `spiSeMe.h`.


### Examples

The `/examples/C++/` directory stores two C++ programs that provide examples of usage of the package functions. The directory also stores a `Makefile` that shows an example of compilation chain on Linux systems.
