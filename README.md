# SpiSeMe

The **Spi**ke **Se**quence mi**me** (SpiSeMe) package provides **C++**, **Matlab** and **Python** functions implementing four different algorithms for generation of surrogates of event (aka spikes) sequences.
The algorithms provided by the package are
* the **JODI** (**JO**int **DI**stribution) surrogate generation algorithm described in L. Ricci, M. Castelluzzo, L. Minati and A. Perinelli _Generation of surrogate event sequences via joint distribution of successive inter-event intervals_, Chaos **29**(12):121102, 2019; [DOI:10.1063/1.5138250](https://doi.org/10.1063/1.5138250);
* the **IAAFT** (**i**terative **a**mplitude **a**djusted **F**ourier **t**ransform) algorithm first proposed by T. Schreiber and A. Schmitz in Phys. Rev. Lett. **77** (1996), 635;
* a generalized randomization algorithm based on simulated annealing (**SA**), first proposed by T. Schreiber in Phys. Rev. Lett. **80** (1998), 2105;
* a dithering method.

In addition, the package also provides a function for the assessment of the autocorrelation of event sequences, a function for the assessment of cross correlation between event sequences, and a function to test the compatibility of the inter-event interval distributions of two event sequences.

A summary of the algorithms can be found in Chapter 2 of the user manual (`/docs/manual.pdf`).


## License

This package is free software. It is distributed under the terms of the GNU General Public License (GPL), version 3.0 - see the LICENSE.txt file for details.


## Authors

- Alessio Perinelli (1), alessio.perinelli@unitn.it
- Michele Castelluzzo (1)
- Ludovico Minati (2,3,4)
- Leonardo Ricci (1,2), leonardo.ricci@unitn.it

(1) Department of Physics, University of Trento, Trento, Italy  
(2) CIMeC, Center for Mind/Brain Sciences, University of Trento, Rovereto, Italy  
(3) Tokyo Tech World Research Hub Initiative (WRHI), Institute of Innovative Research, Tokyo Institute of Technology,
Yokohama 226-8503, Japan  
(4) Complex Systems Theory Department, Institute of Nuclear Physics, Polish Academy of Sciences, 31-342 Kraków, Poland

If the implementation of the JODI algorithm turns out to be useful for your research, please cite our paper:
	L. Ricci, M. Castelluzzo, L. Minati and A. Perinelli _Generation of surrogate event sequences via joint distribution of successive inter-event intervals_, Chaos **29**(12):121102,


## C++
The C++ implementation of the package requires the [GNU Scientific Libraries](https://www.gnu.org/software/gsl/) to properly compile and run. Functions are provided in source files (`.cpp`) and their declaration is reported in the single header file `spiSeMe.h`. The functions documentation is available in Chapter 3 of the user manual `/docs/manual.pdf`.

## Matlab
The Matlab implementation of the package was developed and tested on Matlab version R2017b. Previous versions of Matlab might not be able to run the package. The package does not require any setup. To make the package functions available in Matlab, the source (`.m`) files have to be copied in a directory that is part of Matlab _path_. For further details, please refer to the [official Mathworks website/user guide](https://www.mathworks.com/help/matlab/matlab_env/what-is-the-matlab-search-path.html). The functions documentation is available in Chapter 4 of the user manual `/docs/manual.pdf`.

## Python
The Python implementation requires Python 3.5 and NumPy 1.15 or later. The package might not work properly with older versions of Python or NumPy. The examples included in the package also require the Matplotlib library in order to produce graphic output. Examples were tested with Matplotlib version 3.0.3. Further information can be found in `/code/Python/README.md` The functions documentation is available in Chapter 5 of the user manual `/docs/manual.pdf`, chapter
