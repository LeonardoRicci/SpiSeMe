### Check required binaries ###

ifeq (, $(shell which g++))
$(error "Compiler g++ not found. Modify the Makefile to change compiler.")
endif
COMPILER = g++
CXXFLAGS = -std=c++11
CFLAGS = -lm

ifeq (, $(shell which gsl-config))
$(error "GSL configuration executable `gsl-config` not found: cannot compile the example. Modify the Makefile if your GSL setup does not provide the `gsl-config` executable.")
endif
GSL_CFLAGS = `gsl-config --cflags`
# 	If gsl-config is not available, change with: -I/path/to/gsl/parent/directory
# 	e.g. if GSL headers are in /usr/incldude/gsl:
# GSL_CFLAGS = -I/usr/include
GSL_LINKFLAGS = `gsl-config --libs`
# 	If gsl-config is not available, change with: -L/path/to/directory/of/gsl/binaries -lgsl -lgslcblas
#	e.g. if /usr/lib/x86_64-linux-gnu contains GSL binaries (libgsl.a, etc.):
# GSL_LINKFLAGS = -L/usr/lib/x86_64-linux-gnu -lgsl -lgslcblas

CODE_DIR = ../../code/C++/
# Directory containing the functions source files (.cpp, .h). Change if necessary.


### Recipes ###

.PHONY: all
all: example example_cross_correlation

example: example.cpp
	cd $(CODE_DIR); $(COMPILER) -c spiSeMe_event_autocorrelation.cpp -o ac.o $(CXXFLAGS) $(CFLAGS) $(GSL_CFLAGS)
	cd $(CODE_DIR); $(COMPILER) -c spiSeMe_surrogate_jodi.cpp -o jodi.o $(CXXFLAGS) $(CFLAGS) $(GSL_CFLAGS)
	cd $(CODE_DIR); $(COMPILER) -c spiSeMe_surrogate_iaaft.cpp -o iaaft.o $(CXXFLAGS) $(CFLAGS) $(GSL_CFLAGS)
	cd $(CODE_DIR); $(COMPILER) -c spiSeMe_surrogate_sa.cpp -o sa.o $(CXXFLAGS) $(CFLAGS) $(GSL_CFLAGS)
	cd $(CODE_DIR); $(COMPILER) -c spiSeMe_surrogate_dither.cpp -o dither.o $(CXXFLAGS) $(CFLAGS) $(GSL_CFLAGS)
	mv $(CODE_DIR)*.o .
	$(COMPILER) example.cpp ac.o jodi.o iaaft.o sa.o dither.o $(CXXFLAGS) $(CFLAGS) $(GSL_LINKFLAGS) -I$(CODE_DIR) -o example

example_cross_correlation: example_cross_correlation.cpp
	cd $(CODE_DIR); $(COMPILER) -c spiSeMe_event_cross_correlation.cpp -o xc.o $(CXXFLAGS) $(CFLAGS) $(GSL_CFLAGS)
	cd $(CODE_DIR); $(COMPILER) -c spiSeMe_distribution_test.cpp -o kstest.o $(CXXFLAGS) $(CFLAGS) $(GSL_CFLAGS)
	cd $(CODE_DIR); $(COMPILER) -c spiSeMe_surrogate_jodi.cpp -o jodi.o $(CXXFLAGS) $(CFLAGS) $(GSL_CFLAGS)
	mv $(CODE_DIR)*.o .
	$(COMPILER) example_cross_correlation.cpp xc.o jodi.o kstest.o $(CXXFLAGS) $(CFLAGS) $(GSL_LINKFLAGS) -I$(CODE_DIR) -o example_cross_correlation

.PHONY: clean
clean:
	rm -f example example_cross_correlation
	rm -f ac.o xc.o kstest.o jodi.o iaaft.o sa.o dither.o
	rm -f ac_*.dat surr_iei_*.dat cross_correlation.dat
