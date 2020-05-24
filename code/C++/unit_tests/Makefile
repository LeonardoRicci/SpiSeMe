### Check required binaries ###

ifeq (, $(shell which g++))
$(error "Compiler g++ not found. Modify the Makefile to change compiler. Tests cannot be run.")
endif
COMPILER = g++
CXXFLAGS = -std=c++11
CFLAGS = -lm

ifeq (, $(shell which gsl-config))
$(error "GSL configuration executable `gsl-config` not found: cannot compile the tests. Modify the Makefile if your GSL setup does not provide the `gsl-config` executable.")
endif
GSL_CFLAGS = `gsl-config --cflags`
# 	If gsl-config is not available, change with: -I/path/to/gsl/parent/directory
# 	e.g. if GSL headers are in /usr/incldude/gsl:
# GSL_CFLAGS = -I/usr/include
GSL_LINKFLAGS = `gsl-config --libs`
# 	If gsl-config is not available, change with: -L/path/to/directory/of/gsl/binaries -lgsl -lgslcblas
#	e.g. if /usr/lib/x86_64-linux-gnu contains GSL binaries (libgsl.a, etc.):
# GSL_LINKFLAGS = -L/usr/lib/x86_64-linux-gnu -lgsl -lgslcblas

CODE_DIR = ../
# Directory containing the functions source files (.cpp, .h). Change if necessary.


### Recipes ###

.PHONY: all
all:
	@echo "To perform the tests, type either:"
	@echo "    test-surrogate    [tests all surrogate generating functions]"
	@echo "    test-auxiliary    [tests all auxiliary functions]"

.PHONY: test-surrogate
test-surrogate: test_surrogate_generation
	@echo "\n\nCompiled successfully: now starting tests...\n"
	./test_surrogate_generation
	rm -f test_surrogate_generation

.PHONY: test-auxiliary
test-auxiliary: test_auxiliary_functions
	@echo "\n\nCompiled successfully: now starting tests...\n"
	./test_auxiliary_functions
	rm -f test_auxiliary_functions

test_surrogate_generation: test_surrogate_generation.cpp
	cd $(CODE_DIR); $(COMPILER) -c spiSeMe_surrogate_jodi.cpp -o jodi.o $(CXXFLAGS) $(CFLAGS) $(GSL_CFLAGS)
	cd $(CODE_DIR); $(COMPILER) -c spiSeMe_surrogate_iaaft.cpp -o iaaft.o $(CXXFLAGS) $(CFLAGS) $(GSL_CFLAGS)
	cd $(CODE_DIR); $(COMPILER) -c spiSeMe_surrogate_sa.cpp -o sa.o $(CXXFLAGS) $(CFLAGS) $(GSL_CFLAGS)
	cd $(CODE_DIR); $(COMPILER) -c spiSeMe_surrogate_dither.cpp -o dither.o $(CXXFLAGS) $(CFLAGS) $(GSL_CFLAGS)
	cd $(CODE_DIR); $(COMPILER) -c spiSeMe_event_autocorrelation.cpp -o ac.o $(CXXFLAGS) $(CFLAGS) $(GSL_CFLAGS)
	cd $(CODE_DIR); $(COMPILER) -c spiSeMe_distribution_test.cpp -o kstest.o $(CXXFLAGS) $(CFLAGS) $(GSL_CFLAGS)
	mv $(CODE_DIR)*.o .
	$(COMPILER) test_surrogate_generation.cpp jodi.o iaaft.o sa.o dither.o ac.o kstest.o $(CXXFLAGS) $(CFLAGS) $(GSL_LINKFLAGS) -I$(CODE_DIR) -o test_surrogate_generation
	rm -f *.o

test_auxiliary_functions: test_auxiliary_functions.cpp
	cd $(CODE_DIR); $(COMPILER) -c spiSeMe_event_autocorrelation.cpp -o ac.o $(CXXFLAGS) $(CFLAGS) $(GSL_CFLAGS)
	cd $(CODE_DIR); $(COMPILER) -c spiSeMe_event_cross_correlation.cpp -o xc.o $(CXXFLAGS) $(CFLAGS) $(GSL_CFLAGS)
	cd $(CODE_DIR); $(COMPILER) -c spiSeMe_distribution_test.cpp -o kstest.o $(CXXFLAGS) $(CFLAGS) $(GSL_CFLAGS)
	mv $(CODE_DIR)*.o .
	$(COMPILER) test_auxiliary_functions.cpp ac.o xc.o kstest.o $(CXXFLAGS) $(CFLAGS) $(GSL_LINKFLAGS) -I$(CODE_DIR) -o test_auxiliary_functions
	rm -f *.o

.PHONY: clean
clean:
	rm -f test_surrogate_generation test_auxiliary_functions
	rm -f *.o
	rm -f test*.log
