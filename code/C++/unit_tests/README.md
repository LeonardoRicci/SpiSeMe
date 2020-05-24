Unit tests for the __C++__ implementation of the SpiSeMe package. The chosen unit testing framework is __Catch2__, which can be found [here](https://github.com/catchorg/Catch2), and whose header file `catch.hpp` is included in the present folder. __Catch2__ is licensed under the Boost Software License 1.0.

Two programs for unit testing are provided:
```bash
test_surrogate_generation.cpp
test_auxiliary_functions.cpp
```

A Makefile is provided for a quick compilation and run of the tests. To invoke the tests, use
```bash
make test-surrogate
make test-auxiliary
```

Both tests will check the correct behavior of functions both for invalid and for valid arguments. For this reason, during the tests the console will display error statements coming from the spiSeMe functions intentionally fed with invalid arguments. If the tests succeed, the output will begin and end with lines looking like the following ones:
```bash
===============================================================================
All tests passed (28 assertions in 1 test case)
```

Please note that, because surrogate generation is not a deterministic process, it is impossible to provide a unit testing that checks the numerical output of functions. Nevertheless, wherever possible (JODI, IAAFT, SA), `test_surrogate_generation.cpp` checks that IEI distributions are exactly conserved.
