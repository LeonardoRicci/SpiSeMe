Unit tests for the __Python__ implementation of the SpiSeMe package.

Two scripts for unit testing are provided:
```bash
test_surrogate_generation.py
test_auxiliary_functions.py
```

To invoke the tests, run (in this directory)
```bash
python3 test_surrogate_generation.py
python3 test_auxiliary_functions.py
```

Both tests will check the correct behavior of functions both for invalid and for valid arguments. For this reason, during the tests the console will display error statements coming from the spiSeMe functions intentionally fed with invalid arguments. If the tests succeed, the output will end with lines looking like the following ones:
```bash
----------------------------------------------------------------------
Ran 6 tests in 0.053s

OK
```

Please note that, because surrogate generation is not a deterministic process, it is impossible to provide a unit testing that checks the numerical output of functions. Nevertheless, wherever possible (JODI, IAAFT, SA), `test_surrogate_generation.py` checks that IEI distributions are exactly conserved.
