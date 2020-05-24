Unit tests for the __Matlab__ implementation of the SpiSeMe package.

Two scripts for unit testing are provided:
```bash
surrogate_generation_test.m
auxiliary_functions_test.m
```

To invoke the tests, run (within Matlab, from this directory)
```matlab
results = runtests('surrogate_generation_test.m')
results = runtests('auxiliary_functions_test.m')
```

Both tests will check the correct behavior of functions both for invalid and for valid arguments. If the tests succeed, the output will begin and end with lines looking like the following ones:
```matlab
Running surrogate_generation_test
........
Done surrogate_generation_test

[...]

Totals:
   8 Passed, 0 Failed, 0 Incomplete.
   0.61618 seconds testing time.
```

Please note that, because surrogate generation is not a deterministic process, it is impossible to provide a unit testing that checks the numerical output of functions. Nevertheless, wherever possible (JODI, IAAFT, SA), `surrogate_generation_test.m` checks that IEI distributions are exactly conserved.
