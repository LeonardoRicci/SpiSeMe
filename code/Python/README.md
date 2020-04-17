__Python__ implementation of the SpiSeMe package.

The package provides four surrogate-generating functions:
```python
spiSeMe.spiSeMe_surrogate_jodi
spiSeMe.spiSeMe_surrogate_iaaft
spiSeMe.spiSeMe_surrogate_sa
spiSeMe.spiSeMe_surrogate_dither
```
In addition, the package provides auxiliary functions to analyze event sequences:
```python
spiSeMe.spiSeMe_event_autocorrelation
spiSeMe.spiSeMe_event_cross_correlation
spiSeMe.spiSeMe_distribution_test
```
Function documentation can be found in `/docs/manual.pdf`.

### Setup

The package requires __Python version >= 3.5__ and __Numpy version >= 1.15__.

To install the package, run the following command from within the `/code/Python/` directory:
```Bash
python3 setup.py install --user
```
Depending on the system, it might be necessary to change `python3` with `python`. The `--user` flag might not be necessary, depending on the system. Please refer to the [Python packaging user guide](https://packaging.python.org/tutorials/installing-packages/) for further information on package setup.

From within Python, each module in the package has to be imported separately:
```python
from spiSeMe import <source> as ...
```
Each module provides the function with the same name. For example:
```python
from spiSeMe import spiSeMe_surrogate_jodi as jodi

jodi.spiSeMe_surrogate_jodi(...)
```

### Examples

The `/examples/Python/` directory stores two Python scripts containing examples.

### Uninstall

To remove the package use `pip`:

```Bash
pip3 uninstall SpiSeMe
```
