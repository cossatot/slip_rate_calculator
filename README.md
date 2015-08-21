Slip Rate Calculator
====================

The Slip Rate Calculator is a tool that takes probability distributions for
age and offset distance for geologic (or other) features cut by a fault and
calculates the probability distribution of the fault's slip rate. Typically,
the offset features (called `offset markers`) show *cumulative offset*, i.e.
all the `offset markers` that are present on a fault at a given time interval
experience the same offset amount over that time.

The Slip Rate Calculator can also calculate slip rates that change through
time, either by fitting a piecewise line (with a specified number of segments)
or a cubic spline. The piecewise fitting routine also includes a statistical
test (based on the [Bayesian Information Criterion][bic]) to determine whether
there was a slip rate change at some point in the past.

[bic]: https://en.wikipedia.org/wiki/Bayesian_information_criterion


## Installation
The Slip Rate Calculator is a Python 3 application. It depends heavily on the
Python scientific stack, as well as PyQt4. The easiest way to get everything
running is to install the free [Anaconda Python 3.4][anaconda] distribution. It
doesn't interfere with your system Python, or any other Python versions you
may have on your computer.

[anaconda]: https://store.continuum.io/cshop/anaconda/

Then, the git repository should be cloned:

```
git clone https://github.com/cossatot/slip_rate_calculator.git
```

The Slip Rate Calculator (which does need a more clever name) itself doesn't
need to be installed. Instead, go into the `slip_rate_calculator/app`
directory and run `python SlipRateGUI.pyw` or `./SlipRateGIU.pyw`.

Note: Because the Slip Rate Calculator has an embedded IPython console, it
can't be run using IPython. Plain Python works fine.

## Basic Usage

