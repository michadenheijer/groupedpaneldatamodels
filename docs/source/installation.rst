Installation
================

This package is distributed as pure Python with minimal dependencies, so in most environments a C/Fortran compiler is not required (all heavy dependencies ship pre‑built wheels).

Installation of the package is fairly straightforward for (most) machines with Python installed. As the package is published on the Python Package Index (PyPI), it can easily be installed using the following command.

.. code-block:: bash

    pip install groupedpaneldatamodels

If there are updates available can these be installed using:

.. code-block:: bash

    pip install --upgrade groupedpaneldatamodels

In most cases it is advised to run the second command, as it automatically upgrades the dependencies of your machine if the package requires a newer version of a certain dependencies. To download the latest development version of the package, it can simply be downloaded from Github and be installed as follows:

.. code-block:: bash

    git pull https://github.com/michadenheijer/groupedpaneldatamodels.git
    cd groupedpaneldatamodels
    pip install .

On most Windows machines the ``pip`` command may not be available, thus the alternative ``python -m pip`` may have to be used. In the case that multiple versions of Python are installed ``pip3`` may have to be used.

Supported Python versions
-------------------------

``groupedpaneldatamodels`` requires Python 3.11 or later and is tested on
CPython 3.12.

Dependencies
-------------
The package has the following dependencies:
- `numpy <https://numpy.org/>`_ – for numerical operations
- `scipy <https://scipy.org/>`_ – for statistical functions
- `statsmodels <https://www.statsmodels.org/>`_ – for the summary function
- `scikit-learn <https://scikit-learn.org/>`_
- `skglm <https://contrib.scikit-learn.org/skglm/stable/index.html>` - for the SCAD penalty
- `numba <https://numba.pydata.org/>`_ – for JIT compilation of some functions
- `tqdm <https://tqdm.github.io/>`_ – for progress bars during estimation