Installation
================

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