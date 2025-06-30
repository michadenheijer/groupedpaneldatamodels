.. groupedpaneldatamodels documentation master file, created by
   sphinx-quickstart on Thu Jun 26 00:16:32 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


``groupedpaneldatamodels`` documentation
====================================

.. image:: https://img.shields.io/pypi/v/groupedpaneldatamodels.svg?style=flat
   :target: https://pypi.org/project/groupedpaneldatamodels/
.. image:: https://github.com/michadenheijer/groupedpaneldatamodels/actions/workflows/sphinx.yml/badge.svg
   :target: https://github.com/michadenheijer/groupedpaneldatamodels/actions/workflows/sphinx.yml

``groupedpaneldatamodels`` is the first package to bring Grouped Fixed Effects (GFE) and Grouped Interactive Fixed Effects (GIFE) estimators
(Bonhomme & Manresa 2015; Su, Shi & Phillips 2016; Ando & Bai 2016; Su & Ju, 2018) to Python.
Compared with classic fixed effects, grouping boosts efficiency, uncovers hidden groupings, and scales to large `N` panels.
This package implements some of the most popular GFE and GIFE estimators, analytical or bootstrap standard errors, and automatic selection for hyperparameters such as the number of groupings using Information Criterion.

.. code-block:: python

    import groupedpaneldatamodels as gpdm
    model = gpdm.GroupedFixedEffects(dependent=Y, exog=X, G=3)
    model.fit().summary()     # prints coefficients, group sizes and IC

New here? Start with :doc:`usage` for a quick start guide, then you can dive into the API docs.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   usage
   citation
   modules
   changelog