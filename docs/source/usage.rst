Usage Guide
===========

This guide shows how to work with **groupedpaneldatamodels** after installation.
It follows the usage section of the MSc thesis and focuses only on what an end‑user
needs to run, fit, and inspect the bundled estimators.

.. contents::
   :local:
   :depth: 2


Quick start
-----------

First import the library and prepare the data:

.. code-block:: python

    import numpy as np
    import groupedpaneldatamodels as gpdm

    # TODO replace with some sort of real data set

All estimators expect **3‑D NumPy arrays**:
``dependent`` with shape ``(N, T, 1)`` and ``exog`` with shape ``(N, T, K)``.
Pandas `DataFrame` support may be added later, but is **not** available yet.


Grouped Fixed Effects
---------------------

Two grouped‑fixed‑effects estimators are available:

* **Bonhomme & Manresa (2015)** – clustering approach (default)
* **Su, Shi & Phillips (2016)** – C‑Lasso penalisation

Both are exposed through one class::

    model_gfe = gpdm.GroupedFixedEffects(
        dependent=y,
        exog=x,
        G=3,                       # maximum number of groups
        use_bootstrap=False        # analytical s.e. by default
    )
    model_gfe.fit(max_iter=100, tol=1e-4)
    model_gfe.summary()

Common constructor arguments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
   :widths: 20 15 10 55
   :header-rows: 1

   * - Argument
     - Type
     - Default
     - Purpose
   * - ``dependent``
     - ``ndarray``
     - —
     - 3‑D response array ``(N, T, 1)``
   * - ``exog``
     - ``ndarray``
     - —
     - 3‑D regressor array ``(N, T, K)``
   * - ``G``
     - ``int``
     - —
     - Maximum number of latent groups to consider
   * - ``use_bootstrap``
     - ``bool``
     - ``False``
     - If *True* compute bootstrap s.e.'s after fitting
   * - ``model``
     - ``str``
     - ``"bonhomme_manresa"``
     - Choose ``"bonhomme_manresa"`` or ``"su_shi_phillips"``
   * - ``heterogeneous_beta``
     - ``bool``
     - ``True``
     - Make the slope coefficients group‑specific
   * - ``entity_effects``
     - ``bool``
     - ``False``
     - Include individual fixed effects :math:`\\alpha_i`

Model‑specific fit options
^^^^^^^^^^^^^^^^^^^^^^^^^^

Bonhomme & Manresa
^^^^^^^^^^^^^^^^^^

.. list-table::
   :widths: 25 15 15
   :header-rows: 1

   * - Parameter
     - Type
     - Default
   * - ``n_boot``
     - ``int``
     - ``50``
   * - ``max_iter``
     - ``int``
     - ``10000``
   * - ``tol``
     - ``float``
     - ``1e-6``
   * - ``gfe_iterations``
     - ``int``
     - ``100``
   * - ``enable_vns``
     - ``bool``
     - ``False``

Su, Shi & Phillips
^^^^^^^^^^^^^^^^^^

.. list-table::
   :widths: 25 15 15
   :header-rows: 1

   * - Parameter
     - Type
     - Default
   * - ``kappa``
     - ``float``
     - ``0.1``
   * - ``n_boot``
     - ``int``
     - ``50``
   * - ``max_iter``
     - ``int``
     - ``1000``
   * - ``tol``
     - ``float``
     - ``1e-6``
   * - ``only_bfgs``
     - ``bool``
     - ``True``

Inspecting results
^^^^^^^^^^^^^^^^^^

.. code-block:: python

    beta = model_gfe.params                       # coefficients
    se   = model_gfe.params_standard_errors       # bootstrap or analytical
    ci   = model_gfe.get_confidence_intervals()   # 95 % by default

To see more details about this model, view the :class:`groupedpaneldatamodels.GroupedFixedEffects` documentation.


Grouped Interactive Fixed Effects
---------------------------------

Interactive estimators extend the fixed‑effects ideas with unobserved common
factors.  Two variants are implemented:

* **Ando & Bai (2016)** – clustering (default)
* **Su & Ju (2018)** – C‑Lasso penalisation

Typical workflow::

    model_gife = gpdm.GroupedInteractiveFixedEffects(
        dependent=y,
        exog=x,
        G=3,
        use_bootstrap=False
    )
    model_gife.fit(max_iter=100, tol=1e-4)
    model_gife.summary()

Constructor arguments
^^^^^^^^^^^^^^^^^^^^^

Only the additions relative to the fixed‑effects constructor are shown.

.. list-table::
   :widths: 20 15 20
   :header-rows: 1

   * - Argument
     - Type
     - Default
   * - ``model``
     - ``str``
     - ``"ando_bai"``
   * - ``R``
     - ``int``
     - ``1``
   * - ``GF``
     - ``array``
     - ``[1,…,1]``

Fit‑time options match those in the thesis:

* ``n_boot``, ``max_iter``, ``tol``, ``gife_iterations``
* ``kappa`` and ``gamma`` (SCAD penalty, Ando & Bai)
* ``only_bfgs`` (Su & Ju)

Results are accessed the same way as for the fixed‑effects estimator.

To see more details about this model, view the :class:`groupedpaneldatamodels.GroupedInteractiveFixedEffects` documentation.

Information‑criterion model search
----------------------------------

Use ``grid_search_by_ic`` to choose *G* (or any other hyper‑parameter)
automatically:

.. code-block:: python

    best_model, best_ic = gpdm.grid_search_by_ic(
        gpdm.GroupedFixedEffects,
        param_ranges={"G": [2, 3, 4, 5, 6]},
        init_params={"dependent": y, "exog": x},
        fit_params={"max_iter": 100},
        ic_criterion="BIC"
    )

The helper tests every candidate, computes the requested criterion
(**BIC**, **AIC** or **HQIC**) and returns the best‑scoring fitted model.

To see more details about this function, view the :func:`groupedpaneldatamodels.grid_search_by_ic` documentation.


Further reading
---------------

You can find more details in the :ref:`API documentation <api>`. Additionally, you can read the thesis this package is based on,
which is which is `available as a PDF <_static/thesis.pdf>`_.