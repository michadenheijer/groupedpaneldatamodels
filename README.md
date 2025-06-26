# ```groupedpaneldatamodels```: A Python Library for Grouped Fixed and Interactive Effects Models

```groupedpaneldatamodels``` is an open‑source Python library that implements a collection of *Grouped Panel Data Models* (GPDMs) for econometric research.

`groupedpaneldatamodels` is an Open Source Python library that implements multiple Grouped Panel Data Models (GPDMs) for Econometric research. These models offer a middle ground between fully homogeneous (which are often incorrectly specified) and fully heterogeneous (which are often difficult to estimate) by grouping multiple individuals and assuming the same coeficients for all members of the groupings.

## Features
This package implements the models and algorithms proposed by the following four papers, which each suggest different GPDMS.

* **Grouped Fixed Effects (GFE)**
  * [Bonhomme & Manresa (2015)](https://doi.org/10.3982/ECTA11319) clustering estimator
  * [Su, Shi & Phillips (2016)](https://doi.org/10.3982/ECTA12560) C‑Lasso estimator
* **Grouped Interactive Fixed Effects (GIFE)**
  * [Ando & Bai (2016)](https://doi.org/10.1002/jae.2467) clustering estimator
  * [Su & Ju (2018)](https://doi.org/10.1016/j.jeconom.2018.06.014) C‑Lasso estimator
* Automatic group selection via Information Criteria (BIC, AIC, HQIC).
* Analytical or bootstrap standard errors
* Fast NumPy and JIT-compiled Numba core with optional parallel bootstrap for large panels
* Familiar, `statsmodels`‑like API

## Installation

```bash
pip install groupedpaneldatamodels
# or update
pip install --upgrade groupedpaneldatamodels
```

To grab the bleeding‑edge version:

```bash
git clone https://github.com/michadenheijer/groupedpaneldatamodels.git
cd groupedpaneldatamodels
pip install .
```

---

## Quick start

```python
import numpy as np
import groupedpaneldatamodels as gpdm

# Y  shape (N, T, 1); X  shape (N, T, K)
gfe = gpdm.GroupedFixedEffects(Y, X, G=3, model="bonhomme_manresa")
gfe.fit()
print(gfe.summary())

gife = gpdm.GroupedInteractiveFixedEffects(Y, X, G=3,
                                           model="ando_bai",
                                           GF=[2, 2, 2])            # 2 common factors per grouping
gife.fit()
betas = gife.params["beta"]
```

### Selecting the number of groups

```python
best = gpdm.grid_search_by_ic(
    gpdm.GroupedFixedEffects,
    param_ranges={"G": range(1, 7)},
    init_params={"dependent": Y, "exog": X},
    pit_params={"gife_iterations": 100},
    ic_criterion="BIC"
)
print(best.G)          # optimal group count
```

## Documentation

An API reference with proper installation and guidelines is available at
<https://groupedpaneldatamodels.michadenheijer.com>


## Simulation Study

A simulation study has been done for the Master's thesis creating this package. This thesis has shown that this package
can succesfully reproduce the properties of the underlying estimators and can reduce the RMSE compared to a fully heterogeneous
model when `N` is large and `T` is small.

---

## Citation

Please cite the thesis if you use `groupedpaneldatamodels`:

```bibtex
@mastersthesis{denheijer2025,
    author      = {Micha den Heijer},
    title       = {groupedpaneldatamodels: A Python Library for Grouped Fixed and Interactive Effects Models},
    school      = {Vrije Universiteit Amsterdam},
    year        = {2025},
    month       = {July},
    date        = {2025-07-01},
    url         = {https://groupedpaneldatamodels.michadenheijer.com/_static/thesis.pdf}
}
```

---

## License

Released under the MIT License. See [LICENSE](LICENSE) for details.