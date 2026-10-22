# Heteroscedastic-and-homoscedastic-Shrunken-Centroid-Classifiers-Python


Heteroscedastic nearest shrunken centroid (NSC) classification for
high-dimensional, low-sample-size data such as gene expression. Python port of
the R package **hetNSC**.

| Label | Thresholding | Variance | Reference |
|-------|--------------|----------|-----------|
| STh   | soft  | pooled (homoscedastic) | original PAM (Tibshirani et al. 2002) |
| HTh   | hard  | pooled (homoscedastic) | Sahtout, Wang and Ghimire (2024) |
| OTh   | order | pooled (homoscedastic) | Sahtout, Wang and Ghimire (2024) |
| STh3  | soft  | class-specific | Sahtout and Wang |
| HTh3  | hard  | class-specific | Sahtout and Wang |
| OTh3  | order | class-specific | Sahtout and Wang |

The threshold is chosen by the refined ("deep") cross-validation search of
Sahtout, Wang and Ghimire (2024).

## Installation

```bash
pip install hetnsc-0.1.0-py3-none-any.whl          # from the built wheel
pip install "hetnsc[sklearn,pandas]"                # once on PyPI, with extras
```

The only required dependency is NumPy. scikit-learn (optional) makes
`HNSCClassifier` a full scikit-learn estimator; pandas (optional) makes the
comparison summaries data frames.

## Quick start

Data are **samples in rows, genes in columns** (the scikit-learn convention;
note that the R package uses genes in rows).

```python
import hetnsc as h

d = h.simulate_hetero_data(15, 20, 50, 50, p=1000, a=14, random_state=1)

clf = h.HNSCClassifier("HTh3", random_state=0)     # or method="hard", variance="hetero"
clf.fit(d["X_train"], d["y_train"])
print(clf)                                         # threshold, genes used, CV error
pred = clf.predict(d["X_test"])
print((pred != d["y_test"]).mean())                # test error
clf.predict_proba(d["X_test"])                     # posterior probabilities
clf.selected_genes_                                # column indices of genes used
```

`HNSCClassifier` works with scikit-learn tools such as `cross_val_score`,
`Pipeline` and `GridSearchCV`.

## Main functions

| Python | R (hetNSC) | Purpose |
|--------|-----------|---------|
| `HNSCClassifier(...).fit()` | `hnsc()` | Fit with the threshold chosen by deep CV |
| `.predict()`, `.predict_proba()`, `.decision_function()` | `predict()` | Classes, posteriors, scores |
| `train()` | `hnsc_train()` | Fit over a grid of thresholds |
| `cross_validate()` | `hnsc_cv()` | Class-balanced CV along the grid |
| `select_threshold()` | `hnsc_select()` | Deep CV threshold search |
| `compare_nsc()` | `compare_nsc()` | Repeated-run comparison protocol of the paper; `.summary()`, `.recall_table()` |
| `simulate_hetero_data()` | `simulate_hetero_data()` | Correlated-gene simulation designs (`correlated="B"` or `"both"`) |
| `compute_rho2()` | `compute_rho2()` | Relative heterogeneity / crossover check (Section 2.4) |
| `confusion_matrix()`, `per_class_recall()`, `balanced_accuracy()` | `make_cm()`, ... | Class-level metrics |
| `methods()` | `hnsc_methods()` | Table of the six classifiers |

`examples/simulation_correlated.py` runs the correlated-gene simulation study.

## Agreement with the R package

Given the same data and the same cross-validation folds, the Python and R
packages select the same thresholds and give the same CV errors, predictions
and posterior probabilities (to about 1e-12); `tests/test_against_r.py` checks
this against results recorded from R. Because NumPy and R use different random
number generators, the same seed does **not** give the same folds or simulated
data in the two languages, so results from repeated random runs agree only in
distribution.

One known source of tiny platform differences: with **two classes**,
|d_i1| = |d_i2| holds exactly in theory for every gene, so order thresholding
(OTh, OTh3) breaks these ties according to floating-point rounding, which can
differ between R, NumPy and BLAS libraries.

## Citation

* Sahtout, M. O. and Wang, H. Heteroscedastic nearest shrunken centroid
  classification: Methods and applications to high-dimensional data.
  *Statistical Applications in Genetics and Molecular Biology* (accepted).
* Sahtout, M. O., Wang, H. and Ghimire, S. (2024). Different thresholding
  methods on nearest shrunken centroid algorithm. *Communications in
  Statistics - Simulation and Computation*, 53(3), 1444-1460.
  https://doi.org/10.1080/03610918.2022.2047201
* Tibshirani, R., Hastie, T., Narasimhan, B. and Chu, G. (2002). Diagnosis of
  multiple cancer types by shrunken centroids of gene expression. *PNAS*,
  99(10), 6567-6572.

## License

GPL-2.0-or-later. Parts of the algorithm are adapted from the R package
**pamr** by Hastie, Tibshirani, Narasimhan and Chu.
