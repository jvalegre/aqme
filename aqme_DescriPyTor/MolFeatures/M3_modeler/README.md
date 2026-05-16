# M3 Modeler — Machine Learning for Molecular Modeling

A comprehensive package for developing, evaluating, and applying machine‑learning models to molecular datasets. Supports both **regression** and **classification** with specialized tools for **feature selection**, **model evaluation**, **predictive analysis**, and dataset **balancing/sampling**.

---

## Table of Contents

* [Features](#features)
* [Installation](#installation)
* [Data Expectations](#data-expectations)
* [Quick Start](#quick-start)

  * [Command Line](#command-line)
  * [Python API](#python-api)
* [Cross-Validation Methods](#cross-validation-methods)
* [Core Classes](#core-classes)

  * [LinearRegressionModel](#linearregressionmodel)
  * [ClassificationModel](#classificationmodel)
  * [OrdinalLogisticRegression](#ordinallogisticregression)
* [Feature Selection Functions](#feature-selection-functions)
* [Utility Functions](#utility-functions)
* [Database Integration](#database-integration)
* [Similarity & Stratified Sampling](#similarity--stratified-sampling)
* [Reports & Visualizations](#reports--visualizations)
* [Recipes](#recipes)
* [Troubleshooting](#troubleshooting)
* [FAQ](#faq)
* [License](#license)

---

## Features

**Regression**

* Robust evaluation: R², Q², MAE, RMSD with repeated CV (K‑fold / LOO / repeated holdout).
* Exhaustive feature subset search (bounded by `min_features`/`max_features`).
* Prediction intervals and coefficient statistics.
* Multicollinearity checks (VIF) + diagnostics.
* Resilient to missing and non‑numeric features (coercion + safe dropping).

**Classification**

* Metrics: Accuracy, F1, McFadden’s R².
* Stratified K‑fold (class‑aware) CV; automatic adjustment for small datasets.
* Class balancing via stratified sampling & similarity‑based sampling.

**Common**

* Results persisted to SQLite for reproducibility and incremental runs.
* Compact PDF dashboards (coefficients, CV metrics, VIF, diagnostics; similarity plots before/after sampling).

---


## Data Expectations

You can work in **one‑CSV** or **two‑CSV** mode.

* **One CSV**: First column is a sample identifier (e.g., `Name`). Remaining columns contain numeric **features** and a single **target** column.
* **Two CSVs**: One file for features (rows = samples, columns = descriptors) and one for target/labels. Rows must align or include a join key.

**Tidy tips**

* Non‑numeric feature columns are coerced when possible; all‑NaN or constant columns are dropped.
* If you plan to use `--leave-out`, ensure the identifier column exists and matches your samples.

---

## Quick Start

### Command Line

```bash
python __main__.py model\
  --mode {regression|classification} \
  --features_csv path/to/features.csv \
  --target_csv path/to/targets.csv \
  --y_value EXPERIMENT_TAG \
  --n_jobs -1 \
  --min-features 1 \
  --max-features 5 \
  --top-n 20 \
  --bool-parallel \
  --threshold 0.5 \
  --leave-out SAMPLE_A SAMPLE_B
```

**Arguments**

| Flag                   | Type / Choices                   | Required | Description                                                                            |
| ---------------------- | -------------------------------- | -------: | -------------------------------------------------------------------------------------- |
| `-h`, `--help`         | —                                |       No | Show help and exit.                                                                    |
| `-m`, `--mode`         | `regression` \| `classification` |  **Yes** | Task to run; affects metrics and thresholds.                                           |
| `-f`, `--features_csv` | path                             |  **Yes** | Path to CSV with input features (descriptors).                                         |
| `-t`, `--target_csv`   | path                             |  **Yes** | Path to CSV with targets/labels.                                                       |
| `-y`, `--y_value`      | str                              |  **Yes** | Tag/slug used to name outputs (plots, DB).                                             |
| `-j`, `--n_jobs`       | int                              |       No | CPU cores. `-1` = all. Defaults to `$NSLOTS` or all cores.                             |
| `--min-features`       | int                              |       No | Minimum features per model.                                                            |
| `--max-features`       | int                              |       No | Maximum features per model.                                                            |
| `--top-n`              | int                              |       No | Number of top models to keep/evaluate.                                                 |
| `--bool-parallel`      | flag                             |       No | Enable parallel evaluation (honors `--n_jobs`).                                        |
| `--threshold`          | float                            |       No | Initial quality cutoff. **Regression**: min R². **Classification**: min McFadden’s R². |
| `--leave-out`          | list\[str]                       |       No | Space‑separated sample IDs to exclude from fitting/evaluation.                         |

**Examples**

Regression

```bash
python __main__.py.py \
  --mode regression \
  --features_csv data/x_features.csv \
  --target_csv data/y_targets.csv \
  --y_value yield \
  --min-features 1 --max-features 3 \
  --top-n 10 \
  --threshold 0.6
```

Classification

```bash
python __main__.py model \
  --mode classification \
  --features_csv data/x_features.csv \
  --target_csv data/classes.csv \
  --y_value enantio_select \
  --min-features 1 --max-features 4 \
  --top-n 20 \
  --threshold 0.5 \
  --bool-parallel --n_jobs -1 \
  --leave-out LS1714 LS2008
```

### Python API

Regression

```python
from M3_modeler.modeling import LinearRegressionModel

csvs = {
    "features_csv_filepath": "data/Train.csv",
    "target_csv_filepath": ""  # leave empty in one‑CSV mode
}
model = LinearRegressionModel(csvs, process_method="one csv", y_value="pIC50",
                              min_features_num=1, max_features_num=3)

model.search_models(top_n=10, threshold=0.6, bool_parallel=True)
# Predict on left‑out set (after leave_out_samples(...))
pred, intervals = model.predict_for_leftout(return_interval=True)
```

Classification

```python
from M3_modeler.modeling import ClassificationModel

csvs = {
    "features_csv_filepath": "data/Train.csv",
    "target_csv_filepath": ""
}
clf = ClassificationModel(csvs, process_method="one csv", output_name="class",
                          min_features_num=1, max_features_num=4)

# Optionally downsample class 1 by similarity (drop indices from class=1)
drop_pos = clf.simi_sampler_(class_label=1, compare_with=0, plot=True, sample_size=40)

clf.search_models(top_n=20, mcfadden_threshold=0.5)
```

---

## Cross-Validation Methods

This project provides several cross-validation (CV) strategies designed for honest estimates of generalization performance and to prevent data leakage. All preprocessing that can leak signal (e.g., standardization) is fit only on the training fold and applied to the validation fold.

### Shared Principles

* No leakage: scaling/transformations are fit on the train split of each fold, then applied to the validation split.
* Out-of-fold (OOF) predictions: for each sample, we collect predictions from the fold where that sample was in validation. These OOF predictions are used to compute global metrics (e.g., Q²).
* Aggregation: per-fold metrics are summarized (mean/median and optionally stdev). Global metrics like Q² are computed from OOF predictions.
* Small datasets: when the dataset is very small, the workflow can fall back to Leave-One-Out CV (LOOCV) automatically.
* Left-out set vs CV: any explicit left-out samples you provide are held out entirely and not used in CV. They can be evaluated later with `predict_for_leftout`.

### Regression CV

* **K-Fold CV**
  Split the data into K folds; iterate K times, each time training on K−1 folds and validating on the remaining fold. Collect OOF predictions to compute:

  * **Q²** = 1 − SSE\_OOF / SST, where SSE\_OOF = sum((y − yhat\_OOF)^2) over all samples and SST = sum((y − mean(y))^2).
  * **R², MAE, RMSD**: reported per fold and summarized across folds.

* **Leave-One-Out CV (LOOCV)**
  A special case of K-Fold with K = N (one sample per validation split). Provides low-bias estimates for tiny datasets but is computationally heavier.

* **Repeated Holdout / Shuffle-Split (optional)**
  Multiple random train/validation splits (e.g., 80/20) can be run and aggregated. Useful for stress-testing stability when K-Fold may be brittle on small N.

### Classification CV

* **Stratified K-Fold CV**
  Preserves class proportions in each fold. Useful for imbalanced datasets and multi-class problems. Metrics collected per fold include:

  * Accuracy and F1 (macro or weighted depending on configuration),
  * McFadden’s R² (pseudo-R² derived from per-fold log-likelihoods),
    which are then summarized across folds.

* **Edge cases / small classes**
  The number of folds may be reduced automatically to avoid empty class folds. When extremely small, LOOCV-like behavior is used.

> Tip: For reproducibility, set a fixed random seed where applicable. For highly imbalanced data, consider pairing Stratified K-Fold with the built-in similarity-based sampling or stratified sampling utilities before CV.

## Core Classes

### `LinearRegressionModel`

* Calculates **R², Q², MAE, RMSD** with robust cross‑validation (LOO / K‑fold / repeated holdout).
* Exhaustive **feature subset search** (bounded by `min_features_num`/`max_features_num`).
* **Prediction intervals** and coefficient statistics for interpretability.
* **Multicollinearity** detection via VIF; diagnostics and (optional) pruning.
* Robust to missing values and non‑numeric columns.
* Proper scaling performed **inside each fold** to avoid leakage.

### `ClassificationModel`

* Metrics: **Accuracy, F1, McFadden’s R²** (per‑fold; aggregated).
* **Stratified K‑fold** with class‑aware fold generation; automatic adjustments for small classes.
* **Class balancing** via:

  * Stratified sampling (preserve class ratios).
  * **Similarity‑based sampling** to obtain a diverse subset within a class.

### `OrdinalLogisticRegression`

A scikit‑learn‑style wrapper around `statsmodels`’ `OrderedModel` that supports:

* `fit`, `predict`, and `predict_proba` APIs.
* Different link functions (logit, probit, etc.).
* Ordered categorical targets.

---

## Feature Selection Functions

* `fit_and_evaluate_single_combination_regression(X, y, features, ...)`

  > Evaluate one feature combination for regression; returns metrics, coefficients, and diagnostics.

* `fit_and_evaluate_single_combination_classification(X, y, features, ...)`

  > Evaluate one feature combination for classification; returns metrics (Accuracy, F1, McFadden’s R²), confusion matrix, etc.

* `search_models(top_n, threshold, ...)` (method on both classes)

  > Exhaustive search over feature combinations with caching, pruning via threshold, and optional parallelism.

* `run_single_combo_report(regression_model,['x','y'])` (method on both classes)
  > runs a single model regression with added statistics plots.
---

## Utility Functions

* `check_linear_regression_assumptions(X, y)`

  > Statistical tests: Durbin–Watson, Breusch–Pagan, Shapiro–Wilk, plus Q–Q plots.

* `_sort_results(results_df)`

  > Sorts model results by performance (Q², then R² for regression; McFadden’s R² for classification).

* `_compute_vif(X)`

  > Variance Inflation Factor to flag multicollinearity.

---

## Database Integration

* `insert_result_into_db_regression(model_info, metrics, ...)`
* `insert_result_into_db_classification(model_info, metrics, ...)`
* `load_results_from_db(db_path, table)`

> Results are persisted to a SQLite DB per dataset, enabling incremental runs and reproducible top‑N selection.

---

## Similarity & Stratified Sampling

* `simi_sampler(data, class_label, compare_with=0, sample_size=K, return_kind="pos")`

  > Compute cosine‑like similarity on scaled features; **keep K** evenly spaced over the similarity range within `class_label` and **drop the rest**. Returns indices to **drop** (or names/labels if requested).

* `ClassificationModel.simi_sampler_(class_label, compare_with=0, plot=True, sample_size=None)`

  > Convenience wrapper that runs the sampler and **applies the drop** (only within the target class). Plots “before vs after” with dropped points marked.

* `stratified_sampling_with_plots(data, plot, sample_size)`

  > Create a balanced dataset with optional class distribution plots.

**Tip:** Ensure the sampler uses/reset index (`reset_index(drop=True)`) so positional indices align with your working `X`/`y`.

---

## Reports & Visualizations

Each top model can generate a compact PDF “summary dashboard” including:

* Coefficients table (with signs and magnitudes),
* VIF table and multicollinearity warnings,
* Cross‑validation metrics table (per fold and aggregated),
* For classification: confusion matrix/ROC/PR as configured,
* For sampling: similarity scatter **before vs after** (dropped points highlighted).

> The table renderer (`_safe_table`) auto‑wraps and shrinks to fit, preferring readability and preventing unreadably tiny text.

---

## Recipes

**Balance a specific class by similarity, then search**

```python
clf = ClassificationModel(csvs, process_method="one csv", output_name="class",
                          min_features_num=1, max_features_num=4)
# Drop within class=1, keep 10 most representative across similarity
clf.simi_sampler_(class_label=1, compare_with=0, sample_size=10, plot=True)
clf.search_models(top_n=30, mcfadden_threshold=0.5)
```

**Hold out a list of samples by ID, train on the rest, then predict**

```python
model.leave_out_samples(["Mol_0012", "Mol_0345"], keep_only=False)
model.search_models(top_n=10, threshold=0.65)
pred, intervals = model.predict_for_leftout(return_interval=True)
```

---

## Troubleshooting

* **“dtype='numeric' is not compatible with arrays of bytes/strings”**
  Your routine expects numeric arrays but received strings (e.g., class labels). Coerce features/targets with `pd.to_numeric(..., errors='coerce')` (regression) or skip numeric plots in classification.

* **SHAP error with KernelExplainer**
  Pass a **DataFrame** with real column names to SHAP; ensure `X` is numeric (coerce and drop NaNs) and call with `X` (not just a list of feature names).

* **Similarity sampler dropped nothing**
  If `sample_size >= #rows in class_label`, nothing is dropped. Reduce `sample_size`. Also ensure `class_label` dtype matches (`1` vs `'1'`).

* **Tables in PDF are tiny**
  Use the updated `_safe_table` (targets nearly full axes width; enlarges figure before shrinking text). Give the table a generous subplot or reduce `max_rows`/`max_cols`.

---

## FAQ

**Q: Does `simi_sampler` remove from the compared class or the target class?**
A: It removes **from `class_label`** only. `compare_with` chooses the similarity axis (self‑similarity when 0).

**Q: What does `sample_size` mean in `simi_sampler`?**
A: The number of samples to **keep** within `class_label`, evenly spread across the similarity range; the rest are dropped.

**Q: How do I apply drops only inside one class?**
A: Build a mask: keep all rows from other classes, then either keep or drop the selected positions inside `class_label` only.

---

## License

MIT License (see `LICENSE` for details).
