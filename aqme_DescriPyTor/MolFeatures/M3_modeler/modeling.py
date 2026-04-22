# -*- coding: latin-1 -*-
import cProfile
import os
if not os.environ.get("MPLBACKEND"):
    os.environ["MPLBACKEND"] = "Agg"   # headless, file-only
import pstats
import random
import sqlite3
import sys
import time
import multiprocessing
from itertools import combinations
from tkinter import filedialog, messagebox
from typing import Iterable, List, Sequence, Set, Tuple, Union, Optional

# --- Third-party ---
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import PIL
import scipy.stats as stats
from tqdm import tqdm
from joblib import Parallel, delayed
from scipy.stats import t

# scikit-learn
from sklearn.compose import ColumnTransformer
from sklearn.linear_model import (
    BayesianRidge,
    Lasso,
    LinearRegression,
    LogisticRegression,
    SGDRegressor,
)
from sklearn.utils.validation import check_is_fitted
from sklearn.exceptions import NotFittedError 
from sklearn.metrics import (
    accuracy_score,
    confusion_matrix,
    f1_score,
    make_scorer,
    mean_absolute_error,
    mean_squared_error,
    precision_score,
    r2_score,
    recall_score,
    roc_auc_score,
)
from sklearn.model_selection import (
    KFold,
    LeaveOneOut,
    RepeatedKFold,
    cross_val_predict,
    cross_validate,
    cross_val_score,
    train_test_split,
)
from sklearn.preprocessing import StandardScaler
from statsmodels.miscmodels.ordinal_model import OrderedModel
from statsmodels.stats.outliers_influence import variance_inflation_factor
import pymc as pm
import arviz as az



sys.path.append(os.path.dirname(os.path.abspath(__file__)))
try:
    from plot import *
    from modeling_utils import simi_sampler, stratified_sampling_with_plots, _normalize_combination_to_columns , _parse_tuple_string
    from modeling_utils import *
except:
    from M3_modeler.plot import *
    from M3_modeler.modeling_utils import simi_sampler, stratified_sampling_with_plots, _normalize_combination_to_columns , _parse_tuple_string
    from M3_modeler.modeling_utils import *



def _combo_key(combo: Sequence[str]) -> str:
    """Canonical string key for a combo (order-insensitive)."""
    return ",".join(sorted(map(str, combo)))

def _to_df(maybe_df) -> pd.DataFrame:
    """Coerce list/None/DataFrame to a DataFrame."""
    
    if isinstance(maybe_df, pd.DataFrame):
        return maybe_df
    if maybe_df is None:
        return pd.DataFrame()
    if isinstance(maybe_df, list):
        return pd.DataFrame(maybe_df)
    return pd.DataFrame(maybe_df)

def _ensure_numeric(s: pd.Series) -> pd.Series:
    """Coerce to numeric; invalid -> NaN."""
    return pd.to_numeric(s, errors="coerce")

def _extract_series(df: pd.DataFrame, names: Tuple[str, ...]) -> Optional[pd.Series]:
    """Return the first matching column (case-insensitive), else None."""
    if df.empty:
        return None
    lower = {c.lower(): c for c in df.columns}
    for n in names:
        if n.lower() in lower:
            return df[lower[n.lower()]]
    return None

def _extract_q2(df: pd.DataFrame) -> Optional[pd.Series]:
    """Get Q2 (or q2) as numeric; if only in 'scores' dicts, pull from there."""
    if df.empty:
        return None
    s = _extract_series(df, ("q2", "Q2"))
    if s is not None:
        return _ensure_numeric(s)
    # scores column path
    scores = _extract_series(df, ("scores",))
    if scores is not None:
        s = scores.apply(lambda d: (d or {}).get("q2", (d or {}).get("Q2", np.nan)))
        return _ensure_numeric(s)
    return None

def _extract_r2(df: pd.DataFrame) -> Optional[pd.Series]:
    s = _extract_series(df, ("r2", "R2"))
 
    return _ensure_numeric(s) if s is not None else None

def _is_all_q2_neg_inf(df: pd.DataFrame) -> bool:
    """True iff Q2 exists and every non-NaN value is exactly -inf (NaN => not -inf)."""
    q2 = _extract_q2(df)
    if q2 is None or q2.dropna().empty:
        return False
    return bool((q2.dropna() == -np.inf).all())

def _best_r2(df: pd.DataFrame) -> Optional[float]:
    r2 = _extract_r2(df)
    if r2 is None or r2.dropna().empty:
        return None
    return float(r2.max())

def _sort_results(df: pd.DataFrame) -> pd.DataFrame:
    """
    Sort model results by Q² (descending) if present, then by R².
    Falls back to R²-only sorting if Q² not found.

    Parameters
    ----------
    df : pd.DataFrame
        Results dataframe expected to contain columns like:
        ['id', 'threshold', 'model', 'predictions', 'q2', 'r2', ...]

    Returns
    -------
    pd.DataFrame
        Sorted dataframe (copy), cleaned of non-metric columns.
    """
    if df is None or df.empty:
        print("⚠️ Empty DataFrame received in _sort_results.")
        return pd.DataFrame()

    # drop non-metric columns if they exist
    drop_cols = [c for c in ['id', 'threshold', 'model', 'predictions'] if c in df.columns]
    df = df.drop(columns=drop_cols, errors='ignore').copy()

    # filter out any rows with inf or NaN values in numeric columns
    numeric_cols = df.select_dtypes(include=[float, int]).columns
    if len(numeric_cols) > 0:
        df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=numeric_cols, how="any")

    if df.empty:
        print("⚠️ No valid rows after cleaning in _sort_results.")
        return df

    # --- your existing logic follows ---
    q2 = _extract_q2(df)
    r2 = _extract_r2(df)

    temp = df.copy()
    if q2 is not None:
        temp["_q2__"] = q2
    if r2 is not None:
        temp["_r2__"] = r2

    by = []
    if "_q2__" in temp.columns:
        by.append("_q2__")
    if "_r2__" in temp.columns:
        by.append("_r2__")

    if not by:
        return df

    temp = temp.sort_values(by=by, ascending=[False] * len(by))
    tmp_re = temp.drop(columns=[c for c in ("_q2__", "_r2__") if c in temp.columns])
    
    return tmp_re 



        # Print the count of combinations for each number of features
def fit_and_evaluate_single_combination_classification(
    model,
    combination,
    threshold: float = 0.5,
    return_probabilities: bool = False,
    categorical_max_unique: int = 20,
    categorical_frac: float = 0.05,
):
    import numpy as np
    import pandas as pd

    try:
        colms = _normalize_combination_to_columns(combination)
        selected_features = model.features_df[colms]
    except Exception:
        selected_features = model.features_df[list(combination)]

    X_df = selected_features.copy()
    y = model.target_vector.to_numpy()

    # --- Separate categorical-like integer columns ---
   
    categorical_cols = []
    for c in X_df.columns:
        series = X_df[c].dropna()
        if pd.api.types.is_integer_dtype(series):
            categorical_cols.append(c)
        elif pd.api.types.is_float_dtype(series):
            # Check if all non-NaN values are integers (e.g. 0.0, 1.0, 2.0)
            if np.allclose(series, series.astype(int)):
                categorical_cols.append(c)

    continuous_cols = [c for c in X_df.columns if c not in categorical_cols]

    # --- One-hot encode categorical columns ---
    X_cont = X_df[continuous_cols]

    valid_cats = [c for c in categorical_cols if c in X_df.columns]
    # print(valid_cats,'valid_cats',categorical_cols,'categorical_cols', X_df,'X_df.')
    if valid_cats:
        X_cats = pd.get_dummies(
            X_df,
            columns=valid_cats,
            prefix={c: c for c in categorical_cols},
            prefix_sep="_",
            dummy_na=False
        )
        X_cont = X_df[continuous_cols]
        X_processed = pd.concat([X_cont, X_cats], axis=1)
    else:
        X_processed = X_df


    X = X_processed.to_numpy()

    # --- Fit the model ---
    model.fit(X, y)

    # --- Evaluate the model ---
    evaluation_results, y_pred = model.evaluate(X, y)
    mcf = evaluation_results.get("mcfadden_r2", None)
    acc = evaluation_results.get("accuracy", None)

    # Use accuracy if McFadden's R² is missing or not a finite number
    use_accuracy = (mcf is None) or (isinstance(mcf, float) and not np.isfinite(mcf))

    gate_metric = acc if use_accuracy else mcf
    metric_name = "accuracy" if use_accuracy else "mcfadden_r2"

    if gate_metric is not None and gate_metric > threshold:
        avg_accuracy, avg_f1, avg_r2 = model.cross_validation(X, y, model.n_splits)
        evaluation_results["avg_accuracy"] = avg_accuracy
        evaluation_results["avg_f1_score"] = avg_f1
        evaluation_results["avg_mcfadden_r2"] = avg_r2



    results = {
        "combination": combination,
        "scores": evaluation_results,
        "model": model,
        "predictions": y_pred,
    }

    # --- Optionally return probabilities ---
    if return_probabilities:
        if getattr(model, "ordinal", False):
            cumulative_probs = model.result.predict(exog=X)
            class_probs = np.zeros_like(cumulative_probs)

            class_probs[:, 0] = cumulative_probs[:, 0]
            for i in range(1, cumulative_probs.shape[1]):
                class_probs[:, i] = cumulative_probs[:, i] - cumulative_probs[:, i - 1]
            class_probs[:, -1] = 1 - cumulative_probs[:, -2]

            prob_df = pd.DataFrame(
                class_probs,
                columns=[f"Prob_Class_{i+1}" for i in range(class_probs.shape[1])]
            )
            prob_df["Predicted_Class"] = model.model.predict(X)
            prob_df["True_Class"] = y
            return results, prob_df
        else:
            probabilities = model.model.predict_proba(X)
            prob_df = pd.DataFrame(
                probabilities,
                columns=[f"Prob_Class_{i+1}" for i in range(probabilities.shape[1])]
            )
            prob_df["Predicted_Class"] = model.model.predict(X)
            prob_df["True_Class"] = y
            return results, prob_df

    return results


def fit_and_evaluate_single_combination_regression(model, combination, r2_threshold=0.7):
    try:
        colms = _normalize_combination_to_columns(combination)
        selected_features = model.features_df[colms]
    except:
        selected_features = model.features_df[_parse_tuple_string(combination)] 

    X = selected_features.to_numpy()
    y = model.target_vector.to_numpy()
  
    # Fit the model
    t0 = time.time()
    model.fit(X, y)
    ## check
    model.model.fit(X, y)
    fit_time=time.time()-t0
    # Evaluate the modeld
    t1=time.time()
    model._trained_features = list(combination)
    evaluation_results, y_pred = model.evaluate(X, y)
    eval_time=time.time()-t1
    coefficients,intercepts = model.get_coefficients_from_trained_estimator()
 
    # Check if R-squared is above the threshold
    t3=time.time()
    if evaluation_results['r2'] > r2_threshold:
        q2, mae,rmsd = model.calculate_q2_and_mae(X, y, n_splits=1)
        evaluation_results['Q2'] = q2
        evaluation_results['MAE'] = mae
        evaluation_results['RMSD'] = rmsd
        # print(f'R2:{evaluation_results["r2"]:.3f} Q2: {q2:.3f}, MAE: {mae:.3f}, RMSD: {rmsd:.3f} for combination: {combination}')

    q2_time=time.time()-t3

    result = {
        'combination': combination,
        'scores': evaluation_results,
        'intercept': intercepts,
        'coefficients': coefficients,
        'model': model,
        'predictions': y_pred
    }

    if  'scores' in result:
        try:
            r2 = result['scores'].get('r2', float('-inf'))
            q2 = result['scores'].get('Q2', float('-inf'))
            mae = result['scores'].get('MAE', float('inf'))
            rmsd = result['scores'].get('RMSD', float('inf'))
            csv_path=model.db_path.replace('.db','.csv')
            model = result['model']
            # Insert into DB
            result_dict = {
                'combination': str(combination),
                'r2': r2,
                'q2': q2,
                'mae': mae,
                'rmsd': rmsd,
                'threshold': r2_threshold,
                'model': model,
                'predictions': y_pred,
            }
            # print(type(model),'type of model variable')
            insert_result_into_db_regression(
                db_path=model.db_path,
                combination=combination,
                r2=r2,
                q2=q2,
                mae=mae,
                rmsd=rmsd,
                threshold=r2_threshold,
                csv_path=csv_path,
                model=model,
                predictions=y_pred
            )
            return result_dict
        except Exception as e:
            print(f'failed to insert combination : {combination}')
            print(e)

    return result




class PlotModel:
    def __init__(self, model):
        self.model = model
    
    


class LinearRegressionModel:

    def __init__(
            self, 
            csv_filepaths, 
            process_method='one csv', 
            names_column=None,
            y_value='output', 
            leave_out=None, 
            min_features_num=2, 
            max_features_num=None, 
            n_splits=5, 
            metrics=None, 
            return_coefficients=False, 
            model_type='linear',    # <--- Choose 'linear' or 'lasso'
            alpha=1.0,
            app=None,
            db_path='results', 
            scale=True ,
            seed = 42 ,           # <--- If lasso, this is the regularization strength
            drop_columns = 'id',
    ):

        self.random_state = seed
        self.setup_random_seeds()
        self.drop_columns=drop_columns
        self.csv_filepaths = csv_filepaths
        self.process_method = process_method
        self.y_value = y_value
        self.names_column = names_column
        self.leave_out = leave_out
        self.min_features_num = min_features_num
        self.max_features_num = max_features_num
        self.metrics = metrics if metrics is not None else ['r2', 'neg_mean_absolute_error']
        self.return_coefficients = return_coefficients
        self.model_type = model_type
        self.alpha = alpha
        self.app = app
        self.n_splits = n_splits
        self.scale=scale
        try:
            self.name=os.path.splitext(os.path.basename(self.csv_filepaths.get('features_csv_filepath')))[0]
        except:
            self.name='in_memory_dataset'

        self.paths = prepare_run_dirs(
            base_dir="runs",
            dataset_name=self.name,
            y_value=self.y_value,
            tag=self.model_type  # e.g., 'linear' or 'lasso'
        )

        self.db_path = resolve_db_path(db_path, self.name, self.paths.db)
        create_results_table(self.db_path)
        
        if self.model_type.lower() == 'linear':
            self.model = LinearRegression()
            print('linear model selected')
        elif model_type.lower() == 'lasso':
            self.model = Lasso(alpha=alpha)
            print('lasso model selected')
        else:
            raise ValueError("Invalid model_type. Please specify 'linear' or 'lasso'.")

        try:
            if csv_filepaths:
                if process_method == 'one csv':
                    self.process_features_csv(drop_columns=drop_columns)
                elif process_method == 'two csvs':
                    self.process_features_csv(drop_columns=drop_columns)
                    self.process_target_csv(csv_filepaths.get('target_csv_filepath'))
        except Exception as e:
            self.process_features_csv(drop_columns=drop_columns)
            self.features_df = self.features_df.dropna(axis=1)

        self.scaler = StandardScaler()
        self.original_features_df = self.features_df.copy()
        self.features_df =  pd.DataFrame(self.scaler.fit_transform(self.features_df), columns=self.features_df.columns)
        self.feature_names=self.features_df.columns.tolist()
        
        # self.compute_correlation()
       
        self.leave_out_samples(self.leave_out)
        self.determine_number_of_features()

        

        self.theta       = None
        self.X_b_train   = None
        self.sigma2      = None
        self.XtX_inv     = None
        
        ## need to be preformed on the chosen model

        # self.check_linear_regression_assumptions()
    

    from pathlib import Path


    def load_results(self, sort_key: str = "r2", ascending: bool = False):
        """
        Load all CSV results from the run's db directory into DataFrames,
        optionally sorting them by a metric (default: R²).

        Parameters
        ----------
        sort_key : str
            Column name to sort by (default 'r2').
        ascending : bool
            Sort order. Default False = best (highest values first).

        Returns
        -------
        results_dict : dict[str, pd.DataFrame]
            Mapping {filename (no extension): DataFrame}.
        """
        results_dict = {}

        # Use Path.glob so csv_file is a Path, not a str
        for csv_file in Path(self.paths.db).glob("*.csv"):
            try:
                df = pd.read_csv(csv_file)

                if sort_key in df.columns:
                    df = df.sort_values(by=sort_key, ascending=ascending).reset_index(drop=True)
                else:
                    print(f"⚠️ Column '{sort_key}' not found in {csv_file.name}, returning unsorted.")

                results_dict[csv_file.stem] = df
            except Exception as e:
                print(f"❌ Failed to load {csv_file.name}: {e}")
        # remove threshold, model, predictions columns
        for df in results_dict.values():
            df.drop(columns=['threshold', 'model', 'predictions'], errors='ignore', inplace=True)
        return results_dict



    def setup_random_seeds(self) -> None:
        """
        Set random seeds across Python, NumPy, and optionally PyTorch 
        to ensure reproducible results. Also stores a `self.random_state`
        attribute for use in other methods.
        """
        # Python built-in
        random.seed(self.random_state)

        # NumPy
        np.random.seed(self.random_state)

        # Make Python hash-based operations deterministic
        os.environ["PYTHONHASHSEED"] = str(self.random_state)

    

    def check_linear_regression_assumptions(self):
        """
        Check linear regression assumptions (linearity, homoscedasticity, normality).
        """
        # check if indices align , if not change features_df index to target_vector index
        if not self.features_df.index.equals(self.target_vector.index):
            self.features_df.index = self.target_vector.index
        return check_linear_regression_assumptions(self.features_df, self.target_vector)

    
    def compute_multicollinearity(self, vif_threshold=5.0):
        """
        Compute the Variance Inflation Factor (VIF) for each feature in the dataset.
        """
        # Compute VIF
        vif_results = self._compute_vif(self.features_df)
        # Identify features with high VIF
        high_vif_features = vif_results[vif_results['VIF'] > vif_threshold]
        
        return vif_results
    

    def process_features_csv(
    self,
    drop_columns: Optional[Union[str, List[str]]] = None,
    drop_smiles: bool = True,
) -> dict:
        """
        Load and process a features dataset from:
        - a CSV file path (string),
        - a dict with key 'features_csv_filepath',
        - a single-path list/tuple,
        - or an in-memory pandas DataFrame.

        Initializes:
            - self.molecule_names : list[str]
            - self.target_vector  : pd.Series (indexed by molecule names)
            - self.features_df    : pd.DataFrame (numeric-only, deduplicated, indexed by molecule names)
            - self.features_list  : list[str]
    """
        import os
        from typing import List, Optional, Union
        import numpy as np
        import pandas as pd

        # ---------- Load source ----------
        csv_source = getattr(self, "csv_filepaths", None)
        if csv_source is None:
            raise AttributeError(
                "`self.csv_filepaths` must be set to a CSV path (str), a dict, a pandas DataFrame, "
                "or a single-path list/tuple."
            )

        if isinstance(csv_source, (list, tuple)):
            if len(csv_source) == 1 and isinstance(csv_source[0], str):
                csv_source = csv_source[0]
            else:
                raise TypeError(
                    "`csv_filepaths` list/tuple must contain exactly one string path; "
                    f"got {csv_source!r}"
                )

        names_from_index = False
        path_used = None

        if isinstance(csv_source, pd.DataFrame):
            df = csv_source.copy()
        elif isinstance(csv_source, dict):
            csv_path = csv_source.get("features_csv_filepath")
            if not csv_path or not os.path.exists(csv_path):
                raise FileNotFoundError(f"CSV not found at: {csv_path!r}")
            df = pd.read_csv(csv_path)
            path_used = os.path.abspath(csv_path)
        elif isinstance(csv_source, str):
            if not os.path.exists(csv_source):
                raise FileNotFoundError(f"CSV not found at: {csv_source!r}")
            df = pd.read_csv(csv_source)
            path_used = os.path.abspath(csv_source)
        else:
            raise TypeError("`self.csv_filepaths` must be a str, dict, DataFrame, or single-path list/tuple.")

        if df.empty:
            raise ValueError("The provided dataset is empty.")

        # ---------- Helpers ----------
        def _case_insensitive_find(col_name: str, columns: pd.Index):
            lc = str(col_name).lower()
            for c in columns:
                if str(c).lower() == lc:
                    return c
            return None

        def _resolve_many(names: List[str], columns: pd.Index) -> List[str]:
            resolved = []
            for n in names:
                hit = _case_insensitive_find(n, columns)
                if hit is not None:
                    resolved.append(hit)
            return list(dict.fromkeys(resolved))  # dedupe

        def _first_match(candidates: List[str], columns: pd.Index):
            for c in candidates:
                hit = _case_insensitive_find(c, columns)
                if hit is not None:
                    return hit
            return None

        def _clean_and_make_unique(names: List[str]) -> List[str]:
            # strip, normalize 'nan', and make unique with suffixes: name, name#2, name#3, ...
            base = []
            for x in names:
                s = "" if x is None else str(x)
                s = s.strip()
                if s.lower() == "nan":
                    s = ""
                base.append(s)

            result = []
            seen = {}
            for s in base:
                key = s if s else ""  # allow empty then disambiguate
                if key not in seen:
                    seen[key] = 1
                    result.append(key if key else "row")
                else:
                    seen[key] += 1
                    # append a simple counter
                    tag = f"{key if key else 'row'}#{seen[key]}"
                    result.append(tag)
            return result

        # ---------- Target & name resolution ----------
        y_value = getattr(self, "y_value", None)
        if y_value is None:
            raise AttributeError("`self.y_value` must be set to the target column name.")

        if drop_columns is None and hasattr(self, "drop_columns"):
            drop_columns = getattr(self, "drop_columns")

        if isinstance(drop_columns, str):
            drop_columns = [drop_columns]
        elif drop_columns is None:
            drop_columns = []

        provided_names_col = getattr(self, "names_column", None)
        name_candidates = ["name", "names", "molecule", "molecule_name", "mol", "compound", "sample", "id"]

        # Choose names column (honor provided_names_col first)
        names_col = None
        if provided_names_col:
            names_col = _case_insensitive_find(provided_names_col, df.columns)
            if names_col is None:
                raise KeyError(
                    f"names_column {provided_names_col!r} not found. Available: {list(df.columns)}"
                )
        else:
            if not isinstance(df.index, pd.RangeIndex):
                names_col = None  # use index
                names_from_index = True
            else:
                names_col = _first_match(name_candidates, df.columns) or df.columns[0]

        target_col = _case_insensitive_find(y_value, df.columns)
        if target_col is None:
            raise KeyError(f"Target column {y_value!r} not found. Available: {list(df.columns)}")

        # ---------- Build features frame ----------
        cols_to_drop = {c for c in [names_col, target_col] if c is not None}
        features_df_raw = df.drop(columns=[c for c in cols_to_drop if c in df.columns])
        n_features_total = features_df_raw.shape[1]

        dropped_smiles_like_cols, dropped_user_columns = [], []

        # Drop SMILES-like columns
        if drop_smiles:
            smiles_candidates = ["smiles", "canonical_smiles", "random_smiles"]
            smiles_to_drop = _resolve_many(smiles_candidates, features_df_raw.columns)
            if smiles_to_drop:
                features_df_raw = features_df_raw.drop(columns=smiles_to_drop)
                dropped_smiles_like_cols.extend(smiles_to_drop)

        # Drop user-specified columns (support names or integer positions)
        if drop_columns:
            dc = list(drop_columns)
            by_pos = [c for c in dc if isinstance(c, (int, np.integer))]
            by_name = [c for c in dc if not isinstance(c, (int, np.integer))]

            resolved_user_cols = _resolve_many(by_name, features_df_raw.columns)

            if by_pos:
                pos_names = [
                    features_df_raw.columns[int(i)] for i in by_pos
                    if 0 <= int(i) < features_df_raw.shape[1]
                ]
                resolved_user_cols.extend(pos_names)

            resolved_user_cols = list(dict.fromkeys(resolved_user_cols))  # dedupe
            if resolved_user_cols:
                features_df_raw = features_df_raw.drop(columns=resolved_user_cols)
                dropped_user_columns.extend(resolved_user_cols)

        # Deduplicate columns
        features_df_raw = features_df_raw.loc[:, ~features_df_raw.columns.duplicated()]

        # ---------- Coerce to numeric ----------
        coerced = features_df_raw.apply(pd.to_numeric, errors="coerce")
        all_nan_cols = [c for c in coerced.columns if coerced[c].isna().all()]
        numeric_df = coerced.drop(columns=all_nan_cols, errors="ignore").select_dtypes(include=[np.number])
        dropped_non_numeric = [
            c for c in features_df_raw.columns if c not in numeric_df.columns and c not in all_nan_cols
        ]

        # ---------- Names / index handling (MAIN FIX) ----------
        if names_col is None:
            names_series = pd.Series(df.index, copy=False)
            names_label = "<index>"
        else:
            names_series = pd.Series(df[names_col], copy=False)
            names_label = names_col

        clean_unique_names = _clean_and_make_unique(names_series.tolist())

        # Apply the same name index to both X and y
        numeric_df = numeric_df.copy()
        numeric_df.index = clean_unique_names

        target_series_raw = pd.to_numeric(df[target_col], errors="coerce")
        # If everything became NaN, fallback to original dtype
        target_series = target_series_raw if not target_series_raw.isna().all() else df[target_col]
        target_series = pd.Series(target_series.values, index=clean_unique_names, name=target_col)

        # Optionally drop rows with NaN targets (keeps model flows sane)
        nan_target_mask = target_series.isna()
        n_nan_targets = int(nan_target_mask.sum())
        if n_nan_targets:
            target_series = target_series[~nan_target_mask]
            numeric_df = numeric_df.loc[target_series.index]

        # Final align & checks
        shared = numeric_df.index.intersection(target_series.index)
        numeric_df = numeric_df.loc[shared]
        target_series = target_series.loc[shared]

        # Sync public attributes
        self.molecule_names = list(numeric_df.index)
        self.target_vector = target_series
        self.features_df = numeric_df
        self.features_list = numeric_df.columns.tolist()

        # ---------- Summary ----------
        name_src_note = "from index" if names_col is None else f"from column '{names_col}'"
        summary = {
            "path": path_used,
            "names_column": names_label,
            "target_column": target_col,
            "n_rows": df.shape[0],
            "n_features_total": n_features_total,
            "n_features_kept_numeric": numeric_df.shape[1],
            "dropped_non_numeric_cols": dropped_non_numeric,
            "dropped_all_nan_cols": all_nan_cols,
            "dropped_smiles_like_cols": dropped_smiles_like_cols,
            "dropped_user_columns": dropped_user_columns,
            "dropped_rows_nan_target": n_nan_targets,
            "n_final_rows": numeric_df.shape[0],
        }

        # ---------- Output/log ----------
        msg = (
            f"Processed {'DataFrame' if path_used is None else os.path.basename(path_used)}\n"
            f"Molecule names: {name_src_note} | Target: {target_col}\n"
            f"Rows: {df.shape[0]} → {summary['n_final_rows']} (dropped NaN targets: {n_nan_targets}) | "
            f"Features: {n_features_total} → {numeric_df.shape[1]}"
        )
        for key, lst in {
            "SMILES-like": dropped_smiles_like_cols,
            "User-specified": dropped_user_columns,
            "Non-numeric": dropped_non_numeric,
            "All-NaN": all_nan_cols,
        }.items():
            if lst:
                msg += f"\nDropped {key} cols: {lst}"

        app = getattr(self, "app", None)
        if app:
            try:
                app.show_result(msg)
            except Exception:
                print(msg)
        else:
            print(msg)

        # Final sanity
        assert self.features_df.index.equals(self.target_vector.index), \
            "X and y must share identical indices after processing."

        return summary



    def run_spike_and_slab_selection(self, slab_sd: float = 5.0, n_samples: int = 2000, n_tune: int = 1000,
                                     target_accept: float = 0.95, p_include: float = 0.5, verbose: bool = True):
        """
        Run spike-and-slab Bayesian variable selection using PyMC.
        Stores inclusion probabilities and posterior in self.spike_and_slab_result.
        """
        X = self.features_df.values
        y = self.target_vector.values if hasattr(self.target_vector, 'values') else self.target_vector

        p = X.shape[1]
        features = self.features_df.columns

        with pm.Model() as model:
            # Slab: wide Gaussian for nonzero
            spike = pm.Bernoulli('spike', p=p_include, shape=p)
            betas_slab = pm.Normal('betas_slab', mu=0, sigma=slab_sd, shape=p)
            betas = pm.Deterministic('betas', betas_slab * spike)
            intercept = pm.Normal('intercept', mu=0, sigma=10)
            sigma = pm.HalfNormal('sigma', sigma=1)

            mu = intercept + pm.math.dot(X, betas)
            y_obs = pm.Normal('y_obs', mu=mu, sigma=sigma, observed=y)
            trace = pm.sample(n_samples, tune=n_tune, target_accept=target_accept, chains=4, cores=1, random_seed=self.random_state, progressbar=verbose)

        inclusion_probs = trace.posterior['spike'].mean(dim=("chain", "draw")).values

        self.spike_and_slab_result = {
            'trace': trace,
            'inclusion_probs': inclusion_probs,
            'features': features,
        }

        if verbose:
            print("Spike-and-Slab Feature Inclusion Probabilities:")
            for fname, prob in zip(features, inclusion_probs):
                print(f"{fname:20s}  P(included) = {prob:.3f}")

    def get_selected_features_spike_and_slab(self, threshold: float = 0.5):
        """
        Return features whose spike-and-slab inclusion probability exceeds the threshold.
        """
        if not hasattr(self, 'spike_and_slab_result'):
            raise ValueError("Run run_spike_and_slab_selection() first.")
        inclusion_probs = self.spike_and_slab_result['inclusion_probs']
        features = self.spike_and_slab_result['features']
        selected = [f for f, prob in zip(features, inclusion_probs) if prob >= threshold]
        return selected
    
    


    def compute_correlation(
        self,
        correlation_threshold: float = 0.80,
        vif_threshold: float = None,
    ) -> dict:
        """
        Analyze multicollinearity via correlation and (optionally) VIF, with optional GUI prompts.

        - Computes the full correlation matrix on `self.features_df`.
        - Finds features involved in any pair with |r| > correlation_threshold.
        - (Optional) Visualizes a heatmap of the correlated subset.
        - (Optional) Prunes multicollinearity using VIF (greedy: drop the highest-VIF
        feature until all VIF <= vif_threshold).
        - Updates `self.features_df` and `self.features_list` in-place if pruning occurs.

        Parameters
        ----------
        correlation_threshold : float, default=0.80
            Absolute correlation cutoff to flag features.
        vif_threshold : float, default=5.0
            If not None, apply greedy VIF-based pruning to reduce multicollinearity.
            Set to None to skip VIF pruning.

        Returns
        -------
        dict
            {
            "corr_matrix": pd.DataFrame,
            "high_corr_features": set[str],
            "sub_corr": pd.DataFrame | None,
            "vif_table_before": pd.DataFrame | None,
            "vif_table_after": pd.DataFrame | None,
            "dropped_features": list[str],
            "kept_features": list[str],
            "message": str
            }
        """
        import numpy as np
        import pandas as pd
        import seaborn as sns
        import matplotlib.pyplot as plt
        from tkinter import messagebox
        from statsmodels.stats.outliers_influence import variance_inflation_factor

        def _numeric_df(df: pd.DataFrame) -> pd.DataFrame:
            return df.select_dtypes(include=[np.number]).copy()

        def _find_high_corr_features(corr: pd.DataFrame, thr: float) -> set:
            feats = set()
            cols = corr.columns
            for i in range(len(cols)):
                for j in range(i + 1, len(cols)):
                    c = corr.iloc[i, j]
                    if abs(c) > thr:
                        feats.add(cols[i])
                        feats.add(cols[j])
            return feats

        def _compute_vif_table(df_numeric: pd.DataFrame) -> pd.DataFrame:
            if df_numeric.shape[1] < 2:
                return pd.DataFrame({"feature": df_numeric.columns, "VIF": [np.nan] * df_numeric.shape[1]})
            X = df_numeric.values
            vif_vals = []
            for i in range(df_numeric.shape[1]):
                vif_vals.append(variance_inflation_factor(X, i))
            return pd.DataFrame({"feature": df_numeric.columns, "VIF": vif_vals}).sort_values("VIF", ascending=False)

        app = getattr(self, "app", None)

        # 1) Correlation matrix (numeric only)
        num_df = _numeric_df(self.features_df)
        self.corr_matrix = num_df.corr()

        # 2) High-correlation set (use helper if present, else fallback)
        try:
            high_corr_features = self._get_highly_correlated_features(self.corr_matrix, threshold=correlation_threshold)
            # Ensure it's a set of strings (column names)
            high_corr_features = set(list(high_corr_features))
        except AttributeError:
            high_corr_features = _find_high_corr_features(self.corr_matrix, correlation_threshold)

        msgs = []
        dropped_features: list[str] = []

        if high_corr_features:
            msgs.append(
                f"--- Correlation Report ---\n"
                f"Features with |r| > {correlation_threshold}:\n{sorted(high_corr_features)}"
            )

            # 3) Visualization decision
            visualize = False  # default if no app
            if app is not None:
                visualize = messagebox.askyesno(
                    title="Visualize Correlated Features?",
                    message="Show a heatmap of correlations among the flagged features?"
                )

            sub_corr = None
            if visualize and len(high_corr_features) >= 2:

                sub_corr = self.corr_matrix.loc[sorted(high_corr_features), sorted(high_corr_features)]
                interactive_corr_heatmap(sub_corr, title=f"Correlation Heatmap (|r| > {correlation_threshold})")
                # interactive_corr_heatmap_with_highlights(self.features_df)
    
            else:
                sub_corr = None

            # 4) VIF-based pruning (greedy) if requested
            vif_table_before = None
            vif_table_after = None
            apply_vif = vif_threshold is not None

            if app is not None:
                # Ask the user if they want VIF pruning
                if apply_vif:
                    apply_vif = messagebox.askyesno(
                        title="VIF-based Pruning?",
                        message=f"Apply VIF pruning to reduce multicollinearity (threshold={vif_threshold})?"
                    )

            if apply_vif:
                work_df = _numeric_df(self.features_df)  # numeric subset only
                # Start with the union of high-corr features (operating set).
                # If union has < 2, it's pointless to compute VIF—skip.
                vif_operating_cols = [c for c in work_df.columns if c in high_corr_features]
                if len(vif_operating_cols) < 2:
                    msgs.append("VIF pruning skipped (fewer than 2 correlated numeric features).")
                else:
                    vif_table_before = _compute_vif_table(work_df[vif_operating_cols])
                    msgs.append(f"Initial VIF (top 10):\n{vif_table_before.head(10).to_string(index=False)}")

                    # Greedy drop the highest-VIF feature until all VIF <= threshold
                    keep_cols = vif_operating_cols[:]
                    while True:
                        vif_tab = _compute_vif_table(work_df[keep_cols])
                        max_vif = float(vif_tab["VIF"].max())
                        if np.isnan(max_vif) or max_vif <= vif_threshold or len(keep_cols) <= 1:
                            vif_table_after = vif_tab
                            break
                        # Drop the feature with the highest VIF
                        to_drop = vif_tab.iloc[0]["feature"]
                        keep_cols.remove(to_drop)
                        dropped_features.append(to_drop)

            else:
                msgs.append("VIF pruning skipped by user/preference.")
                vif_table_before = None
                vif_table_after = None

            # Report via app if available
            final_msg = "\n".join(msgs)
            if app is not None:
                app.show_result(final_msg)
            else:
                print(final_msg)

            return {
                "corr_matrix": self.corr_matrix,
                "high_corr_features": high_corr_features,
                "sub_corr": sub_corr,
                "vif_table_before": vif_table_before,
                "vif_table_after": vif_table_after,
                "dropped_features": dropped_features,
                "kept_features": getattr(self, "features_list", list(self.features_df.columns)),
                "message": final_msg,
            }

        # No high correlations found
        msg = "No feature pairs exceeded the correlation threshold."
        if app is not None:
            app.show_result(msg)
        else:
            print(msg)

        return {
            "corr_matrix": self.corr_matrix,
            "high_corr_features": set(),
            "sub_corr": None,
            "vif_table_before": None,
            "vif_table_after": None,
            "dropped_features": [],
            "kept_features": getattr(self, "features_list", list(self.features_df.columns)),
            "message": msg,
        }



    def _compute_vif(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Compute the Variance Inflation Factor for each numeric column in df.
        """
        # 1) Select only numeric columns
        df_num = df.select_dtypes(include=[np.number]).copy()

        # 2) Handle missing or infinite values
        df_num.replace([np.inf, -np.inf], np.nan, inplace=True)
        df_num.dropna(axis=1, how='any', inplace=True)
        df_num.dropna(axis=0, how='any', inplace=True)  # ensure no NaNs in rows
       
        # 3) Check matrix rank for multicollinearity warnings
        rank = np.linalg.matrix_rank(df_num.values)
        if rank < df_num.shape[1]:
            print(f"Warning: Linear dependence detected. Matrix rank: {rank} < {df_num.shape[1]}")

        # 4) Standardize data before VIF calculation
        scaler = StandardScaler()
        df_scaled = pd.DataFrame(
            scaler.fit_transform(df_num),
            columns=df_num.columns,
            index=df_num.index
        )

        # 5) Compute VIF for each variable
        vif_data = {
            "variable": [],
            "VIF": []
        }
        for i, col in enumerate(df_scaled.columns):
            vif_val = variance_inflation_factor(df_scaled.values, i)
            vif_data["variable"].append(col)
            vif_data["VIF"].append(vif_val)

        vif_df = pd.DataFrame(vif_data)
        return vif_df


    
    def _get_highly_correlated_features(self, corr_matrix, threshold=0.8):
        """
        Identify any features whose pairwise correlation is above the threshold.
        Returns a set of the implicated feature names.
        """
        corr_matrix_abs = corr_matrix.abs()
        columns = corr_matrix_abs.columns
        high_corr_features = set()

        # We only need to look at upper triangular part to avoid duplication
        for i in range(len(columns)):
            for j in range(i + 1, len(columns)):
                if corr_matrix_abs.iloc[i, j] > threshold:
                    # Add both columns in the pair
                    high_corr_features.add(columns[i])
                    high_corr_features.add(columns[j])
        
        return high_corr_features


    def process_target_csv(self, csv_filepath):
        target_vector_unordered = pd.read_csv(csv_filepath)[self.y_value]
        
        self.target_vector = target_vector_unordered.loc[self.molecule_names]
        
    def leave_out_samples(self, leave_out=None, keep_only=False):
        """
        Remove or retain specific rows from features_df/target_vector based on indices or molecule names.
        Stores the removed (or retained) samples into:
            - self.leftout_samples
            - self.leftout_target_vector
            - self.molecule_names_predict

        Parameters:
            leave_out (int, str, list[int|str]): Indices or molecule names to leave out (or to keep if keep_only=True).
            keep_only (bool): If True, only the given indices/names are kept, all others are left out.
        """
       
        if leave_out is None:
            return

        # Normalize input to a list
        if isinstance(leave_out, (int, np.integer, str)):
            leave_out = [leave_out]
        else:
            leave_out = list(leave_out)

        # Determine indices from mixed inputs (names or indices)
        indices = []
        name_to_index = {name: idx for idx, name in enumerate(self.molecule_names)}

        for item in leave_out:
            if isinstance(item, (int, np.integer)):
                indices.append(int(item))

            elif isinstance(item, str):
                # try exact match first
                if item in name_to_index:
                    indices.append(name_to_index[item])
                else:
                    # try to parse as int
                    try:
                        idx_val = int(item)
                        indices.append(idx_val)
                    except ValueError:
                        raise ValueError(
                            f"Leave-out entry '{item}' is not a known molecule name "
                            f"and cannot be parsed as an integer index."
                        )

            else:
                raise ValueError("Items in leave_out must be integers or strings.")

        selected_features = self.features_df.iloc[indices].copy()
        selected_target = self.target_vector.iloc[indices].copy()
        selected_names = [self.molecule_names[i] for i in indices]

        # Ensure we have the full set of columns in the right order
        if hasattr(self, 'feature_names'):
            selected_features = selected_features.reindex(
                columns=self.feature_names,
                fill_value=0
            )

        # locate index labels for dropping
        feat_labels = self.features_df.index.take(indices)

        # --- decide how to address the target_vector (labels vs positions) ---
        use_label_drop_for_target = False
        targ_labels = None
        if self.target_vector.index.equals(self.features_df.index):
            # same index as features → safe to drop by the same labels
            targ_labels = feat_labels
            use_label_drop_for_target = True
        elif hasattr(self, "molecule_names") and len(self.molecule_names) == len(self.target_vector):
            # try mapping via molecule_names if they mirror the target order
            cand = pd.Index(self.molecule_names).take(indices)
            if pd.Index(cand).isin(self.target_vector.index).all():
                targ_labels = cand
                use_label_drop_for_target = True
        # else: we will drop target by positions

        # cache for downstream / debugging
        self.idx_labels = feat_labels
        self.indices_predict = indices

        if keep_only:
            # Everything NOT in indices becomes the "left out" set
            self.leftout_samples = self.features_df.drop(index=feat_labels).copy()
            if use_label_drop_for_target:
                self.leftout_target_vector = self.target_vector.drop(index=targ_labels).copy()
            else:
                self.leftout_target_vector = self.target_vector.drop(self.target_vector.index[indices]).copy()
            self.molecule_names_predict = [n for i, n in enumerate(self.molecule_names) if i not in indices]

            # Keep only the selected rows as the main set (already sliced earlier)
            self.features_df   = selected_features
            self.target_vector = selected_target
            self.molecule_names = selected_names
        else:
            # The selected rows become the "left out" set
            self.leftout_samples        = selected_features
            self.leftout_target_vector  = selected_target
            self.molecule_names_predict = selected_names

            # Drop them from the main set
            self.features_df = self.features_df.drop(index=feat_labels)
            if use_label_drop_for_target:
                self.target_vector = self.target_vector.drop(index=targ_labels)
            else:
                self.target_vector = self.target_vector.drop(self.target_vector.index[indices])
            self.molecule_names = [n for i, n in enumerate(self.molecule_names) if i not in indices]








        

    def determine_number_of_features(self):
        total_features_num = self.features_df.shape[0]
        self.max_features_num = set_max_features_limit(total_features_num, self.max_features_num)
  


    def get_feature_combinations(self):
       
        self.features_combinations = list(get_feature_combinations(self.features_list, self.min_features_num, self.max_features_num))
   

    def _add_leftout_mae_to_results(self, results: pd.DataFrame) -> pd.DataFrame:
        """
        Compute MAE of left-out predictions for each model in results.
        Fits a fresh model on training data for each combination and evaluates on left-out samples.
        
        Parameters
        ----------
        results : pd.DataFrame
            Results dataframe with 'combination' column.
            
        Returns
        -------
        pd.DataFrame
            Updated results with 'leftout_mae' column added.
        """
        if results.empty or not hasattr(self, 'leftout_samples') or self.leftout_samples is None:
            return results
        
        results = results.copy()
        leftout_mae_values = []
        
        for idx, row in results.iterrows():
            try:
                combination = row.get('combination')
                
                if combination is None:
                    leftout_mae_values.append(np.nan)
                    continue
                
                # Extract features for both training and left-out samples
                try:
                    colms = _normalize_combination_to_columns(combination)
                    train_X = self.features_df[colms]
                    leftout_X = self.leftout_samples[colms]
                except:
                    colms = _parse_tuple_string(combination)
                    train_X = self.features_df[colms]
                    leftout_X = self.leftout_samples[colms]
                
                # Get target values (ensure they're numpy arrays)
                train_y = self.target_vector.values if hasattr(self.target_vector, 'values') else np.asarray(self.target_vector)
                leftout_y = self.leftout_target_vector.values if hasattr(self.leftout_target_vector, 'values') else np.asarray(self.leftout_target_vector)
                
                # Ensure X data is numpy arrays with correct shape
                train_X_arr = train_X.values if hasattr(train_X, 'values') else np.asarray(train_X)
                leftout_X_arr = leftout_X.values if hasattr(leftout_X, 'values') else np.asarray(leftout_X)
                
                # Validation: check that shapes align
                if train_X_arr.shape[0] != train_y.shape[0]:
                    raise ValueError(
                        f"Training data mismatch: X has {train_X_arr.shape[0]} samples, "
                        f"y has {train_y.shape[0]} samples"
                    )
                if leftout_X_arr.shape[0] != leftout_y.shape[0]:
                    raise ValueError(
                        f"Left-out data mismatch: X has {leftout_X_arr.shape[0]} samples, "
                        f"y has {leftout_y.shape[0]} samples"
                    )
                
                # Create and fit a fresh model for this combination
                if self.model_type.lower() == 'linear':
                    temp_model = LinearRegression()
                elif self.model_type.lower() == 'lasso':
                    temp_model = Lasso(alpha=self.alpha)
                else:
                    raise ValueError(f"Unknown model_type: {self.model_type}")
                
                # Fit on training data
                temp_model.fit(train_X_arr, train_y)
                
                # Predict on left-out samples
                leftout_predictions = temp_model.predict(leftout_X_arr)
                
                # Calculate MAE
                mae = mean_absolute_error(leftout_y, leftout_predictions)
                leftout_mae_values.append(mae)
                
            except Exception as e:
                print(f"⚠️ Failed to compute left-out MAE for combination {row.get('combination')}: {e}")
                leftout_mae_values.append(np.nan)
        
        results['leftout_mae'] = leftout_mae_values
        return results

    # analyze_shap_values(model, X, feature_names=None, target_name="output", n_top_features=10):

    def plot_shap_values(self, X, feature_names=None, target_name="output", n_top_features=10, plot=True, dir=None):
        """
        Plot SHAP values for the model's predictions.
        
        Args:
            X (np.ndarray): Feature matrix.
            feature_names (list, optional): Names of the features.
            target_name (str, optional): Name of the target variable.
            n_top_features (int, optional): Number of top features to display.
        """
        model= self.model
        
        sample_names = self.molecule_names if hasattr(model, 'molecule_names') else None
        save_dir =self.paths.png_dir if hasattr(model, 'paths') and hasattr(self.paths, 'png_dir') else dir
        print(f"Saving SHAP plots to: {save_dir}, sample names: {sample_names}")
        analyze_shap_values(model, X, feature_names=feature_names, sample_names=sample_names, target_name=target_name, n_top_features=n_top_features, plot=plot, dir=save_dir)

    def calculate_q2_and_mae(
        self,
        X,
        y,
        n_splits=None,
        test_size=0.1,
        random_state=42,
        n_iterations=5,
        exclude_leftout: bool = True,
        plot: bool = False,
    ):
        """
        Cross-validation metrics: Q², MAE, RMSD.
        Optionally plots out-of-fold predictions instead of printing per-fold stats.
        """
        import numpy as np
        import matplotlib.pyplot as plt
        from sklearn.model_selection import train_test_split, LeaveOneOut, KFold
        from sklearn.pipeline import make_pipeline
        from sklearn.preprocessing import StandardScaler
        from sklearn.base import clone
        from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error

        X = np.asarray(X)
        y = np.asarray(y).ravel()

        # Exclude reserved molecules
        if exclude_leftout and hasattr(self, "molecule_names_predict") and hasattr(self, "molecule_names"):
            try:
                keep_mask = np.array(
                    [name not in set(self.molecule_names_predict) for name in self.molecule_names],
                    dtype=bool,
                )
                if keep_mask.shape[0] == X.shape[0]:
                    X, y = X[keep_mask], y[keep_mask]
            except Exception:
                pass

        n = X.shape[0]
        if X.ndim != 2 or y.ndim != 1 or y.shape[0] != n:
            raise ValueError("X must be 2D and y must be 1D with matching rows")

        def _fit_predict(train_idx, test_idx):
            pipe = make_pipeline(StandardScaler(), clone(self.model))
            pipe.fit(X[train_idx], y[train_idx])
            return pipe.predict(X[test_idx])

        q2_list, mae_list, rmsd_list = [], [], []
        all_preds, all_true, all_folds = [], [], []
        n_iter = max(1, int(n_iterations))
        base_seed = None if random_state is None else int(random_state)

        # -----------------
        # HOLDOUT
        # -----------------
        if n_splits is None or n_splits == 0:
            for it in range(n_iter):
                rs = None if base_seed is None else base_seed + it
                tr, te = train_test_split(np.arange(n), test_size=test_size, shuffle=True, random_state=rs)
                y_hat = _fit_predict(tr, te)
                q2 = float(r2_score(y[te], y_hat)) if len(te) >= 2 else float("nan")
                mae = float(mean_absolute_error(y[te], y_hat))
                rmsd = float(np.sqrt(mean_squared_error(y[te], y_hat)))
                q2_list.append(q2); mae_list.append(mae); rmsd_list.append(rmsd)
                all_preds.extend(y_hat); all_true.extend(y[te]); all_folds.extend([it]*len(te))

        # -----------------
        # LOO
        # -----------------
        elif n_splits == 1:
            loo = LeaveOneOut()
            y_pred_oof = np.empty_like(y, dtype=float)
            for tr, te in loo.split(X):
                y_pred_oof[te] = _fit_predict(tr, te)
            q2 = float(r2_score(y, y_pred_oof))
            mae = float(mean_absolute_error(y, y_pred_oof))
            rmsd = float(np.sqrt(mean_squared_error(y, y_pred_oof)))
            q2_list = [q2]; mae_list = [mae]; rmsd_list = [rmsd]
            all_preds = list(y_pred_oof); all_true = list(y); all_folds = list(range(n))

        # -----------------
        # K-FOLD
        # -----------------
        else:
            k = int(min(max(2, n // 2), int(n_splits)))
            for it in range(n_iter):
                rs = None if base_seed is None else base_seed + it
                kf = KFold(n_splits=k, shuffle=True, random_state=rs)
                y_pred_oof = np.empty_like(y, dtype=float)
                for fold, (tr, te) in enumerate(kf.split(X)):
                    y_hat = _fit_predict(tr, te)
                    y_pred_oof[te] = y_hat
                    all_preds.extend(y_hat); all_true.extend(y[te]); all_folds.extend([fold]*len(te))
                q2 = float(r2_score(y, y_pred_oof))
                mae = float(mean_absolute_error(y, y_pred_oof))
                rmsd = float(np.sqrt(mean_squared_error(y, y_pred_oof)))
                q2_list.append(q2); mae_list.append(mae); rmsd_list.append(rmsd)

        # --- Final metrics ---
        q2 = float(np.nanmean(q2_list))
        mae = float(np.mean(mae_list))
        rmsd = float(np.mean(rmsd_list))

        # -----------------
        # Plot diagnostics
        # -----------------
        if plot and len(all_preds) > 0:
            all_preds = np.array(all_preds)
            all_true = np.array(all_true)
            all_folds = np.array(all_folds)
            plt.figure(figsize=(6,6))
            scatter = plt.scatter(all_true, all_preds, c=all_folds, cmap="tab10", s=60, alpha=0.8, edgecolor="k")
            lims = [min(all_true.min(), all_preds.min())-0.5, max(all_true.max(), all_preds.max())+0.5]
            plt.plot(lims, lims, "r--", lw=2, label="Ideal")
            plt.xlabel("True values")
            plt.ylabel("Predicted values (OOF)")
            plt.title(f"CV predictions\nQ²={q2:.3f}, MAE={mae:.3f}, RMSD={rmsd:.3f}")
            plt.legend()
            plt.colorbar(scatter, label="Fold")
            plt.tight_layout()
            plt.show()

        return q2, mae, rmsd




    def fit(self, X, y, alpha=1e-5):
        """
        Train linear regression with bias and L2 regularization,
        and precompute everything for predict().
        """
        # 1) build training design matrix with bias
        X_b = np.c_[np.ones((X.shape[0], 1)), X]
        self.X_b_train = X_b

        # 2) normal eqn with regularization (no penalty on bias term)
        A = X_b.T @ X_b
        I = np.eye(A.shape[0])
        I[0, 0] = 0
        self.theta = np.linalg.inv(A + alpha * I) @ (X_b.T @ y)

        # 3) residual variance σ² = SSE / (n – p)
        residuals = y - X_b.dot(self.theta)
        n, p = X_b.shape
        self.sigma2 = (residuals**2).sum() / (n - p)

        # 4) store (XᵀX)⁻¹ for later interval widths
        self.XtX_inv = np.linalg.inv(X_b.T @ X_b)

        return self

    def predict(self, X, return_interval=False, cl=0.95):
        """
        Make point predictions, and optionally return (lower, upper)
        prediction intervals at confidence level cl.
        """
        # 1) build new design matrix
        X_b = np.c_[np.ones((X.shape[0], 1)), X]

        # 2) point predictions
        yhat = X_b.dot(self.theta)

        if not return_interval:
            return yhat
        residuals = self.target_vector - X_b.dot(self.theta)  # Assuming target_vector and features_df are stored
        self.residual_variance = np.sum(residuals ** 2) / (X_b.shape[0] - X_b.shape[1])
        # Calculate and store the covariance matrix of the coefficients
        self.model_covariance_matrix = self.residual_variance * np.linalg.inv(X_b.T @ X_b)
        # 3) t‐value for two-sided interval
        df   = self.X_b_train.shape[0] - self.X_b_train.shape[1]
        tval = t.ppf(1 - (1 - cl) / 2, df)

        # 4) standard errors: sqrt(σ² * (1 + xᵀ (XᵀX)⁻¹ x))
        #    vectorized for all rows
        #    (X_b @ XtX_inv) has shape [n_pred, p]; elementwise * X_b then sum across axis=1
        inside = 1 + np.sum((X_b @ self.XtX_inv) * X_b, axis=1)
        se_pred = np.sqrt(self.sigma2 * inside)

        # 5) intervals
        lower = yhat - tval * se_pred
        upper = yhat + tval * se_pred
        return yhat, lower, upper



    def predict_for_leftout(
    self,
    X: pd.DataFrame,
    y: Optional[object] = None,
    X_train: Optional[pd.DataFrame] = None,
    y_train: Optional[object] = None,
    *,
    calc_interval: bool = False,
    confidence_level: float = 0.95,
):
        """
        Predict on left-out samples using an already-fitted or freshly refitted model.

        Parameters
        ----------
        X : pd.DataFrame or np.ndarray
            Left-out feature matrix.
        y : array-like, optional
            Left-out target values.
        X_train, y_train : array-like, optional
            If provided, the model will be (re)fit on this data before predicting.
        calc_interval : bool, default False
            If True and the model/predict wrapper supports it, compute prediction intervals.
        confidence_level : float, default 0.95
            Confidence level for prediction intervals (if supported).

        Returns
        -------
        preds or (preds, errors) or (preds, lower, upper) or (preds, lower, upper, errors)
        """

        # ----------------------- 1) Ensure fitted (or fit now) -----------------------
        must_refit = (X_train is not None) and (y_train is not None)
        if must_refit:
            Xt = X_train.to_numpy() if hasattr(X_train, "to_numpy") else np.asarray(X_train)
            yt = np.asarray(y_train).ravel()
            self.model.fit(Xt, yt)          # fit now
            # Optionally keep for future checks
            self.X_, self.y_ = Xt, yt
        else:
            try:
                check_is_fitted(self.model)
            except Exception:
                # try cached training data, else error
                if hasattr(self, "X_") and hasattr(self, "y_"):
                    self.model.fit(self.X_, self.y_)
                else:
                    raise NotFittedError(
                        "This estimator is not fitted yet. Provide X_train,y_train or fit before calling predict_for_leftout."
                    )

        # ----------------------- 2) Align columns/schema -----------------------
        X_is_df = isinstance(X, pd.DataFrame)
        if X_is_df and isinstance(X_train, pd.DataFrame):
            # reindex to training columns; missing -> NaN
            X_aligned = X.reindex(columns=X_train.columns)
        elif X_is_df and hasattr(self, "features_df") and isinstance(self.features_df, pd.DataFrame):
            X_aligned = X.reindex(columns=self.features_df.columns)
        else:
            X_aligned = X

        # handle missing columns (NaNs after reindex)
        if isinstance(X_aligned, pd.DataFrame):
            if X_aligned.isna().any().any():
                # Simple, safe default: fill NaNs with 0.0
                # (swap to column means if you store them, e.g., self.feature_means_)
                X_aligned = X_aligned.fillna(0.0)
            X_arr = X_aligned.to_numpy()
        else:
            X_arr = np.asarray(X_aligned)

        # Ensure 2D
        if X_arr.ndim == 1:
            X_arr = X_arr.reshape(-1, 1)

        # ----------------------- 3) Apply scaler if present -----------------------
        if hasattr(self, "scaler_") and self.scaler_ is not None:
            try:
                X_arr = self.scaler_.transform(X_arr)
            except Exception:
                # if scaler fails (shape mismatch), skip silently or raise as you prefer
                pass

        # ----------------------- 4) Predict (with or without intervals) -----------------------
        # If your model class has a wrapper predict(X, return_interval=..., cl=...),
        # prefer that; else fall back to sklearn's .predict
        preds = lower = upper = None

        if calc_interval and hasattr(self, "predict"):
            # Your wrapper signature uses 'cl' in earlier code; keep both for safety
            try:
                preds, lower, upper = self.predict(X_arr, return_interval=True, cl=confidence_level)
            except TypeError:
                preds, lower, upper = self.predict(X_arr, return_interval=True, confidence_level=confidence_level)
        elif hasattr(self, "predict"):
            # If your wrapper exists but intervals not requested
            out = self.predict(X_arr, return_interval=False)
            preds = out[0] if isinstance(out, (tuple, list)) else out
        else:
            # plain sklearn estimator
            preds = self.model.predict(X_arr)

        preds = np.asarray(preds).ravel()

        # ----------------------- 5) Errors (if y given) -----------------------
        errors = None
        if y is not None:
            y_arr = y.values.ravel() if isinstance(y, (pd.Series, pd.DataFrame)) else np.asarray(y).ravel()
            if preds.shape != y_arr.shape:
                raise ValueError(f"Prediction shape {preds.shape} != y shape {y_arr.shape}")
            errors = y_arr - preds

        # ----------------------- 6) Return -----------------------
        if calc_interval and (lower is not None) and (upper is not None):
            lower = np.asarray(lower).ravel()
            upper = np.asarray(upper).ravel()
            return (preds, lower, upper, errors) if errors is not None else (preds, lower, upper)
        else:
            return (preds, errors) if errors is not None else preds

        
    def evaluate(self, X, y):
        
        ## must use model.predict() to get the predictions
        y_pred = self.model.predict(X)
        results = {}
        if 'r2' in self.metrics:
            results['r2'] = r2_score(y, y_pred)
        return results , y_pred

    def cross_validate(self, X, y):
        n_splits=self.n_splits 
        cv = KFold(n_splits=n_splits, shuffle=True, random_state=42)
        y_pred = cross_val_predict(self.model, X, y, cv=cv)
        scores = cross_validate(self.model, X, y, cv=cv, scoring=self.metrics)
        return y_pred, scores
    
    def get_coefficients_from_trained_estimator(self):
        intercept = self.theta[0]
        coefficients = self.theta[1:]  
        return intercept, coefficients
    
    def get_covariance_matrix(self,features):
        intercept = self.theta[0]
        coefficients = self.theta[1:]
        cov_matrix = self.model_covariance_matrix
        std_errors = np.sqrt(np.diag(cov_matrix))
        t_values = coefficients / std_errors[1:]  # Ignoring intercept for t-value
        degrees_of_freedom = len(self.target_vector) - len(coefficients) - 1
        p_values = 2 * (1 - stats.t.cdf(np.abs(t_values), df=degrees_of_freedom))
        # Include intercept values for std_error, t_value, and p_value
        intercept_std_error = std_errors[0]
        intercept_t_value = intercept / intercept_std_error
        intercept_p_value = 2 * (1 - stats.t.cdf(np.abs(intercept_t_value), df=degrees_of_freedom))
        coefficient_table = {
            'Estimate': np.concatenate([[intercept], coefficients]),
            'Std. Error': std_errors,
            't value': np.concatenate([[intercept_t_value], t_values]),
            'p value': np.concatenate([[intercept_p_value], p_values])
        }

        # Convert to DataFrame
        coef_df = pd.DataFrame(coefficient_table)
        coef_df.index = ['(Intercept)'] + features  # Assuming feature_names is a list of your feature names

        return coef_df



    def search_models(
        self,
        top_n: int = 20,
        n_jobs: int = -1,
        threshold: float = 0.7,
        bool_parallel: bool = False,
        required_features: Optional[Iterable[str]] = None,
        min_models_to_keep: int = 5,
        threshold_step: float = 0.05,
        min_threshold: float = 0.2,
    ):
        """
        Fit and evaluate feature combinations with an optional set of required features.
        Threshold is relaxed automatically if all Q2 are -inf.

        Returns
        -------
        pd.DataFrame
            Top results, sorted by q2 then r2, limited to top_n.
        """
        app = self.app

        # Generate combinations up-front
        self.get_feature_combinations()
        all_combos: List[Tuple[str, ...]] = list(self.features_combinations or [])
        req: Set[str] = set(required_features or [])

        # Decide parallel
        cpu_count = multiprocessing.cpu_count()
        effective_jobs = 1 if (cpu_count == 1 or not bool_parallel) else (n_jobs if n_jobs != -1 else cpu_count)
        print(f"Using {effective_jobs} jobs for evaluation. Found {cpu_count} cores.")

        # Load existing results (avoid re-doing work)
        try:
            existing_results = _to_df(load_results_from_db(self.db_path))
            print(f"Loaded {len(existing_results)} existing results from DB.")
            # normalize combination keys
            if "combination" in existing_results.columns:
                done_combos = set(map(str, existing_results["combination"].tolist()))
            else:
                done_combos = set()
        except Exception:
            existing_results = pd.DataFrame()
            done_combos = set()

        def _filter_new_combos() -> List[Tuple[str, ...]]:
            """Not in DB and contains all required features."""
            out = []
            for combo in all_combos:
                key = str(combo)
                if key in done_combos:
                    continue
                if req and not req.issubset(set(combo)):
                    continue
                out.append(combo)
            return out

        def _evaluate_block(threshold: float) -> pd.DataFrame:
            """
            Evaluate only combos not yet done & satisfying required features.
            If parallel, assume each call writes to DB; then reload.
            """
            results_df = existing_results.copy()

            combos_to_run = _filter_new_combos()
            print(f"Combos to run: {len(combos_to_run)}, done_combos: {len(done_combos)}")
            if not combos_to_run:
                print(f"No new combinations to evaluate at threshold {threshold:.3f}.")
                return results_df

            print(f"Evaluating {len(combos_to_run)} new combos with R2 >= {threshold:.3f}...")
            if effective_jobs == 1:
                new_results = []
                for combo in tqdm(combos_to_run, desc=f"Threshold {threshold:.3f} (single-core)"):
                    try:
                        res = fit_and_evaluate_single_combination_regression(self, combo, threshold)
                        new_results.append(res)
                    except Exception:
                        # swallow or log; keep going
                        pass
                new_df = _to_df(new_results)
                results_df = pd.concat([results_df, new_df], ignore_index=True) if not new_df.empty else results_df
               
            else:
                Parallel(n_jobs=effective_jobs)(
                    delayed(fit_and_evaluate_single_combination_regression)(self, combo, threshold)
                    for combo in tqdm(combos_to_run, desc=f"Threshold {threshold:.3f} (parallel)")
                )
                # reload DB to include new results
                results_df = _to_df(load_results_from_db(self.db_path))
          
            return results_df

        # ---- initial pass -------------------------------------------------------
        results = _evaluate_block(threshold)
        results = _sort_results(results)
      

        # If everything's -inf on Q2, relax threshold and retry (bounded loop)
        attempts = 0
        threshold = float(threshold)
        while _is_all_q2_neg_inf(results) and threshold > min_threshold:
            print("All Q2 values are -inf, lowering R2 threshold and retrying...")
            best = _best_r2(results)
            # choose the next threshold: either step down from best or by a fixed step
            if best is not None:
                threshold = max(min_threshold, best - 0.15)
            else:
                threshold = max(min_threshold, threshold - threshold_step)
            print(f"New threshold: {threshold:.3f}")
            results = _evaluate_block(threshold)
            results = _sort_results(results)
            attempts += 1
            if attempts > 5:
                # avoid infinite loops
                break

        # keep at least min_models_to_keep even if top_n is small
        results = results.head(max(top_n, min_models_to_keep))

        # Show table
        if not results.empty:
            print_models_regression_table(results, app, self)

        self.results=results

        return results



import numpy as np
from sklearn.base import BaseEstimator, ClassifierMixin
from statsmodels.miscmodels.ordinal_model import OrderedModel

class OrdinalLogisticRegression(BaseEstimator, ClassifierMixin):
    def __init__(self, distr='logit'):
        """
        Custom wrapper for OrderedModel from statsmodels.
        
        Parameters:
        ------------
        distr: str, optional (default='logit')
            The distribution function for the ordinal logistic regression.
            Available options: 'logit', 'probit', 'loglog', 'cloglog', 'cauchit'.
        """
        self.distr = distr
        self.model = None
        self.result = None
        self.classes_ = None

    def fit(self, X, y):
        self.classes_ = np.unique(y)
        self.model = OrderedModel(y, X, distr=self.distr)
        self.result = self.model.fit(method='bfgs', maxiter=1000, disp=False)
        return self

    def predict(self, X):
        # Get predicted probabilities for each class
        class_probs = self.result.predict(exog=X, which='prob')
        # Get the index of the class with the highest probability
        y_pred_indices = class_probs.argmax(axis=1)
        # Map indices to class labels
        y_pred = self.classes_[y_pred_indices]
        return y_pred



    def predict_proba(self, X):
        # Get predicted probabilities for each class
        class_probs = self.result.predict(exog=X, which='prob')
      
        return class_probs


    def get_params(self, deep=True):
        """
        Get parameters for this estimator.
        
        Returns:
        ------------
        params: dict
            Dictionary of parameters.
        """
        return {"distr": self.distr}

    def set_params(self, **params):
        """
        Set the parameters of this estimator.
        
        Parameters:
        ------------
        params: dict
            Dictionary of parameters to set.
            
        Returns:
        ------------
        self: object
            Returns self.
        """
        for key, value in params.items():
            setattr(self, key, value)
        return self


class ClassificationModel:

    def __init__(self, csv_filepaths, process_method='one csv', y_value='class', names_column=None, leave_out=None, min_features_num=2, max_features_num=None, n_splits=5, metrics=None, return_coefficients=False, ordinal=False, exclude_columns=None, app=None,db_path='results'):
        self.csv_filepaths = csv_filepaths
        self.process_method = process_method
        self.y_value = y_value
        self.leave_out = leave_out
        self.names_column = names_column
        self.min_features_num = min_features_num
        self.max_features_num = max_features_num
        self.metrics = metrics if metrics is not None else ['accuracy', 'precision', 'recall', 'f1', 'roc_auc', 'mcfadden_r2']
        self.return_coefficients = return_coefficients
        self.n_splits = n_splits
        self.ordinal = ordinal
        self.app=app
        name=name = os.path.splitext(os.path.basename(self.csv_filepaths.get('features_csv_filepath')))[0]
        
        self.setup_random_seeds()
        self.paths = prepare_run_dirs(
            base_dir="runs",
            dataset_name=name,
            y_value=self.y_value,
            tag="cls"
        )
        self.db_path = resolve_db_path(db_path, name, self.paths.db)
        
        create_results_table_classification(self.db_path)
        if csv_filepaths:
      
            if process_method == 'one csv':
                self.process_features_csv()
            elif process_method == 'two csvs':
                self.process_features_csv()
                self.process_target_csv(csv_filepaths.get('target_csv_filepath'))
            self.compute_correlation()
            self.scaler = StandardScaler()
            if exclude_columns is None:
                exclude_columns = []

            # Identify columns to scale
            # --- Identify columns to scale: only non-integers (floats/continuous) ---
            continuous_cols = [
                col for col in self.features_df.columns
                if col not in (exclude_columns or [])
                and not pd.api.types.is_integer_dtype(self.features_df[col])
            ]

            # Columns to leave untouched
            other_cols = [c for c in self.features_df.columns if c not in continuous_cols]

            print(f"Scaling only continuous (non-integer) columns: {continuous_cols}")
            print(f"Leaving integer columns untouched: {other_cols}")

            # Apply scaling only to continuous columns
            self.scaler = ColumnTransformer(
                transformers=[
                    ("num", StandardScaler(), continuous_cols),
                    ("passthrough", "passthrough", other_cols),
                ],
                remainder="drop"
            )

            # Fit and transform the data
            scaled_array = self.scaler.fit_transform(self.features_df)
            new_columns = continuous_cols + other_cols
            self.features_df = pd.DataFrame(scaled_array, columns=new_columns, index=self.features_df.index)

            self.leave_out_samples(self.leave_out, keep_only=False)
            self.determine_number_of_features()
            # self.get_feature_combinations()
            

        if ordinal:
            self.model = OrdinalLogisticRegression()
        else:    
            self.model = LogisticRegression(solver='lbfgs', random_state=42)


    def setup_random_seeds(self) -> None:
        """
        Set random seeds across Python, NumPy, and optionally PyTorch 
        to ensure reproducible results. Also stores a `self.random_state`
        attribute for use in other methods.
        """
        self.random_state = 42
        # Python built-in
        random.seed(self.random_state)

        # NumPy
        np.random.seed(self.random_state)

        # Make Python hash-based operations deterministic
        os.environ["PYTHONHASHSEED"] = str(self.random_state)

    def load_results(self, sort_key: str = "mcfadden_r2", ascending: bool = False):
        """
        Load all CSV results from the run's tables directory into DataFrames,
        optionally sorting them by a metric (default: McFadden's R²).

        Parameters
        ----------
        sort_key : str
            Column name to sort by (default 'mcfadden_r2').
        ascending : bool
            Sort order. Default False = best (highest values first).

        Returns
        -------
        results_dict : dict[str, pd.DataFrame]
            Mapping {filename (no extension): DataFrame}.
        """
        results_dict = {}
        for csv_file in os.listdir(self.paths.db):
            if not csv_file.endswith(".csv"):
             
                continue
            print(f"Loading results from {csv_file}...")
            try:
                full_path = self.paths.db / csv_file
                print(f"  Full path: {full_path}")
                df = pd.read_csv(full_path)
                print(f"  {len(df)} rows, columns: {list(df.columns)}, head:\n{df.head(2)}")
                if sort_key in df.columns:
                    df = df.sort_values(by=sort_key, ascending=ascending).reset_index(drop=True)
                else:
                    print(f"⚠️ Column '{sort_key}' not found in {csv_file.name}, returning unsorted.")
                results_dict[csv_file.stem] = df
            except Exception as e:
                print(f"❌ Failed to load {csv_file}: {e}")

        for df in results_dict.values():
            df.drop(columns=['threshold', 'model', 'predictions'], errors='ignore', inplace=True)

        return results_dict



    def compute_correlation(self, correlation_threshold=0.8, vif_threshold=None):
        app = self.app
        self.corr_matrix = self.features_df.corr()
        high_corr_features = self._get_highly_correlated_features(
            self.corr_matrix, threshold=correlation_threshold
        )
        if high_corr_features:
            # Report findings
            msg = (
                f"\n--- Correlation Report ---\n"
                f"Features with correlation above {correlation_threshold}:\n"
                f"{list(high_corr_features)}\n"
            )
            if app:
                app.show_result(msg)
            print(msg)

            # === Always visualize if no app ===
            visualize = True
            if app:
                visualize = messagebox.askyesno(
                    title="Visualize Correlated Features?",
                    message="Would you like to see a heatmap of the correlation among these features?"
                )

            if visualize:
                # interactive_corr_heatmap_with_highlights(self.features_df)
                sub_corr = self.corr_matrix.loc[sorted(high_corr_features), sorted(high_corr_features)]
                interactive_corr_heatmap(sub_corr, title=f"Correlation Heatmap (|r| > {correlation_threshold})")

            # === Ask to drop correlated features ===
            drop_corr = False
            if app:
                drop_corr = messagebox.askyesno(
                    title="Drop Correlated Features?",
                    message=(
                        f"Features above correlation {correlation_threshold}:\n"
                        f"{list(high_corr_features)}\n\n"
                        "Do you want to randomly drop some of these correlated features?"
                    )
                )

            if drop_corr:
                count_to_drop = len(high_corr_features) // 2
                features_to_drop = random.sample(list(high_corr_features), k=count_to_drop)

                drop_msg = (
                    f"\nRandomly selected {count_to_drop} features to drop:\n"
                    f"{features_to_drop}\n"
                )
                if app:
                    app.show_result(drop_msg)
                print(drop_msg)

                self.features_df.drop(columns=features_to_drop, inplace=True)
                self.features_list = self.features_df.columns.tolist()

                remaining_msg = f"Remaining features: {self.features_list}\n"
                if app:
                    app.show_result(remaining_msg)
                print(remaining_msg)
            else:
                no_drop_msg = "\nCorrelated features were not dropped.\n"
                if app:
                    app.show_result(no_drop_msg)
                print(no_drop_msg)

        else:
            msg = "\nNo features exceeded the correlation threshold.\n"
            if app:
                app.show_result(msg)
            print(msg)


    def simi_sampler_(self, class_label, compare_with=0, sample_size=None, plot=False):
        """
        Perform similarity sampling to balance classes in the dataset.
        """
        data = pd.concat([self.features_df, self.target_vector], axis=1)

        # If sample_size not provided, use the minimum class count
        if sample_size is None:
            sample_size = data["class"].astype("category").value_counts().min()
            print(f"Sample size not provided. Using minimum class count: {sample_size}")
        self.indices = simi_sampler(data, class_label, compare_with, sample_size, plot)

        return self.indices
        # # Rebuild X and y without the sampled indices
        # X = self.features_df.reset_index(drop=True)
        # y = self.target_vector.reset_index(drop=True)

        # mask = np.ones(len(X), dtype=bool)
        # mask[np.array(self.indices, dtype=int)] = False

        # self.features_df = X.iloc[mask]
        # self.target_vector = y.iloc[mask]

        # # Update molecule_names to exclude indices
        # self.molecule_names = [
        #     self.molecule_names[i] for i in range(len(self.molecule_names)) if i not in self.indices
        # ]


        
        
    def stratified_sampling_(self,plot=False, sample_size=None):
        data=pd.concat([self.features_df,self.target_vector],axis=1)
        df = stratified_sampling_with_plots(data, plot, sample_size)
        self.features_df = df.drop(columns=[self.y_value])
        self.target_vector = df[self.y_value]

    def get_feature_combinations(self):
        self.features_combinations = list(get_feature_combinations(self.features_list, self.min_features_num, self.max_features_num))


    def determine_number_of_features(self):
        total_features_num = self.features_df.shape[0]
        self.max_features_num = set_max_features_limit(total_features_num, self.max_features_num)
    

    def leave_out_samples(self, leave_out=None, keep_only=False):
        """
        Remove or retain specific rows from features_df/target_vector based on indices or molecule names.
        Stores the removed (or retained) samples into:
            - self.leftout_samples
            - self.leftout_target_vector
            - self.molecule_names_predict

        Parameters:
            leave_out (int, str, list[int|str]): Indices or molecule names to leave out (or to keep if keep_only=True).
            keep_only (bool): If True, only the given indices/names are kept, all others are left out.
        """
        print('Debugging leave')
        if leave_out is None:
            return

        # Normalize input to a list
        if isinstance(leave_out, (int, np.integer, str)):
            leave_out = [leave_out]
        else:
            leave_out = list(leave_out)

        # Determine indices from mixed inputs (names or indices)
        indices = []
        name_to_index = {name: idx for idx, name in enumerate(self.molecule_names)}

        for item in leave_out:
            if isinstance(item, (int, np.integer)):
                indices.append(int(item))

            elif isinstance(item, str):
                # try exact match first
                if item in name_to_index:
                    indices.append(name_to_index[item])
                else:
                    # try to parse as int
                    try:
                        idx_val = int(item)
                        indices.append(idx_val)
                    except ValueError:
                        raise ValueError(
                            f"Leave-out entry '{item}' is not a known molecule name "
                            f"and cannot be parsed as an integer index."
                        )

            else:
                raise ValueError("Items in leave_out must be integers or strings.")

        # --- stash the relevant data ---
        selected_features = self.features_df.iloc[indices].copy()
        selected_target   = self.target_vector.iloc[indices].copy()
        selected_names    = [self.molecule_names[i] for i in indices]

        # — ensure we have the full 106 columns in the right order —
        if hasattr(self, 'feature_names'):
            selected_features = selected_features.reindex(
                columns=self.feature_names,
                fill_value=0
            )

        # locate index labels for dropping
        idx_labels = self.features_df.index[indices]

        if keep_only:
            # Everything *not* in indices becomes the "left out" set
            self.leftout_samples       = self.features_df.drop(index=idx_labels).copy()
            self.leftout_target_vector = self.target_vector.drop(index=idx_labels).copy()
            self.molecule_names_predict = [
                n for i, n in enumerate(self.molecule_names) if i not in indices
            ]

            # The new main set is just our selected rows (already reindexed)
            self.features_df     = selected_features
            self.target_vector   = selected_target
            self.molecule_names  = selected_names
        else:
            # The selected rows become the "left out" set (already reindexed)
            self.leftout_samples       = selected_features
            self.leftout_target_vector = selected_target
            self.molecule_names_predict = selected_names

            # Drop them from the main set
            self.features_df    = self.features_df.drop(index=idx_labels)
            self.target_vector  = self.target_vector.drop(index=idx_labels)
            self.molecule_names = [
                n for i, n in enumerate(self.molecule_names) if i not in indices
            ]

        # debug prints
        print(f"Left out samples:   {self.molecule_names_predict}")
        print(f"Remaining samples:  {self.molecule_names}")
        
    def process_features_csv(
    self,
    drop_columns: Optional[Union[str, List[str]]] = None,
    drop_smiles: bool = True,
) -> dict:
        """
        Load and process a features dataset from:
        - a CSV file path (string),
        - a dict with key 'features_csv_filepath',
        - a single-path list/tuple,
        - or an in-memory pandas DataFrame.

        Initializes:
            - self.molecule_names : list[str]
            - self.target_vector  : pd.Series
            - self.features_df    : pd.DataFrame (numeric-only, deduplicated)
            - self.features_list  : list[str]

        Parameters
        ----------
        drop_columns : str or list[str|int], optional
            Additional columns to drop (case-insensitive for names; ints treated as positions).
        drop_smiles : bool, default=True
            Drop SMILES-like columns if found.

        Returns
        -------
        dict
            Summary of processing steps.
        """
        import os
        from typing import List, Optional, Union
        import numpy as np
        import pandas as pd

        # ---------- Load source ----------
        csv_source = getattr(self, "csv_filepaths", None)
        if csv_source is None:
            raise AttributeError(
                "`self.csv_filepaths` must be set to a CSV path (str), a dict, a pandas DataFrame, "
                "or a single-path list/tuple."
            )

        # Normalize a single-path list/tuple into a string
        if isinstance(csv_source, (list, tuple)):
            if len(csv_source) == 1 and isinstance(csv_source[0], str):
                csv_source = csv_source[0]
            else:
                raise TypeError(
                    "`csv_filepaths` list/tuple must contain exactly one string path; "
                    f"got {csv_source!r}"
                )

        path_used = None

        if isinstance(csv_source, pd.DataFrame):
            df = csv_source.copy()
            # we will still honor self.names_column if provided
        elif isinstance(csv_source, dict):
            csv_path = csv_source.get("features_csv_filepath")
            if not csv_path or not os.path.exists(csv_path):
                raise FileNotFoundError(f"CSV not found at: {csv_path!r}")
            df = pd.read_csv(csv_path)
            path_used = os.path.abspath(csv_path)
        elif isinstance(csv_source, str):
            if not os.path.exists(csv_source):
                raise FileNotFoundError(f"CSV not found at: {csv_source!r}")
            df = pd.read_csv(csv_source)
            path_used = os.path.abspath(csv_source)
        else:
            raise TypeError("`self.csv_filepaths` must be a str, dict, DataFrame, or single-path list/tuple.")

        if df.empty:
            raise ValueError("The provided dataset is empty.")

        # ---------- Helpers ----------
        def _case_insensitive_find(col_name: str, columns: pd.Index):
            lc = str(col_name).lower()
            for c in columns:
                if str(c).lower() == lc:
                    return c
            return None

        def _resolve_many(names: List[str], columns: pd.Index) -> List[str]:
            resolved = []
            for n in names:
                hit = _case_insensitive_find(n, columns)
                if hit is not None:
                    resolved.append(hit)
            return list(dict.fromkeys(resolved))  # dedupe

        def _first_match(candidates: List[str], columns: pd.Index):
            for c in candidates:
                hit = _case_insensitive_find(c, columns)
                if hit is not None:
                    return hit
            return None

        # ---------- Target & name resolution ----------
        y_value = getattr(self, "y_value", None)
        if y_value is None:
            raise AttributeError("`self.y_value` must be set to the target column name.")

        # Allow caller-provided drop_columns or self.drop_columns
        if drop_columns is None and hasattr(self, "drop_columns"):
            drop_columns = getattr(self, "drop_columns")

        if isinstance(drop_columns, str):
            drop_columns = [drop_columns]
        elif drop_columns is None:
            drop_columns = []

        provided_names_col = getattr(self, "names_column", None)
        name_candidates = ["name", "names", "molecule", "molecule_name", "mol", "compound", "sample", "id"]

        # Choose names column (honor provided_names_col first)
        names_col = None
        if provided_names_col:
            names_col = _case_insensitive_find(provided_names_col, df.columns)
            if names_col is None:
                raise KeyError(
                    f"names_column {provided_names_col!r} not found. Available: {list(df.columns)}"
                )
        else:
            # Use a meaningful index if not a default RangeIndex; otherwise try heuristics then fallback to first column
            if not isinstance(df.index, pd.RangeIndex):
                names_col = None  # will use index
            else:
                names_col = _first_match(name_candidates, df.columns) or df.columns[0]

        target_col = _case_insensitive_find(y_value, df.columns)
        if target_col is None:
            raise KeyError(f"Target column {y_value!r} not found. Available: {list(df.columns)}")

        # ---------- Build features frame ----------
        cols_to_drop = {c for c in [names_col, target_col] if c is not None}
        features_df_raw = df.drop(columns=[c for c in cols_to_drop if c in df.columns])
        n_features_total = features_df_raw.shape[1]

        dropped_smiles_like_cols, dropped_user_columns = [], []

        # Drop SMILES-like columns
        if drop_smiles:
            smiles_candidates = ["smiles", "canonical_smiles", "random_smiles"]
            smiles_to_drop = _resolve_many(smiles_candidates, features_df_raw.columns)
            if smiles_to_drop:
                features_df_raw = features_df_raw.drop(columns=smiles_to_drop)
                dropped_smiles_like_cols.extend(smiles_to_drop)

        # Drop user-specified columns (support names or integer positions)
        if drop_columns:
            dc = list(drop_columns)
            by_pos = [c for c in dc if isinstance(c, (int, np.integer))]
            by_name = [c for c in dc if not isinstance(c, (int, np.integer))]

            resolved_user_cols = _resolve_many(by_name, features_df_raw.columns)

            if by_pos:
                pos_names = [
                    features_df_raw.columns[i] for i in by_pos
                    if 0 <= int(i) < features_df_raw.shape[1]
                ]
                resolved_user_cols.extend(pos_names)

            resolved_user_cols = list(dict.fromkeys(resolved_user_cols))  # dedupe
            if resolved_user_cols:
                features_df_raw = features_df_raw.drop(columns=resolved_user_cols)
                dropped_user_columns.extend(resolved_user_cols)

        # Deduplicate columns
        features_df_raw = features_df_raw.loc[:, ~features_df_raw.columns.duplicated()]

        # ---------- Coerce to numeric ----------
        coerced = features_df_raw.apply(pd.to_numeric, errors="coerce")
        all_nan_cols = [c for c in coerced.columns if coerced[c].isna().all()]
        numeric_df = coerced.drop(columns=all_nan_cols, errors="ignore").select_dtypes(include=[np.number])
        dropped_non_numeric = [
            c for c in features_df_raw.columns if c not in numeric_df.columns and c not in all_nan_cols
        ]

        # ---------- Assign outputs ----------
   

        if names_col is None:
            # Convert Index -> Series so we can use Series string ops safely
            names_series = pd.Series(df.index, copy=False)
            names_label = "<index>"
        else:
            names_series = df[names_col]
            names_label = names_col

        # Clean names to strings and remove NaNs cleanly
        names_series = pd.Series(names_series, copy=False)
        names_series = names_series.where(names_series.notna(), "")
        self.molecule_names = (
            names_series.astype(str)
            .str.replace(r'^\s*nan\s*$', '', regex=True)
            .tolist()
        )

        target_series_numeric = pd.to_numeric(df[target_col], errors="coerce")
        self.target_vector = target_series_numeric if not target_series_numeric.isna().all() else df[target_col]
        self.features_df = numeric_df
        self.features_list = numeric_df.columns.tolist()

        # ---------- Summary ----------
        name_src_note = "from index" if names_col is None else f"from column '{names_col}'"
        summary = {
            "path": path_used,
            "names_column": names_label,
            "target_column": target_col,
            "n_rows": df.shape[0],
            "n_features_total": n_features_total,
            "n_features_kept_numeric": numeric_df.shape[1],
            "dropped_non_numeric_cols": dropped_non_numeric,
            "dropped_all_nan_cols": all_nan_cols,
            "dropped_smiles_like_cols": dropped_smiles_like_cols,
            "dropped_user_columns": dropped_user_columns,
        }

        # ---------- Output/log ----------
        msg = (
            f"Processed {'DataFrame' if path_used is None else os.path.basename(path_used)}\n"
            f"Molecule names: {name_src_note} | Target: {target_col}\n"
            f"Rows: {df.shape[0]} | Features: {n_features_total} → {numeric_df.shape[1]}"
        )
        for key, lst in {
            "SMILES-like": dropped_smiles_like_cols,
            "User-specified": dropped_user_columns,
            "Non-numeric": dropped_non_numeric,
            "All-NaN": all_nan_cols,
        }.items():
            if lst:
                msg += f"\nDropped {key} cols: {lst}"

        app = getattr(self, "app", None)
        if app:
            try:
                app.show_result(msg)
            except Exception:
                print(msg)
        else:
            print(msg)

        return summary



    
    def compute_multicollinearity(self,correlation_threshold=0.8):
        app=self.app
        self.corr_matrix = self.features_df.corr()

        # Identify highly-correlated features above correlation_threshold
        high_corr_features = self._get_highly_correlated_features(
            self.corr_matrix, threshold=correlation_threshold
        )


        if high_corr_features:
            # Show correlation report
            app.show_result(f"\n--- Correlation Report ---\n")
            app.show_result(
                f"Features with correlation above {correlation_threshold}:\n"
                f"{list(high_corr_features)}\n"
            )
            visualize_corr = messagebox.askyesno(
            title="Visualize Correlated Features?",
            message=(
                "Would you like to see a heatmap of the correlation among these features?"
            )
        )
            if visualize_corr:
                # Subset the correlation matrix for the correlated features only
                sub_corr = self.corr_matrix.loc[list(high_corr_features), list(high_corr_features)]
                
                # Create a heatmap with Seaborn
                fig, ax = plt.subplots(figsize=(6, 5))
                sns.heatmap(sub_corr, annot=False, cmap='coolwarm', square=True, ax=ax)
                ax.set_title(f"Correlation Heatmap (>{correlation_threshold})")
                plt.tight_layout()
                plt.show()

            # Ask user if they want to drop them (yes/no)
            drop_corr = messagebox.askyesno(
                title="Drop Correlated Features?",
                message=(
                    f"Features above correlation {correlation_threshold}:\n"
                    f"{list(high_corr_features)}\n\n"
                    "Do you want to randomly drop some of these correlated features?"
                )
            )
            if drop_corr:
                # Decide how many to drop. (Here: half the set, randomly)
                count_to_drop = len(high_corr_features) // 2
                features_to_drop = random.sample(list(high_corr_features), k=count_to_drop)

                app.show_result(f"\nRandomly selected {count_to_drop} features to drop:")
                app.show_result(f"{features_to_drop}\n")

                # Remove from DataFrame
                self.features_df.drop(columns=features_to_drop, inplace=True)
                self.features_list = self.features_df.columns.tolist()

                app.show_result(f"Remaining features: {self.features_list}\n")
            else:
                app.show_result("\nCorrelated features were not dropped.\n")
        else:
            app.show_result("\nNo features exceeded the correlation threshold.\n")

       

    def _compute_vif(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Compute Variance Inflation Factor (VIF) for each numeric column in df.
        Skips computation (returns NaN VIFs) if fewer than 2 usable features remain.
        """
        import numpy as np
        import pandas as pd
        from sklearn.preprocessing import StandardScaler
        from statsmodels.stats.outliers_influence import variance_inflation_factor

        # 1) Keep numeric, drop inf/NaN rows
        X = df.select_dtypes(include=[np.number]).copy()
        X = X.replace([np.inf, -np.inf], np.nan).dropna(axis=0)

        # 2) Drop zero-variance columns
        zero_var = X.columns[X.nunique() <= 1].tolist()
        if zero_var:
            X = X.drop(columns=zero_var)

        # 3) Guard: need at least two columns to compute VIF
        if X.shape[1] < 2:
            msg = (
                f"VIF not computed: need >=2 numeric, non-constant columns; "
                f"found {X.shape[1]}. "
                f"Dropped zero-variance: {zero_var}" if zero_var else ""
            )
            print(msg)
            return pd.DataFrame({"variables": X.columns, "VIF": [np.nan] * X.shape[1]})

        # 4) Standardize
        scaler = StandardScaler()
        Xs = pd.DataFrame(scaler.fit_transform(X), columns=X.columns, index=X.index)

        # 5) Optional: quick check for rank deficiency
        rank = np.linalg.matrix_rank(Xs.values)
        if rank < Xs.shape[1]:
            print(f"Warning: linear dependence detected (rank {rank} < {Xs.shape[1]}). "
                "Some VIFs may be very large or infinite.")

        # 6) Compute VIFs safely
        vifs = []
        for i in range(Xs.shape[1]):
            try:
                v = variance_inflation_factor(Xs.values, i)
            except Exception:
                v = np.inf
            vifs.append(v)

        vif_df = pd.DataFrame({"variables": Xs.columns, "VIF": vifs}).sort_values("VIF", ascending=False)
        print("VIF summary:\n", vif_df)
        return vif_df



    def _get_highly_correlated_features(self, corr_matrix, threshold=0.8):
        """
        Identify any features whose pairwise correlation is above the threshold.
        Returns a set of the implicated feature names.
        """
        corr_matrix_abs = corr_matrix.abs()
        columns = corr_matrix_abs.columns
        high_corr_features = set()

        # We only need to look at upper triangular part to avoid duplication
        for i in range(len(columns)):
            for j in range(i + 1, len(columns)):
                if corr_matrix_abs.iloc[i, j] > threshold:
                    # Add both columns in the pair
                    high_corr_features.add(columns[i])
                    high_corr_features.add(columns[j])
        
        return high_corr_features
        

    def process_target_csv(self, csv_filepath):
        target_vector_unordered = pd.read_csv(csv_filepath)[self.y_value]
        self.target_vector = target_vector_unordered.loc[self.molecule_names]

    def fit(self, X, y):

        self.model.fit(X, y)

    def predict(self, X):

        return self.model.predict(X)

    def compute_multicollinearity(self, vif_threshold=5.0):
        """
        Compute the Variance Inflation Factor (VIF) for each feature in the dataset.
        """
        # Compute VIF
        vif_results = self._compute_vif(self.features_df)
        
        # Identify features with high VIF
        high_vif_features = vif_results[vif_results['VIF'] > vif_threshold]
        
        return vif_results
    
    def calculate_vif(self):
        X = self.features_df
        vif_data = pd.DataFrame()
        vif_data["feature"] = X.columns
        vif_data["VIF"] = [variance_inflation_factor(X.values, i) for i in range(X.shape[1])]
        return vif_data

    def evaluate(self, X, y):
        # Evaluate the classifier using different metrics
        y_pred = self.predict(X)
        accuracy = accuracy_score(y, y_pred)
        precision = precision_score(y, y_pred,average='weighted', zero_division=0)
        recall = recall_score(y, y_pred,average='weighted', zero_division=0)
        f1 = f1_score(y, y_pred,average='weighted', zero_division=0)
        # auc = roc_auc_score(y, self.model.predict_proba(X)[:, 1])
        mcfadden_r2_var = self.mcfadden_r2(X, y)
        results = {
            'accuracy': accuracy,
            'precision': precision,
            'recall': recall,
            'f1_score': f1,
            'mcfadden_r2': mcfadden_r2_var
            #'auc': auc
        }
        return results , y_pred



    def cross_validation(self, X, y, n_splits=5):
        
        # Initialize lists to store metrics
        accuracy_list = []
        f1_list = []
        mcfadden_r2_list = []

        # Use appropriate KFold splitter
        kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)

        # Loop over the folds manually
        for train_index, test_index in kf.split(X):
            X_train, X_test = X[train_index], X[test_index]
            y_train, y_test = y[train_index], y[test_index]

            # Fit the model on the training data
            if self.ordinal:
               
                result = self.model.fit(X_train, y_train)
                # Predict probabilities and classes
                y_pred = result.predict(X_test)
                # y_pred = np.argmax(y_pred_prob, axis=1)
                mcfadden_r2= self.mcfadden_r2(X_test, y_test)
                mcfadden_r2_list.append(mcfadden_r2)

            else:
                
                result = self.model.fit(X_train, y_train)
                y_pred = result.predict(X_test)
                mcfadden_r2 = self.mcfadden_r2(X_test, y_test)
                mcfadden_r2_list.append(mcfadden_r2)

            # Compute accuracy and F1-score
            accuracy = accuracy_score(y_test, y_pred)
            f1 = f1_score(y_test, y_pred, average='weighted', zero_division=0)

            accuracy_list.append(accuracy)
            f1_list.append(f1)

        # Calculate mean values
        avg_accuracy = np.mean(accuracy_list)
        avg_f1 = np.mean(f1_list)
        avg_mcfadden_r2 = np.mean(mcfadden_r2_list)

        return avg_accuracy, avg_f1, avg_mcfadden_r2

    def cross_validation(self, X, y, n_splits: int = 5, verbose: bool = True):
        """
        Stratified CV for classification with robust McFadden R^2:
        - StratifiedKFold with n_splits capped by min class count
        - Accuracy, weighted F1
        - McFadden R^2 computed stably:
            * sklearn: from predict_proba on the test fold
            * ordinal (statsmodels OrderedModel): from fitted fold probabilities (no refit)
            * skipped (NaN) if the test fold has <2 classes
        """
        import numpy as np
        import pandas as pd
        from sklearn.model_selection import StratifiedKFold
        from sklearn.metrics import accuracy_score, f1_score

        def _log(msg: str):
            if verbose:
                print(msg)
                app = getattr(self, "app", None)
                if app:
                    try:
                        app.show_result(str(msg))
                    except Exception:
                        pass

        # --- Make arrays; drop rows with non-finite values ---
        X_arr = np.asarray(X)
        y_arr = np.asarray(y)
        mask = np.isfinite(X_arr).all(axis=1) & np.isfinite(y_arr)
        X_arr, y_arr = X_arr[mask], y_arr[mask]

        # --- Cap splits by minimum class count to avoid single-class folds ---
        classes, y_enc = np.unique(y_arr, return_inverse=True)
        if classes.size < 2:
            _log("[CV] Single-class data; returning NaNs.")
            return float("nan"), float("nan"), float("nan")

        min_per_class = np.bincount(y_enc).min()
        eff_splits = max(2, min(n_splits, int(min_per_class)))
        if eff_splits != n_splits:
            _log(f"[CV] Adjusted n_splits: requested={n_splits}, effective={eff_splits} (min per-class count = {min_per_class})")

        skf = StratifiedKFold(n_splits=eff_splits, shuffle=True, random_state=42)

        accs, f1s, r2s = [], [], []

        for fold, (tr, te) in enumerate(skf.split(X_arr, y_enc), start=1):
            X_train, X_test = X_arr[tr], X_arr[te]
            y_train, y_test = y_arr[tr], y_arr[te]

            # ---- Fit model on the training split ----
            result = self.model.fit(X_train, y_train)

            # ---- Predict labels for metrics ----
            if getattr(self, "ordinal", False):
                # statsmodels OrderedModel: result.predict often returns probabilities by default
                y_pred_raw = result.predict(X_test)
                # If it's a probability matrix, argmax -> predicted class index; map back to labels from y_train
                if hasattr(y_pred_raw, "ndim") and y_pred_raw.ndim == 2:
                    pred_idx = np.asarray(y_pred_raw).argmax(axis=1)
                    train_cats = pd.Categorical(y_train)
                    y_pred = np.asarray([train_cats.categories[i] for i in pred_idx])
                else:
                    # assume it already returned class labels
                    y_pred = np.asarray(y_pred_raw)
            else:
                # sklearn-style predictor
                y_pred = np.asarray(result.predict(X_test))

            # ---- Accuracy & F1 (weighted) ----
            accs.append(accuracy_score(y_test, y_pred))
            f1s.append(f1_score(y_test, y_pred, average="weighted", zero_division=0))

            # ---- McFadden R^2 on the test split (stable computation) ----
            # Skip if test split lacks at least 2 classes
            if np.unique(y_test).size < 2:
                r2s.append(np.nan)
                _log(f"[CV][fold {fold}] Skipping McFadden R^2 (single-class test split).")
                continue

            try:
                if getattr(self, "ordinal", False):
                    # Try to get per-class probabilities from the fitted OrderedModel
                    try:
                        probs = result.model.predict(result.params, X_test, which="prob")  # (n,k)
                    except Exception:
                        # Some versions return probs via result.predict with kw
                        probs = result.predict(X_test, which="prob")
                    probs = np.asarray(probs)
                    k = probs.shape[1]

                    # Map y_test to indices consistent with training categories
                    train_cats = pd.Categorical(y_train)
                    class_index = {cat: i for i, cat in enumerate(train_cats.categories)}
                    y_idx = np.array([class_index[yy] for yy in y_test], dtype=int)

                    eps = 1e-15
                    p_model = np.clip(probs[np.arange(len(y_idx)), y_idx], eps, 1 - eps)
                    LL_model = float(np.sum(np.log(p_model)))

                    counts = np.bincount(y_idx, minlength=k).astype(float)
                    class_probs = counts / counts.sum()
                    p_null = np.clip(class_probs[y_idx], eps, 1 - eps)
                    LL_null = float(np.sum(np.log(p_null)))

                    r2 = 1.0 - (LL_model / LL_null)
                    r2s.append(r2)
                else:
                    # sklearn path: compute from predicted probabilities (no refit on test fold)
                    if hasattr(result, "predict_proba"):
                        proba = np.asarray(result.predict_proba(X_test))  # (n,k)
                        cls = getattr(result, "classes_", getattr(self.model, "classes_", None))
                        if cls is None:
                            raise AttributeError("Classifier lacks classes_ for aligning probabilities.")

                        idx_map = {c: i for i, c in enumerate(cls)}
                        y_idx = np.array([idx_map[yy] for yy in y_test], dtype=int)

                        eps = 1e-15
                        p_model = np.clip(proba[np.arange(len(y_idx)), y_idx], eps, 1 - eps)
                        LL_model = float(np.sum(np.log(p_model)))

                        k = proba.shape[1]
                        counts = np.bincount(y_idx, minlength=k).astype(float)
                        class_probs = counts / counts.sum()
                        p_null = np.clip(class_probs[y_idx], eps, 1 - eps)
                        LL_null = float(np.sum(np.log(p_null)))

                        r2 = 1.0 - (LL_model / LL_null)
                        r2s.append(r2)
                    else:
                        # No probabilities available -> cannot compute McFadden reliably
                        r2s.append(np.nan)
                        _log(f"[CV][fold {fold}] Classifier has no predict_proba; McFadden R^2 set to NaN.")
            except Exception as e:
                _log(f"[CV][fold {fold}] McFadden R^2 computation failed: {e}")
                r2s.append(np.nan)

        # --- Averages (ignore NaNs for R^2) ---
        avg_accuracy = float(np.mean(accs)) if accs else float("nan")
        avg_f1       = float(np.mean(f1s))  if f1s  else float("nan")
        avg_r2       = float(np.nanmean(r2s)) if r2s else float("nan")

        return avg_accuracy, avg_f1, avg_r2


    def search_models(
        self,
        top_n: int = 50,
        n_jobs: int = 1,
        threshold: float = 0.5,
        bool_parallel: bool = False,
        required_features: Optional[Iterable[str]] = None,
        min_models_to_keep: int = 5,
    ):
        """
        Fit and evaluate classification feature combinations with an optional set of required features.
        Results are cached in a SQLite database and merged with any existing results.

        Parameters
        ----------
        top_n : int
            Number of top models to return.
        n_jobs : int
            Number of parallel jobs (-1 = use all cores).
        threshold : float
            Minimum McFadden R² threshold for evaluation.
        bool_parallel : bool
            Whether to use joblib Parallel for evaluation.
        required_features : list[str], optional
            Require all combos to include these features.
        min_models_to_keep : int
            Always keep at least this many models, even if top_n is smaller.

        Returns
        -------
        pd.DataFrame
            Top results sorted by McFadden R².
        """
        app = self.app

        # Generate all candidate feature combinations
        self.get_feature_combinations()
        all_combos: List[Tuple[str, ...]] = list(self.features_combinations or [])
        req: Set[str] = set(required_features or [])

        # Decide on parallel execution
        cpu_count = multiprocessing.cpu_count()
        effective_jobs = (
            1
            if (cpu_count == 1 or not bool_parallel)
            else (n_jobs if n_jobs != -1 else cpu_count)
        )
        print(f"Using {effective_jobs} jobs for evaluation. Found {cpu_count} cores.")

        # Load existing results
        try:
            existing_results = load_results_from_db(self.db_path, table="classification_results")
            print(f"Loaded {len(existing_results)} existing results from DB.")
            done_combos = set(map(str, existing_results["combination"].tolist()))
        except Exception:
            existing_results = pd.DataFrame()
            done_combos = set()
            print("[INFO] No existing results found, starting fresh.")

        def _filter_new_combos() -> List[Tuple[str, ...]]:
            """Filter combos: not in DB, must contain required features if any."""
            out = []
            for combo in all_combos:
                key = str(combo)
                if key in done_combos:
                    continue
                if req and not req.issubset(set(combo)):
                    continue
                out.append(combo)
            return out

        def _evaluate_block(threshold: float) -> pd.DataFrame:
            """Evaluate unseen combos at given threshold."""
            results_df = existing_results.copy()
            combos_to_run = _filter_new_combos()
            print(f"Combos to run: {len(combos_to_run)}, already done: {len(done_combos)}")

            if not combos_to_run:
                print(f"No new combos left at threshold {threshold:.3f}.")
                return results_df

            print(f"Evaluating {len(combos_to_run)} new combos with McFadden R² >= {threshold:.3f}...")
            if effective_jobs == 1:
                new_results = []
                for combo in tqdm(combos_to_run, desc=f"Threshold {threshold:.3f} (single-core)"):
                    try:
                        print(f"[INFO] Evaluating combo: {combo}")
                        result = fit_and_evaluate_single_combination_classification(
                            self, combo, threshold=threshold
                        )
                        insert_result_into_db_classification(
                            self.db_path,
                            combo,
                            result["scores"],
                            threshold,
                            csv_path="classification_results.csv",
                        )
                        new_results.append(result)
                    except Exception as e:
                        print(f"[WARN] Combo {combo} failed: {e}")
                if new_results:
                    results_df = pd.concat([results_df, pd.DataFrame(new_results)], ignore_index=True)
            else:
                Parallel(n_jobs=effective_jobs)(
                    delayed(fit_and_evaluate_single_combination_classification)(
                        self, combo, threshold
                    )
                    for combo in tqdm(combos_to_run, desc=f"Threshold {threshold:.3f} (parallel)")
                )
                results_df = load_results_from_db(self.db_path, table="classification_results")

            return results_df

        # ---- Run evaluation ------------------------------------------------
        results = _evaluate_block(threshold)

        if results.empty:
            print("[INFO] No models were successfully evaluated.")
            return results

        # Sort results (by McFadden R² descending)
        results = results.sort_values(by="mcfadden_r2", ascending=False)

        # Keep at least min_models_to_keep
        results = results.head(max(top_n, min_models_to_keep))

        # Show results
        if not results.empty:
            print_models_classification_table(results, app, self)
            self.combinations_list = results["combination"].tolist()

        # ---- Left-out evaluation ------------------------------------------
        if self.leave_out:
            X = self.leftout_samples.to_numpy()
            y = self.leftout_target_vector.to_numpy()
            self.fit(X, y)
            preds = self.predict(X)
            df_lo = pd.DataFrame(
                {
                    "sample_name": self.molecule_names_predict,
                    "true": y.ravel(),
                    "predicted": np.array(preds).ravel(),
                }
            )
            if app:
                app.show_result("\n\nPredictions on left-out samples\n\n")
                app.show_result(df_lo.to_markdown(tablefmt="pipe", index=False))
            else:
                print(df_lo.to_markdown(tablefmt="pipe", index=False))
        self.results=results

        return results


    
    def iterative_cross_validation(self, combination, n_iter=10, n_splits=5):
        """
        Perform iterative cross-validation on a selected model and feature combination.

        Parameters:
        ------------
        model: statsmodels model object
            The fitted model to be evaluated.

        combination: tuple or list
            The feature combination used in the model.

        n_iter: int, optional (default=10)
            Number of iterations for cross-validation.

        n_splits: int, optional (default=5)
            Number of folds in each cross-validation iteration.

        Returns:
        ------------
        results: dict
            A dictionary containing:
            - 'overall_avg_accuracy': Overall average accuracy across all iterations.
            - 'overall_avg_f1_score': Overall average F1 score across all iterations.
            - 'overall_avg_mcfadden_r2': Overall average McFadden's R-squared across all iterations.
            - 'iteration_metrics': List of metrics for each iteration.
        """
        # Initialize lists to store overall metrics
        overall_accuracy_list = []
        overall_f1_list = []
        overall_mcfadden_r2_list = []
        iteration_metrics = []

        for i in range(n_iter):
            # For each iteration, perform cross-validation with a different random state
            random_state = 42 + i  # Change the random state for each iteration
            kf = KFold(n_splits=n_splits, shuffle=True, random_state=random_state)

            # Initialize lists for this iteration
            accuracy_list = []
            f1_list = []
            mcfadden_r2_list = []
            try:
                selected_features = self.features_df[list(combination)]
            except:
                selected_features = self.features_df[_parse_tuple_string(combination)]

            X = selected_features.to_numpy()
            y = self.target_vector.to_numpy()
   
            for train_index, test_index in kf.split(X):
                X_train, X_test = X[train_index], X[test_index]
                y_train, y_test = y[train_index], y[test_index]

                if self.ordinal:
                
                    result = self.model.fit(X_train, y_train)
                
                    # Predict probabilities and classes
                    y_pred = result.predict(X_test)
                    mcfadden_r2 = self.mcfadden_r2(X_test, y_test)
                else:
                    result = self.model.fit(X_train, y_train)
             
                    y_pred = result.predict(X_test)
                    mcfadden_r2 = self.mcfadden_r2(X_test, y_test)
                    mcfadden_r2_list.append(mcfadden_r2)


                # Compute accuracy and F1-score
                accuracy_cv = accuracy_score(y_test, y_pred)
                f1_cv = f1_score(y_test, y_pred, average='weighted', zero_division=0)

                accuracy_list.append(accuracy_cv)
                f1_list.append(f1_cv)

            # Calculate average metrics for this iteration
            avg_accuracy = np.mean(accuracy_list)
            avg_f1 = np.mean(f1_list)
            avg_mcfadden_r2 = np.mean(mcfadden_r2_list)

            # Store the metrics for this iteration
            iteration_metrics.append({
                'iteration': i + 1,
                'avg_accuracy': avg_accuracy,
                'avg_f1_score': avg_f1,
                'avg_mcfadden_r2': avg_mcfadden_r2
            })

            # Add to overall metrics
            overall_accuracy_list.append(avg_accuracy)
            overall_f1_list.append(avg_f1)
            overall_mcfadden_r2_list.append(avg_mcfadden_r2)

        # Calculate overall average metrics
        overall_avg_accuracy = np.mean(overall_accuracy_list)
        overall_avg_f1 = np.mean(overall_f1_list)
        overall_avg_mcfadden_r2 = np.mean(overall_mcfadden_r2_list)

        results = {
            'overall_avg_accuracy': overall_avg_accuracy,
            'overall_avg_f1_score': overall_avg_f1,
            'overall_avg_mcfadden_r2': overall_avg_mcfadden_r2,
            'iteration_metrics': iteration_metrics
        }

        return results


    def mcfadden_r2(self, X, y, verbose: bool = True, alpha: float = 1e-6):
        """
        Robust McFadden's R^2 using statsmodels with separation handling.
        R2 = 1 - (LL_model / LL_null)

        - Binary: Logit
        - Multiclass: MNLogit
        - Fallback to ridge penalty (fit_regularized) on separation/singularity
        - LL_null is computed analytically from class frequencies (stable)
        """
        import warnings
        import numpy as np
        import pandas as pd
        import statsmodels.api as sm
        from statsmodels.tools.sm_exceptions import PerfectSeparationWarning
        from numpy.linalg import LinAlgError

        # ---------------------- helpers ----------------------
        def _clean_rows(Xdf: pd.DataFrame, yser: pd.Series):
            Xdf = pd.DataFrame(Xdf)
            yser = pd.Series(yser)

            # y may be non-numeric; use notna instead of np.isfinite
            mask_y = yser.notna()

            # Only check finiteness on numeric X columns
            numX = Xdf.select_dtypes(include=[np.number])
            mask_x = np.ones(len(Xdf), dtype=bool)
            if not numX.empty:
                mask_x &= np.isfinite(numX).all(axis=1).to_numpy()

            mask = mask_y & mask_x
            return Xdf.loc[mask].copy(), yser.loc[mask].copy()

        def _encode_X(Xdf: pd.DataFrame) -> pd.DataFrame:
            # booleans -> ints
            for c in Xdf.columns:
                if Xdf[c].dtype == bool:
                    Xdf[c] = Xdf[c].astype(int)

            # non-numeric -> one-hot with drop_first (avoid dummy trap)
            nonnum = [c for c in Xdf.columns if not np.issubdtype(Xdf[c].dtype, np.number)]
            if nonnum:
                Xdf = pd.get_dummies(Xdf, columns=nonnum, drop_first=True)

            # drop zero-variance & exact duplicates
            zero_var = [c for c in Xdf.columns if Xdf[c].nunique(dropna=False) <= 1]
            if zero_var and verbose:
                print(f"[R2] Dropping zero-variance columns: {zero_var}")
            if zero_var:
                Xdf = Xdf.drop(columns=zero_var)

            # drop duplicate columns
            dup_mask = Xdf.T.duplicated()
            if dup_mask.any() and verbose:
                dups = list(Xdf.columns[dup_mask])
                print(f"[R2] Dropping duplicate columns: {dups}")
            Xdf = Xdf.loc[:, ~dup_mask]

            return Xdf

        def _repair_singularity(Xf: pd.DataFrame, label: str, tol: float = 1e-8) -> Optional[pd.DataFrame]:
            # Drop zero-variance and duplicates
            Xf = Xf.loc[:, Xf.nunique(dropna=False) > 1]
            Xf = Xf.loc[:, ~Xf.T.duplicated()]

            # Approx collinearity check via correlation
            corr = Xf.corr().abs()
            upper = corr.where(np.triu(np.ones(corr.shape), k=1).astype(bool))

            to_drop = [col for col in upper.columns if any(upper[col] > (1 - tol))]
            if to_drop:
                print(f"[INFO] {label}: dropping highly collinear columns: {to_drop}")
                Xf = Xf.drop(columns=to_drop)

            rank = np.linalg.matrix_rank(Xf.values, tol=tol)
            print(f"[OK] {label}: reduced to {Xf.shape[1]} cols, rank {rank}")
            return Xf


        def _fit_full_binary(Xf, ybin):
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always", category=PerfectSeparationWarning)
                try:
                    res = sm.Logit(ybin, Xf).fit(disp=False, method="newton", maxiter=2000)
                except (LinAlgError, np.linalg.LinAlgError, ValueError):
                    res = None
                need_fallback = res is None or any(
                    issubclass(getattr(wi, "category", type(None)), PerfectSeparationWarning) for wi in w
                )
            if need_fallback:
                if verbose:
                    print("[R2] Falling back to regularized Logit (L2)")
                res = sm.Logit(ybin, Xf).fit_regularized(alpha=alpha, L1_wt=0.0, maxiter=2000, disp=False)
            return res

        def _fit_full_mnlogit(Xf, ycodes):
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always", category=PerfectSeparationWarning)
                try:
                    res = sm.MNLogit(ycodes, Xf).fit(disp=False, method="newton", maxiter=2000)
                except Exception:
                    res = None
                need_fallback = res is None or any(
                    issubclass(getattr(wi, "category", type(None)), PerfectSeparationWarning) for wi in w
                )
            if need_fallback:
                if verbose:
                    print("[R2] Falling back to regularized MNLogit (L2)")
                res = sm.MNLogit(ycodes, Xf).fit_regularized(alpha=alpha, L1_wt=0.0, maxiter=2000, disp=False)
            return res

        # ---------------------- pipeline ----------------------
        X, y = _clean_rows(X, y)
        X = _encode_X(X)

        n, p = X.shape
        if n == 0 or p == 0:
            if verbose:
                print("[R2] No usable rows/features after cleaning.")
            return float("nan")

        # Add intercept
        X_full = sm.add_constant(X, has_constant="add")

        cats = pd.Categorical(y)
        k = len(cats.categories)

        # ---------------------- fit full model ----------------------
        try:
            if k == 2:
                # Binary
                y_bin = (y == cats.categories[1]).astype(int)

                # Repair singularity if needed
                X_try = _repair_singularity(X_full, label="binary") 
                if X_try is None or X_try.empty:
                    return float("nan")
                else:
                    X_full = X_try
                res_full = _fit_full_binary(X_try, y_bin)
                LL_model = float(res_full.llf)
                counts = np.bincount(y_bin, minlength=2).astype(float)

            else:
                # Multiclass
                y_codes = cats.codes  # 0..k-1

                X_try = _repair_singularity(X_full, label="multiclass")
                if X_try is None or X_try.empty:
                    return float("nan")
                else:
                    X_full = X_try
                res_full = _fit_full_mnlogit(X_try, y_codes)
                LL_model = float(res_full.llf)
                counts = np.bincount(y_codes, minlength=k).astype(float)

        except Exception as e:
            if verbose:
                print(f"[WARN] Fit failed: {e}")
            return float("nan")

        # ---------------------- null model ----------------------
        nobs = counts.sum()
        if nobs <= 0 or (counts > 0).sum() < 2:
            return float("nan")

        probs = counts / nobs
        # LL_null = sum_c n_c * log(p_c), for appearing classes only
        LL_null = float(np.sum(counts[counts > 0] * np.log(probs[counts > 0])))

        # Guard against degenerate null likelihood
        if not np.isfinite(LL_null) or np.isclose(LL_null, 0.0):
            LL_null = 1e-15

        # ---------------------- McFadden R^2 ----------------------
        r2 = 1.0 - (LL_model / LL_null)
        return r2



    def _compute_log_likelihood(self, y_indices, y_prob):
        """
        Compute the log-likelihood of the model.

        Parameters:
        ------------
        y_indices: array-like, shape (n_samples,)
            True labels encoded as integers starting from 0.

        y_prob: array-like, shape (n_samples, n_classes)
            Predicted probabilities for each class.

        Returns:
        ------------
        log_likelihood: float
            The log-likelihood value.
        """
        # Clip probabilities to avoid log(0)
        eps = 1e-15
        y_prob = np.clip(y_prob, eps, 1 - eps)

        # Get the probabilities of the true classes
        prob_true_class = y_prob[np.arange(len(y_indices)), y_indices]

        # Compute log-likelihood
        log_likelihood = np.sum(np.log(prob_true_class))
        return log_likelihood









