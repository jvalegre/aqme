import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.linear_model import LogisticRegression
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from adjustText import adjust_text
from itertools import combinations
import sqlite3
import os 
from typing import  List , Literal ,Union, Optional
import statsmodels.api as sm
from statsmodels.stats.stattools import durbin_watson
from statsmodels.stats.diagnostic import het_breuschpagan
from scipy.stats import shapiro, probplot

def _add_drop_status(simi_table: pd.DataFrame, class_label, keep_pos, label_col_hint="Label"):
    df = simi_table.copy()
    # pick a label column for annotations
    if label_col_hint not in df.columns:
        if "Name" in df.columns:
            df["Label"] = df["Name"].astype(str)
        else:
            df["Label"] = df.index.astype(str)
    # status: other / kept / dropped
    is_target = (df["class"] == class_label)
    kept_mask = is_target & df["_pos_"].isin(keep_pos)
    drop_mask = is_target & ~df["_pos_"].isin(keep_pos)
    df["status"] = np.where(drop_mask, "dropped", np.where(kept_mask, "kept", "other"))
    return df


def plot_similarity_scatter_with_status(df: pd.DataFrame, x_col: str, title: str,
                                        annotate_dropped: bool = True, annotate_kept: bool = False):
    # clean
    dd = df[[x_col, "class", "status", "Label"]].copy()
    dd = dd.dropna(subset=[x_col, "class"])
    dd = dd[np.isfinite(dd[x_col])]
    if dd.empty:
        print("[plot] nothing to plot.")
        return

    fig, ax = plt.subplots(figsize=(10, 5))

    # OTHER classes (low emphasis)
    other = dd[dd["status"] == "other"]
    if not other.empty:
        ax.scatter(other[x_col], other["class"], s=24, alpha=0.25, color="0.6", label="other")

    # KEPT (target class)
    kept = dd[dd["status"] == "kept"]
    if not kept.empty:
        ax.scatter(kept[x_col], kept["class"], s=40, alpha=0.95, label=f"kept (n={len(kept)})",
                   edgecolors="k", linewidths=0.6)
        if annotate_kept:
            for _, r in kept.iterrows():
                ax.text(r[x_col], r["class"], str(r["Label"]), fontsize=8, va="bottom", ha="left")

    # DROPPED (target class) — highlight strongly
    dropped = dd[dd["status"] == "dropped"]
    if not dropped.empty:
        ax.scatter(dropped[x_col], dropped["class"], s=64, marker="x", linewidths=1.8,
                   color="crimson", label=f"dropped (n={len(dropped)})")
        if annotate_dropped:
            for _, r in dropped.iterrows():
                ax.text(r[x_col], r["class"], str(r["Label"]), fontsize=8, color="crimson",
                        va="bottom", ha="left")

    ax.set_title(title)
    ax.set_xlabel(x_col)
    ax.set_ylabel("class")
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    plt.show()



from typing import Optional, Literal, List
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler


def simi_sampler(
    data: pd.DataFrame,
    class_label,
    compare_with: Optional[str] = "self",
    sample_size: Optional[int] = None,
    plot: bool = False,
    return_kind: Literal["pos", "label", "name"] = "pos",
) -> List:
    """
    Select a subset of rows from `data` for the specified `class_label`,
    spreading selections across the chosen similarity range, and
    RETURN THE INDICES TO DROP for that class.

    Parameters
    ----------
    data : pd.DataFrame
        Must contain 'class' column and numeric feature columns.
    class_label : str | int
        The class whose samples will be downsampled.
    compare_with : {"self", class_name}, default="self"
        Which similarity axis to use.
        - "self" → similarity to own class centroid.
        - class_name → similarity to another class centroid.
    sample_size : int, optional
        Number of samples to KEEP. Defaults to min(class counts).
    plot : bool, default=False
        If True, scatter plots before/after selection.
    return_kind : {"pos","label","name"}
        Format of indices to return:
        - "pos": positional indices (0..n-1, for .iloc)
        - "label": original DataFrame index values
        - "name": molecule names (or fallback strings)

    Returns
    -------
    List
        Indices (according to `return_kind`) to DROP.
    """
    if "class" not in data.columns:
        raise KeyError("data must contain a 'class' column.")

    # keep original index and pos
    df = data.copy()
    df["_orig_index_"] = df.index
    df.reset_index(drop=True, inplace=True)
    df["_pos_"] = np.arange(len(df), dtype=int)

    # choose a label column
    label_col = None
    for cand in ("Name", "molecule"):
        if cand in df.columns:
            label_col = cand
            break
    if label_col is None:
        df["Name"] = df["_orig_index_"].astype(str)
        label_col = "Name"

    # feature columns
    meta_cols = {"class", "flag", "Sample", "_orig_index_", "_pos_", "Name"}
    var_cols = [c for c in df.columns if c not in meta_cols]

    # numeric features only
    sampler_data = df[var_cols].apply(pd.to_numeric, errors="coerce")
    sampler_data = sampler_data.dropna(axis=1, how="all")
    sampler_data = sampler_data.fillna(sampler_data.mean())

    # drop zero-variance cols
    sampler_data = sampler_data.loc[:, sampler_data.nunique(dropna=False) > 1]
    if sampler_data.empty:
        raise ValueError("No usable numeric feature columns after preprocessing.")

    # scale
    scaler = StandardScaler()
    sampler_scaled = pd.DataFrame(
        scaler.fit_transform(sampler_data),
        columns=sampler_data.columns,
        index=df.index,
    )

    # compute class centroids
    unique_classes = df["class"].unique()
    class_vec, class_mag = {}, {}
    for cls in unique_classes:
        m = sampler_scaled[df["class"] == cls].mean(axis=0).values
        class_vec[cls] = m
        class_mag[cls] = np.linalg.norm(m)

    # similarity column
    if compare_with == "self":
        # similarity to own centroid
        sim = []
        for idx, row in sampler_scaled.iterrows():
            c = df.at[idx, "class"]
            v, m = class_vec[c], class_mag[c]
            rn = np.linalg.norm(row.values)
            sim.append(0.0 if (m == 0.0 or rn == 0.0) else float(np.dot(v, row) / (m * rn)))
        df["_sim_"] = sim
    else:
        # similarity to another class centroid
        if compare_with not in class_vec:
            raise KeyError(f"compare_with={compare_with} not in {list(class_vec.keys())}")
        v, m = class_vec[compare_with], class_mag[compare_with]
        sim = []
        for _, row in sampler_scaled.iterrows():
            rn = np.linalg.norm(row.values)
            sim.append(0.0 if (m == 0.0 or rn == 0.0) else float(np.dot(v, row) / (m * rn)))
        df["_sim_"] = sim

    # target class only
    simi_class = df.loc[df["class"] == class_label, "_sim_"]
    if simi_class.empty:
        raise ValueError(f"No rows found for class_label={class_label}")

    n_class = len(simi_class)
    k = min(sample_size if sample_size is not None else n_class, n_class)

    # positions of all rows in class
    all_pos_in_class = df.loc[df["class"] == class_label, "_pos_"].tolist()

    # pick k evenly spaced along similarity axis
    if k >= n_class:
        keep_pos = all_pos_in_class
    else:
        steps = np.linspace(simi_class.min(), simi_class.max(), num=k)
        sim_series = simi_class.copy()
        keep_pos = []
        for step in steps:
            idx = (sim_series - step).abs().idxmin()
            keep_pos.append(int(df.at[idx, "_pos_"]))
            sim_series = sim_series.drop(idx)

    drop_pos = [p for p in all_pos_in_class if p not in set(keep_pos)]

    # optionally plot (only if plotting helpers exist)
    if plot:
        try:
            import matplotlib.pyplot as plt

            plt.scatter(df["_sim_"], df["class"], c="gray", alpha=0.5)
            plt.scatter(
                df.loc[df["_pos_"].isin(keep_pos), "_sim_"],
                df.loc[df["_pos_"].isin(keep_pos), "class"],
                c="green",
                label="kept",
            )
            plt.scatter(
                df.loc[df["_pos_"].isin(drop_pos), "_sim_"],
                df.loc[df["_pos_"].isin(drop_pos), "class"],
                c="red",
                label="dropped",
            )
            plt.title(f"Similarity distribution for class {class_label}")
            plt.legend()
            plt.show()
        except Exception:
            print("[simi_sampler] Plotting skipped (matplotlib not available).")

    # return DROP indices
    if return_kind == "pos":
        return drop_pos
    elif return_kind == "label":
        return df.loc[drop_pos, "_orig_index_"].tolist()
    elif return_kind == "name":
        return df.loc[drop_pos, label_col].astype(str).tolist()
    else:
        raise ValueError("return_kind must be one of {'pos','label','name'}")




def stratified_sampling(df, group, size):
    """
    Performs stratified sampling on a DataFrame based on the grouping variable.

    Parameters:
    - df (pd.DataFrame): The DataFrame to sample from.
    - group (str): The name of the column to group by for stratification.
    - size (int, float, dict): Desired sample size:
                               - If int, it specifies a fixed number of samples per group.
                               - If float < 1, it specifies the proportion to sample from each group.
                               - If dict, it specifies the exact number of samples for each group.

    Returns:
    - pd.DataFrame: A stratified sample of the original DataFrame.
    """
    
    # Group by the specified column
    grouped = df.groupby(group)

    # Initialize an empty DataFrame to store the stratified sample
    stratified_sample = pd.DataFrame()

    # Iterate through each group in the DataFrame
    for name, group_data in grouped:
        group_size = len(group_data)

        # Determine the number of samples based on the type of `size`
        if isinstance(size, int):
            n_samples = min(size, group_size)  # Take all if group size is smaller than `size`
        elif isinstance(size, float) and size < 1:
            n_samples = max(1, int(size * group_size))  # Take a proportion
        elif isinstance(size, dict):
            n_samples = min(size.get(name, 0), group_size)  # Take from dict, but ensure we don’t exceed group size
        else:
            raise ValueError("size must be an integer, a float < 1, or a dictionary")

        # Sample from the group
        sampled_data = group_data.sample(n=n_samples, replace=False)
        stratified_sample = pd.concat([stratified_sample, sampled_data])

    # Reset the index of the sampled DataFrame
    stratified_sample = stratified_sample.reset_index(drop=True)
    
    return stratified_sample


def plot_distribution_scatter(data, group, x_column='index', title='Class Distribution'):
    """
    Plots the class distribution using a scatter plot, with text annotations and adjusted labels.

    Parameters:
    - data (pd.DataFrame): The DataFrame whose distribution to plot.
    - group (str): The name of the column representing classes/groups.
    - x_column (str): Column to use for x-axis (default is 'index' for basic class scatter).
    - title (str): The title of the plot.
    """
    # Ensure the index column exists for x-axis if not provided
    if x_column == 'index':
        data = data.reset_index()

    # Prepare plot
    plt.figure(figsize=(10, 5))
    
    # Scatterplot with class hue
    sns.scatterplot(
        x=data[x_column],
        y=data[group],
        hue=data[group]
    )
    
    # Add text labels to each point, adjusting for overlaps
    texts = []
    for cls in data[group].unique():
        class_data = data[data[group] == cls]
    
        for i in range(class_data.shape[0]):
            texts.append(
                plt.text(
                    class_data[x_column].iloc[i],
                    class_data[group].iloc[i],
                    class_data.iloc[:,1].iloc[i],  # Use index or add a 'Name' column if needed
                    fontsize=9
                )
            )
    
    # Adjust text positions to avoid overlap
    adjust_text(texts, arrowprops=dict(arrowstyle='->', color='red'))
    
    # Set title, labels, and legend
    plt.title(title)
    plt.xlabel(f"{x_column}")
    plt.ylabel(f"{group}")
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))  # Legend on the right side
    
    # Show the plot
    plt.show()


def stratified_sampling_with_plots(df, size ,group='class' , plot=True):
    """
    Performs stratified sampling and plots the distribution before and after sampling.

    Parameters:
    - df (pd.DataFrame): The DataFrame to sample from.
    - group (str): The name of the column to group by for stratification.
    - size (int, float, dict): Desired sample size (same as stratified_sampling function).
    
    Returns:
    - pd.DataFrame: A stratified sample of the original DataFrame.
    """
    if plot:
    # Plot distribution before sampling
        plot_distribution_scatter(df, group, title="Class Distribution Before Stratified Sampling")

    # Perform stratified sampling
    stratified_sample = stratified_sampling(df, group, size)
    
    if plot:
    # Plot distribution after sampling
        plot_distribution_scatter(stratified_sample, group, title="Class Distribution After Stratified Sampling")
    
    return stratified_sample



def create_results_table_classification(db_path='results.db'):
    """Create the classification_results table if it does not already exist."""
    db_exists = os.path.isfile(db_path)

    if db_exists:
        print(f"Database already exists at: {db_path}")
    else:
        print(f"Database does not exist. It will be created at: {db_path}")

    # Create table
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    c.execute('''
        CREATE TABLE IF NOT EXISTS classification_results (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            combination TEXT,
            accuracy REAL,
            precision REAL,
            recall REAL,
            f1_score REAL,
            mcfadden_r2 REAL,
            avg_mcfadden_r2 REAL,
            avg_accuracy REAL,
            avg_f1_score REAL,
            threshold REAL
        );
    ''')
    print("Table 'classification_results' has been ensured to exist.")
    
    conn.commit()
    conn.close()

def insert_result_into_db_classification(db_path, combination, results, threshold, csv_path='classification_results.csv'):
    """
    Insert classification metrics into the SQLite database and append to a CSV file.

    Args:
        db_path (str): Path to SQLite database.
        combination (str): Feature combination used.
        results (dict): Dictionary with keys: 'accuracy', 'precision', 'recall', 'f1_score', 'mcfadden_r2'.
        threshold (float): Threshold used.
        csv_path (str): Path to CSV file for backup/logging.
    """
    accuracy = results.get('accuracy')
    precision = results.get('precision')
    recall = results.get('recall')
    f1 = results.get('f1_score')
    mcfadden_r2 = results.get('mcfadden_r2')
    avg_mcfadden_r2 = results.get('avg_mcfadden_r2')
    avg_accuracy = results.get('avg_accuracy')
    avg_f1_score = results.get('avg_f1_score')
    # Insert into SQLite
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    c.execute('''
        INSERT INTO classification_results (
            combination, accuracy, precision, recall, f1_score, mcfadden_r2, threshold, avg_accuracy, avg_f1_score, avg_mcfadden_r2
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?);
    ''', (str(combination), accuracy, precision, recall, f1, mcfadden_r2, threshold, avg_accuracy, avg_f1_score, avg_mcfadden_r2))
    conn.commit()
    conn.close()

    # Append to CSV
    result_dict = {
        'combination': [str(combination)],
        'accuracy': [accuracy],
        'precision': [precision],
        'recall': [recall],
        'f1_score': [f1],
        'mcfadden_r2': [mcfadden_r2],
        'avg_mcfadden_r2': [avg_mcfadden_r2],
        'avg_accuracy': [avg_accuracy],
        'avg_f1_score': [avg_f1_score],
        'threshold': [threshold]
    }

    result_df = pd.DataFrame(result_dict)
   
    if not os.path.isfile(csv_path):
        result_df.to_csv(csv_path, index=False, mode='w')
    else:
        result_df.to_csv(csv_path, index=False, mode='a', header=False)



def create_results_table(db_path='results.db'):
    """Create the regression_results table if it does not already exist."""
    # Check if DB file already exists
    db_exists = os.path.isfile(db_path)
    
    if db_exists:
        print(f"Database already exists at: {db_path}")
    else:
        print(f"Database does not exist. It will be created at: {db_path}")

    # Connect and create table if needed
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    c.execute('''
        CREATE TABLE IF NOT EXISTS regression_results (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            combination TEXT,
            r2 REAL,
            q2 REAL,
            mae REAL,
            rmsd REAL,

            threshold REAL,
            model TEXT,
            predictions TEXT
        );
    ''')
    print("Table 'regression_results' has been ensured to exist.")
    
    conn.commit()
    conn.close()


def insert_result_into_db_regression(db_path, combination, r2, q2, mae, rmsd, threshold, model, predictions, csv_path='results.csv'):
    """
    Insert one row of results into the SQLite database and append to a CSV file.
    
    Args:
        db_path (str): Path to the SQLite database.
        combination (str): Feature combination.
        formula (str): Model formula.
        r2 (float): R-squared value.
        q2 (float): Q-squared value.
        mae (float): Mean Absolute Error.
        rmsd (float): Root Mean Squared Deviation.
        threshold (float): Threshold used.
        csv_path (str): Path to the CSV file.
    """
    # Insert into SQLite database
    # print(f'Inserting results for combination: {combination} | R2: {r2} | Q2: {q2}')
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    c.execute('''
        INSERT INTO regression_results (combination, r2, q2, mae, rmsd, threshold, model, predictions)
        VALUES (?, ?, ?, ?, ?, ?, ?, ?);
    ''', (str(combination), r2, q2, mae, rmsd, threshold, str(model), str(predictions)))
    conn.commit()
    conn.close()

    # Prepare data for CSV
    result_dict = {
        'combination': [str(combination)],
        'r2': [r2],
        'q2': [q2],
        'mae': [mae],
        'rmsd': [rmsd],
        'threshold': [threshold],
        'model': [model],
        'predictions': [predictions]
    }
  
    result_df = pd.DataFrame(result_dict)
    
    # Check if CSV exists; if not, write header
    if not os.path.isfile(csv_path):
        result_df.to_csv(csv_path, index=False, mode='w')
     
    else:
        result_df.to_csv(csv_path, index=False, mode='a', header=False)
       


def load_results_from_db(db_path, table='regression_results'):
    """
    Load the entire results table from the SQLite database.

    Args:
        db_path (str): Path to the SQLite database.
        table   (str): Table name ('regression_results' or 'classification_results').

    Returns:
        List[dict]: One dict per row, with column names as keys.
    """
    conn = sqlite3.connect(db_path)
    df = pd.read_sql_query(f"SELECT * FROM {table}", conn)
  
    conn.close()
    return df

def sample(x, size=None, replace=False, prob=None, random_state=None):
    """
    Draw random samples from a population.

    Parameters:
        x (int or sequence): 
            - If an integer, represents the range 1 to `x` inclusive.
            - If a sequence (list, tuple, numpy array), samples are drawn from the elements of the sequence.
        size (int, optional): 
            Number of samples to draw. 
            - If `None`, defaults to the length of `x` (only applicable when `replace=False`).
        replace (bool, optional): 
            Whether the sampling is with replacement. Defaults to `False`.
        prob (list or numpy.ndarray, optional): 
            A list of probabilities associated with each element in `x`. Must be the same length as `x` if provided.
        random_state (int, numpy.random.Generator, or None, optional): 
            Seed or random number generator for reproducibility.

    Returns:
        list: A list of sampled elements.
    """
    # Handle random state
    rng = None
    if isinstance(random_state, int):
        rng = np.random.default_rng(random_state)
    elif isinstance(random_state, np.random.Generator):
        rng = random_state
    elif random_state is not None:
        raise ValueError("random_state must be an int, numpy.random.Generator, or None.")
    
    # If x is an integer, interpret as 1 to x inclusive
    if isinstance(x, int):
        if x < 1:
            raise ValueError("When 'x' is an integer, it must be greater than or equal to 1.")
        population = list(range(1, x + 1))
    elif isinstance(x, (list, tuple, np.ndarray, pd.Series)):
        population = list(x)
    else:
        raise TypeError("x must be an integer or a sequence (list, tuple, numpy array, or pandas Series).")
    
    # Determine default size
    if size is None:
        if replace:
            size = len(population)
        else:
            size = len(population)
    
    # Validate size
    if not replace and size > len(population):
        raise ValueError("Cannot take a larger sample than population when 'replace' is False.")
    
    # Handle probabilities
    if prob is not None:
        if len(prob) != len(population):
            raise ValueError("'prob' must be the same length as 'x'.")
        prob = np.array(prob)
        if not np.isclose(prob.sum(), 1):
            raise ValueError("The sum of 'prob' must be 1.")
    
    # Perform sampling
    if rng is not None:
        sampled = rng.choice(population, size=size, replace=replace, p=prob)
    else:
        sampled = np.random.choice(population, size=size, replace=replace, p=prob)
    
    return sampled.tolist()

def assign_folds_no_empty(n_samples, n_folds, random_state=None):
    """
    Assign each data point to a fold ensuring that no fold is empty.

    Parameters:
        n_samples (int): Number of data points.
        n_folds (int): Number of folds.
        random_state (int, optional): Seed for reproducibility.

    Returns:
        list: List of fold assignments (1-based indexing).
    """
    if n_folds > n_samples:
        raise ValueError("Number of folds cannot exceed number of samples.")

    # Assign one unique data point to each fold
    initial_assignments = sample(n_samples, size=n_folds, replace=False, random_state=random_state)
    # initial_assignments are unique data points indices (1-based)

    # Initialize all assignments to 0
    fold_assignments = [0] * n_samples

    # Assign each fold its unique data point
    for fold_num, idx in enumerate(initial_assignments, start=1):
        fold_assignments[idx - 1] = fold_num  # Convert to 1-based indexing

    # Assign the remaining data points to any fold using the sample function
    remaining_indices = [i for i in range(1, n_samples + 1) if fold_assignments[i - 1] == 0]
    if remaining_indices:
        # Define uniform probabilities for simplicity; modify if needed
        fold_probs = [1.0 / n_folds] * n_folds
        sampled_folds = sample(n_folds, size=len(remaining_indices), replace=True, prob=fold_probs, random_state=random_state)
        for idx, fold_num in zip(remaining_indices, sampled_folds):
            fold_assignments[idx - 1] = fold_num

    return fold_assignments


# --- your main function, refactored -----------------------------------------

def r_squared(pred, obs, formula="corr", na_rm=False):
    """
    Compute R-squared between observed and predicted values.

    Args:
        pred (array-like): Predicted values.
        obs (array-like): Observed (actual) values.
        formula (str): Method for R² calculation: "corr" (default) or "traditional".
        na_rm (bool): If True, remove NaN values before computation.

    Returns:
        float: R-squared value.
    """
    # Convert inputs to numpy arrays for easier computation
    pred = np.array(pred)
    obs = np.array(obs)
    
    # Handle missing values (NaN) if na_rm is True
    if na_rm:
        mask = ~np.isnan(pred) & ~np.isnan(obs)
        pred = pred[mask]
        obs = obs[mask]

    n = len(pred)

    if formula == "corr":
        # Correlation-based R²
        corr_matrix = np.corrcoef(obs, pred)
        r_squared_value = corr_matrix[0, 1] ** 2

    elif formula == "traditional":
        # Traditional R²: 1 - (SS_res / SS_tot)
        ss_res = np.sum((obs - pred) ** 2)
        ss_tot = (n - 1) * np.var(obs, ddof=1)  # Sample variance
        r_squared_value = 1 - (ss_res / ss_tot)

    else:
        raise ValueError("Invalid formula type. Choose 'corr' or 'traditional'.")

    return r_squared_value


def set_max_features_limit(total_features_num, max_features_num=None):
    if max_features_num is None:
        max_features_num = int(total_features_num / 5.0)
    return max_features_num

def get_feature_combinations(features, min_features_num=2, max_features_num=None):
    max_features_num = set_max_features_limit(len(features), max_features_num)
    total_combinations = 0  # Counter to track the total number of combinations generated
    for current_features_num in range(min_features_num, max_features_num + 1):
        count_combinations = 0  # Counter for each specific number of features
        for combo in combinations(features, current_features_num):
            count_combinations += 1
            total_combinations += 1
            yield combo

def _parse_tuple_string(s: str):

    return [x.strip(" '") for x in s.strip("()").split(",")]

def _normalize_combination_to_columns(combination) -> list[str]:
    """Return a clean list of column names from tuple/list/str without using ast.
    Handles items like 'Dist(1, 2)' that contain internal commas."""
    # Already a sequence of names
    if isinstance(combination, (list, tuple, pd.Index, np.ndarray)):
        return [str(x).strip().strip("'\"") for x in combination]

    if not isinstance(combination, str):
        # Fallback: just cast to string
        s = str(combination)
    else:
        s = combination

    s = s.strip()
    # Strip outer parentheses if it's a tuple-like string
    if s.startswith("(") and s.endswith(")"):
        s = s[1:-1].strip()

    # If there are no commas at top level, it's a single feature name
    # (e.g., "dip_x" or "'dip_x'")
    # We still run the splitter below, but it will return one token.
    items = []
    buf = []
    depth = 0
    in_sq = False
    in_dq = False

    def flush():
        token = "".join(buf).strip().strip("'\"")
        if token:
            items.append(token)

    for ch in s:
        if ch == "'" and not in_dq:
            in_sq = not in_sq
            buf.append(ch)
            continue
        if ch == '"' and not in_sq:
            in_dq = not in_dq
            buf.append(ch)
            continue

        if not in_sq and not in_dq:
            if ch in "([{":
                depth += 1
            elif ch in ")]}":
                depth -= 1
            elif ch == "," and depth == 0:
                flush()
                buf = []
                continue

        buf.append(ch)

    flush()
    return items or ([] if s == "" else [s.strip().strip("'\"")])

from dataclasses import dataclass
from pathlib import Path
from datetime import datetime
import os

@dataclass
class RunPaths:
    root: Path
    db: Path
    pdf: Path
    figs: Path
    tables: Path
    exports: Path
    logs: Path

def prepare_run_dirs(
    base_dir="runs",
    dataset_name="dataset",
    y_value=None,
    tag=None,
    with_date=True,
) -> RunPaths:
    """
    Prepare (or reuse) run directories for storing results.

    Parameters
    ----------
    base_dir : str
        Root directory for all runs.
    dataset_name : str
        Dataset identifier.
    y_value : str or None
        Target variable name (optional).
    tag : str or None
        Extra tag for distinguishing runs (optional).
    with_date : bool
        If True, append today's date (YYYYMMDD).

    Returns
    -------
    RunPaths
        A namedtuple containing paths for root, db, pdf, figs, tables, exports, logs.
    """
    base_dir = Path(base_dir)
    parts = [dataset_name]
    if y_value:
        parts.append(str(y_value))
    if tag:
        parts.append(str(tag))
    if with_date:
        parts.append(datetime.now().strftime("%Y%m%d"))  # date only

    run_name = "_".join(parts)
    root = base_dir / run_name

    # Subdirectories
    db     = root / "db"
    pdf    = root / "pdf"
    figs   = root / "figs"
    tables = root / "tables"
    exports= root / "exports"
    logs   = root / "logs"

    # If directory already exists, reuse paths
    if root.exists():
        print(f"Reusing existing run directory: {root}")
    else:
        print(f"Creating new run directory: {root}")
        for p in [db, pdf, figs, tables, exports, logs]:
            p.mkdir(parents=True, exist_ok=True)

    # Ensure headless plotting
    os.environ.setdefault("MPLBACKEND", "Agg")

    return RunPaths(root, db, pdf, figs, tables, exports, logs)


def resolve_db_path(db_path_arg: Union[str, Path],
                    dataset_name: str,
                    db_dir: Optional[Path]) -> str:

    """
    If db_path_arg ends with `.db`, treat it as a full path.
    Else, treat it as a prefix and build `<prefix>_<dataset_name>.db`.
    If db_dir is provided, place the DB inside that directory.
    """
    db_path_arg = str(db_path_arg)
    if db_path_arg.lower().endswith(".db"):
        dbfile = Path(db_path_arg).name
    else:
        dbfile = f"{db_path_arg}_{dataset_name}.db"

    if db_dir is not None:
        return str((db_dir / dbfile).resolve())
    return str(Path(dbfile).resolve())

import os
import numpy as np
import pandas as pd
import statsmodels.api as sm
import matplotlib.pyplot as plt

from statsmodels.stats.diagnostic import (
    het_breuschpagan, het_white, het_goldfeldquandt,
    acorr_breusch_godfrey, acorr_ljungbox, linear_rainbow, linear_harvey_collier
)
from statsmodels.stats.stattools import durbin_watson, jarque_bera
from statsmodels.stats.diagnostic import normal_ad
from statsmodels.stats.outliers_influence import variance_inflation_factor
from scipy.stats import shapiro, probplot

def check_linear_regression_assumptions(X, y, dir=None, plot=False, k_best=10):
    """
    Extended diagnostics for OLS assumptions and influence.
    Returns a dict of test results and (optionally) writes plots to 'dir'.
    """
    # Ensure DataFrame for consistent column labels
    X = pd.DataFrame(X).copy()
    Xc = sm.add_constant(X)
    model = sm.OLS(y, Xc).fit()
    resid = model.resid
    fitted = model.fittedvalues

    results = {}

    # === Linearity / specification ===
    try:
        reset = sm.stats.diagnostic.linear_reset(model, power=2, use_f=True)
        results["RESET_F"], results["RESET_p"] = float(reset.fvalue), float(reset.pvalue)
    except Exception:
        results["RESET_F"], results["RESET_p"] = np.nan, np.nan

    try:
        rb = linear_rainbow(model)
        results["Rainbow_F"], results["Rainbow_p"] = float(rb[0]), float(rb[1])
    except Exception:
        results["Rainbow_F"], results["Rainbow_p"] = np.nan, np.nan

    try:
        hc = linear_harvey_collier(model)
        results["HarveyCollier_t"], results["HarveyCollier_p"] = float(hc[0]), float(hc[1])
    except Exception:
        results["HarveyCollier_t"], results["HarveyCollier_p"] = np.nan, np.nan

    # === Normality ===
    try:
        jb_stat, jb_p, _, _ = jarque_bera(resid)
        results["JarqueBera_stat"], results["JarqueBera_p"] = float(jb_stat), float(jb_p)
    except Exception:
        results["JarqueBera_stat"], results["JarqueBera_p"] = np.nan, np.nan

    try:
        ad_stat, ad_p = normal_ad(resid)
        results["AndersonDarling_stat"], results["AndersonDarling_p"] = float(ad_stat), float(ad_p)
    except Exception:
        results["AndersonDarling_stat"], results["AndersonDarling_p"] = np.nan, np.nan

    try:
        sh_stat, sh_p = shapiro(resid)
        results["Shapiro_stat"], results["Shapiro_p"] = float(sh_stat), float(sh_p)
    except Exception:
        results["Shapiro_stat"], results["Shapiro_p"] = np.nan, np.nan

    # === Homoscedasticity ===
    try:
        bp = het_breuschpagan(resid, model.model.exog)
        # (Lagrange multiplier stat, p-value, f-value, f p-value)
        results["BP_LM"], results["BP_p"], results["BP_F"], results["BP_F_p"] = map(float, bp)
    except Exception:
        results["BP_LM"] = results["BP_p"] = results["BP_F"] = results["BP_F_p"] = np.nan

    try:
        w = het_white(resid, model.model.exog)
        results["White_LM"], results["White_p"], results["White_F"], results["White_F_p"] = map(float, w)
    except Exception:
        results["White_LM"] = results["White_p"] = results["White_F"] = results["White_F_p"] = np.nan

    try:
        gq = het_goldfeldquandt(y, model.model.exog, alternative="two-sided")
        # returns (F, pvalue, split)
        results["GoldfeldQuandt_F"], results["GoldfeldQuandt_p"], results["GoldfeldQuandt_split"] = float(gq[0]), float(gq[1]), int(gq[2])
    except Exception:
        results["GoldfeldQuandt_F"] = results["GoldfeldQuandt_p"] = np.nan
        results["GoldfeldQuandt_split"] = None

    # === Autocorrelation ===
    try:
        results["DurbinWatson"] = float(durbin_watson(resid))
    except Exception:
        results["DurbinWatson"] = np.nan

    try:
        bg = acorr_breusch_godfrey(model, nlags=1)  # increase nlags if needed
        # (LM stat, p-value, F stat, F p-value)
        results["BG_LM"], results["BG_p"], results["BG_F"], results["BG_F_p"] = map(float, bg)
    except Exception:
        results["BG_LM"] = results["BG_p"] = results["BG_F"] = results["BG_F_p"] = np.nan

    try:
        lb = acorr_ljungbox(resid, lags=[5], return_df=True)
        results["LjungBox_Q(5)"] = float(lb["lb_stat"].iloc[0])
        results["LjungBox_p(5)"] = float(lb["lb_pvalue"].iloc[0])
    except Exception:
        results["LjungBox_Q(5)"] = results["LjungBox_p(5)"] = np.nan

    try:
        XtX = Xc.values.T @ Xc.values
        results["ConditionNumber"] = float(np.linalg.cond(XtX))
    except Exception:
        results["ConditionNumber"] = np.nan

    # === Influence & outliers ===
    try:
        infl = model.get_influence()
        summ = infl.summary_frame()
        # Cook's D, leverage (hat_diag), studentized residuals (student_resid), DFFITS, DFBETAS...
        influence_df = summ[[
            "cooks_d", "hat_diag", "student_resid", "dffits_internal", "cov_ratio"
        ]].copy()
        # DFBETAS as array -> add max |DFBETA| per obs for quick flagging
        dfbetas = infl.dfbetas
        influence_df["max_abs_dfbeta"] = np.abs(dfbetas).max(axis=1)

        # Heuristics for flags
        n, p = Xc.shape[0], Xc.shape[1] - 1
        cooks_thresh = 4 / (n - p - 1) if n > p + 1 else np.inf
        leverage_thresh = 2 * (p + 1) / n if n else np.inf

        influence_df["flag_cooks"] = influence_df["cooks_d"] > cooks_thresh
        influence_df["flag_leverage"] = influence_df["hat_diag"] > leverage_thresh
        influence_df["flag_student_resid"] = np.abs(influence_df["student_resid"]) > 3
        influence_df["flag_dffits"] = np.abs(influence_df["dffits_internal"]) > 2*np.sqrt((p+1)/n) if n else False
        influence_df["flag_dfbetas"] = influence_df["max_abs_dfbeta"] > 2/np.sqrt(n) if n else False

        results["influence_table"] = influence_df
        results["influential_indices"] = influence_df[
            influence_df.filter(like="flag_").any(axis=1)
        ].index.tolist()
    except Exception:
        results["influence_table"] = None
        results["influential_indices"] = []

    # === Plots (optional) ===
    if plot:
        if dir and not os.path.isdir(dir):
            os.makedirs(dir, exist_ok=True)

        # Residuals vs fitted (scale-location style)
        plt.figure()
        plt.scatter(fitted, resid, alpha=0.7)
        plt.axhline(0, ls="--")
        plt.xlabel("Fitted values")
        plt.ylabel("Residuals")
        plt.title("Residuals vs Fitted")
        if dir: plt.savefig(os.path.join(dir, "resid_vs_fitted.png"), dpi=200)
        # plt.show()

        # Q-Q plot
        plt.figure()
        probplot(resid, dist="norm", plot=plt)
        plt.title("Q-Q plot of residuals")
        if dir: plt.savefig(os.path.join(dir, "qq_plot_residuals.png"), dpi=200)
        # plt.show()

        # Scale-location plot (|studentized residuals| vs fitted)
        try:
            infl = model.get_influence()
            stud = infl.resid_studentized_internal
            plt.figure()
            plt.scatter(fitted, np.sqrt(np.abs(stud)), alpha=0.7)
            plt.xlabel("Fitted values")
            plt.ylabel("sqrt(|Studentized residuals|)")
            plt.title("Scale-Location")
            if dir: plt.savefig(os.path.join(dir, "scale_location.png"), dpi=200)
            # plt.show()
        except Exception:
            pass

        # Cook’s distance stem plot
        try:
            cooks = results["influence_table"]["cooks_d"].values
            plt.figure()
            markerline, stemlines, baseline = plt.stem(range(len(cooks)), cooks, use_line_collection=True)
            plt.setp(markerline, markersize=4)
            plt.xlabel("Observation")
            plt.ylabel("Cook's distance")
            plt.title("Influence: Cook's Distance")
            if dir: plt.savefig(os.path.join(dir, "cooks_distance.png"), dpi=200)
            # plt.show()
        except Exception:
            pass

    # === Robust SEs (optional summary) ===
    try:
        robust = model.get_robustcov_results(cov_type="HC3")
        results["robust_summary"] = robust.summary().as_text()
    except Exception:
        results["robust_summary"] = None
    # Peek key outputs:
    print("RESET p:", results["RESET_p"])
    print("White p:", results["White_p"])
    print("BG (autocorr) p:", results["BG_p"])
    print("Influential indices:", results["influential_indices"])
    return results, model


# def check_linear_regression_assumptions(X,y,dir=None,plot=False):
#     # Load data
#     # Fit linear regression model
#     model = sm.OLS(y, sm.add_constant(X)).fit()
#     residuals = model.resid
#     predictions = model.predict(sm.add_constant(X))

#     if plot:
#         print("\n----- Independence of Errors (Durbin-Watson) -----")
#         dw_stat = durbin_watson(residuals)
#         print(f"Durbin-Watson statistic: {dw_stat:.3f}")
#         if 1.5 < dw_stat < 2.5:
#             print("✅ No autocorrelation detected.")
#         else:
#             print("⚠️ Possible autocorrelation in residuals.")

#         print("\n----- Homoscedasticity (Breusch-Pagan Test) -----")
#         bp_test = het_breuschpagan(residuals, model.model.exog)
#         p_value_bp = bp_test[1]
#         print(f"Breusch-Pagan p-value: {p_value_bp:.3f}")
#         if p_value_bp > 0.05:
#             print("✅ Homoscedasticity assumed (good).")
#         else:
#             print("⚠️ Heteroscedasticity detected (bad).")

#         print("\n----- Normality of Errors (Shapiro-Wilk Test) -----")
#         shapiro_stat, shapiro_p = shapiro(residuals)
#         print(f"Shapiro-Wilk p-value: {shapiro_p:.3f}")
#         if shapiro_p > 0.05:
#             print("✅ Residuals appear normally distributed.")
#         else:
#             print("⚠️ Residuals may not be normally distributed.")

#         print("\n----- Normality of Errors (Q-Q Plot) -----")
#         plt.figure()
#         probplot(residuals, dist="norm", plot=plt)
#         plt.title('Q-Q plot of residuals')
#         plt.show()
#         if dir:
#             plt.savefig(os.path.join(dir, 'qq_plot_residuals.png'))
#     else:
#         return