import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix, precision_score, recall_score, accuracy_score, r2_score, mean_absolute_error, mean_squared_error, f1_score
from tkinter import filedialog, messagebox
from tkinter.simpledialog import askstring
import statsmodels.api as sm
import tkinter as tk
import traceback
from matplotlib.colors import ListedColormap
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch
from adjustText import adjust_text
from tkinter import ttk
import sys 
from pathlib import Path
import os 
from IPython import get_ipython
import math
import plotly.graph_objects as go
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from typing import Optional, List
from matplotlib.colors import ListedColormap
from sklearn import metrics
from sklearn.metrics import confusion_matrix,f1_score
from sklearn.tree import DecisionTreeClassifier
import textwrap
from matplotlib.table import Table
from datetime import datetime
import warnings
warnings.filterwarnings("ignore", message="Glyph.*missing from font")
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
import matplotlib 
matplotlib.rcParams["font.family"] = "Apotos"
try:
    from .modeling import fit_and_evaluate_single_combination_regression , fit_and_evaluate_single_combination_classification
    from .modeling_utils import _normalize_combination_to_columns, check_linear_regression_assumptions
except ImportError as e:
    from modeling import fit_and_evaluate_single_combination_regression , fit_and_evaluate_single_combination_classification
    from modeling_utils import _normalize_combination_to_columns, check_linear_regression_assumptions

def show_table_window(title, df):
    """
    Creates a new Tkinter window to display the DataFrame in a Treeview widget.

    Parameters:
        title (str): The title of the window.
        df (pd.DataFrame): The DataFrame to display.
    """
    # Create a new top-level window
    window = tk.Toplevel()
    window.title(title)

    # Set window size (optional)
    window.geometry("800x600")

    # Create a frame for the Treeview and scrollbar
    frame = ttk.Frame(window)
    frame.pack(fill=tk.BOTH, expand=True)

    # Create the Treeview
    tree = ttk.Treeview(frame, columns=list(df.columns), show='headings')

    # Define headings and column properties
    for col in df.columns:
        tree.heading(col, text=col)
        # Optionally, set column width and alignment
        tree.column(col, anchor=tk.CENTER, width=100)

    # Insert data into the Treeview
    for _, row in df.iterrows():
        tree.insert('', tk.END, values=list(row))

    # Add a vertical scrollbar
    scrollbar = ttk.Scrollbar(frame, orient=tk.VERTICAL, command=tree.yview)
    tree.configure(yscroll=scrollbar.set)
    scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
    tree.pack(side=tk.LEFT, fill=tk.BOTH, expand=True)

    # Optional: Add a horizontal scrollbar
    scrollbar_x = ttk.Scrollbar(frame, orient=tk.HORIZONTAL, command=tree.xview)
    tree.configure(xscroll=scrollbar_x.set)
    scrollbar_x.pack(side=tk.BOTTOM, fill=tk.X)


def get_valid_integer(prompt, default_value):
            while True:
                value_str = askstring("Input", prompt)
                try:
                    return int(value_str)
                except (ValueError, TypeError):
                    # If the input is not valid, return the default value
                    print(f"Using default value: {default_value}")
                    return default_value
                

def build_regression_equation(formula, coefficients, r_squared):
    """
    Build a regression equation string with proper LaTeX formatting.
    """
    import re

    intercept = getattr(coefficients, "iloc", coefficients)[0]  # Intercept
    feature_coeffs = coefficients[1:]

    def escape_latex(name: str) -> str:
        """Escape special LaTeX characters in feature names."""
        name = str(name)
        name = name.replace("_", r"\_")
        name = name.replace("%", r"\%")
        name = name.replace("Δ", r"\Delta ")
        # Wrap in \mathrm{} so text stays readable
        return r"\mathrm{" + name + "}"

    safe_formula = [escape_latex(name) for name in formula]

    equation_terms = []
    for i, coef in enumerate(feature_coeffs):
        sign = "+" if coef >= 0 else "-"
        equation_terms.append(f" {sign} {abs(coef):.2f}·{safe_formula[i]}")

    # Add intercept at the end
    sign_intercept = "+" if intercept >= 0 else "-"
    equation_terms.append(f" {sign_intercept} {abs(intercept):.2f}")

    # Build LaTeX equation
    equation = (
        r"$y =" + "".join(equation_terms).strip() + r"$" + "\n" 
    )

    return equation


## might change in the future to plot confidence intervals as dotted lines calculated from the covariance matrix
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from adjustText import adjust_text

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from adjustText import adjust_text
import matplotlib.patheffects as pe  

def generate_q2_scatter_plot(
    y,
    y_pred,
    labels,
    folds_df,
    formula,
    coefficients,
    r=None,
    lower_bound=None,
    upper_bound=None,
    figsize=(4.6, 4.6),
    width_scale=1.0,
    height_scale=1.0,
    equal_aspect=True,
    fontsize=13,
    scatter_color="#0b0e0b",
    band_color="cadetblue",
    identity_color="#1f77b4",
    palette="deep",
    dpi=300,
    plot=True,
    dir=None,
    margin_frac=0.10,
    ligand_types=None,
    type_markers=None,
    color_by_labels=False,
    show_type_legend=True,
    compact=True,
    label_fraction=1.0,
    label_strategy="residual",
    random_state=42,
    leftout_mask=None,
    leftout_pred_df=None,
    leftout_label_col="Molecule",
    leftout_meas_col="Actual",
    leftout_pred_col="Predicted",
    leftout_index_col=None,
    leftout_marker="X",
    leftout_size=50,
    leftout_edgecolor="k",
    leftout_facecolor=None,
    leftout_alpha=0.98,
    marker_size=50,
    marker_edgewidth=1.0,
    label_max=10,
    label_min_abs_residual=None,
    label_fontsize=13,
    leftout_label_fontsize=13,
    show_metrics=True,
    footer_y=0.015,
    auto_expand_figure=False,
    max_adjust_iterations=500,
    show_confidence_interval=False,
    confidence_level=0.95,
    confidence_alpha=0.22,
    confidence_on="fit",
    show_identity_line=False,
    identity_lw=1.0,
    fit_line_color="black",
    fit_line_style="--",
    fit_line_width=1.5,
    save_name="q2_scatter_plot.png",
):
    import os
    import math
    import numpy as np
    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt
    import matplotlib.patheffects as pe
    from adjustText import adjust_text
    from matplotlib.lines import Line2D

    try:
        from scipy.stats import t as student_t
        has_scipy = True
    except Exception:
        has_scipy = False

    def _build_eqn(formula, coefficients, corr):
        try:
            coeffs = np.asarray(coefficients, dtype=float).ravel()
            n = len(coeffs)
        except Exception:
            return f"R = {corr:.2f}"

        if isinstance(formula, str):
            features = [formula]
        else:
            try:
                features = list(formula)
            except Exception:
                features = [f"x{i+1}" for i in range(n)]

        if n == 0:
            return f"R = {corr:.2f}"

        if n == 1:
            return rf"$\Delta\Delta G^{{\ddagger}}$ = {coeffs[0]:.3g}{features[0]}"

        if n == 2:
            a, b = coeffs
            name = features[0] if len(features) >= 1 else "x"
            sign = "+" if b >= 0 else "-"
            return rf"$\Delta\Delta G^{{\ddagger}}$ = {a:.3g}{name} {sign} {abs(b):.3g}"

        feat_names = features[: n - 1] if len(features) >= n - 1 else [f"x{i+1}" for i in range(n - 1)]
        terms = [f"{c:.3g}{f}" for c, f in zip(coeffs[:-1], feat_names)]
        intercept = coeffs[-1]
        sign = "+" if intercept >= 0 else "-"
        eqn = " + ".join(terms) + f" {sign} {abs(intercept):.3g}"
        return rf"$\Delta\Delta G^{{\ddagger}}$ = {eqn}"

    def _t_crit(conf_level, dof):
        alpha = 1.0 - conf_level
        if dof <= 0:
            return 1.96
        if has_scipy:
            return float(student_t.ppf(1.0 - alpha / 2.0, dof))
        return 1.96

    def _compute_fit_and_ci(x, y_obs, x_grid, slope, intercept, conf_level=0.95, mode="fit"):
        x = np.asarray(x, dtype=float)
        y_obs = np.asarray(y_obs, dtype=float)
        x_grid = np.asarray(x_grid, dtype=float)

        y_hat_obs = slope * x + intercept
        resid = y_obs - y_hat_obs
        n = len(x)

        if n < 3:
            y_fit = slope * x_grid + intercept
            return y_fit, None, None

        x_bar = np.mean(x)
        sxx = np.sum((x - x_bar) ** 2)

        if sxx <= 0:
            y_fit = slope * x_grid + intercept
            return y_fit, None, None

        dof = n - 2
        rss = np.sum(resid ** 2)
        mse = rss / dof if dof > 0 else 0.0
        s = np.sqrt(max(mse, 0.0))
        tcrit = _t_crit(conf_level, dof)

        y_fit = slope * x_grid + intercept
        base = (1.0 / n) + ((x_grid - x_bar) ** 2) / sxx

        mode = str(mode).lower()
        if mode in {"prediction", "pred"}:
            se = s * np.sqrt(1.0 + base)
        else:
            se = s * np.sqrt(base)

        lower = y_fit - tcrit * se
        upper = y_fit + tcrit * se
        return y_fit, lower, upper

    rng = np.random.default_rng(random_state)

    y = np.asarray(y).ravel().astype(float)
    y_pred = np.asarray(y_pred).ravel().astype(float)
    labels = np.asarray(labels).astype(str)

    if not (len(y) == len(y_pred) == len(labels)):
        raise ValueError("y, y_pred, and labels must have the same length.")

    data = pd.DataFrame({
        "Measured": y,
        "Predicted": y_pred,
        "Labels": labels
    })
    data["Residual"] = data["Predicted"] - data["Measured"]

    if ligand_types is not None:
        ligand_types = np.asarray(ligand_types).astype(str)
        if len(ligand_types) != len(data):
            raise ValueError("ligand_types must have the same length as y and y_pred.")
        data["Type"] = ligand_types
        type_levels = pd.Categorical(data["Type"]).categories.tolist()
    else:
        type_levels = []

    if ligand_types is not None:
        if type_markers is None:
            default_cycle = ["o", "s", "D", "^", "v", "P", "X", "*", "h", "<", ">"]
            if len(type_levels) > len(default_cycle):
                k = int(np.ceil(len(type_levels) / len(default_cycle)))
                default_cycle = (default_cycle * k)[:len(type_levels)]
            type_markers = {t: m for t, m in zip(type_levels, default_cycle)}
        else:
            missing = [t for t in type_levels if t not in type_markers]
            if missing:
                raise ValueError(f"Missing marker shapes for types: {missing}")

    sns.set_theme(style="white", palette=palette)

    n_labels_total = len(data)
    longest_label = max((len(str(x)) for x in labels), default=8)

    fw, fh = figsize
    fw *= width_scale
    fh *= height_scale

    if auto_expand_figure:
        fw += min(8.0, 0.12 * n_labels_total + 0.018 * longest_label * max(n_labels_total, 1))
        fh += min(8.0, 0.10 * n_labels_total + 0.012 * longest_label * max(n_labels_total, 1))
        fw = max(fw, 10.0)
        fh = max(fh, 10.0)

    fig, ax = plt.subplots(figsize=(fw, fh), dpi=dpi)

    scatter_kwargs = dict(
        x="Measured",
        y="Predicted",
        s=marker_size,
        edgecolor="white",
        linewidth=marker_edgewidth,
        zorder=3,
        ax=ax
    )

    if ligand_types is not None and color_by_labels:
        sns.scatterplot(
            data=data,
            hue="Type",
            style="Type",
            markers=type_markers,
            palette=palette,
            **scatter_kwargs
        )
        show_type_legend = True
    else:
        sns.scatterplot(
            data=data,
            color=scatter_color,
            legend=False,
            **scatter_kwargs
        )
        show_type_legend = False

    lo_df = None
    texts_lo_all = []

    if leftout_pred_df is not None and len(leftout_pred_df):
        lo_df = leftout_pred_df.copy()
        y_meas = pd.to_numeric(lo_df.get(leftout_meas_col, lo_df.get("Measured")), errors="coerce")
        y_pr = pd.to_numeric(lo_df.get(leftout_pred_col, lo_df.get("Predicted")), errors="coerce")
        lo_df = lo_df.assign(Measured=y_meas, Predicted=y_pr).dropna(subset=["Measured", "Predicted"])

        ax.scatter(
            lo_df["Measured"],
            lo_df["Predicted"],
            marker=leftout_marker,
            s=leftout_size,
            facecolors=leftout_facecolor or "blue",
            edgecolors=leftout_edgecolor,
            linewidths=1.2,
            alpha=leftout_alpha,
            zorder=5,
            label="External Validation"
        )

    try:
        corr = float(r) if r is not None else float(np.corrcoef(y, y_pred)[0, 1])
    except Exception:
        corr = np.nan

    eqn = _build_eqn(formula, coefficients, corr)

    lo = float(np.nanmin([np.nanmin(y), np.nanmin(y_pred)]))
    hi = float(np.nanmax([np.nanmax(y), np.nanmax(y_pred)]))

    if lo_df is not None and not lo_df.empty:
        lo = min(lo, float(np.nanmin(lo_df[["Measured", "Predicted"]].to_numpy())))
        hi = max(hi, float(np.nanmax(lo_df[["Measured", "Predicted"]].to_numpy())))

    if lower_bound is not None:
        lo = min(lo, float(lower_bound))
    if upper_bound is not None:
        hi = max(hi, float(upper_bound))

    span = hi - lo
    pad = margin_frac * span if span > 0 else 1.0
    lo_p, hi_p = lo - pad, hi + pad

    xmin = np.floor(lo_p * 2) / 2
    xmax = np.ceil(hi_p * 2) / 2

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(xmin, xmax)

    # identical ticks
    ticks = np.arange(xmin, xmax + 0.001, 0.5)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)

    # identical formatting (same digits)
    from matplotlib.ticker import FormatStrFormatter
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

    if equal_aspect:
        ax.set_aspect("equal", adjustable="box")

    coeffs_flat = None
    try:
        if coefficients is not None:
            coeffs_flat = np.asarray(coefficients, dtype=float).ravel()
    except Exception:
        coeffs_flat = None

    try:
        if coeffs_flat is not None and len(coeffs_flat) == 2:
            slope, intercept = float(coeffs_flat[0]), float(coeffs_flat[1])
        else:
            slope, intercept = np.polyfit(y, y_pred, 1)
    except Exception:
        slope, intercept = np.polyfit(y, y_pred, 1)

    if show_identity_line:
        xx_id = np.linspace(*ax.get_xlim(), 200)
        ax.plot(
            xx_id,
            xx_id,
            color=identity_color,
            lw=identity_lw,
            linestyle="-",
            alpha=0.8,
            zorder=1
        )

    xx = np.linspace(*ax.get_xlim(), 300)

    if show_confidence_interval:
        ci_mode = "fit" if str(confidence_on).lower() not in {"prediction", "pred"} else "prediction"
        yy_fit, yy_lo, yy_hi = _compute_fit_and_ci(
            x=y,
            y_obs=y_pred,
            x_grid=xx,
            slope=slope,
            intercept=intercept,
            conf_level=confidence_level,
            mode=ci_mode,
        )
        if yy_lo is not None and yy_hi is not None:
            ax.fill_between(
                xx,
                yy_lo,
                yy_hi,
                color=band_color,
                alpha=confidence_alpha,
                zorder=1.5,
                linewidth=0
            )
    else:
        yy_fit = slope * xx + intercept

    ax.plot(
        xx,
        yy_fit,
        color=fit_line_color,
        lw=fit_line_width,
        linestyle=fit_line_style,
        zorder=2
    )

    n_total = len(data)

    if label_strategy == "residual":
        order = np.argsort(np.abs(data["Residual"].values))[::-1]
    else:
        order = rng.permutation(n_total)

    if label_min_abs_residual is not None:
        order = [i for i in order if abs(data["Residual"].iat[i]) >= label_min_abs_residual]

    if label_max is None:
        idx = np.array(order, dtype=int)
    else:
        n_default = max(1, int(np.ceil(label_fraction * n_total)))
        n_to_label = min(label_max if label_max is not None else n_default, n_total)
        idx = np.array(order[:n_to_label], dtype=int)

    if label_max is None and label_fraction >= 1.0:
        idx = np.arange(n_total, dtype=int)

    texts = []
    initial_dx = 0.03  * (hi_p - lo_p)
    initial_dy = 0.03  * (hi_p - lo_p)

    for k, (_, row) in enumerate(data.iloc[idx].iterrows()):
        x0 = float(row["Measured"])
        y0 = float(row["Predicted"])

        angle = 2.0 * math.pi * (k / max(len(idx), 1))
        x_text = x0 + 1.3 * initial_dx * math.cos(angle)
        y_text = y0 + 1.3 * initial_dy * math.sin(angle)

        t = ax.text(
            x_text,
            y_text,
            str(row["Labels"]),
            fontsize=label_fontsize,
            ha="center",
            va="center",
            color="black",
            fontweight="normal",
            clip_on=False,
            zorder=6,
            path_effects=[
                pe.withStroke(linewidth=2.8, foreground="white", alpha=1.0)
            ]
        )
        texts.append(t)

    if lo_df is not None and not lo_df.empty:
        for k, (_, row) in enumerate(lo_df.iterrows()):
            x0 = float(row["Measured"])
            y0 = float(row["Predicted"])

            angle = 2.0 * math.pi * (k / max(len(lo_df), 1))
            x_text = x0 + 1.2 * initial_dx * math.cos(angle)
            y_text = y0 + 1.2 * initial_dy * math.sin(angle)

            t = ax.text(
                x_text,
                y_text,
                str(row.get(leftout_label_col, "")),
                fontsize=leftout_label_fontsize,
                ha="center",
                va="center",
                color="blue",
                fontweight="normal",
                clip_on=False,
                zorder=7,
                path_effects=[
                    pe.withStroke(linewidth=2.8, foreground="white", alpha=1.0)
                ]
            )
            texts_lo_all.append(t)

    all_x = data["Measured"].to_numpy()
    all_y = data["Predicted"].to_numpy()
    if lo_df is not None and not lo_df.empty:
        all_x = np.concatenate([all_x, lo_df["Measured"].to_numpy()])
        all_y = np.concatenate([all_y, lo_df["Predicted"].to_numpy()])

    try:
        adjust_text(
            texts + texts_lo_all,
            x=all_x,
            y=all_y,
            ax=ax,
            expand_text=(1.25, 1.35),
            expand_points=(1.15, 1.25),
            force_text=(1.2, 1.5),
            force_points=(0.8, 1.0),
            force_pull=(0.08, 0.10),
            only_move={"points": "y", "text": "xy"},
            arrowprops=dict(
                arrowstyle="-",
                color="gray",
                lw=0.5,
                alpha=0.7
            ),
            ensure_inside_axes=False,
            expand_axes=True,
            lim=max_adjust_iterations
        )
    except Exception:
        pass

    ax.margins(x=0.18, y=0.18)

    ax.set_xlabel(
    r"Measured $\Delta\Delta G^{\ddagger}$",
    fontsize=fontsize + 1,
    fontweight="normal"
    )

    ax.set_ylabel(
        r"Predicted $\Delta\Delta G^{\ddagger}$",
        fontsize=fontsize + 1,
        fontweight="normal"
    )
    # ax.set_title("Predicted vs Measured", fontsize=fontsize + 3, fontweight="normal")

    if compact:
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.tick_params(
            axis="both",
            which="both",
            length=4,
            width=1.0,
            labelsize=max(9, fontsize - 1)
        )

    handles = [
        Line2D(
            [0], [0],
            marker="o",
            color="w",
            markerfacecolor=scatter_color,
            markeredgecolor="w",
            markeredgewidth=marker_edgewidth,
            markersize=np.sqrt(marker_size),
            label="Training Set"
        )
    ]
    labels_legend = ["Training Set"]

    if leftout_pred_df is not None and len(leftout_pred_df):
        handles.append(
            Line2D(
                [0], [0],
                marker=leftout_marker,
                color="w",
                markerfacecolor=leftout_facecolor or "blue",
                markeredgecolor=leftout_edgecolor,
                markeredgewidth=1.0,
                markersize=np.sqrt(leftout_size),
                label="External Validation"
            )
        )
        labels_legend.append("External Validation")

    if show_confidence_interval:
        ci_text = f"{int(round(confidence_level * 100))}% CI"
        if str(confidence_on).lower() in {"prediction", "pred"}:
            ci_text = f"{int(round(confidence_level * 100))}% Prediction Interval"

        handles.append(
            Line2D(
                [0], [0],
                color=band_color,
                lw=8,
                alpha=confidence_alpha,
                label=ci_text
            )
        )
        labels_legend.append(ci_text)

    ax.legend(
        handles,
        labels_legend,
        loc="lower right",
        frameon=True,
        facecolor="white",
        edgecolor="black",
        framealpha=0.95,
        fontsize=max(9, fontsize),
        markerscale=1.0,
        handlelength=1.3,
        handletextpad=0.6,
        borderpad=0.5,
        labelspacing=0.35
    )

    if show_metrics:
        if folds_df is not None and not getattr(folds_df, "empty", True):
            q = folds_df.iloc[0]
            q_txt = (
            f"R$^2$ = {corr:.2f}"
            f"\n3-fold Q$^2$: {q.get('Q2_3_Fold', np.nan):.2f}"
            f"\n5-fold Q$^2$: {q.get('Q2_5_Fold', np.nan):.2f}"
            f"\nLOOCV Q$^2$: {q.get('Q2_LOOCV', np.nan):.2f}"
        )
        else:
            q_txt = f"R$^2$ = {corr:.2f}"

        ax.text(
            0.02,
            0.98,
            q_txt,
            transform=ax.transAxes,
            fontsize=fontsize,
            fontweight="normal",
            va="top",
            ha="left",
            bbox=dict(
                facecolor="white",
                edgecolor="black",
                alpha=0.97,
                boxstyle="round,pad=0.35"
            ),
            zorder=8
        )

    fig.subplots_adjust(
        left=0.12,
        right=0.96,
        bottom=max(0.2, footer_y + 0.07),
        top=0.93
    )

    fig.text(
        0.5,
        footer_y,
        eqn,
        ha="center",
        va="center",
        fontsize=max(9, fontsize - 1),
        fontweight="normal",
        bbox=dict(
            facecolor="white",
            edgecolor="black",
            alpha=0.97,
            boxstyle="round,pad=0.35"
        )
    )

    if dir is not None:
        os.makedirs(dir, exist_ok=True)
        outpath = os.path.join(dir, save_name)
        fig.savefig(outpath, dpi=dpi, bbox_inches="tight")

    if plot:
        plt.show()

    return fig, ax





def plot_probabilities(probabilities_df, sample_names):
    df = probabilities_df.copy()
    
    # Rename if needed
    if 'prediction' in df.columns:
        df.rename(columns={'prediction': 'Predicted_Class'}, inplace=True)
    if 'True_Class' in df.columns:
        df.rename(columns={'True_Class': 'Actual_Class'}, inplace=True)
    
    # Identify probability columns by prefix
    prob_cols = [col for col in df.columns if col.startswith('Prob_Class_')]
    if not prob_cols:
        raise ValueError("No probability columns found (expecting columns like 'Prob_Class_1', ...)")
    
    # Ensure classes are integers for building column names
    df['Actual_Class'] = df['Actual_Class'].astype(int)
    df['Predicted_Class'] = df['Predicted_Class'].astype(str)
    
    # Compute rankings across the probability columns
    rankings = df[prob_cols].rank(axis=1, ascending=False, method='min')
    
    # Helper to get the rank of the true class
    def _true_rank(row):
        prob_col = f"Prob_Class_{row['Actual_Class']}"
        if prob_col not in rankings.columns:
            raise KeyError(f"Expected column {prob_col} in probabilities, got {prob_cols}")
        return int(rankings.at[row.name, prob_col])
    
    df['Rank'] = df.apply(_true_rank, axis=1)
    
    # Color‐code by rank (1=green, 2=yellow, 3=red, >3=gray)
    color_map = {1: 'green', 2: 'yellow', 3: 'red'}
    df['Color_Code'] = df['Rank'].map(color_map).fillna('gray')
    
    # Build sample labels from your list
    if len(sample_names) != len(df):
        raise ValueError("sample_names length must match number of rows in probabilities_df")
    df['Label'] = [
        f"{name} (Pred: {pred}, Actual: {act})"
        for name, pred, act in zip(
            sample_names,
            df['Predicted_Class'].astype(str),
            df['Actual_Class'].astype(str)
        )
    ]
    
    # Plot heatmap
    plt.figure(figsize=(12, 10))
    sns.heatmap(
        df[prob_cols].astype(float),
        cmap='Blues',
        annot=True,
        fmt=".2f",
        cbar_kws={'label': 'Probability'}
    )

    ax = plt.gca()
    # Y‐ticks: use your sample labels
    ax.set_yticks(np.arange(0.5, len(df), 1))
    ax.set_yticklabels(df['Label'], rotation=0, fontsize=10)
    for label, color in zip(ax.get_yticklabels(), df['Color_Code']):
        label.set_color(color)
    
    # X‐ticks: strip the prefix for readability
    classes = [col.replace('Prob_Class_', '') for col in prob_cols]
    ax.set_xticklabels(classes, rotation=45, ha='right')
    
    plt.title('Probability Heatmap with Prediction vs. Actual')
    plt.xlabel('Class')
    plt.ylabel('Samples')
    plt.tight_layout()
    plt.show()




def plot_enhanced_confusion_matrix(cm, classes, precision, recall, accuracy, figsize=(10, 10), annot_fontsize=12):
    """
    Plot an enhanced confusion matrix with precision, recall, and accuracy.
    
    The confusion matrix will include an extra row and column for total samples and precision.
    
    Parameters:
    - cm: Confusion matrix (as a 2D array)
    - classes: List of class names
    - precision: List of precision values for each class
    - recall: List of recall values for each class
    - accuracy: Overall accuracy value
    - figsize: Size of the figure
    - annot_fontsize: Font size of the annotations inside the matrix
    """
    # Number of classes
    n = len(classes)

    # Expand the confusion matrix to 4x4
    cm_expanded = np.zeros((n + 1, n + 1))

    # Fill the confusion matrix values into the expanded matrix
    cm_expanded[:n, :n] = cm

    # Calculate row totals (true class totals) and column totals (predicted class totals)
    row_totals = cm.sum(axis=1)  # Row totals for the right column
    col_totals = cm.sum(axis=0)  # Column totals for the bottom row

    # Fill the last row with column totals and precision percentages
    cm_expanded[n, :n] = col_totals

    # Fill the last column with row totals and percentages of the total dataset
    cm_expanded[:n, n] = row_totals

    # Add overall total and accuracy in the bottom-right corner
    cm_expanded[n, n] = cm.sum()

    # Create a DataFrame for plotting
    classes_with_total = classes + ['Total']
    cm_df = pd.DataFrame(cm_expanded, index=classes_with_total, columns=classes_with_total)

    # Create annotations for the confusion matrix cells
    annot = np.empty_like(cm_expanded).astype(str)

    # Fill in the confusion matrix cells with counts and percentages (diagonal and off-diagonal)
    for i in range(n):
        for j in range(n):
            annot[i, j] = f"{int(cm_expanded[i, j])}\n({cm_expanded[i, j] / row_totals[i] * 100:.1f}%)" if row_totals[i] != 0 else f"{int(cm_expanded[i, j])}\n(0%)"

    # Fill the last row with precision percentages
    for i in range(n):
        annot[n, i] = f"{int(col_totals[i])}\n({precision[i] * 100:.1f}%)"

    # Fill the last column with row totals and percentage of total dataset
    for i in range(n):
        annot[i, n] = f"{int(row_totals[i])}\n({row_totals[i] / cm.sum() * 100:.1f}%)"

    # Add the total and accuracy in the bottom-right corner
    annot[n, n] = f"Total: {int(cm.sum())}\nAcc: {accuracy * 100:.1f}%"

    # Plot the expanded confusion matrix
    plt.figure(figsize=figsize)

    # Create a mask for diagonal cells (true positives)
    mask_diag = np.zeros_like(cm_expanded, dtype=bool)
    np.fill_diagonal(mask_diag, True)

    # Create a mask for the bottom-right cell (total and accuracy)
    mask_acc = np.zeros_like(cm_expanded, dtype=bool)
    mask_acc[-1, -1] = True

    # Mask for off-diagonal cells
    mask_false = np.ones_like(cm_expanded, dtype=bool)
    np.fill_diagonal(mask_false, False)
    mask_false[-1, :] = False  # Do not mask the totals row
    mask_false[:, -1] = False  # Do not mask the totals column

    # Mask for total row (yellow gradient)
    mask_total_row = np.zeros_like(cm_expanded, dtype=bool)
    mask_total_row[-1, :-1] = True  # The entire last row except the bottom-right corner

    # Mask for total column (grey gradient)
    mask_total_column = np.zeros_like(cm_expanded, dtype=bool)
    mask_total_column[:-1, -1] = True  # The entire last column except the bottom-right corner

    # Plot true positive cells (diagonal) using green gradient (Greens cmap)
    sns.heatmap(cm_df, annot=annot, annot_kws={"size": annot_fontsize}, fmt='', mask=~mask_diag, cmap="Greens", vmin=0, vmax=np.max(cm), cbar=False, square=True, linewidths=1, linecolor='black')

    # Plot false predictions (off-diagonal) using red gradient (Reds cmap)
    sns.heatmap(cm_df, annot=annot, annot_kws={"size": annot_fontsize}, fmt='', mask=~mask_false, cmap="Reds", vmin=0, vmax=np.max(cm), cbar=False, square=True, linewidths=1, linecolor='black')

    # Plot total row (yellow gradient)
    sns.heatmap(cm_df, annot=annot, annot_kws={"size": annot_fontsize}, fmt='', mask=~mask_total_row, cmap="YlOrBr", vmin=0, vmax=np.max(cm), cbar=False, square=True, linewidths=1, linecolor='black')

    # Plot total column (grey gradient)
    sns.heatmap(cm_df, annot=annot, annot_kws={"size": annot_fontsize}, fmt='', mask=~mask_total_column, cmap="Greys", vmin=0, vmax=np.max(cm), cbar=False, square=True, linewidths=1, linecolor='black')

    # Plot total samples and precision/accuracy using blue gradient (Blues cmap)
    sns.heatmap(cm_df, annot=annot, annot_kws={"size": annot_fontsize}, fmt='', mask=~mask_acc, cmap="Blues", vmin=0, vmax=np.max(cm), cbar=False, square=True, linewidths=1, linecolor='black')

    # Add a color legend for each section
    
    legend_elements = [Patch(facecolor='green', edgecolor='black', label='True Positives (Diagonal)'),
                       Patch(facecolor='red', edgecolor='black', label='False Positives/Negatives'),
                       Patch(facecolor='yellow', edgecolor='black', label='Precision'),
                       Patch(facecolor='grey', edgecolor='black', label='Total Column (Grey)'),
                       Patch(facecolor='blue', edgecolor='black', label='Overall Total and Accuracy')]
    
    plt.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1.35, 1), title='Legend')

    plt.title('Enhanced Confusion Matrix with Totals, Precision, and Accuracy', fontsize=16)
    plt.xlabel('Predicted Class')
    plt.ylabel('True Class')
    plt.tight_layout()
    plt.show()





def print_models_classification_table(results , app=None , model=None):
    results_tmp = results.reset_index(drop=True).copy()
    formulas=results_tmp['combination']
    accuracy=results_tmp['accuracy']
    precision=results_tmp['precision']
    recall=results_tmp['recall']
    f1=results_tmp['f1_score']
    mcfaden=results_tmp['mcfadden_r2']
    model_ids=[i for i in range(len(results_tmp))]
    avg_accuracy=results_tmp['avg_accuracy']
    avg_f1=results_tmp['avg_f1_score']
    avg_mcfaden=results_tmp['avg_mcfadden_r2']

    df = pd.DataFrame({
        'formula': formulas,
        'accuracy': accuracy,
        'precision': precision,
        'recall': recall,
        'f1': f1,
        'mcfadden_r2': mcfaden,
        'avg_mcfadden_r2': avg_mcfaden,
        'avg_accuracy': avg_accuracy,
        'avg_f1': avg_f1,
        'Model_id': model_ids
    })
    df.sort_values(by='avg_accuracy', ascending=False, inplace=True)
    # Set the index to range from 1 to n (1-based indexing)
    df.index = range(1, len(df) + 1)
    if app:
        app.show_result(df.to_markdown(index=False, tablefmt="pipe"))
        messagebox.showinfo('Models List: ', df.to_markdown(index=False, tablefmt="pipe"))  
        
    else:
        print(df.to_markdown(index=False, tablefmt="pipe"))
        

    try:
        df.to_csv('models_classification_table.csv', index=False)
    except:
        print('could not save the table')

    
    while True:
        if app:
            selected_model = get_valid_integer('Select a model number: default is 0', 0)
        else:
            try:
                selected_model = int(input("Select a model number (or -1 to exit): "))
            except ValueError:
                print("Invalid input. Please enter a number.")
                continue
            
        if selected_model == -1:
            print("Exiting model selection.")
            break


        _, probablities_df = fit_and_evaluate_single_combination_classification(model, formulas[selected_model], return_probabilities=True)
        colms=_normalize_combination_to_columns(formulas[selected_model])
        X=model.features_df[colms]
        # x=pd.DataFrame(X, columns=formulas[selected_model])
        vif_df = model._compute_vif(X)
        samples_names=model.molecule_names
        plot_probabilities(probablities_df, samples_names)
        print_models_vif_table(vif_df)
        # Print the confusion matrix
        y_pred = model.predict(model.features_df[colms].to_numpy())
        y_true = model.target_vector.to_numpy()
        cm = confusion_matrix(y_true, y_pred)
        print("\nConfusion Matrix\n")
        class_names = np.unique(y_true)
        class_names = [f'Class_{i}' for i in class_names]

        model_precision = precision_score(y_true, y_pred, average=None)
        model_recall = recall_score(y_true, y_pred, average=None)
        model_accuracy = accuracy_score(y_true, y_pred)

        # Save text file with the results
        with open('classification_results.txt', 'a') as f:
            f.write(f"Models List\n\n{df.to_markdown(index=False, tablefmt='pipe')}\n\n Probabilities\n\n{probablities_df.to_markdown(tablefmt='pipe')}\n\nConfusion Matrix\n\n{cm}\n\nPrecision\n\n{model_precision}\n\nRecall\n\n{model_recall}\n\nAccuracy\n\n{model_accuracy}\n\n")
            print('Results saved to classification_results.txt in {}'.format(os.getcwd()))

        plot_enhanced_confusion_matrix(cm, class_names, model_precision, model_recall, model_accuracy)

        # Ask the user if they want to select another model or exit
        if not app:
            cont = input("Do you want to select another model? (y/n): ").strip().lower()
            if cont != 'y':
                print("Exiting model selection.")
                break
        else:
            cont=messagebox.askyesno('Continue','Do you want to select another model?')
            if not cont:
                break



def print_models_vif_table(results, app=None):
    if app:
        app.show_result('\n\n\n')
        app.show_result('VIF Table\n')
        app.show_result('---\n')
        app.show_result(results.to_markdown(index=False, tablefmt="pipe"))
    else:
        print('\n\n\n')
        print('VIF Table\n')
        print('---\n')
        print(results.to_markdown(index=False, tablefmt="pipe"))

def _capture_new_figs(run_fn):
    """
    Run a plotting function and return the list of NEW matplotlib Figure objects it created.
    Works by diffing figure numbers before/after.
    """
    before = set(plt.get_fignums())
    run_fn()
    after = set(plt.get_fignums())
    new_nums = sorted(after - before)
    return [plt.figure(num) for num in new_nums]

# ---- helpers ---------------------------------------------------------------
from matplotlib.gridspec import GridSpec

def _nice_table(ax, df, title=None, fontsize=9):
    ax.axis('off')
    if title:
        ax.set_title(title, pad=8, fontsize=11, fontweight='normal')
    tbl = ax.table(cellText=df.values, colLabels=df.columns, loc='center', cellLoc='center')
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(fontsize)
    # Adaptive column widths
    ncols = len(df.columns)
    for (row, col), cell in tbl.get_celld().items():
        # header row
        if row == 0:
            cell.set_text_props(weight='normal')
            cell.set_height(0.08)
        else:
            cell.set_height(0.06)
        cell.set_edgecolor('#DDDDDD')
        if col < ncols:
            cell.set_linewidth(0.6)
    # squeeze to panel
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)

def _page_header(fig, subtitle):
    fig.suptitle(subtitle, y=0.98, fontsize=13, fontweight='normal')
    

def _normalize_combination_to_columns(combination) -> list[str]:
    """
    Robustly convert a 'combination' (tuple/list/string) into a list of feature names.
    Handles tuple-like strings with internal commas such as "('Dist(1, 2)', 'dip_x')".
    """
    if isinstance(combination, (list, tuple, pd.Index, np.ndarray)):
        return [str(x).strip().strip("'\"") for x in combination]

    if not isinstance(combination, str):
        s = str(combination)
    else:
        s = combination
    s = s.strip()

    # Remove outer parentheses if tuple-like
    if s.startswith("(") and s.endswith(")"):
        s = s[1:-1].strip()

    # Split on top-level commas only (ignore commas inside (), [], {}, or quotes)
    items, buf, depth = [], [], 0
    in_sq = in_dq = False

    def flush():
        token = "".join(buf).strip().strip("'\"")
        if token:
            items.append(token)

    for ch in s:
        if ch == "'" and not in_dq:
            in_sq = not in_sq
            buf.append(ch); continue
        if ch == '"' and not in_sq:
            in_dq = not in_dq
            buf.append(ch); continue

        if not in_sq and not in_dq:
            if ch in "([{": depth += 1
            elif ch in ")]}": depth -= 1
            elif ch == "," and depth == 0:
                flush(); buf = []
                continue
        buf.append(ch)
    flush()

    # Single bare name case
    return items or ([s.strip().strip("'\"")] if s else [])

def _stringify_formula(features: list[str]) -> str:
    """Pretty print a features list as 'f1 + f2 + ...'."""
    return " + ".join(map(str, features))

def _footer(fig, left_text: str = "", page_num: Optional[int] = None):
    """Add a tiny footer with optional page number."""
    if page_num is not None:
        fig.text(0.99, 0.01, f"p. {page_num}", ha="right", va="bottom", fontsize=8, alpha=0.6)
    if left_text:
        fig.text(0.01, 0.01, left_text, ha="left", va="bottom", fontsize=8, alpha=0.6)


def _safe_table(
    ax,
    df: "pd.DataFrame",
    title: str,
    *,
    max_rows: int = 30,
    max_cols: int = 12,
    wrap_width: int = 36,          # wrap long strings
    floatfmt: str = "{:.4g}",
    start_fontsize: int = 11,      # larger default
    min_fontsize: int = 8,         # don't go below readable
    cell_scale: float = 1.0,       # no pre-shrink
    max_shrink_steps: int = 20,    # softer shrinking
    show_index: bool = False,
    zebra: bool = True,
    header_facecolor: str = "#ECECEC",
    header_text_weight: str = "normal",
    header_text_align: str = "center",
    # NEW: sizing strategy
    auto_resize_figure: bool = True,
    min_figsize: Optional[tuple[float, float]] = None,        # (w,h) inches; None -> current fig size
    max_figsize: Optional[tuple[float, float]] = (8.2, 10.5), # A4-ish portrait area inside margins
    target_col_total: float = 0.98,          # use (almost) full axes width
    debug: bool = False,
):
    """
    Render df as a readable table inside an axes. Prefers readability:
    1) limit rows/cols, 2) wrap text, 3) grow figure up to max_figsize,
    4) only then shrink font/scale a bit. Falls back to monospace if still too big.
    """
    ax.set_title(title, fontsize=12, fontweight="normal", pad=8)
    ax.axis("off")

    # ---- Early exits
    if df is None:
        ax.text(0.5, 0.5, "No data (None)", ha="center", va="center", fontsize=10, alpha=0.6)
        return
    if not isinstance(df, pd.DataFrame):
        ax.text(0.5, 0.5, f"Unsupported type: {type(df)}", ha="center", va="center", fontsize=10, alpha=0.6)
        return
    if df.shape[0] == 0 or df.shape[1] == 0:
        ax.text(0.5, 0.5, "No rows/columns to display", ha="center", va="center", fontsize=10, alpha=0.6)
        return

    # ---- 0) Subset
    df_disp = df.copy()
    if not show_index:
        df_disp = df_disp.reset_index(drop=True)
    df_disp = df_disp.iloc[:max_rows, :max_cols]

    # ---- 1) Format + wrap
    def _fmt_cell(x):
        if pd.isna(x):
            return ""
        if isinstance(x, (float, np.floating, int, np.integer)):
            try:
                return floatfmt.format(x)
            except Exception:
                try:
                    return f"{float(x):.4g}"
                except Exception:
                    return str(x)
        s = str(x)
        if wrap_width and len(s) > wrap_width:
            s = textwrap.fill(s, width=wrap_width)
        return s

    formatted = df_disp.map(_fmt_cell)
    col_labels = [textwrap.fill(str(c), width=wrap_width) if wrap_width else str(c) for c in df_disp.columns]

    # ---- 2) Estimate relative column widths (downweight wrapping)
    def _max_len(col_series, col_label):
        lens = [len(str(v)) for v in col_series]
        mx = max(lens + [len(col_label), 1])
        mx = mx / math.sqrt(max(1, math.ceil(mx / max(1, wrap_width))))
        return mx

    col_lengths = np.array(
        [_max_len(formatted.iloc[:, j], col_labels[j]) for j in range(formatted.shape[1])],
        dtype=float
    )
    if not np.isfinite(col_lengths).all() or np.all(col_lengths <= 0):
        col_lengths = np.ones_like(col_lengths)
    rel_widths = col_lengths / col_lengths.sum()

    # ---- 3) Build table
    tbl = ax.table(cellText=formatted.values, colLabels=col_labels, loc="center")
    tbl.auto_set_font_size(False)
    fs = start_fontsize
    tbl.set_fontsize(fs)
    tbl.scale(cell_scale, cell_scale)

    # Header style
    for (r, c), cell in tbl.get_celld().items():
        if r == 0:
            cell.set_facecolor(header_facecolor)
            cell.set_text_props(weight=header_text_weight, ha=header_text_align, va="center")

    # Fit columns to (almost) full axes width — NOT >1 like before
    widths = (rel_widths / rel_widths.sum()) * target_col_total
    ncols = formatted.shape[1]
    for j in range(ncols):
        for (r, c), cell in tbl.get_celld().items():
            if c == j:
                cell.set_width(widths[j])

    # Align numbers on the right
    num_like = {j for j in range(ncols) if pd.api.types.is_numeric_dtype(df_disp.dtypes.iloc[j])}
    for (r, c), cell in tbl.get_celld().items():
        if r == 0:
            continue
        if c in num_like:
            cell.set_text_props(ha="right", va="center")
        else:
            cell.set_text_props(ha="left", va="center")

    # Zebra stripes on data rows
    if zebra:
        nrows = formatted.shape[0]
        for (r, c), cell in tbl.get_celld().items():
            if r > 0 and (r % 2 == 0) and r <= nrows:
                cell.set_facecolor("#F5F5F5")

    fig = ax.get_figure()

    # ---- Helper: estimate required figure size (inches) for readability
    def _estimate_needed_inches(nrows, ncols, col_lengths, fontsize):
        # 1 pt = 1/72 inch; rough width per char ≈ 0.6*fontsize pt
        char_w_in = 0.6 * fontsize / 72.0
        row_h_in = 1.35 * fontsize / 72.0
        # padding ~ 4 chars per column
        col_in = char_w_in * (col_lengths + 4)
        width = float(col_in.sum()) + 0.6  # margins
        height = row_h_in * (nrows + 1) + 0.8  # + header + padding
        return width, height

    # ---- 4) Prefer growing the figure (up to max_figsize) before shrinking
    if auto_resize_figure:
        cur_w, cur_h = fig.get_size_inches()
        if min_figsize is None:
            min_figsize = (cur_w, cur_h)
        need_w, need_h = _estimate_needed_inches(formatted.shape[0], formatted.shape[1], col_lengths, fs)
        new_w = max(min_figsize[0], need_w)
        new_h = max(min_figsize[1], need_h)
        if max_figsize is not None:
            new_w = min(new_w, max_figsize[0])
            new_h = min(new_h, max_figsize[1])
        # Only enlarge (don’t shrink here)
        if new_w > cur_w or new_h > cur_h:
            fig.set_size_inches(new_w, new_h, forward=True)

    # ---- 5) Shrink-to-fit (gentle)
    def _fits():
        fig.canvas.draw()
        ren = fig.canvas.get_renderer()
        tbbox = tbl.get_window_extent(renderer=ren)
        abbox = ax.get_window_extent(renderer=ren)
        margin = 2  # px
        return (tbbox.width <= abbox.width - margin) and (tbbox.height <= abbox.height - margin)

    steps = 0
    while not _fits() and steps < max_shrink_steps:
        steps += 1
        if fs > min_fontsize:
            fs -= 1
            tbl.set_fontsize(fs)
        else:
            # gentle scaling if we already hit min font
            tbl.scale(0.98, 0.97)

    # ---- 6) Final fallback: monospace dump
    if not _fits():
        try:
            tbl.remove()
        except Exception:
            pass
        preview = df_disp.copy()
        with pd.option_context("display.max_rows", max_rows, "display.max_columns", max_cols, "display.width", 1000):
            text = preview.to_string(index=False)
        ax.text(0.0, 1.0, text, ha="left", va="top", fontsize=9, family="monospace")
        note = []
        if len(df) > max_rows: note.append(f"showing first {max_rows}/{len(df)} rows")
        if df.shape[1] > max_cols: note.append(f"first {max_cols}/{df.shape[1]} cols")
        if note:
            ax.text(0.5, -0.06, " · ".join(note), transform=ax.transAxes,
                    ha="center", va="top", fontsize=8, alpha=0.6)
        if debug:
            ax.text(0.5, -0.10, "FALLBACK: monospace dump (table too large)", transform=ax.transAxes,
                    ha="center", va="top", fontsize=8, color="crimson")
        return

    # ---- 7) Truncation note when it *does* fit
    notes = []
    if len(df) > max_rows:
        notes.append(f"showing first {max_rows} rows of {len(df)}")
    if df.shape[1] > max_cols:
        notes.append(f"showing first {max_cols} columns of {df.shape[1]}")
    if notes:
        ax.text(0.5, -0.08, " · ".join(notes), transform=ax.transAxes,
                ha="center", va="top", fontsize=8, alpha=0.6)


def _try(func, default=None, note: Optional[str] = None):
    """Call func() and return default on failure; optionally print a note."""
    try:
        return func()
    except Exception as e:
        if note:
            print(f"[PDF] {note}: {e}")
            # print on what line the error occurred
            print(f"[PDF] Error occurred on line {e.__traceback__.tb_lineno}")
        return default

def _save_top5_pdf_regression(results: pd.DataFrame, model, pdf_path: str = "top_models_regression_report.pdf", k: int = 5):
    """
    Build a multi-page PDF summarizing the top-k models by R² (default 5).
    Structure:
      • Cover page with ranked table (Formula | R² | #Features | #Samples)
      • For each model:
          Page 1: Coefficients, VIF, CV metrics, meta
          Page 2: Q² scatter plot(s)
          Page 3+: SHAP plots (if available)
    """
    # --- Validate and pick top-k ---
    if "combination" not in results.columns:
        raise KeyError("results must have a 'combination' column")
    r2_col = "r2" if "r2" in results.columns else None
    if r2_col is None:
        raise KeyError("results must have an 'r2' column for ranking")

    res = results.copy()
    res = res.loc[res["combination"].notna() & np.isfinite(res[r2_col])]

    if len(res) == 0:
        raise ValueError("No valid rows to rank.")

    # Prepare parsed features and a tidy printable formula
    parsed_features = res["combination"].apply(_normalize_combination_to_columns)
    pretty_formula = parsed_features.apply(_stringify_formula)

    # Rank
    order = np.argsort(-res[r2_col].values)
    top_idx = order[: min(k, len(order))]

    # --- Cover page data ---
    cover_df = pd.DataFrame({
        "Rank": np.arange(1, len(top_idx) + 1),
        "R² (train)": res.iloc[top_idx][r2_col].round(3).values,
        "#Features": parsed_features.iloc[top_idx].apply(len).values,
        "#Samples": [len(model.target_vector)] * len(top_idx),
        "Formula": pretty_formula.iloc[top_idx].values
    })

    page_counter = 0
    with PdfPages(pdf_path) as pdf:
        # ---------------- COVER PAGE ----------------
        fig = plt.figure(figsize=(11.5, 8.2))
        try:
            _page_header(fig, f"Top {len(top_idx)} Models (ranked by R²)")
        except Exception:
            fig.suptitle(f"Top {len(top_idx)} Models (ranked by R²)", y=0.98, fontsize=16, fontweight="normal")

        gs = GridSpec(1, 1, figure=fig)
        ax = fig.add_subplot(gs[0, 0])
        _safe_table(ax, cover_df, title="Overview")
        page_counter += 1
        _footer(fig, left_text=f"MolFeatures report", page_num=page_counter)
        fig.tight_layout(rect=[0, 0, 1, 0.95])
        pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

        # ------------- PER-MODEL PAGES -------------
        for rank, idx in enumerate(top_idx, start=1):
            feats = parsed_features.iloc[idx]
            formula_str = pretty_formula.iloc[idx]
            r2_val = float(res.iloc[idx][r2_col])

            # --- Extract and fit ---
            # Skip if features missing
            missing = [f for f in feats if f not in model.features_df.columns]
            if missing:
                print(f"[PDF] SKIP model rank {rank}: missing features {missing}")
                continue

            X_df = model.features_df[feats]
            X = X_df.to_numpy()
            y = np.asarray(model.target_vector)

            # Fit & predict with intervals
            _ = model.fit(X, y)
            pred, lwr_band, upr_band = _try(lambda: model.predict(X, return_interval=True),
                                            default=(model.predict(X), None, None),
                                            note="predict with intervals failed")

            # Coefficients & VIF
            coef_df = _try(lambda: model.get_covariance_matrix(feats), default=None, note="get_covariance_matrix failed")
            if coef_df is not None:
                # pick a sensible numeric coefficient column for table highlight
                num_cols = coef_df.select_dtypes(include="number").columns
                if len(num_cols) and "Estimate" in coef_df.columns:
                    coef_df = coef_df[["Estimate"] + [c for c in coef_df.columns if c != "Estimate"]]
                elif len(num_cols):
                    lead = [num_cols[0]] + [c for c in coef_df.columns if c != num_cols[0]]
                    coef_df = coef_df[lead]
            vif_df = _try(lambda: model._compute_vif(X_df), default=None, note="_compute_vif failed")

            # CV metrics
            Q2_3, MAE_3, RMSD_3 = _try(lambda: model.calculate_q2_and_mae(X, y, n_splits=3,plot=False),
                                       default=(np.nan, np.nan, np.nan), note="Q2/MAE/RMSD (3-fold) failed")
            Q2_5, MAE_5, RMSD_5 = _try(lambda: model.calculate_q2_and_mae(X, y, n_splits=5,plot=False),
                                       default=(np.nan, np.nan, np.nan), note="Q2/MAE/RMSD (5-fold) failed")
            Q2_loo, MAE_loo, RMSD_loo = _try(lambda: model.calculate_q2_and_mae(X, y, n_splits=1,plot=False),
                                             default=(np.nan, np.nan, np.nan), note="Q2/MAE/RMSD (LOOCV) failed")

            folds_df = pd.DataFrame({
                'Q2_3_Fold': [Q2_3], 'MAE (3)': [MAE_3], 'RMSD (3)': [RMSD_3],
                'Q2_5_Fold': [Q2_5], 'MAE (5)': [MAE_5], 'RMSD (5)': [RMSD_5],
                'Q2_LOOCV': [Q2_loo], 'MAE (LOO)': [MAE_loo], 'RMSD (LOO)': [RMSD_loo],
            }).round(4)
            
            # -------- PAGE 1: Summary dashboard --------
            fig1 = plt.figure(figsize=(11.5, 8.2))
            title = f"Rank #{rank} | R²={r2_val:.3f} | Features={len(feats)} | Samples={len(y)}"
            try:
                _page_header(fig1, f"{title}\n{formula_str}")
            except Exception:
                fig1.suptitle(f"{title}\n{formula_str}", y=0.98, fontsize=14, fontweight="normal")

            # Adapt the height ratios to give more space to long coefficient tables
            coef_rows = 0 if coef_df is None else getattr(coef_df, "shape", (0, 0))[0]
            top_ratio = float(np.clip(0.9 + 0.02 * min(coef_rows, 35), 0.9, 1.6))
            gs = GridSpec(
                2, 3, figure=fig1,
                height_ratios=[top_ratio, max(0.8, 2.1 - top_ratio)],
                width_ratios=[1.3, 1.0, 1.0]
            )

            ax_coef = fig1.add_subplot(gs[0, :])
            ax_vif  = fig1.add_subplot(gs[1, 0])
            ax_cv   = fig1.add_subplot(gs[1, 1])

            def _render_cell_table(ax, df, title):
                """Choose readable table parameters from the subplot size, then draw."""
                ax.figure.canvas.draw()  # ensure renderer exists for size calc
                bbox = ax.get_window_extent()
                w_px, h_px = float(bbox.width), float(bbox.height)
                ncols = 0 if df is None else max(1, df.shape[1])

                # Heuristics: per-column width -> wrap width; height -> max rows; height -> font size
                per_col_px = w_px / ncols
                if per_col_px < 90:       wrap = 16
                elif per_col_px < 140:    wrap = 24
                else:                     wrap = 36

                # Roughly ~20 px per row + header; keep within sensible bounds
                max_rows = max(8, min(60, int(h_px / 20) - 1))

                start_fs = 12 if h_px > 320 else 10
                min_fs   = 8

                _safe_table(
                    ax, df, title,
                    max_rows=max_rows,
                    max_cols=min(ncols, 12),
                    wrap_width=wrap,
                    start_fontsize=start_fs,
                    min_fontsize=min_fs,
                    cell_scale=1.0,
                    max_shrink_steps=16,
                    auto_resize_figure=False,   # don't resize the whole page layout
                    target_col_total=0.96,      # use most of the subplot width
                    zebra=True
                )

            _render_cell_table(ax_coef, coef_df, "Coefficients")
            _render_cell_table(ax_vif,  vif_df,  "VIF")
            _render_cell_table(ax_cv,   folds_df.round(3), "Cross-validation metrics")
            # Meta box
            ax_meta = fig1.add_subplot(gs[1, 2]); ax_meta.axis("off")
            ax_meta.set_title("Model Summary", pad=8, fontsize=11, fontweight="normal")
            meta_txt = [
                f"Formula: {formula_str}",
                f"Train R²: {r2_val:.3f}",
                f"Q² (3/5/LOO): {Q2_3:.3f} / {Q2_5:.3f} / {Q2_loo:.3f}",
            ]
            ax_meta.text(0.02, 0.98, "\n".join(meta_txt), va="top", ha="left", fontsize=10)

            page_counter += 1
            _footer(fig1, left_text="MolFeatures • Summary", page_num=page_counter)
            fig1.tight_layout(rect=[0, 0, 1, 0.95])
            pdf.savefig(fig1, bbox_inches="tight"); plt.close(fig1)

            # -------- PAGE 2: Pred vs Obs / Q² scatter (your function) --------
            def _run_q2_plot():
                coef_col = None
                if coef_df is not None:
                    # choose a numeric column to pass (your function expects 'Estimate' historically)
                    if "Estimate" in coef_df.columns:
                        coef_col = coef_df["Estimate"]
                    else:
                        num_cols = coef_df.select_dtypes(include="number").columns
                        if len(num_cols): coef_col = coef_df[num_cols[0]]
                fig_q, _ = generate_q2_scatter_plot(
                    y=y, y_pred=pred, labels=getattr(model, "molecule_names", None),
                    folds_df=folds_df, formula=feats, coefficients=coef_col if coef_col is not None else pd.Series(dtype=float),
                    r=r2_val, plot=False
                )

            q2_figs = _try(lambda: _capture_new_figs(_run_q2_plot), default=[], note="Q² scatter generation failed")
            if q2_figs:
                for fig in q2_figs:
                    page_counter += 1
                    _footer(fig, left_text="MolFeatures • Q² scatter", page_num=page_counter)
                    pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)
            else:
                # Provide a fallback page
                fig_fallback = plt.figure(figsize=(11.5, 8.2))
                try:
                    _page_header(fig_fallback, "Q² scatter")
                except Exception:
                    fig_fallback.suptitle("Q² scatter", y=0.98, fontsize=14, fontweight="normal")
                ax = fig_fallback.add_subplot(111); ax.axis("off")
                ax.text(0.5, 0.5, "Q² scatter plot unavailable.", ha="center", va="center", fontsize=11, alpha=0.7)
                page_counter += 1
                _footer(fig_fallback, left_text="MolFeatures • Q² scatter", page_num=page_counter)
                pdf.savefig(fig_fallback, bbox_inches="tight"); plt.close(fig_fallback)

            molecule_names = model.molecule_names
            # -------- PAGE 3+: SHAP analysis (optional) --------
            shap_res = _try(
                lambda: analyze_shap_values(
                    model=model,
                    X=model.features_df[feats],
                    feature_names=feats,
                    sample_names=molecule_names,
                    target_name=getattr(model, "output_name", "target"),
                    n_top_features=min(10, len(feats)),
                    plot=True
                ),
                default=None, note="SHAP analysis failed"
            )
            if shap_res and isinstance(shap_res, dict):
                fig = shap_res.get('figure',[])
                fig.text(0.01, 0.01, f"Rank #{rank} | SHAP", fontsize=8, ha="left", va="bottom", alpha=0.7)
                page_counter += 1
                _footer(fig, left_text="MolFeatures • SHAP", page_num=page_counter)
                pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)
                
      
            # fig = _try(lambda: plot_feature_vs_target(feature_values, target_vector, feature_name=feature, point_labels=molecule_names, figsize=(adjusted_figsize[0]/n_cols*2, adjusted_figsize[1]/n_rows*1.5)), default=None, note="Feature vs Target plot failed")
            if fig:
                page_counter += 1
                _footer(fig, left_text="MolFeatures • Feature vs Target", page_num=page_counter)
                pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

    print(f"[PDF] Saved top-{len(top_idx)} models report to: {pdf_path}")

def _parse_tuple_string(s: str):
    # "('L_11-6', 'buried_volume')" -> ['L_11-6','buried_volume']
    return [x.strip(" '") for x in s.strip("()").split(",")]

def print_models_regression_table(results, app=None ,model=None):
   
    formulas=results['combination'].values
    r_squared=results['r2'].values
    q_squared=results['q2'].values
    mae=results['mae'].values
    model_ids=[i for i in range(len(results))]


    # Create a DataFrame from the inputs
    df = pd.DataFrame({
        'formula': formulas,
        'R.sq': r_squared,
        'Q.sq': q_squared,
        'MAE': mae,
        'Model_id': model_ids
    })

    df = df.sort_values(by='Q.sq', ascending=False)
    df.index = range(1, len(df) + 1)
    print(df.head().to_markdown(index=False, tablefmt="pipe"))
    try:
        pdf_path = model.paths.pdf / f"{model.name}_top_models_report.pdf"
        # _save_top5_pdf_regression(results,model, pdf_path=pdf_path)
    except Exception as e:
        print(f"[PDF] Skipping top-5 export due to error: {e}")

    while True:

        if app:
            print('debug')
        else:
            # print(df.head().to_markdown(index=False, tablefmt="pipe"))
            try:
                selected_model = int(input("Select a model number (or -1 to exit): "))
            except ValueError:
                print("Invalid input. Please enter a number.")
                continue
            
        if selected_model == -1:
            print("Exiting model selection.")
            break


        s = formulas[selected_model]
        features = _normalize_combination_to_columns(s)
   
        X = model.features_df[features]
        vif_df = model._compute_vif(X)
       
        X = model.features_df[features].to_numpy()
        y = model.target_vector.to_numpy()
        model.fit(X, y)
        pred, lwr, upr = model.predict(X, return_interval=True)
        coef_df = model.get_covariance_matrix(features)
        leftout_pred_df, r2_leftout, mae_leftout = _compute_leftout(model, features)
        x_min, y_min = y.min(), y.min()
        x_max, y_max = pred.max(), pred.max()
        padding_x = (x_max - x_min) * 0.05
        padding_y = (y_max - y_min) * 0.05
        lwr = [x_min - padding_x, y_min - padding_y]
        upr = [x_max + padding_x, y_max + padding_y]

        # Debug statements
       

        # Validate bounds
        if not (isinstance(lwr, (list, tuple, np.ndarray)) and len(lwr) == 2 and np.all(np.isfinite(lwr))):
            print("Invalid lower_bound. Using default.")
            lwr = None
        if not (isinstance(upr, (list, tuple, np.ndarray)) and len(upr) == 2 and np.all(np.isfinite(upr))):
            print("Invalid upper_bound. Using default.")
            upr = None
        
        if app:
            app.show_result('\nModel Coefficients\n')
            app.show_result(coef_df.to_markdown(tablefmt="pipe"))
            print_models_vif_table(vif_df, app)
        else:
            print("\nModel Coefficients\n")
            print(coef_df.to_markdown(tablefmt="pipe"))
            print(f"\nSelected Model: {formulas[selected_model]}\n")
            print_models_vif_table(vif_df)
        
        Q2_3, MAE_3 , rmsd_3 = model.calculate_q2_and_mae(X, y, n_splits=3)
        Q2_5, MAE_5, rmsd_5 = model.calculate_q2_and_mae(X, y, n_splits=5)
        ## LOOCV
        Q2_loo, MAE_loo , rmsd_loo = model.calculate_q2_and_mae(X, y, n_splits=1)
        
        if app:
            app.show_result(f'\n\n Model Picked: {selected_model}_{formulas[selected_model]}\n')
            app.show_result(pd.DataFrame({'Q2_3_Fold': [Q2_3], 'MAE': [MAE_3], 'RMSD':[rmsd_3]}).to_markdown(tablefmt="pipe", index=False))
            app.show_result(pd.DataFrame({'Q2_5_Fold':[Q2_5], 'MAE': [MAE_5], 'RMSD':[rmsd_5]}).to_markdown(tablefmt="pipe", index=False))
            app.show_result(pd.DataFrame({'Q2_LOOCV':[Q2_loo], 'MAE': [MAE_loo], 'RMSD':[rmsd_loo]}).to_markdown(tablefmt="pipe", index=False))
        else:
            print("\n3-fold CV\n")
            print(pd.DataFrame({'Q2_3_Fold': [Q2_3], 'MAE': [MAE_3]}).to_markdown(tablefmt="pipe", index=False))
            print("\n5-fold CV\n")
            print(pd.DataFrame({'Q2_5_Fold':[Q2_5], 'MAE': [MAE_5]}).to_markdown(tablefmt="pipe", index=False))
            print("\nLOOCV\n")
            print(pd.DataFrame({'Q2_LOOCV':[Q2_loo], 'MAE': [MAE_loo]}).to_markdown(tablefmt="pipe", index=False))
        
        txt_path = model.paths.logs / 'regression_results.txt'
        with open(txt_path, 'a') as f:
            f.write(
                f"Models list {df.to_markdown(index=False, tablefmt='pipe')} \n\n"
                f"Model Coefficients\n\n{coef_df.to_markdown(tablefmt='pipe')}\n\n"
                f"3-fold CV\n\n{pd.DataFrame({'Q2': [Q2_3], 'MAE': [MAE_3]}).to_markdown(tablefmt='pipe', index=False)}\n\n"
                f"5-fold CV\n\n{pd.DataFrame({'Q2':[Q2_5], 'MAE': [MAE_5]}).to_markdown(tablefmt='pipe', index=False)}\n\n"
            )
            print('Results saved to regression_results.txt in {}'.format(os.getcwd()))

            # ✅ Build folds_df with unique column names
            folds_df = pd.DataFrame({
                'Q2_3_Fold': [Q2_3], 'MAE_3': [MAE_3], 'RMSD_3': [rmsd_3],
                'Q2_5_Fold': [Q2_5], 'MAE_5': [MAE_5], 'RMSD_5': [rmsd_5],
                'Q2_LOOCV': [Q2_loo], 'MAE_LOO': [MAE_loo], 'RMSD_LOO': [rmsd_loo]
            })

            r = r_squared[selected_model]

            # ✅ Generate scatter plot safely (don’t auto-show inside the function)
            lb = None if lwr is None else float(np.nanmin(lwr))
            ub = None if upr is None else float(np.nanmax(upr))
           
            _, ax_q2 = generate_q2_scatter_plot(
                y=y,
                y_pred=pred,
                labels=model.molecule_names,
                folds_df=folds_df,
                formula=features,
                coefficients=coef_df['Estimate'],
                r=r,
                lower_bound=lb,                 # <- keyword, scalar
                upper_bound=ub,                 # <- keyword, scalar
                plot=False,
                dir=model.paths.figs,
                figsize=(10, 5),                # wider
                equal_aspect=False,             # let it expand horizontally
                fontsize=10,                    # bigger axes text
                label_fontsize=8,               # smaller point labels
                leftout_pred_df=leftout_pred_df
            )

            ip = get_ipython()
            is_colab = "google.colab" in str(ip).lower() if ip else False

            # # --- show figure non-blocking ---
            # if ip is not None:  # Notebook / Colab
            #     from IPython.display import display
            #     display(fig_q2)       # render inline immediately
            #     plt.pause(0.5)        # <-- short pause lets Colab finish rendering
            # else:                     # Terminal / script
            #     plt.show(block=False)
            #     plt.pause(0.5)

            # --- handle prompt logic ---
            if not app:
                if is_colab:
                    # Colab can't handle `input()` interactively — fallback
                    print("\n[Colab detected] Skipping interactive input prompt.")
                    cont = "n"  # or auto-continue if you prefer: cont = "y"
                else:
                    cont = input("Do you want to select another model? (y/n): ").strip().lower()

                if cont != "y":
                    print("Exiting model selection.")
                    break
            else:
                cont = messagebox.askyesno("Continue", "Do you want to select another model?")
                if not cont:
                    break



import re
import unicodedata
import hashlib


def run_model_sanity_checks(
    model, X, y, features, Q2_loo=None, folds_df=None, n_runs=100, random_state=42,
    pdf=None, show=False, png_dir=None
):
    """
    Run Y-randomization, one-hot test, global/per-feature shuffling.
    Saves each figure to PDF (if provided) and PNG (if png_dir provided).
    """
    # ----- normalize output dir ONCE and create it -----
    png_dir = Path(png_dir) if png_dir is not None else None
    
    if png_dir is not None:
        png_dir.mkdir(parents=True, exist_ok=True)

    print(f"[SANITY] Running sanity checks on model '{model.name}' with {X.shape[1]} features and {X.shape[0]} samples..., path={png_dir}")

    try:
        rng = np.random.default_rng(random_state)
        if folds_df is None:
            folds_df = pd.DataFrame()

        # --- Baseline (real) metrics on the given X,y (5-fold CV) ---
        _, MAE_real, RMSD_real = model.calculate_q2_and_mae(X, y, n_splits=5)

        # --- Y-randomization (distributions of MAE/RMSD) ---
        random_mae, random_rmsd = [], []
        for _ in range(n_runs):
            y_perm = rng.permutation(y)
            _, mae_p, rmsd_p = model.calculate_q2_and_mae(X, y_perm, n_splits=5)
            random_mae.append(mae_p)
            random_rmsd.append(rmsd_p)

        folds_df["MAE_random_mean"]  = [float(np.mean(random_mae))]
        folds_df["MAE_random_std"]   = [float(np.std(random_mae))]
        folds_df["RMSD_random_mean"] = [float(np.mean(random_rmsd))]
        folds_df["RMSD_random_std"]  = [float(np.std(random_rmsd))]

        # --- One-hot encoding test ---
        onehot_df = pd.get_dummies(model.molecule_names, prefix="mol")
        X_onehot  = onehot_df.to_numpy()
        y_true    = model.target_vector.to_numpy()
        _, MAE_oh, RMSD_oh = model.calculate_q2_and_mae(X_onehot, y_true, n_splits=5)
        folds_df["MAE_onehot"]  = [MAE_oh]
        folds_df["RMSD_onehot"] = [RMSD_oh]

        # --- Global X-shuffling ---
        X_shuffled = np.apply_along_axis(rng.permutation, 0, X.copy())
        _, MAE_xshuf, RMSD_xshuf = model.calculate_q2_and_mae(X_shuffled, y, n_splits=5)
        folds_df["MAE_global_shuffle"]  = [MAE_xshuf]
        folds_df["RMSD_global_shuffle"] = [RMSD_xshuf]

        # --- Per-feature shuffling ---
        mae_per_feat, rmsd_per_feat = [], []
        for j in range(X.shape[1]):
            X_tmp = X.copy()
            rng.shuffle(X_tmp[:, j])
            _, mae_j, rmsd_j = model.calculate_q2_and_mae(X_tmp, y, n_splits=5)
            mae_per_feat.append(mae_j)
            rmsd_per_feat.append(rmsd_j)

        # ---------- helper (no shadowing; dir already exists) ----------
        def _save_current_fig(stem):
            fig = plt.gcf()
            plt.show()
            if pdf is not None:
                pdf.savefig(fig)
            if png_dir is not None:
                fig_path = png_dir / f"{stem}.png"
                fig.savefig(fig_path, dpi=300, bbox_inches="tight")
                print(f"[FIG] Saved {stem} plot to: {fig_path}")
            plt.close(fig)

        # =======================
        # Figures
        # =======================

        # 1) Y-randomization MAE histogram
        plt.figure(figsize=(8, 5))
        plt.hist(random_mae, bins=12, color="lightblue", edgecolor="black", alpha=0.7, label="Y-random MAE")
        plt.axvline(MAE_real,   color="red",   lw=2, label=f"Real MAE = {MAE_real:.3f}")
        plt.axvline(np.mean(random_mae), color="blue",  lw=2, ls="--", label=f"Mean Y-rand MAE = {np.mean(random_mae):.3f}")
        plt.axvline(MAE_oh,     color="green", lw=2, ls="--", label=f"One-hot MAE = {MAE_oh:.3f}")
        plt.axvline(MAE_xshuf,  color="purple",lw=2, ls="--", label=f"Global X-shuffle MAE = {MAE_xshuf:.3f}")
        plt.xlabel("MAE (5-fold CV)"); plt.ylabel("Count"); plt.title("Validation (lower is better): MAE"); plt.legend(); plt.tight_layout()
        if show: plt.show()
        _save_current_fig("sanity_checks_validation_MAE")

        # 2) Y-randomization RMSD histogram
        plt.figure(figsize=(8, 5))
        plt.hist(random_rmsd, bins=12, color="lightblue", edgecolor="black", alpha=0.7, label="Y-random RMSD")
        plt.axvline(RMSD_real,   color="red",   lw=2, label=f"Real RMSD = {RMSD_real:.3f}")
        plt.axvline(np.mean(random_rmsd), color="blue",  lw=2, ls="--", label=f"Mean Y-rand RMSD = {np.mean(random_rmsd):.3f}")
        plt.axvline(RMSD_oh,     color="green", lw=2, ls="--", label=f"One-hot RMSD = {RMSD_oh:.3f}")
        plt.axvline(RMSD_xshuf,  color="purple",lw=2, ls="--", label=f"Global X-shuffle RMSD = {RMSD_xshuf:.3f}")
        plt.xlabel("RMSD (5-fold CV)"); plt.ylabel("Count"); plt.title("Validation (lower is better): RMSD"); plt.legend(); plt.tight_layout()
        if show: plt.show()
        _save_current_fig("sanity_checks_validation_RMSD")

        # 3) Per-feature shuffle bar (MAE)
        plt.figure(figsize=(1.5*X.shape[1], 5))
        plt.bar(range(len(features)), mae_per_feat, color="orange", edgecolor="black", alpha=0.8)
        plt.axhline(MAE_real, color="red", ls="--", label=f"Real MAE = {MAE_real:.3f}")
        plt.xticks(range(len(features)), features, rotation=30, ha="right")
        plt.ylabel("MAE (5-fold CV)"); plt.title("Per-feature shuffle test: MAE (lower is better)"); plt.legend(); plt.tight_layout()
        if show: plt.show()
        _save_current_fig("sanity_checks_per_feature_shuffle_MAE")

        # 4) Per-feature shuffle bar (RMSD)
        plt.figure(figsize=(1.5*X.shape[1], 5))
        plt.bar(range(len(features)), rmsd_per_feat, color="orange", edgecolor="black", alpha=0.8)
        plt.axhline(RMSD_real, color="red", ls="--", label=f"Real RMSD = {RMSD_real:.3f}")
        plt.xticks(range(len(features)), features, rotation=30, ha="right")
        plt.ylabel("RMSD (5-fold CV)"); plt.title("Per-feature shuffle test: RMSD (lower is better)"); plt.legend(); plt.tight_layout()
        if show: plt.show()
        _save_current_fig("sanity_checks_per_feature_shuffle_RMSD")

        # Package per-feature metrics
        per_feature_metrics = pd.DataFrame({
            "feature": list(features),
            "MAE_permute": mae_per_feat,
            "RMSD_permute": rmsd_per_feat
        })

        # Store real metrics too
        folds_df["MAE_real"]  = [MAE_real]
        folds_df["RMSD_real"] = [RMSD_real]

        return {"folds_df": folds_df, "per_feature_metrics": per_feature_metrics}

    except Exception as e:
        print("Error during sanity checks:", e)
        return {
            "folds_df": folds_df if folds_df is not None else pd.DataFrame(),
            "per_feature_metrics": pd.DataFrame()
        }
    

def _slugify(s: str, maxlen: int = 120) -> str:
    s = str(s)
    s = unicodedata.normalize("NFKD", s).encode("ascii", "ignore").decode("ascii")
    s = s.replace(" ", "_")
    s = re.sub(r"[^\w\-\.]+", "_", s)  # keep word chars, dash, dot, underscore
    s = re.sub(r"_+", "_", s).strip("_.")
    if len(s) <= maxlen:
        return s
    # truncate but keep extension if present
    root, ext = os.path.splitext(s)
    keep = maxlen - len(ext) - 9  # leave room for _trunc and hash
    h = hashlib.sha1(s.encode()).hexdigest()[:8]
    return f"{root[:max(1, keep)]}_trunc_{h}{ext}"

def _shorten_if_needed(path: str, stem: str, ext: str) -> str:
    """
    If path is too long for Windows, shorten filename using a hash.
    """
    try_path = os.path.join(path, f"{stem}{ext}")
    # Heuristic: keep paths under ~240 chars to be safe on Windows
    if len(try_path) < 240:
        return try_path
    h = hashlib.sha1(stem.encode()).hexdigest()[:10]
    short_stem = (stem[:80] + f"_{h}") if len(stem) > 100 else f"{stem}_{h}"
    return os.path.join(path, f"{short_stem}{ext}")

def draw_table_on_ax(ax, df, title, max_rows=12, fontsize=8, round_digits=None):
    """Draw a compact table (head of df) on a given axes."""
    ax.axis("off")
    ax.set_title(title, fontsize=12, pad=6)
    if df is None or df.empty:
        ax.text(0.5, 0.5, "(no data)", ha="center", va="center", fontsize=fontsize+1)
        return

    df_show = df.copy()
    if round_digits is not None:
        with np.errstate(all='ignore'):
            df_show = df_show.round(round_digits)
    df_show = df_show.head(max_rows)

    tbl = ax.table(
        cellText=df_show.values,
        colLabels=[str(c) for c in df_show.columns],
        loc="center",
        cellLoc="left",
        rowLoc="center",
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(fontsize)
    tbl.scale(1.0, 1.15)

    # Style header row
    cells = tbl.get_celld()
    max_row = max(r for (r, c) in cells.keys())
    max_col = max(c for (r, c) in cells.keys())
    for c in range(0, max_col+1):
        cell = cells.get((0, c))
        if cell:
            cell.set_facecolor("#E8EEF6")
            cell.set_edgecolor("#CCCCCC")
            cell.set_text_props(weight="normal")
    # Zebra stripes
    for r in range(1, max_row+1):
        color = "#F8F9FB" if (r % 2 == 0) else "white"
        for c in range(0, max_col+1):
            cell = cells.get((r, c))
            if cell:
                cell.set_facecolor(color)


            

def save_fig_both(fig, pdf: PdfPages, outdir: str, stem: str, dpi: int = 300, tight: bool = True):
    """
    Save figure to the open PdfPages and as a PNG in outdir, robust to Windows path limits.
    """
    os.makedirs(outdir, exist_ok=True)

    # sometimes tables + tight_layout cause warnings; allow opt-out per-call
    try:
        if tight:
            # Use constrained_layout if available to avoid tight_layout warning on GridSpec tables
            fig.set_constrained_layout(True)
    except Exception:
        pass

    # Save to PDF first (does not hit Windows MAX_PATH for filenames)
    pdf.savefig(fig)

    # Build safe PNG path
    safe_stem = _slugify(stem, maxlen=120)
    png_path = _shorten_if_needed(outdir, safe_stem, ".png")

    # Save PNG with fallback strategies
    try:
        fig.canvas.draw()  # ensure render
        fig.savefig(png_path, dpi=dpi, bbox_inches="tight" if tight else None)
    except (FileNotFoundError, OSError) as e:
        # Fallback: shorten further + disable tight bbox
        try:
            very_short = _shorten_if_needed(outdir, hashlib.md5(safe_stem.encode()).hexdigest()[:12], ".png")
            fig.savefig(very_short, dpi=dpi)
        except Exception as e2:
            print(f"[save_fig_both] Failed to save PNG even after shortening: {e2}")
    finally:
        plt.close(fig)

def save_current_fig_both(pdf: PdfPages, outdir: str, stem: str, dpi: int = 300):
    """
    Save the current matplotlib figure both into PDF and as PNG.
    """
    fig = plt.gcf()
    save_fig_both(fig, pdf, outdir, stem, dpi=dpi)


def _extract_Xy(model, features):
    X = model.features_df[features].to_numpy()
    y = model.target_vector.to_numpy()
    return X, y

def _fit_and_predict(model, features):
    result = fit_and_evaluate_single_combination_regression(model, features)
    pred = result["predictions"]
    return pred, result

def _compute_cv_metrics(model, X, y):
    Q2_3, MAE_3, rmsd_3   = model.calculate_q2_and_mae(X, y, n_splits=3, plot=False)
    Q2_5, MAE_5, rmsd_5   = model.calculate_q2_and_mae(X, y, n_splits=5, plot=False)
    Q2_loo, MAE_loo, rmsd_loo = model.calculate_q2_and_mae(X, y, n_splits=1, plot=False)
    folds_df = pd.DataFrame({
        'Q2_3_Fold': [Q2_3], 'MAE_3': [MAE_3], 'RMSD_3': [rmsd_3],
        'Q2_5_Fold': [Q2_5], 'MAE_5': [MAE_5], 'RMSD_5': [rmsd_5],
        'Q2_LOOCV': [Q2_loo], 'MAE_LOOCV': [MAE_loo], 'RMSD_LOOCV': [rmsd_loo]
    })
    return folds_df, Q2_loo

def _compute_leftout(model, features):
    import numpy as np
    import pandas as pd
    from sklearn.metrics import r2_score, mean_absolute_error

    # ---- Guard: no left-out set ----
    if getattr(model, "leftout_samples", None) is None or len(model.leftout_samples) == 0:
        return None, None, None

    # ---- Training data for this feature set (as DataFrame/Series) ----
    X_train = model.features_df[features]
    y_train = model.target_vector

    # ---- Align left-out columns to training schema ----
    # (reindex adds missing columns as NaN; model.predict_for_leftout should impute/handle or you can fillna here)
    X_left = model.leftout_samples.reindex(columns=features)
    y_left = None
    if getattr(model, "leftout_target_vector", None) is not None:
        y_left = np.asarray(model.leftout_target_vector).ravel()

    # ---- Predict (pass X_train/y_train so the method can fit if needed) ----
    pred_out = model.predict_for_leftout(
        X_left,
        y=y_left,
        X_train=X_train,
        y_train=y_train,
        calc_interval=False
    )

    # Unpack predictions
    if isinstance(pred_out, tuple):
        y_pred_left = np.asarray(pred_out[0]).ravel()
    else:
        y_pred_left = np.asarray(pred_out).ravel()

    # ---- Names handling (robust to length mismatches) ----
    if hasattr(model, "molecule_names_predict") and model.molecule_names_predict is not None:
        names = list(model.molecule_names_predict)
    else:
        names = X_left.index.astype(str).tolist()

    if len(names) != len(y_pred_left):
        names = [str(n) for n in range(len(y_pred_left))]  # fallback safe indexing

    # ---- Build result frame ----
    if y_left is not None and len(y_left) == len(y_pred_left):
        df = pd.DataFrame({
            "Molecule": names,
            "Actual": y_left,
            "Predicted": y_pred_left
        })
        df["Error in %"] = np.where(
            df["Actual"] != 0,
            np.abs(df["Actual"] - df["Predicted"]) / np.abs(df["Actual"]) * 100.0,
            np.nan
        )
        r2_leftout = r2_score(y_left, y_pred_left)
        mae_leftout = mean_absolute_error(y_left, y_pred_left)
    else:
        # No targets → just predictions
        df = pd.DataFrame({
            "Molecule": names,
            "Predicted": y_pred_left
        })
        df["Error in %"] = np.nan
        r2_leftout, mae_leftout = None, None

    return df, r2_leftout, mae_leftout


def _compute_margins():
    x_ideal = np.linspace(-3.0, 0.6, 100)
    margin_lines = {
        "identity": x_ideal,
        "+0.25":    x_ideal + 0.25,
        "-0.25":    x_ideal - 0.25,
        "+0.50":    x_ideal + 0.50,
        "-0.50":    x_ideal - 0.50
    }
    return x_ideal, margin_lines

def _compute_vif(model, features):
    try:
        return model._compute_vif(model.features_df[features])
    except Exception:
        return pd.DataFrame()

def _in_sample_stats(y, pred):
    from sklearn.metrics import mean_absolute_error
    r2_in = float(np.corrcoef(y, pred)[0, 1]**2)
    mae = float(mean_absolute_error(y, pred))
    return r2_in, mae

def _coeffs_and_intervals(model, X, features):
    try:
        y_i, upr, lwr = model.predict(X, return_interval=True)
    except Exception:
        y_i, upr, lwr = None, None, None
    try:
        coef_df = model.get_covariance_matrix(features)
    except Exception:
        coef_df = pd.DataFrame()
    return (y_i, upr, lwr), coef_df

# ----------------------------- paths & pdf -----------------------------------

def _prepare_paths(model, features, pdf_name=None):
    figs_dir = getattr(getattr(model, "paths", None), "figs", ".")
    os.makedirs(figs_dir, exist_ok=True)

    feat_slug = _slugify("_".join(features)) or "model"
    base_name = f"report_{feat_slug}"
    pdf_path = os.path.join(figs_dir, f"{base_name}.pdf")
    
    png_dir = os.path.join(figs_dir, base_name)
    try:
        os.makedirs(png_dir, exist_ok=True)
        model.paths.png_dir = png_dir
    except:
        print('Name is too long creating dir with default name: Model_Res')
        png_dir = os.path.join(figs_dir, 'Model_Res')
        os.makedirs(png_dir, exist_ok=True)
        model.paths.png_dir = png_dir
    return figs_dir, base_name, pdf_path, png_dir

# ----------------------------- page builders ---------------------------------

def _add_summary_page(pdf, png_dir, base_name, folds_df, vif_df, coef_df, leftout_pred_df):
    fig = plt.figure(figsize=(11.0, 8.5))
    gs = fig.add_gridspec(nrows=3, ncols=3, height_ratios=[1.0, 1.0, 0.8], hspace=0.6, wspace=0.4)

    # Row 0: CV metrics
    ax_cv3 = fig.add_subplot(gs[0, 0])
    draw_table_on_ax(ax_cv3, folds_df[["Q2_3_Fold", "MAE_3", "RMSD_3"]], "3-fold CV Metrics", max_rows=12, round_digits=3)

    ax_cv5 = fig.add_subplot(gs[0, 1])
    draw_table_on_ax(ax_cv5, folds_df[["Q2_5_Fold", "MAE_5", "RMSD_5"]], "5-fold CV Metrics", max_rows=12, round_digits=3)

    ax_loo = fig.add_subplot(gs[0, 2])
    draw_table_on_ax(ax_loo, folds_df[["Q2_LOOCV", "MAE_LOOCV", "RMSD_LOOCV"]], "LOOCV Metrics", max_rows=12, round_digits=3)

    # Row 1: VIF + Coef
    ax_vif = fig.add_subplot(gs[1, 0:2])
    vif_for_pdf = vif_df.copy()
    if not vif_for_pdf.empty and "VIF" in vif_for_pdf.columns:
        vif_for_pdf = vif_for_pdf.sort_values("VIF", ascending=False)
    draw_table_on_ax(ax_vif, vif_for_pdf, "Variance Inflation Factors", max_rows=12, round_digits=2)

    ax_coef = fig.add_subplot(gs[1, 2])
    coef_for_pdf = coef_df.copy()
    preferred = [c for c in ["Estimate", "Std. Error", "t value", "Pr(>|t|)", "p-value"] if c in coef_for_pdf.columns]
    if preferred:
        rest = [c for c in coef_for_pdf.columns if c not in preferred]
        coef_for_pdf = coef_for_pdf[preferred + rest]
    draw_table_on_ax(ax_coef, coef_for_pdf, "Coefficient Estimates", max_rows=12, round_digits=3)

    # Row 2: Leftout table
    ax_leftout = fig.add_subplot(gs[2, :])
    if leftout_pred_df is not None:
        draw_table_on_ax(ax_leftout, leftout_pred_df, "Left-out Predictions", max_rows=12, round_digits=3)
    else:
        ax_leftout.axis("off")
        ax_leftout.set_title("Left-out Predictions (none)", fontsize=12, pad=6)

    # plt.tight_layout()
    save_fig_both(fig, pdf, png_dir, f"{base_name}__tables_summary")

def _add_q2_scatter_page(pdf, png_dir, base_name, y, pred, names, folds_df, features, coef_df, r2_in, lig_types=None,leftout_pred_df=None):
    print("Calling generate_q2_scatter_plot...")
    print('predicitions,',leftout_pred_df)
    est = coef_df['Estimate'] if 'Estimate' in coef_df.columns else None
    fig_q2, ax_q2 = generate_q2_scatter_plot(y, pred, names, folds_df, features, est, r2_in, plot=True, dir=None, ligand_types=lig_types,leftout_pred_df=leftout_pred_df)
    plt.show()
    save_fig_both(fig_q2, pdf, png_dir, f"{base_name}__03_pred_vs_meas")

import mplcursors
def _add_violin_page(pdf, png_dir, base_name, model, features):
    print("Generating violin plots for selected features...")
    fig, ax = plt.subplots(figsize=(1.5 * len(features), 6))
    sns.violinplot(data=model.features_df[features], inner="point", cut=0, linewidth=1, ax=ax)
    ax.set_title("Distribution of Selected Features", fontsize=14)
    ax.set_ylabel("Feature Value", fontsize=12)
    ax.set_xlabel("Feature", fontsize=12)
    plt.xticks(rotation=30, ha="right")
    # plt.tight_layout()
    plt.show()
    save_fig_both(fig, pdf, png_dir, f"{base_name}__01_violin")

def _add_shap_page(pdf, png_dir, base_name, model, features):
    try:
     
        shap_res = _try(
                lambda: analyze_shap_values(
                    model=model,
                    X=model.features_df[features],
                    feature_names=features,
                    sample_names=model.molecule_names,
                    target_name=getattr(model, "output_name", "target"),
                    n_top_features=min(10, len(features)),
                    plot=True,
                    dir=png_dir
                ),
                default=None, note="SHAP analysis failed"
            )
        fig = shap_res.get("figure", [])
        plt.show()
        if fig is None:
            save_current_fig_both(pdf, png_dir, f"{base_name}__04_shap")
        else:
            save_fig_both(fig, pdf, png_dir, f"{base_name}__04_shap")
    except Exception as e:
        print("Error generating SHAP plot:", e)

def _add_threshold_pages(pdf, png_dir, base_name, model, features):
    _, figs = threshold_analysis_plot(
        model.target_vector,
        model.features_df[features],
        cutoff=("percentile", 30, "low_is_positive"),
        data_labels=model.molecule_names
    )
    for i, f in enumerate(figs):
        save_fig_both(f, pdf, png_dir, f"{base_name}__02_threshold_analysis_{i+1}")
        plt.close(f)

def _run_sanity(model, X, y, features, Q2_loo, pdf, folds_df, png_dir):
    sanity_out = run_model_sanity_checks(
        model, X, y, features, Q2_loo,
        folds_df=folds_df, n_runs=30, random_state=42,
        pdf=pdf, show=False, png_dir=png_dir
    )
    return sanity_out.get("folds_df", folds_df)

# ----------------------------- public wrapper --------------------------------

def run_single_combo_report(model, features, app=None, pdf_name=None, lig_types=None):
    """
    Smaller, testable steps; relies on pre-existing helpers:
    - save_fig_both, _slugify, generate_q2_scatter_plot, draw_table_on_ax,
      run_model_sanity_checks, save_current_fig_both, threshold_analysis_plot
    """
    # 1) data & fit
    X, y = _extract_Xy(model, features)
    pred, _ = _fit_and_predict(model, features)

    # 2) metrics
    folds_df, Q2_loo = _compute_cv_metrics(model, X, y)
    leftout_pred_df, r2_leftout, mae_leftout = _compute_leftout(model, features)

    # 3) aux calcs
    _ = _compute_margins()  # if you want to pass on later
    vif_df = _compute_vif(model, features)
    r2_in, mae_in = _in_sample_stats(y, pred)
    (_, _, _), coef_df = _coeffs_and_intervals(model, X, features)

    # 4) paths
    figs_dir, base_name, pdf_path, png_dir = _prepare_paths(model, features, pdf_name)
    print(f"Saving report to: {pdf_path} and PNGs to: {png_dir}")

    # 5) build pdf
    with PdfPages(pdf_path) as pdf:
        # Q2 scatter first (cover)
        try:
            
            _add_q2_scatter_page(pdf, png_dir, base_name, y, pred, model.molecule_names, folds_df, features, coef_df, r2_in, lig_types=lig_types, leftout_pred_df=leftout_pred_df)
        except Exception as e:
            print("Error in Q2 scatter page:", e)

        # Summary tables page
        try:
            _add_summary_page(pdf, png_dir, base_name, folds_df, vif_df, coef_df, leftout_pred_df)
        except Exception as e:
            print("Error in summary tables page:", e)

        # Sanity checks (updates folds_df possibly)
        try:
            folds_df = _run_sanity(model, X, y, features, Q2_loo, pdf, folds_df, png_dir)
        except Exception as e:
            print("Error during sanity checks:", e)

        # Violin plots
        try:
            _add_violin_page(pdf, png_dir, base_name, model, features)
        except Exception as e:
            print("Error generating violin plots:", e)

        # # SHAP
        # _add_shap_page(pdf, png_dir, base_name, model, features)

        # Threshold pages
        try:
            _add_threshold_pages(pdf, png_dir, base_name, model, features)
        except Exception as e:
            print("Error generating threshold analysis plot:", e)

        try:
            results, _ = check_linear_regression_assumptions(X, y, dir="diag_plots", plot=False)
        except Exception as e:
            print("Error generating regression diagnostics plots:", e)

    results_dict = {
        "pdf_path": pdf_path,
        "png_dir": png_dir,
        "folds_df": folds_df,
        "vif_df": vif_df,
        "coef_df": coef_df,
        "in_sample": {"R2": r2_in, "MAE": mae_in},
        "leftout": {"df": leftout_pred_df, "R2": r2_leftout, "MAE": mae_leftout},
    }
    return results_dict




def plot_feature_vs_target(feature_values, y_values, feature_name, y_name="Target", point_labels=None, figsize=(10, 6), dir=None):
    """
    Plot a single feature against the target variable with optional point labels,
    with x-axis = y_value (target), y-axis = feature_value (feature).
    """
    feature_values = np.asarray(feature_values)
    y_values = np.asarray(y_values)
    
    fig, ax = plt.subplots(figsize=figsize)

    # Scatter: x is y_values (target), y is feature_values
    scatter = ax.scatter(y_values, feature_values, s=70, alpha=0.7, edgecolor='w', linewidth=1)
    
    # Regression line (same axes)
    sns.regplot(x=y_values, y=feature_values, scatter=False, ci=95, line_kws={'color':'red'}, ax=ax)
    
    # Correlation coefficient
    corr_coef = np.corrcoef(y_values, feature_values)[0, 1]
    r_squared = corr_coef**2
    
    # Axis labels (x: target, y: feature)
    ax.set_xlabel(y_name, fontsize=14)
    ax.set_ylabel(feature_name, fontsize=14)
    ax.set_title(f"{feature_name} vs {y_name}\nPearson r = {corr_coef:.3f}, R² = {r_squared:.3f}", fontsize=16)
    ax.grid(True, linestyle='--', alpha=0.7)
    
    # Point labels
    if point_labels is not None:
        texts = []
        for i, label in enumerate(point_labels):
            texts.append(ax.annotate(label, (y_values[i], feature_values[i]),
                                     fontsize=10, ha='right', va='bottom'))
        try:
            adjust_text(texts, arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))
        except ImportError:
            print("adjustText package not found. Labels may overlap.")
    
    plt.tight_layout()
    if dir is not None:
        plt.savefig(os.path.join(dir, f'feature_vs_target_{feature_name}.png'), dpi=300, bbox_inches='tight')
    return fig

def plot_all_features_vs_target(features_df, target_vector, molecule_names=None, figsize=(12, 10), dir=None):
    """
    Plot each feature against the target variable in a grid of subplots.
    
    Parameters:
    -----------
    features_df : pandas.DataFrame
        DataFrame containing all features to plot.
    target_vector : array-like
        The target variable values.
    molecule_names : array-like, optional
        Names to use as point labels.
    figsize : tuple, optional
        Base figure size that will be adjusted based on number of features.
    
    Returns:
    --------
    figs : list
        List of created figure objects.
    """
    features = features_df.columns
    n_features = len(features)
    figs = []
    
    # Calculate grid dimensions
    n_cols = min(3, n_features)
    n_rows = (n_features + n_cols - 1) // n_cols
    
    # Adjust figsize based on the number of subplots
    adjusted_figsize = (figsize[0] * n_cols / 3, figsize[1] * n_rows / 3)

    for i, feature in enumerate(features):

        feature_values = features_df[feature].values
        fig = plot_feature_vs_target(feature_values, target_vector, 
                                    feature_name=feature, 
                                    point_labels=molecule_names,
                                    figsize=(adjusted_figsize[0]/n_cols*2, adjusted_figsize[1]/n_rows*1.5))
        figs.append(fig)
        
        # Save each plot
        plt.savefig(f'feature_plot_{feature}.png', dpi=300, bbox_inches='tight')
    plt.show()
    plt.close('all')  # Close all figures to save memory
    return figs


from typing import Optional, Dict, Any, List, Union
import warnings

import shap
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import warnings
from typing import Any, Dict, List, Optional, Union

def analyze_shap_values(
    model,
    X: Union[pd.DataFrame, np.ndarray, List[List[float]]],
    y: Optional[Union[pd.Series, np.ndarray]] = None,
    feature_names: Optional[List[str]] = None,
    sample_names: Optional[List[str]] = None,
    target_name: str = "output",
    n_top_features: int = 10,
    max_background: int = 100,
    plot: bool = True,
    dir: Optional[str] = None,
    outlier_z: float = 2.5,
    label_outliers: bool = True,
    max_labels_per_feature: int = 10,
) -> Dict[str, Any]:
    """
    Compute SHAP and build static Matplotlib beeswarm-like plots.
    Outlier points (large |SHAP| values) are annotated with sample names.
    """
    model = model.model
    # ---------------------- X handling & names ----------------------
    X_is_df = hasattr(X, "columns")
    if X_is_df:
        X_df = X.copy()
        feature_names = list(X_df.columns)
        X_arr = X_df.values
    else:
        X_arr = np.asarray(X)
        if feature_names is None:
            feature_names = [f"Feature {i}" for i in range(X_arr.shape[1])]
        X_df = pd.DataFrame(X_arr, columns=feature_names)

    n, d = X_arr.shape

    if sample_names is None:
        sample_names = [f"s{i}" for i in range(n)]
    else:
        # force correct length and string type
        sample_names = list(map(str, sample_names))
        if len(sample_names) != n:
            warnings.warn(f"sample_names length {len(sample_names)} != n={n}; regenerating default names.")
            sample_names = [f"s{i}" for i in range(n)]

    # ---------------------- Build explainer ------------------------
    if shap is None:
        warnings.warn("SHAP not available; returning without plots.")
        return {"explainer": None, "shap_values": None, "importance_df": None, "top_features": None}

    try:
        # Try model-aware Explainer first (tree, linear, deep, etc.)
        explainer = shap.Explainer(model, X_df)
        shap_vals = explainer(X_df)
        shap_array = np.array(shap_vals.values)
    except Exception:
        warnings.warn("Falling back to KernelExplainer (may be slow).")
        bg_idx = np.random.choice(n, size=min(max_background, n), replace=False)
        bg = X_arr[bg_idx]
        explainer = shap.KernelExplainer(lambda z: np.asarray(model.predict(z)).ravel(), bg)
        shap_array = np.array(explainer.shap_values(X_arr))

    # ---------------------- Normalize to 2D ------------------------
    # Handle multi-output/class shapes robustly
    if isinstance(shap_array, list):
        # e.g., list per class -> stack [C, N, D] and reduce classes by mean |.| 
        shap_array = np.stack(shap_array, axis=0)
        shap_2d = np.mean(np.abs(shap_array), axis=0)  # [N, D]
    elif shap_array.ndim == 3:
        # common shapes: [N, D, C] or [C, N, D]
        if shap_array.shape[0] == n:
            shap_2d = np.mean(np.abs(shap_array), axis=-1)  # [N, D]
        else:
            shap_2d = np.mean(np.abs(shap_array), axis=0)   # [N, D]
    else:
        # [N, D]
        shap_2d = shap_array

    # ---------------------- Importance ranking ---------------------
    mean_abs = np.mean(np.abs(shap_2d), axis=0)
    importance_df = (
        pd.DataFrame({"Feature": feature_names, "MeanAbsSHAP": mean_abs})
        .sort_values("MeanAbsSHAP", ascending=False, ignore_index=True)
    )
    top_features = importance_df["Feature"].head(n_top_features).tolist()

    # ---------------------- Plotting --------------------------------
    fig = None
    if plot:
        # --- take only top 3 features ---
        top3 = top_features
        fig, ax = plt.subplots(figsize=(8, max(3.5, 0.9 * len(top3))))
        ypos = np.arange(len(top3))

        try:
            from adjustText import adjust_text
            _have_adjust = True
        except Exception:
            _have_adjust = False
            adjust_text = None

        # color by absolute SHAP impact
        import matplotlib as mpl
        vmax = np.nanmax(np.abs(shap_2d)) if np.size(shap_2d) else 1.0
        norm = mpl.colors.Normalize(vmin=0.0, vmax=vmax)
        cmap = mpl.cm.viridis
        sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)

        texts_all = []
        arrows_all = []
        for i, feat in enumerate(top3):
            j = feature_names.index(feat)
            vals = shap_2d[:, j]

            # jitter for the beeswarm row
            jitter = (np.random.rand(len(vals)) - 0.5) * 0.15
            yrow = np.full_like(vals, ypos[i], dtype=float) + jitter

            # scatter colored by |SHAP|
            colors = sm.to_rgba(np.abs(vals))
            ax.scatter(vals, yrow, alpha=0.75, s=22, c=colors, edgecolor="none")

            n_labels = 5
            abs_sorted_idx = np.argsort(-np.abs(vals))  # sort descending by |SHAP|
            top_idx = abs_sorted_idx[:n_labels]

            for rank, idx in enumerate(top_idx):
                xt, yt = float(vals[idx]), float(yrow[idx])
                color = "tab:red" if vals[idx] > 0 else "tab:blue"

                t = ax.text(
                    xt, yt + 0.22, str(sample_names[idx]),
                    fontsize=8, color=color, ha="center", va="bottom",
                    bbox=dict(facecolor="white", alpha=0.85, edgecolor="none", pad=0.6),
                    zorder=5
                )
                texts_all.append(t)

                # arrow from label to point
                arr = ax.annotate(
                    "", xy=(xt, yt), xytext=(xt, yt + 0.18),
                    arrowprops=dict(arrowstyle="-", color="gray", lw=0.5)
                )
                arrows_all.append(arr)

        # axis / labels
        ax.axvline(0, color="black", linestyle="--", lw=1)
        ax.set_yticks(ypos)
        ax.set_yticklabels(top3)
        ax.set_xlabel("SHAP value (impact on model output)")
        ax.set_title(f"SHAP beeswarm — top 3 features ({target_name})")

        # colorbar for |impact|
        cbar = plt.colorbar(sm, ax=ax, pad=0.02)
        cbar.set_label(r"|SHAP impact|")

        plt.tight_layout()

        # Reduce overlap if adjustText is available
        if _have_adjust and texts_all:
            try:
                adjust_text(
                    texts_all,
                    ax=ax,
                    only_move={'points': 'y', 'text': 'y', 'objects': 'y'},
                    arrowprops=dict(arrowstyle="-", color="gray", lw=0.5),
                )
            except Exception:
                pass

        if dir is not None:
            os.makedirs(dir, exist_ok=True)
            fig.savefig(os.path.join(dir, "shap_beeswarm_top3.png"), dpi=300)

    return {
        "explainer": explainer,
        "shap_values": shap_2d,
        "importance_df": importance_df,
        "top_features": top_features,
        "figure": fig,
    }



def univariate_threshold_analysis(X: pd.DataFrame, y: pd.Series, thresholds_per_feature=100, plot_top_n=5):
    results = []
    feature_curves = {}

    for feature in X.columns:
        x_vals = X[feature].values
        thresholds = np.linspace(np.min(x_vals), np.max(x_vals), thresholds_per_feature)

        best_result = {
            'Feature': feature,
            'Best Threshold': None,
            'Accuracy': 0,
            'F1 Score': 0,
            'Direction': None
        }

        f1_scores_greater = []
        f1_scores_less = []

        for thresh in thresholds:
            # 'greater' direction
            y_pred_greater = (x_vals > thresh).astype(int)
            f1_greater = f1_score(y, y_pred_greater)
            f1_scores_greater.append(f1_greater)
            if f1_greater > best_result['F1 Score']:
                best_result.update({
                    'Best Threshold': thresh,
                    'Accuracy': accuracy_score(y, y_pred_greater),
                    'F1 Score': f1_greater,
                    'Direction': 'greater'
                })

            # 'less' direction
            y_pred_less = (x_vals < thresh).astype(int)
            f1_less = f1_score(y, y_pred_less)
            f1_scores_less.append(f1_less)
            if f1_less > best_result['F1 Score']:
                best_result.update({
                    'Best Threshold': thresh,
                    'Accuracy': accuracy_score(y, y_pred_less),
                    'F1 Score': f1_less,
                    'Direction': 'less'
                })

        results.append(best_result)
        feature_curves[feature] = {
            'thresholds': thresholds,
            'f1_greater': f1_scores_greater,
            'f1_less': f1_scores_less
        }

    result_df = pd.DataFrame(results).sort_values(by='F1 Score', ascending=False).reset_index(drop=True)

    # Plot top N features
    top_features = result_df.head(plot_top_n)['Feature']
    for feature in top_features:
        curve = feature_curves[feature]
        plt.figure(figsize=(8, 4))
        plt.plot(curve['thresholds'], curve['f1_greater'], label='greater than threshold', color='blue')
        plt.plot(curve['thresholds'], curve['f1_less'], label='less than threshold', color='red')
        plt.axvline(result_df[result_df['Feature'] == feature]['Best Threshold'].values[0], 
                    color='black', linestyle='--', label='Best threshold')
        plt.title(f'F1 Score vs. Threshold for {feature}')
        plt.xlabel('Threshold')
        plt.ylabel('F1 Score')
        plt.legend()
        plt.tight_layout()
        plt.show()

    return result_df

def ddG_to_ee_and_class(ddg_values, temperature=298.15, threshold_ee=35, return_df=True):
    """
    Converts ΔΔG values (kcal/mol) to enantiomeric excess (%ee), then binarizes.

    Parameters
    ----------
    ddg_values : array-like
        ΔΔG values in kcal/mol.
    temperature : float, optional
        Temperature in Kelvin. Default is 298.15 K.
    threshold_ee : float, optional
        Threshold (in %ee) to binarize outcome. Default is 50% ee.
    return_df : bool, optional
        If True, returns DataFrame with ΔΔG, %ee, and binary label. Else, returns only labels.

    Returns
    -------
    pd.DataFrame or np.ndarray
        If return_df=True, DataFrame with columns: ΔΔG, %ee, Binary Class.
        If return_df=False, just the binary class labels.
    """
    R = 0.0019872041  # kcal/mol·K
    ddg_values = np.asarray(ddg_values).astype(float)

    # Compute enantiomeric excess using tanh form
    ee = np.tanh(-ddg_values / (2 * R * temperature)) * 100  # in %

    # Binary classification: is abs(%ee) >= threshold?
    threshold_ee = float(threshold_ee)
    binary_class = (np.abs(ee) >= threshold_ee)

    if return_df:
        return pd.DataFrame({
            'ΔΔG (kcal/mol)': ddg_values,
            'Predicted ee (%)': ee,
            'Binary Class': binary_class
        })
    else:
        return binary_class
    
    
try:
    from ipywidgets import interact, FloatSlider, VBox
    _HAS_WIDGETS = True
except Exception:
    _HAS_WIDGETS = False


def interactive_corr_heatmap(
    data_or_corr: pd.DataFrame,
    *,
    initial_threshold: float = 0.9,
    title: str = "Correlation Heatmap",
    width: int = 700,
    height: int = 650,
    sort_by_strength: bool = True,
    dir=None
):
    """
    Show an interactive correlation heatmap with an adjustable |r| threshold.
    - Displays full corr (light gray) and highlights cells with |r| >= threshold in color.
    - Values are constrained to [-1, 1].
    - Works with either a raw DataFrame (computes df.corr()) or a correlation DataFrame.

    Parameters
    ----------
    data_or_corr : DataFrame
        Raw data (numeric) or a correlation matrix (square, symmetric, index=columns).
    initial_threshold : float
        Initial absolute correlation threshold to highlight.
    title : str
        Plot title.
    width, height : int
        Figure size in pixels.
    sort_by_strength : bool
        If True, reorder features by total |corr| strength for a cleaner view.
    """
    # ---- compute/validate corr
    corr = data_or_corr.corr() if not is_corr_matrix(data_or_corr) else data_or_corr.copy()
    corr = corr.clip(lower=-1, upper=1)  # keep within [-1, 1]
    labels = corr.columns.tolist()

    if sort_by_strength:
        order = corr.abs().sum().sort_values(ascending=False).index
        corr = corr.loc[order, order]
        labels = corr.index.tolist()

    # ---- base (light gray) layer
    def make_base_heatmap():
        return go.Heatmap(
            z=corr.values,
            x=labels, y=labels,
            colorscale=[(0.0, "#e0e0e0"), (1.0, "#7f7f7f")],  # gray background
            zmin=-1, zmax=1,
            showscale=False,
            hovertemplate="x=%{x}<br>y=%{y}<br>r=%{z:.3f}<extra></extra>",
            name="all",
            opacity=0.5,
        )

    # ---- highlighted layer for |r| >= threshold
    def make_masked_layer(th: float):
        z = corr.values.copy()
        mask = np.abs(z) < th
        z_masked = z.astype(float)
        z_masked[mask] = np.nan  # hide cells under threshold
        return go.Heatmap(
            z=z_masked,
            x=labels, y=labels,
            colorscale="RdBu",
            reversescale=True,
            zmin=-1, zmax=1,
            colorbar=dict(title="r"),
            hovertemplate="x=%{x}<br>y=%{y}<br>r=%{z:.3f}<extra></extra>",
            name=f"|r|≥{th:.2f}",
        )

    def make_figure(th: float) -> go.Figure:
        fig = go.Figure(data=[make_base_heatmap(), make_masked_layer(th)])
        fig.update_layout(
            title=f"{title} (|r| ≥ {th:.2f})",
            xaxis=dict(scaleanchor="y", constrain="domain", tickangle=45),
            yaxis=dict(autorange="reversed"),
            width=width, height=height,
            margin=dict(l=60, r=20, t=60, b=60),
        )
        return fig

    # ---- Jupyter: live slider; non-Jupyter: a figure with preset slider steps
    if _HAS_WIDGETS:
        slider = FloatSlider(
            value=float(initial_threshold),
            min=0.0, max=1.0, step=0.01,
            description="|r| ≥",
            continuous_update=False,
            readout_format=".2f",
        )

        @interact(threshold=slider)
        def _show(threshold):
            fig = make_figure(threshold)
            fig.show()

        return  # ipywidgets handles display

    # Fallback: Plotly steps (works anywhere, not as smooth but interactive)
    steps = []
    thresholds = [round(t, 2) for t in np.linspace(0.0, 1.0, 21)]
    fig = go.Figure()
    fig.add_trace(make_base_heatmap())             # trace 0
    fig.add_trace(make_masked_layer(initial_threshold))  # trace 1

    for t in thresholds:
        steps.append(dict(
            method="update",
            args=[
                {"visible": [True, True]},
                {"title": f"{title} (|r| ≥ {t:.2f})"}
            ],
            label=f"{t:.2f}",
        ))
    fig.update_layout(
        updatemenus=[dict(
            type="buttons",
            direction="right",
            x=0.5, xanchor="center", 
            y=1.12, yanchor="top",
            buttons=[
                dict(
                    label="Reset",
                    method="update",
                    args=[{}, {"title": title}]
                )
            ],
            showactive=False
        )],

        # Slider stays horizontal but positioned beside the plot (right side)
        sliders=[dict(
            active=int(initial_threshold * 20),
            steps=steps,
            x=1.05, xanchor="left",      # move to the right of the plot
            y=0.5,  yanchor="middle",    # vertically centered
            len=0.6,                     # slider length along the y-direction
            pad=dict(t=0, b=0),
            currentvalue=dict(prefix="|r| ≥ ", visible=True),
        )],

        title=f"{title} (|r| ≥ {initial_threshold:.2f})",
        xaxis=dict(
            scaleanchor="y",
            constrain="domain",
            tickangle=45,
        ),
        yaxis=dict(
            autorange="reversed",
        ),
        width=width, height=height,
        margin=dict(l=60, r=100, t=80, b=60),  # extra right margin for slider
    )
    fig.show()


def is_corr_matrix(df: pd.DataFrame) -> bool:
    """Heuristic: square, symmetric, and index==columns."""
    if not isinstance(df, pd.DataFrame):
        return False
    if df.shape[0] != df.shape[1]:
        return False
    if not df.index.equals(df.columns):
        return False
    # symmetry check with tolerance
    return np.allclose(df.values, df.values.T, atol=1e-12, equal_nan=True)


def interactive_corr_heatmap_with_highlights(df, initial_threshold=0.9, highlight_thresh=0.95, cell_size=40, min_size=400, max_size=1200):
    """
    Interactive Plotly heatmap of the correlation matrix with a slider for the minimum absolute correlation.
    Features that are in a highly correlated pair (|corr| > highlight_thresh) are marked with an asterisk (*).
    """
    corr = df.corr()
    order = corr.abs().sum().sort_values(ascending=False).index
    corr_sorted = corr.loc[order, order]
    n = len(order)
    size = min(max(n * cell_size, min_size), max_size)

    # Find features to highlight
    highly_corr_features = set()
    for i, feat_i in enumerate(order):
        for j, feat_j in enumerate(order):
            if i != j and abs(corr_sorted.iloc[i, j]) > highlight_thresh:
                highly_corr_features.add(feat_i)
                highly_corr_features.add(feat_j)
    # Add asterisk to highly correlated features
    labels = [
        f"{name}{'' if name in highly_corr_features else ''}"
        for name in order
    ]

    def plot(threshold):
        # Mask correlations below threshold and zeros
        corr_display = corr_sorted.where(corr_sorted.abs() >= threshold)
        corr_display = corr_display.where(corr_display != 0)

        fig = go.Figure(
            data=go.Heatmap(
                z=corr_display.values,
                x=labels,
                y=labels,
                colorscale='RdBu',
                zmid=0,
                colorbar=dict(title='Correlation'),
                hovertemplate='Feature 1: %{y}<br>Feature 2: %{x}<br>Correlation: %{z:.3f}<extra></extra>'
            )
        )
        fig.update_layout(
            title=f'Correlation Heatmap (|corr| ≥ {threshold})',
            xaxis_nticks=n,
            yaxis_nticks=n,
            autosize=False,
            width=size,
            height=size,
        )
        fig.show()

    interact(plot, threshold=FloatSlider(min=0, max=1, step=0.01, value=initial_threshold, description='Threshold'))

    
def threshold_analysis_plot(
    y_value,
    X_dict,
    types=None,
    feature_units=None,
    marker_shapes=None,
    target_names=None,
    cutoff="median",   # "median", "mean", float threshold, percentile 0–100, or tuple ("percentile", p, direction)
    data_labels=None,
    plot=True
):
    """
    Perform threshold-based classification analysis for selected features and optionally plot results.

    Parameters
    ----------
    y_value : pd.Series or np.ndarray
        Output values (e.g., yield, ΔΔG).
    X_dict : dict[str, array-like]
        Mapping {feature_name: values} of selected features.
    types : array-like or None
        Optional ligand types for plotting.
    feature_units : list[str] or None
        Units for axis labels.
    marker_shapes : dict or None
        Mapping from ligand_type to marker symbols.
    target_names : list or None
        Names for binary classes. Defaults to ["0","1"].
    cutoff : str | float | int | tuple
        Cutoff specification:
          - "median" → median of y
          - "mean" → mean of y
          - int/float 0–100 → percentile (high_is_positive)
          - float > 1 → absolute cutoff
          - tuple ("percentile", p, "low_is_positive"|"high_is_positive")
    plot : bool
        If True, generate contour + scatter plots.

    Returns
    -------
    results : list of dict
        Metrics and thresholds per feature.
    """

    # --- Compute cutoff ---
    direction = "high_is_positive"
    if isinstance(cutoff, tuple):
        if cutoff[0] != "percentile":
            raise ValueError("tuple cutoff must be ('percentile', p, direction)")
        p, direction = cutoff[1], cutoff[2]
        y_cut = np.percentile(y_value, p)
    elif cutoff == "median":
        y_cut = np.median(y_value)
    elif cutoff == "mean":
        y_cut = np.mean(y_value)
    elif isinstance(cutoff, (int, float)) and 0 <= cutoff <= 100:
        y_cut = np.percentile(y_value, cutoff)
    elif isinstance(cutoff, (int, float)):
        y_cut = float(cutoff)
    else:
        raise ValueError("cutoff must be 'median','mean', percentile, numeric, or tuple form")

    print(f"Using cutoff value = {y_cut:.3f}")

    # --- Default values ---
    if target_names is None:
        target_names = ["0", "1"]
    if feature_units is None:
        feature_units = ["" for _ in X_dict]

    results = []
    figs = []
    for i, (feature, X_all) in enumerate(X_dict.items()):
        y_all = np.asarray(y_value)

        # Optional ligand_type
        ligand_types = types if types is not None and len(types) == len(y_all) else None

        # Binarize by cutoff
        if direction == "low_is_positive":
            y_class = np.array([1 if v <= y_cut else 0 for v in y_all])
        else:  # high_is_positive
            y_class = np.array([1 if v >= y_cut else 0 for v in y_all])

        print(f"\nFeature: {feature}")
        print(f"Min={np.min(X_all):.2f}, Max={np.max(X_all):.2f}, Unique={len(np.unique(X_all))}")
        print(f"Class balance: {np.bincount(y_class)}")

        if len(np.unique(y_class)) < 2 or len(np.unique(X_all)) < 2:
            print("Not enough variation for classification.")
            continue

        # Fit decision stump
        dt = DecisionTreeClassifier(max_depth=1, class_weight="balanced", random_state=0)
        dt.fit(np.array(X_all).reshape(-1, 1), y_class)

        thr = dt.tree_.threshold[0]
        y_pred = dt.predict(np.array(X_all).reshape(-1, 1))

        acc = metrics.accuracy_score(y_class, y_pred)
        f1 = f1_score(y_class, y_pred)
        rec = metrics.recall_score(y_class, y_pred)

        print(f"Decision threshold on feature {feature}: {thr:.2f}")
        print(f"Accuracy={acc:.2f}, F1={f1:.2f}, Recall={rec:.2f}")
        print(metrics.classification_report(y_class, y_pred, target_names=target_names))

        results.append({
            "feature": feature,
            "decision_threshold": thr,
            "y_cut": y_cut,
            "direction": direction,
            "accuracy": acc,
            "f1": f1,
            "recall": rec
        })

        # Plot if requested
        if plot:
            x_min, x_max = np.min(X_all), np.max(X_all)
            y_min, y_max = np.min(y_all), np.max(y_all)
            dx, dy = x_max - x_min, y_max - y_min
            xx, yy = np.meshgrid(
                np.linspace(x_min - 0.05*dx, x_max + 0.05*dx, 200),
                np.linspace(y_min - 0.05*dy, y_max + 0.05*dy, 200)
            )
            Z = dt.predict(xx.ravel().reshape(-1, 1)).reshape(xx.shape)

            cMap_background = ListedColormap(["white", "#dfe0e0"])
            cMap_points = ListedColormap(["#c00", "#4caf50"])

            fig, ax = plt.subplots(figsize=(7, 6))   # <-- own the figure
            ax.contourf(xx, yy, Z, cmap=cMap_background, alpha=0.8)

            

            ax.axvline(x=thr, color="gray", linestyle=":", lw=1)
            ax.axhline(y=y_cut, color="gray", linestyle="--", lw=1)
            ax.text(
                1.02, 0.9, f"x-thr = {thr:.2f}",
                transform=ax.transAxes,        # both x,y in axes coords
                ha="left", va="center",
                fontsize=10, color="black",
                bbox=dict(facecolor="white", alpha=0.7, edgecolor="none", pad=1.0),
                clip_on=False
            )

            ax.text(
                1.02, 0.85, f"y-cut = {y_cut:.2f}",
                transform=ax.transAxes,        # both x,y in axes coords
                ha="left", va="center",
                fontsize=10, color="black",
                bbox=dict(facecolor="white", alpha=0.7, edgecolor="none", pad=1.0),
                clip_on=False
            )


            # Scatter points
            if ligand_types is not None:
                for lt in np.unique(ligand_types):
                    mask = ligand_types == lt
                    ax.scatter(
                        np.array(X_all)[mask], np.array(y_all)[mask],
                        c=y_class[mask], cmap=cMap_points, alpha=0.9, edgecolor="black", s=100,
                        marker=marker_shapes.get(lt, "o") if marker_shapes else "o", label=str(lt)
                    )
                ax.legend(title="Ligand Type", fontsize=12)
            else:
                ax.scatter(
                    X_all, y_all,
                    c=y_class, cmap=cMap_points, alpha=0.9, edgecolor="black", s=100, marker="o"
                )

            if data_labels is not None and len(data_labels) == len(y_all):
                texts = []
                for xi, yi, label in zip(X_all, y_all, data_labels):
                    texts.append(
                        ax.text(
                            xi, yi, str(label),
                            fontsize=8, color="black", ha="right", va="bottom"
                        )
                    )
                adjust_text(texts, arrowprops=dict(arrowstyle="->", color='gray', lw=0.5))

            ax.set_xlabel(f"{feature} {f'({feature_units[i]})' if feature_units[i] else ''}", fontsize=14)
            ax.set_ylabel("Y value", fontsize=14)
            ax.set_title(f"Threshold analysis: {feature}", fontsize=15)
            plt.show()

            # Prevent display in non-interactive environments
            figs.append(fig)   # <-- keep the open figure

            
    return results, figs


def y_randomization_check(model_class, X, y, n_runs=50, cv=5, random_state=42, **model_kwargs):
    """
    Perform Y-randomization (response permutation) test on a regression/classification model.

    Parameters
    ----------
    model_class : sklearn-like estimator class (e.g., LinearRegression, DecisionTreeClassifier)
    X : array-like, shape (n_samples, n_features)
        Feature matrix
    y : array-like, shape (n_samples,)
        Target values
    n_runs : int
        Number of randomization runs
    cv : int
        Cross-validation folds
    random_state : int
        Reproducibility
    model_kwargs : dict
        Extra keyword args passed to the model_class constructor

    Returns
    -------
    dict
        {
          "real_scores": [Q², MAE, RMSD],
          "random_scores": list of dicts for each run
        }
    """
    from sklearn.model_selection import cross_validate
    rng = np.random.default_rng(random_state)

    # --- Train with real y ---
    model = model_class(**model_kwargs)
    cv_results = cross_validate(
        model, X, y,
        cv=cv,
        scoring=("r2", "neg_mean_absolute_error", "neg_root_mean_squared_error"),
        return_train_score=False
    )
    real_scores = {
        "Q2": np.mean(cv_results["test_r2"]),
        "MAE": -np.mean(cv_results["test_neg_mean_absolute_error"]),
        "RMSD": -np.mean(cv_results["test_neg_root_mean_squared_error"])
    }

    # --- Y-randomization ---
    random_scores = []
    for run in range(n_runs):
        y_perm = rng.permutation(y)
        model = model_class(**model_kwargs)
        cv_results = cross_validate(
            model, X, y_perm,
            cv=cv,
            scoring=("r2", "neg_mean_absolute_error", "neg_root_mean_squared_error"),
            return_train_score=False
        )
        random_scores.append({
            "Q2": np.mean(cv_results["test_r2"]),
            "MAE": -np.mean(cv_results["test_neg_mean_absolute_error"]),
            "RMSD": -np.mean(cv_results["test_neg_root_mean_squared_error"])
        })

    return {"real_scores": real_scores, "random_scores": random_scores}