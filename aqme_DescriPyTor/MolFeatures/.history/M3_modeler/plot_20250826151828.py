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
import shap
from matplotlib.patches import Patch
from adjustText import adjust_text
from tkinter import ttk
import sys 
import os 
import math
import plotly.graph_objects as go
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from typing import Optional, List

import textwrap
from matplotlib.table import Table
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
try:
    from .modeling import fit_and_evaluate_single_combination_regression , fit_and_evaluate_single_combination_classification
    from .modeling_utils import _normalize_combination_to_columns
except ImportError as e:
    from modeling import fit_and_evaluate_single_combination_regression , fit_and_evaluate_single_combination_classification
    from modeling_utils import _normalize_combination_to_columns

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
                



def set_q2_plot_settings(ax, lower_bound, upper_bound, fontsize=15):
    bounds_array = np.array([lower_bound, upper_bound])
    ax.plot(bounds_array, bounds_array, 'k--', linewidth=2)  # black dashed line
    ax.set_xlabel('Measured', fontsize=fontsize)  # Assuming 'Measured' is the label you want
    ax.set_ylabel('Predicted', fontsize=fontsize)
    ax.set_ylim(bounds_array)
    ax.set_xlim(bounds_array)
    ax.grid(True)  # Adding a grid

def build_regression_equation(formula, coefficients, r_squared):
    """
    Build a regression equation string with proper LaTeX formatting.
    """
    intercept = getattr(coefficients, "iloc", coefficients)[0]  # Intercept
    feature_coeffs = coefficients[1:]

    # Escape underscores in feature names
    safe_formula = [str(name).replace("_", r"\_") for name in formula]

    equation_terms = []

    for i, coef in enumerate(feature_coeffs):
        sign = "+" if coef >= 0 else "-"
        equation_terms.append(f" {sign} {abs(coef):.2f}·{safe_formula[i]}")

    # Add intercept at the end
    sign_intercept = "+" if intercept >= 0 else "-"
    equation_terms.append(f" {sign_intercept} {abs(intercept):.2f}")

    # Build the LaTeX equation
    equation = f'$y = {"".join(equation_terms).strip()}$\n$R^2 = {r_squared:.2f}$'

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
    figsize=(16, 6),
    fontsize=12,
    scatter_color='#2ca02c',
    band_color='cadetblue',
    identity_color='#1f77b4',
    palette='deep',
    dpi=300,
    plot=True
):
    """
    Plots Predicted vs Measured with a smooth 90% confidence/prediction band
    around the regression line using seaborn's regplot, plus actual regression line,
    point labels, and Q² metrics.
    """
    y = np.asarray(y)
    y_pred = np.asarray(y_pred)
    labels = np.asarray(labels)

    data = pd.DataFrame({
    'Measured':  np.asarray(y, dtype=float).ravel(),
    'Predicted': np.asarray(y_pred, dtype=float).ravel(),
    'Labels':    np.asarray(labels, dtype=str).ravel()
})

    sns.set_theme(style='whitegrid', palette=palette)
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)

    # Scatter plot
    sns.scatterplot(
        data=data,
        x='Measured',
        y='Predicted',
        hue='Labels',
        palette=palette,
        edgecolor='w',
        s=100,
        ax=ax,
        legend=False
    )

    # Plot CI band only, no line (regplot)
    sns.regplot(
        data=data,
        x='Measured',
        y='Predicted',
        scatter=False,
        ci=90,
        line_kws={'color': band_color, 'linewidth': 1},
        ax=ax
    )

    # Range for plotting regression line
    mn = min(data['Measured'].min(), data['Predicted'].min())
    mx = max(data['Measured'].max(), data['Predicted'].max())
    x_ideal = np.linspace(mn, mx, 100)

    # Plot your regression line, dotted and clearly visible
    if coefficients is not None and len(coefficients) == 2:
        a, b = coefficients
        y_reg = a * x_ideal + b
        ax.plot(
            x_ideal, y_reg,
            linestyle='-',      # solid for visibility (change to ':' if you prefer)
            color='black',      # stands out
            linewidth=2.5,
            label='Regression Line'
        )

    # Set plot limits so the line is always visible
    ax.set_xlim(mn, mx)
    ax.set_ylim(mn, mx)

    # Equation & Pearson r
    corr = r if r is not None else np.corrcoef(y, y_pred)[0, 1]
    eqn = build_regression_equation(formula, coefficients, corr)
    ax.text(
        0.05, 0.95,
        f"{eqn}\nPearson r = {corr:.2f}",
        transform=ax.transAxes,
        fontsize=fontsize,
        va='top',
        bbox=dict(facecolor='white', alpha=0.8)
    )

    # Q² metrics
    if folds_df is not None and not folds_df.empty:
        q = folds_df.iloc[0]
        q_txt = (
            f"3-fold Q²: {q['Q2_3_Fold']:.2f}\n"
            f"5-fold Q²: {q['Q2_5_Fold']:.2f}\n"
            f"LOOCV Q²: {q['Q2_LOOCV']:.2f}"
        )
        ax.text(
            0.05, 0.80, q_txt,
            transform=ax.transAxes,
            fontsize=fontsize, va='top',
            bbox=dict(facecolor='white', alpha=0.8)
        )

    # Add point labels (optional, but can clutter if many points)
    texts = []
    for _, row in data.iterrows():
        texts.append(
            ax.text(row['Measured'], row['Predicted'], row['Labels'],
                    fontsize=fontsize-2, ha='center', va='bottom', color='gray')
        )
    adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle='-', color='gray', lw=0.5))

    ax.set_xlabel('Measured', fontsize=fontsize+2)
    ax.set_ylabel('Predicted', fontsize=fontsize+2)
    ax.set_title('Predicted vs Measured', fontsize=fontsize+4)
    # Show only your regression line in the legend
    handles, labels_ = ax.get_legend_handles_labels()
    reg_handles = [h for h, l in zip(handles, labels_) if l == 'Regression Line']
    reg_labels = [l for l in labels_ if l == 'Regression Line']
    valid = [lbl for lbl in labels if lbl and not lbl.startswith('_')]
    if reg_handles:
        ax.legend(reg_handles, reg_labels, loc='lower right', frameon=True)
    elif handles and valid:
        ax.legend().remove()  # Hide legend if only points
    else:
        leg = ax.get_legend()
    if leg:
        leg.remove()
    
    # plt.savefig(f'model_plot_{formula}.png', dpi=dpi)
    if plot:  
        plt.show()
   
    return fig



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
        ax.set_title(title, pad=8, fontsize=11, fontweight='bold')
    tbl = ax.table(cellText=df.values, colLabels=df.columns, loc='center', cellLoc='center')
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(fontsize)
    # Adaptive column widths
    ncols = len(df.columns)
    for (row, col), cell in tbl.get_celld().items():
        # header row
        if row == 0:
            cell.set_text_props(weight='bold')
            cell.set_height(0.08)
        else:
            cell.set_height(0.06)
        cell.set_edgecolor('#DDDDDD')
        if col < ncols:
            cell.set_linewidth(0.6)
    # squeeze to panel
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)

def _page_header(fig, subtitle):
    fig.suptitle(subtitle, y=0.98, fontsize=13, fontweight='bold')
    

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
    header_text_weight: str = "bold",
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
    ax.set_title(title, fontsize=12, fontweight="bold", pad=8)
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

def _save_top5_pdf(results: pd.DataFrame, model, pdf_path: str = "top_models_report.pdf", k: int = 5):
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
            fig.suptitle(f"Top {len(top_idx)} Models (ranked by R²)", y=0.98, fontsize=16, fontweight="bold")

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
            Q2_3, MAE_3, RMSD_3 = _try(lambda: model.calculate_q2_and_mae(X, y, n_splits=3),
                                       default=(np.nan, np.nan, np.nan), note="Q2/MAE/RMSD (3-fold) failed")
            Q2_5, MAE_5, RMSD_5 = _try(lambda: model.calculate_q2_and_mae(X, y, n_splits=5),
                                       default=(np.nan, np.nan, np.nan), note="Q2/MAE/RMSD (5-fold) failed")
            Q2_loo, MAE_loo, RMSD_loo = _try(lambda: model.calculate_q2_and_mae(X, y, n_splits=1),
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
                fig1.suptitle(f"{title}\n{formula_str}", y=0.98, fontsize=14, fontweight="bold")

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
            ax_meta.set_title("Model Summary", pad=8, fontsize=11, fontweight="bold")
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
                _ = generate_q2_scatter_plot(
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
                    fig_fallback.suptitle("Q² scatter", y=0.98, fontsize=14, fontweight="bold")
                ax = fig_fallback.add_subplot(111); ax.axis("off")
                ax.text(0.5, 0.5, "Q² scatter plot unavailable.", ha="center", va="center", fontsize=11, alpha=0.7)
                page_counter += 1
                _footer(fig_fallback, left_text="MolFeatures • Q² scatter", page_num=page_counter)
                pdf.savefig(fig_fallback, bbox_inches="tight"); plt.close(fig_fallback)

            # -------- PAGE 3+: SHAP analysis (optional) --------
            shap_res = _try(
                lambda: analyze_shap_values(
                    model=model,
                    X=model.features_df[feats],
                    feature_names=feats,
                    target_name=getattr(model, "output_name", "target"),
                    n_top_features=min(10, len(feats)),
                    plot=True
                ),
                default=None, note="SHAP analysis failed"
            )
            if shap_res and isinstance(shap_res, dict):
            
                for fig in shap_res.get("figures", []):
                    fig.text(0.01, 0.01, f"Rank #{rank} | SHAP", fontsize=8, ha="left", va="bottom", alpha=0.7)
                    page_counter += 1
                    _footer(fig, left_text="MolFeatures • SHAP", page_num=page_counter)
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
    try:
        _save_top5_pdf(results,model, pdf_path=f"{model.name}_top_models_report.pdf")
    except Exception as e:
        print(f"[PDF] Skipping top-5 export due to error: {e}")

    while True:

        if app:
            messagebox.showinfo('Models List:',df.to_markdown(index=False, tablefmt="pipe"))
            print(df.to_markdown(index=False, tablefmt="pipe"))
            selected_model = get_valid_integer('Select a model number: default is 0', 0)
            show_table_window('Models List:',df)
        else:
            print(df.to_markdown(index=False, tablefmt="pipe"))
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
        
        # Create a text file with the results
        with open('regression_results.txt', 'a') as f:
            f.write(f"Models list {df.to_markdown(index=False, tablefmt='pipe')} \n\n Model Coefficients\n\n{coef_df.to_markdown(tablefmt='pipe')}\n\n3-fold CV\n\n{pd.DataFrame({'Q2': [Q2_3], 'MAE': [MAE_3]}).to_markdown(tablefmt='pipe', index=False)}\n\n5-fold CV\n\n{pd.DataFrame({'Q2':[Q2_5], 'MAE': [MAE_5]}).to_markdown(tablefmt='pipe', index=False)}\n\n")
            print('Results saved to regression_results.txt in {}'.format(os.getcwd()))
        ## make a 3 5 loocv table to plot
        folds_df=pd.DataFrame({'Q2_3_Fold': [Q2_3], 'MAE': [MAE_3],'RMSD':[rmsd_3],'Q2_5_Fold':[Q2_5], 'MAE': [MAE_5],'RMSD':[rmsd_5],'Q2_LOOCV':[Q2_loo], 'MAE': [MAE_loo],'RMSD':[rmsd_loo]})
        r=r_squared[selected_model]
        # Generate and display the Q2 scatter plot
        
        _ = generate_q2_scatter_plot(y, pred, model.molecule_names,folds_df ,features,coef_df['Estimate'] ,r,X, lwr, upr, plot=False)

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



def generate_and_display_single_combination_plot(model, features, app=None):
    """
    Computes extra calculations (fitting, predictions, CV metrics, coefficient estimates, 
    and axis bounds) and then generates the Q2 scatter plot.

    Parameters:
        model: The regression model instance.
        features (list): List of feature names (columns in model.features_df) to use.
        app (optional): An application interface to display results (if provided).

    Returns:
        The return value of generate_q2_scatter_plot.
    """

    
    # Extract features and target values
    try:
        
        print("Extracting features from model.features_df...")
        X = model.features_df[features].to_numpy()
        y = model.target_vector.to_numpy()
       
    except Exception as e:
        print("Error extracting features:", e)
        return
    try:
        result = fit_and_evaluate_single_combination_regression(model, features)
        pred = result['predictions']
       
    except Exception as e:
        print("Error during model fitting/prediction:", e)
        return
    # Retrieve coefficient estimates
    # Compute cross-validation metrics
    try:
        print("Calculating cross-validation metrics for 3-fold CV...")
        Q2_3, MAE_3, rmsd_3 = model.calculate_q2_and_mae(X, y, n_splits=3)
        print("3-fold CV metrics: Q2: {}, MAE: {}, RMSD: {}".format(Q2_3, MAE_3, rmsd_3))
        
        print("Calculating cross-validation metrics for 5-fold CV...")
        Q2_5, MAE_5, rmsd_5 = model.calculate_q2_and_mae(X, y, n_splits=5)
        print("5-fold CV metrics: Q2: {}, MAE: {}, RMSD: {}".format(Q2_5, MAE_5, rmsd_5))
        
        print("Calculating cross-validation metrics for LOOCV...")
        Q2_loo, MAE_loo, rmsd_loo = model.calculate_q2_and_mae(X, y, n_splits=1)
        print("LOOCV metrics: Q2: {}, MAE: {}, RMSD: {}".format(Q2_loo, MAE_loo, rmsd_loo))


        leftout_pred = None
        # print models attributes and params
        
        try:
          
            if model.leftout_samples is not None and len(model.leftout_samples) > 0:
                print("Calculating left-out samples prediction and metrics...")
                X_left = model.leftout_samples[features]  # DataFrame shape (n_leftout, 4)
                X_left = X_left.reindex()
                y_left = model.leftout_target_vector  # Series shape (n_leftout,)
                # 3) call your predictor; let it add constant & reorder itself
                try:
                    leftout_pred = model.predict_for_leftout(X_left, y=y_left, calc_interval=False)
                    print("Left-out samples prediction completed.", leftout_pred)
                    if isinstance(leftout_pred, tuple):
                        y_pred = np.array(leftout_pred[0]).ravel()  # Use only predictions, not errors
                    else:
                        y_pred = np.array(leftout_pred).ravel()
                        
                    print(f"Successfully predicted left-out samples: {leftout_pred}")
                except Exception as e:
                    print("Error predicting left-out samples:", e)
                    leftout_pred = None  # Ensure variable exists if exception occurs

                if leftout_pred is not None:
                    y_true = np.array(model.leftout_target_vector).ravel()
                    names  = list(model.molecule_names_predict)

                    # 2) build DataFrame
                    prediction_df = pd.DataFrame({
                        'Molecule':   names,
                        'Actual':     y_true,
                        'Predicted':  y_pred
                    })

                    # 3) compute absolute percent error (always positive, avoids division by zero)
                    prediction_df['Error in %'] = np.where(
                        prediction_df['Actual'] != 0,
                        np.abs(prediction_df['Actual'] - prediction_df['Predicted']) / np.abs(prediction_df['Actual']) * 100,
                        np.nan
                    )

                    print(prediction_df)
                    # Calculate and print R2 as well
                    r2_leftout = r2_score(y_true, y_pred)
                    mae_leftout = mean_absolute_error(y_true, y_pred)
                    print(f"R² for left-out predictions: {r2_leftout:.4f}")
                    print(f"MAE for left-out predictions: {mae_leftout:.4f}")

                        
                else:
                    print("No left-out predictions available; skipping result table.")
        except Exception as e:
            print("Error:", e)
            print("No left-out samples available; skipping result table.")
            pass

        # Prepare a folds DataFrame with CV results
        folds_df = pd.DataFrame({
            'Q2_3_Fold': [Q2_3],
            'MAE_3': [MAE_3],
            'RMSD_3': [rmsd_3],
            'Q2_5_Fold': [Q2_5],
            'MAE_5': [MAE_5],
            'RMSD_5': [rmsd_5],
            'Q2_LOOCV': [Q2_loo],
            'MAE_LOOCV': [MAE_loo],
            'RMSD_LOOCV': [rmsd_loo]
        })
        print("Folds DataFrame prepared:")
        # print(folds_df)
    except Exception as e:
        
        print(traceback.format_exc())
        print("Error calculating cross-validation metrics:", e)
        return

    # Compute axis bounds for plotting
    try:
        print("Calculating fixed margin lines and axis bounds…")
        # Set axis bounds
        plot_scale_start = -3.0
        plot_scale_end = 0.6
        x_ideal = np.linspace(plot_scale_start, plot_scale_end, 100)

        # Fixed margin lines (identity ±0.25, ±0.50)
        margin_lines = {
            "identity":     x_ideal,                           # y = x
            "+0.25":        x_ideal + 0.25,                    # y = x + 0.25
            "-0.25":        x_ideal - 0.25,                    # y = x - 0.25
            "+0.50":        x_ideal + 0.50,                    # y = x + 0.50
            "-0.50":        x_ideal - 0.50                     # y = x - 0.50
        }

    except Exception as e:
        print("Error calculating margin lines:", e)
        return

    try:
        vif_df = model._compute_vif(model.features_df[features])
        print_models_vif_table(vif_df)
    except Exception as e:
        print("Error calculating VIF:", e)
        return

    # Compute R^2 (or r-squared) as a measure of correlation
    try:
        print("Calculating R^2 value...")
        r = np.corrcoef(y, pred)[0, 1]**2
        mae = mean_absolute_error(y, pred)
        print(f"R^2 value: {r:.4f}, MAE: {mae:.4f}")
      
    except Exception as e:
        print("Error calculating R^2:", e)
        return
    
    try:
        y_i,upr,lwr = model.predict(X, return_interval=True)
        print("Retrieving coefficient estimates...")
        coef_df = model.get_covariance_matrix(features)
        print("Coefficient estimates retrieved:")
        print(coef_df.head())
    except Exception as e:
        print("Error retrieving coefficient estimates:", e)
        return
    
    try:
        print("Calling generate_q2_scatter_plot with computed values...")
        # Remove pi_lower and pi_upper if you are not calculating prediction intervals
        plot_output = generate_q2_scatter_plot(
            y, pred, model.molecule_names, folds_df, features, coef_df['Estimate'], r,plot=True
        )
        X=model.features_df[features]
        shap_plot = model.plot_shap_values(X, plot=True)
        print("Plot generated successfully.")
    except Exception as e:
        print("Error in generate_q2_scatter_plot:", e)
        return

    return


def plot_feature_vs_target(feature_values, y_values, feature_name, y_name="Target", point_labels=None, figsize=(10, 6)):
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
    return fig

def plot_all_features_vs_target(features_df, target_vector, molecule_names=None, figsize=(12, 10)):
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



def analyze_shap_values(model, X, feature_names=None, target_name="output",
                        n_top_features=10, plot=False):
    """
    Compute SHAP analysis and RETURN results (figures optional).
    Guaranteed to return explicit Matplotlib Figure handles when plot=True.
    """
    import numpy as np
    import pandas as pd
    import shap
    import matplotlib.pyplot as plt

    # ---- feature names
    if hasattr(X, "columns"):  # DataFrame
        feature_names = X.columns.tolist()
        X = X.values  # keep the array form
    elif hasattr(X, "shape"):  # NumPy array
        feature_names = [f"Feature {i}" for i in range(X.shape[1])]
    elif isinstance(X, list):
        if len(X) > 0 and isinstance(X[0], (list, tuple)):
            feature_names = [f"Feature {i}" for i in range(len(X[0]))]
        else:  # flat list
            feature_names = ["Feature 0"]
            X = np.array(X).reshape(-1, 1)  # force it into 2D
    else:
        raise ValueError(f"Unsupported type for X: {type(X)}")

    results = {}

    # ---- explainer
    try:
        explainer = shap.Explainer(model)
    except Exception:
        background = shap.sample(X, min(50, X.shape[0])) if isinstance(X, pd.DataFrame) else X
        background = background if isinstance(background, np.ndarray) else background.values
        explainer = shap.KernelExplainer(model.predict, background)

    # ---- SHAP values
    shap_values = explainer.shap_values(X)
    if isinstance(shap_values, list):
        if len(shap_values) > 1:
            mean_shap = np.abs(np.array(shap_values)).mean(axis=0)
            results['shap_values'] = shap_values
        else:
            mean_shap = np.abs(shap_values[0])
            results['shap_values'] = shap_values[0]
    else:
        mean_shap = np.abs(shap_values)
        results['shap_values'] = shap_values
    results['mean_shap'] = mean_shap

    # ---- importance
    feature_importance = pd.DataFrame({
        'Feature': feature_names,
        'Mean Absolute SHAP Value': np.mean(np.abs(mean_shap), axis=0),
        'Max Absolute SHAP Value':  np.max(np.abs(mean_shap), axis=0)
    }).sort_values('Mean Absolute SHAP Value', ascending=False)
    results['feature_importance'] = feature_importance
    results['top_features'] = feature_importance['Feature'].head(n_top_features).tolist()

    # ---- figures (explicit creation)
    results['figures'] = []
    results['fig_summary'] = None

    if plot:
        # Create a NEW figure explicitly and draw into it
        fig = plt.figure(figsize=(9, 6))
        ax = fig.add_subplot(111)

        # shap.summary_plot draws on the current figure; show=False prevents blocking
        # Use the raw matrix (values) if X is a DataFrame
        X_mat = X.values if isinstance(X, pd.DataFrame) else X
        shap.summary_plot(
            results['shap_values'],
            X_mat,
            feature_names=feature_names,
            show=False,
            max_display=n_top_features
        )

        # Tighten and hand back the figure
        try:
            fig.tight_layout()
        except Exception:
            pass

        results['fig_summary'] = fig
        results['figures'].append(fig)
        
    return results




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
            x=0.5, xanchor="center", y=1.12, yanchor="top",
            buttons=[
                dict(
                    label="Reset",
                    method="update",
                    args=[{}, {"title": title}]
                )
            ],
            showactive=False
        )],
        sliders=[dict(
            active=int(initial_threshold*20),
            steps=steps,
            x=0.5, xanchor="center", y=-0.08, len=0.8,
            currentvalue=dict(prefix="|r| ≥ ", visible=True)
        )],
        title=f"{title} (|r| ≥ {initial_threshold:.2f})",
        xaxis=dict(scaleanchor="y", constrain="domain", tickangle=45),
        yaxis=dict(autorange="reversed"),
        width=width, height=height,
        margin=dict(l=60, r=20, t=60, b=80),
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