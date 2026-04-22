import sys
import os
if not os.environ.get("MPLBACKEND"):
    os.environ["MPLBACKEND"] = "Agg"
import gc
import time
import psutil
import argparse
from typing import Optional, List
try:
    from tkinter import Tk, Frame, Label, Button, Entry, StringVar, OptionMenu, Toplevel, filedialog, Text, Scrollbar, Checkbutton, IntVar, Canvas, LEFT, SOLID,END
    import customtkinter  # Assuming you have this library
    
    pkg_dir = os.path.dirname(os.path.abspath(__file__))
    # 2) Compute your project root (the parent of MolFeatures/)
    project_root = os.path.dirname(pkg_dir)
    project_root = os.path.join(project_root, 'MolFeatures')
    print(f'Starting at {project_root}')
    # 3) If it’s not already on sys.path, insert it at the front
    if project_root not in sys.path:
        sys.path.insert(0, project_root)
        
    import matplotlib
    matplotlib.use(os.environ["MPLBACKEND"])
    from tkinter.simpledialog import askstring
    import csv
    import pandas as pd
    import re
    import shutil
    import subprocess
    import warnings
    from tkinter import filedialog, messagebox
    import warnings
    from typing import Dict
    from tkinter import ttk
    from PIL import Image , ImageTk
    from datetime import datetime
    import argparse
    import webbrowser
    import seaborn as sns
    import json
    import inspect
    import traceback
    import urllib.request

    try:
        from .utils import help_functions, file_handlers
        from .M2_data_extractor.data_extractor import Molecules
        from .M1_pre_calculations.main import Module1Handler
        # from .Mol_align.renumbering import batch_renumbering
        from .M2_data_extractor.feather_extractor import logs_to_feather
        from .M3_modeler.modeling import ClassificationModel, LinearRegressionModel
        from .M2_data_extractor.cube_sterimol import cube_many
        from .M2_data_extractor.sterimol_standalone import Molecules_xyz
    except ImportError:

        import utils.help_functions as help_functions
        import utils.file_handlers as file_handlers
        from M2_data_extractor.data_extractor import Molecules
        from M1_pre_calculations.main import Module1Handler
        # from Mol_align.renumbering import batch_renumbering
        from M2_data_extractor.feather_extractor import logs_to_feather
        from M3_modeler.modeling import ClassificationModel, LinearRegressionModel
        from M2_data_extractor.cube_sterimol import cube_many
        from M2_data_extractor.sterimol_standalone import Molecules_xyz
        




except ImportError or ModuleNotFoundError as e:
    print(f"An error occurred: {e}, Run install_packages.py script to install the required packages.")


def load_answers_json(path = None):

        if path is not None:
            file_path = path
        else:
            file_path = filedialog.askopenfilename(
                defaultextension=".json",
                filetypes=[("JSON files", "*.json"), ("All files", "*.*")]
            )
        if not file_path:
            return {}
        with open(file_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
        print(f'printing loaded data from json file: {data}')
        return data


def get_local_setup_version():
    """
    Look for a line like `version='0.9009'` in the sibling setup.py
    and return the version string. Falls back to 'unknown'.
    """
    here = os.path.abspath(os.path.dirname(__file__))
    setup_path = os.path.join(here, 'setup.py')
    try:
        version_pattern = re.compile(r"version\s*=\s*['\"]([^'\"]+)['\"]")
        with open(setup_path, encoding='utf-8') as f:
            for line in f:
                match = version_pattern.search(line)
                if match:
                    return match.group(1)
    except FileNotFoundError:
        pass
    
    return "unknown"

__version__ = get_local_setup_version()

# Function to get the current date
def get_current_date():
    return datetime.now().strftime("%Y-%m-%d")

def show_results(message):
        messagebox.showinfo("Results", message)
            
class ToolTip(object):
    def __init__(self, widget):
        self.widget = widget
        self.tipwindow = None
        self.id = None
        self.x = self.y = 0

    def showtip(self, text):
        "Display text in tooltip window"
        self.text = text
        if self.tipwindow or not self.text:
            return
        x, y, _, _ = self.widget.bbox("insert")
        x = x + self.widget.winfo_rootx() + 25
        y = y + self.widget.winfo_rooty() + 25
        self.tipwindow = tw = Toplevel(self.widget)
        tw.wm_overrideredirect(True)
        tw.wm_geometry(f"+{x}+{y}")
        label = Label(tw, text=self.text, justify=LEFT,
                         background="#ffffe0", relief=SOLID, borderwidth=1,
                         font=("tahoma", "8", "normal"))
        label.pack(ipadx=1)

    def hidetip(self):
        tw = self.tipwindow
        self.tipwindow = None
        if tw:
            tw.destroy()

def createToolTip(widget, text):
    toolTip = ToolTip(widget)
    def enter(event):
        toolTip.showtip(text)
    def leave(event):
        toolTip.hidetip()
    widget.bind('<Enter>', enter)
    widget.bind('<Leave>', leave)

def convert_to_custom_nested_list(input_str):
    """
    Converts a comma-separated string into a nested list based on the format:
    - "1,2,3,4,5,6" -> [[1,2,3,4], 5, 6]
    - "1,2,3,4,5,6 2,3,1,5,6,7" -> [[[1,2,3,4], 5, 6], [[2,3,1,5], 6, 7]]
    """
   
    split_by_space = input_str.split(' ')  # Split by space for multiple sections

    def process_sublist(sublist_str):
        """Convert a single comma-separated list to the required nested structure."""
        elements = list(map(int, sublist_str.split(',')))  # Convert to list of integers
        if len(elements) > 2:  # Ensure we separate all but the last two elements
            return [elements[:-2]] + elements[-2:]
        return elements  # If fewer than 3 elements, return as-is

    # Process each segment and decide if it's a nested list or a single flat list
    if len(split_by_space) == 1:
        # Single segment, no spaces
        return process_sublist(split_by_space[0])
    else:
        # Multiple segments separated by spaces
        nested_list = []
        for sublist_str in split_by_space:
            nested_list.append(process_sublist(sublist_str))
        return nested_list


def convert_to_list_or_nested_list(input_str):
    
    split_by_space = input_str.split(' ')

    # If there are no spaces, return a flat list
    if len(split_by_space) == 1:
        try:
            return list(map(int, split_by_space[0].split(',')))
        except:
            return list(map(int, split_by_space[0]))
    
    # Otherwise, return a nested list
    nested_list = []
    for sublist_str in split_by_space:
        sublist = list(map(int, sublist_str.split(',')))
        nested_list.append(sublist)

    if nested_list:
        return nested_list
    else:
        return None
   

class MoleculeApp:
    def __init__(self, master):
        
        self.master = master
        master.title("Molecule Data Extractor - version {}".format(__version__))
        self.current_file_path = os.path.abspath(__file__)
        # Get the directory of the current file
        self.current_directory = os.path.dirname(self.current_file_path)
        os.chdir(self.current_directory)
        


        self.sidebar_frame_left = customtkinter.CTkFrame(master, width=100,height=250, corner_radius=0)
        self.sidebar_frame_left.grid(row=0, column=0, rowspan=4, sticky="nsew")

        self.sidebar_frame_right = customtkinter.CTkFrame(master, width=140,height=250, corner_radius=0)
        self.sidebar_frame_right.grid(row=0, column=3, rowspan=4, sticky="nsew")
        

        
        self.output_text = Text(master, wrap="word", height=30, width=100)
        self.v_scrollbar = Scrollbar(master, command=self.output_text.yview)
        # self.h_scrollbar = Scrollbar(master, orient='horizontal', command=self.output_text.xview)

        self.output_text.config(yscrollcommand=self.v_scrollbar.set)  # , xscrollcommand=self.h_scrollbar.set

        # Place the Text widget and scrollbars in the grid
        self.output_text.grid(row=0, column=1, rowspan=4, sticky="nsew")
        self.v_scrollbar.grid(row=0, column=2, rowspan=4, sticky='ns')
        # self.h_scrollbar.grid(row=4, column=1, columnspan=2, sticky='ew')

        self.output_text.bind(
            "<Enter>",
            lambda e: e.widget.focus_set()
        )

        # 2) Bind wheel events to all Text widgets via the Text *class*:
        self.master.bind_class(
            "Text",
            "<MouseWheel>",
            lambda e: e.widget.yview_scroll(int(-1*(e.delta/120)), "units")
        )
        # Linux variants, if you care:
        self.master.bind_class("Text", "<Button-4>", lambda e: e.widget.yview_scroll(-1, "units"))
        self.master.bind_class("Text", "<Button-5>", lambda e: e.widget.yview_scroll( 1, "units"))
        self.show_result(f"Current directory: {self.current_directory}\n List of files: {os.listdir()}\n")
        self.print_description()
        
        #label for parameters
        self.param_description = Label(master, text="")
        self.param_description.grid(row=4, column=1, sticky='w')

        # Entry for parameters
        self.param_entry = Entry(master,width=50)
        self.param_entry.grid(row=5, column=1, sticky='w')

        # Submit button
        self.submit_button = Button(master, text="Submit", command=self.activate_method)
        self.submit_button.grid(row=5, column=0, sticky='e')

        self.label = customtkinter.CTkLabel(self.sidebar_frame_left, text="Choose Directory to Load Feather files:")
        self.label.grid(row=0, column=0, padx=20, pady=10)

        self.folder_path = StringVar()
        self.browse_button = customtkinter.CTkButton(self.sidebar_frame_left, text="Browse Feather Directory", command=self.browse_directory)
        self.browse_button.grid(row=1, column=0, padx=20, pady=10)
        createToolTip(self.browse_button, "Choose the directory where the Feather files are located to initialize them as Molecule objects.")
        
       
        self.method_var = StringVar(master)
        self.method_var.set("Choose a method")  # Default value
        self.method_var.trace_add("write", lambda *args: self.open_param_window())

        self.method_menu = OptionMenu(self.sidebar_frame_left, self.method_var, "Windows Command",
                                    "get_sterimol_dict", "get_npa_dict", "get_stretch_dict", "get_ring_dict",
                                    "get_dipole_dict", "get_bond_angle_dict", "get_bond_length_dict",
                                    "get_charge_dict", "get_charge_diff_dict", "get_bending_dict")
        self.method_menu.grid(row=3, column=0, padx=20, pady=10)
        
        # StringVar for dropdown menu selection
        self.file_handler_var = StringVar(master)
        self.file_handler_var.set("File Handler")  # Default value
    
        # Dropdown menu for file handling options
        self.file_handler_menu = OptionMenu(self.sidebar_frame_left, self.file_handler_var, "Smiles to XYZ", "Create com Files", "Log to Feather")
        self.file_handler_menu.grid(row=3, column=1, padx=20, pady=10)
        self.file_handler_var.trace_add("write", lambda *args: self.handle_file_action())

        # Separate button for Visualization
        self.visualize_button = customtkinter.CTkButton(self.sidebar_frame_left, text="Visualize Molecules", command=self.filter_molecules_vis)
        self.visualize_button.grid(row=5, column=0, padx=20, pady=10)
        createToolTip(self.visualize_button, "Select Molecules to visualize.")

        self.model_button = customtkinter.CTkButton(self.sidebar_frame_left, text="Model Data", command=self.run_model_in_directory)
        self.model_button.grid(row=4, column=1, padx=20, pady=10)
        createToolTip(self.model_button, "Run a model in a specified directory using provided CSV filepaths.\n Choose between classification and regression \n Choose between 2-4 features and provide a features with target CSV file.")

        # Separate button for Export Data
        self.export_button = customtkinter.CTkButton(self.sidebar_frame_left, text="Extract DataFrames", command=self.export_data)
        self.export_button.grid(row=6, column=0, padx=20, pady=10)
        self.molecules = None  # Placeholder for Molecules object
        createToolTip(self.export_button, "Extract DataFrames from the Molecules object and save them to CSV files in a directory for each molecule.")
        
        self.export_button = customtkinter.CTkButton(self.sidebar_frame_left, text="Extract xyz files", command=self.export_xyz)
        self.export_button.grid(row=7, column=0, padx=20, pady=10)
        createToolTip(self.export_button, "Extract  xyz files from the Molecules object to a directory called xyz_files.")

        self.filter_molecules_button = customtkinter.CTkButton(self.sidebar_frame_left, text="Filter Molecules", command=self.filter_molecules)
        self.filter_molecules_button.grid(row=1, column=1, padx=20, pady=10)
        createToolTip(self.filter_molecules_button, "Open a window to select specific molecules out of the Loaded.")

        self.comp_set_button = customtkinter.CTkButton(self.sidebar_frame_left, text="Extract Features", command=self.open_question_window)
        self.comp_set_button.grid(row=4, column=0, padx=20, pady=10)
        createToolTip(self.comp_set_button, "Open a window to input parameters for the Extract Features.\n you can load a previous input file to save time.")
         
        if self.molecules is not None:
            self.check_vars = [IntVar(value=1) for _ in self.molecules_names]
        else:
            self.check_vars = []

        # save text button
        self.save_txt_button = Button(master, text="Save Output Text", command=self.save_text)
        self.save_txt_button.grid(row=5, column=1, sticky='e')
        createToolTip(self.save_txt_button, "Save the output text from the box to a .txt file.")
        self.master.report_callback_exception = self._show_error

        restart_button = customtkinter.CTkButton(
            master=self.sidebar_frame_right,
            text="Restart App",
            command=self.restart_app
        )
        restart_button.pack(pady=10)

    def _show_error(self, exc, val, tb):
        """Display a traceback in a pop-up, then let the app continue."""
        err = "".join(traceback.format_exception(exc, val, tb))
        messagebox.showerror("Unexpected Error", err)

    
    def restart_app(self):
        # 1) Path to the Python executable
        python = sys.executable

        # 2) Compute the project root (one level above MolFeatures/)
        #    __file__ lives in MolFeatures/<some_module>.py
        project_root = os.path.dirname(os.path.dirname(__file__))

        # 3) Switch cwd to the project root so -m MolFeatures works
        os.chdir(project_root)

        # 4) Re-invoke exactly: python -m MolFeatures gui
        os.execv(python, [python, "-m", "MolFeatures", "gui"])
            
    def print_description(self):
        # Path to description.txt file
        txt_path = 'description.txt'
        try:
            with open(txt_path, 'r') as txt_file:
                string = txt_file.read()
                # Show the description text
                self.show_result(string)
                
                # Insert a clickable HTML link
                self.output_text.insert(END, "\nFor more information, visit: ")
                self.output_text.insert(END, "Click Here", "link")

                # Bind the link to open a web browser
                self.output_text.tag_config("link", foreground="blue", underline=1)
                self.output_text.tag_bind("link", "<Button-1>", lambda e: webbrowser.open_new("https://github.com/edenspec2/LabCode"))

        except FileNotFoundError:
            self.show_result("Description file not found.")

    

    def validate_integer(self, P):
            # Validation function to allow only integer input
            if P.isdigit() or P == "":
                return True
            return False
    
    def leave_out_molecules(self):
        """
        Opens a new Tkinter window with a list of molecules as checkboxes.
        Returns the list of indices for any molecules the user checks.
        """
        names = self.molecules_csv_names  # Assuming this is a list of molecule names
     
        # Create a new Toplevel window
        new_window = Toplevel(self.master)
        new_window.title("Select Molecules to Remove")

        # We'll store the (index, IntVar) for each molecule in a list
        var_list = []
        for i, name in enumerate(names):
            var = IntVar()
            chk = Checkbutton(new_window, text=name, variable=var)
            chk.pack(anchor="w", padx=10, pady=2)
            var_list.append((i, var))

        # We'll store selected indices here after the user clicks 'Confirm'
        selected_indices = []

        def confirm_selection():
            """Collect all checked molecules and close the window."""
            for idx, variable in var_list:
                if variable.get() == 1:
                    selected_indices.append(idx)
            
            self.leave_out_indices = selected_indices

            # 2) Also store as a comma-separated string in the StringVar
            if selected_indices:
                csv_str = ",".join(str(x) for x in selected_indices)
                self.leave_out_mols.set(csv_str)
            else:
                self.leave_out_mols.set("None")
            new_window.destroy()  # Closes the Toplevel window
            

        # Add a 'Confirm' button at the bottom
        confirm_button = Button(new_window, text="Confirm", command=confirm_selection)
        confirm_button.pack(pady=10)

        # This forces the function to wait until the user closes the window,
        # allowing us to return the final list of selected indices.
        new_window.wait_window()

    
    def run_model_in_directory(self):
        new_window = Toplevel(self.master)
        new_window.title("Run Model")
        
        
        
        vcmd = (new_window.register(self.validate_integer), '%P')
        # Model Type Selection
        model_type_label = Label(new_window, text="Select Model Type:")
        model_type_label.grid(row=0, column=0, padx=10, pady=10)
        
        model_type_var = StringVar(new_window)
        model_type_dropdown = OptionMenu(new_window, model_type_var, 'classification', 'linear_regression')
        model_type_dropdown.grid(row=0, column=1, padx=10, pady=10)
        # 2) SUBTYPE (only relevant if linear_regression is selected): ordinary vs. lasso

        regression_subtype_label = Label(new_window, text="Regression Subtype:")
        regression_subtype_var = StringVar(new_window, value='ordinary')
        regression_subtype_dropdown = OptionMenu(new_window, regression_subtype_var, 'ordinary', 'lasso')

        # 3) ALPHA ENTRY (only for lasso)
        lasso_alpha_label = Label(new_window, text="Lasso Alpha:")
        lasso_alpha_var = StringVar(new_window, value='0.1')
        lasso_alpha_entry = Entry(new_window, textvariable=lasso_alpha_var, width=5)

        def update_model_type(*args):
            """
            Show or hide the regression subtype and alpha entry
            depending on whether model_type_var is 'linear_regression'.
            """
            if model_type_var.get() == 'linear_regression':
                # Show the regression subtype widgets
                regression_subtype_label.grid(row=0, column=2, padx=10, pady=10)
                regression_subtype_dropdown.grid(row=0, column=3, padx=10, pady=10)

                # Show or hide alpha label/entry based on current subtype
                if regression_subtype_var.get() == 'lasso':
                    lasso_alpha_label.grid(row=0, column=4, padx=10, pady=10)
                    lasso_alpha_entry.grid(row=0, column=5, padx=10, pady=10)
                else:
                    lasso_alpha_label.grid_remove()
                    lasso_alpha_entry.grid_remove()
            else:
                # Hide them all if 'classification'
                regression_subtype_label.grid_remove()
                regression_subtype_dropdown.grid_remove()
                lasso_alpha_label.grid_remove()
                lasso_alpha_entry.grid_remove()

        def update_subtype(*args):
            """
            Show or hide lasso alpha depending on whether 
            the subtype is 'lasso'.
            """
            # Only show alpha if we're in linear_regression *and* 'lasso'
            if model_type_var.get() == 'linear_regression' and regression_subtype_var.get() == 'lasso':
                lasso_alpha_label.grid(row=0, column=4, padx=10, pady=10)
                lasso_alpha_entry.grid(row=0, column=5, padx=10, pady=10)
            else:
                lasso_alpha_label.grid_remove()
                lasso_alpha_entry.grid_remove()

        # Trace both variables so UI updates any time the user changes a dropdown
        model_type_var.trace_add("write", update_model_type)
        regression_subtype_var.trace_add("write", update_subtype)

        # Initialize UI state
        update_model_type()
        update_subtype()


        # Feature CSV Selection
        def select_features_csv():
            filepath = filedialog.askopenfilename(defaultextension=".csv",
                                                  filetypes=[("CSV files", "*.csv"), ("All files", "*.*")])
            feature_csv_entry.delete(0, END)
            feature_csv_entry.insert(0, filepath)
            self.molecules_csv_names = pd.read_csv(filepath)['Unnamed: 0'].tolist()
            

        features_csv_label = Label(new_window, text="Select Features CSV:")
        features_csv_label.grid(row=1, column=0, padx=10, pady=10)

        feature_csv_entry = Entry(new_window, width=50)
        feature_csv_entry.grid(row=1, column=1, padx=10, pady=10)

        feature_csv_button = Button(new_window, text="Browse...", command=select_features_csv)
        feature_csv_button.grid(row=1, column=2, padx=10, pady=10)

        # Target CSV Selection (optional)
        def select_target_csv():
            filepath = filedialog.askopenfilename(defaultextension=".csv",
                                                  filetypes=[("CSV files", "*.csv"), ("All files", "*.*")])
            target_csv_entry.delete(0, END)
            target_csv_entry.insert(0, filepath)

        target_csv_label = Label(new_window, text="Select Target CSV (Optional):")
        target_csv_label.grid(row=2, column=0, padx=10, pady=10)

        target_csv_entry = Entry(new_window, width=50)
        target_csv_entry.grid(row=2, column=1, padx=10, pady=10)

        target_csv_button = Button(new_window, text="Browse...", command=select_target_csv)
        target_csv_button.grid(row=2, column=2, padx=10, pady=10)

        leave_out_button=Button(new_window, text="Show Molecules Numbers", command=self.leave_out_molecules)
        leave_out_button.grid(row=3, column=2, padx=10, pady=10)
        Label(new_window, text="Leave Out Molecules:").grid(row=3, column=0)
        self.leave_out_mols=StringVar(new_window,value='None')
        Entry(new_window,textvariable=self.leave_out_mols, validate='key').grid(row=3, column=1)
        # self.leave_out_mols=convert_to_list_or_nested_list(self.leave_out_mols)

        Label(new_window, text="Min Features:").grid(row=4, column=0)
        self.min_features_var = StringVar(new_window, value='2')
        Entry(new_window, textvariable=self.min_features_var, validate='key', validatecommand=vcmd).grid(row=4, column=1)

        Label(new_window, text="Max Features:").grid(row=5, column=0)
        self.max_features_var = StringVar(new_window, value='None')
        Entry(new_window, textvariable=self.max_features_var, validate='key').grid(row=5, column=1)

        Label(new_window, text="Top N Models:").grid(row=6, column=0)
        self.top_n_var = StringVar(new_window, value='50')
        Entry(new_window, textvariable=self.top_n_var, validate='key', validatecommand=vcmd).grid(row=6, column=1)

        Label(new_window, text="Threshold:").grid(row=7, column=0)
        self.threshold_var = StringVar(new_window, value='0.80')
        Entry(new_window, textvariable=self.threshold_var, validate='key', validatecommand=vcmd).grid(row=7, column=1)


        # Additional Options and Run Button
        def run_model():
            # Here you can collect all inputs and run the model based on selected options
            model_type = model_type_var.get()
            features_csv = feature_csv_entry.get()
            target_csv = target_csv_entry.get() if target_csv_entry.get() else None
            min_features_num = int(self.min_features_var.get())
            max_features_num = int(self.max_features_var.get()) if self.max_features_var.get() != 'None' else None
            top_n = int(self.top_n_var.get())
            threshold = float(self.threshold_var.get())
            leave_out_indices = self.leave_out_mols.get() if self.leave_out_mols.get() !='None' else None
            if leave_out_indices:
                leave_out_indices = self.leave_out_indices #convert_to_list_or_nested_list(leave_out_indices)
                
            
            if model_type == 'classification':
                classification = ClassificationModel({'features_csv_filepath': features_csv, 'target_csv_filepath': target_csv}, process_method='one csv', output_name='class', leave_out=leave_out_indices, min_features_num=min_features_num, max_features_num=max_features_num, metrics=None, return_coefficients=False,app=self)
                classification_results = classification.search_models(top_n=top_n, accuracy_threshold=threshold)
            elif model_type =='linear_regression':
               
                linear_regression = LinearRegressionModel({'features_csv_filepath': features_csv, 'target_csv_filepath': target_csv}, process_method='one csv', output_name='output', leave_out=leave_out_indices, min_features_num=min_features_num, max_features_num=max_features_num, metrics=None, return_coefficients=False,app=self)
                regression_results = linear_regression.search_models(top_n=top_n, initial_r2_threshold=threshold)

        run_button = Button(new_window, text="Run Model", command=run_model)
        run_button.grid(row=8, column=1, padx=10, pady=20)

    def handle_file_action(self, *args):
        selected_action = self.file_handler_var.get()
        if selected_action == "Smiles to XYZ":
            self.smiles_to_xyz_files()
        elif selected_action == "Create com Files":
            self.open_com_window()
        elif selected_action == "Log to Feather":
            self.log_to_feather()

    def log_to_feather(self):
        directory = filedialog.askdirectory()
        string_report=logs_to_feather(directory)
        self.show_result(f"Log to Feather Report: {string_report}")

    def renumber_directory(self):
        # Ask user if they want to create a new directory
        create_new_dir = messagebox.askyesno("Choose Directory", "Do you want to create a new directory for XYZ files?")
        
        if create_new_dir:
            # Let the user choose a location and name for the new directory
            new_dir_path = filedialog.asksaveasfilename(title="Select location for new directory",
                                                        filetypes=[('All Files', '*.*')])
            if new_dir_path:
                os.makedirs(new_dir_path, exist_ok=True)
                directory = new_dir_path
                os.chdir(directory)
                try:
                    [mol.write_xyz_file() for mol in self.molecules.molecules]
                except AttributeError:
                    self.show_result(f"Failed to write XYZ files to {directory} ...")
                
            else:
                return  # User cancelled the action
        else:
            # Let the user select an existing directory
            directory = filedialog.askdirectory()
            os.chdir(directory)
            if not directory:
                return  # User cancelled the action

        
        string_report,_ = batch_renumbering(directory)
        self.show_result(f"Renumbering Report: {string_report}")

    
    def smiles_to_xyz_files(self):
        # Initialize a Module1Handler object
        file_path = filedialog.askopenfilename(defaultextension=".csv",
                                       filetypes=[("Excel files", "*.csv"),
                                                  ("All files", "*.*")])

        module_handler = Module1Handler(file_path)
        os.chdir(module_handler.working_dir)
        help_functions.smiles_to_xyz_files(module_handler.smiles_list, module_handler.names_list, new_dir=True)
        
        
    def save_text(self):
        # dir_path = filedialog.askdirectory()
        text_name = filedialog.asksaveasfilename(defaultextension=".txt",
                                                filetypes=[("Text files", "*.txt"),
                                                            ("All files", "*.*")])
        dir_path = text_name.replace(text_name.split('/')[-1], '')
        if dir_path:
            os.chdir(dir_path)
            self.show_result(f" text saved at {dir_path}")
            with open(text_name, 'w') as f:
                f.write(self.output_text.get(1.0, "end-1c")) 
                f.close()
                                                   
    def morfeus_visualize(self):
        ## open a window to enter indices
        def get_indices():
            indices=askstring("Input", "Enter the indices for the atoms to visualize.")
            indices = convert_to_list_or_nested_list(indices)
            return indices
        Indices = get_indices()
        
        if self.molecules is not None:
            self.molecules.visualize_smallest_molecule_morfeus(Indices)

        
    def build_question_interface(self, parent, questions, loaded_entries=None):
        """
        Build the question interface in the parent window/frame.
        :param parent: Tkinter parent widget
        :param questions: List of question strings
        :param loaded_entries: Dict of question->answer (optional, for reloading)
        :return: Dict of question->Entry widgets
        """
        loaded_entries = loaded_entries or {}
        entry_widgets = {}

        # 1. Canvas + Scrollbar layout
        canvas = Canvas(parent, borderwidth=0, highlightthickness=0, width=800, height=600)
        scrollbar = Scrollbar(parent, orient='vertical', command=canvas.yview)
        canvas.configure(yscrollcommand=scrollbar.set)
        canvas.pack(side='left', fill='both', expand=True)
        scrollbar.pack(side='right', fill='y')

        container_frame = Frame(canvas)
        canvas.create_window((0, 0), window=container_frame, anchor="nw")

        def _on_configure(event):
            canvas.configure(scrollregion=canvas.bbox("all"))
        container_frame.bind("<Configure>", _on_configure)

        canvas.focus_set()
        self._bind_mouse_wheel(canvas)

        # 2. Top control buttons
        Button(container_frame, text="Visualize Basic Structure", command=self.visualize_smallest_molecule).pack(padx=10, pady=5)
        Button(container_frame, text="Choose Parameters", command=self.open_parameter_window).pack(pady=5)
        self.chosen_parameters = Label(container_frame, text=f"Chosen Parameters: {self.parameters}")
        self.chosen_parameters.pack(pady=5)

        # 3. Question widgets
        for i, question in enumerate(questions):
            if any(re.search(r, question, re.IGNORECASE) for r in ["threshold", "geometric center"]):
                continue

            frame_q = Frame(container_frame)
            frame_q.pack(pady=5, fill='x')
            label = Label(frame_q, text=question, wraplength=400)
            label.pack(side="left", padx=5)
            entry = Entry(frame_q, width=30)
            entry.pack(side="left", padx=5)

            # Map question types to handler functions
            self._handle_special_question_types(
                question, frame_q, entry_widgets, loaded_entries, entry
            )

            # Default loading of entry value
            if question in loaded_entries and not question.startswith("NPA manipulation"):
                entry.delete(0, 'end')
                entry.insert(0, loaded_entries[question])

            entry_widgets[question] = entry

        # 4. Bottom control buttons
        bottom_frame = Frame(container_frame)
        bottom_frame.pack(pady=10)
        Button(bottom_frame, text="Submit", command=lambda: self.submit_answers(entry_widgets, self.parameters)).pack(side='left', padx=5)
        Button(bottom_frame, text="Save input", command=lambda: self.save_input_json(entry_widgets)).pack(side='left', padx=5)
        Button(bottom_frame, text="Save output", command=lambda: self.submit_answers(entry_widgets, self.parameters, save_as=True)).pack(side='left', padx=5)
        Button(container_frame, text="Load input", command=lambda: self.on_load_answers(questions)).pack(pady=5)

        return entry_widgets

    def _visualize_dipole(self):
            ## open a window to enter indices
        def get_indices():
            indices=askstring("Input", "Enter the indices for the dipole atoms transformation.")
            indices = convert_to_list_or_nested_list(indices)
            return indices
        Indices = get_indices()
        
        if self.molecules is not None:
            # print(self.molecules,len(self.molecules.molecules))
            self.molecules.molecules[0].get_dipole_gaussian_df(Indices,visualize_bool=True)


    # Utility: Bind mouse wheel scrolling
    def _bind_mouse_wheel(self, canvas):
        self.master.bind_all("<MouseWheel>", lambda e: canvas.yview_scroll(int(-1 * (e.delta / 120)), "units"))
        self.master.bind_all("<Button-4>", lambda e: canvas.yview_scroll(-1, "units"))  # Linux scroll up
        self.master.bind_all("<Button-5>", lambda e: canvas.yview_scroll(1, "units"))   # Linux scroll down

    # Utility: Special question handlers
    def _handle_special_question_types(self, question, frame_q, entry_widgets, loaded_entries, entry):
        """Detect question type by prefix and add extra widgets as needed."""
        if question.startswith("Dipole atoms"):
            self._add_center_atoms_entry(frame_q, entry_widgets, loaded_entries)
            Button(frame_q, text="Show", command=self._visualize_dipole).pack(side="left", padx=5)
        elif question.startswith("NPA manipulation"):
            self._add_sub_atoms_entry(frame_q, entry_widgets, loaded_entries)
            self._load_entry(entry, loaded_entries, question)
        elif question.startswith("Ring Vibration atoms"):
            Button(frame_q, text="Show", command=lambda: self.open_image(r"pictures\rings.png")).pack(side="left", padx=5)
        elif question.startswith("Stretching Vibration atoms"):
            self._add_threshold_entry("Stretch Threshold", frame_q, entry_widgets, loaded_entries, default=1600)
        elif question.startswith("Bending Vibration atoms"):
            self._add_threshold_entry("Bend Threshold", frame_q, entry_widgets, loaded_entries, default=1600)
        elif question.startswith("Sterimol atoms"):
            Button(frame_q, text="Show", command=self.morfeus_visualize).pack(side="left", padx=5)
            # add a button to drop 
    
    def _drop_atoms_sterimol(self, frame, entry_widgets, loaded_entries):
        Label(frame, text="Drop Atoms:").pack(side="left", padx=5)
        e = Entry(frame, width=10)
        e.pack(side="left", padx=5)
        entry_widgets["Drop_Atoms"] = e
        self._load_entry(e, loaded_entries, "Drop_Atoms")

    def _add_center_atoms_entry(self, frame, entry_widgets, loaded_entries):
        Label(frame, text="Center_Atoms:").pack(side="left", padx=5)
        e = Entry(frame, width=10)
        e.pack(side="left", padx=5)
        entry_widgets["Center_Atoms"] = e
        self._load_entry(e, loaded_entries, "Center_Atoms")

    def _add_sub_atoms_entry(self, frame, entry_widgets, loaded_entries):
        Label(frame, text="Sub-Atoms:").pack(side="left", padx=5)
        e = Entry(frame, width=10)
        e.pack(side="left", padx=5)
        entry_widgets["Sub-Atoms"] = e
        self._load_entry(e, loaded_entries, "Sub_Atoms")

    def _add_threshold_entry(self, key, frame, entry_widgets, loaded_entries, default):
        Label(frame, text="Threshold:").pack(side="left", padx=5)
        e = Entry(frame, width=10)
        e.pack(side="left", padx=5)
        entry_widgets[key] = e
        if key in loaded_entries:
            e.delete(0, 'end')
            e.insert(0, loaded_entries[key])
        else:
            e.insert(0, default)

    def _load_entry(self, entry, loaded_entries, key):
        if key in loaded_entries:
            entry.delete(0, 'end')
            entry.insert(0, loaded_entries[key])



    def on_load_answers(self, questions):
        """
        Loads answers from a file, closes the old question window,
        and rebuilds it with loaded answers.
        """
        loaded_entries = load_answers_json()
        if not loaded_entries:
            return  # If user cancels or no entries found, do nothing

        # Destroy the existing window
        if self.new_window:
            self.new_window.destroy()

        # Create a fresh Toplevel
        self.new_window = Toplevel(self.master)
        self.new_window.title("Questions")
        print(f'loading new window with input')
        # Build the interface with the loaded answers
        self.build_question_interface(self.new_window, questions, loaded_entries=loaded_entries)


    

    def save_input_json(self, entry_widgets):
        """
        Save the current entry widget values to a JSON file.
        The output file will be human-editable and easy to reload.
        """
        file_path = filedialog.asksaveasfilename(
            defaultextension=".json",
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")]
        )
        if not file_path:
            return  # User cancelled

        data = {}
        for question, entry in entry_widgets.items():
            value = entry.get()
            # Try to convert to int or float if possible
            try:
                value = int(value)
            except ValueError:
                try:
                    value = float(value)
                except ValueError:
                    pass  # Keep as string if not a number
            data[question] = value

        with open(file_path, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=4, ensure_ascii=False)

        print(f"Input saved to {file_path}")

    def save_input(self, entry_widgets, save_as=False):
        """
        Saves user input (entry_widgets) to a text file if 'save_as' is True.
        """
        answers = {}
        for question, entry in entry_widgets.items():
            try:
                value = entry.get() if hasattr(entry, 'get') else entry
                
                # Check if the entry is an integer greater than 100
                if value.isdigit() and int(value) > 100:
                    answers[question] = int(value)
                else:
                    # Convert the value using convert_to_list_or_nested_list
                    answers[question] = convert_to_list_or_nested_list(value)
                    
            except AttributeError:
                # Handle direct string or int entries
                if isinstance(entry, int) and entry > 100:
                    answers[question] = entry
                else:
                    answers[question] = convert_to_list_or_nested_list(entry)

        if save_as:
            file_path = filedialog.asksaveasfilename(
                defaultextension=".txt",
                filetypes=[("Text files", "*.txt"), ("All files", "*.*")]
            )
            if file_path:
                with open(file_path, 'w') as f:
                    for question, answer in answers.items():
                        f.write(f"{question}\n{answer}\n\n")

    def submit_answers(self, entry_widgets, parameters, save_as=False):
        """
        Gathers user input from entry_widgets, applies parameters,
        extracts features, and optionally saves results to a file.
        """
        answers = {}
        
        for question, entry in entry_widgets.items():
            key = question.split()[0].lower()  # Use first word as key, lowercase
            try:
                answers[key] = entry.get()
            except AttributeError:
                answers[key] = entry
        

        radii = self.parameters['Radii'] 
        iso= self.parameters['Isotropic'] if 'Isotropic' in self.parameters else False
        
        # Example: calling a method from self.molecules (stub for demonstration)
        comp_set = self.molecules.get_molecules_features_set(answers, parameters=self.parameters)
        self.show_result(f"Extracted Features: {comp_set}")

        # For demonstration, just print
        print("Submitted answers:", answers, "\n\n")
        print("Using parameters:", parameters)

        # If save_as, save both text and CSV
        if save_as:
           
            file_path_csv = filedialog.asksaveasfilename(
                defaultextension=".csv",
                filetypes=[("csv files", "*.csv"), ("All files", "*.*")]
            )
            if file_path_csv:
                comp_set.to_csv(file_path_csv, index=True)

    def open_parameter_window(self):
        """
        Opens a small window to adjust self.parameters.
        """
        window = Toplevel(self.master)
        window.title("Parameters")
        window.grab_set()

        frame = Frame(window)
        frame.pack(pady=5)

        # Dipole Mode
      
        # Radii
        var_radii = StringVar(frame)
        var_radii.set(self.parameters.get('Radii', 'bondi'))
        var_radii.trace_add("write", lambda *args: apply_parameters())

        # Isotropic
        var_iso = StringVar(frame)
        var_iso.set(str(self.parameters.get('Isotropic', False)))
        var_iso.trace_add("write", lambda *args: apply_parameters())

        # Menu Widgets
     
        OptionMenu(frame, var_radii, 'bondi', 'CPK', 'Pyykko').grid(row=0, column=1, padx=5)
        OptionMenu(frame, var_iso, 'True', 'False').grid(row=0, column=2, padx=5)

        def apply_parameters():
            
        
            self.parameters['Radii'] = var_radii.get()
            self.parameters['Isotropic'] = (var_iso.get() == 'True')
            if hasattr(self, 'chosen_parameters'):
                self.chosen_parameters.config(text=f"Chosen Parameters: {self.parameters}")

            print('parameters',self.parameters)
            print('chosen',self.chosen_parameters)
            

        Button(frame, text="Apply", command=window.destroy).grid(row=0, column=3, padx=5)

    def open_question_window(self):
        if not hasattr(self, 'parameters'):
            self.parameters = {'Radii': 'bondi', 'Isotropic': False}

        questions = [
                "Ring Vibration atoms - by order -> Pick primary atom and para to it: \n example: 13,17",
                "Stretching Vibration atoms- enter bonded atom pairs: \n example: 1,2 4,5",
                "Strech Threshold - default is 1600 for carbonyl stretch",
                "Bending Vibration atoms - enter atom pairs that have a common atom: \n example: 4,7",
                "Bend Threshold - default is 1600 ",
                "Center Atoms Dipole - indices to move Geometric Center: \n example 1,2,3,4,5,6 - move to Ring Center",
                "Dipole atoms - indices for coordination transformation: \n example: 4,5,6 - origin, y-axis, new xy plane",
                "Sub-Atoms NPA - Insert atoms to show NPA: \n example: 1,2,3,4,5,6",
                "NPA atoms - Insert atoms to show NPA: \n example: 1,2,4",
                "charges values - Insert atoms to show charge: \n example: 1,2,3,4",
                "charge_difference - Insert atoms to show charge difference: \n example: 1,2 3,4",
                "Sterimol atoms - Primary axis along: \n example: 7,8",
                "Bond length - Atom pairs to calculate difference: \n example: 1,2 4,5",
                "Bond Angle - Insert a list of atom triads/quartets for which you wish to have angles/dihedrals: \n example: 1,3,4 5,6,7,4"
            ]
        
        # Create a new Toplevel
        self.new_window = Toplevel(self.master)
        self.new_window.title("Questions")

        # Build the interface with no previously loaded answers
        self.build_question_interface(self.new_window, questions, loaded_entries=None)
                
        

    def open_com_window(self):
        options_window = Toplevel(self.master)
        options_window.title("Conversion Options")
        options_window.grab_set()  # Make the window moda
        gaussian_options = {
            'functionals': [
                'HF','b97d3', 'B3LYP', 'PBE', 'M062x', 'CAM-B3LYP', 'MP2', 'CCSD'
            ],
            'basis_sets': [
                'def2tzvp', 'STO-3G', '3-21G', '6-31G', '6-31G(d) int=sg1', '6-31G(d,p)', '6-31G(2df,p)',
                '6-31+G(d,p)', '6-311G(d,p)', '6-311+G(d,p)', '6-311++G(d,p)', '6-311++G(2d,p)', '6-311++G(3df,2p)'
            ],
            'tasks': [
                'sp', 'opt'] # , 'freq', 'nmr', 'ts', 'irc', 'opt+freq', 'opt+freq+nmr'
        }


        # Create OptionMenus for functional, basis_set, and task
        self.functional_var = StringVar(value='B3LYP')
        Label(options_window, text='Functional:').pack()
        OptionMenu(options_window, self.functional_var, *gaussian_options['functionals']).pack()

        self.basisset_var = StringVar(value='6-31G(d)')
        Label(options_window, text='Basis Set:').pack()
        OptionMenu(options_window, self.basisset_var, *gaussian_options['basis_sets']).pack()

        self.task_var = StringVar(value='opt')
        Label(options_window, text='Task:').pack()
        OptionMenu(options_window, self.task_var, *gaussian_options['tasks']).pack()

        # Parameters to be entered by the user

        self.charge_var = StringVar(value='0 1')
        self.charge_var = StringVar(value='n')
        self.title_var = StringVar(value='title')
        # Create labels and entry widgets for each parameter
        Label(options_window, text='Spin & Charge:').pack()
        Entry(options_window, textvariable=self.charge_var).pack()
        Label(options_window, text='Title:').pack()
        Entry(options_window, textvariable=self.title_var).pack()

        self.freeze_var = StringVar(value='')
        Label(options_window, text='Freeze:').pack()
        Entry(options_window, textvariable=self.freeze_var).pack()


        # Add button to select directory and execute conversion
        customtkinter.CTkButton(options_window, text="Select Directory and Convert", command=self.convert_xyz_to_com).pack()
        customtkinter.CTkButton(options_window, text="Create New Directory", command=self.create_new_directory).pack()

    def convert_xyz_to_com(self):
        folder_selected = filedialog.askdirectory(title="Select folder with XYZ files")
        if not folder_selected:
            return

        output_folder = filedialog.askdirectory(title="Select output folder for .com files")
        if not output_folder:
            return

        os.chdir(folder_selected)
        xyz_files = [f for f in os.listdir(folder_selected) if f.endswith('.xyz')]
        if not xyz_files:
            self.show_result("No .xyz files found in the selected directory.")
            return

        n_success, n_fail = 0, 0

        for filename in xyz_files:
            try:
                # Fetch user-selected options from UI variables
                functional = self.functional_var.get()
                basis = self.basisset_var.get()
                charge = self.charge_var.get()
                title = self.title_var.get()
                task = self.task_var.get()
                freeze = self.freeze_var.get() if hasattr(self, "freeze_var") else None
         

                # Your actual file writing function (adapt as needed)
                file_handlers.write_gaussian_file(
                    filename, functional, basis, charge, title, task, freeze
                )

                com_filename = filename.replace('.xyz', '.com')
                output_path = os.path.join(output_folder, com_filename)

                shutil.move(com_filename, output_path)
                self.show_result(f"Converted and moved: {com_filename}")
                n_success += 1
            except Exception as e:
                self.show_result(f"Failed for {filename}: {e}")
                n_fail += 1

        self.show_result(f"Done! {n_success} files converted, {n_fail} failed.")
        messagebox.showinfo("Conversion complete", f"{n_success} files converted successfully.\n{n_fail} files failed.")


            


    def get_answers(self):
        for question, entry in self.answers.items():
            self.show_result(f"{question}: {entry.get()}")

    def filter_molecules(self):
        self.new_window = Toplevel(self.master)
        self.new_window.title("Filter Molecules")

        canvas = Canvas(self.new_window)
        scrollbar = Scrollbar(self.new_window, orient='vertical', command=canvas.yview)
        scrollbar.pack(side='right', fill='y')
        scrollbar.bind("<MouseWheel>", lambda event: canvas.yview_scroll(int(-1 * (event.delta / 120)), "units"))
        canvas.pack(side='left', fill='both', expand=True)
        canvas.configure(yscrollcommand=scrollbar.set)

        frame = Frame(canvas)
        canvas_frame = canvas.create_window((0, 0), window=frame, anchor='nw')

        self.check_vars = [IntVar(value=1) for _ in self.molecules.old_molecules_names]
        for index, molecule in enumerate(self.molecules.old_molecules_names):
            Checkbutton(frame, text=molecule, variable=self.check_vars[index]).pack(anchor='w')

        Button(frame, text="Submit", command=self.get_selected_molecules).pack()
        Button(frame, text="Uncheck", command=self.uncheck_all_boxes).pack()
        Button(frame, text="Check", command=self.check_all_boxes).pack()

        frame.update_idletasks()
        canvas.config(scrollregion=canvas.bbox('all'))
        # allow scrooling with scrollwheel
        canvas.bind_all("<MouseWheel>", lambda event: canvas.yview_scroll(int(-1 * (event.delta / 120)), "units"))

    def check_all_boxes(self):
        for var in self.check_vars:
            var.set(1)

    def uncheck_all_boxes(self):
        for var in self.check_vars:
            var.set(0)

    def get_selected_molecules(self):
        self.molecules.molecules_names = self.molecules.old_molecules_names
        self.molecules.molecules = self.molecules.old_molecules
        selected_indices = [i for i, var in enumerate(self.check_vars) if var.get() == 1]
        self.show_result(f"Selected indices: {selected_indices}")
        self.new_window.destroy()
        self.molecules.filter_molecules(selected_indices)
        self.show_result(f"Initializing Molecules: {self.molecules.molecules_names}")

    def choose_directory(self):
        folder_selected = filedialog.askdirectory()
        if folder_selected:
            os.chdir(folder_selected)
            self.show_result(f"Working directory changed to {folder_selected}")


    def create_new_directory(self):
        folder_selected = filedialog.askdirectory()
        os.chdir(folder_selected)
        folder_name = filedialog.asksaveasfilename(title="Enter a Name")
        self.new_directory_path = os.path.join(folder_selected, folder_name)
        if folder_name:
            try:
                os.makedirs(folder_name)
                self.show_result(f"Directory {folder_name} created.")
            except FileExistsError:
                self.show_result(f"Directory {folder_name} already exists.")
            except Exception as e:
                self.show_result(f"An error occurred: {e}")
        os.chdir(folder_name)
        self.show_result(f"Working directory changed to {folder_name}")

    def browse_directory(self):
        
    
        folder_selected = filedialog.askdirectory(initialdir=self.current_directory)

        if folder_selected:
            self.folder_path.set(folder_selected)
            self.initialize_molecules()

    def initialize_molecules(self):
        add_mols=False
        directory = self.folder_path.get()
        os.chdir(directory)
        files_list = os.listdir(directory)
        feather_files = [file for file in files_list if file.endswith('.feather')]
        if directory:
            if len(feather_files) == 0:
                dir_list=os.listdir()
                try:
                    os.mkdir('feather_files')
                except FileExistsError:
                    pass
                for dir in dir_list:
                    os.chdir(dir)     
                    feather_file = [file for file in os.listdir() if (file.endswith('.feather') and file.split('-')[0]=='xyz_files')][0]
                    try:
                        shutil.copy(feather_file, directory + '/feather_files')
                    except shutil.SameFileError:    
                        pass
                    os.chdir('..')
                if hasattr(self, 'molecules') and getattr(self.molecules, 'molecules', None):
                    previous_molecules = self.molecules.molecules
                    previous_molecules_names = self.molecules.molecules_names
                    add_mols = messagebox.askyesno("Add", "Add to existing Molecules set ?", parent=self.master)

                self.molecules = Molecules(directory+'/feather_files') # , renumber=True
                if add_mols:
                    print(f"previous_molecules: {previous_molecules}")  # Debugging
                    self.molecules.molecules.extend(previous_molecules)
                    self.molecules.molecules_names.extend(previous_molecules_names)
                self.show_result(f"Molecules initialized with directory: {self.molecules.molecules_names}\n")
                self.show_result(f'Failed to load Molecules: {self.molecules.failed_molecules}\n')
                self.show_result(f"Initializing Molecules with directory: {directory}\n")  # Debugging
                os.chdir('..')
            else:

                if hasattr(self, 'molecules') and getattr(self.molecules, 'molecules', None):
                    previous_molecules = self.molecules.molecules
                    previous_molecules_names = self.molecules.molecules_names
                    add_mols = messagebox.askyesno("Add", "Add to existing Molecules set ?", parent=self.master)

                print(f"Initializing Molecules with directory: {directory}")  # Debugging
                self.show_result(f"Initializing Molecules with directory: {directory}\n")  
                self.molecules = Molecules(directory) # , renumber=True
                if add_mols:
                    print(f"previous_molecules: {previous_molecules}")
                    self.molecules.molecules.extend(previous_molecules)
                    self.molecules.molecules_names.extend(previous_molecules_names)

                self.show_result(f"Molecules initialized : {self.molecules.molecules_names}\n")
                self.show_result(f'Failed to load Molecules: {self.molecules.failed_molecules}\n')
                self.show_result(f"Initializing Molecules with directory: {directory}\n") 
                
            

    def open_param_window(self):
        
        selected_method = self.method_var.get()
        description_text = f"Enter parameters for {selected_method}:"
        self.param_description.config(text=description_text)

        if selected_method == "Windows Command":
            self.show_result(f"Use as Command Line:\n")
        elif selected_method == "get_sterimol_dict":
            self.show_result(f"Method: {(self.molecules.molecules[0].get_sterimol.__doc__)}\n)")
        elif selected_method == "get_npa_dict":
            self.show_result(f"Method: {(self.molecules.molecules[0].get_npa_df.__doc__)}\n)")
        elif selected_method == "get_stretch_dict":
            self.show_result(f"Method: {self.molecules.molecules[0].get_stretch_vibration.__doc__}\n)")
        elif selected_method == "get_ring_dict":
            self.show_result(f"Method: {self.molecules.molecules[0].get_ring_vibrations.__doc__}\n)")
        elif selected_method == "get_dipole_dict":
            self.show_result(f"Method: {self.molecules.molecules[0].get_dipole_gaussian_df.__doc__}\n)")
        elif selected_method == "get_bond_angle_dict":
            self.show_result(f"Method: {self.molecules.molecules[0].get_bond_angle.__doc__}\n)")
        elif selected_method == "get_bond_length_dict":
            self.show_result(f"Method: {self.molecules.molecules[0].get_bond_length.__doc__}\n)")
        elif selected_method == "get_charge_dict":
            self.show_result(f"Method: {self.molecules.molecules[0].get_charge_df.__doc__}\n)")
        elif selected_method == "get_charge_diff_dict":
            self.show_result(f"Method: {self.molecules.molecules[0].get_charge_diff_df.__doc__}\n)")
        elif selected_method == "get_bending_dict":
            self.show_result(f"Method: {self.molecules.molecules[0].get_bend_vibration.__doc__}\n)")
        
    def show_result(self, result):
        # Update Text widget instead of creating a new Toplevel window
        self.output_text.insert('end', str(result) + '\n')
        self.output_text.see('end')  # Auto-scroll to the end

    def activate_method(self):

        method = self.method_var.get()
        params = self.param_entry.get()

        # Now you have the method and parameters, you can activate the method
        if method == "Windows Command":
            self.use_command_line(params)
            
        elif method == "get_sterimol_dict":
            self.get_sterimol(params)
        elif method == "get_npa_dict":
            self.get_npa(params)
        elif method == "get_stretch_dict":
            self.get_stretch(params)
        elif method == "get_ring_dict":
            self.get_ring(params)
        elif method == "get_dipole_dict":
            self.get_dipole(params)
        elif method == "get_bond_angle_dict":
            self.get_bond_angle(params) 
        elif method == "get_bond_length_dict":
            self.get_bond_length(params)
        elif method == "get_charge_dict":
            self.get_charge(params)
        elif method == "get_charge_diff_dict":
            self.get_charge_diff(params)
        elif method == "get_bending_dict":
            self.get_bending(params)
        elif method == "get_molecules_comp_set":
            self.get_molecules_features_set()

    def use_command_line(self, params):
        try:
            # Execute the command and capture the output
            result = subprocess.run(params, shell=True, check=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            
            # Print the standard output of the command
            self.show_result(f"Output: {result.stdout}\n")
            

            # Optionally, print the standard error if there is any
            if result.stderr:
                self.show_result(f"Errors: {result.stderr}\n")
                
        except subprocess.CalledProcessError as e:
            # This block will run if the command exits with a non-zero status
            self.show_result(f"An error occurred: {e}")

    def get_sterimol(self,base_atoms_str):
        base_atoms = convert_to_list_or_nested_list(base_atoms_str)
        sterimol_data = self.molecules.get_sterimol_dict(base_atoms)
        self.show_result(f"Sterimol values:\n {sterimol_data}\n")

    def get_npa(self,base_atoms_str):
        if base_atoms_str:
            base_atoms = convert_to_list_or_nested_list(base_atoms_str)
            npa_data = self.molecules.get_npa_dict(base_atoms)
            self.show_result(f"NPA Charges:\n {npa_data}\n")

    def get_stretch(self,base_atoms_str):
        if base_atoms_str:
            base_atoms = convert_to_list_or_nested_list(base_atoms_str)
            stretch_data = self.molecules.get_stretch_vibration_dict(base_atoms)
            self.show_result(f"Stretch Vibration:\n {stretch_data}\n")

    def get_ring(self,base_atoms_str):
        if base_atoms_str:
            base_atoms = convert_to_list_or_nested_list(base_atoms_str)
            ring_data = self.molecules.get_ring_vibration_dict(base_atoms)
            self.show_result(f"Ring Vibrations:\n {ring_data}\n")

    def get_dipole(self,base_atoms_str):

        if base_atoms_str:
            base_atoms = convert_to_custom_nested_list(base_atoms_str)
            dipole_data = self.molecules.get_dipole_dict(base_atoms)
            self.show_result(f"Dipole Moment:\n {dipole_data}\n")
    
    def get_bond_angle(self,base_atoms_str):
        if base_atoms_str:
            base_atoms = convert_to_list_or_nested_list(base_atoms_str)
            bond_angle_data = self.molecules.get_bond_angle_dict(base_atoms)
            self.show_result(f"Bond Angles:\n {bond_angle_data}\n")

    def get_bond_length(self,base_atoms_str):
        if base_atoms_str:
            base_atoms = convert_to_list_or_nested_list(base_atoms_str)
            bond_length_data = self.molecules.get_bond_length_dict(base_atoms)
            self.show_result(f"Bond Lengths:\n {bond_length_data}\n")

    def get_charge(self,base_atoms_str):
        if base_atoms_str:
            single_numbers = [int(num) for num in base_atoms_str.split(',')]
            charge_data = self.molecules.get_charge_df_dict(single_numbers)
            self.show_result(f"Charge Analysis:\n {charge_data}\n")
    
    def get_charge_diff(self,base_atoms_str):
        if base_atoms_str:
            base_atoms = convert_to_list_or_nested_list(base_atoms_str)
            charge_diff_data = self.molecules.get_charge_diff_dict(base_atoms)
            self.show_result(f"Charge Differences:\n {charge_diff_data}\n")

    def get_bending(self,base_atoms_str):
        if base_atoms_str:
            base_atoms = convert_to_list_or_nested_list(base_atoms_str)
            bending_data = self.molecules.get_bend_vibration_dict(base_atoms)
            self.show_result(f"Bending Vibrations:\n {bending_data}\n")
    # TODO enabe to choose which mols to visualize
            
    def filter_molecules_vis(self):
        self.new_window = Toplevel(self.master)
        self.new_window.title("Filter Molecules")

        canvas = Canvas(self.new_window)
        scrollbar = Scrollbar(self.new_window, orient='vertical', command=canvas.yview)
        scrollbar.pack(side='right', fill='y')
        scrollbar.bind("<MouseWheel>", lambda event: event.yview_scroll(int(-1 * (event.delta / 120)), "units"))
        canvas.pack(side='left', fill='both', expand=True)
        canvas.configure(yscrollcommand=scrollbar.set)
        frame = Frame(canvas)
        canvas_frame = canvas.create_window((0, 0), window=frame, anchor='nw')

        self.check_vars = [IntVar(value=0) for _ in self.molecules.old_molecules_names]
        for index, molecule in enumerate(self.molecules.old_molecules_names):
            Checkbutton(frame, text=molecule, variable=self.check_vars[index]).pack(anchor='w')

        Button(frame, text="Visualize", command=self.visualize).pack()
        Button(frame, text="Uncheck", command=self.uncheck_all_boxes).pack()
        Button(frame, text="Check", command=self.check_all_boxes).pack()

        frame.update_idletasks()
        canvas.config(scrollregion=canvas.bbox('all'))
        # allow scrooling with scrollwheel
        canvas.bind_all("<MouseWheel>", lambda event: event.widget.yview_scroll(int(-1*(event.delta/120)), "units"))



    def visualize(self):
        
        selected_indices = [i for i, var in enumerate(self.check_vars) if var.get() == 1]
        self.show_result(f"Selected indices: {selected_indices}")
        self.new_window.destroy()
        self.molecules.visualize_molecules(selected_indices)

   
    def visualize_smallest_molecule(self):
        if self.molecules:
            self.molecules.visualize_smallest_molecule()
            

    def open_image(self,image_path):
        new_window = Toplevel(self.master)
        new_window.title("Image Display")
        # Load and resize the image
        image = Image.open(image_path)
        # Resize the image to desired dimensions, e.g., (width, height)
        photo = ImageTk.PhotoImage(image)
        # Display the image
        label = Label(new_window, image=photo)
        label.image = photo  # Keep a reference!
        label.grid(row=0, column=0)
        # Bind a mouse click event to the label
        def close_window(event=None):
            new_window.destroy()
        label.bind("<Button-1>", close_window)
        new_window.protocol("WM_DELETE_WINDOW", close_window)

    def export_data(self):
        self.molecules.extract_all_dfs()
        self.show_result(f"DataFrames extracted.")
    
    def export_xyz(self):
        self.molecules.export_all_xyz()
        self.show_result(f"XYZ files exported.")


def get_function_args(function):
    sig = inspect.signature(function)
    return sig.parameters

def prompt_for_args(method):
    sig = inspect.signature(method)
    values = []
    for param in sig.parameters.values():
        param_name = param.name.replace("_", " ").capitalize()
        if param.default != inspect.Parameter.empty:
            
            prompt = f"Enter value for {param_name} (default: {param.default}): "
        else:
            prompt = f"Enter value for {param_name}: "

        user_input = input(prompt) #.strip()
        
        if user_input == "":
            # Use default value if input is empty
            values.append(param.default)
        else:
            try:
                # values.append(ast.literal_eval(user_input))
                values.append(user_input)
            except (ValueError, SyntaxError):

                values.append(user_input)
 
    return values



def resolve_n_jobs(user_n_jobs) -> int:
    """Pick n_jobs from CLI, NSLOTS, or all logical cores."""
    if user_n_jobs is not None:
        return int(user_n_jobs)
    env_slots = os.environ.get("NSLOTS")
    if env_slots:
        try:
            return int(env_slots)
        except ValueError:
            pass
    return int(psutil.cpu_count(logical=True))

def add_modeling_path():
    """Add the modeling path to sys.path if given (or if a sensible default exists)."""
    default_path = os.path.join(os.path.dirname(__file__), "modeling")
    path = default_path
    if os.path.isdir(path) and path not in sys.path:
        sys.path.append(path)

def safe_chdir(workdir, features_csv):
    """
    Set working directory:
    - If --workdir given, chdir there.
    - Else chdir to the directory containing features_csv.
    """
    if workdir and os.path.isdir(workdir):
        os.chdir(workdir)
        return
    if features_csv:
        abs_feat = os.path.abspath(features_csv)
        base = os.path.dirname(abs_feat)
        if os.path.isdir(base):
            os.chdir(base)

def run_experiment(
    mode: str,
    features_csv: str,
    target_csv: str = "",
    y_value: str = "output",
    n_jobs: int = 1,
    min_features_num: int = 2,
    max_features_num: int = 4,
    threshold: float = 0.70,
    top_n: int = 50,
    leave_out: Optional[List[str]] = None,
    process_method: str = "one csv",
    bool_parallel: bool = False,
) -> float:
    """
    Run regression or classification and return elapsed seconds.
    - For regression, calls: fit_and_evaluate_combinations(...) if available, else search_models(...)
    - For classification, calls: search_models(...)
    """


    csv_filepaths = {
        "features_csv_filepath": features_csv,
        "target_csv_filepath": target_csv or ""
    }

    leave_out = list(leave_out or [])

    # Instantiate correct model class
    if mode.lower() == "regression":
        model = LinearRegressionModel(
            csv_filepaths,
            process_method=process_method,
            y_value=y_value,
            leave_out=leave_out,
            min_features_num=min_features_num,
            max_features_num=max_features_num,
           
        )
    elif mode.lower() == "classification":
        # Assumes you have a ClassificationModel class with similar init signature
        model = ClassificationModel(
            csv_filepaths,
            process_method=process_method,
            y_value=y_value,
            leave_out=leave_out,
            min_features_num=min_features_num,
            max_features_num=max_features_num,
    
        )
    else:
        raise ValueError("mode must be 'regression' or 'classification'")

    # Decide parallel cores
    eff_n_jobs = n_jobs if bool_parallel else 1
    print(f"\n[run_experiment] mode={mode} | n_jobs={eff_n_jobs} (requested={n_jobs}, parallel={bool_parallel})")
    print(f"[run_experiment] features={features_csv} | target={target_csv or '(none)'} | y_value={y_value}")
    print(f"[run_experiment] min/max features: {min_features_num}/{max_features_num}\n")

    # Run
    start = time.time()
    if mode.lower() == "regression":
        if hasattr(model, "search_models"):
            model.search_models(
                top_n=top_n,
                n_jobs=eff_n_jobs,
                threshold=threshold,
                bool_parallel=bool_parallel
            )
        else:
            raise AttributeError("Regression model lacks search_models.")
    else:  # classification
        if hasattr(model, "search_models"):
            model.search_models(
                top_n=top_n,
                n_jobs=eff_n_jobs,
                threshold=threshold,
                bool_parallel=bool_parallel
            )
        else:
            raise AttributeError("Classification model lacks search_models().")

    elapsed = time.time() - start
    print(f"\n[run_experiment] Completed in {elapsed:.2f} s\n")
    gc.collect()
    return elapsed

def run_gui_app():
    print("Running GUI app...")
    root = Tk()
    app = MoleculeApp(root)
    root.mainloop()
    # Your code to launch the GUI app goes here

def run_feature_extraction(input_file, output_file, verbose=False):
    print(f"Running feature extraction...")
    print(f"Input file: {input_file}")
    print(f"Output file: {output_file}")
    if verbose:
        print(f"Verbose mode is ON")
    

def load_molecules(molecules_dir_name, renumber=False):
    return Molecules(molecules_dir_name, renumber=renumber)

def interactive_cli(molecules):
    exclude_methods = ['get_molecules_features_set', 'visualize_smallest_molecule_morfeus','visualize_smallest_molecule', 'filter_molecules', 'renumber_molecules']
    methods = [method for method in dir(molecules) if callable(getattr(molecules, method)) and not method.startswith("__") and method not in exclude_methods]
    
    while True:
        print("\nAvailable methods:")
        for method in methods:
            print(f"- {method}")
        command = input("\nEnter method to invoke (or 'help method_name' to see documentation, or 'exit' to quit): ")
        if command.lower() == 'exit':
            break
        
        if command.startswith('help '):
            method_name = command.split(' ')[1]
            if hasattr(molecules, method_name):
                method = getattr(molecules, method_name)
                help(method)
            else:
                print(f"Unknown method: {method_name}")
            continue
        
        if not hasattr(molecules, command):
            print(f"Unknown method: {command}")
            continue
        method = getattr(molecules, command)
        user_values = prompt_for_args(method)
        try:
            user_values[0]=convert_to_list_or_nested_list(user_values[0])
            result = method(*user_values)
            # If the method returns a dictionary of DataFrames, print them
            if isinstance(result, dict):
                for name, df in result.items():
                    print(f"\nResults for {name}:\n{df}\n")
            elif isinstance(result, pd.DataFrame):
                print(f"\n{result}\n")
            elif result is not None:
                print(result)
        except Exception as e:
            print(f"Error invoking {command}: {e}")
            continue

def interactive_modeling(csv_path):
    model_bool=input("Do you want to model the data? (yes/no): ").strip().lower()


def main():
    parser = argparse.ArgumentParser(description="MolFeatures Package")
    subparsers = parser.add_subparsers(dest="command")

    # Subcommand for running the GUI app
    gui_parser = subparsers.add_parser("gui", help="Run the GUI app")
    interactive_parser = subparsers.add_parser("interactive", help="Start interactive CLI for cmd line operations")
    conver_parser = subparsers.add_parser("logs_to_feather", help="Convert log files to feather files")
    cube_parser = subparsers.add_parser("cube", help="Calculates cube sterimol from cube files")
    sterimol_parser = subparsers.add_parser("sterimol", help="Calculate sterimol values from xyz files")
    install_parser = subparsers.add_parser("install", help="Install required packages")
    model_parser = subparsers.add_parser("model", help="Run regression or classification")
    model_parser.add_argument("-m", "--mode", choices=["regression", "classification"], required=True,
                              help="Which task to run.")
    model_parser.add_argument("-f", "--features_csv", required=True, help="Path to features CSV.")
    model_parser.add_argument("-t", "--target_csv", default="", help="Path to target/labels CSV.")
    model_parser.add_argument("-y", "--y_value", default="output", help="Base name for output files.")
    model_parser.add_argument("-j", "--n_jobs", type=int,
                              help="Cores to use (-1 = all). If omitted, uses NSLOTS or all logical cores.")
    model_parser.add_argument("--min-features", type=int, default=4, help="Minimum features per model.")
    model_parser.add_argument("--max-features", type=int, default=4, help="Maximum features per model.")
    model_parser.add_argument("--top-n", type=int, default=50, help="How many top models to keep/evaluate.")
    model_parser.add_argument("--bool-parallel", action="store_false", help="Enable parallel evaluation.")
    model_parser.add_argument("--threshold", type=float, default=0.70, help="Initial threshold (regression(R2)/classification(mcfadden_R2)).")
    model_parser.add_argument("--leave-out", nargs="*", default=[], help="Space-separated list of samples to leave out.")
   
    args = parser.parse_args()

    feature_extraction = subparsers.add_parser("feature_extraction", help="Run feature extraction - complete set - from input file")
    feature_extraction.add_argument("-i", "--input", required=True, help="Path to input file.")
    feature_extraction.add_argument("-o", "--output", required=True, help="Path to output file.")
    feature_extraction.add_argument("--verbose", action="store_true", help="Enable verbose output.")

    if args.command == "gui":
        run_gui_app()
        
    elif args.command == "logs_to_feather":
        log_dir = input("Enter the path to the log files directory: ")
        logs_to_feather(log_dir)
        

    elif args.command == "cube":
        while True:
            cube_file_path = input("Enter the full path to directory with cube files: ")
            base_atoms = input("Enter the atom indices: ")
            base_atoms=convert_to_list_or_nested_list(base_atoms)
            sterimol=cube_many(cube_file_path, base_atoms)
            ## save to csv
            # df=pd.DataFrame(sterimol.sterimol_df)
            sterimol.sterimol_df.to_csv('cube_sterimol.csv')

            another_dir = input("Do you want to input another directory? (yes/no): ").strip().lower()
            if another_dir != 'yes':
                break
    elif args.command == "sterimol":
        while True:
            xyz_dir = input("Enter the path to the xyz files directory: ")
            base_atoms = input("Enter the atom indices: ")
            params = inspect.signature(Molecules_xyz.get_sterimol_df)
            radii=list(params.parameters.values())[2].default
            radius_input = input(f"Enter the radii (default: {radii}): ")
            if radius_input=='':
                radius_input=radii
            base_atoms=convert_to_list_or_nested_list(base_atoms)
            sterimol=Molecules_xyz(xyz_dir)
            sterimol_df=sterimol.get_sterimol_df(base_atoms, radii=radius_input)
            os.chdir(xyz_dir)
            sterimol_df.to_csv('xyz_sterimol.csv')
            print(sterimol_df.head())
            print(f"Sterimol values saved to xyz_sterimol.csv in {xyz_dir}")

            another_dir = input("Do you want to input another directory? (yes/no): ").strip().lower()
            if another_dir != 'yes':
                break

            
    elif args.command == "model":
        # Prepare environment/import path & working dir
        n_jobs = resolve_n_jobs(args.n_jobs)
        if n_jobs == -1:
            n_jobs = int(psutil.cpu_count(logical=True))
        add_modeling_path()
        # safe_chdir(args.workdir, args.features_csv)

        elapsed = run_experiment(
            mode=args.mode,
            features_csv=args.features_csv,
            target_csv=args.target_csv,
            y_value=args.y_value,
            n_jobs=n_jobs,
            min_features_num=args.min_features,
            max_features_num=args.max_features,
            threshold=args.threshold,
            top_n=args.top_n,
            leave_out=args.leave_out,
            bool_parallel=args.bool_parallel,
        )
        print(f"[model] Total execution time: {elapsed:.2f} seconds")





if __name__ == "__main__":
    main()






