#!/usr/bin/env python3

import subprocess
import sys
import os 

# List of packages to check/install
packages = [
    "pandas==2.1.4",
    "rdkit==2023.9.3",
    "python-igraph==0.11.3",
    "XlsxWriter==3.1.9",
    "ipywidgets==8.1.1",
    "pyarrow==14.0.2",
    "plotly==5.18.0",
    "customtkinter==5.2.1",
    "chardet==5.2.0",
    "matplotlib==3.8.2",
    "rmsd==1.5.1",
    "networkx==3.1",
    "dash==2.14.2",
    "pyvista==0.43.3",
    "pyvistaqt==0.11.0",
    "morfeus-ml==0.7.2",
    "scikit-learn==1.3.2",
    "seaborn==0.13.2",
    "Pillow==10.3.0",   # instead of "pillow"
    "scipy==1.11.4",
    "tqdm==4.66.1",
    "statsmodels==0.14.1",
    "adjustText==1.2.0",
    "multiprocess==0.70.17",
    "shap==0.48.0",
    "pymc==5.12.0",
    "tabulate==0.9.0",
    "ipykernel==6.28.0"
]

# (Optional) Map any package name to its import name if they differ
# For example, "PIL" -> "PIL" or "python-igraph" -> "igraph"
import_name_map = {
    "pillow": "PIL",
    "python-igraph": "igraph",
    "morfeus-ml": "morfeus"
    # Add more if needed
}


def base_package_name(pkg: str) -> str:
    """
    Return the base package name without any version specifier.
    E.g. 'chardet==5.2.0' -> 'chardet'
         'numpy>=1.23'   -> 'numpy'
    """
    for sep in ["==", ">=", "<=", "~=", ">", "<"]:
        if sep in pkg:
            return pkg.split(sep, 1)[0].strip()
    return pkg.strip()


def install_package(package_name):
    """
    Installs the given package using pip in the current Python environment.
    """
    print(f"\nAttempting to install '{package_name}' ...")
    result = subprocess.run(
        [sys.executable, "-m", "pip", "install", package_name],
        capture_output=True,
        text=True
    )
    if result.returncode == 0:
        print(f"✅ Installed successfully: {package_name}")
    else:
        print(f"❌ Error installing {package_name}: {result.stderr}")

def check_import(package_name):
    """
    Attempts to import a package. If unsuccessful, installs the package and re-attempts.
    """
    # Handle cases where pip name ≠ import name
    dist_name = base_package_name(package_name)
    print(dist_name,'heyy')
    mod_name = import_name_map.get(dist_name, package_name)

    try:
        __import__(mod_name)
        print(f"✔ Imported successfully: {dist_name} (import as '{mod_name}')")
    except ImportError:
        install_package(dist_name)
        # Retry import
        try:
            __import__(mod_name)
            # print(f"✔ Imported successfully after install: {dist_name} (import as '{mod_name}')")
        except ImportError as e:
            # print(f"❌ Failed to import {dist_name} even after installation. Error: {e}")
            pass

def print_environment_info():
    """
    Prints information about the Python environment.
    """
    print("\n🔎 Environment Info")
    print(f"Python executable : {sys.executable}")
    print(f"Python version    : {sys.version}")
    print(f"Working directory : {os.getcwd()}")
    print(f"PATH              : {os.environ.get('PATH', '').split(os.pathsep)[0]} ...")

def main():
    """
    Checks (and if necessary installs) all packages, then tries re-importing them.
    """
    print_environment_info()
    print("\nChecking required packages...\n")
    for pkg in packages:
        check_import(pkg)
    print("\n🎉 All done. If no errors were reported above, your environment is ready.")

if __name__ == "__main__":
    main()