from setuptools import setup, find_packages
"""
Setup script for the MolFeatures package.
This script handles the installation and distribution of the MolFeatures package.
It configures metadata such as package name, version, author information, and
dependencies. It also specifies which non-Python files should be included in the
distribution.
The package uses setuptools to find and include all Python packages automatically
through the find_packages() function.
Attributes:
    name (str): The name of the package as it will appear in PyPI.
    version (str): The current version of the package.
    packages (list): Automatically discovered Python packages.
    description (str): Short description of the package.
    long_description (str): Detailed description loaded from README.md file.
    long_description_content_type (str): Format of the long description (markdown).
    author (str): Package creator's name.
    author_email (str): Contact email address for the package author.
    url (str): Repository URL where the package source is hosted.
    install_requires (list): List of package dependencies.
    package_data (dict): Non-Python files to include in the package.
    include_package_data (bool): Whether to include additional data files.
Example:
    To install the package in development mode:
    ```
    pip install -e .
    ```
    To build a distribution package:
    ```
    python setup.py sdist bdist_wheel
    ```
"""

with open('README.md', 'r', encoding='utf-8') as fh:
    long_description = fh.read()
    print(long_description)

setup(
    name='MolFeatures',
    version='0.10005',
    packages=find_packages(),
    description='Your package description',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Eden Specktor',
    author_email='edenpsec@post.bgu.ac.il',
    url='https://github.com/edenspec2/Automation_code-main',
    install_requires=[],
    package_data={
        # If 'feather_example' is a Python package with __init__.py
        'MolFeatures': ['cube_example/*', 'logs_example/*','feather_example/*','Workshop_Example_Data/*','pictures/*', 'description.txt', 'README.md', 'requirements.txt','setup.py','Practical_Notebook_Modeling.ipynb','Practical_Notebook_Features.ipynb'],
    },
    include_package_data=True,
)

