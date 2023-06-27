from setuptools import setup, find_packages
version = "1.5.1"
setup(
    name="aqme",
    packages=find_packages(exclude=["tests"]),
    package_data={"aqme": ["templates/*"]},
    version=version,
    license="MIT",
    description="Automated Quantum Mechanical Environments",
    long_description="Documentation in Read The Docs: https://aqme.readthedocs.io",
    long_description_content_type="text/markdown",
    author="Shree Sowndarya S. V., Juan V. Alegre Requena",
    author_email="svss@colostate.edu, jvalegre@unizar.es",
    keywords=[
        "workflows",
        "computational chemistry",
        "conformational sampling",
        "cheminformatics",
        "quantum mechanics",
        "DFT",
        "automated",
    ],
    url="https://github.com/jvalegre/aqme",
    download_url=f"https://github.com/jvalegre/aqme/archive/refs/tags/{version}.tar.gz",
    classifiers=[
        "Development Status :: 5 - Production/Stable",  # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
        "Intended Audience :: Developers",  # Define that your audience are developers
        "Topic :: Software Development :: Build Tools",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.7",  # Specify which python versions you want to support
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
    ],
    install_requires=[
        "PyYAML",
        "pandas>=2.0.2",
        "progress",
        "ase",
        "numpy",
        "cclib",
        "matplotlib",
        "seaborn",
        "cffi",
        "goodvibes",
        "py3Dmol",
        "ipywidgets",
        "dbstep"
    ],
    python_requires=">=3.0",
    include_package_data=True,
)
