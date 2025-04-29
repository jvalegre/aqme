from setuptools import setup, find_packages
version = "1.7.2"
setup(
    name="aqme",
    packages=find_packages(exclude=["tests"]),
    package_data={"aqme": ["templates/*"]},
    version=version,
    license="MIT",
    description="Automated Quantum Mechanical Environments",
    long_description="Documentation in Read The Docs: https://aqme.readthedocs.io",
    long_description_content_type="text/markdown",
    author="Juan V. Alegre Requena",
    author_email="jv.alegre@csic.es",
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
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
    ],
    install_requires=[
        "PyYAML",
        "pandas==2.2.3",
        "progress",
        "numpy==1.26.4",
        "rdkit==2024.3.3",
        "cclib==1.7.2",
        "cffi",
        "morfeus-ml==0.7.2"
    ],
    python_requires=">=3.10",
    include_package_data=True,
)
