from setuptools import setup, find_packages
version = "2.0.0"
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
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "Programming Language :: Python :: 3.14",
    ],
    install_requires=[
        "PyYAML==6.0.3",
        "pandas==2.3.3",
        "progress==1.6.1",
        "numpy==2.3.4",
        "rdkit==2025.9.1",
        "cclib==1.8.1",
        "cffi==2.0.0",
        "morfeus-ml==0.8.0",
    ],
    python_requires=">=3.11",
    include_package_data=True,
)
