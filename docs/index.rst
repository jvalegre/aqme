=====================================
Welcome to pyCONFORT's documentation!
=====================================

This user guide focuses on the applications of pyCONFORT for creation of
conformers and provides wrapper for  analysis of Gaussian output files.

We provide a bunch examples depicting the uses of pyCONFORT in
organic molecules, metal complexes extended upto crystals.

All guides follow the a structure which allows the use of command
line for simple tasks and .yaml files for more detailed tasks.

Inputs structures for different sections are provided with examples and the respective outputs!


Please report any bugs or missing details at the mailing list or open an issue at `github`_.


.. _github: https://github.com/jvalegre/pyCONFORT

.. toctree::
	 :maxdepth: 3
	 :caption: Quickstart

	 Quickstart/setup
	 Quickstart/requirements
	 Quickstart/basic
	 Quickstart/development

.. toctree::
	:maxdepth: 3
	:caption: Conformer Generation

	Conformer Generation/Organic Molecules
	Conformer Generation/Metal Complex
	Conformer Generation/Crystals

.. toctree::
	:maxdepth: 3
	:caption: Templates

	Conformer Generation with Templates/Linear
	Conformer Generation with Templates/Trigonal Planar
	Conformer Generation with Templates/Square Planar
	Conformer Generation with Templates/Square Pyrimidal

.. toctree::
	:maxdepth: 3
	:caption: Methods

	Methods for Conformer Generation/RDKit Only
	Methods for Conformer Generation/RDKit + xTB
	Methods for Conformer Generation/RDKit + ANI1ccx

.. toctree::
	:maxdepth: 3
	:caption: Gaussian Input Files

	Generation of Gaussian Input Files/Options for input parameters


.. toctree::
	:maxdepth: 3
	:caption: Misc

	Misc/versions
	Misc/license
	Misc/help
