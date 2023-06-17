.. _versions:

========
Versions
========

Version 1.5.1 [`url <https://github.com/jvalegre/aqme/releases/tag/1.5.1>`__]
   -  Added the verbose option. When verbose=False, no DAT and CSV files with summaries are printed
   -  Fixed bug for using constraints with SMILES that are not mapped
   -  Changed directory to create AQME-ROBERT databases to working dir (to adapt for full ML 
      workflows in ROBERT)
   -  Added SMARTS pattern recogniton to qdescp_atoms (the option now takes lists of atoms as well 
      as lists of functional groups)
   -  Fix minor bugs in QCORR-ORCA
   -  Added the geometry filter (geom option) to CSEARCH
   -  Added geom and complex_type as options in CSV columns
   -  Fixed QDESCP when using csv_name with columns that are not code_name and SMILES
   -  Added more RDKit descriptors in QDESCP with Descriptors.CalcMolDescriptors() (only for 
      rdkit>=2023)

Version 1.5.0 [`url <https://github.com/jvalegre/aqme/releases/tag/1.5.0>`__]
   -  A second ConstrainedEmbed() function with a core with no Hs was added to avoid
      RDKit embedding issues that show up in tricky cases
   -  Fixed the --charge and --mult options when using xyz/sdf/pdf/gjf/com files as inputs in 
      CSEARCH (xyz fixed in QPREP as well)
   -  Options low_check and e_threshold_qprep were added to QPREP (create inputs only for n 
      conformers of below a certain energy threshold)
   -  Option nodup_check was added to QCORR (turns off the duplicate filter)
   -  DBSTEP buried volume added in QDESCP with the qdescp_atoms option
   -  Atomic properties of a single atom type were added in QDESCP with the qdescp_atoms option
   -  Creation of databases for AQME-ROBERT workflows with the --robert option (True by default)

Version 1.4.7 [`url <https://github.com/jvalegre/aqme/releases/tag/1.4.7>`__]
   -  QCORR is compatible with directories that contain a "." symbol  
   -  QCORR is compatible with ORCA calcs (it doesn't detect spin contamination yet)
   -  QCORR includes hessian calculations for calcs with extra imaginary frequencies by default 
      (new option to control this: im_freq_input)
   -  QCORR tries to fix SCF convergence issues in ORCA by adding the SlowConv keyword
   -  qm_end option is added after the genecp section in QPREP
   -  Fixed a bug when using the destination option in CMIN-xTB

Version 1.4.6 [`url <https://github.com/jvalegre/aqme/releases/tag/1.4.6>`__]
   -  The files and input options are compatible with partial PATHs, full PATHs, and direct names 
      from command lines and Jupyter Notebooks  
   -  The SUMM option was fixed in CSEARCH  
   -  The files and input options now tolerate PATHs that contain directories with "." characters

Version 1.4.5 [`url <https://github.com/jvalegre/aqme/releases/tag/1.4.5>`__]
   -  Suffix/prefix options work in CSEARCH, CMIN and QPREP  
   -  Automatic recognition of metals with the auto_metal_atom option  
   -  In QPREP, if qm_input starts with "p ", the Gaussian inputs starts with "#p"  
   -  CSEARCH-CREST updates the CREST outfile as the program calculates (not at the end only)  

Version 1.4.4 [`url <https://github.com/jvalegre/aqme/releases/tag/1.4.4>`__]
   -  When using a CSV as input, the user can specify charge and mult for each species by 
      using the charge/mult columns  
   -  QCORR now detects duplicates including the successful calculations from previous runs  
   -  Fixed an error in full_check from QCORR when using genecp  
   -  Admits lists in command lines specified as ["X"], "[X]" and '["X"]'  

Version 1.4.3 [`url <https://github.com/jvalegre/aqme/releases/tag/1.4.3>`__]
   -  Return metal into RDKit mol object when using the metal_atoms option with CSEARCH-CREST  
   -  Doubles bonds do not add extra charges in metal complexes when using the automated charge 
      calculation from SMILES  
   -  Deprotonated SiR3 groups add -1 charge to metal complexes when using the automated charge 
      calculation from SMILES  

Version 1.4.2 [`url <https://github.com/jvalegre/aqme/releases/tag/1.4.2>`__]
   -  Fixed an error that raised when using CSEARCH-CREST with organic molecules  
   -  Adding more information printed when running CSEARCH  
   -  Updated README with citations from external programs  
   -  Fixed a bug during filtering of xTB conformers in CMIN (using kcal/mol instead of Hartree
      in the filters now)  
   -  Writing CSEARCH-CREST conformers in kcal/mol instead of Hartrees  
   -  Templates are not active when using metals with different number of ligands 
      (i.e. if complex_type='linear' and Cu2+/CuL2 are used simultaneously)  
   -  Fixed squarepyramidal templates  

Version 1.4.1 [`url <https://github.com/jvalegre/aqme/releases/tag/1.4.1>`__]
   -  Changed the way xTB works in CMIN. Before, it worked through xtb-python, but in this 
      version xtb is called through the xTB external command. This change speeds up the 
      calculations and avoids problems for people that do not have xtb-python installed.  
   -  Fixed some bugs in the PATHs when using AQME through command lines  
   -  Updated information printed in QDESCP  
   -  Adding more error prints when no program or files are specified  

Version 1.4.0 [`url <https://github.com/jvalegre/aqme/releases/tag/1.4.0>`__]
   -  Fixed a bug in the automated charge and multiplicity detector for metal complexes  
   -  Adapted CREST workflows to work with metal templates  
   -  Refactored utils and rearrange files to meet code analyzer standards  
   -  The mol object that CREST uses as input now comes from the RDKit 
      conformer generator (otherwise, metal templates aren't applied and 
      stereochemistry information might be lost)  

Version 1.3.1 [`url <https://github.com/jvalegre/aqme/releases/tag/1.3.1>`__]
   -  Workflows were updated  
   -  Small fixes in CREST when using constraints  
   -  Readme was updated  
   -  GoodVibes added in installation requirements  

Version 1.3.0 [`url <https://github.com/jvalegre/aqme/releases/tag/1.3.0>`__]
   -  Publication version  

Version 1.2.0 [`url <https://github.com/jvalegre/aqme/releases/tag/1.2.0>`__]
   -  This version improves how AQME reads PATHs from arguments to make the program more robust  

Version 1.1.0 [`url <https://github.com/jvalegre/aqme/releases/tag/1.1.0>`__]
   -  Fixes pip install issue coming from older versions  

Version 1.0.0 [`url <https://github.com/jvalegre/aqme/releases/tag/1.0.0>`__]
   -  First official version of AQME ready to generate publication-quality results  
