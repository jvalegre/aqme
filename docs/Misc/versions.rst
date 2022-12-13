.. _versions:

========
Versions
========

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
