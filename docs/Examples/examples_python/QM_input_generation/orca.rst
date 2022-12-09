.. |mol_3d| image:: ../../images/Et_rdkit.png
   :width: 300

Generate ORCA Inputs
====================

For these examples we are going to assume that we have a folder named 'sdf_files'
that contains a single file 'ethane.sdf' with a single conformer in .sdf format 
whose orca input file we want to generate. As you might have guessed in this 
specific example we will be working with Ethane. 

|mol_3d|

The sdf file contents are as follows: 

.. literalinclude:: ../../chemfiles/ethane.sdf


.. note:: 
   
   The following code will also work for multiple conformers and/or molecules. 


.. note:: 

   aqme supports various formats for providing the geometries of the conformers.
   If we want to use a format to specify the molecule that does not contain 
   3D coordinates we will need to generate them beforehand, please see the 
   :doc:`Conformer Search <../conformer_search>` section.


First we start importing the required modules. 

.. code:: python

   from pathlib import Path
   from aqme.qprep import qprep

Next we list all the files whose orca input we want. 

.. code:: python

    sdf_files = [str(filepath) for filepath in Path('sdf_files').glob('*.sdf'))]

Now we are going to specify the ORCA calculation. 

.. code:: python 

    ORCA_SP = r'''
    m06 def2qzvpp
    %cpcm
    smd true
    SMDsolvent "CH2Cl2"
    end'''.lstrip()

Now we proceed to generate the orca input files. 

.. code:: python

   qprep(files=sdf_files, 
         qm_input=ORCA_SP, 
         suffix='m06-basic',
         program='orca', 
         mem='16GB', 
         nprocs=8)

With this we have generated a new folder named QCALC that contains the file 
'ethane_conf_1_m06-basic.inp' with the following contents:

.. literalinclude:: ../../chemfiles/ethane_basic.inp


Enforce Charge and Multiplicity
-------------------------------

If we had wanted to specify the charge and multiplicity we just need to add the 
appropriate keywords. 

.. code:: python

    qprep(files=sdf_files, 
          charge=0, 
          mult=3,
          qm_input=ORCA_SP,
          suffix='m06-reduced',
          program='orca', 
          mem='16GB', 
          nprocs=8)

Will lead to the creation of the file 'ethane_conf_1_wb97xd-triplet.com' with the
following contents: 

.. literalinclude:: ../../chemfiles/ethane_triplet.inp




