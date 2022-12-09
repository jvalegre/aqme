.. |mol_3d| image:: ../../images/Et_rdkit.png
   :width: 300

Generate ORCA Inputs
====================

.. warning::

   This section of the documentation is in construction so currently only a copy
   of the python example is displayed

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

As we will be generating qm inputs we will use the :code:`--qprep` module. 

We include the suffix that we want to append to the base name of the generated
files :code:`--suffix m06-basic`

We specify the files whose orca input we want :code:`--files "sdf_files/*.sdf"`

We include the number of processors :code:`--nprocs 8` and memory 
:code:`--mem 16GB` for the calculations 

And we specify the command line of the orca calculation: 

.. code:: shell 

   --qm_input "m06 def2qzvpp
    %cpcm
    smd true
    SMDsolvent \"CH2Cl2\"
    end"

.. note:: 

   the \"CH2Cl2\" in this case will ensure that the generated file contains :code:`"CH2Cl2"`

Our final command line will look like: 

.. code:: shell 

   python -m aqme --qprep --suffix m06-basic --files "sdf_files/*.sdf" --nprocs 8 --mem 16GB --qm_input "m06 def2qzvpp
   %cpcm
   smd true
   SMDsolvent \"CH2Cl2\"
   end"

With this we have generated a new folder named QCALC that contains the file 
'ethane_conf_1_m06-basic.inp' with the following contents:

.. literalinclude:: ../../chemfiles/ethane_basic.inp


Enforce Charge and Multiplicity
-------------------------------

If we had wanted to specify the charge and multiplicity we just need to add the 
appropriate keywords. 

.. code:: shell 

   python -m aqme --qprep --charge 0 --mult --3 --suffix "m06-triplet" --files "sdf_files/*.sdf" --nprocs 8 --mem 16GB --qm_input "m06 def2qzvpp
   %cpcm
   smd true
   SMDsolvent \"CH2Cl2\"
   end"

Will lead to the creation of the file 'ethane_conf_1_wb97xd-triplet.com' with the
following contents: 

.. literalinclude:: ../../chemfiles/ethane_triplet.inp




