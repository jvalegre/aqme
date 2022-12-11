.. |mol_3d| image:: ../../images/Et_rdkit.png
   :width: 300

Generate Gaussian Inputs
========================

For these examples we are going to assume that we have a folder named 'sdf_files'
that contains a single file 'ethane.sdf' with a single conformer in .sdf format 
whose gaussian input file we want to generate. As you might have guessed in this 
specific example we will be working with Ethane. 

|mol_3d|

The sdf file contents are as follows: 

.. literalinclude:: ../../chemfiles/ethane.sdf


.. note:: 
   
   The following code will also work for multiple conformers and/or molecules. 

.. note:: 

   AQME supports various formats for providing the geometries of the conformers.
   If we want to use a format to specify the molecule that does not contain 
   3D coordinates we will need to generate them beforehand, please see the 
   :doc:`Conformer Search <../conformer_search>` section.

We will be using the QPREP module :code:`--qprep`

We indicate the files whose gaussian input we want :code:`--files "sdf_files/*.sdf"` 

.. warning:: 

   Please notice that shell wildcard arguments need to be provided as strings.
   :code:`--files "sdf_files/*.sdf"` should be provided instead of 
   :code:`--files sdf_files/*.sdf`. This feature might change in future to 
   follow the usual conventions. 

We include the suffix that will be appended to the base name of the files 
:code:`--suffix "wb97xd-basic"`

We indicate the details of the gaussian calculation 
:code:`--qm_input "wb97xd/def2qzvpp scrf=(smd,solvent=acetonitrile)"`

And the number of processors :code:`--nprocs 8` and memory :code:`--mem 16GB` for
the gaussian calculations. 

Our final command will now look like this: 

.. code:: shell

   python -m aqme --qprep --files "sdf_files/*.sdf" --qm_input "wb97xd/def2qzvpp scrf=(smd,solvent=acetonitrile)" --suffix wb97xd-basic --program gaussian --mem 16GB --nprocs 8

With this we have generated a new folder named QCALC that contains the file 
'ethane_conf_1_wb97xd-basic.com' with the following contents:

.. highlight:: none 

.. literalinclude:: ../../chemfiles/ethane_basic.com

.. highlight:: default
   

Enforce Charge and Multiplicity
-------------------------------

If we had wanted to specify the charge and multiplicity we just need to add the 
appropriate keywords :code:`--charge 0 --mult 3`. 

Leading to the full command line (note that we have changed the suffix so the 
file name that we will generate will be different)

.. code:: shell

   python -m aqme --qprep --suffix wb97xd-triplet --charge 0 --mult 3 --files "sdf_files/*.sdf" --qm_input "wb97xd/def2qzvpp scrf=(smd,solvent=acetonitrile)" --program gaussian --mem 16GB --nprocs 8

Which will create the file 'ethane_conf_1_wb97xd-triplet.com' with the
following contents: 

.. highlight:: none

.. literalinclude:: ../../chemfiles/ethane_ch_mult.com

.. highlight:: default

Include the gen or genecp section
---------------------------------

If we want to include a genecp with automatic detection of the atoms in the 
molecule we have to add some extra keywords:

We need to specify the atoms whose ECP we are setting :code:`--gen_atoms "['C']"`

As well as the ECP :code:`--bs_gen def2svp`

The basis set for the atoms that will not use the ECP :code:`--bs_nogen "6-31G*"`

Finally we need to also substitute the basis set by ``genecp`` in the qm_input 
parameter :code:`--qm_input "wb97xd/genecp scrf=(smd,solvent=acetonitrile)"`

And we end up with the following command line: 

.. code:: shell

   python -m aqme --qprep --suffix "wb97xd-genecp" --gen_atoms "['C']" --bs_gen def2svp --bs_nogen "6-31G*" --files "sdf_files/*.sdf" --qm_input "wb97xd/genecp scrf=(smd,solvent=acetonitrile)" --program gaussian --mem 16GB --nprocs 8

Which will lead to the creation of the file 'ethane_conf_1_wb97xd-genecp.com' with the
following contents: 

.. highlight:: none

.. literalinclude:: ../../chemfiles/ethane_genecp.com

.. highlight:: default

If we instead do not want to include the ECP section, or in other words we want 
to use the ``gen`` instead of ``genecp`` we only need to substitute it in the 
``qm_input`` parameter. AQME will automatically recognize it and write 
the input file accordingly: 

.. code:: shell

   python -m aqme --qprep --suffix "wb97xd-gen" --gen_atoms "['C']" --bs_gen def2svp --bs_nogen "6-31G*" --files "sdf_files/*.sdf" --qm_input "wb97xd/gen scrf=(smd,solvent=acetonitrile)" --program gaussian --mem 16GB --nprocs 8

Which will lead to the creation of the file 'ethane_conf_1_wb97xd-gen.com' with the
following contents: 

.. highlight:: none

.. literalinclude:: ../../chemfiles/ethane_gen.com

.. highlight:: default


Include instructions after the geometry section
-----------------------------------------------

Finally if we want to specify some extra instructions after the geometry which 
are required for some commands such as the modredundant optimization option or 
for nbo6 calculations. Here we use the NBO as an example. 

We will only need to add a single extra command to include such instructions 
:code:`--qm_end "$nbo bndidx $end"` 

.. warning:: 

   In linux-based systems the :code:`$` needs to be escaped so the previous 
   option would need to be typed as :code:`--qm_end "\$nbo bndidx \$end"` 
   instead.  

But as we are using an NBO calculation as example we also need to use a gaussian
command line that is appropriate for the calculation 
:code:`--qm_input "pop=(nbo6read,savenbos) wb97xd/def2svp"`

Our final command line will look like: 

.. code:: shell

   python -m aqme --qprep --suffix wb97xd-nbo --qm_end "$nbo bndidx $end" --qm_input "pop=(nbo6read,savenbos) wb97xd/def2svp" --files "sdf_files/*.sdf"  --program gaussian --mem 16GB --nprocs 8

Which will lead to the creation of the file 'ethane_conf_1_wb97xd-nbo.com' with the
following contents: 

.. highlight:: none

.. literalinclude:: ../../chemfiles/ethane_nbo.com

.. highlight:: default
