TS involving a trimolecular complex generated from SMILES
---------------------------------------------------------

In the following example we will: 

1) Use the atom ordering of the provided SMILES to set up the constraints.
2) Do a constrained conformational search (using CREST) to generate various 
   conformers  of the SN2 transition state involving Water, Cloride ion and 
   2-Fluoro-2-methylpropane.  
3) Generate ORCA input files for each of the conformers

Step 1: creating SMILES with predefined atom numbers and setting constrains
...........................................................................

.. code:: python

    import glob
    from rdkit import Chem
    from aqme.csearch import csearch
    from aqme.qprep import qprep
    
    smi = 'O.FC(C)(C)C.[Cl-]'
    mol = Chem.MolFromSmiles(smi)
    mol = Chem.AddHs(mol)
    for i,atom in enumerate(mol.GetAtoms()):
        atom.SetAtomMapNum(i+3) 
    # mapped SMILES to use in CSEARCH
    smi_new = Chem.MolToSmiles(mol)
    
    print(smi_new)
    mol

Based on the numbers above, we choose the constraints for the TS

.. code:: python

    # 1) Bond between atoms 4 and 5 with a distance of 1.8 A
    # 2) Bond between atoms 5 and 9 with a distance of 1.8 A
    constraits_dist = [[4,5,1.8],[5,9,1.8]]
    
    # 3) Angle between atoms 4, 5 and 9 of 180 degrees
    constraits_angle = [[4,5,9,180]]

Step 2: constrained CSEARCH conformational sampling (creates SDF files)
.......................................................................

.. code:: python

    # run CSEARCH conformational sampling, specifying:
    # 1) Mapped SMILES string (smi=smi_new)
    # 2) Name for the output SDF files (name='TS-example')
    # 3) RDKit sampling (program='crest')
    # 4) Include CREGEN post-analysis (cregen=True)
    # 5) Specify that this a TS calculation (ts_complex=True)
    # 6) Define distance constraints (constraints_dist=constraits_dist)
    # 7) Define angle constraints (constraints_angle=constraits_angle)
    csearch(smi=smi_new,
            name='TS-example',program='crest',cregen=True,crest_nci=True,
            constraints_dist=constraits_dist,constraints_angle=constraits_angle)

Step 3: Writing Gaussian input files with the SDF files obtained from CSEARCH
.............................................................................

.. code:: python

    # set SDF filenames and directory where the new com files will be created
    sdf_rdkit_files = glob.glob(f'CSEARCH/*.sdf')
    
    # run QPREP input files generator, with:
    # 1) Files to convert (files=sdf_rdkit_files)
    # 2) QM program for the input (program='gaussian')
    # 3) Keyword line for the Gaussian inputs (qm_input='wb97xd/6-31+G* opt=(calcfc,ts,noeigen) freq')
    # 4) Memory to use in the calculations (mem='24GB')
    # 5) Processors to use in the calcs (nprocs=8)
    qprep(files=sdf_rdkit_files,program='gaussian',
            qm_input='wb97xd/6-31+G* opt=(calcfc,ts,noeigen) freq',mem='24GB',nprocs=8)
