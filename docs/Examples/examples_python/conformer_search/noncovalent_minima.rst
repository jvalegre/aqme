Noncovalent bimolecular complex generated from SMILES
-----------------------------------------------------

In the following example we will: 

1) Generate various conformers of the isopentane-water complex using CREST
2) Generate Gaussian input files for each of the conformers

Step 1: CSEARCH conformational sampling (creates SDF files)
...........................................................

.. code:: python

    import glob                                                                                                                                                     
    from aqme.csearch import csearch
    from aqme.qprep import qprep
    
    name = 'isopent-water-complex'
    smi = 'CCC(C)C.O'
    
    # run CSEARCH conformational sampling, specifying:
    # 1) SMILES string (smi=smi)
    # 2) Name for the output SDF files (name=name)
    # 3) CREST sampling (program='crest')
    # 4) Additional CREST keywords (crest_keywords='--nci')
    # 5) Include CREGEN post-analysis (cregen=True)
    # 6) Additional CREGEN keywords (cregen_keywords='--ewin 3')
    csearch(smi=smi,
            name=name,program='crest',crest_keywords='--nci',
            cregen=True,cregen_keywords='--ewin 3')

Step 2: Writing Gaussian input files with the sdf obtained from CSEARCH
.......................................................................

.. code:: python

    # set SDF filenames and directory where the new com files will be created
    sdf_rdkit_files = glob.glob(f'CSEARCH/*.sdf')
    
    # run QPREP input files generator, with:
    # 1) Files to convert (files=file)
    # 2) QM program for the input (program='gaussian')
    # 3) Keyword line for the Gaussian inputs (qm_input='wb97xd/6-31+G* opt freq')
    # 4) Memory to use in the calculations (mem='24GB')
    # 5) Processors to use in the calcs (nprocs=8)
    qprep(files=sdf_rdkit_files,program='gaussian',
            qm_input='wb97xd/6-31+G* opt freq',mem='24GB',nprocs=8)

