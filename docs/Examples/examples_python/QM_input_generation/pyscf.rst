Using AQME to set up a PySCF calculation
========================================

Although support for PySCF is not yet included here we provide an example on how
to use the tools in AQME to run a PySCF job from json files. 

.. code:: python

    # run a PySCF single-point calculation  
    # with the DeepMind 21 (DM21) functional
    # using json files as input geometries 
    
    import os
    from aqme.utils import cclib_atoms_coords
    import json
    import glob
    
    # read json files
    initial_dir = os.getcwd()
    w_dir_main = os.getcwd()+'/QCORR_1/successful_QM_outputs/json_files'
    os.chdir(w_dir_main)
    json_files = glob.glob('*.json')
    
    # run the PySCF calculations
    for file in json_files:
        with open(file) as json_file:
            cclib_data = json.load(json_file)
    
        atom_types,cartesians = cclib_atoms_coords(cclib_data)
    
        coord_input = ''
        for i,atom in enumerate(atom_types):
            if i != 0:
                coord_input += ' '
            coord_input += atom+' '
                
            for j,cart in enumerate(cartesians[i]):
                coord_input += str(cart)
                if j != 2:
                    coord_input += ' '
                else:
                    if i != len(atom_types)-1:
                        coord_input += ';'
        
        charge = cclib_data['properties']['charge']
        mult = cclib_data['properties']['multiplicity']
        spin = mult-1
        basis = 'ccpvdz'
    
        # creates mol object for the calculations
        mol = gto.M(atom=coord_input, basis=basis, charge=charge, spin=spin)
        mol.output = f'./{file.split(".")[0]}.log'
        mol.verbose = 3
        mol.build()
    
        # runs the PySCF calculation
        if spin == 0:
            energy = mol.RKS().run(chkfile = 'expt0.chk', _numint = dm21.NeuralNumInt(dm21.Functional.DM21),
                                conv_tol = 1E-6, conv_tol_grad = 1E-3)
        else:
            energy = mol.UKS().run(chkfile = 'expt0.chk', _numint = dm21.NeuralNumInt(dm21.Functional.DM21),
                                conv_tol = 1E-6, conv_tol_grad = 1E-3) 
        
        # print results in the LOG file specified in mol.output                       
        energy.dump_scf_summary()
        energy.analyze()
        energy.spin_square()
    
    os.chdir(initial_dir)
