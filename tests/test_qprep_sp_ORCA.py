#!/usr/bin/env python

######################################################.
# 		        Testing with pytest: 	             #
#     Conformer generation of organic molecules      #
######################################################.

import os
import pytest
from definitions_testing import conf_gen

# saves the working directory
path_organic = os.getcwd()

# tests for individual organic molecules and metal complexes
@pytest.mark.parametrize("folder, smiles, params_file, charge_orca, multiplicity_orca, param, job_type",
[
    # tests for ORCA in QPREP
    ('Organic_molecules', 'pentane.smi', 'params_test28.yaml', 0, 1, False, 'QPREP'), # ref test for ORCA (all default options)
    ('Organic_molecules', 'pentane.smi', 'params_test29.yaml', 0, 1, 'cpcm_input', 'QPREP'), # ref test for ORCA
    ('Organic_molecules', 'pentane.smi', 'params_test30.yaml', 0, 1, 'orca_scf_iters', 'QPREP'), # ref test for ORCA
    ('Organic_molecules', 'pentane.smi', 'params_test31.yaml', 0, 1, 'mdci_orca', 'QPREP'), # ref test for ORCA
    ('Organic_molecules', 'pentane.smi', 'params_test32.yaml', 0, 1, 'print_mini_orca', 'QPREP'), # ref test for ORCA
    ('Organic_molecules', 'pentane.smi', 'params_test33.yaml', 0, 1, 'set_input_line', 'QPREP'), # ref test for ORCA
    ('Organic_molecules', 'pentane.smi', 'params_test34.yaml', 0, 1, 'solvent_CPCM', 'QPREP'), # ref test for ORCA
    ('Organic_molecules', 'pentane.smi', 'params_test35.yaml', 0, 1, 'solvent_SMD', 'QPREP'), # ref test for ORCA
    ('Organic_molecules', 'pentane.smi', 'params_test36.yaml', 0, 1, 'nprocs', 'QPREP'), # ref test for ORCA
    ('Organic_molecules', 'pentane.smi', 'params_test37.yaml', 0, 1, 'mem', 'QPREP'), # ref test for ORCA

    ('Metal_complexes', 'Ir_hexacoord.smi', 'params_Ir_test6.yaml', 1, 1, 'QPREP'), # test with metals and genecp for ORCA
])


	parser.add_argument("--aux_atoms_orca",default=[], help="List of atoms included in the aux part when using multiple basis sets in ORCA", dest="aux_atoms_orca", type=str, nargs='?')
	parser.add_argument("--aux_basis_set_genecp_atoms",default=[], help="Auxiliary basis set for genecp/gen in ORCA", dest="aux_basis_set_genecp_atoms", type=str, nargs='?')
	parser.add_argument("--aux_fit_genecp_atoms",default=[], help="Fitting for the auxiliary basis set in ORCA (i.e. ['def2-TZVPP/C'])", dest="aux_fit_genecp_atoms", type=str, nargs='?')
	parser.add_argument("--cpcm_input",default=[], help="Additional lines for ORCA input files in the cpcm section", dest="cpcm_input", type=str, nargs='?')
	parser.add_argument("--orca_scf_iters",default=[], help="Number of SCF iterations in ORCA", dest="orca_scf_iters", type=str, nargs='?')
	parser.add_argument("--mdci_orca",default='None', help="mdci section in ORCA", dest="mdci_orca", type=str, nargs='?')
	parser.add_argument("--print_mini_orca",action="store_true",default=True, help="Option to print 'mini' (reduced outputs) in ORCA")
	parser.add_argument("--set_input_line",  help="(i) keywords used in Gaussian input files (overiding opt and freq) or (ii) additional keywords for the ORCA input line", default="None", dest="set_input_line")
	parser.add_argument("--solvent_model",  help="Type of solvent model in Gaussian and ORCA", default="gas_phase", dest="solvent_model", type=str)
	parser.add_argument("--solvent_name",  help="Name of the solvent in Gaussian and ORCA", default="Acetonitrile", dest="solvent_name", type=str)
	parser.add_argument("--nprocs", help="Number of processors for the DFT calculations", default=24, type=int, dest="nprocs")
	parser.add_argument("--mem", help="Memory for the DFT calculations (i) Gaussian: total memory; (ii) ORCA: memory per processor", default="96GB", type=str, dest="mem")

def test_confgen_organic(folder, smiles, params_file, n_confs_organic, prefilter_confs_rdkit_organic_organic, filter_confs_rdkit_organic, E_confs, charge_organic, multiplicity_organic, dihedral, xTB_ANI):
    # runs the program with the different tests
    cmd_organic = ['python', '-m', 'pyconfort', '--varfile', params_file]

    test_init_rdkit_confs_organic,test_prefilter_rdkit_confs_organic,test_filter_rdkit_confs_organic,round_confs_organic,test_round_confs_organic,test_charge_organic,test_unique_confs_organic,_,charge_organic_com,multiplicity_com_organic = conf_gen(path_organic, precision_organic, cmd_organic, folder, smiles, E_confs, dihedral, xTB_ANI, metal=False, template=False)

    # the assert statements are placed here, otherwise pytest doesn't explain the AssertionError
    # first, dicard tests 8 and 9 since they are designed to fail
    if n_confs_organic != 'nan':
        # dihedral vs no dihedral scans
        if not dihedral:
            assert str(n_confs_organic) == str(test_init_rdkit_confs_organic[0])
            assert str(prefilter_confs_rdkit_organic_organic) == str(test_prefilter_rdkit_confs_organic[0])
            assert str(filter_confs_rdkit_organic) == str(test_filter_rdkit_confs_organic[0])
        else:
            assert str(n_confs_organic) == str(test_init_rdkit_confs_organic[0])
            # I use the filter_confs_rdkit_organic variable to assert for unique confs in dihedral scan
            assert str(filter_confs_rdkit_organic) == str(test_unique_confs_organic[0])

        assert str(round_confs_organic) == str(test_round_confs_organic)
        assert str(charge_organic) == str(test_charge_organic[0])

        # make sure the COM files have the right charge_organic and multiplicity
        assert str(charge_organic_com) == str(charge_organic)
        assert str(multiplicity_com_organic) == str(multiplicity_organic)

    elif params_file == 'params_test8.yaml' or params_file == 'params_test9.yaml':
        assert str(test_filter_rdkit_confs_organic) == 'nan'
        assert str(test_round_confs_organic) == 'nan'

    else:
        assert 3 ==2
