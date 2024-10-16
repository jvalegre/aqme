#!/usr/bin/env python

######################################################.
# 		        Testing with pytest: 	             #
#                 CSEARCH module                     #
######################################################.

import os
import pytest
import glob
from aqme.csearch import csearch
import rdkit
import shutil

# saves the working directory
w_dir_main = os.getcwd()
csearch_methods_dir = w_dir_main + "/tests/csearch_methods"
csearch_rdkit_summ_dir = w_dir_main + "/tests/csearch_rdkit_summ"
csearch_fullmonte_dir = w_dir_main + "/tests/csearch_fullmonte"
csearch_crest_dir = w_dir_main + "/tests/csearch_crest"
csearch_others_dir = w_dir_main + "/tests/csearch_others"
csearch_input_dir = w_dir_main + "/tests/csearch_input"
csearch_varfile_dir = w_dir_main + "/tests/csearch_varfile"

if not os.path.exists(csearch_methods_dir):
    os.mkdir(csearch_methods_dir)
if not os.path.exists(csearch_rdkit_summ_dir):
    os.mkdir(csearch_rdkit_summ_dir)
if not os.path.exists(csearch_fullmonte_dir):
    os.mkdir(csearch_fullmonte_dir)
if not os.path.exists(csearch_crest_dir):
    os.mkdir(csearch_crest_dir)
if not os.path.exists(csearch_others_dir):
    os.mkdir(csearch_others_dir)
if not os.path.exists(csearch_input_dir):
    os.mkdir(csearch_input_dir)
if not os.path.exists(csearch_varfile_dir):
    os.mkdir(csearch_varfile_dir)

# tests for varfile
@pytest.mark.parametrize(
    "varfile, nameinvarfile, output_nummols",
    [
        # tests for conformer generation with RDKit
        ("params.yaml", "pentane_varfile", 2),
    ],
)
def test_csearch_varfile(varfile, nameinvarfile, output_nummols):
    os.chdir(csearch_varfile_dir)
    # runs the program with the different tests
    csearch(w_dir_main=csearch_varfile_dir, varfile=varfile)

    # tests here
    file = str("CSEARCH/" + nameinvarfile + "_rdkit" + ".sdf")
    mol1 = rdkit.Chem.SDMolSupplier(file, removeHs=False)
    assert len(mol1) == output_nummols
    os.chdir(w_dir_main)


# tests for input types
@pytest.mark.parametrize(
    "program, input, output_nummols",
    [
        # tests for conformer generation with RDKit
        ("rdkit", "pentane.smi", [2, 4]),
        ("rdkit", "pentane.csv", [2, 4]),
        ("rdkit", "Cu.csv", [1, 1, 1]),
        ("rdkit", "blank_smi.csv", [1, 1]), # check blank cells
        ("rdkit", "unique_smi.csv", [1]), # check that only unique SMILES are used
        ("rdkit", "partial_path", [2, 4]), # checks partial PATHs
        ("rdkit", "file_name", [2, 4]), # checks file_name
        ("rdkit", "molecules.cdx", [4, 2]),
        ("rdkit", "pentane_gjf.gjf", 4),
        ("rdkit", "pentane_com.com", 4),
        ("rdkit", "pentane_xyz.xyz", 4),
        ("rdkit", "pentane_sdf.sdf", 4),
        ("rdkit", "pentane_mol.mol", 4),
        ("rdkit", "pentane_pdb.pdb", 4),
        # ("rdkit", "pentane_mol2.mol2", None ), # not working currently in rdkit
    ],
)
def test_csearch_input_parameters(program, input, output_nummols):
    
    # runs the program with the different tests
    os.chdir(w_dir_main)
    if input == "partial_path":
        input = "tests/csearch_input/pentane.csv"
        csearch(destination=f'{csearch_input_dir}/CSEARCH', program=program, input=input)
        input = "pentane.csv"
        os.chdir(csearch_input_dir)
    elif input == "file_name":
        input = "pentane.csv"
        os.chdir(csearch_input_dir)
        csearch(destination=f'{csearch_input_dir}/CSEARCH', program=program, input=input)
    elif input == "Cu.csv":
        csearch(destination=f'{csearch_input_dir}/CSEARCH', program=program, input=f'{csearch_input_dir}/{input}', sample=5)
        os.chdir(csearch_input_dir)
    else:
        csearch(destination=f'{csearch_input_dir}/CSEARCH', program=program, input=f'{csearch_input_dir}/{input}')
        os.chdir(csearch_input_dir)
    
    # tests here
    if input in ["pentane.smi", "pentane.csv"]:
        file1 = f'{csearch_input_dir}/CSEARCH/butane_{input.split(".")[1]}_{program}.sdf'
        file2 = f'{csearch_input_dir}/CSEARCH/pentane_{input.split(".")[1]}_{program}.sdf'

        with rdkit.Chem.SDMolSupplier(file1, removeHs=False) as mol1:
            assert len(mol1) == output_nummols[0]
        with rdkit.Chem.SDMolSupplier(file2, removeHs=False) as mol2:
            assert len(mol2) == output_nummols[1]
        os.remove(file1)
        os.remove(file2)

    elif input == "Cu.csv":
        file1 = f'{csearch_input_dir}/CSEARCH/cu1_{program}.sdf'
        file2 = f'{csearch_input_dir}/CSEARCH/cu2_{program}.sdf'
        file3 = f'{csearch_input_dir}/CSEARCH/cu3_{program}.sdf'

        with rdkit.Chem.SDMolSupplier(file1, removeHs=False) as mol1:
            assert len(mol1) == output_nummols[0]
        with rdkit.Chem.SDMolSupplier(file2, removeHs=False) as mol2:
            assert len(mol2) == output_nummols[1]
        with rdkit.Chem.SDMolSupplier(file3, removeHs=False) as mol3:
            assert len(mol3) == output_nummols[2]
        os.remove(file1)
        os.remove(file2)
        os.remove(file3)

    elif input == "blank_smi.csv":
        file1 = f'{csearch_input_dir}/CSEARCH/Me_{program}.sdf'
        file2 = f'{csearch_input_dir}/CSEARCH/Et_{program}.sdf'

        with rdkit.Chem.SDMolSupplier(file1, removeHs=False) as mol1:
            assert len(mol1) == output_nummols[0]
        with rdkit.Chem.SDMolSupplier(file2, removeHs=False) as mol2:
            assert len(mol2) == output_nummols[1]
        os.remove(file1)
        os.remove(file2)

    elif input == "unique_smi.csv":
        file1 = f'{csearch_input_dir}/CSEARCH/Me1_{program}.sdf'
        file2 = f'{csearch_input_dir}/CSEARCH/Me2_{program}.sdf'

        assert os.path.exists(file1)
        assert not os.path.exists(file2)
        os.remove(file1)

        # ensure that the warning is shown in the DAT file
        file_dat = str(w_dir_main+f"/CSEARCH_data.dat")
        outfile = open(file_dat, "r")
        outlines_dat = outfile.readlines()
        outfile.close()

        find_warn_dup = False
        for line in outlines_dat:
            if 'x  SMILES "C" used in Me2 was already used with a different code_name!' in line:
                find_warn_dup = True
        assert find_warn_dup

    elif input in ["molecules.cdx"]:
        file1 = f'{csearch_input_dir}/CSEARCH/molecules_0_{program}.sdf'
        file2 = f'{csearch_input_dir}/CSEARCH/molecules_1_{program}.sdf'
        mol1 = rdkit.Chem.SDMolSupplier(file1, removeHs=False)
        mol2 = rdkit.Chem.SDMolSupplier(file2, removeHs=False)
        assert len(mol1) == output_nummols[0]
        assert len(mol2) == output_nummols[1]
    else:
        file = f'{csearch_input_dir}/CSEARCH/{input.split(".")[0]}_{program}.sdf'
        mols = rdkit.Chem.SDMolSupplier(file, removeHs=False)
        assert len(mols) == output_nummols
    os.chdir(w_dir_main)


# tests for parameters of SUMM
@pytest.mark.parametrize(
    "program, smi, name, charge, mult, ang_summ, output_nummols",
    [
        ("summ", "CCCCC", "pentane_summ", 3, 4, 120, 4),
    ],
)
def test_csearch_summ_parameters(
    program,
    smi,
    name,
    charge,
    mult,
    ang_summ,
    output_nummols,
):
    os.chdir(csearch_rdkit_summ_dir)
    # runs the program with the different tests
    csearch(
        program=program,
        smi=smi,
        name=name,
        charge=charge,
        mult=mult,
        degree=ang_summ,
    )

    # tests here
    file = str("CSEARCH/" + name + "_" + program + ".sdf")
    mols = rdkit.Chem.SDMolSupplier(file, removeHs=False)
    assert len(mols) == output_nummols
    assert charge == int(mols[0].GetProp("Real charge"))
    assert mult == int(mols[0].GetProp("Mult"))
    os.chdir(w_dir_main)


# tests for parameters of CREST
@pytest.mark.parametrize(
    "program, smi, name, cregen, cregen_keywords, crest_keywords, charge, mult, output_nummols",
    [
        # tests for conformer generation with CREST and CREGEN
        (
            "crest",
            "C",
            "methane",
            True,
            "--ethr 1 --rthr 0.5 --bthr 0.7 --ewin 4 --cluster",
            None,
            0,
            1,
            1,
        ),
        # tests for conformer generation with multiple confs
        (
            "crest",
            "CCCC",
            "butane",
            False,
            None,
            None,
            0,
            1,
            3,
        ),
        # tests for crest_keywords
        (
            "crest",
            "C",
            "methane_solvent",
            True,
            "--ethr 1 --rthr 0.5 --bthr 0.7 --ewin 4 --cluster",
            "--alpb benzene",
            0,
            1,
            1,
        ),
        # tests for n of processors
        (
            "crest",
            "C",
            "methane_nprocs",
            True,
            "--ethr 1 --rthr 0.5 --bthr 0.7 --ewin 4 --cluster",
            None,
            0,
            1,
            1,
        ),
        # test for charge and mult
        (
            "crest",
            "C",
            "methane_charged",
            True,
            "--ethr 1 --rthr 0.5 --bthr 0.7 --ewin 4 --cluster",
            None,
            1,
            2,
            1,
        ),
    ],
)
def test_csearch_crest_parameters(
    program,
    smi,
    name,
    cregen,
    cregen_keywords,
    crest_keywords,
    charge,
    mult, 
    output_nummols,
):
    os.chdir(csearch_crest_dir)

    # runs the program with the different tests
    if name == 'methane_nprocs':
        csearch(
            w_dir_main=csearch_crest_dir,
            program=program,
            smi=smi,
            name=name,
            cregen=cregen,
            cregen_keywords=cregen_keywords,
            crest_keywords=crest_keywords,
            charge=charge,
            mult=mult,
            nprocs=14
        )
    elif name == 'butane':
        csearch(
            w_dir_main=csearch_crest_dir,
            program=program,
            smi=smi,
            name=name
        )
    elif name == 'methane_solvent':
        csearch(
            w_dir_main=csearch_crest_dir,
            program=program,
            smi=smi,
            name=name,
            cregen=cregen,
            cregen_keywords=cregen_keywords,
            crest_keywords=crest_keywords,
            charge=charge,
            mult=mult,
            xtb_keywords='--alpb benzene'
        )        
    else:
        csearch(
            w_dir_main=csearch_crest_dir,
            program=program,
            smi=smi,
            name=name,
            cregen=cregen,
            cregen_keywords=cregen_keywords,
            crest_keywords=crest_keywords,
            charge=charge,
            mult=mult
        )

    if name != 'butane':
        file_cluster = str(
            csearch_crest_dir
            + "/CSEARCH/"
            + program
            + "_xyz/"
            + name + "_"
            + "crest_clustered.xyz"
        )
        assert os.path.exists(file_cluster)
    file = str("CSEARCH/" + name + "_" + program + ".sdf")
    mols = rdkit.Chem.SDMolSupplier(file, removeHs=False)
    assert charge == int(mols[0].GetProp("Real charge"))
    assert mult == int(mols[0].GetProp("Mult"))
    if name == 'butane': # CREST sometimes gives 2 conformers and other times 3
        assert len(mols) >= 1
    else:
        assert len(mols) == output_nummols
    os.chdir(w_dir_main)

    file_crest = str(csearch_crest_dir+f"/CSEARCH/crest_xyz/{name}_crest.out")
    outfile = open(file_crest, "r")
    outlines_crest = outfile.readlines()
    outfile.close()
    if name == 'methane_charged':
        for line in outlines_crest:
            if line.startswith(' > crest'):
                # check if charge and mult are correct in CREST
                assert line.find('--chrg 1 --uhf 1') > -1
                break

        file_xtb2 = str(csearch_crest_dir+"/CSEARCH/crest_xyz/methane_charged_crest_xtb1.out")
        outfile_xtb2 = open(file_xtb2, "r")
        outlines_xtb2 = outfile_xtb2.readlines()
        outfile.close()
        for _,line in enumerate(outlines_xtb2):
            if line.find('program call') > -1:
                # check if charge and mult are correct in xTB
                assert line.find('-c 1 --uhf 1') > -1
    elif name == 'methane_solvent':
        for line in outlines_crest:
            if line.startswith(' > crest'):
                # check if crest_keywords are correct in CREST
                assert line.find('--alpb benzene') > -1
                break
        file_xtb2 = str(csearch_crest_dir+"/CSEARCH/crest_xyz/methane_solvent_crest_xtb1.out")
        outfile_xtb2 = open(file_xtb2, "r")
        outlines_xtb2 = outfile_xtb2.readlines()
        outfile.close()
        for _,line in enumerate(outlines_xtb2):
            if line.find('program call') > -1:
                # check if xtb_keywords are correct in CREST
                assert line.find('--alpb benzene') > -1
    elif name == 'methane_nprocs':
        for line in outlines_crest:
            if line.startswith(' > crest'):
                # check if n of procs are correct in CREST
                assert line.find('-T 14') > -1
                break
        file_xtb2 = str(csearch_crest_dir+"/CSEARCH/crest_xyz/methane_nprocs_crest_xtb1.out")
        outfile_xtb2 = open(file_xtb2, "r")
        outlines_xtb2 = outfile_xtb2.readlines()
        outfile.close()
        for _,line in enumerate(outlines_xtb2):
            if line.find('program call') > -1:
                # check if xtb_keywords are correct in CREST
                assert line.find('-P 14') > -1


# tests for parameters of csearch fullmonte
@pytest.mark.parametrize(
    "program, smi, name, charge, mult, ewin_fullmonte, ewin_sample_fullmonte, nsteps_fullmonte, nrot_fullmonte, ang_fullmonte, output_nummols",
    [
        ("fullmonte", "CCCCC", "pentane_fullmonte", 3, 4, 12, 3, 200, 4, 10, 4),
    ],
)
def test_csearch_fullmonte_parameters(
    program,
    smi,
    name,
    charge,
    mult,
    ewin_fullmonte,
    ewin_sample_fullmonte,
    nsteps_fullmonte,
    nrot_fullmonte,
    ang_fullmonte,
    output_nummols,
):
    os.chdir(csearch_fullmonte_dir)
    # runs the program with the different tests
    csearch(
        w_dir_main=csearch_fullmonte_dir,
        program=program,
        smi=smi,
        name=name,
        charge=charge,
        mult=mult,
        ewin_fullmonte=ewin_fullmonte,
        ewin_sample_fullmonte=ewin_sample_fullmonte,
        nsteps_fullmonte=nsteps_fullmonte,
        nrot_fullmonte=nrot_fullmonte,
        ang_fullmonte=ang_fullmonte,
    )

    # tests here
    file = str("CSEARCH/" + name + "_" + program + ".sdf")
    mols = rdkit.Chem.SDMolSupplier(file, removeHs=False)
    assert len(mols) == output_nummols
    assert charge == int(mols[0].GetProp("Real charge"))
    assert mult == int(mols[0].GetProp("Mult"))
    os.chdir(w_dir_main)


# tests for parameters of csearch rdkit
@pytest.mark.parametrize(
    "program, smi, name, charge, mult, sample, opt_steps_rdkit, heavyonly, ewin_csearch, initial_energy_threshold, energy_threshold, rms_threshold, output_nummols ",
    [
        # tests for conformer generation with RDKit
        (
            "rdkit",
            "CCCCC",
            "pentane_rdkit",
            3,
            4,
            "auto",
            1000,
            False,
            10,
            0.0001,
            2,
            0.2,
            4,
        ),
        (
            "rdkit",
            "CC[CH]CC",
            "radical_rdkit",
            0,
            2,
            15,
            10,
            True,
            10,
            0.000001,
            0.6,
            0.3,
            24,
        ),
        (
            "rdkit",
            "C[NH2+]CC",
            "charged_rdkit",
            1,
            1,
            "auto",
            10,
            True,
            10,
            0.0001,
            4,
            0.6,
            3,
        ),
    ],
)
def test_csearch_rdkit_parameters(
    program,
    smi,
    name,
    charge,
    mult,
    sample,
    opt_steps_rdkit,
    heavyonly,
    ewin_csearch,
    initial_energy_threshold,
    energy_threshold,
    rms_threshold,
    output_nummols,
):
    os.chdir(csearch_rdkit_summ_dir)
    # runs the program with the different tests
    if name != 'radical_rdkit':
        csearch(
            w_dir_main=csearch_rdkit_summ_dir,
            program=program,
            smi=smi,
            name=name,
            charge=charge,
            mult=mult,
            opt_steps_rdkit=opt_steps_rdkit,
            heavyonly=heavyonly,
            ewin_csearch=ewin_csearch,
            initial_energy_threshold=initial_energy_threshold,
            energy_threshold=energy_threshold,
            rms_threshold=rms_threshold,
        )
    else:
        csearch(
            w_dir_main=csearch_rdkit_summ_dir,
            program=program,
            smi=smi,
            name=name,
            charge=charge,
            mult=mult,
            sample=sample,
            opt_steps_rdkit=opt_steps_rdkit,
            heavyonly=heavyonly,
            ewin_csearch=ewin_csearch,
            initial_energy_threshold=initial_energy_threshold,
            energy_threshold=energy_threshold,
            rms_threshold=rms_threshold,
            pytest_testing=True
        )

    # tests here
    if name != 'radical_rdkit':
        file = str("CSEARCH/" + name + "_" + program + ".sdf")
    else: # checks all the conformers initially generated
        file = str("CSEARCH/" + name + "_" + program + "_all_confs.sdf")
    mols = rdkit.Chem.SDMolSupplier(file, removeHs=False)
    assert len(mols) == output_nummols
    assert charge == int(mols[0].GetProp("Real charge"))
    assert mult == int(mols[0].GetProp("Mult"))
    # check if all the energies are sorted
    mol_energies = []
    for mol in mols:
        mol_energies.append(float(mol.GetProp("Energy")))
    assert mol_energies == sorted(mol_energies)

    if name == 'radical_rdkit': # test clusterized conformers
        initial_three_E = []
        for mol in mols:
            initial_three_E.append(float(mol.GetProp("Energy")))
            if len(initial_three_E) == 3:
                break
        file = str("CSEARCH/" + name + "_" + program + ".sdf")
        mols = rdkit.Chem.SDMolSupplier(file, removeHs=False)
        final_three_E = []
        for mol in mols:
            final_three_E.append(float(mol.GetProp("Energy")))
            if len(final_three_E) == 3:
                break
        assert len(mols) == sample
        assert charge == int(mols[0].GetProp("Real charge"))
        assert mult == int(mols[0].GetProp("Mult"))
        # check if the first five energies are included
        assert initial_three_E == final_three_E
        # check if all the energies are sorted
        mol_energies = []
        for mol in mols:
            mol_energies.append(float(mol.GetProp("Energy")))
        assert mol_energies == sorted(mol_energies)
        
        # check if the DAT file contains the Butina print
        dat_rdkit = str(csearch_rdkit_summ_dir+f"/CSEARCH_data.dat")
        datfile = open(dat_rdkit, "r")
        datfile_rdkit = datfile.readlines()
        datfile.close()
        butina_found = False
        for line in datfile_rdkit:
            if 'conformers using a combination of energies and Butina RMS-based clustering' in line:
                butina_found = True
        assert butina_found

    # check that H atoms are included
    outfile = open(file,"r")
    outlines = outfile.readlines()
    outfile.close()
    Hatoms_found = False
    for line in outlines:
        if "H   0" in line:
            Hatoms_found = True
    assert Hatoms_found    
    os.chdir(w_dir_main)


# tests for individual organic molecules and metal complexes with different types of csearch methods
@pytest.mark.parametrize(
    "program, smi, name, complex, metal_complex, complex_type, constraints_dist, constraints_angle, constraints_dihedral, charge, mult, crest_keywords, destination, output_nummols",
    [
        # tests for conformer generation with RDKit, SUMM, FullMonte and CREST
        (
            "rdkit",
            "CCCCC",
            "pentane_RD",
            False,
            False,
            None,
            [],
            [],
            [],
            0,
            1,
            None,
            False,
            4,
        ),
        # (
        #     "summ",
        #     "CCCCC",
        #     "pentane",
        #     False,
        #     False,
        #     None,
        #     [],
        #     [],
        #     [],
        #     0,
        #     1,
        #     None,
        #     False,
        #     4,
        # ),
        (
            "fullmonte",
            "CCCCC",
            "pentane_FM",
            False,
            False,
            None,
            [],
            [],
            [],
            0,
            1,
            None,
            False,
            4,
        ),
        # metal atoms
        (
            "rdkit",
            "I[Pd]([PH3])(F)Cl",
            "Pd_metal_only",
            False,
            True,
            None,
            [],
            [],
            [],
            -1,
            1,
            None,
            False,
            1
        ),
        # multiple templates
        (
            "rdkit",
            "I[Pd]([PH3])(F)Cl",
            "Pd_complex",
            False,
            True,
            "squareplanar",
            [],
            [],
            [],
            -1,
            1,
            None,
            False,
            1
        ),
        (
            "rdkit",
            "N#C[Cu](C#N)C#N",
            "Cu_trigonal",
            False,
            True,
            "trigonalplanar",
            [],
            [],
            [],
            -2,
            1,
            None,
            False,
            1
        ),
        (
            "rdkit",
            "[O-][V](Cl)(Cl)(Cl)Cl",
            "V_squarepyramidal",
            False,
            True,
            "squarepyramidal",
            [],
            [],
            [],
            -2,
            1,
            None,
            False,
            1
        ),
        (
            "rdkit",
            "rule_IrSP.csv",
            "rule_IrSP",
            False,
            True,
            "Ir_squareplanar",
            [],
            [],
            [],
            0,
            1,
            None,
            False,
            3
        ),

        # compatibility of CREST with metal complexes and templates
        (
            "crest",
            "[NH3+][Ag][NH3+]",
            "Ag_complex_crest",
            False,
            True,
            "linear",
            [],
            [],
            [],
            1,
            1,
            None,
            False,
            1,
        ),
        (
            "crest",
            "[NH3][Ag][NH3]",
            "Ag_metal_crest",
            False,
            True,
            None,
            [],
            [],
            [],
            1,
            1,
            None,
            False,
            1,
        ),
        # compatibility of CREST with destination
        (
            "crest",
            "CC",
            "ethane",
            False,
            False,
            None,
            [],
            [],
            [],
            0,
            1,
            None,
            True,
            1,
        ),
        # without --nci, the resulting sdf has more than 400 conformers
        (
            "crest",
            "C.O",
            "nci_keyword",
            True,
            False,
            None,
            [],
            [],
            [],
            0,
            1,
            '--nci --cbonds 0.5',
            False,
            None,
        ),
        (
            "crest",
            "C.O",
            "cregen_clustering",
            True,
            False,
            None,
            [],
            [],
            [],
            0,
            1,
            '--nci --cbonds 0.5',
            False,
            None,
        ),
        (
            "crest",
            "C.O",
            "sample_keyword",
            True,
            False,
            None,
            [],
            [],
            [],
            0,
            1,
            '--nci --cbonds 0.5',
            False,
            None,
        ),
        (
            "crest",
            "[Cl-:9].[F:4][C:5]([H:15])([H:16])[H:17]",
            "ts",
            True,
            False,
            None,
            [[4, 5, 1.8], [5, 9, 1.8]],
            [[4, 5, 9, 180]],
            [],
            -1,
            1,
            None,
            False,
            1,
        ),
    ],
)
def test_csearch_methods(
    program,
    smi,
    name,
    complex,
    metal_complex,
    complex_type,
    constraints_dist,
    constraints_angle,
    constraints_dihedral,
    charge,
    mult,
    crest_keywords,
    destination,
    output_nummols,
):
    os.chdir(csearch_methods_dir)
    # runs the program with the different tests
    if destination:
        csearch(
            w_dir_main=csearch_methods_dir,
            destination=csearch_methods_dir+'/Et_sdf_files',
            program=program,
            smi=smi,
            name=name
        )
        
    elif not complex and not metal_complex:
        csearch(
            w_dir_main=csearch_methods_dir,
            program=program,
            smi=smi,
            name=name,
        )

    elif metal_complex is True:
        if name in ['Pd_metal_only','Ag_metal_crest']:
            csearch(
                w_dir_main=csearch_methods_dir,
                program=program,
                smi=smi,
                name=name,
                mult=mult,
                charge=charge,
                sample=10
            )
        elif name == 'rule_IrSP':
            csearch(
                program=program,
                input='rule_IrSP.csv',
                sample=10
            )
        else:
            csearch(
                w_dir_main=csearch_methods_dir,
                program=program,
                smi=smi,
                name=name,
                charge=charge,
                complex_type=complex_type,
                mult=mult,
                sample=10
            )

    elif complex is True:
        if name == 'nci_keyword':
            csearch(
                w_dir_main=csearch_methods_dir,
                program=program,
                smi=smi,
                name=name,
                crest_keywords=crest_keywords,
                auto_cluster=False # to keep track that crest options like --nci work
            )
        elif name == 'sample_keyword':
            csearch(
                w_dir_main=csearch_methods_dir,
                program=program,
                smi=smi,
                name=name,
                crest_keywords=crest_keywords,
                constraints_dist=constraints_dist,
                constraints_angle=constraints_angle,
                constraints_dihedral=constraints_dihedral,
                sample=17,
            )
        elif name == "cregen_clustering":
            csearch(
                w_dir_main=csearch_methods_dir,
                program=program,
                smi=smi,
                name=name,
                crest_keywords=crest_keywords,
                constraints_dist=constraints_dist,
                constraints_angle=constraints_angle,
                constraints_dihedral=constraints_dihedral,
                pytest_testing=True
            )
        else:
            csearch(
                w_dir_main=csearch_methods_dir,
                program=program,
                smi=smi,
                name=name,
                crest_keywords=crest_keywords,
                constraints_dist=constraints_dist,
                constraints_angle=constraints_angle,
                constraints_dihedral=constraints_dihedral,
            )

    if destination:
        file = str(csearch_methods_dir+"/Et_sdf_files/" + name + "_" + program + ".sdf")
    elif metal_complex is False or name in ['Ag_complex_crest','Cu_trigonal','Pd_metal_only','Ag_metal_crest']:
        file = str(csearch_methods_dir+"/CSEARCH/" + name + "_" + program + ".sdf")
    else:
        file = str(
            csearch_methods_dir
            + "/CSEARCH/"
            + name
            + "_1_"
            + program
            + ".sdf"
        )

    if program == 'crest':
        if name != 'ethane':
            file_crest = str(csearch_methods_dir+f"/CSEARCH/crest_xyz/{name}_crest.out")
            sdf_crest = str(csearch_methods_dir+f"/CSEARCH/{name}_crest.sdf")
        else:
            file_crest = str(csearch_methods_dir+f"/Et_sdf_files/crest_xyz/{name}_crest.out")
        outfile = open(file_crest, "r")
        outlines_crest = outfile.readlines()
        outfile.close()
    if name == 'rule_IrSP':
        # only the NAME_2_rdkit.sdf file passes the rule. For consistency, file_2 is file (since it exists)
        for suffix in ['A','B']:
            file_0 = str(csearch_methods_dir+"/CSEARCH/" + name + "_" + suffix + "_0_" + program + ".sdf")
            file_1 = str(csearch_methods_dir+"/CSEARCH/" + name + "_" + suffix + "_1_" + program + ".sdf")
            file = str(csearch_methods_dir+"/CSEARCH/" + name + "_" + suffix + "_2_" + program + ".sdf")
            assert not os.path.exists(file_0)
            assert not os.path.exists(file_1)
            assert os.path.exists(file)
    else:
        assert os.path.exists(file)
    mols = rdkit.Chem.SDMolSupplier(file, removeHs=False, sanitize=False)
    assert charge == int(mols[-1].GetProp("Real charge"))
    assert mult == int(mols[-1].GetProp("Mult"))

    # check that the metal is added back to the RDKit mol objects
    metal_found,iodine_found = False,False
    if name in ['Pd_complex','Pd_metal_only']:
        outfile = open(file, "r")
        outlines_sdf = outfile.readlines()
        outfile.close()
        for line in outlines_sdf:
            if ' Pd ' in line:
                metal_found = True
                break
        assert metal_found

    if name in ['Ag_complex_crest','Ag_metal_crest']:
        outfile = open(file, "r")
        outlines_sdf = outfile.readlines()
        outfile.close()
        for line in outlines_sdf:
            if ' Ag ' in line:
                metal_found = True
            if ' I ' in line:
                iodine_found = True

        assert metal_found
        assert not iodine_found

        # check that the metal atom is already added during the CREST calculation
        xyz_crest = str(csearch_methods_dir+f"/CSEARCH/crest_xyz/{name}_crest_xtb1.xyz")
        xyzfile = open(xyz_crest, "r")
        xyzlines_crest = xyzfile.readlines()
        xyzfile.close()
        xyz_metal,xyz_iodine = False,False
        for line in xyzlines_crest:
            if 'Ag  ' in line:
                xyz_metal = True
            if 'I  ' in line:
                xyz_iodine = True

        assert xyz_metal
        assert not xyz_iodine
    # if name == 'nci':
    #     assert len(mols) > 350
    #     # check that the conformers are sorted by their energy
    #     outfile_sdf = open(sdf_crest, "r")
    #     outlines_sdf = outfile_sdf.readlines()
    #     outfile_sdf.close()
    #     sdf_energies = []
    #     for i,line in enumerate(outlines_sdf):
    #         if line.startswith('>  <Energy>'):
    #             # check if the crest_keywords option works
    #             sdf_energies.append(float(outlines_sdf[i+1].split()[0]))
    #     assert len(sdf_energies) > 10
    #     assert sdf_energies == sorted(sdf_energies)

    # the n of conformers decreases when --nci is used
    elif name == 'nci_keyword':
        assert 100 < len(mols) < 250 # the number isn't exact, but it's between 100 and 250
        outfile = open(file_crest, "r")
        outlines_crest = outfile.readlines()
        outfile.close()
        for line in outlines_crest:
            if line.startswith(' > crest'):
                # check if the crest_keywords option works
                assert line.find('--nci --cbonds 0.5') > -1
                # check if there are no --cinp for constrained calcs
                assert line.find('--cinp') == -1
                # check if charge and mult are correct
                assert line.find('--chrg 0 --uhf 0') > -1
                break
        # check that the conformers are sorted by their number (avoid going from 1 to 10 instead of to 2)
        outfile_sdf = open(sdf_crest, "r")
        outlines_sdf = outfile_sdf.readlines()
        outfile_sdf.close()
        sdf_energies = []
        for i,line in enumerate(outlines_sdf):
            if line.startswith('>  <Energy>'):
                # check if the crest_keywords option works
                sdf_energies.append(float(outlines_sdf[i+1].split()[0]))
        assert len(sdf_energies) > 10
        assert sdf_energies == sorted(sdf_energies)

    elif name == 'cregen_clustering':
        assert len(mols) == 25
        dat_crest = str(csearch_methods_dir+f"/CSEARCH_data.dat")
        datfile = open(dat_crest, "r")
        datfile_crest = datfile.readlines()
        datfile.close()
        cregen_found = False
        for line in datfile_crest:
            if 'Starting CREGEN sorting' in line:
                cregen_found = True
        assert cregen_found
        # check that the conformers are sorted by their energ
        outfile_sdf = open(sdf_crest, "r")
        outlines_sdf = outfile_sdf.readlines()
        outfile_sdf.close()
        sdf_energies = []
        for i,line in enumerate(outlines_sdf):
            if line.startswith('>  <Energy>'):
                # check if the crest_keywords option works
                sdf_energies.append(float(outlines_sdf[i+1].split()[0]))
        assert len(sdf_energies) > 10
        assert sdf_energies == sorted(sdf_energies)

        # check if the first five most stable energies are included
        file = str("CSEARCH/" + name + "_" + program + "_all_confs.sdf")
        mols = rdkit.Chem.SDMolSupplier(file, removeHs=False)
        initial_five_E = []
        for mol in mols:
            initial_five_E.append(float(mol.GetProp("Energy")))
            if len(initial_five_E) == 5:
                break

        file = str("CSEARCH/" + name + "_" + program + ".sdf")
        mols = rdkit.Chem.SDMolSupplier(file, removeHs=False)
        final_five_E = []
        for mol in mols:
            final_five_E.append(float(mol.GetProp("Energy")))
            if len(final_five_E) == 5:
                break

        assert initial_five_E == final_five_E


    elif name == 'sample_keyword':
        assert len(mols) == 17
        dat_crest = str(csearch_methods_dir+f"/CSEARCH_data.dat")
        datfile = open(dat_crest, "r")
        datfile_crest = datfile.readlines()
        datfile.close()
        butina_found, cregen_clust_found = False,False
        for line in datfile_crest:
            if 'conformers using a combination of energies and Butina RMS-based clustering' in line:
                butina_found = True
            if 'Starting CREGEN sorting' in line:
                cregen_clust_found = True
        assert butina_found
        assert cregen_clust_found
        # check that the conformers are sorted by their number (avoid going from 1 to 10 instead of to 2)
        outfile_sdf = open(sdf_crest, "r")
        outlines_sdf = outfile_sdf.readlines()
        outfile_sdf.close()
        sdf_energies = []
        for i,line in enumerate(outlines_sdf):
            if line.startswith('>  <Energy>'):
                # check if the crest_keywords option works
                sdf_energies.append(float(outlines_sdf[i+1].split()[0]))
        assert len(sdf_energies) == 17
        assert sdf_energies == sorted(sdf_energies)
    
    elif name == 'ethane':
        assert len(mols) >= 1
    elif name == 'ts':
        assert len(mols) == 1
        constrain_line_found = False
        for i,line in enumerate(outlines_crest):
            if line.startswith(' > crest'):
                # check if the the constraints are active in the CREST opt
                assert line.find('--cinp .xcontrol.sample') > -1
                # check if charge and mult are correct
                assert line.find('--chrg -1 --uhf 0') > -1
            if line.startswith('> $constrain'):
                constrain_line_found = True
                assert outlines_crest[i+1].find('distance: 2,3,1.8') > -1
                assert outlines_crest[i+2].find('distance: 3,1,1.8') > -1
                assert outlines_crest[i+3].find('angle: 2,3,1,180.0') > -1
                assert outlines_crest[i+4].find('force constant=0.5') > -1
                assert outlines_crest[i+5].find('reference=coord.ref') > -1
                assert outlines_crest[i+6].find('$metadyn') > -1
                assert outlines_crest[i+7].find('atoms: 4-6') > -1
            if line.find('Generating MTD length') > -1:
                break
        assert constrain_line_found

        file_xtb1 = str(csearch_methods_dir+"/CSEARCH/crest_xyz/ts_crest_xtb1.out")
        outfile_xtb1 = open(file_xtb1, "r")
        outlines_xtb1 = outfile_xtb1.readlines()
        outfile.close()
        constrain_found_xtb = False
        for i,line in enumerate(outlines_xtb1):
            if line.find('program call') > -1:
                # check if charge and mult are correct in xTB
                assert line.find('-c -1 --uhf 0') > -1
            elif line.find('constraining bond 2 3 to    1.8') > -1:
                constrain_found_xtb = True
                # 5 bond constraints should be set in xtb2
                assert outlines_xtb1[i+6].find('molecular fragmentation') > -1
            elif line.find('G F N 2 - x T B') > -1:
                break
        assert constrain_found_xtb
                
        file_xtb2 = str(csearch_methods_dir+"/CSEARCH/crest_xyz/ts_crest_xtb2.out")
        outfile_xtb2 = open(file_xtb2, "r")
        outlines_xtb2 = outfile_xtb2.readlines()
        outfile.close()
        constrain_found_xtb = False
        for i,line in enumerate(outlines_xtb2):
            if line.find('program call') > -1:
                # check if charge and mult are correct in xTB
                assert line.find('-c -1 --uhf 0') > -1
            elif line.find('constraining bond 2 3 to    1.8') > -1:
                constrain_found_xtb = True
                # only 3 constraints should be set in xtb2
                assert outlines_xtb2[i+4].find('molecular fragmentation') > -1
            elif line.find('G F N 2 - x T B') > -1:
                break
        assert constrain_found_xtb

    elif name != 'rule_IrSP':
        assert len(mols) == output_nummols
    os.chdir(w_dir_main)


# tests for removing foler
@pytest.mark.parametrize(
    "folder_list, file_list",
    [
        # tests for conformer generation with RDKit
        (
            ["tests/csearch_methods/CSEARCH","tests/csearch_rdkit_summ/CSEARCH","tests/csearch_fullmonte/CSEARCH","tests/csearch_crest/CSEARCH","tests/csearch_input/CSEARCH","tests/csearch_varfile/CSEARCH"],
            ["tests/csearch_methods/CSEARCH*","tests/csearch_rdkit_summ/CSEARCH*","tests/csearch_fullmonte/CSEARCH*","tests/csearch_crest/CSEARCH*","tests/csearch_input/CSEARCH*","tests/csearch_varfile/CSEARCH*"]
        ),
    ],
)
def test_remove(folder_list, file_list):
    for i,folder in enumerate(folder_list):
        shutil.rmtree(w_dir_main + "/" + folder)
        for f in glob.glob(file_list[i]):
            os.remove(f)
    os.chdir(w_dir_main)
