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
    os.chdir(csearch_input_dir)
    # runs the program with the different tests
    csearch(w_dir_main=csearch_input_dir, program=program, input=input)

    # tests here
    if input in ["pentane.smi", "pentane.csv"]:
        file1 = str(
            "CSEARCH/"
            + "butane_"
            + input.split(".")[1]
            + "_"
            + program
            + ".sdf"
        )
        file2 = str(
            "CSEARCH/"
            + "pentane_"
            + input.split(".")[1]
            + "_"
            + program
            + ".sdf"
        )
        mol1 = rdkit.Chem.SDMolSupplier(file1, removeHs=False)
        mol2 = rdkit.Chem.SDMolSupplier(file2, removeHs=False)
        assert len(mol1) == output_nummols[0]
        assert len(mol2) == output_nummols[1]
    elif input in ["molecules.cdx"]:
        file1 = str("CSEARCH/" + "molecules_0_" + program + ".sdf")
        file2 = str("CSEARCH/" + "molecules_1_" + program + ".sdf")
        mol1 = rdkit.Chem.SDMolSupplier(file1, removeHs=False)
        mol2 = rdkit.Chem.SDMolSupplier(file2, removeHs=False)
        assert len(mol1) == output_nummols[0]
        assert len(mol2) == output_nummols[1]
    elif input in ["pentane_sdf.sdf"]:
        file = str(
            "CSEARCH/" + input.split(".")[0] + "_" + program + ".sdf"
        )
        mols = rdkit.Chem.SDMolSupplier(file, removeHs=False)
        assert len(mols) == output_nummols
    else:
        file = str(
            "CSEARCH/" + input.split(".")[0] + "_" + program + ".sdf"
        )
        mols = rdkit.Chem.SDMolSupplier(file, removeHs=False)
        assert len(mols) == output_nummols
    os.chdir(w_dir_main)


# # tests for parameters of csearch random initialzation
# @pytest.mark.parametrize(
#     "program, smi, name, max_matches_rmsd , max_mol_wt , ff, degree, output, max_torsions, prefix, output_nummols ",
#     [
#         # tests for conformer generation with RDKit
#         ("summ", "CCCCC", "pentane", 500, 200, "MMFF", 30, ".sdf", 20, "mol", 4),
#     ],
# )
# def test_csearch_others_parameters(
#     program,
#     smi,
#     name,
#     max_matches_rmsd,
#     max_mol_wt,
#     ff,
#     degree,
#     output,
#     max_torsions,
#     prefix,
#     output_nummols,
# ):
#     os.chdir(csearch_others_dir)
#     # runs the program with the different tests
#     csearch(
#         w_dir_main=csearch_others_dir,
#         program=program,
#         smi=smi,
#         name=name,
#         max_matches_rmsd=max_matches_rmsd,
#         max_mol_wt=max_mol_wt,
#         ff=ff,
#         degree=degree,
#         output=output,
#         max_torsions=max_torsions,
#         prefix=prefix,
#     )

#     # tests here
#     file = str(
#         "CSEARCH/" + prefix + "_" + name + "_" + program + ".sdf"
#     )
#     mols = rdkit.Chem.SDMolSupplier(file, removeHs=False)
#     assert len(mols) == output_nummols
#     os.chdir(w_dir_main)


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

    # tests here
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
    assert len(mols) == output_nummols
    os.chdir(w_dir_main)

    file_crest = str(csearch_crest_dir+f"/CSEARCH/crest_xyz/{name}.out")
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
        # tests for conformer generation with RDKit
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
        output_nummols=output_nummols,
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
        # (
        #     "summ",
        #     "CCCCC",
        #     "pentane_summ",
        #     0,
        #     1,
        #     "auto",
        #     100,
        #     True,
        #     40,
        #     0.0001,
        #     0.2,
        #     0.1,
        #     5,
        # ),
        (
            "rdkit",
            "CC[CH]CC",
            "radical_rdkit",
            0,
            2,
            "auto",
            10,
            True,
            10,
            0.000001,
            0.6,
            0.3,
            41,
        ),
        # (
        #     "summ",
        #     "CC[CH]CC",
        #     "radical_summ",
        #     0,
        #     2,
        #     "auto",
        #     500,
        #     True,
        #     2,
        #     0.0001,
        #     0.2,
        #     0.2,
        #     3,
        # ),
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
            5,
        ),
        # (
        #     "summ",
        #     "C[NH2+]CC",
        #     "charged_summ",
        #     1,
        #     0,
        #     "auto",
        #     1000,
        #     False,
        #     10,
        #     0.0001,
        #     0.2,
        #     0.2,
        #     2,
        # ),
    ],
)
def test_csearch_rdkit_summ_parameters(
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
        output_nummols=output_nummols,
    )

    # tests here
    file = str("CSEARCH/" + name + "_" + program + ".sdf")
    mols = rdkit.Chem.SDMolSupplier(file, removeHs=False)
    assert len(mols) == output_nummols
    assert charge == int(mols[0].GetProp("Real charge"))
    assert mult == int(mols[0].GetProp("Mult"))
    os.chdir(w_dir_main)


# tests for individual organic molecules and metal complexes with different types of csearch methods
@pytest.mark.parametrize(
    "program, smi, name, complex, metal_complex, metal, metal_oxi, complex_type, constraints_dist, constraints_angle, constraints_dihedral, charge, mult, crest_keywords, destination, output_nummols",
    [
        # tests for conformer generation with RDKit, SUMM, FullMonte and CREST
        (
            "rdkit",
            "CCCCC",
            "pentane_RD",
            False,
            False,
            None,
            None,
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
        #     None,
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
            None,
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
        (
            "rdkit",
            "I[Pd]([PH3+])(F)Cl",
            "Pd_complex",
            False,
            True,
            ["Pd"],
            [2],
            "squareplanar",
            None,
            None,
            None,
            -1,
            1,
            None,
            False,
            1
        ),
        # compatibility of CREST with metal complexes
        (
            "crest",
            "[NH3+][Ag][NH3+]",
            "Ag_complex_crest",
            False,
            True,
            ["Ag"],
            [1],
            "linear",
            None,
            None,
            None,
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
            None,
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
        (
            "crest",
            "C.O",
            "nci",
            True,
            False,
            None,
            None,
            None,
            [],
            [],
            [],
            0,
            1,
            None,
            False,
            None,
        ),
        (
            "crest",
            "C.O",
            "nci_keyword",
            True,
            False,
            None,
            None,
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
            None,
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
    metal,
    metal_oxi,
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
        csearch(
            w_dir_main=csearch_methods_dir,
            program=program,
            smi=smi,
            name=name,
            metal_atoms=metal,
            metal_oxi=metal_oxi,
            complex_type=complex_type,
        )

    elif complex is True:
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
    elif metal_complex is False or name == 'Ag_complex_crest':
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
        file_crest = str(csearch_methods_dir+f"/CSEARCH/crest_xyz/{name}.out")
        outfile = open(file_crest, "r")
        outlines_crest = outfile.readlines()
        outfile.close()
    assert os.path.exists(file)
    mols = rdkit.Chem.SDMolSupplier(file, removeHs=False, sanitize=False)
    assert charge == int(mols[-1].GetProp("Real charge"))
    assert mult == int(mols[-1].GetProp("Mult"))
    if name == 'nci':
        assert len(mols) > 350
    # the n of conformers decreases when --nci is used
    elif name == 'nci_keyword':
        assert len(mols) < 300
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
                assert outlines_crest[i+7].find('atoms: 4,5,6') > -1
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

    else:
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
