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
    file = str("CSEARCH/" + "rdkit/" + "pentane_varfile" + "_rdkit" + ".sdf")
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
            + program
            + "/"
            + "butane_"
            + input.split(".")[1]
            + "_"
            + program
            + ".sdf"
        )
        file2 = str(
            "CSEARCH/"
            + program
            + "/"
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
        file1 = str("CSEARCH/" + program + "/" + "molecules_0_" + program + ".sdf")
        file2 = str("CSEARCH/" + program + "/" + "molecules_1_" + program + ".sdf")
        mol1 = rdkit.Chem.SDMolSupplier(file1, removeHs=False)
        mol2 = rdkit.Chem.SDMolSupplier(file2, removeHs=False)
        assert len(mol1) == output_nummols[0]
        assert len(mol2) == output_nummols[1]
    elif input in ["pentane_sdf.sdf"]:
        file = str(
            "CSEARCH/" + program + "/" + input.split(".")[0] + "_" + program + ".sdf"
        )
        mols = rdkit.Chem.SDMolSupplier(file, removeHs=False)
        assert len(mols) == output_nummols
    else:
        file = str(
            "CSEARCH/" + program + "/" + input.split(".")[0] + "_" + program + ".sdf"
        )
        mols = rdkit.Chem.SDMolSupplier(file, removeHs=False)
        assert len(mols) == output_nummols
    os.chdir(w_dir_main)


# tests for parameters of csearch random initialzation
@pytest.mark.parametrize(
    "program, smi, name, max_matches_rmsd , max_mol_wt , ff, degree, verbose, output, max_torsions, prefix, output_nummols ",
    [
        # tests for conformer generation with RDKit
        ("summ", "CCCCC", "pentane", 500, 200, "MMFF", 30, True, ".sdf", 20, "mol", 4),
    ],
)
def test_csearch_others_parameters(
    program,
    smi,
    name,
    max_matches_rmsd,
    max_mol_wt,
    ff,
    degree,
    verbose,
    output,
    max_torsions,
    prefix,
    output_nummols,
):
    os.chdir(csearch_others_dir)
    # runs the program with the different tests
    csearch(
        w_dir_main=csearch_others_dir,
        program=program,
        smi=smi,
        name=name,
        max_matches_rmsd=max_matches_rmsd,
        max_mol_wt=max_mol_wt,
        ff=ff,
        degree=degree,
        verbose=verbose,
        output=output,
        max_torsions=max_torsions,
        prefix=prefix,
    )

    # tests here
    file = str(
        "CSEARCH/" + program + "/" + prefix + "_" + name + "_" + program + ".sdf"
    )
    mols = rdkit.Chem.SDMolSupplier(file, removeHs=False)
    assert len(mols) == output_nummols
    os.chdir(w_dir_main)


# tests for parameters of csearch
@pytest.mark.parametrize(
    "program, smi, name, cregen, cregen_keywords, crest_keywords, output_nummols",
    [
        # tests for conformer generation with RDKit
        (
            "crest",
            "C",
            "methane",
            True,
            "--ethr 1 --rthr 0.5 --bthr 0.7 --ewin 4 --cbonds 0.8 --cluster",
            "--alpb benzene",
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
    output_nummols,
):
    os.chdir(csearch_crest_dir)
    # runs the program with the different tests
    csearch(
        w_dir_main=csearch_crest_dir,
        program=program,
        smi=smi,
        name=name,
        cregen=cregen,
        cregen_keywords=cregen_keywords,
    )

    # tests here
    file_cluster = str(
        csearch_crest_dir
        + "/CSEARCH/"
        + program
        + "_xyz/"
        + name
        + "_"
        + program
        + "/crest_clustered.xyz"
    )
    assert os.path.exists(file_cluster)
    file = str("CSEARCH/" + program + "/" + name + "_" + program + ".sdf")
    mols = rdkit.Chem.SDMolSupplier(file, removeHs=False)
    assert len(mols) == output_nummols
    os.chdir(w_dir_main)


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
    file = str("CSEARCH/" + program + "/" + name + "_" + program + ".sdf")
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
            "summ",
            "CCCCC",
            "pentane_summ",
            0,
            1,
            "auto",
            100,
            True,
            40,
            0.0001,
            0.2,
            0.1,
            5,
        ),
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
        (
            "summ",
            "CC[CH]CC",
            "radical_summ",
            0,
            2,
            "auto",
            500,
            True,
            2,
            0.0001,
            0.2,
            0.2,
            3,
        ),
        (
            "rdkit",
            "C[NH2+]CC",
            "charged_rdkit",
            1,
            0,
            "auto",
            10,
            True,
            10,
            0.0001,
            4,
            0.6,
            5,
        ),
        (
            "summ",
            "C[NH2+]CC",
            "charged_summ",
            1,
            0,
            "auto",
            1000,
            False,
            10,
            0.0001,
            0.2,
            0.2,
            2,
        ),
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
    file = str("CSEARCH/" + program + "/" + name + "_" + program + ".sdf")
    mols = rdkit.Chem.SDMolSupplier(file, removeHs=False)
    assert len(mols) == output_nummols
    assert charge == int(mols[0].GetProp("Real charge"))
    assert mult == int(mols[0].GetProp("Mult"))
    os.chdir(w_dir_main)


# tests for individual organic molecules and metal complexes with different types of csearch methods
@pytest.mark.parametrize(
    "program, smi, name, complex, metal_complex, metal, metal_oxi, complex_type, constraints_dist, constraints_angle, constraints_dihedral, cregen, output_nummols",
    [
        # tests for conformer generation with RDKit, xTB and ANI1
        (
            "rdkit",
            "CCCCC",
            "pentane",
            False,
            False,
            None,
            None,
            None,
            [],
            [],
            [],
            False,
            4,
        ),
        (
            "summ",
            "CCCCC",
            "pentane",
            False,
            False,
            None,
            None,
            None,
            [],
            [],
            [],
            False,
            4,
        ),
        (
            "fullmonte",
            "CCCCC",
            "pentane",
            False,
            False,
            None,
            None,
            None,
            [],
            [],
            [],
            False,
            4,
        ),
        (
            "crest",
            "C",
            "methane",
            False,
            False,
            None,
            None,
            None,
            [],
            [],
            [],
            True,
            1,
        ),
        # (
        #     "rdkit",
        #     "I[Pd](Cl)([PH3+])[N+]1=CC=CC=C1",
        #     "Pd_complex",
        #     False,
        #     True,
        #     ["Pd"],
        #     [2],
        #     "squareplanar",
        #     None,
        #     None,
        #     None,
        # False,
        # None
        # ),
        # (
        #     "summ",
        #     "I[Pd](Cl)([PH3+])[N+]1=CC=CC=C1",
        #     "Pd_complex",
        #     False,
        #     True,
        #     ["Pd"],
        #     [2],
        #     "squareplanar",
        #     [],
        #     [],
        #     [],
        # False,
        # None
        # ),
        # (
        #     "fullmonte",
        #     "I[Pd](Cl)([PH3+])[N+]1=CC=CC=C1",
        #     "Pd_complex",
        #     False,
        #     True,
        #     ["Pd"],
        #     [2],
        #     "squareplanar",
        #     [],
        #     [],
        #     [],
        # False
        # None
        # ),
        # (
        #     "crest",
        #     "I[Pd](Cl)([PH3+])[N+]1=CC=CC=C1",
        #     "Pd_complex",
        #     False,
        #     True,
        #     ["Pd"],
        #     [2],
        #     "squareplanar",
        #     [],
        #     [],
        #     [],
        # True
        # None
        # ),
        # (
        #     "crest",
        #     "C.O",
        #     "nci",
        #     True,
        #     False,
        #     None,
        #     None,
        #     None,
        #     [],
        #     [],
        #     [],
        #     True,
        #     436,
        # ),  # CHECK THIS TEST again. working now
        # (
        #     "crest",
        #     "[Cl-:9].[F:4][C:5]([C:6]([H:12])([H:13])[H:14])([C:7]([H:15])([H:16])[H:17])[C:8]([H:18])([H:19])[H:20].[O:3]([H:10])[H:11]",
        #     "ts",
        #     True,
        #     False,
        #     None,
        #     None,
        #     None,
        #     [[4, 5, 1.8], [5, 9, 1.8]],
        #     [[4, 5, 9, 180]],
        #     [],
        #     True,
        #     1,
        # ),
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
    cregen,
    output_nummols,
):
    os.chdir(csearch_methods_dir)
    # runs the program with the different tests
    if not complex and not metal_complex:
        csearch(
            w_dir_main=csearch_methods_dir,
            program=program,
            smi=smi,
            name=name,
            cregen=cregen,
        )

    if metal_complex is True:
        csearch(
            w_dir_main=csearch_methods_dir,
            program=program,
            smi=smi,
            name=name,
            metal_complex=metal_complex,
            metal=metal,
            metal_oxi=metal_oxi,
            complex_type=complex_type,
            cregen=cregen,
        )

    if complex is True:
        csearch(
            w_dir_main=csearch_methods_dir,
            program=program,
            smi=smi,
            name=name,
            complex=complex,
            cregen=cregen,
            constraints_dist=constraints_dist,
            constraints_angle=constraints_angle,
            constraints_dihedral=constraints_dihedral,
        )

    if metal_complex is False:
        file = str("CSEARCH/" + program + "/" + name + "_" + program + ".sdf")
    else:
        file = str(
            csearch_cmin_dir
            + "/CSEARCH/"
            + program
            + "/"
            + name
            + "_1_"
            + program
            + ".sdf"
        )

    assert os.path.exists(file)
    mols = rdkit.Chem.SDMolSupplier(file, removeHs=False)
    assert len(mols) == output_nummols
    os.chdir(w_dir_main)


# tests for removing foler
@pytest.mark.parametrize(
    "remove, folder, file",
    [
        # tests for conformer generation with RDKit
        (
            True,
            "tests/csearch_methods/CSEARCH",
            "tests/csearch_methods/CSEARCH*",
        ),
        (
            True,
            "tests/csearch_rdkit_summ/CSEARCH",
            "tests/csearch_rdkit_summ/CSEARCH*",
        ),
        (
            True,
            "tests/csearch_fullmonte/CSEARCH",
            "tests/csearch_fullmonte/CSEARCH*",
        ),
        (True, "tests/csearch_crest/CSEARCH", "tests/csearch_crest/CSEARCH*"),
        (True, "tests/csearch_others/CSEARCH", "tests/csearch_others/CSEARCH*"),
        (True, "tests/csearch_input/CSEARCH", "tests/csearch_input/CSEARCH*"),
        (
            True,
            "tests/csearch_varfile/CSEARCH",
            "tests/csearch_varfile/CSEARCH*",
        ),
    ],
)
def test_remove(remove, folder, file):
    # os.chdir(w_dir_main)
    shutil.rmtree(w_dir_main + "/" + folder)
    for f in glob.glob(file):
        os.remove(f)
