=====================
Descriptor Generation
=====================

.. tabs::

   .. tab:: 1. Organics

      .. tabs::

         .. tab:: 1.1. Basic

            **Input CSV:**

            +---------+-----------+
            | SMILES  | code_name |
            +=========+===========+
            | C       | mol_1     |
            +---------+-----------+
            | CC      | mol_2     |
            +---------+-----------+
            | CCC     | mol_3     |
            +---------+-----------+
            | CCCC    | mol_4     |
            +---------+-----------+

            **Python (Jupyter Notebook):**

            .. code:: python

               from aqme.qdescp import qdescp

               qdescp(
                   input="FILENAME.csv"
               )

            **Outputs:**

            - CSV files with descriptors (standard: ``AQME-ROBERT_interpret_FILENAME.csv``)
            - No common substructure → only molecular descriptors generated


         .. tab:: 1.2. Auto detection of common substructure

            **Input CSV:**

            +---------+-----------+
            | SMILES  | code_name |
            +=========+===========+
            | P       | mol_1     |
            +---------+-----------+
            | PC      | mol_2     |
            +---------+-----------+
            | P(C)C   | mol_3     |
            +---------+-----------+
            | P(C)(C)C| mol_4     |
            +---------+-----------+

            **Python (Jupyter Notebook):**

            .. code:: python

               from aqme.qdescp import qdescp

               qdescp(
                   input="FILENAME.csv"
               )

            **Outputs:**

            - CSV files with descriptors
            - Common substructure detected (P atom)
            - Molecular and atomic descriptors generated


         .. tab:: 1.3. Custom SMARTS for atomic descriptors

            **Input CSV:**

            +-----------+-----------+
            | SMILES    | code_name |
            +===========+===========+
            | C(O)N     | mol_1     |
            +-----------+-----------+
            | CC(O)N    | mol_2     |
            +-----------+-----------+
            | CCC(O)N   | mol_3     |
            +-----------+-----------+
            | CCCC(O)N  | mol_4     |
            +-----------+-----------+

            **Python (Jupyter Notebook):**

            .. code:: python

               from aqme.qdescp import qdescp

               qdescp(
                   input="FILENAME.csv",
                   qdescp_atoms=["O"]
               )

            **Outputs:**

            - Only O atoms used for atomic descriptors
            - Molecular descriptors always included

            .. note::

               Other SMARTS patterns can be used (e.g., ``C=O``, ``C#N``). The pattern must appear once in all molecules.


         .. tab:: 1.4. Atom indices for atomic descriptors

            **Input CSV:**

            +--------------------------------+-----------+
            | SMILES                         | code_name |
            +================================+===========+
            | [P:1]([H])([H])([H])           | mol_1     |
            +--------------------------------+-----------+
            | [P:1]([H])([H])C               | mol_2     |
            +--------------------------------+-----------+
            | [P:1]([H])(C)C                 | mol_3     |
            +--------------------------------+-----------+
            | [P:1](C)(C)C                   | mol_4     |
            +--------------------------------+-----------+

            **Python (Jupyter Notebook):**

            .. code:: python

               from aqme.qdescp import qdescp

               qdescp(
                   input="FILENAME.csv",
                   qdescp_atoms=["1"]
               )

            **Outputs:**

            - Descriptors only for atoms labeled with index ``1``
            - Molecular descriptors always included


         .. tab:: 1.5. Combining sets of two molecules

            **Input CSV:**

            +-----------+------------------+-----------------+
            | code_name | SMILES_phosph    | SMILES_alk      |
            +===========+==================+=================+
            | mol_1     | P                | C               |
            +-----------+------------------+-----------------+
            | mol_2     | PC               | CC              |
            +-----------+------------------+-----------------+
            | mol_3     | P(C)C            | CCC             |
            +-----------+------------------+-----------------+
            | mol_4     | P(C)(C)C         | CCCC            |
            +-----------+------------------+-----------------+

            **Python (Jupyter Notebook):**

            .. code:: python

               from aqme.qdescp import qdescp

               qdescp(
                   input="FILENAME.csv"
               )

            **Outputs:**

            - Descriptors generated independently for each component and combined into a single descriptor array


   .. tab:: 2. Metal complexes

      .. tabs::

         .. tab:: 2.1. Basic

            **Input CSV:**

            +-----------+------------------------------------------------------+--------+------+
            | code_name | SMILES                                               | charge | mult |
            +===========+======================================================+========+======+
            | mol_1     | [H][N+]1([H])[Cu][N+]([H])([H])CC1                   | 2      | 2    |
            +-----------+------------------------------------------------------+--------+------+
            | mol_2     | [H][N+]1([H])[Cu]N([H])CC1                           | 1      | 2    |
            +-----------+------------------------------------------------------+--------+------+
            | mol_3     | [H]N1[Cu]N([H])CC1                                   | 0      | 2    |
            +-----------+------------------------------------------------------+--------+------+
            | mol_4     | C[N+]1(C)[Cu][N+](C)(C)C=C1                          | 2      | 2    |
            +-----------+------------------------------------------------------+--------+------+

            **Python (Jupyter Notebook):**

            .. code:: python

               from aqme.qdescp import qdescp

               qdescp(
                   input="FILENAME.csv"
               )

            **Outputs:**

            - Atomic descriptors for Cu automatically generated
            - Molecular descriptors always included
            - Charges and multiplicity required


         .. tab:: 2.2. Forcing squareplanar geometry

            **Input CSV:**

            +-----------+------------------------------------------------+--------+------+----------------+--------------------------------------+
            | code_name | SMILES                                         | charge | mult | complex_type   | geom                                 |
            +===========+================================================+========+======+================+======================================+
            | mol_1     | Cl[Pd]([PH3+])(Cl)[NH3+]                       | 0      | 1    | squareplanar   | "['[Cl][Pd][Cl]',180]"               |
            +-----------+------------------------------------------------+--------+------+----------------+--------------------------------------+
            | mol_2     | Cl[Pd](Cl)([P+](C)(C)C)[NH3+]                  | 0      | 1    | squareplanar   | "['[Cl][Pd][Cl]',180]"               |
            +-----------+------------------------------------------------+--------+------+----------------+--------------------------------------+
            | mol_3     | Cl[Pd]([PH3+])(Cl)[N+](C)(C)C                  | 0      | 1    | squareplanar   | "['[Cl][Pd][Cl]',180]"               |
            +-----------+------------------------------------------------+--------+------+----------------+--------------------------------------+
            | mol_4     | Cl[Pd](Cl)([P+](C)(C)C)[N+](C)(C)C             | 0      | 1    | squareplanar   | "['[Cl][Pd][Cl]',180]"               |
            +-----------+------------------------------------------------+--------+------+----------------+--------------------------------------+

            **Python (Jupyter Notebook):**

            .. code:: python

               from aqme.qdescp import qdescp

               qdescp(
                   input="FILENAME.csv"
               )

            **Outputs:**

            - Atomic descriptors for Pd, N and P (detected automatically)
            - Molecular descriptors always included
            - Charges and multiplicity required
            - Squareplanar geometry enforced via the ``complex_type`` column
            - Two chloride ligands in trans (Cl-Pd-Cl angle at 180 degrees) enforced via the ``geom`` column


         .. tab:: 2.3. Custom SMARTS for atomic descriptors

            **Input CSV:**

            +-----------+---------------------------------------------+--------+------+
            | code_name | SMILES                                      | charge | mult |
            +===========+=============================================+========+======+
            | mol_1     | O=CC1=CN(C)[Cu]N1C                          | 0      | 2    |
            +-----------+---------------------------------------------+--------+------+
            | mol_2     | O=CC1=C(C)N(C)[Cu]N1C                       | 0      | 2    |
            +-----------+---------------------------------------------+--------+------+
            | mol_3     | O=CC1=C(N([Cu]N1C(C)(C)C)C)C                | 0      | 2    |
            +-----------+---------------------------------------------+--------+------+
            | mol_4     | O=C(C)C1=C(C)N([Cu]N1C)C                    | 0      | 2    |
            +-----------+---------------------------------------------+--------+------+

            **Python (Jupyter Notebook):**

            .. code:: python

               from aqme.qdescp import qdescp

               qdescp(
                   input="FILENAME.csv",
                   qdescp_atoms=["C=O"]
               )

            **Outputs:**

            - Atomic descriptors only for C and O atoms from carbonyl groups
            - Molecular descriptors always included


         .. tab:: 2.4. Atom indices for atomic descriptors

            **Input CSV:**

            +-----------+--------------------------------------------------------------+--------+------+
            | code_name | SMILES                                                       | charge | mult |
            +===========+==============================================================+========+======+
            | mol_1     | [H][N+]1([Cu][N:1](C=C1)[H])[H]                              | 1      | 2    |
            +-----------+--------------------------------------------------------------+--------+------+
            | mol_2     | [H][N+]1([H])[Cu][N:1]([H])CC1                               | 1      | 2    |
            +-----------+--------------------------------------------------------------+--------+------+
            | mol_3     | [H][N+]1([Cu][N:1](C(C)=C1)[H])[H]                           | 1      | 2    |
            +-----------+--------------------------------------------------------------+--------+------+
            | mol_4     | [H][N+]1([Cu][N:1](C(C)=C1C)[H])[H]                          | 1      | 2    |
            +-----------+--------------------------------------------------------------+--------+------+

            **Python (Jupyter Notebook):**

            .. code:: python

               from aqme.qdescp import qdescp

               qdescp(
                   input="FILENAME.csv",
                   qdescp_atoms=["1"]
               )

            **Outputs:**

            - Atomic descriptors only for indexed atoms
            - Molecular descriptors always included


   .. tab:: 3. From SDF/XYZ

      **Python (Jupyter Notebook):**

      .. code:: python

         from aqme.qdescp import qdescp

         qdescp(
             files="*.sdf"
         )

      **Outputs:**

      - CSV descriptor files generated from 3D structures

      .. note::

         Useful when SMILES are not suitable (e.g., axial chirality or noncovalent complexes).