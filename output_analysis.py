import sys, os, glob, shutil

'''
This part reads the files and retrieve some info. Also, if the calculations
finished with imaginary frequencies or Error termination, it will generate
new input files to resubmit. When the calculations end with imag freq,
the new input file is a displaced geometry obtained from the imag freq. When
the calcs end with Error termination, the new inout file corresponds to the
point with minimum gradient in the optimization that failed.

How this works:
- Place all the output files to analyze in folder XXX
'''

def output_analyzer(log_files, w_dir, lot, bs, chk, nprocs, mem, input):

    # Atom IDs
    periodictable = ["","H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr",
        "Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl",
        "Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Uub","Uut","Uuq","Uup","Uuh","Uus","Uuo"]

    #made it global for all functions
    rms = 10000
    #defined the variable stop_rms, standor
    stop_rms = 0
    standor = 0
    NATOMS =0

    for file in log_files:
        outfile = open(file,"r")
        outlines = outfile.readlines()
        ATOMTYPES, CARTESIANS = [],[]
        FREQS, REDMASS, FORCECONST, NORMALMODE = [],[],[],[]; IM_FREQS = 0
        freqs_so_far = 0
        TERMINATION = "unfinished"
        for i in range(0,len(outlines)):
            # Get the name of the compound (specified in the title)
            if outlines[i].find('Symbolic Z-matrix:') > -1:
                name = outlines[i-2].split()[0]
            # Determine the kind of job termination
            if outlines[i].find("Normal termination") > -1:
                TERMINATION = "normal"
            elif outlines[i].find("Error termination") > -1:
                TERMINATION = "error"
            # Determine charge and multiplicity
            if outlines[i].find("Charge = ") > -1:
                CHARGE = int(outlines[i].split()[2])
                MULT = int(outlines[i].split()[5].rstrip("\n"))

        for i in range(0,len(outlines)):
            if TERMINATION == "normal":
                # Sets where the final coordinates are inside the file
                if outlines[i].find("Input orientation") > -1: standor = i
                if outlines[i].find("Standard orientation") > -1: standor = i
                if outlines[i].find("Distance matrix") > -1 or outlines[i].find("Rotational constants") >-1:
                    if outlines[i-1].find("-------") > -1:
                        NATOMS = i-standor-6
                # Get the frequencies and identifies negative frequencies
                if outlines[i].find(" Frequencies -- ") > -1:
                    nfreqs = len(outlines[i].split())
                    for j in range(2, nfreqs):
                        FREQS.append(float(outlines[i].split()[j]))
                        NORMALMODE.append([])
                        if float(outlines[i].split()[j]) < 0.0: IM_FREQS += 1
                    for j in range(3, nfreqs+1): REDMASS.append(float(outlines[i+1].split()[j]))
                    for j in range(3, nfreqs+1): FORCECONST.append(float(outlines[i+2].split()[j]))
                    for j in range(0,NATOMS):
                        for k in range(0, nfreqs-2):
                            NORMALMODE[(freqs_so_far + k)].append([float(outlines[i+5+j].split()[3*k+2]), float(outlines[i+5+j].split()[3*k+3]), float(outlines[i+5+j].split()[3*k+4])])
                    freqs_so_far = freqs_so_far + nfreqs - 2
            if TERMINATION != "normal":
                if outlines[i].find('Cartesian Forces:  Max') > -1:
                    if float(outlines[i].split()[5]) < rms:
                        rms = float(outlines[i].split()[5])
                        stop_rms = i

        if TERMINATION == "normal":
            # Get the coordinates for jobs that finished well with and without imag. freqs
            try: standor
            except NameError: pass
            else:
                for i in range (standor+5,standor+5+NATOMS):
                    massno = int(outlines[i].split()[1])
                    if massno < len(periodictable):
                        atom_symbol = periodictable[massno]
                    else: atom_symbol = "XX"
                    ATOMTYPES.append(atom_symbol)
                    CARTESIANS.append([float(outlines[i].split()[3]), float(outlines[i].split()[4]), float(outlines[i].split()[5])])

        if TERMINATION != "normal":
            # Get the coordinates for jobs that did not finished or finished with an error
            for i in range(0,stop_rms):
                # Sets where the final coordinates are inside the file
                if outlines[i].find("Input orientation") > -1: standor = i
                if outlines[i].find("Standard orientation") > -1: standor = i
                if outlines[i].find("Distance matrix") > -1 or outlines[i].find("Rotational constants") >-1:
                    if outlines[i-1].find("-------") > -1:
                        NATOMS = i-standor-6
            for i in range (standor+5,standor+5+NATOMS):
                massno = int(outlines[i].split()[1])
                if massno < len(periodictable):
                    atom_symbol = periodictable[massno]
                else: atom_symbol = "XX"
                ATOMTYPES.append(atom_symbol)
                CARTESIANS.append([float(outlines[i].split()[3]), float(outlines[i].split()[4]), float(outlines[i].split()[5])])

        # This part fixes jobs with imaginary freqs
        if IM_FREQS > 0:
            # Multiplies the imaginary normal mode vector by this amount (from -1 to 1).
            amplitude = 0.2 # default in pyQRC
            shift = []

            # Save the original Cartesian coordinates before they are altered
            orig_carts = []
            for atom in range(0,NATOMS):
                orig_carts.append([CARTESIANS[atom][0], CARTESIANS[atom][1], CARTESIANS[atom][2]])

            # could get rid of atomic units here, if zpe_rat definition is changed
            for mode, wn in enumerate(FREQS):
                # Either moves along any and all imaginary freqs, or a specific mode requested by the user
                if FREQS[mode] < 0.0:
                    shift.append(amplitude)
                else: shift.append(0.0)

            # The starting geometry is displaced along the each normal mode according to the random shift
                for atom in range(0,NATOMS):
                    for coord in range(0,3):
                        CARTESIANS[atom][coord] = CARTESIANS[atom][coord] + NORMALMODE[mode][atom][coord] * shift[mode]
        outfile.close()

        # This part places the calculations in different folders depending on the type of
        # termination and number of imag. freqs
        source = w_dir+file

        if IM_FREQS > 0:
            destination = w_dir+'Imaginary frequencies/'+file
            try:
                os.makedirs(destination)
                shutil.move(source, destination)
            except OSError:
                if  os.path.isdir(destination):
                    pass
                else:
                    raise

        if IM_FREQS == 0 and TERMINATION == "normal":
            destination = w_dir+'Finished/'+file
            try:
                os.makedirs(destination)
                shutil.move(source, destination)
            except OSError:
                if  os.path.isdir(destination):
                    pass
                else:
                    raise

        if IM_FREQS == 0 and TERMINATION == "error":
            destination = w_dir+'Failed_Error/'+file
            try:
                os.makedirs(destination)
                shutil.move(source, destination)
            except OSError:
                if  os.path.isdir(destination):
                    pass
                else:
                    raise

        if IM_FREQS == 0 and TERMINATION == "unfinished":
            destination = w_dir+'Failed_Unfinished/'+file
            try:
                os.makedirs(destination)
                shutil.move(source, destination)
            except OSError:
                if  os.path.isdir(destination):
                    pass
                else:
                    raise

        if IM_FREQS > 0 or TERMINATION != "normal":
            # Settings for the com files
            n_of_processors = '24'
            memory = '96GB'
            functional = lot
            basis_set = bs
            basis_set_I = 'LANL2DZ'
            # Options for genecp
            ecp_list,ecp_I = [],False
            possible_atoms = ['N', 'P', 'As', 'C', 'Si', 'Ge', 'B', 'H', 'S', 'O', 'Se', 'F', 'Br', 'Cl', 'I']
            for i in range(len(ATOMTYPES)):
                if ATOMTYPES[i] not in ecp_list and ATOMTYPES[i] in possible_atoms:
                    ecp_list.append(ATOMTYPES[i])
                if ATOMTYPES[i] == 'I':
                   ecp_I = True
            if ecp_I == False:
                genecp = 'gen'
            if ecp_I == True:
                genecp = 'genecp'
            solvent = ''
            keywords_opt = functional+'/'+genecp+' '+solvent+' opt=(calcfc,maxstep=5) freq=noraman '

            print(os.getcwd())

            fileout = open(file.split(".")[0]+'.com', "w")
            fileout.write("%mem="+memory+"\n")
            fileout.write("%nprocshared="+n_of_processors+"\n")
            fileout.write("# "+keywords_opt+"\n")
            fileout.write("\n")
            fileout.write(name+"\n")
            fileout.write("\n")
            fileout.write(str(CHARGE)+' '+str(MULT)+'\n')
            for atom in range(0,NATOMS):
                fileout.write('{0:>2} {1:12.8f} {2:12.8f} {3:12.8f}'.format(ATOMTYPES[atom], CARTESIANS[atom][0],  CARTESIANS[atom][1],  CARTESIANS[atom][2]))
                fileout.write("\n")
            fileout.write("\n")
            for i in range(len(ecp_list)):
                if ecp_list[i] != 'I':
                    fileout.write(ecp_list[i]+' ')
            fileout.write('0\n')
            fileout.write(basis_set+'\n')
            fileout.write('****\n')
            if ecp_I == False:
                fileout.write('\n')
            else:
                fileout.write('I     0\n')
                fileout.write(basis_set_I+'\n')
                fileout.write('****\n\n')
                fileout.write('I 0\n')
                fileout.write(basis_set_I+'\n\n')
            fileout.close()
