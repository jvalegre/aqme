import os
import subprocess
import shutil
import shlex
import re
from collections import namedtuple
import glob
from pathlib import Path

# ensure a proper import of pybel for openbabel v2 and v3
try:
	import pybel
except ImportError:
	from openbabel import pybel # for openbabel>=3.0.0


class TurbomoleOutput(object):
    """
    Represents the output of a turbomole calculation in a specified folder.
    
    Parameters
    ----------
    folder : str
        A valid path to the folder where the turbomole calculation has 
        been carried out or is in the process.
    T : float, optional
        Temperature in K for thermochemistry calculations, by default 298.15
    P : float, optional
        Pressure in MPa for thermochemistry calculations, by default 0.1
    scale_factor : float, optional
        scaling factor for thermochemistry calculations between 0 and 1, 
        by default 1.0
    qh_thermo : bool, optional
        If True the a qh rrho formalism will be used instead of the rrho 
        formalism to calculate the thermochemistry, by default False
    fcutoffs : tuple, optional
        lower and upper limits frequency thresholds to the qh calculations in 
        cm^-1, by default (100,150)

    Attributes
    ----------
    calc_type
    calc_state
    cosmo
    geometry
    energy
    frequencies
    rotational_constants
    zpe
    enthalpy
    entropy
    gibbs
    """
    def __init__(self,folder,T=298.15,P=0.1,scale_factor=1.0,qh_thermo=False,
                 fcutoffs=(100,150)):
        self.folder = Path(folder)
        self.calc_type = self.guess_calculation_type()
        self.calc_state = self.guess_calculation_state()
        self.cosmo = self.uses_cosmo()
        self.geometry = self.parse_geometry()
        self.energy = self.parse_energy()
        self.T = T # Temperature in K
        self.P = P # pressure in MPa
        self.scale_factor = scale_factor # Scaling factor for frequencies 
        self.qh_thermo = qh_thermo       # If the qh free energy is to be calculated 
        self.fcutoffs = fcutoffs         # Cutoffs for frequencies for qh in cm^-1
        self.frequencies = self.parse_frequencies()
        self.rotational_constants = self.parse_rotational_constants()
        self._zpe = None
        self._enthalpy = None
        self._entropy = None
        self._gibbs = None
        self.calc_thermo()
    def guess_calculation_type(self):
        """
        Guesses the type of calculation using the existence/nonexistence of 
        certain files. 

        Returns
        -------
        str
            4 possible outcomes: 'sp','opt','freq', 'opt+freq'
        """
        calc_type = []
        file = self.folder / 'ridft.out'
        if file.exists():
            return 'sp'
        file = self.folder / 'jobex.out'
        if file.exists():
            calc_type.append('opt') 
        file = self.folder / 'vibspectrum'
        if file.exists():
            calc_type.append('freq')
        return '+'.join(calc_type)
    def guess_calculation_state(self):
        """
        Guesses the current state of the calculation. 

        Returns
        -------
        str
            1 of the following states: 'running','converged','failed'
        """
        if self.calc_type == 'sp': 
            return self._guess_calculation_state_sp()
        else:
            return self._guess_calculation_state_optfreq()
    def _guess_calculation_state_sp(self):
        """
        Guesses the current state of the calculation depending on the last line
        of the ridft.out file.
        """
        file = self.folder / 'ridft.out'
        with open(file,'r') as F: 
            lines = [line.strip() for line in F if line.strip()]

        for line in lines: 
            if 'ridft did not converge!' in line: 
                return 'failed'
        
        if 'ridft ended normally' in lines[-1]:
            return 'converged'
        
        return 'running' 
    def _guess_calculation_state_optfreq(self):
        """
        Guesses the current state of the calculation depending on the existence
        of the 'GEO_OP_XXXX' file. 

        Returns
        -------
        str
            3 accepted states: 'running','converged','failed'
        """
        running_file = self.folder / 'GEO_OPT_RUNNING'
        if running_file.exists():
            # Check for out of time
            # if out of time -> return 'failed'
            # else:
            return 'running'
        for item in ['converged','failed']: 
            file = self.folder / f'GEO_OPT_{item.upper()}'
            if file.exists():
                return item
        else:
            raise RuntimeError('No state was found for the calculation')

    def parse_geometry(self):
        """
        Finds the last geometry from de coord file and returns its representation
        as a .xyz file. 

        Returns
        -------
        str
            contents of the xyz version of the coord file without title.
        """
        file = str(self.folder/'coord')
        mol = next(pybel.readfile('tmol',file))
        return mol.write('xyz')
    def parse_energy(self):
        """
        Reads the potential energy from the appropiate location depending on the
        solvation and type of calculation. 

        Returns
        -------
        float
            Final potential energy
        """
        pattern = r'\|\s*?total\s*?energy\s*?\=\s*?(\-+[0-9]*\.[0-9]*\s*?\|)'
        if self.calc_type == 'sp': 
            file = Path('ridft.out')
        else:
            file = Path('job.last')
        with open(self.folder / file,'r') as F: 
            txt = F.read()
        if self.cosmo: 
            pattern = r'Total\senergy\s\+\sOC\scorr\.\s\=\s*?(\-+[0-9]*\.[0-9]*)\n'
        else:
            pattern = r'\|\s*?total\s*?energy\s*?\=\s*?(\-+[0-9]*\.[0-9]*)\s*?\|'
        regex = re.compile(pattern)
        matches = regex.findall(txt)
        if matches: 
            return float(matches[-1])
    def uses_cosmo(self):
        """
        Checks if the calculation uses cosmo solvation. Returns True if it finds
        the 'out.cosmo' file.
        """
        file = self.folder / 'out.cosmo'
        if file.exists(): 
            return True
        # Check in the control file
        regex = re.compile('cosmo') 
        with open(self.folder/'control','r') as F: 
            txt = F.read()
        match = regex.findall(txt)
        return bool(match)
    def calc_thermo(self,executable=None):
        """
        Runs the appropiate executable to calculate the thermochemistry to 
        calculate the zpe, H, S and G at T and P and parses the output to 
        retrieve said magnitudes. 

        Parameters
        ----------
        executable : str
            Absolute path to the executable script in case that it is not 
            accesible with shutils.which. Grimme's "thermo" script for qh and
            Turbomole's "freeh" for non-qh.
        """
        if self.qh_thermo:
            zpe,enthalpy,entropy,gibbs = self._calc_thermo_qh(executable)
        else:
            zpe,enthalpy,entropy,gibbs = self._calc_thermo_noqh(executable)
        self._zpe = zpe
        self._enthalpy = enthalpy
        self._entropy = entropy
        self._gibbs = gibbs
    def _calc_thermo_noqh(self,executable=None):
        """
        Runs the 'freeh' command line utility provided that the it is accesible
        to calculate the zpe, H, S and G at T and parses the output to retrieve 
        said magnitudes. 

        Parameters
        ----------
        executable : str
            Absolute path to the executable script in case that it is not 
            accesible with shutils.which.

        Returns
        -------
        zpe,enthalpy,entropy,gibbs
            Energy corrections to the potential energy in hartree if the output 
            of freeh is in kJ/mol, otherwise the units are the same as the ones 
            provided by the freeh script. (Returns (None,)*4 if the freeh 
            command is not accesible)
        """
        freeh = shutil.which('freeh')  # turbomole's executable
        if freeh is None and executable is None: 
            return None, None, None, None
        elif executable is not None: 
            freeh = executable
        # Run the freeh command
        cwd = os.path.abspath(os.getcwd())
        os.chdir(self.folder)
        tstart = tend = self.T
        pstart = pend = self.P
        keywords = f'tstart={tstart} tend={tend} numt=1 pstart={pstart} pend={pend} nump=1'
        echo_cmd = f'echo "\n1\n{keywords}\n*\n"'
        cmd = freeh
        p1 = subprocess.Popen(shlex.split(echo_cmd), stdout=subprocess.PIPE)
        p2 = subprocess.Popen(shlex.split(cmd), stdin=p1.stdout, 
                              stdout=subprocess.PIPE,stderr=subprocess.DEVNULL)
        p1.stdout.close()
        txt = [line.decode('utf-8') for line in p2.stdout]
        p2.stdout.close()
        os.chdir(cwd)
        # Here do the parsing of txt
        _iter = txt.__iter__()
        for line in _iter:
            if 'your wishes are' in line: 
                    break
        # read zpe
        for line in _iter: 
            if 'zpe=' in line:
                _, zpe, zpe_unit = line.strip().split()
                zpe = float(zpe)
                break
        else: 
            zpe = None
            zpe_unit = None
        for line in _iter:
            if 'entropy' in line:
                items = line.strip().split()
                index_e = items.index('entropy')
                index_g = items.index('chem.pot.')
                _ = next(_iter)
                _ = next(_iter)
                line = next(_iter)
                entropy = float(line.strip().split()[index_e])
                gibbs = float(line.strip().split()[index_g])
                break
        else:
            entropy = gibbs = None
        for line in _iter: 
            if 'enthalpy' in line:
                items = line.strip().split()
                index = items.index('enthalpy')
                _ = next(_iter)
                line = next(_iter)
                enthalpy= float(line.strip().split()[index])
                break
        else:
            enthalpy = None
        # We assume congruent units of all properties: 
        if zpe_unit == 'kJ/mol':
            factor = 2625.5 #kJ/mol -> a.u.
            zpe = zpe/factor
            enthalpy = enthalpy/factor
            entropy = entropy/factor
            gibbs = entropy/factor
        return zpe,enthalpy,entropy,gibbs
    def _calc_thermo_qh(self,executable=None):
        """
        Runs the 'thermo' script (developed by Grimme's) provided that the it 
        accesible to calculate the zpe, H, S and G at T and parses the output to
        retrieve said magnitudes. 

        Parameters
        ----------
        executable : str
            Absolute path to the executable script in case that it is not 
            accesible with shutils.which

        Returns
        -------
        zpe,enthalpy,entropy,gibbs
            Energies in kcal/mol if the output of freeh is in kJ/mol, otherwise
            the units are the same as the ones provided by the freeh script. 
            (Returns (None,)*4 if the freeh command is not accesible)
        """
        thermo = shutil.which('thermo') # Grimme's script
        if thermo is None and executable is None: 
            return None, None, None, None
        elif executable is not None: 
            thermo = executable
        cwd = os.path.abspath(os.getcwd())
        os.chdir(self.folder)
        cutoff1,cutoff2 = self.fcutoffs
        cmd = f'{thermo} {cutoff1} {cutoff2} {self.T} {self.scale_factor}\n'
        p1 = subprocess.Popen(shlex.split(cmd), 
                              stdout=subprocess.PIPE,
                              stderr=subprocess.DEVNULL)
        zpe = enthalpy = entropy = gibbs = None
        for bytes in p1.stdout:
            line = bytes.decode('utf-8')
            if 'ZPVE' in line: 
                zpe = float(line.strip().split()[1])
            if 'H(T)' in line and 'H(0)' not in line: 
                enthalpy = float(line.strip().split()[1])
            if 'T*S' in line and len(line.split())==4: 
                entropy = float(line.strip().split()[1])/self.T
            if 'G(T)' in line: 
                gibbs = float(line.strip().split()[1])
        p1.stdout.close()
        p1.wait()
        if Path('.H298').exists():
            os.remove('.H298')
        if Path('.G298').exists():
            os.remove('.G298')
        os.chdir(cwd)
        return zpe, enthalpy, entropy, gibbs
    def parse_frequencies(self):
        """
        Reads the frequency mode and wavenumber if the turbomole calculation 
        has frequencies.  

        Returns
        -------
        list
            list of tuples where the first element is the number of the 
            vibrational mode and the second element is the wavenumber in cm^-1 
        """
        Frequency = namedtuple('frequency','mode wavenumber')
        if 'freq' not in self.calc_type:
            return []
        with open(self.folder/'vibspectrum','r') as F: 
            _iter = F.__iter__()
            line = next(_iter)
            while "$vibrational spectrum" not in line: 
                line = next(_iter)
            line = next(_iter)
            while line.startswith('#'):
                line = next(_iter)
            lines = []
            while '$end' not in line: 
                lines.append(line.strip().split())
                line = next(_iter)
        frequencies = [Frequency(int(line[0]),float(line[2])) for line in lines if len(line)==6]
        return frequencies
    def parse_rotational_constants(self):
        """
        Reads the rotational constants in MHz and returns them in GHz.

        Returns
        -------
        list
            A three element list with the rotational constants around the principal
            moments of inertia in GHz. 
        """
        if self.cosmo: 
            file = self.folder / 'numforce/aoforce.out'
        else:
            file = self.folder / 'aoforce.out'
        if not file.exists(): 
            return []
        pattern = r'.*?b.*?\:(.*)\s*\(MHz\)'
        with open(file,'r') as F: 
            txt = F.read() 
        regex = re.compile(pattern)
        match = regex.findall(txt)
        if match: 
            line = match[-1]
            items = line.split()
            return [float(item)/1000.0 for item in items]
        return []

    @property
    def zpe(self):
        if self._zpe is None: 
            return None
        return self.energy + self._zpe
    @property
    def enthalpy(self):
        if self._enthalpy is None: 
            return None
        return self.energy + self._enthalpy
    @property
    def entropy(self):
        if self._entropy is None: 
            return None
        return self.energy + self._entropy
    @property
    def gibbs(self):
        if self._gibbs is None: 
            return None
        return self.energy + self._gibbs

class TurbomoleInput(object):
    """
    Turbomole Input Generator class. This class provides a simplified interface 
    for automating the usage of the 'define' command of Turbomole. For more 
    details on the options check turbomole define's manual.

    Parameters
    ----------
    folder : str
        Valid path to the folder where the input is going to be generated.
        If the folder does not exist it will be created. If parent folders 
        do not exist it will raise an error.  
    functional : str, optional
        dft functional to use or 'hf' or 'hf-3c', by default 'TPSS'. Ensure that
        the functional is written as in Turbomole. 
    basis : str or dict, optional
        basis set to use for all atoms, by default 'def2-SVP'.
    auxbasis : str or dict, optional
        auxiliar basis set, by default it attempts to use the same as the 
        selected basis set.
    ecp : dict, optional
        dictionary of atom -> ecp for the calculation.
    dispersion : str, optional
        Whether to include dft dispersion and which kind, by default 'off'. 
    charge : int, optional
        Formal Charge, by default 0
    multiplicity : int, optional
        Overall multiplicity 's','d','t' are also accepted as arguments but 
        integer argument is preferred, by default 1
    grid : str, optional
        DFT integration grid, by default 'm4'
    epsilon : float, optional
        Dielectric constant of the cosmo solvation or 'gas' or 'infinity', 
        by default 'gas'
    maxcore : int, optional
        Maximum memory per core in MB, by default 200
    ricore : int, optional
        Maximum memory per core in MB for the RI formalism, by default 200
    disable_ired : bool, optional
        If True it will not generate internal redundant coordinates, 
        by default False
    cavity : str, optional
        Solvent cavity definition for cosmo, by default 'fine'
    title : str, optional
        title of the calculation, by default ''
    """
    def __init__(self,folder, functional='TPSS',basis='def2-SVP', 
                 auxbasis=None,ecp=None,dispersion='off',charge=0,
                 multiplicity=1,grid='m4',epsilon='gas',maxcore=200,
                 ricore=200,disable_ired=False,cavity='none',title=''):
        self.folder = Path(folder)
        self.ofile = self.folder/f'{self.folder.stem}.sh' 
        self.cwd = Path(os.path.abspath(os.getcwd()))

        self.title = title

        self.functional = functional
        self.basis = basis
        self.auxbasis = auxbasis
        self.ecp = ecp
        self.dispersion = dispersion
        self.charge = charge
        # initialize multiplicity
        self.multiplicity = multiplicity
        self.check_multiplicity()
        self.grid = grid
        self.epsilon = epsilon
        self.maxcore = maxcore
        self.ricore = ricore
        self.disable_ired = disable_ired
        self.cavity = cavity

        self.atoms = []
        self.check_atom_list()
        self.check_auxbasis_availability()
    
    def __str__(self):
        output = f'\nInformation for directory\n\n{self.folder}\n\n'
        output += f'Functional: {self.functional}\n'
        output += f'Basis set: {self.basis}\n'
        output += f'Auxilary basis set: {self.auxbasis}\n'
        if self.dispersion == 'on': 
            disp = 'd3(0)'
        else: 
            disp = self.dispersion
        output += f'Dispersion: {disp}\n'
        output += f'Charge: {self.charge}\n'
        output += f'Multiplicity: {self.multiplicity}\n'
        output += f'ricore: {self.ricore}\n'
        output += f'maxcore: {self.maxcore}\n'
        output += f'Grid: {self.grid}\n'
        #Infinity, gas or a number a.bcd
        if not self.cavity :
            output += f"COSMO settings: {self.epsilon} with {self.cavity}\n"
        else :
            output += f"COSMO settings: {self.epsilon}\n"

        return output

    # Initialization methods
    def check_multiplicity(self):
        """
        Checks the 'multiplicity' attribute and ensures proper values for the 
        'multiplicity', 'is_closed_shell' and 'unpaired_electrons' attributes. 
        """
        # Ensure that multiplicity is a integer
        m = self.multiplicity
        if m == 's': 
            multiplicity = 1
        elif m == 'd':
            multiplicity = 2
        elif m == 't': 
            multiplicity = 3
        else:
            multiplicity = int(m)
        # Calculate check for restricted/unrestricted and number of unpaired 
        # electrons
        if multiplicity != 1:
            is_closed_shell = False 
            unpaired_electrons = multiplicity-1
        else:
            is_closed_shell = True
            unpaired_electrons = 0
        # Assign the new object attributes
        self.multiplicity = multiplicity
        self.is_closed_shell = is_closed_shell
        self.unpaired_electrons = unpaired_electrons
    def check_atom_list(self):
        """
        Gets a list of the different atoms present in the molecule in the 'atoms'
        attribute.
        """
        coord_file = self.folder/'coord' 
        if coord_file.exists():
            mol = next(pybel.readfile('tmol',str(coord_file)))
            lines = mol.write('xyz').split('\n')
        else: 
            xyzfiles = glob.glob(str(self.folder)+'/*.xyz')
            if xyzfiles: 
                with open(xyzfiles[0],'r') as F: 
                    lines = F.read().split('\n')
            else: 
                return

        atoms = set([line.strip().split()[0] for line in lines[2:] if line.strip()])
        self.atoms = [element.lower() for element in list(atoms)]
    def check_auxbasis_availability(self):
        """
        Checks availability of the auxiliar basis sets in the jbasen and jkbasen
        directories in of turbomole. If the auxiliar basis set is not specified
        it assumes that the auxiliar basis set = basis set. 

        Raises
        ------
        ValueError
            If a basis set for a certain set of atoms is not specified.
        """
        if self.auxbasis is None and self.basis == 'minix':  
            self.auxbasis = 'universal'
        elif self.auxbasis is None: 
            self.auxbasis = self.basis

        if type(self.auxbasis) == str:
            self._check_auxbasis_availability_all() 
        else:
            self._check_auxbasis_availability_per_atom()

    def _check_auxbasis_availability_all(self):
        """
        Checks availability of the auxiliar basis sets in the jbasen and jkbasen
        directories in of turbomole. If the auxiliar basis set is not specified
        it assumes that the auxiliar basis set = basis set. Assumes both the 
        basis and auxiliar basis are specified for "all" atoms. 

        Raises
        ------
        ValueError
            If a basis set for a certain set of atoms is not specified.
        """
        turbodir = Path(os.environ['TURBODIR'])
        self.use_jkbasis_as_jbasis = False
        found_list  = []
        for atom in self.atoms: 
            jbasen = turbodir/f'jbasen/{atom}'
            jkbasen = turbodir/f'jkbasen/{atom}'
            basis_name = f'{atom} {self.auxbasis}'
            found_jbasis = self.find_basis_in_file(basis_name,jbasen)
            found_jkbasis = self.find_basis_in_file(basis_name,jkbasen)
            found_list.append((atom,found_jbasis,found_jkbasis))
       
        errors = []
        for atom,found_jbasis,found_jkbasis in found_list:
            if not found_jbasis and not found_jkbasis:
                errors.append((atom,self.auxbasis))
            elif not found_jbasis:
                self.use_jkbasis_as_jbasis = True

        if errors:
            msg = ','.join([f"'{atom} {basis}'" for atom,basis in errors])
            raise ValueError(f'Basis sets [{msg}] are not found in jkbasen nor jbasen')
    def _check_auxbasis_availability_per_atom(self):
        """
        Checks availability of the auxiliar basis sets in the jbasen and jkbasen
        directories in of turbomole. If the auxiliar basis set is not specified
        it assumes that the auxiliar basis set = basis set. Assumes both the 
        basis and auxiliar basis are specified per each atom as a dictionary. 

        Raises
        ------
        ValueError
            If a basis set for a certain set of atoms is not specified.
        """
        turbodir = Path(os.environ['TURBODIR'])
        self.use_jkbasis_as_jbasis = False
        found_list  = []
        for atom,auxbasis in self.auxbasis.items():
            sym = atom.lower()
            jbasen = turbodir/f'jbasen/{sym}'
            jkbasen = turbodir/f'jkbasen/{sym}'
            basis_name = f'{sym} {auxbasis}'
            found_jbasis = found_jkbasis = False
            if jbasen.exists():
                found_jbasis = self.find_basis_in_file(basis_name,jbasen)
            if jkbasen.exists():
                found_jkbasis = self.find_basis_in_file(basis_name,jkbasen)
            found_list.append((sym,auxbasis,found_jbasis,found_jkbasis))
       
        errors = []
        for sym,auxbasis,found_jbasis,found_jkbasis in found_list:
            if not found_jbasis and not found_jkbasis:
                errors.append((sym,auxbasis))
            elif not found_jbasis:
                self.use_jkbasis_as_jbasis = True

        if errors:
            msg = ','.join([f"'{atom} {basis}'" for atom,basis in errors])
            raise ValueError(f'Basis sets [{msg}] are not found in jkbasen nor jbasen')

    @staticmethod
    def find_basis_in_file(basis,file):
        """
        Attempts to find a basis set in the specified file. 

        Parameters
        ----------
        basis : str
            basis set in the format '[sym] [basis name]' where the symbol is 
            lowercased. 
        file : Path
            Path to the file where the basis is to be searched. 

        Returns
        -------
        bool
            True if the basis is found within the file.
        """
        assert file.exists()
        with open(file,'r') as F: 
            for line in F: 
                if basis in line:
                    return True        
        return False
    
    # Generation methods
    def write_input_command(self,script):
        """
        Writes a file named 'input_command' with a line of text that contains 
        the necessary options used for the generation of the input file. 
        It is legacy code for input generation scripts.  

        Parameters
        ----------
        script : str
            name of the python script used to generate the turbomole input. If 
            none it uses this file's name. 
        """
        if script is None: 
            script = __file__
        options = [('f',self.functional),
                   ('bs',self.basis),
                   ('d',self.dispersion),
                   ('c',self.charge),
                   ('u',self.multiplicity),
                   ('rm',self.ricore),
                   ('m',self.maxcore),
                   ('g',self.grid),
                   ('s',self.epsilon)]
        if not self.disable_ired: 
            options.append(('ired',''))
        if self.cavity == 'fine': 
            options.append(('fine',''))
        
        input_str = f'#python {script} '
        input_str += ' '.join([f'-{k} {val}' for k,val in options])
        with open(self.folder/'input_command','w') as F: 
            F.write(input_str)
    def create_coord(self):
        """
        Ensures the existence of a 'coord' file. Uses pybel to transform from 
        xyz to tmol format in case that a 'coord' file does not exist but a .xyz
        does exist.
        In case of several .xyz files (not recommended), it will use the first 
        one in alphabetical order. 

        Raises
        ------
        Run
            If no xyz or coord file is found in the folder it will stop 
            completely. 
        """
        coord_file = self.folder/'coord' 
        if coord_file.exists():
            return 
        #print("File coord was not found in the folder.")
        xyzfiles = glob.glob(str(self.folder)+'/*.xyz')
        if xyzfiles:
            #print("File coord will be generated from an .xyz-file.")
            mol = next(pybel.readfile('xyz',xyzfiles[0]))
            mol.write('tmol',str(self.folder/'coord'))
        else:
            msg = "There is no .xyz-file either in the folder. Exiting the script."
            raise RuntimeError(msg)
    
    # define's interfacing methods
    def write_geometry_menu(self):
        """
        Encapsulates the interactions with 'define' related with the 
        geometry menu and returns a string with those interactions.
        """
        output = 'a coord\n'         # Add coordinates file 'coord'
        if not self.disable_ired:    # 
            output += '*\nno\n'
        else:
            output += 'ired\n*\n'
        return output
    def write_atom_attr_menu(self):
        """
        Encapsulates the interactions with 'define' related with the 
        Atom Attributes menu (Essentially, Isotope definition and Basis set 
        definition) and returns a string with those interactions.
        """
        if type(self.basis) == str: 
            return f'b all {self.basis}\n*\n' # Asign the basis set to all atoms and exit
        output = ''
        for key,val in self.basis.items(): 
            output += f'b "{key}" {val}\n'        # Add basis set per atom
        if self.ecp is not None: 
            ecps = [(k,val) for k,val in self.ecp.items()]
            output += 'ecp\n'                     # Enter ecp menu
            if ecps:
                for key,val in ecps:
                    output += f'"{key}" {val}\n'  # Add ecps
            output += '\n'                        # Exit the ecp menu
        output += '*\n'                           # Exit the atom attribute menu
        return output
    def write_occupation_MO_menu(self):
        """
        Encapsulates the interactions with 'define' related with the 
        Occupation and Molecular Orbitals menu and returns a string with those 
        interactions.
        """
        output = ''
        if 'cu' in self.atoms: 
            # RAUL: Maybe there is an extra newline between eht and charge
            output = f'eht\n\n\n{self.charge}\n'
        else:
            output = f'eht\n\n{self.charge}\n'

        if not self.is_closed_shell and self.unpaired_electrons != 0:
            output += f'n\nu {self.unpaired_electrons}\n*\n\n'
        elif not self.is_closed_shell:
            output += f'n\n{self.unpaired_electrons}\n*\n\n'
        else:
            output += f'y\n'

        return output
    
    def write_dft_submenu(self):
        """
        Encapsulates the interactions with 'define' related with the 
        DFT submenu within the General Menu and returns a string with those 
        interactions.
        """
        if self.functional in ['hf','hf-3c']:
            return ''
        output = 'dft\n'                         # Enter the DFT menu
        output += 'on\n'                         # Enable DFT
        output += f'func {self.functional}\n'    # Set the functional
        output += f'grid {self.grid}\n'          # Set the integration grid
        output += f'\n'                          # End the DFT menu
        return output
    def write_auxbasis_submenu(self):
        """
        Encapsulates the interactions with 'define' related with the 
        Auxiliary Basis sets submenu within the RI menu and returns a 
        string with those interactions.
        """
        if type(self.auxbasis) == str: 
            basis_definition = ''
            if self.basis != self.auxbasis: 
                for atom in self.atoms:
                    basis_definition += f'nick {atom} {self.auxbasis}\n'
            basis_definition += f'b all {self.auxbasis}\n'
        else:
            basis_definition = ''
            if self.basis != self.auxbasis: 
                for atom,auxbasis in self.basis.items():
                    basis_definition += f'nick {atom} {auxbasis}\n'
            for atom,auxbasis in self.basis.items():
                basis_definition += f'b "{atom}" {auxbasis}\n'
        return f'jbas\n{basis_definition}*\n'
    def write_ri_submenu(self):
        """
        Encapsulates the interactions with 'define' related with the 
        RI submenu within the General menu and returns a string with those 
        interactions.
        """
        # TODO LEGACY CODE, this comment should be removed when we are 
        # sure that it is definetly not usefull. 
        ## hybrid_or_HF = ['b2-plyp','b3-lyp','tpssh','pw6b95',
        ##                 'pbe0','m06-2x','bhlyp','hf','hf-3c']
        ## is_hybrid_or_HF = self.functional in hybrid_or_HF

        if self.use_jkbasis_as_jbasis:
            output  = 'rijk\n'    # Enter the RI-JK-HF menu
            output += 'on\n'      # Activate it
            #output += f'm {self.ricore}\n'   # Set the maximum memory per core in MB
            output += '\n'        # Exit the RI-JK-HF menu
            return output

        auxbasis_menu = self.write_auxbasis_submenu()
        output  = 'ri\n'                 # Enter the RI menu
        output += 'on\n'                 # Activate it 
        output += f'm {self.ricore}\n'   # Set the maximum memory per core in MB
        output += f'{auxbasis_menu}'     # Enter and exit the auxbasis submenu
        output += '\n'                   # Exit the RI menu 
        return output    
    def write_multipole_menu(self):
        """
        Encapsulates the interactions with 'define' related with the 
        Multipole submenu within the General menu and returns a string with 
        those interactions.
        """
        output = 'marij\n'   # Enter the multipole RI menu
        output += '\n'       # Accept defaults and exit
        return output
    def write_dispersion_menu(self):
        """
        Encapsulates the interactions with 'define' related with the 
        Dispersion submenu within the General menu and returns a string with 
        those interactions.
        """
        output = 'dsp\n'                    # Enter the dispersion menu
        output += f'{self.dispersion}\n'    # Set the dispersion
        output += '\n'                      # End the menu
        return output
    def write_ricc_menu(self):
        """
        Encapsulates the interactions with 'define' related with the 
        RICC submenu within the General menu and returns a string with those 
        interactions.
        """
        output = 'cc\n'                     # Enter the ricc menu
        output += f'memory {self.maxcore}'  # Set the maximum memory per core
        output += '*\n'                     # End the ricc menu
        return f'cc\nmemory {self.maxcore}\n*\n'
    def write_scf_menu(self):
        """
        Encapsulates the interactions with 'define' related with the 
        SCF submenu within the General menu and returns a string with those 
        interactions.
        """
        output  = 'scf\n'   # Enter the scf menu
        output += 'iter\n'  # Select the max iterations property for editing
        output += '300\n'   # Set it to 300 iterations
        output += '\n'      # End the scf menu 
        return output
    def write_general_menu(self):
        """
        Encapsulates the interactions with 'define' related with the 
        the General menu and centralizes all the interactions of its submenus. 
        returns a string with all interactions.
        """
        dft_menu = self.write_dft_submenu()
        ri_menu = self.write_ri_submenu()
        multipole_menu = self.write_multipole_menu()
        dispersion_menu = self.write_dispersion_menu()
        ricc_menu = self.write_ricc_menu()
        scf_menu = self.write_scf_menu()
        
        general_menu = ''
        general_menu += dft_menu
        general_menu += ri_menu
        general_menu += multipole_menu
        general_menu += dispersion_menu
        general_menu += ricc_menu
        general_menu += scf_menu
        general_menu += '*\n'

        return general_menu

    def generate_cosmoprep_string(self):
        """
        Centralizes the generation of all the interactions with the 'cosmoprep'
        command and returns them in the appropiate format as a single string.
        """
        if self.epsilon == 'gas': 
            return ''                # No-solvation
        output  = '\n'*12
        output += 'r all b\n'
        output += '*\n'
        output += 'out.cosmo\n'
        output += 'n\n'
        output += '\n'
        return output
    def generate_define_string(self):
        """
        Centralizes the generation of all the interactions with the 'define'
        command and returns them in the appropiate format as a single string. 
        """
        molecular_geom = self.write_geometry_menu()
        atom_attrib = self.write_atom_attr_menu()
        ocupation_MO = self.write_occupation_MO_menu()
        general_menu = self.write_general_menu()

        exec_str = f'\n{self.title}\n'
        exec_str += f'{molecular_geom}'
        exec_str += f'{atom_attrib}'
        exec_str += f'{ocupation_MO}'
        
        exec_str += general_menu
        exec_str += '*\n'
        return exec_str
    
    # Post-generation modification methods
    def modify_generated_control_file(self):
        """
        Modifies the generated control file: 
        - Adds a scftol of 1d-16
        - Adds the cosmo_isorad keyword when needed
        - Ensures that if the jkbasis were used, they are treated as jbasis and 
          removes the rik keyword
        - Ensures proper dispersion treatment with hf-3c.  
        """
        with open(self.folder/'control','r') as F:
            txt = F.read()
        
        head = f'$scftol 1d-16\n'
        if self.cavity == 'fine':
            head += '$cosmo_isorad\n'
        
        txt = head + txt
        
        if self.use_jkbasis_as_jbasis: 
            newtxt = txt.replace('jkbas','jbas')
            lines = [line for line in newtxt.split('\n') if '$rik' not in line]
            txt = '\n'.join(lines)

        if self.functional == 'hf-3c': 
            txt = txt.replace('disp3 bj','disp3 bj func hf-3c')
        
        with open(self.folder/'control','w') as F: 
            F.write(txt)
    def modify_generated_auxbasis_file(self):
        """
        Modifies the generated auxbasis file:
        - If no jkbases were used, it does nothing
        - Otherwise ensures that any jkbasis is treated as jbasis 
        """
        if not self.use_jkbasis_as_jbasis:
            return
        
        with open(self.folder/'auxbasis','r') as F:
            txt = F.read()

        newtxt = txt.replace('jkbas','jbas')
        
        with open(self.folder/'auxbasis','w') as F: 
            F.write(newtxt)
    def modify_generated_files(self):
        """
        Applies custom post-modifications to generated files. Currently only 
        the control and auxbasis files might be modified. 
        """
        self.modify_generated_control_file()
        self.modify_generated_auxbasis_file()

    # Main API
    def generate(self):
        """
        Ensures the existence of the coord file, and generates the input 
        interacting with 'define' and 'cosmoprep' commands. Requires that 
        'define' and 'cosmoprep' are accessible at the current path.  
        """
        self.create_coord()
        define_str = self.generate_define_string()
        cosmo_str = self.generate_cosmoprep_string()
        # the define command requires to be started in the folder of the 
        # calculation
        os.chdir(self.folder)
        # Run the define and the cosmoprep processes
        p_echo1 = subprocess.Popen(shlex.split(f'echo "{define_str}"'),
                                              stdout=subprocess.PIPE)
        p_define = subprocess.Popen(shlex.split('define'),
                                    stdin=p_echo1.stdout, 
                                    stdout=subprocess.PIPE,
                                    stderr=subprocess.DEVNULL)
        p_echo1.wait()
        p_define.wait()
        p_echo1.stdout.close()
        p_define.stdout.close()

        if cosmo_str:
            p_echo2 = subprocess.Popen(shlex.split(f'echo "{cosmo_str}"'),
                                                  stdout=subprocess.PIPE)
            p_cosmoprep = subprocess.Popen(shlex.split('cosmoprep'),
                                           stdin=p_echo2.stdout, 
                                           stdout=subprocess.PIPE,
                                           stderr=subprocess.DEVNULL)
            p_echo2.wait()
            p_echo2.stdout.close()
            p_cosmoprep.wait()
            p_cosmoprep.stdout.close()
        
        # After the define and cosmoprep are finished return to the original 
        # directory 
        os.chdir(self.cwd)

        # Apply the post-modifications to the generated files (will look for the
        # files in the folder, so no need to be in the same directory).
        self.modify_generated_files()
    

def read_basiset_file(filepath):
    """
    Reads a file containing information on basis sets, auxiliar basis sets and 
    ECPs. The file follows the following format: 
    $basis
    sym1 basis1
    sym2 basis2
    ...
    $auxbasis
    sym1 auxbasis1
    sym2 auxbasis2
    ...
    $ecp
    sym1 ecp1
    sym2 ecp2
    ...
    $end

    Parameters
    ----------
    filepath : str or Path
        a valid path to a file with the previously specified format.

    Returns
    -------
    dict or None
        basis, auxbasis, ecp
    """
    regex = r'\$([^\$]*)'
    regex = re.compile(regex)
    with open(filepath,'r') as F: 
        txt = F.read()
    keywords = dict()
    kwblocks = regex.findall(txt)
    for block in kwblocks: 
        lines = [line for line in block.split('\n') if line]
        kw_line = lines[0]
        keyword = kw_line.split(' ')[0]
        if keyword == 'end':
            continue
        keywords[keyword] = dict()
        for line in lines[1:]:
            sym,basis = line.split(' ',1)
            keywords[keyword][sym.lower()] = basis.strip()
    basis = keywords.get('basis',None)
    auxbasis = keywords.get('auxbasis',None)
    ecp = keywords.get('ecp',None)
    return basis,auxbasis,ecp

def create_tmol_input_opt_from_args(base_folder,args):
    folder_name = ''
    folder = f'{base_folder}/{folder_name}'
    
    kwargs = dict()
    kwargs['functional'] = args.tmfunctional
    if args.tmbasis is not None: 
        kwargs['basis'] = args.tmbasis
    else: 
        basis,auxbasis,ecp = read_basiset_file(args.tmbasisfile)
        kwargs['basis'] = basis
        kwargs['auxbasis'] = auxbasis
        kwargs['ecp'] = ecp
        
    kwargs['dispersion'] = args.tmdispersion
    kwargs['charge'] = 0
    kwargs['multiplicity'] = 1
    kwargs['grid'] = args.tmgrid
    kwargs['epsilon'] = args.tmepsilon
    kwargs['maxcore'] = args.tmmaxcore
    kwargs['ricore'] = args.tmricore
    kwargs['cavity'] = args.tmcavity
    kwargs['title'] = ''
    tmol_input = TurbomoleInput(folder,**kwargs)
    tmol_input.generate()

def create_tmol_input_SP_from_args(base_folder,args):
    folder_name = ''
    folder = f'{base_folder}/{folder_name}'
    
    kwargs = dict()
    kwargs['functional'] = args.tmfunctionalsp
    if args.tmbasissp is not None: 
        kwargs['basis'] = args.tmbasissp
    else: 
        basis,auxbasis,ecp = read_basiset_file(args.tmbasisfilesp)
        kwargs['basis'] = basis
        kwargs['auxbasis'] = auxbasis
        kwargs['ecp'] = ecp
        
    kwargs['dispersion'] = args.tmdispersion
    kwargs['charge'] = 0
    kwargs['multiplicity'] = 1
    kwargs['grid'] = args.tmgrid
    kwargs['epsilon'] = args.tmepsilon
    kwargs['maxcore'] = args.tmmaxcore
    kwargs['ricore'] = args.tmricore
    kwargs['cavity'] = args.tmcavity
    kwargs['title'] = ''
    tmol_input = TurbomoleInput(folder,**kwargs)
    tmol_input.generate()