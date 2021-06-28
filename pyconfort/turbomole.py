import os
import subprocess
import shutil
import shlex
import re
from collections import namedtuple
from pathlib import Path

# ensure a proper import of pybel for openbabel v2 and v3
try:
    import pybel
except ImportError: 
    from openbabel import pybel


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
        if 'ridft ended normally' in lines[-1]:
            return 'converged'
        
        for line in lines: 
            if 'ridft did not converge!' in line: 
                return 'failed'
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
    def __init__(self):
        pass
        # geometry
        # control file 
        # basis file 
        # auxbasis file
        # coord file 
        # mos file ? 
    def write_files(self,folder=None):
        """
        Writes and creates (If it does not exist) the input files in the 
        specified folder. Otherwise it creates the files in the current working 
        directory. 

        Parameters
        ----------
        folder : str, optional
            A valid path to a folder. The folder may not exist, but its parents 
            should exist, by default creates the files in the current directory.
        """
        pass

    def write_coord_file(self): 
        pass
    
    def write_basis_file(self):
        pass

    def write_auxbasis_file(self):
        pass

    def coord_file(self):
        pass
    
    def mos_file(self):
        pass