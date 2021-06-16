import os
import subprocess
import re
from collections import namedtuple
from pathlib import Path

# ensure a proper import of pybel for openbabel v2 and v3
try:
    import pybel
except ImportError: 
    from openbabel import pybel


class TurbomoleOutput(object): 
    def __init__(self,folder,T=298.15,scale_factor=1,qh_thermo=False,
                 fcutoffs=(100,150)):
        self.folder = Path(folder)
        self.calc_type = self.guess_calculation_type()
        self.calc_state = self.guess_calculation_state()
        self.cosmo = self.uses_cosmo()
        self.geometry = self.parse_geometry()
        self.energy = self.parse_energy()
        self.freeh_log = ''
        self.T = T # Temperature in K
        self.scale_factor = scale_factor # Scaling factor for frequencies 
        self.qh_thermo = qh_thermo       # If the qh free energy is to be calculated 
        self.fcutoffs = fcutoffs         # Cutoffs for frequencies for qh in cm^-1
        self.frequencies = self.parse_frequencies()
        self.rotational_constants = self.parse_rotational_constants()
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
            return self._guess_calculation_state_optfreq()
        else:
            return self._guess_calculation_state_sp()
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
        file = self.folder / 'out.cosmo'
        return file.exists()
    def parse_thermo(self):
        if self.qh_thermo:
            self._parse_thermo_qh()
        else:
            self._parse_thermo_qh()
    def _parse_thermo_qh(self):
        cwd = os.path.abspath(os.getcwd())
        os.chdir(self.folder)
        cmd = ['./freeh','>','temp.temp']
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        for line in p.stdout:
            self.freeh_log += line
        p.wait()
        # Here do the parsing of temp.temp
        os.remove('temp.temp')
        os.chdir(cwd)
    def _parse_thermo_noqh(self):
        cwd = os.path.abspath(os.getcwd())
        os.chdir(self.folder)
        cutoff1,cutoff2 = self.fcutoffs
        with open('.thermorc','w') as F: 
            F.write(f'{cutoff1} {cutoff2} {self.T} {self.scale_factor}\n')
        cmd = ['./grimme','>','temp.temp']
        p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
        for line in p.stdout:
            self.freeh_log += line
        p.wait()
        # Here do the parsing of temp.temp
        os.remove('temp.temp')
        os.chdir(cwd)
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