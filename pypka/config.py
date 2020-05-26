import os
from constants import *
from psutil import cpu_count
from pprint import pformat
from numpy import arange
from copy import copy

class Config:
    @classmethod
    def storeParams(cls, titration_obj, log, debug, parameters):
        cls.log = log
        cls.debug = debug
        cls.titration = titration_obj
        cls.pypka_params = PypKaConfig(log)
        cls.delphi_params = DelPhiConfig(log)
        cls.mc_params = MCConfig(log)
        cls.parallel_params = ParallelConfig()

        cls.filter_params(parameters)

        return cls.pypka_params, cls.delphi_params, cls.mc_params, cls.parallel_params

    @classmethod
    def loadParams(cls, parameters):
        pypka_params, delphi_params, mc_params, parallel_params = parameters
        cls.pypka_params = pypka_params
        cls.delphi_params = delphi_params
        cls.mc_params = mc_params
        cls.parallel_params = parallel_params

    @classmethod
    def getConfigObj(cls, name):
        if name in cls.__dict__:
            return cls.__dict__
        elif name in cls.pypka_params:
            return cls.pypka_params
        elif name in cls.delphi_params:
            return cls.delphi_params
        elif name in cls.mc_params:
            return cls.mc_params
        else:
            return None

    @classmethod
    def filter_params(cls, parameters):
        """Checks the validity of the input parameters

        Args:
            parameters (dict): input parameters
        """
        # Define all parameters
        for param_name, param_value in parameters.items():
            if param_name.startswith('sites') or param_name in IGNORED_PARAMS:
                continue
            config_obj = cls.getConfigObj(param_name)
            if config_obj:
                config_obj[param_name] = param_value
            else:
                info = '{} is not a valid parameter.'.format(param_name)
                cls.log.report_warning(info, stdout=True)

        # Check if mandatory parameters were defined
        for param_name in MANDATORY_PARAMS:
            config_obj = cls.getConfigObj(param_name)
            if not (config_obj and config_obj[param_name] is not None):
                cls.log.raise_required_param_error(param_name)

        cls.pypka_params.set_structure_extension()
        cls.pypka_params.set_ncpus()
        cls.pypka_params.set_radii_charges_paths()
        cls.pypka_params.readTermini()

        cls.mc_params.set_pH_values(parameters)

        if cls.pypka_params['structure_output']:
            cls.pypka_params.set_structure_output(cls.mc_params['pHmin'],
                                                   cls.mc_params['pHmax'])

        if cls.delphi_params['pbc_dim'] == 2:
            cls.delphi_params.set_nonlinear_params(cls.pypka_params, parameters)

        if 'lipid_definition' in parameters:
            cls.pypka_params.define_lipids(parameters['lipid_definition'])


class ParametersDict:
    input_conversion = {'structure': 'f_in',
                        'grid_fill': 'perfil',
                        'pbc_dimensions': 'pbc_dim',
                        'convergence': 'maxc',
                        'output': 'f_out',
                        'logfile': 'f_log',
                        'titration_output': 'f_prot_out',
                        'clean': 'clean_pdb'}

    def __init__(self, log):
        self.log = log
        self.input_special_conditions = {}
        self.input_type = {}

    def convert_param_name(self, name):
        if name in self.input_conversion.keys():
            return self.input_conversion[name]
        return name

    def check_param_type(self, param_name, param_value, param_type, msg=''):
        """Checks if param_value is of type param_type

        If param_value is not of param_type, it tries to convert into the desired type.
        An error is raised if the conversion is not possible.

        Args:
            param_name (str): name of the parameter
            param_value: value of the parameter
            param_type (type): desired type for the parameter
            msg (str): error help message

        Returns:
            param_type: param_value is returned with the type param_type
        """
        if param_type is bool and not isinstance(param_value, bool):
            param_value = param_value.lower()
            if param_value in ('false', 'no'):
                param_value = False
            elif param_value in ('true', 'yes'):
                param_value = True

        if not isinstance(param_value, param_type):
            try:
                param_value = param_type(param_value)
            except ValueError:
                self.log.raise_input_param_error(param_name, param_type, msg)
        return param_value

    def check_conditions(self, param_name, param_value):
        if param_name in self.input_type:
            param_type = self.input_type[param_name]
            param_value = self.check_param_type(param_name, param_value, param_type)
        else:
            info = 'parameter {} is not being checked for type. ' \
                    'Please warn the developement team.'.format(param_name)
            self.log.report_warning(info, stdout=True)
        if param_name in self.input_special_conditions:
            condition = self.input_special_conditions[param_name]
            if condition == '>0':
                if  param_value <= 0:
                    self.log.raise_input_param_error(param_name,
                                                     'greater than zero.', '')
            elif param_value not in condition:
                    self.log.raise_input_param_error(param_name,
                                                     'in {}'.format(condition), '')

        return param_value

    def __getitem__(self, name):
        name = self.convert_param_name(name)
        if name not in self.__dict__:
            raise Exception('{} not in {}'.format(name, self.name))
        return self.__dict__[name]

    def __setitem__(self, name, key):
        name = self.convert_param_name(name)
        key = self.check_conditions(name, key)

        self.__dict__[name] = key

    def __contains__(self, item):
        item = self.convert_param_name(item)
        return item in self.__dict__

    def __str__(self):
        self_dict = {}
        for key, value in self.__dict__.items():
            if key not in self.not_to_print and key not in ('name', 'not_to_print', 'log'):
                self_dict[key] = value
        return '# {} Parameters\n{}\n'.format(self.name, pformat(self_dict))


class PypKaConfig(ParametersDict):
    """Configuration parameters
    """
    def __init__(self, log):
        super().__init__(log)

        self.name = 'PypKa'
        self.not_to_print = ['tmpsites', 'temp', 'input_special_conditions', 'input_type', 'NTR_atoms', 'CTR_atoms', 'box']

        self.tmpsites    = {}
        self.pid         = os.getpid()
        self.debug       = False
        self.ncpus       = None
        self.temp        = 298

        self.CpHMD_mode = False

        # Paths
        self.file_dir   = os.path.dirname(os.path.abspath(__file__))
        self.script_dir = os.path.dirname(__file__)
        self.pdb2pqr    = "{0}/pdb2pqr/pdb2pqr.py".format(self.file_dir)
        self.userff     = "{0}/pdb2pqr/dat/GROMOS.DAT".format(self.file_dir)
        self.usernames  = "{0}/pdb2pqr/dat/GROMOS.names".format(self.file_dir)

        # File Naming
        self.f_in               = None
        self.f_in_extension     = None
        self.f_out              = None
        self.f_prot_out         = None
        self.f_structure_out    = None

        # Output File
        self.structure_output   = None
        self.f_structure_out_pH = None
        self.ff_structure_out   = None


        # Force Field
        self.f_crg = None
        self.f_siz = None
        self.ffID  = 'G54A7'
        self.ff_family = 'GROMOS'
        self.NTR_atoms = None
        self.CTR_atoms = None
        self.LIPIDS = {}

        # Preprocessing parameters
        self.ffinput    = 'GROMOS'
        self.clean_pdb  = True
        self.keep_ions  = False
        self.ser_thr_titration = True

        self.cutoff     = -1
        self.slice      = 0.05
        self.box        = []

        # Parameters Validity
        self.input_special_conditions = {
            'temp': '>0',
            'ffID': ('G54A7', 'CHARMM36m'),
            'ffinput': ('GROMOS', 'AMBER', 'CHARMM'),
            'ff_structure_out': ('gromos_cph', 'amber')
        }
        self.input_type = {
            'debug'             : bool,
            'ncpus'             : int,
            'temp'              : float,
            'f_in'              : str,
            'f_out'             : str,
            'f_prot_out'        : str,
            'f_structure_out'   : str,
            'f_structure_out_pH': float,
            'ff_structure_out'  : str,
            'structure_output'  : str,
            'ffID'              : str,
            'ff_family'         : str,
            'ffinput'           : str,
            'cutoff'            : float,
            'slice'             : float,
            'clean_pdb'         : bool,
            'keep_ions'         : bool,
            'ser_thr_titration' : bool,
            'f_crg'             : str,
            'f_siz'             : str,
            'box'               : list,
            'CpHMD_mode'        : bool
        }

    def set_structure_extension(self):
        structure = self['structure']
        f_in_parts = structure.split('.')
        if len(f_in_parts) <= 1:
            self.log.raise_input_param_error('structure',
                'a string containing a file extension.',
                'Ex: structure.pdb or structure.gro')

        extension = f_in_parts[-1].lower().replace('pqr', 'pdb')
        if extension not in ('gro', 'pdb'):
            self.log.raise_input_param_error('structure',
                'a string containing a valid file extension.',
                'Ex: structure.pdb or structure.gro or structure.pqr')

        self.f_in_extension = extension

    def setBox(self, box):
        self['box'] = box

    def set_ncpus(self):
        ncpus = self['ncpus']
        if ncpus == -1:
            self['ncpus'] = cpu_count(logical=False)
        elif ncpus < 1:
            self['ncpus'] = 1

    def set_radii_charges_paths(self):
        file_path = os.path.join(self['script_dir'], self['ffID'])
        if not self['f_crg']:
            self['f_crg'] = '{}/DataBaseT.crg'.format(file_path)
        if not self['f_siz']:
            self['f_siz'] = '{}/DataBaseT.siz'.format(file_path)

        ffID = Config.pypka_params['ffID'].lower()
        if 'charmm36m' in ffID:
            self['ff_family'] = 'CHARMM'
        elif 'g54a7' in ffID:
            self['ff_family'] = 'GROMOS'
        else:
            raise Exception('Forcefield {0} not supported'.format(ffID))


    def define_lipids(self, lipids):
        for lipid in lipids:
            resname = lipids[lipid]
            self.LIPIDS[lipid] = resname

    def set_structure_output(self, pHmin, pHmax):
        error_raise = False
        structure_output = self['structure_output'].split(',')
        msg = 'CLI Example: "structure_output": ("structure.pdb", 7, "amber")\n '\
              'API Example: structure_output = structure.pdb, 7, amber'
        if len(structure_output) != 3:
            error_msg = 'a tuple containing a filename, the desired pH value and force field naming scheme.'
            self.log.raise_input_param_error('structure_output', error_msg, msg)


        outfilename = structure_output[0].strip('()"\' ')
        pH =  structure_output[1].strip('()"\' ')
        ff_out = structure_output[2].strip('()"\' ').lower()

        pH = self.check_param_type('structure_output_pH', pH, float, msg)

        self['f_structure_out'] = outfilename
        self['f_structure_out_pH'] = pH
        self['ff_structure_out'] = ff_out

        if pH < pHmin or pH > pHmax:
            message = 'in range [pHmin, pHmax].'
            self.log.raise_input_param_error('structure_output', message, '')

    def readTermini(self):
        script_dir = self['script_dir']
        ffID = self['ffID']
        NTR_atoms = []
        ntr_fname = "{}/{}/sts/NTRtau1.st".format(script_dir, ffID)
        with open(ntr_fname) as f:
            for line in f:
                parts = line.split()
                if len(parts) > 1:
                    NTR_atoms.append(parts[1].strip())
        CTR_atoms = []
        ctr_fname = "{}/{}/sts/CTRtau1.st".format(script_dir, ffID)
        with open(ctr_fname) as f:
            for line in f:
                parts = line.split()
                if len(parts) > 1:
                    CTR_atoms.append(parts[1].strip())
        self.NTR_atoms = NTR_atoms
        self.CTR_atoms = CTR_atoms

    def redefine_f_in(self, new_f_in):
        self['f_in'] = new_f_in
        self.set_structure_extension()


class DelPhiConfig(ParametersDict):
    """DelPhi configuration parameters
    """
    def __init__(self, log):
        super().__init__(log)

        self.name = 'DelPhi'
        self.not_to_print = ['p_atpos', 'p_rad3', 'p_chrgv4', 'atinf', 'p_iatmed',
                             'delphimol', 'input_type', 'input_special_conditions', 'lookup_atoms']

        self.perfil     = 0.9
        self.gsize      = 81
        self.scaleP     = 1
        self.scaleM     = 4
        self.precision  = 'single'
        self.ionicstr   = 0.1
        self.bndcon     = 3
        self.maxc       = 0.01
        self.nlit       = 500
        self.nonit      = 0
        self.relfac     = 0.75
        self.relpar     = 0.75
        self.pbx        = False
        self.pby        = False
        self.epssol     = 80.0
        self.pbc_dim    = 0
        self.epsin      = 20.0

        self.lookup_atoms = None

        # -1 NanoShaper off
        # 0 connolly surface
        # 1 skin
        # 2 blobby
        # 3 mesh
        # 4 msms
        self.nanoshaper = -1

        self.p_atpos   = None
        self.p_rad3    = None
        self.p_chrgv4  = None
        self.atinf     = None
        self.p_iatmed  = None
        self.delphimol = None

        self.input_type = {
            'perfil'    : float,
            'gsize'     : int,
            'scaleP'    : float,
            'scaleM'    : float,
            'precision' : str,
            'ionicstr'  : float,
            'bndcon'    : int,
            'maxc'      : float,
            'nlit'      : int,
            'nonit'     : int,
            'relfac'    : float,
            'relpar'    : float,
            'pbx'       : bool,
            'pby'       : bool,
            'epssol'    : float,
            'pbc_dim'   : int,
            'epsin'     : float,
            'nanoshaper': int
        }

        self.input_special_conditions = {
            'scaleP'   : '>0',
            'scaleM'   : '>0',
            'maxc'     : '>0',
            'gsize'    : '>0',
            'perfil'   : '>0',
            'ionicstr' : '>0',
            'nlit'     : '>0',
            'epssol'   : '>0',
            'epsin'    : '>0',
            'bndcon'   : (1, 2, 3, 4),
            'precision': ('single', 'double'),
            'pbc_dim'  : (0, 2),
            'nanoshaper': (0, 1, 2, 3, 4, -1)
        }

    def set_nonlinear_params(self, pypka_config, input_params):
        if self['relfac'] != 0.2 and \
            'relfac' not in input_params:
            self['relfac'] = 0.2
        if self['nonit'] != 5 and \
            'nonit' not in input_params:
            self['nonit'] = 5
        if pypka_config['cutoff'] == -1 and \
            'cutoff' not in input_params:
            pypka_config['cutoff'] = 5

    def redefineScale(self):
        scaleP = (self['gsize'] - 1) / (Config.pypka_params['box'][0])
        scaleM = int(4 / scaleP + 0.5) * scaleP

        self['scaleP'] = scaleP
        self['scaleM'] = scaleM

    def store_run_params(self, delphimol):
        if delphimol != 'reload':
            self.delphimol = delphimol
        else:
            delphimol = self.delphimol
        self.p_atpos  = copy(delphimol.get_atpos())
        self.p_rad3   = copy(delphimol.get_rad3())
        self.p_chrgv4 = copy(delphimol.get_chrgv4())
        self.atinf    = copy(delphimol.get_atinf())
        self.p_iatmed = copy(delphimol.get_iatmed())


class MCConfig(ParametersDict):
    """Monte Carlo configuration parameters
    """
    def __init__(self, log):
        super().__init__(log)

        self.name = 'Monte Carlo'
        self.not_to_print = ['input_type', 'input_special_conditions']
        self.pHmin      = 0
        self.pHmax      = 14
        self.pHstep     = 0.25
        self.seed       = 1234567
        self.couple_min = 2.0
        self.mcsteps    = 200000
        self.eqsteps    = 1000

        self.pH_values = []

        self.input_special_conditions = {
            'pHstep' : '>0',
            'mcsteps': '>0',
            'eqsteps': '>0'
        }
        self.input_type = {
            'pHmin'     : float,
            'pHmax'     : float,
            'pHstep'    : float,
            'seed'      : int,
            'couple_min': float,
            'mcsteps'   : int,
            'eqsteps'   : int
        }

    def set_pH_values(self, parameters):
        if 'pH' in parameters:
            pH = parameters['pH']
        else:
            pH = '{}-{}'.format(self.pHmin, self.pHmax)

        pH_error_info = 'pH can be a single value or a range. '\
                        'As default pH is set to [0, 14]\n'\
                        'API Example: pH = 7 or pH = 5,8\n'\
                        'CLI Example: "pH": 7 or "pH": "5,18"'

        pH_parts = pH.split(',')
        if len(pH_parts) == 1:
            pH_parts = pH.split('-')
        if len(pH_parts) == 2:
            self.pHmin = self.check_param_type('pHmin', pH_parts[0], float)
            self.pHmax = self.check_param_type('pHmax', pH_parts[1], float)

            diff = self.pHmax - self.pHmin
            if diff < 0:
                self.log.raise_input_param_error('pH', 'correctly defined.',
                                                 pH_error_info)
            if diff % self.pHstep != 0:
                self.pHmax += self.pHstep
            self.pH_values = arange(self.pHmin, self.pHmax + 0.001, self.pHstep)
            self.pHmax = float(self.pH_values[-1])
        elif len(pH_parts) == 1:
            pH = self.check_param_type('pH', pH, float)
            self.pH_values = (pH)
            self.pHmin = pH
            self.pHmax = pH
        else:
            self.log.raise_input_param_error('pH', 'correctly defined.',
                                             pH_error_info)


class ParallelConfig:
    def __init__(self):
        self.njobs      = None
        self.pb_time    = None
        self.total_jobs = None

        self.all_tautomers_order = None
        self.site_interactions = None
        self.interactions_look = None
        self.all_sites = None

        self.npossible_states = None
        self.possible_states_g = None
        self.possible_states_occ = None
        self.interactions = None
