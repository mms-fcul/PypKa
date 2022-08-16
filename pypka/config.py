import os
from copy import copy
from pprint import pformat
import logging
from numpy import arange
from psutil import cpu_count

from pypka.constants import *

logger = logging.getLogger(__name__)


class Config:
    @classmethod
    def storeParams(cls, titration_obj, debug, parameters, sites):
        cls.debug = debug
        cls.titration = titration_obj
        cls.pypka_params = PypKaConfig(sites)
        cls.delphi_params = DelPhiConfig()
        cls.mc_params = MCConfig()
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
            if (
                param_name.startswith("sites")
                or param_name.startswith("fixed_sites")
                or param_name in IGNORED_PARAMS
            ):
                continue
            config_obj = cls.getConfigObj(param_name)
            if config_obj:
                config_obj[param_name] = param_value
            else:
                info = "{} is not a valid parameter.".format(param_name)
                logger.warning(info)

        # Check if mandatory parameters were defined
        for param_name in MANDATORY_PARAMS:
            config_obj = cls.getConfigObj(param_name)
            if not (config_obj and config_obj[param_name] is not None):
                raise IOError(
                    'Required input parameter "{0}" '
                    "is not defined.".format(param_name)
                )

        cls.pypka_params.set_structure_extension()
        cls.pypka_params.set_ncpus()
        cls.pypka_params.set_radii_charges_paths()
        cls.pypka_params.readTermini()

        cls.mc_params.set_pH_values(parameters)

        if cls.pypka_params["structure_output"]:
            cls.pypka_params.set_structure_output(cls.mc_params["pH_values"])

        if cls.pypka_params["isoelectric_point"] and cls.pypka_params["sites"] != "all":
            warn_message = (
                "The isoelectric point can only be calculated when titrating all sites."
            )
            raise Exception(warn_message)

        if cls.delphi_params["pbc_dim"] == 2:
            cls.delphi_params.set_nonlinear_params(cls.pypka_params, parameters)

        if "lipid_definition" in parameters:
            cls.pypka_params.define_lipids(parameters["lipid_definition"])


class ParametersDict:
    input_conversion = {
        "structure": "f_in",
        "grid_fill": "perfil",
        "pbc_dimensions": "pbc_dim",
        "convergence": "maxc",
        "output": "f_out",
        "logfile": "f_log",
        "titration_output": "f_prot_out",
        "clean": "clean_pdb",
    }

    def __init__(self):
        self.input_special_conditions = {}
        self.input_type = {}

    def convert_param_name(self, name):
        if name in self.input_conversion.keys():
            return self.input_conversion[name]
        return name

    def check_param_type(self, param_name, param_value, param_type, msg=""):
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
            if param_value in ("false", "no"):
                param_value = False
            elif param_value in ("true", "yes"):
                param_value = True

        if not isinstance(param_value, param_type):
            try:
                param_value = param_type(param_value)
            except ValueError:
                raise ValueError(
                    'Input parameter "{0}" is not {1}\n '
                    "{2}".format(param_name, param_type, msg)
                )
        return param_value

    def check_conditions(self, param_name, param_value):
        if param_name in self.input_type:
            param_type = self.input_type[param_name]
            param_value = self.check_param_type(param_name, param_value, param_type)
        else:
            logger.warning(
                "parameter {} is not being checked for type. "
                "Please warn the development team.".format(param_name)
            )

        if param_name in self.input_special_conditions:
            condition = self.input_special_conditions[param_name]
            if condition == ">0":
                if param_value <= 0:
                    raise ValueError(
                        'Input parameter "{0}" is not greater than zero.'.format(
                            param_name
                        )
                    )
            elif param_value not in condition:
                raise ValueError(
                    'Input parameter "{0}" is not in {1}'.format(param_name, condition)
                )

        return param_value

    def __getitem__(self, name):
        name = self.convert_param_name(name)
        if name not in self.__dict__:
            raise Exception("{} not in {}".format(name, self.name))
        return self.__dict__[name]

    def __setitem__(self, name, key):
        name = self.convert_param_name(name)
        key = self.check_conditions(name, key)

        self.__dict__[name] = key

    def __contains__(self, item):
        item = self.convert_param_name(item)
        return item in self.__dict__

    def __str__(self):
        self_dict = self.get_clean_params()
        return "# {} Parameters\n{}\n".format(self.name, pformat(self_dict))

    def get_clean_params(self):
        self_dict = {}
        for key, value in self.__dict__.items():
            if key not in self.not_to_print and key not in (
                "name",
                "not_to_print",
                "log",
            ):
                self_dict[key] = value
        return self_dict


class PypKaConfig(ParametersDict):
    """Configuration parameters"""

    def __init__(self, sites):
        super().__init__()

        self.name = "PypKa"
        self.not_to_print = [
            "tmpsites",
            "temp",
            "input_special_conditions",
            "input_type",
            "NTR_atoms",
            "CTR_atoms",
            "box",
        ]

        self.sites = sites
        self.tmpsites = {}
        self.pid = os.getpid()
        self.debug = False
        self.ncpus = None
        self.temp = 298

        self.CpHMD_mode = False

        # Paths
        self.file_dir = os.path.dirname(os.path.abspath(__file__))
        self.ffs_dir = os.path.dirname(__file__)

        # File Naming
        self.f_in = None
        self.f_in_extension = None
        self.f_out = None
        self.f_prot_out = None
        self.f_structure_out = None

        # Output File
        self.structure_output = None
        self.coupled_sites_output = None
        self.f_structure_out_pH = None
        self.ff_structure_out = None
        self.save_pdb = None
        self.save_mc_energies = None
        self.isoelectric_point = None

        # Force Field
        self.f_crg = None
        self.f_siz = None
        self.ffID = "G54A7"
        self.ff_family = "GROMOS"
        self.sts = "sts"
        self.NTR_atoms = None
        self.CTR_atoms = None
        self.LIPIDS = {}

        # Preprocessing parameters
        self.ffinput = "GROMOS"
        self.clean_pdb = True
        self.pdb2pqr_inputfile = "input_clean.pdb"
        self.pdb2pqr_h_opt = True
        self.remove_hs = True
        self.keep_ions = False
        self.ser_thr_titration = True

        self.cutoff = -1
        self.slice = 0.05
        self.box = []

        # Parameters Validity
        self.input_special_conditions = {
            "temp": ">0",
            "ffID": ("G54A7", "CHARMM36m"),
            "ffinput": ("GROMOS", "AMBER", "CHARMM"),
            "ff_structure_out": ("gromos_cph", "amber"),
        }
        self.input_type = {
            "debug": bool,
            "ncpus": int,
            "temp": float,
            "f_in": str,
            "f_out": str,
            "f_prot_out": str,
            "f_structure_out": str,
            "f_structure_out_pH": float,
            "ff_structure_out": str,
            "structure_output": str,
            "coupled_sites_output": str,
            "isoelectric_point": bool,
            "ffID": str,
            "ff_family": str,
            "ffinput": str,
            "cutoff": float,
            "slice": float,
            "clean_pdb": bool,
            "pdb2pqr_h_opt": bool,
            "remove_hs": bool,
            "keep_ions": bool,
            "ser_thr_titration": bool,
            "f_crg": str,
            "f_siz": str,
            "box": list,
            "CpHMD_mode": bool,
            "save_pdb": str,
            "save_mc_energies": str,
            "ffs_dir": str,
            "sts": str,
        }

    def set_structure_extension(self):
        structure = self["structure"]
        f_in_parts = structure.split(".")
        if len(f_in_parts) <= 1:
            raise ValueError(
                'Input parameter "structure" is not a string containing a file extension.\n '
                "Ex: structure.pdb or structure.gro or structure.pqr"
            )

        extension = f_in_parts[-1].lower().replace("pqr", "pdb")
        if extension not in ("gro", "pdb"):
            raise ValueError(
                'Input parameter "structure" is not a string containing a valid file extension.\n '
                "Ex: structure.pdb or structure.gro or structure.pqr"
            )

        self.f_in_extension = extension

    def setBox(self, box):
        self["box"] = box

    def set_ncpus(self):
        ncpus = self["ncpus"]
        if ncpus == -1:
            self["ncpus"] = cpu_count(logical=False)
        elif ncpus < 1:
            self["ncpus"] = 1

    def set_radii_charges_paths(self):
        file_path = os.path.join(self["ffs_dir"], self["ffID"])
        if not self["f_crg"]:
            self["f_crg"] = "{}/DataBaseT.crg".format(file_path)
        if not self["f_siz"]:
            self["f_siz"] = "{}/DataBaseT.siz".format(file_path)

        ffID = Config.pypka_params["ffID"].lower()
        if "charmm36m" in ffID:
            self["ff_family"] = "CHARMM"
            # self["ffinput"] = "CHARMM"
            self["ser_thr_titration"] = False
        elif "g54a7" in ffID:
            self["ff_family"] = "GROMOS"
        else:
            raise Exception("Forcefield {0} not supported".format(ffID))

        if not os.path.isfile(self["f_crg"]):
            raise Exception(
                "DelPhi charges database file {} does not exist.".format(self["f_crg"])
            )
        elif not os.path.isfile(self["f_siz"]):
            raise Exception(
                "DelPhi radii database file {} does not exist.".format(self["f_siz"])
            )

    def define_lipids(self, lipids):
        for lipid in lipids:
            resname = lipids[lipid]
            self.LIPIDS[lipid] = resname

    def set_structure_output(self, pH_values):
        structure_output = self["structure_output"].split(",")
        msg = (
            'CLI Example: "structure_output": ("structure.pdb", 7, "amber")\n '
            "API Example: structure_output = structure.pdb, 7, amber"
        )
        if len(structure_output) != 3:
            raise ValueError(
                'Input parameter "structure_output" is not a tuple containing a filename, '
                "the desired pH value and force field naming scheme.\n{0}".format(msg)
            )

        outfilename = structure_output[0].strip("()\"' ")
        pH = structure_output[1].strip("()\"' ")
        ff_out = structure_output[2].strip("()\"' ").lower()

        pH = self.check_param_type("structure_output_pH", pH, float, msg)

        self["f_structure_out"] = outfilename
        self["f_structure_out_pH"] = pH
        self["ff_structure_out"] = ff_out

        if pH not in pH_values:
            pH_values.append(pH)
            pH_values.sort()
            # raise Exception(
            #    f"pH value specified in structure_output ({pH}) must be compatible with pH values being calculated ({pH_values})."
            # )

    def readTermini(self):
        self.sts_path = "{0}/{1}/{2}/".format(self.ffs_dir, self.ffID, self.sts)

        NTR_atoms = []
        ntr_fname = "{}/NTRtau1.st".format(self.sts_path)
        with open(ntr_fname) as f:
            for line in f:
                parts = line.split()
                if len(parts) > 1:
                    NTR_atoms.append(parts[1].strip())
        CTR_atoms = []
        ctr_fname = "{}/CTRtau1.st".format(self.sts_path)
        with open(ctr_fname) as f:
            for line in f:
                parts = line.split()
                if len(parts) > 1:
                    CTR_atoms.append(parts[1].strip())
        self.NTR_atoms = NTR_atoms
        self.CTR_atoms = CTR_atoms

    def redefine_f_in(self, new_f_in):
        self["f_in"] = new_f_in
        self.set_structure_extension()


class DelPhiConfig(ParametersDict):
    """DelPhi configuration parameters"""

    def __init__(self):
        super().__init__()

        self.name = "DelPhi"
        self.not_to_print = [
            "p_atpos",
            "p_rad3",
            "p_chrgv4",
            "atinf",
            "p_iatmed",
            "delphimol",
            "input_type",
            "input_special_conditions",
            "lookup_atoms",
        ]

        self.perfil = 0.9
        self.gsize = 81
        self.scaleP = 1
        self.scaleM = 4
        self.precision = "single"
        self.ionicstr = 0.1
        self.bndcon = 3
        self.maxc = 0.01
        self.nlit = 500
        self.nonit = 0
        self.relfac = 0.75
        self.relpar = 0.75
        self.pbx = False
        self.pby = False
        self.epssol = 80.0
        self.pbc_dim = 0
        self.epsin = 20.0

        self.lookup_atoms = None

        # -1 NanoShaper off
        # 0 connolly surface
        # 1 skin
        # 2 blobby
        # 3 mesh
        # 4 msms
        self.nanoshaper = -1

        self.p_atpos = None
        self.p_rad3 = None
        self.p_chrgv4 = None
        self.atinf = None
        self.p_iatmed = None
        self.delphimol = None

        self.input_type = {
            "perfil": float,
            "gsize": int,
            "scaleP": float,
            "scaleM": float,
            "precision": str,
            "ionicstr": float,
            "bndcon": int,
            "maxc": float,
            "nlit": int,
            "nonit": int,
            "relfac": float,
            "relpar": float,
            "pbx": bool,
            "pby": bool,
            "epssol": float,
            "pbc_dim": int,
            "epsin": float,
            "nanoshaper": int,
        }

        self.input_special_conditions = {
            "scaleP": ">0",
            "scaleM": ">0",
            "maxc": ">0",
            "gsize": ">0",
            "perfil": ">0",
            "ionicstr": ">0",
            "nlit": ">0",
            "epssol": ">0",
            "epsin": ">0",
            "bndcon": (1, 2, 3, 4),
            "precision": ("single", "double"),
            "pbc_dim": (0, 2),
            "nanoshaper": (0, 1, 2, 3, 4, -1),
        }

    def set_nonlinear_params(self, pypka_config, input_params):
        if self["relfac"] != 0.2 and "relfac" not in input_params:
            self["relfac"] = 0.2
        if self["nonit"] != 5 and "nonit" not in input_params:
            self["nonit"] = 5
        if pypka_config["cutoff"] == -1 and "cutoff" not in input_params:
            pypka_config["cutoff"] = 5

    def redefineScale(self):
        scaleP = (self["gsize"] - 1) / (Config.pypka_params["box"][0])
        scaleM = int(4 / scaleP + 0.5) * scaleP

        self["scaleP"] = scaleP
        self["scaleM"] = scaleM

    def store_run_params(self, delphimol):
        if delphimol != "reload":
            self.delphimol = delphimol
        else:
            delphimol = self.delphimol
        self.p_atpos = copy(delphimol.get_atpos())
        self.p_rad3 = copy(delphimol.get_rad3())
        self.p_chrgv4 = copy(delphimol.get_chrgv4())
        self.atinf = copy(delphimol.get_atinf())
        self.p_iatmed = copy(delphimol.get_iatmed())


class MCConfig(ParametersDict):
    """Monte Carlo configuration parameters"""

    def __init__(self):
        super().__init__()

        self.name = "Monte Carlo"
        self.not_to_print = ["input_type", "input_special_conditions"]
        self.pHmin = 0
        self.pHmax = 14
        self.pHstep = 0.25
        self.seed = 1234567
        self.couple_min = 2.0
        self.mcsteps = 200000
        self.eqsteps = 1000

        self.pH_values = []

        self.input_special_conditions = {
            "pHstep": ">0",
            "mcsteps": ">0",
            "eqsteps": ">0",
        }
        self.input_type = {
            "pHmin": float,
            "pHmax": float,
            "pHstep": float,
            "seed": int,
            "couple_min": float,
            "mcsteps": int,
            "eqsteps": int,
        }

    def set_pH_values(self, parameters):
        if "pH" in parameters:
            pH = parameters["pH"]
        else:
            pH = "{}-{}".format(self.pHmin, self.pHmax)

        pH_error_info = (
            "pH can be a single value or a range. "
            "As default pH is set to [0, 14]\n"
            "API Example: pH = 7 or pH = 5,8\n"
            'CLI Example: "pH": 7 or "pH": "5,18"'
        )

        pH_parts = pH.split(",")
        if len(pH_parts) == 1:
            pH_parts = pH.split("-")
        if len(pH_parts) == 2:
            self.pHmin = self.check_param_type("pHmin", pH_parts[0], float)
            self.pHmax = self.check_param_type("pHmax", pH_parts[1], float)

            diff = self.pHmax - self.pHmin
            if diff < 0:
                raise ValueError(
                    'Input parameter "pH" is not correctly defined.\n '
                    "{2}".format(pH_error_info)
                )
            if diff % self.pHstep != 0:
                self.pHmax += self.pHstep

            pH_values = arange(self.pHmin, self.pHmax + 0.001, self.pHstep).tolist()
            ndecimals = str(self.pHstep)[::-1].find(".")
            self.pH_values = [round(pH, ndecimals + 1) for pH in pH_values]
            self.pHmax = float(self.pH_values[-1])
        elif len(pH_parts) == 1:
            pH = self.check_param_type("pH", pH, float)
            self.pH_values = [pH]
            self.pHmin = pH
            self.pHmax = pH
        else:
            raise ValueError(
                'Input parameter "pH" is not correctly defined.\n '
                "{2}".format(pH_error_info)
            )


class ParallelConfig:
    def __init__(self):
        self.njobs = None
        self.pb_time = None
        self.total_jobs = None

        self.all_tautomers_order = None
        self.site_interactions = None
        self.interactions_look = None
        self.all_sites = None

        self.npossible_states = None
        self.possible_states_g = None
        self.possible_states_occ = None
        self.interactions = None
