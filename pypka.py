import argparse
from ctypes import *
from multiprocessing import Pool
import cProfile
import os

import time

import delphi1
import delphi2

__author__  = "Pedro Reis"
__version__ = "0.3"

__email__   = "pdreis@fc.ul.pt"
__status__  = "Development"


"""
features to be implemented in the future are marked as TODO
"""


parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=""" 
Object-Oriented Script to Calculate the pKint of each site as well as the pairwise energies
Requires:
DelPhi2Py module installation
Nanoshaper and cppSolver are optional libraries

Objects:
DelPhiParams stores the DelPhi input parameters like a .prm file

TitratingMolecule is the molecule which has more than one Site 

Site 

Tautomer



Example: 
python delphiT.py $run_n $cpun ../TMP_auxD1.sites ../TMP_delphi1_Add.pqr 
DataBaseT.crg DataBaseT.siz cent_z-corrected ${dimension} -prm $gsizeP 
$scaleP $epsin $epssol $ionicstr $bndcon $maxc $linit $gsizeM $scaleM 
-temp $temp --debug -boxsize $boxsizex10
""")
parser.add_argument('sysname', help='System name', action='store')
parser.add_argument('ncpus', help='Number of CPU cores to be used',
                    action='store')
parser.add_argument('sites', help='.sites file location',
                    action='store')
parser.add_argument('pqr', help='.pqr file location')
parser.add_argument('crg', help='charge(CRG) database file location')
parser.add_argument('siz', help='radii(SIZ) database')
parser.add_argument('cent', help='file with the center of the sites')
parser.add_argument('dimensions', help='0 for isotropic systems and '
                    '2 for anisotropic systems', type=int,
                    choices=[0, 2])
parser.add_argument('-prm', help='DelPhi input parameters: '
                    'gsizeP scaleP epsin epssol ionicstr '
                    'bndcon maxc linit gsizeM scaleM',
                    required=True, nargs=10)
parser.add_argument('-temp', help='temperature in Kelvin',
                    required=True)
parser.add_argument('-boxsize', help='size of the box in Angstrom, '
                    'should only be used with dimensions==2',
                    required=False)
parser.add_argument('--debug', help='activation of the debug mode '
                    'to print extra information', action='store_true')
args = parser.parse_args()

def checkParsedInput():
    """Checks the input arguments parsed are of the correct type"""
    def wrongType(argument, complaint, explanation):
        raise IOError('Input variable {0} is not {1}.\n '
                      '{2}.'.format(argument, complaint, explanation))
    
    sysname = args.sysname
    try:
        ncpus   = int(args.ncpus)
        if ncpus < 1:
            raise
    except:
        wrongType('$ncpus', 'an integer bigger than zero',
                  'The number of CPU threads to use has to '
                  'be an integer greater than zero')
    f_sites = args.sites
    f_pqr   = args.pqr
    f_crg   = args.crg
    f_siz   = args.siz
    f_cent  = args.cent
    ndim    = int(args.dimensions)
    prm_params = args.prm
    boxsize = args.boxsize
    try:
        temp    = float(args.temp)
        if temp < 0:
            raise
    except:
        wrongType('$temp', 'a float bigger than zero',
                  'Temperature (in Kelvin) has to be bigger than zero')
        

    kBoltz = 5.98435e-6 ; # e^2/(Angstrom*K)
    log10 = 2.302585092994046

    return (sysname, ncpus, f_sites, f_pqr, f_crg, f_siz, f_cent,
            ndim, prm_params, temp, kBoltz, log10, boxsize)

def startPoolProcesses(targetFunction, iterable_job_arguments_list,
                       ncpus, assign='distributed'):
    """Start a pool of ncpus processes that execute the targetFunction
    with the argument iterable_job_arguments_list

    Args:
      targetFunction (function): function to be execute in each process
      iterable_job_arguments_list (iterable): object enumerating 
        the arguments to be passed to the targetFunction.
        when using assing='ordered', it can not be a generator
      ncpus  (int): number of processes to launch
      assign (str): mode to split the jobs among processes
        Example: iterable_job_arguments_list = [0, 1, 2, 3, 4], ncpus = 2
        'distributed' p0 = [0, 2, 4] p1 = [1, 3]
        'ordered' p0 = [0, 1, 2] p1 = [3, 4]
    
    Ensures:
      unpacked_results (list): list of lists containing the returned 
        results from the executed jobs. the order is dependent 
        only on $assign and $ncpus but not on the execution time.
    """      
    pool = Pool(processes=ncpus)
    jobs = [[] for i in range(ncpus)]
    results = []
    i = -1
    if assign == 'distributed': 
        for tautomer in iterable_job_arguments_list:
            i += 1
            jobs[i % ncpus].append(i)
    elif assign == 'ordered':
        max_njobs = len(iterable_job_arguments_list) / ncpus
        ncores_extra_load = len(iterable_job_arguments_list) % ncpus
        core = 0
        core_jobs = 0
        extra_load = False
        for tautomer in iterable_job_arguments_list:
            i += 1
            
            if core_jobs == max_njobs:
                if ncores_extra_load > core and not extra_load:
                    core_jobs -= 1
                    extra_load = True
                else:
                    core_jobs = 0
                    core += 1
                    extra_load = False
            jobs[core].append(i)
            core_jobs += 1
    for job in jobs:        
        result = pool.apply_async(targetFunction, args=(job, ))
        results.append(result)

    pool.close()
    pool.join()

    unpacked_results = []
    for results_percore in results:
        unpacked_results.append(results_percore.get())

    return unpacked_results

def runDelPhiSims(job_list):
    """Run sets of simulations of tautomers included in job_list

    Args:
      job_list (list): list of tautomer numbers to be analysed.
        The use of numbers is needed since class objects 
        can not be pickled and therefore can not be parsed 
        in multiprocessing.

    Ensures:
      results (list)
        tauname (str):   the tautomer name
        sitenum (int):   site id number
        esolvM  (float): the Molecule solvation energy
        sitpotM (list):  the site potential in the Molecule
          the index corresponds to the atoms in the site
        esolvS  (float): the Site solvation energy
        sitpotS (list):  the site potential in the Site
          the index corresponds to the atoms in the site
    """
    if args.debug:
        # each process will run in its own directory
        os.system('mkdir -p core_{0}'.format('_'.join([str(i) for i in job_list])))
        os.chdir('core_{0}'.format('_'.join([str(i) for i in job_list])))

    results = []
    for tau_number in job_list:
        tauname, sitenum, esolvM, sitpotM, esolvS, sitpotS = calcPotential(tau_number)
        results.append([tauname, sitenum, esolvM, sitpotM, esolvS, sitpotS])

    return results


def calcPotential(taut):
    """Run two DelPhi simulations: one for the tautomer and 
    other for the same tautomer when all other sites are neutral.

    Args:
      taut (int): number of the tautomer

    Ensures:     
      tauname (str):   the tautomer name
      sitenum (int):   site id number
      esolvM  (float): the Molecule solvation energy
      sitpotM (list):  the site potential in the Molecule
        the index corresponds to the atoms in the site
      esolvS  (float): the Site solvation energy
      sitpotS (list):  the site potential in the Site
        the index corresponds to the atoms in the site
    """
    # get the tautomer object from the tautomer number
    # transformation needed due to multiprocessing
    taut = tit_mole.getTautomerNumber(taut)
    
    # Whole Molecule with zero charge except for the tautomer being evaluated
    e_solvationM, sitpotM = taut.CalcPotentialTitratingMolecule()

    if args.debug:
        print 'finished Whole Molecule', taut._name, e_solvationM
    
    # Single Site Tautomer
    e_solvationS, sitpotS = taut.CalcPotentialTautomer()

    if args.debug:
        print 'finished Tautomer', taut._name, e_solvationM

    return (taut._name, taut.getSiteResNumber(), e_solvationM,
            sitpotM, e_solvationS, sitpotS)
    
def runInteractionCalcs(job_list):
    """Run interaction calculations between 
    sites in job_list

    Args:
      job_list (list)
        interaction_number (int): interaction number to be analysed.
          The use of numbers is needed since class objects 
          can not be pickled and therefore can not be parsed 
          in multiprocessing.

    Ensures:
      results (list)
        interaction_energies (str): .dat formatted energy 
          interactions between all tautomers of two sites
    """
    results = []
    for interaction_number in job_list:
        interaction_energies = tit_mole.calcInteractionNumber(interaction_number)
        results.append(interaction_energies)

    return results


class DelPhiParams:
    """DelPhi input parameters (.prm file) 
    
    Dictionary Class
    Their attributes are accessed directly outside the class
    """
    def __init__(self, prm_params, f_crg, f_siz, ndim):
        # TODO #
        # extend number of input parameters to include all/most
        # DelPhi input parameters
        #######################################################
        self._prm_params = prm_params

        self.igrid_P   = int(prm_params[0])
        self.scale_P   = float(prm_params[1])
        self.igrid_M   = int(prm_params[8])
        self.scale_M   = float(prm_params[9])

        self.repsin  = float(prm_params[2])
        self.repsout = float(prm_params[3])
        self.conc    = float(prm_params[4])
        self.ibctyp  = int(prm_params[5])
        self.res2    = float(prm_params[6])
        self.nlit    = int(prm_params[7])

        self.acent  = []
        self.in_crg = f_crg
        self.in_siz = f_siz
        self.in_pdb = 'delphi_in_stmod.pdb'

        self.radprb = 1.4
        self.energy = ['s', 'c']
        self.site   = ['a', 'q', 'p']
        self.scale1 = self.scale_P
        self.in_crg_len = len(self.in_crg)
        self.in_siz_len = len(self.in_siz)
        self.in_pdb_len = len(self.in_pdb)

        self.rmaxdim = 999.999

        self.precision = 'single'        
        self.perfil  = 0
        self.isurftype = -1
        self.parallel = False

        # Defines precision
        if self.precision == 'double':
            self.float_type = c_double
        elif self.precision == 'single':
            self.float_type = c_float
        else:
            raise IOError('Unknown precision definition {0}. '
                          'It should be either "double" or "single"')

        
        # Membrane specific
        if ndim == 2:
            self.relfac = 0.20
            self.relpar = 0.750
            self.nonit  = 5
            self.fcrg   = True
            self.pbx    = True
            self.pby    = True
            self.pbz    = False # is ignored
        else:
            self.relfac = 0.0
            self.relpar = 0.0
            self.nonit  = 0
            self.fcrg   = False
            self.pbx    = False
            self.pby    = False
            self.pbz    = False # is ignored
        
    def __str__(self):
        """
        Outputs the parameters used
        """
        out = """
DelPhi Parameters

py_igrid_P   = {}
py_scale_P   = {}
py_igrid_M   = {}
py_scale_M   = {}

py_repsin  = {}
py_repsout = {}
py_conc    = {}
py_ibctyp  = {}
py_res2    = {}
py_nlit    = {}

py_acent  = {}
py_in_crg = {}
py_in_siz = {}
py_in_pdb = {}

py_radprb  = {}
py_energy  = {}
py_site    = {}
py_scale1  = {}
py_in_crg_len = {}
py_in_siz_len = {}
py_in_pdb_len = {}
        """.format(self.igrid_P, self.scale_P,
                   self.igrid_M, self.scale_M,
                   self.repsin, self.repsout, self.conc,
                   self.ibctyp, self.res2, self.nlit,
                   self.acent, self.in_crg, self.in_siz,
                   self.in_pdb, self.radprb, self.energy,
                   self.site, self.scale1,
                   self.in_crg_len, self.in_siz_len,
                   self.in_pdb_len)
        return out

def truncate(f, n):
    """Truncates a float f to n decimal places without rounding"""
    s = '{}'.format(f)
    if 'e' in s or 'E' in s:
        return '{0:.{1}f}'.format(f, n)
    i, p, d = s.partition('.')
    return float('.'.join([i, (d+'0'*n)[:n]]))

class TitratingMolecule:
    """Molecule with more than one titrable sites
    """
    def __init__(self):
        """
        # DelPhi Parameters
        self._delphi_refparams (list): input DelPhi parameters

        # Sites
        self._sites (dict): key is residue number 
                            value is site object instance

        self._sites_order (list): values are site object instances
                                  order respects .sites

        # Atoms
        self._natoms (int): number of atoms in the molecule
        self._atoms (dict): key is atom id number
                            value is atom name
        self._atoms_array_position (dict):
            key is atom number
            value is atom position in pqr
            Example: if the first atom in pqr is number 541
                     then self._atoms_array_position[541] = 1

        # Offset of the first site
        This is only used because the membrane part is not finished
        and uses a pqr input of the first site
        self._InputPQRoffSet (list): 
            x, y, z coordinates of the offset in the input pqr

        # Site Interactions
        self._site_interactions (list): 
            contains tuples with pairs of site objects
        """
        self._delphi_refparams = ''
        self._sites = {}
        self._sites_order = []
        self._natoms = -1
        self._atoms = {}        
        self._atoms_array_position = {}        
        self._InputPQRoffSet = []
        self._site_interactions = []

    # Set Methods
    def addAtom(self, aname, anumb, position):
        self._atoms[anumb] = aname
        self._atoms_array_position[anumb] = position

    def loadParams(self, params):
        self._delphi_refparams = params

    def resetDelPhiData(self, params):
        """Resets all DelPhi Input data structures"""
        # internal DelPhi DataStructures
        # defined as c_type arrays and
        # passed to DelPhi as pointers
        self._py_atpos   = params.float_type * 3 * self._natoms
        self._p_atpos    = self._py_atpos()
        self._py_i_atpos = addressof(self._p_atpos)
        
        self._py_rad3   = params.float_type * self._natoms
        self._p_rad3    = self._py_rad3()
        self._py_i_rad3 = addressof(self._p_rad3)
        
        self._py_chrgv4   = params.float_type * self._natoms
        self._p_chrgv4    = self._py_chrgv4()
        self._py_i_chrgv4 = addressof(self._p_chrgv4)
        
        self._py_atinf   = (c_char * 15  * self._natoms)()
        self._py_i_atinf = addressof(self._py_atinf)
        
        self._nmedia = 1
        self._nobject = 1
        self._len_medeps = self._nmedia + self._nobject
        self._py_medeps  = params.float_type * self._len_medeps
        self._p_medeps    = self._py_medeps()
        self._py_i_medeps = addressof(self._p_medeps)
        
        self._len_iatmed  = self._natoms + 1
        self._py_iatmed   = c_int * self._len_iatmed
        self._p_iatmed    = self._py_iatmed()
        self._py_i_iatmmed = addressof(self._p_iatmed)
        
        self._py_dataobject   = (c_char * 96 * self._nobject * 2)()
        self._py_i_dataobject = addressof(self._py_dataobject)
        
    # Get Methods
    def getArrayPosition(self, atom_id):
        """Returns the index of the atom in DelPhi position array
        Disclamer: used only in testing"""
        return self._atoms_array_position[atom_id]
    
    def getAtomsList(self):
        """Returns a list of atoms details 
        (id, instance, atom index in DelPhi data structures)
        sorted by atom id number"""
        atom_list = []
        for atom in sorted(self._atoms.keys()):
            atom_list.append((self._atoms[atom], atom,
                              self._atoms_array_position[atom]))
        return atom_list

    def getDelPhiParams(self):
        return self._delphi_refparams
    
    def getTautomerNumber(self, t_index):
        """Uncodes the tautomer index t_index needed for multiprocessing
        and returns the corresponding tautomer instance"""
        i = -1
        for site in self._sites_order:
            i += 1
            if i == t_index:
                return site._ref_tautomer           
            for tautomer in sorted(site._tautomers.values()):
                i += 1
                if i == t_index:
                    return tautomer

    def getTautomerInstance(self, tautname, site_resnum):
        """Return the tautomer instance named tautname 
        existent in the site of the residue number site_resnum"""
        site = self._sites[site_resnum]
        if tautname in site._tautomers.keys():
            return site._tautomers[tautname]
        elif tautname == site._ref_tautomer._name:
            return site._ref_tautomer
        raise Exception('Something is very wrong!!!')

    def getNAtoms(self):
        """Return number of atoms of the Site"""
        return self._natoms
    
    def getTautNAtoms(self, taut_name):
        """Return number of atoms of the tautomer named taut_name
        Disclamer: it was used for testing"""
        for site in self._sites.values():
            for tautomer in site.getTautomers():
                if tautomer._name == taut_name:
                    return tautomer._natoms

    def getPQROffset(self):
        return self._InputPQRoffSet

    # Iter Methods
    def iterAtoms(self):
        """Generator that iterates through all atoms details (name, id,
        position) in the Site"""
        for atom in sorted(self._atoms.keys()):
            yield (self._atoms[atom], atom,
                   self._atoms_array_position[atom])

    def iterAllSitesTautomers(self):
        """Generator that iterates through all Tautomer instances.
        The iteration is sorted by site and within each site the first
        to be yielded is the reference tautomer and then the rest of
        the tautomers by order"""
        for site in self._sites_order:
            yield site._ref_tautomer                
            for tautomer in sorted(site._tautomers.values()):
                yield tautomer

    def iterNonRefSitesTautomers(self):
        """Generator that iterates through all Tautomer instances except
        the reference ones"""
        for site in self._sites_order:
            for tautomer in sorted(site._tautomers.values()):
                yield tautomer

    # Printing Methods
    def printAllSites(self):
        """Prints all Site names"""
        for site in self._sites_order:
            print site.getName()

    def printAllTautomers(self):
        """Prints all Tautomer details"""
        for site in self._sites_order:
            print site.getName()
            for tautomer in site.iterTautomers():
                print tautomer
            print site._ref_tautomer
            
    # Input Files Manipulation Methods
    def loadSites(self, f_sites):
	"""Opens .sites file 
        Stores info about the sites names and numbers
        Fetches info about reference sites and charge sets of the
        different tautomers
        
        Args:
          f_sites (str): .sites file name
        """
	with open(f_sites) as f:
	    for line in f:
	        cols = line.split()
	        resnum = cols[0]
	        cols = line.split()
	        res_tauts = cols[1:]
                
                sID = Site(resnum, self)
                self._sites[resnum] = sID
                self._sites_order.append(sID)
                
                sID.setTautomers(res_tauts)

        for site in self._sites.values():            
            site.addReferenceTautomer()
            site.addChargeSets()

    def createRefPDB(self, f_pqr):
        """Transforms pqr into a pdb 
        where the residues are in their standard state.
	The standard state is the last tautomer is the site file
	because tau1 in the site file is TAU0 in the pdb
        """
	natoms = 0
	with open(f_pqr) as f:
	    new_pdb_content = ''
	    for line in f:
	        if line[0:6] == 'ATOM  ':
	            natoms += 1
	            resnumb = line[22:26]
	            aname = line[12:16].replace(' ', '')
                    anumb = int(line[6:11])
                    aposition = natoms - 1
                    self.addAtom(aname, anumb, aposition)
	            if resnumb in self._sites.keys():
                        # change res name to reference tautomer
	                ref_tau_name = self._sites[resnumb].getRefTautomerName()
	                new_pdb_content += line[:17] + ref_tau_name + line[20:54]+'\n'
                        # add atom to corresponding site
                        self._sites[resnumb].addAtom(aname, anumb)                                             
	            else:
	                new_pdb_content += line[:54]+'\n'

	with open('delphi_in_stmod.pdb', 'w') as f_new:
	    f_new.write(new_pdb_content)

        self._natoms = natoms
	
    def readCenters(self, f_cent, boxsize=False):
	"""Opens cent file 
        Stores centers of sites
	"""        
	with open(f_cent) as f:
            nsites = 0
	    for line in f:
                nsites += 1
	        cols = line.split()
                res_number = cols[0]
                x = float(cols[1]) * 10
                y = float(cols[2]) * 10
                z = float(cols[3]) * 10
                if nsites == 1:
                    x = truncate(x, 3)
                    y = truncate(y, 3)
                    z = truncate(z, 3)
                    self._InputPQRoffSet = [x, y, z]
                if res_number in self._sites:
                    if boxsize:
                        self._sites[res_number].addCenter([x, y, z], boxsize=boxsize)
                    else:
                        self._sites[res_number].addCenter([x, y, z])
                else:
                    raise Exception('Something is wrong... Site in Centers '
                                    'not defined in cent file')

    def readDelPhiInputFiles(self):
        """Calls delphi routine to read PDB, CRG, SIZ
        Stores info about:
            self._p_atpos:  (x, y, z) position
            self._p_rad3:   atom radius
            self._p_chrgv4: atom charge
            self._py_atinf: atom type, residue info
        """
        params = self._delphi_refparams
        self.resetDelPhiData(params)

        py_rmaxdim = delphi1.delphi(params.igrid_P, params.scale_P,
                                    params.repsin, params.repsout,
                                    params.radprb, params.conc,
                                    params.ibctyp, params.res2,
                                    params.nlit, params.acent,
                                    params.energy, params.site,
                                    params.in_pdb, params.in_crg,
                                    params.in_siz, params.in_pdb_len,
                                    params.in_crg_len,
                                    params.in_siz_len, self._natoms,
                                    self._nobject, self._py_i_atpos,
                                    self._py_i_rad3,
                                    self._py_i_chrgv4,
                                    self._py_i_atinf,
                                    self._py_i_medeps,
                                    self._py_i_iatmmed,
                                    self._py_i_dataobject,
                                    params.rmaxdim)
        ### TODO ###
        # The py_rmaxdim return is only made so it is easy to
        # incorporate the perfil calculation in the future
        #
        # py_rmaxdim can be used to calculate igrid from perfil
        #if py_perfil != 0:
        #    py_igrid = int(py_scale * 100 / py_perfil * py_rmaxdim)
        #    if py_igrid % 2 == 0:
        #        py_igrid += 1
        ############################################################

        if args.debug:
            print '    x        y        z     radius  charge       atinf'
            for i in range(self._natoms):
                print '{0:8.3f} {1:8.3f} {2:8.3f} {3:7.3f} {4:7.3f} ' \
                    '"{5}"'.format(self._p_atpos[i][0], self._p_atpos[i][1],
                                   self._p_atpos[i][2], self._p_rad3[i], self._p_chrgv4[i],
                                   self._py_atinf[i].value)
    
    def calcSiteInteractionsNonParallel(self):
        """Calculates the pairwise interaction energies 
        and writes them in a formatted .dat file
        All interactions are calculated in serial
        """
        self.writeDatHeader()
        counter = 0
        to_write = ''
        iterAtomsList = self.getAtomsList()
        for site1 in self._sites_order[:-1]:
            counter += 1
            for site2 in self._sites_order[counter:]:
                for taut1 in site1.iterOrderedTautomers():
                    site_atom_list = taut1._site.getAtomNumbersList()
                    for taut2 in site2.iterOrderedTautomers():                        
                        if taut1 == taut2:
                            interaction = 0.0
                        else:
                            interaction = taut1.calcInteractionWith(taut2,
                                                                    site_atom_list,
                                                                    iterAtomsList)
        
                        to_write += self.convertIntoDatFormat(taut1, taut2, interaction)
        with open('{0}.dat'.format(sysname), 'a') as f_new:
            f_new.write(to_write)

    def calcSiteInteractionsParallel(self, ncpus):
        """Calculates the pairwise interaction energies 
        and writes them in a formatted .dat file
        Interactions are calculated using a pool of processes

        Args:
          ncpus (int): number of cpus to be used
        """        
        self.writeDatHeader()
        counter = 0
        for site1 in self._sites_order[:-1]:
            counter += 1
            for site2 in self._sites_order[counter:]:
                self._site_interactions.append((site1, site2))

        ncpus = min(len(self._site_interactions), ncpus)
        results = startPoolProcesses(runInteractionCalcs,
                                     self._site_interactions, ncpus, assign='ordered')

        with open('{0}.dat'.format(sysname), 'a') as f_new:
            to_write = ''
            for formatted_interactions in results:
                if not formatted_interactions: # if the core was not assigned anything
                    continue
                for formatted_interaction in formatted_interactions:
                    to_write += formatted_interaction

            f_new.write(to_write)

    def calcInteractionNumber(self, inumber):
        """Calculates the pairwise interaction energies
        related to the two sites specified in inumber
        
        Args:
          inumber (int): site interaction pair code

        Ensures: 
          to_write (str): site interaction energies formatted
        to be written in .dat file
        """
        to_write = ''
        sites = self._site_interactions[inumber]
        site1 = sites[0]
        site2 = sites[1]

        iterAtomsList = self.getAtomsList()

        for taut1 in site1.iterOrderedTautomers():
            site_atom_list = taut1._site.getAtomNumbersList()
            for taut2 in site2.iterOrderedTautomers():                        
                if taut1 == taut2:
                    interaction = 0.0
                else:
                    interaction = taut1.calcInteractionWith(taut2,
                                                            site_atom_list,
                                                            iterAtomsList)

                to_write += self.convertIntoDatFormat(taut1, taut2, interaction)
        
        return to_write

    def convertIntoDatFormat(self, tau1, tau2, interaction):
        """Returns a .dat format interaction line
        
        Args:
          tau1 (Tautomer): first tautomer of the pair
          tau2 (Tautomer): second tautomer of the pair
          interaction (float): interaction energy of the tautomers pair

        Ensures:
          line (str): .dat formatted interaction line
        """
        site1_index = self._sites_order.index(tau1._site)
        tau1_index = tau1._name[-1]
        tau1_dat = '{0:>3} {1:>2}'.format(site1_index, tau1_index)
        site2_index = self._sites_order.index(tau2._site)
        tau2_index = tau2._name[-1]
        tau2_dat = '{0:>3} {1:>2}'.format(site2_index, tau2_index)
        line = '{0}   {1}   {2:13.6e}\n'.format(tau1_dat, tau2_dat, interaction)
        return line

    def writeDatHeader(self):
        """Writes pKint energies in .dat file header"""
        to_write = '{0}\n'.format(len(self._sites))
        for site in self._sites_order:        
            to_write += '{0:3s}-{1:<7}{2:>2}  P  *\n'.format(site._res_name,
                                                             site._res_number,
                                                             len(site._tautomers) + 1)
            if site._type == 'c':
                tau_prot_state = 0
                ref_prot_state = 1
            elif site._type == 'a':
                tau_prot_state = 1
                ref_prot_state = 0
                
            for tautomer in site.iterOrderedTautomersWithoutRef():
                to_write += '{0:1d} {1:13.6e}\n'.format(tau_prot_state, tautomer._dg)
            to_write += '{0:1d}  0.000000e+00\n'.format(ref_prot_state)
        with open('{0}.dat'.format(sysname), 'w') as f_new:
            f_new.write(to_write)


    def calcpKint(self, unpacked_results):
        """Calculation the pKint of all tautomers
        """
        i = -1
        if args.debug:        
            print '############ results ############'                        
        for tautomer in tit_mole.iterAllSitesTautomers():
            i += 1
            core_index = i % ncpus
            job_index = i / ncpus
            result = unpacked_results[core_index][job_index]
        
            tautname    = result[0]
            tautresnumb = result[1]
        
            esolvationM = result[2]
            sitpotM     = result[3]
            esolvationS = result[4]
            sitpotS     = result[5]
                
            tautomer = tit_mole.getTautomerInstance(tautname, tautresnumb)
        
            
            tautomer.saveDelPhiResults(esolvationS, sitpotS, esolvationM,
                                       sitpotM)
            if args.debug:
                print '### new tautomer ###'
                print (i, core_index, job_index, tautname,
                       tautomer._name, esolvationM, esolvationS)
                print (tautomer._name, tautomer._esolvationS,
                       len(tautomer._sitpotS), tautomer._esolvationM,
                       len(tautomer._sitpotM))

            tautomer.calcBackEnergy()
            if not tautomer.isRefTautomer():
                tautomer.calcpKint()
                if args.debug:
                    print ('pkint', tautomer._name, tautomer._dg,
                           id(tautomer))


class Site:
    """Titrable Site with more than one Tautomer"""
    def __init__(self, res_number, molecule):
        """
        Args: 
          res_number (int): number of the site in the .pdb 
                            and .sites file
          molecule (TitratingMolecule): molecule to
                                        which the Site belongs

        # Site Details
        self._res_name (str): residue name
        self._type (str): site type ('c' or 'a')
                          'c' if site is cationic
                          'a' if site is anionic
        self._center (list) : Site geometric center
            if dimension is 0, then center is defined in cent file
            elif dimension is 2, then center (z) is defined in cent file
                                 and center (x, y) are defined in boxsize variable
        self._center_original (list): Site geometric center on 
                                      cent file when dimension is 2

        # Site Tautomers
        self._tautomers (dict): tautomers of Site (except reference tautomer)
                                key is tautomer name
                                value is tautomer instance
        self._ref_tautomer (Tautomer): reference tautomer instance

        # Site Atoms
        self._atoms (dict): existing atoms in site
                            key is atom id name
                            value is atom id number
        """
        self._molecule = molecule
        self._res_number = res_number
        self._res_name = ''
        self._atoms  = {}        
        self._tautomers = {}
        self._ref_tautomer = ''
        self._center = []
        self._type = ''

    # Set Methods
    def setTautomers(self, res_tauts):
        """Adds Tautomers from res_tauts to self._tautomers

        Args:
            res_tauts (list): tautomers from sites file
        """
        for tautomer in res_tauts:
            self._res_name = tautomer[:3]
            correct_number = str(int(tautomer[-1]) - 1)
            tautomer = tautomer[:2] + correct_number

            tID = Tautomer(tautomer, self, self._molecule)
            self._tautomers[tautomer] = tID

    def addReferenceTautomer(self):
        """Gets last tautomer from .sites file adds one and saves it as the
        reference tautomer"""
        last_tautomer = max(sorted(self._tautomers.keys()))
        correct_number = str(int(last_tautomer[-1]) + 1)
        ref_tautomer = last_tautomer[:2] + correct_number
        tID = Tautomer(ref_tautomer, self, self._molecule)
        self._ref_tautomer = tID


    def addChargeSets(self):
        """Stores the charge set of each existing tautomer 
        for the present Site"""
        for tautomer in self._tautomers.values():
            tautomer.loadChargeSet(self._res_name,
                                   self._ref_tautomer)
            

    def addAtom(self, aname, anumb):
        """
        Args:
            aname (str): atom name
            anumb (int): atom id number
        """
        self._atoms[aname] = anumb

    def addCenter(self, center, boxsize=False):        
        x = center[0]
        y = center[1]
        z = center[2]
        if boxsize:
            if ndim != 2:
                raise Exception('ERROR: The original center is only '
                                'needed for 2-dimensional calculation')
            self._center_original = [x, y, z]
            x = float(boxsize) / 2
            y = float(boxsize) / 2        
        self._center = [x, y, z]

    def setType(self, stype):
        self._type = stype

    # Get Methods
    def getName(self):
        return self._res_name

    def getTautomers(self):
        """Returns list of all tautomers instances 
        except the tautomers of reference"""
        return self._tautomers.values()

    def getAtomNumbersList(self):
        return self._atoms.values()

    def getRefTautomerName(self):
        return self._ref_tautomer.getName()

    def getRefTautomer(self):
        return self._ref_tautomer

    def getMolecule(self):
        return self._molecule

    def getType(self):
        return self._type
        
    def getRefProtState(self):
        if self._type == 'c':
            reftau_prot_state = 1
        elif self._type == 'a':
            reftau_prot_state = -1
        return reftau_prot_state

    def getCenterOriginal(self):
        if ndim != 2:
            raise Exception('ERROR: The original center is only '
                            'needed for 2-dimensional calculation')
        return self._center_original

    def getCenter(self):
        return self._center

    def getResNumber(self):
        return self._res_number

    # Iter Methods
    def iterTautomers(self):
        for i in self._tautomers.values():
            yield i

    def iterOrderedTautomersWithoutRef(self):
        for i in sorted(self._tautomers.keys()):
            yield self._tautomers[i]

    def iterOrderedTautomers(self):
        for i in sorted(self._tautomers.keys()):
            yield self._tautomers[i]
        yield self._ref_tautomer

    # Print Methods
    def __str__(self):
        tautomers = ' '.join(sorted(self._tautomers.keys()))
        natoms = len(self._atoms)
        center = [ round(i, 2) for i in self._center ]
        out = 'Site Number -> {0:5s}   '\
              'Tautomers -> {1:30}  '\
              'Reference -> {2:5} '\
              'NAtoms -> {3:6} '\
              'Center -> {4}'.format(self._res_number,tautomers,
                                     self._ref_tautomer.getName(),
                                     natoms, center)
        return out


class Tautomer:
    """Tautomers share the same site and atoms
    Tautomers have different charge sets for the same atoms
    """
    def __init__(self, name, site, molecule):
        """
        Args:
            name (str): name of the Tautomer
            site (Site): belonging site
            molecule (TitratingMolecule): belonging molecule

        # Inheritance
        self._molecule
        self._site

        # Tautomer Details
        self._name
        self._charge_set (dict): charge set of the tautomer
                                 key is atom name str
                                 value is charge float
        self._natoms (int): number of atoms of the Tautomer
        Redundant to self._site._natoms as it must be the same 
        between all Tautomers belonging to the same Site.

        # Tautomer Energy Details
        self._esolvation (float): solvation energy
        self._e_back (float): background interaction energy
        self._dg (float): pKint energy
        """
        self._molecule = molecule
        self._site = site
        self._name = name
        self._esolvation = ''
        self._charge_set = {}
        self._natoms = 0
        self._e_back = 0.0
        self._dg = 0.0

    # Set Methods
    def loadChargeSet(self, res_name, ref_tautomer):
        """Reads .st file related to Tautomer with residue name res_name
        Stores charge set related to both the Tautomer and the
        Reference Tautomer

        .st file named TYRtau1.st has charge set of TY0 and reference TY2
                       TYRtau2.st has charge set of TY1 and reference TY2
        """
        tau_number = int(self._name[-1]) + 1
        fname = '{0}tau{1}.st'.format(res_name, tau_number)
        with open(fname) as f:
            nline = 0
            charge_set1 = {}
            charge_set2 = {}
            floats_1 = []
            floats_2 = []
            for line in f:
                nline += 1
                cols = line.split()
                if nline > 1:
                    atom_name = cols[1]
                    charge_set1[atom_name] = float(cols[-2])
                    charge_set2[atom_name] = float(cols[-1])
                else:
                    self._pKmod = float(line)                    
                    
        if sum(charge_set1.values()) < 0.001 and sum(charge_set2.values()) > -1.001:
            if args.debug:
                print fname, 'anionic'
            self._charge_set = charge_set1
            ref_tautomer._charge_set = charge_set2
            self._site.setType('a')            
        elif sum(charge_set1.values()) > 0.99 and sum(charge_set2.values()) < 0.001:
            if args.debug:
                print fname, 'cationic'
            self._charge_set = charge_set2
            ref_tautomer._charge_set = charge_set1
            self._site.setType('c')
        self._natoms = len(self._charge_set)
        ref_tautomer._natoms = len(self._charge_set)
        if args.debug:
            print self._charge_set
            print ref_tautomer._charge_set
            print 'finished reading', fname

    def saveDelPhiResults(self, esolvationS, sitpotS, esolvationM, sitpotM):
        self._esolvationS  =  esolvationS
        self._sitpotS      =  sitpotS
        self._esolvationM  =  esolvationM
        self._sitpotM      =  sitpotM

    def defineFocusingPrepParams(self, params, turn_on=False):
        """Set focusing parameters for a DelPhi calculations.  
        Focusing can be turned on or shutdown by turn_on trigger variable
        
        Args:
            params (DelPhiParams): object containing all DelPhi parameters
            turn_on (boolean): manages the focusing parameters
                               if turn_on, then focusing parameters are on
                               if not turn_on, then focusing parameters are off

        """
        molecule = self._molecule
        self._esolvation   = 999.999
        if turn_on:
            self._py_in_frc    = 'self'            
            
            self._py_sitpot    = params.float_type * molecule._natoms
            self._p_sitpot     = self._py_sitpot()
            self._py_i_sitpot  = addressof(self._p_sitpot)
            
            self._py_out_phi   = True
            
            self._len_phimap   = params.igrid_P * params.igrid_P * params.igrid_P
            self._py_phimap4   = c_float * self._len_phimap
            self._p_phimap4    = self._py_phimap4()
            self._py_i_phimap4 = addressof(self._p_phimap4)
        else:
            self._py_in_frc    = ''            
            
            self._py_sitpot    = params.float_type * self._natoms
            self._p_sitpot     = self._py_sitpot()
            self._py_i_sitpot  = addressof(self._p_sitpot)
            
            self._py_out_phi   = False
            
            self._len_phimap   = 0
            self._py_phimap4   = c_float * self._len_phimap
            self._p_phimap4    = self._py_phimap4()
            self._py_i_phimap4 = addressof(self._p_phimap4)

    def setCommonDetails(self, params):
        """Set DelPhi parameters common to all DelPhi runs"""
        molecule = self._molecule
	self._nmedia = molecule._nmedia
	self._nobject = molecule._nobject
	self._len_medeps = self._nmedia + self._nobject
	self._py_medeps  = params.float_type * self._len_medeps
	self._p_medeps    = self._py_medeps()
	self._py_i_medeps = addressof(self._p_medeps)
        memmove(self._py_i_medeps, molecule._py_i_medeps, sizeof(params.float_type) * self._len_medeps)

	self._len_iatmed  = molecule._natoms + 1
	self._py_iatmed   = c_int * self._len_iatmed
	self._p_iatmed    = self._py_iatmed()
	self._py_i_iatmmed = addressof(self._p_iatmed)
        memmove(self._py_i_iatmmed, molecule._py_i_iatmmed, sizeof(c_int) * self._len_iatmed)

	self._py_dataobject   = (c_char * 96 * self._nobject * 2)()
	self._py_i_dataobject = addressof(self._py_dataobject)
        memmove(self._py_i_dataobject, molecule._py_i_dataobject, sizeof(c_char) * 96 * self._nobject * 2)
        
    def setDetailsFromWholeMolecule(self, params):
        """Set DelPhi parameters to run a calculation of a whole molecule 
        (all sites neutral, except one)"""
        molecule = self._molecule
        self._py_atpos   = params.float_type * 3 * molecule._natoms
	self._p_atpos    = self._py_atpos()
	self._py_i_atpos = addressof(self._p_atpos)

        memmove(self._py_i_atpos, molecule._py_i_atpos, sizeof(params.float_type) * 3 * molecule._natoms)

	self._py_rad3   = params.float_type * molecule._natoms
	self._p_rad3    = self._py_rad3()
	self._py_i_rad3 = addressof(self._p_rad3)

        memmove(self._py_i_rad3, molecule._py_i_rad3, sizeof(params.float_type) * molecule._natoms)

	self._py_chrgv4   = params.float_type * molecule._natoms
	self._p_chrgv4    = self._py_chrgv4()
	self._py_i_chrgv4 = addressof(self._p_chrgv4)

        memmove(self._py_i_chrgv4, molecule._py_i_chrgv4, sizeof(params.float_type) * molecule._natoms)

	self._py_atinf   = (c_char * 15  * molecule._natoms)()
	self._py_i_atinf = addressof(self._py_atinf)

        memmove(self._py_i_atinf, molecule._py_i_atinf, sizeof(c_char) * 15 * molecule._natoms)

        if ndim == 2:                
            half_box_xy = self._site.getCenter()[0]
            site_center = self._site.getCenterOriginal()
            pqr_offset = self._site._molecule.getPQROffset()

            # TODO: has to be corrected after the input file change from TMP_delphi1_Add
            offset_x = -(half_box_xy - pqr_offset[0]) + (half_box_xy - site_center[0])
            offset_y = -(half_box_xy - pqr_offset[1]) + (half_box_xy - site_center[1])

        for atom_name, atom_id, atom_position in molecule.iterAtoms():
            # quick fix, should be done only once per site
            if ndim == 2:                
                self._p_atpos[atom_position][0] += offset_x
                self._p_atpos[atom_position][1] += offset_y
            if atom_id in self._site.getAtomNumbersList():
                self._p_chrgv4[atom_position] = self.getCharge(atom_name)

            else:
                self._p_chrgv4[atom_position] = 0.0

            if args.debug:
                print self._name, 'starting'
                print '    x        y        z     radius  charge       atinf'
        
                print '{0:8.3f} {1:8.3f} {2:8.3f} {3:7.3f} {4:7.3f} ' \
                    '"{5}"'.format(self._p_atpos[atom_position][0], self._p_atpos[atom_position][1],
                                   self._p_atpos[atom_position][2], self._p_rad3[atom_position], self._p_chrgv4[atom_position],
                                   self._py_atinf[atom_position].value)

    def setDetailsFromTautomer(self, params):
        """Set DelPhi parameters to run a calculation 
        of a single site tautomer"""
        molecule = self._molecule
        self._py_atpos   = params.float_type * 3 * self._natoms
	self._p_atpos    = self._py_atpos()
	self._py_i_atpos = addressof(self._p_atpos)
        
	self._py_rad3   = params.float_type * self._natoms
	self._p_rad3    = self._py_rad3()
	self._py_i_rad3 = addressof(self._p_rad3)

	self._py_chrgv4   = params.float_type * self._natoms
	self._p_chrgv4    = self._py_chrgv4()
	self._py_i_chrgv4 = addressof(self._p_chrgv4)

	self._py_atinf   = (c_char * 15  * self._natoms)()
	self._py_i_atinf = addressof(self._py_atinf)

        if ndim == 2:
            half_box_xy = self._site.getCenter()[0]
            site_center = self._site.getCenterOriginal()
            pqr_offset = self._site._molecule.getPQROffset()

            offset_x = -(half_box_xy - pqr_offset[0]) + (half_box_xy - site_center[0])
            offset_y = -(half_box_xy - pqr_offset[1]) + (half_box_xy - site_center[1])

        site_atom_position = -1
        for atom_name, atom_id, atom_position in molecule.iterAtoms():
            if atom_id in self._site.getAtomNumbersList():
                site_atom_position += 1
                self._p_atpos[site_atom_position] = molecule._p_atpos[atom_position]
                self._p_rad3[site_atom_position] = molecule._p_rad3[atom_position]
                self._p_chrgv4[site_atom_position] = self.getCharge(atom_name)
                self._py_atinf[site_atom_position].value = molecule._py_atinf[atom_position].value
                # quick fix, should be being done only once per site
                if ndim == 2:                
                    self._p_atpos[site_atom_position][0] += offset_x
                    self._p_atpos[site_atom_position][1] += offset_y

                if args.debug:
                    print (self._p_atpos[site_atom_position][0],
                           self._p_atpos[site_atom_position][1],
                           self._p_atpos[site_atom_position][2],
                           self._p_rad3[site_atom_position],
                           self._p_chrgv4[site_atom_position],
                           self._py_atinf[site_atom_position].value)
    
    # Get Methods
    def getSiteResNumber(self):
        return self._site.getResNumber()

    def getCharge(self, atom_name):
        return self._charge_set[atom_name]

    def getName(self):
        return self._name

    # Print Methods
    def __str__(self):
        out = self._name + '\n'
        for i in self._charge_set.keys():
            out += '{0:>7.3f} {1}\n'.format(self._charge_set[i], i)
        return out

    # Assertion Methods
    def isRefTautomer(self):
        if self == self._site._ref_tautomer:
            return True
        else:
            False

    # Calculation Methods
    def CalcPotentialTautomer(self):
        """Run DelPhi simulation of single site tautomer
        
        Ensures:
            self._esolvation (float): tautomer solvation energy
            self._p_sitpot (list): potential on site atoms
        """
        if args.debug:
            start = time.clock()
            print self._name
            print self._charge_set
        molecule = self._molecule
        params = molecule.getDelPhiParams()
        self.setCommonDetails(params)
        self.setDetailsFromTautomer(params)
        self.defineFocusingPrepParams(params)

        if args.debug:
            t0 = time.clock() - start
        self._esolvation = delphi2.delphi(params.igrid_M,
                                          params.scale_M,
                                          params.repsin,
                                          params.repsout,
                                          params.radprb, params.conc,
                                          params.ibctyp, params.res2,
                                          500, self._site._center,
                                          params.energy, params.site,
                                          0, 0.0,
                                          0.0, False,
                                          False,
                                          params.in_pdb,
                                          params.in_crg,
                                          params.in_siz,
                                          self._py_in_frc,
                                          params.in_pdb_len,
                                          params.in_crg_len,
                                          params.in_siz_len,
                                          self._natoms,
                                          molecule._nmedia,
                                          molecule._nobject,
                                          self._py_i_atpos,
                                          self._py_i_rad3,
                                          self._py_i_chrgv4,
                                          self._py_i_atinf,
                                          self._py_i_medeps,
                                          self._py_i_iatmmed,
                                          self._py_i_dataobject,
                                          self._py_i_phimap4,
                                          params.scale1,
                                          self._py_out_phi,
                                          self._py_i_sitpot,
                                          self._esolvation,
                                          params.isurftype,
                                          params.parallel)
        if args.debug:
            t1 = time.clock() - start
            filename = '{0}_{1}.profl'.format(self._name, self.getSiteResNumber())
            with open(filename, 'a') as f_new:
                f_new.write('time -> {0:10} {1:10}\n'.format(t0, t1))
                f_new.write('{} {} {} {} {}\n'
                            '{} {} {} {} {}\n'
                            '{} {} {} {} {}\n'
                            '{} {} {} {} {}\n'
                            '{} {} {} {} {}\n'
                            '{} {} {} {} {}\n'
                            '{} {} {} {} {}\n'
                            '{} {} {} {}\n'.format(params.igrid_M, params.scale_M,
                                                   params.repsin, params.repsout,
                                                   params.radprb, params.conc,
                                                   params.ibctyp, params.res2, 500,
                                                   self._site._center, params.energy,
                                                   params.site, 0, 0.0, 0.0, False,
                                                   False, '', '',
                                                   '', self._py_in_frc,
                                                   0,
                                                   0,
                                                   0, self._natoms,
                                                   molecule._nmedia, molecule._nobject,
                                                   self._py_i_atpos, self._py_i_rad3,
                                                   self._py_i_chrgv4, self._py_i_atinf,
                                                   self._py_i_medeps,
                                                   self._py_i_iatmmed,
                                                   self._py_i_dataobject,
                                                   self._py_i_phimap4, params.scale1,
                                                   self._py_out_phi, self._py_i_sitpot,
                                                   self._esolvation))
            
            with open('tmp_site', 'w') as f_new:
                text = str(self._esolvation) + str(self._name) + '\n'
                for i in range(self._natoms):            
                    text += '{0:8.3f} {1:8.3f} {2:8.3f} {3:7.3f} {4:7.3f} ' \
                            '"{5}" {6:10.3f}\n'.format(self._p_atpos[i][0],
                                                       self._p_atpos[i][1],
                                                       self._p_atpos[i][2],
                                                       self._p_rad3[i],
                                                       self._p_chrgv4[i],
                                                       self._py_atinf[i].value,
                                                       self._p_sitpot[i])
                f_new.write(text)

        return self._esolvation, self._p_sitpot[:]

    def CalcPotentialTitratingMolecule(self):
        """Run DelPhi simulation of the site tautomer 
        within the whole molecule

        Ensures:
            self._esolvation (float): tautomer solvation energy
            self._p_sitpot (list): potential on site atoms
        """
        if args.debug:
            start = time.clock()
        molecule = self._molecule
        params = molecule.getDelPhiParams()
        self.setCommonDetails(params)
        self.setDetailsFromWholeMolecule(params)
        self.defineFocusingPrepParams(params, turn_on=True)

        if args.debug:
            t0 = time.clock() - start
        self._esolvation = delphi2.delphi(params.igrid_P,
                                          params.scale_P,
                                          params.repsin,
                                          params.repsout,
                                          params.radprb, params.conc,
                                          params.ibctyp, params.res2,
                                          params.nlit, self._site._center,
                                          params.energy, params.site,
                                          params.nonit, params.relfac,
                                          params.relpar, params.pbx,
                                          params.pby,
                                          params.in_pdb,
                                          params.in_crg,
                                          params.in_siz,
                                          self._py_in_frc,
                                          params.in_pdb_len,
                                          params.in_crg_len,
                                          params.in_siz_len,
                                          molecule._natoms,
                                          molecule._nmedia,
                                          molecule._nobject,
                                          self._py_i_atpos,
                                          self._py_i_rad3,
                                          self._py_i_chrgv4,
                                          self._py_i_atinf,
                                          self._py_i_medeps,
                                          self._py_i_iatmmed,
                                          self._py_i_dataobject,
                                          self._py_i_phimap4,
                                          params.scale1,
                                          self._py_out_phi,
                                          self._py_i_sitpot,
                                          self._esolvation,
                                          params.isurftype,
                                          params.parallel)

        if args.debug:
            t1 = time.clock() - start
        self._py_out_phi = False
        
        self._esolvation = delphi2.delphi(params.igrid_P,
                                          params.scale_M,
                                          params.repsin,
                                          params.repsout,
                                          params.radprb, params.conc,
                                          3, params.res2,
                                          500, self._site._center,
                                          params.energy, params.site,
                                          0, 0.0,
                                          0.0, False,
                                          False,
                                          params.in_pdb,
                                          params.in_crg,
                                          params.in_siz,
                                          self._py_in_frc,
                                          params.in_pdb_len,
                                          params.in_crg_len,
                                          params.in_siz_len,
                                          molecule._natoms,
                                          molecule._nmedia,
                                          molecule._nobject,
                                          self._py_i_atpos,
                                          self._py_i_rad3,
                                          self._py_i_chrgv4,
                                          self._py_i_atinf,
                                          self._py_i_medeps,
                                          self._py_i_iatmmed,
                                          self._py_i_dataobject,
                                          self._py_i_phimap4,
                                          params.scale1,
                                          self._py_out_phi,
                                          self._py_i_sitpot,
                                          self._esolvation,
                                          params.isurftype,
                                          params.parallel)
        if args.debug:
            t2 = time.clock() - start
            filename = '{0}_{1}.profl'.format(self._name, self.getSiteResNumber())
            with open(filename, 'a') as f_new:
                f_new.write('time -> {0:10}     {1:10}    {2:10}\n'.format(t0, t1, t2))
                f_new.write('{} {} {} {} {}\n'
                            '{} {} {} {} {}\n'
                            '{} {} {} {} {}\n'
                            '{} {} {} {} {}\n'
                            '{} {} {} {} {}\n'
                            '{} {} {} {} {}\n'
                            '{} {} {} {} {}\n'
                            '{} {} {} {}\n\n'.format(params.igrid_P,
                                              params.scale_P,
                                              params.repsin,
                                              params.repsout,
                                              params.radprb, params.conc,
                                              params.ibctyp, params.res2,
                                              params.nlit, self._site._center,
                                              params.energy, params.site,
                                              params.nonit, params.relfac,
                                              params.relpar, params.pbx,
                                              params.pby,
                                              params.in_pdb,
                                              params.in_crg,
                                              params.in_siz,
                                              self._py_in_frc,
                                              params.in_pdb_len,
                                              params.in_crg_len,
                                              params.in_siz_len,
                                              molecule._natoms,
                                              molecule._nmedia,
                                              molecule._nobject,
                                              self._py_i_atpos,
                                              self._py_i_rad3,
                                              self._py_i_chrgv4,
                                              self._py_i_atinf,
                                              self._py_i_medeps,
                                              self._py_i_iatmmed,
                                              self._py_i_dataobject,
                                              self._py_i_phimap4,
                                              params.scale1,
                                              self._py_out_phi,
                                              self._py_i_sitpot,
                                              self._esolvation,
                                              params.isurftype,
                                              params.parallel))
                f_new.write('{} {} {} {} {}\n'
                            '{} {} {} {} {}\n'
                            '{} {} {} {} {}\n'
                            '{} {} {} {} {}\n'
                            '{} {} {} {} {}\n'
                            '{} {} {} {} {}\n'
                            '{} {} {} {} {}\n'
                            '{} {} {} {}\n'.format(params.igrid_P,
                                              params.scale_M,
                                              params.repsin,
                                              params.repsout,
                                              params.radprb, params.conc,
                                              3, params.res2,
                                              500, self._site._center,
                                              params.energy, params.site,
                                              0, 0.0,
                                              0.0, False,
                                              False,
                                              params.in_pdb,
                                              params.in_crg,
                                              params.in_siz,
                                              self._py_in_frc,
                                              params.in_pdb_len,
                                              params.in_crg_len,
                                              params.in_siz_len,
                                              molecule._natoms,
                                              molecule._nmedia,
                                              molecule._nobject,
                                              self._py_i_atpos,
                                              self._py_i_rad3,
                                              self._py_i_chrgv4,
                                              self._py_i_atinf,
                                              self._py_i_medeps,
                                              self._py_i_iatmmed,
                                              self._py_i_dataobject,
                                              self._py_i_phimap4,
                                              params.scale1,
                                              self._py_out_phi,
                                              self._py_i_sitpot,
                                              self._esolvation))
            
                                
            print self._esolvation, self._name
            print '    x        y        z     radius  charge       atinf'
            with open('tmp', 'w') as f_new:
                text = str(self._esolvation) + str(self._name) + '\n'
                for i in range(molecule.getNAtoms()):            
                    text += '{0:8.3f} {1:8.3f} {2:8.3f} {3:7.3f} {4:7.3f} ' \
                            '"{5}" {6:10.3f}\n'.format(self._p_atpos[i][0],
                                                       self._p_atpos[i][1],
                                                       self._p_atpos[i][2],
                                                       self._p_rad3[i],
                                                       self._p_chrgv4[i],
                                                       self._py_atinf[i].value,
                                                       self._p_sitpot[i])
                f_new.write(text)

        return self._esolvation, self._p_sitpot[:]

    def calcBackEnergy(self):
        """Calculates background energy contribution"""
        if args.debug:
            print self._name, 'background energy start'
        molecule = self._molecule
        for atom_name, atom_id, atom_position in molecule.iterAtoms():
            if atom_id not in self._site.getAtomNumbersList():
                if args.debug:
                    print self._sitpotM[atom_position]
                self._e_back += molecule._p_chrgv4[atom_position] * self._sitpotM[atom_position]
        if args.debug:
            print self._name, self._e_back
            print 'e_back finished'

    def calcpKint(self):
        """Calculates the pKint of the tautomer"""
        ref_taut = self._site._ref_tautomer
        
        dG_solvationM = ref_taut._esolvationM - self._esolvationM
        dG_solvationS = ref_taut._esolvationS - self._esolvationS
        dG_back       = ref_taut._e_back - self._e_back
        
        dG_solvationM /= log10
        dG_solvationS /= log10 
        dG_back       /= log10

        if self._site.getType() == 'a':
            pKint = self._pKmod + (dG_solvationM - dG_solvationS + dG_back)
            chargediff = -1
        elif self._site.getType() == 'c':
            pKint = self._pKmod - (dG_solvationM - dG_solvationS + dG_back)
            chargediff = 1            
        else:
            raise Exception('Site files were poorly interpreted')

        dg = pKint * log10 * kBoltz * temp * chargediff
        if args.debug:
            print 'pKint ->', self._name, pKint, dg
        self._dg = dg

    def calcInteractionWith(self, tautomer2, site_atom_list, iterAtomsList):
        """Calculates the interaction energy 
        between self tautomer and tautomer2

        Args:
            tautomer2 (Tautomer)
            site_atom_list (list): atom numbers belonging to the site
            iterAtomsList (list): details on the titrating molecule atoms
        """
        molecule = self._molecule
        interaction = 0.0
        tau2_ref = tautomer2._site._ref_tautomer
        for atom_name, atom_id, atom_position in iterAtomsList:
            if atom_id in site_atom_list:
                charge_ref = molecule._p_chrgv4[atom_position]
                charge_tau = self.getCharge(atom_name)
                charge = charge_ref - charge_tau
                
                potential_ref = tau2_ref._sitpotM[atom_position]
                potential_tau2 = tautomer2._sitpotM[atom_position]
                potential = potential_ref - potential_tau2
                
                interaction += charge * potential

        site1_chrgtyp = self._site.getRefProtState()
        site2_chrgtyp = tautomer2._site.getRefProtState()

        dG_interaction = site1_chrgtyp * site2_chrgtyp
        dG_interaction *= abs(interaction * kBoltz * temp)
        
        return dG_interaction


if __name__ == '__main__':
    # Declaration of global variables
    (sysname, ncpus, f_sites, f_pqr, f_crg, f_siz, f_cent, ndim,
     prm_params, temp, kBoltz, log10, boxsize) = checkParsedInput()

    # Creating instance of TitratingMolecule
    tit_mole = TitratingMolecule()

    # Reading .sites and .st files
    tit_mole.loadSites(f_sites)

    # Creating .pdb from input .pqr
    # In .pdb residues are in their standard state
    tit_mole.createRefPDB(f_pqr)

    # Reading cent file
    # Dealing with offset of the .pqr input file in ndim==2
    if ndim == 0:
        tit_mole.readCenters(f_cent)
    elif ndim == 2:
        tit_mole.readCenters(f_cent, boxsize=boxsize)

    # Calling original DelPhi parameters
    delphi_params = DelPhiParams(prm_params, f_crg, f_siz, ndim)
    
    if args.debug:
        tit_mole.printAllSites()    
        tit_mole.printAllTautomers()
        print delphi_params

    # Loads original delphi_params as TitratingMolecule attributes
    tit_mole.loadParams(delphi_params)

    # Creates DelPhi data structures by
    # read DelPhi input files (.pdb, .crg, .siz)
    tit_mole.readDelPhiInputFiles()

    # Runs DelPhi simulations for all tautomers
    results = startPoolProcesses(runDelPhiSims,
                                 tit_mole.iterAllSitesTautomers(),
                                 ncpus)
    # Calculates the pKint of all tautomers
    tit_mole.calcpKint(results)

    # Non parallel version testing
    # tit_mole.calcSiteInteractionsNonParallel()
    
    # Profiling effort
    # cProfile.runctx('tit_mole.calcSiteInteractionsNonParallel()', globals(), locals(), 'profile-calcSiteInteraction.out')

    # Calculates sites interaction energies and write .dat file
    tit_mole.calcSiteInteractionsParallel(ncpus)
    
    print 'exited successfully'
