#! /usr/bin/python

"""
A python API and CLI to perform pKa calculations on peptides,
proteins or lipid bilayers.
"""

import config
from cli import checkParsedInput, readSettings, inputParametersFilter
from cleaning import inputPDBCheck, cleanPDB
from formats import convertTermini, pdb2gro
import log

from delphi4py.delphi4py import DelPhi

from molecule import Molecule
from concurrency import startPoolProcesses, runDelPhiSims
import os

__author__ = "Pedro Reis"
__version__ = "0.4"

__email__ = "pdreis@fc.ul.pt"
__status__ = "Development"


def getTitrableSites(pdb):
    """Gets the all titrable sites from the pdb

    Returns an input structure to the Titration class containing all
    titrable sites found in the pdb file. 

    Args:
        pdb (str): The filename a PDB file
    
    Returns: 
        A dict mapping all titrable sites found in the pdb file

    """
    config.f_in = pdb
    config.params['ffID'] = 'G54A7'

    tit_mol = Molecule()

    tit_mol.makeSimpleSites()

    sites = tit_mol.getSitesOrdered()
    sites_keys = []
    for site in sites:
        sitename = site.getName()
        sitenumber = site.getResNumber()
        sitenumber = convertTermini(sitenumber)
        if sitename in ('NTR', 'CTR'):
            if sitename == 'NTR':
                sites_keys.append('{0}N'.format(sitenumber))
            else:
                sites_keys.append('{0}C'.format(sitenumber))
        else:
            sites_keys.append(str(sitenumber))

    return sites_keys


class Titration(object):
    """Main pypka class

    Serves as a wrapper to all other classes

    Args:
        parameters
        sites
        debug
        datfile

    Attributes:
        _pKas: A dict with the calculated pKa values.

    """
    def __init__(self, parameters, sites='all', debug=None, datfile=None):
        """
        self._pKas = {}
        """
        self._pKas = {}

        if sites:  # None if from CLI
            parameters['sites_A'] = sites
            if sites != 'all':
                config.sites = sites
            else:
                config.sites = {}

        if datfile:
            config.f_dat = datfile
            parameters = readSettings(datfile) 

        config.debug = debug

        # Check Input Variables Validity
        inputParametersFilter(parameters)

        print('Start Preprocessing')
        self.preprocessing()
        self.processDelPhiParams()

        print('Start PB Calculations')

        self.DelPhiLaunch()
        print('API exited successfully')

    def preprocessing(self):
        # Creating instance of TitratingMolecule
        config.tit_mole = Molecule()

        if config.f_in_extension == 'gro':
            groname = config.f_in
            if len(config.sites) > 0:
                chains_length, chains_res = inputPDBCheck(config.f_in, config.sites)
                config.tit_mole.loadSites(chains_length, chains_res)
            else:
                log.inputVariableError('sites',
                                       'defined.',
                                       'When using a .gro file format input is used, '
                                       'sites needs to be defined')
        if config.f_in_extension == 'pdb':
            # Reading .st files
            # If the titrable residues are defined
            if len(config.sites) > 0:
                chains_length, chains_res = inputPDBCheck(config.f_in,
                                                          config.sites)
                config.tit_mole.loadSites(chains_length, chains_res)
            # If the titrable residues are not defined and
            # the input pdb file is incomplete
            elif config.params['clean_pdb']:
                chains_res, sites = config.tit_mole.makeSimpleSites()
            # If the titrable residues are not defined and
            # the input pdb is not missing any atoms
            else:
                chains_res, sites = config.tit_mole.makeSimpleSites()
                config.tit_mole.deleteAllSites()
                config.tit_mole.makeSites(sites=sites)

            # Creates a .pdb input for DelPhi
            # where residues are in their standard state
            if config.params['clean_pdb']:
                inputpqr = 'clean.pqr'
                outputpqr = 'cleaned_tau.pqr'
                sites = config.tit_mole.getSites()
                site_numb_n_ref = {}
                for site in sites:
                    site_numb_n_ref[site] = sites[site].getRefTautomerName()
                
                cleanPDB(config.f_in, config.pdb2pqr, chains_res,
                         inputpqr, outputpqr, site_numb_n_ref)
                config.tit_mole.deleteAllSites()
                config.tit_mole.makeSites(useTMPgro=True, sites=list(sites.keys()))
            else:
                pdb2gro(config.f_in, 'TMP.gro', config.tit_mole.box,
                        config.sites)
            groname = 'TMP.gro'

        config.tit_mole.readGROFile(groname)
        if not config.debug and config.f_in_extension == 'pdb' and os.path.isfile('TMP.gro'):
            os.remove('TMP.gro')

    def processDelPhiParams(self):
        # Storing DelPhi parameters and Creates DelPhi data structures
        if config.debug:
            logfile = 'LOG_readFiles'
        else:
            logfile = config.f_log

        delphimol = DelPhi(config.f_crg, config.f_siz, 'delphi_in_stmod.pdb',
                           config.tit_mole.getNAtoms(),
                           config.params['gsize'],
                           config.params['scaleM'],
                           config.params['precision'],
                           epsin=config.params['epsin'],
                           epsout=config.params['epssol'],
                           conc=config.params['ionicstr'],
                           ibctyp=config.params['bndcon'],
                           res2=config.params['maxc'],
                           nlit=config.params['nlit'],
                           nonit=config.params['nonit'],
                           relfac=0.0,
                           relpar=0.0,
                           pbx=config.params['pbx'],
                           pby=config.params['pby'],
                           debug=config.debug, outputfile=logfile)
        if not config.debug:
            os.remove('delphi_in_stmod.pdb')

        log.checkDelPhiErrors(logfile, 'readFiles')
        
        # Loads delphi4py object as TitratingMolecule attributes
        config.tit_mole.loadDelPhiParams(delphimol)

        if config.debug:
            config.tit_mole.printAllSites()
            config.tit_mole.printAllTautomers()
            print(delphimol)

    def DelPhiLaunch(self):
        if len(config.tit_mole.getSites()) < 1:
            raise Exception('At least one site has to be correctly defined.')

        # Runs DelPhi simulations for all tautomers
        results = startPoolProcesses(runDelPhiSims,
                                     config.tit_mole.iterAllSitesTautomers(),
                                     config.params['ncpus'])

        # Calculates the pKint of all tautomers
        config.tit_mole.calcpKint(results)

        # Calculates sites interaction energies and write .dat file
        config.tit_mole.calcSiteInteractionsParallel(config.params['ncpus'])

        #  Monte Carlo sampling
        pKas, pmeans = config.tit_mole.runMC()

        sites = config.tit_mole.getSitesOrdered()

        c = -1
        for i in pKas:
            c += 1
            site = sites[c]._res_number
            pK = i[0]
            if pK != 100.0:
                self._pKas[site] = pK
            else:
                self._pKas[site] = '-'
        self._pmeans = pmeans

    def getAverageProt(self, site, pH):
        """Calculates the average protonation of a site at a given pH value

        Args:
            site (str): Residue number of the site with the suffix 'N'
            or 'C' to indicate a N-ter or C-ter, respectively.

            pH (float): pH value 
        
        Returns:  
            A float of the average protonation of the site at the
            selected pH value
        """
        pKa = self[site]
        if isinstance(pKa, str):
            return 'pk Not In Range'
        average_prot = 10 ** (pKa - pH) / (1 + 10 ** (pKa - pH))
        return average_prot

    def getProtState(self, site, pH):
        """Indicates the most probable protonation state of a site at a given
        pH value

        Args:
            site (str): Residue number of the site with the suffix 'N'
            or 'C' to indicate a N-ter or C-ter, respectively.

            pH (float): pH value 

        Returns: 

            A tuple with two elements. Example: ('protonated', 0.9)

            The first element is a string indicating the most probable
            state of the site. It can be either 'protonated',
            'deprotonated' or 'undefined'. The 'undefined' state is
            prompted if the average protonation state is between 0.1 and 0.9.
          
            The second element is a float of the average protonation of the site

        """
        state = 'undefined'
        average_prot = self.getAverageProt(site, pH)

        if isinstance(average_prot, str):
            return state, average_prot

        if average_prot > 0.9:
            state = 'protonated'
        elif average_prot < 0.1:
            state = 'deprotonated'

        return state, average_prot

    def __iter__(self):
        self._iterpKas = []
        for site in config.tit_mole.getSitesOrdered():
            self._iterpKas.append(site.getResNumber())
        self._iternumb = -1
        self._itermax = len(self._iterpKas)
        return self

    def __next__(self):
        self._iternumb += 1
        if self._iternumb < self._itermax:
            site = self._iterpKas[self._iternumb]
            if site > config.terminal_offset:
                if site - config.terminal_offset == config.tit_mole._NTR:
                    site = 'NTR'
                elif site - config.terminal_offset == config.tit_mole._CTR:
                    site = 'CTR'
                else:
                    raise Exception('Something is terribly wrong')
            return site
        else:
            raise StopIteration

    def correct_site_numb(self, numb):
        if numb == 'NTR':
            numb = config.tit_mole._NTR + config.terminal_offset
        elif numb == 'CTR':
            numb = config.tit_mole._CTR + config.terminal_offset
        if isinstance(numb, str):
            try:
                numb = int(numb)
            except ValueError:
                raise Exception('Unknown site')
        return numb

    def getParameters(self):
        """Get the parameters used in the calculations
        """
        return config.tit_mole.getDelPhi()
    
    def __getitem__(self, numb):
        numb = self.correct_site_numb(numb)
        return self._pKas[numb]

    def __str__(self):
        output = '  Site   Name      pK'
        sites = config.tit_mole.getSites()

        for site in self._pKas:
            sitename = sites[site].getName()
            pk = self._pKas[site]
            if pk != '-':
                pk = round(pk, 2)
            else:
                pk = 'Not In Range'
            site = convertTermini(site)
            output += '\n{0:>6}    {1:3}    {2:>}'.format(site,
                                                          sitename, pk)
        return output


def CLI():
    # Assignment of global variables
    parameters, debug = checkParsedInput()
    Titration(parameters, sites=None, debug=debug)
    print('CLI exited successfully')


if __name__ == "__main__":
    CLI()
