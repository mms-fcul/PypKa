#! /usr/bin/python3

"""
A python API and CLI to perform pKa calculations on peptides,
proteins or lipid bilayers.
"""

from config import Config
from constants import *
from checksites import identify_tit_sites, check_sites_integrity, make_delphi_inputfile, get_chains_from_file
from cli import check_cli_args
from cleaning import inputPDBCheck, cleanPDB
from formats import convertTermini, gro2pdb, pdb2gro, read_pdb_line, new_pdb_line
from molecule import Molecule
from concurrency import startPoolProcesses, runDelPhiSims, runInteractionCalcs, runMCCalcs
import log
from ffconverter import AMBER_protomers, GROMOS_protomers, gromos2amber, mainchain_Hs

from delphi4py.delphi4py import DelPhi4py

import numpy as np
import os

__author__ = "Pedro Reis"
__version__ = "0.4"

__email__ = "pdreis@fc.ul.pt"
__status__ = "Development"


def getTitrableSites(pdb, ser_thr_titration=True, debug=False):
    """Gets the all titrable sites from the pdb

    Returns an input structure to the Titration class containing all
    titrable sites found in the pdb file.

    Args:
        pdb (str): The filename a PDB file

    Returns:
        A dict mapping all titrable sites found in the pdb file

    """

    parameters = {'structure': pdb,
                  'ncpus': 1,
                  'epsin': 15,
                  'ser_thr_titration': ser_thr_titration}
    Config.storeParams('', log.Log(), debug, parameters)

    chains = get_chains_from_file(pdb)
    sites = {chain: 'all' for chain in chains}

    molecules = {}
    for chain, site_list in sites.items():
        molecules[chain] = Molecule(chain, site_list)

    chains_res = identify_tit_sites(molecules, instanciate_sites=False)

    out_sites = {chain: [] for chain in chains_res.keys()}

    for chain, sites in chains_res.items():
        for resnumb, resname in sites.items():
            out_site = resnumb
            if resname == 'NTR':
                out_site += 'N'
            elif resname == 'CTR':
                out_site += 'C'

            out_sites[chain].append(out_site)

    return out_sites, chains_res


class Titration:
    """Main PypKa class

    Attributes:
        params (Config): method parameters
        molecules (dict): titrating molecules ordered by chain
    """

    def __init__(self, parameters, sites='all', debug=False, run='all'):
        """Runs the pKa prediction

        Args:
            parameters (dict): input parameters to overwrite the default configuration.
            sites (dict, optional): tuple titrable residue numbers indexed by chain.
                                    Defaults to 'all' which will titrate all titrable residues.
            debug (boolean, optional): debug mode switch. Defaults to False.
        """

        self.molecules = {}
        self.__parameters = Config.storeParams(self, log.Log(), debug, parameters)
        self.pKas = {}

        print('Start Preprocessing')
        self.preprocessing(sites)
        self.processDelPhiParams()

        if run == 'preprocess':
            return

        print('Start PB Calculations')
        self.DelPhiLaunch()

        if run == 'PB':
            return

        # Calculates sites interaction energies and write .dat file
        self.calcSiteInteractionsParallel()

        #  Monte Carlo sampling
        self.run_mc()

        if Config.pypka_params['f_structure_out']:
            self.writeOutputStructure()

        print('Results')
        print(self)

        print('API exited successfully')

    def preprocessing(self, sites):
        def create_tit_sites(chains_res, TMPpdb=False):
            for chain, molecule in self.molecules.items():
                molecule.deleteAllSites()
            check_sites_integrity(self.molecules, chains_res, useTMPpdb=TMPpdb)

        if Config.pypka_params['f_in_extension'] == 'gro':
            groname = Config.pypka_params['f_in']
            f_in = 'TMP.pdb'
            gro2pdb(groname, f_in, save_box=True)
            Config.pypka_params.redefine_f_in(f_in)


        # Creating instance of TitratingMolecule
        automatic_sites = False
        if sites == 'all':
            f_in = Config.pypka_params['f_in']
            chains = get_chains_from_file(f_in)
            sites = {chain: ['all'] for chain in chains}

        for chain, site_list in sites.items():
            if site_list == ['all']:
                site_list = 'all'
            else:
                site_list = [str(site) for site in site_list]
                sites[chain] = site_list
            self.molecules[chain] = Molecule(chain, site_list)
            if site_list == 'all':
                automatic_sites = True

        if Config.pypka_params['f_in_extension'] == 'pdb':
            f_in = Config.pypka_params['f_in']
            # Reading .st files
            # If the titrable residues are defined
            if not automatic_sites:
                clean_pdb = Config.pypka_params['clean_pdb']
                chains_length, chains_res = inputPDBCheck(f_in,
                                                          sites, clean_pdb)

                for chain, molecule in self.molecules.items():
                    molecule.loadSites(chains_res[chain])

            # If the titrable residues are not defined and
            # the input pdb file is incomplete
            elif Config.pypka_params['clean_pdb']:
                chains_res = identify_tit_sites(self.molecules)

            # If the titrable residues are not defined and
            # the input pdb is not missing any atoms
            else:
                chains_res = identify_tit_sites(self.molecules)
                create_tit_sites(chains_res)

            # Creates a .pdb input for DelPhi
            # where residues are in their standard state
            if Config.pypka_params['clean_pdb']:
                inputpqr = 'clean.pqr'
                outputpqr = 'cleaned_tau.pqr'
                #sites = {}
                #site_numb_n_ref = {}
                #for chain, molecule in self.molecules.items():
                #    sites[chain] = molecule.getSites()
                #    site_numb_n_ref[chain] = {}
                #    for site in sites[chain]:
                #        tau_ref_name = sites[chain][site].getRefTautomerName()
                #        site_numb_n_ref[chain][site] = tau_ref_name

                cleanPDB(self.molecules, chains_res, inputpqr, outputpqr)
                create_tit_sites(chains_res, TMPpdb=True)
                f_in = 'TMP.pdb'
        else:
            f_format = Config.pypka_params['f_in_extension']
            raise Exception("{0} file format is not currently supported".format(f_format))


        f_out = 'delphi_in_stmod.pdb'
        make_delphi_inputfile(f_in, f_out, self.molecules)
        if not Config.debug and os.path.isfile('TMP.pdb'):
            os.remove('TMP.pdb')

    def get_total_atoms(self):
        total_atoms = 0
        for molecule in self.molecules.values():
            total_atoms += molecule.getNAtoms()
        return total_atoms

    def get_all_sites(self, get_list=False):
        sites = {}
        if get_list:
            sites = []
        for chain, molecule in self.molecules.items():
            chain_sites = molecule.getSitesOrdered()
            if get_list:
                #if isinstance(chain_sites, dict):
                #    chain_sites = chain_sites.values()
                sites += chain_sites
            else:
                sites[chain] = chain_sites
        return sites

    def get_molecule_sites(self):
        sites = {}
        for chain, molecule in self.molecules.items():
            chain_sites = molecule.getSites()
            sites[chain] = chain_sites
        return sites

    def get_total_sites(self):
        total = 0
        for sites in self.get_all_sites().values():
            total += len(sites)
        return total

    def processDelPhiParams(self):
        # Storing DelPhi parameters and Creates DelPhi data structures
        logfile = 'LOG_readFiles'

        delphimol = DelPhi4py(Config.pypka_params['f_crg'],
                              Config.pypka_params['f_siz'],
                              'delphi_in_stmod.pdb',
                              Config.delphi_params['gsize'],
                              Config.delphi_params['scaleM'],
                              Config.delphi_params['precision'],
                              epsin=Config.delphi_params['epsin'],
                              epsout=Config.delphi_params['epssol'],
                              conc=Config.delphi_params['ionicstr'],
                              ibctyp=Config.delphi_params['bndcon'],
                              res2=Config.delphi_params['maxc'],
                              nlit=Config.delphi_params['nlit'],
                              nonit=Config.delphi_params['nonit'],
                              relfac=0.0,
                              relpar=0.0,
                              pbx=Config.delphi_params['pbx'],
                              pby=Config.delphi_params['pby'],
                              isurftype=Config.delphi_params['nanoshaper'],
                              debug=Config.debug, outputfile=logfile)

        if Config.pypka_params['f_structure_out']:
            with open('delphi_in_stmod.pdb') as f:
                self.delphi_input_content = f.readlines()
        if not Config.debug:
            os.remove('delphi_in_stmod.pdb')

        log.checkDelPhiErrors(logfile, 'readFiles')

        lookup_atoms = {}
        for molecule in Config.titration.molecules.values():
            for atom_name, atom_id, atom_position in molecule.iterAtoms():
                lookup_atoms[atom_position] = (atom_id, atom_name)

        # Stores DelPhi4py object and attributes
        Config.delphi_params.store_run_params(delphimol)
        Config.delphi_params.lookup_atoms = lookup_atoms

        if Config.debug:
            for chain, molecule in self.molecules.items():
                print('### CHAIN {chain} ###'.format(chain=chain))
                molecule.printAllSites()
                molecule.printAllTautomers()
            print(delphimol)

    def iterAllSitesTautomers(self):
        """Generator that iterates through all Tautomer instances.
        The iteration is sorted by site and within each site the first
        to be yielded is the reference tautomer and then the rest of
        the tautomers by order"""
        for molecule in self.molecules.values():
            for site in molecule.sites_order:
                yield site.ref_tautomer
                for tautomer in sorted(site.tautomers.values()):
                    yield tautomer


    def calcpKint(self, unpacked_results):
        """Calculation the pKint of all tautomers
        """
        i = -1
        if Config.debug:
            print('############ results ############')
        pkints = ''
        contributions = ''
        for tautomer in self.iterAllSitesTautomers():
            i += 1
            core_index = i % Config.pypka_params['ncpus']
            job_index = int(i / Config.pypka_params['ncpus'])
            result = unpacked_results[core_index][job_index]

            tautname    = result[0]
            tautresnumb = result[1]

            esolvationM = result[2]
            sitpotM     = result[3]
            esolvationS = result[4]
            sitpotS     = result[5]

            tautomer = Config.parallel_params.all_tautomers_order[i]

            tautomer.saveDelPhiResults(esolvationS, sitpotS, esolvationM,
                                       sitpotM)
            if Config.debug:
                print('### new tautomer ###')
                print((i, core_index, job_index, tautname,
                       tautomer.name, esolvationM, esolvationS))
                print((tautomer.name, tautomer.esolvationS,
                       len(tautomer.sitpotS), tautomer.esolvationM,
                       len(tautomer.sitpotM)))

            tautomer.calcBackEnergy()
            if not tautomer.isRefTautomer():
                tautomer.calcpKint()
                if Config.debug:
                    print(('pkint', tautomer.name, tautomer.dg,
                           id(tautomer)))
                    pkints += '{} {} {}\n'.format(tautomer.name,
                                                  tautomer.dg, tautomer.pKint)
                    contributions += '{}{} {} {} {} {}\n'.format(tautomer.site.res_number,
                                                                 tautomer.name,
                                                                 tautomer.dG_solvationM,
                                                                 tautomer.dG_solvationS,
                                                                 tautomer.dG_solvationM - tautomer.dG_solvationS,
                                                                 tautomer.dG_back)
        if Config.debug:
            with open('pkint', 'w') as f_new1, \
                 open('contributions', 'w') as f_new2:
                f_new1.write(pkints)
                f_new2.write(contributions)


    def calcSiteInteractionsParallel(self):
        """Calculates the pairwise interaction energies
        and writes them in a formatted .dat file
        Interactions are calculated using a pool of processes

        Args:
          ncpus (int): number of cpus to be used
        """
        def writeDatHeader(sites):
            """Writes pKint energies in .dat file header"""
            to_write = '{0}\n'.format(len(sites))
            for site in sites:
                to_write += '{0:3s}-{1:<7}{2:>2}  P  *\n'.format(site.res_name,
                                                                 site.res_number,
                                                                 len(site.tautomers) + 1)
                if site.type == 'c':
                    tau_prot_state = 0
                    ref_prot_state = 1
                elif site.type == 'a':
                    tau_prot_state = 1
                    ref_prot_state = 0

                for tautomer in site.iterOrderedTautomersWithoutRef():
                    to_write += '{0:1d} {1:13.6e}\n'.format(tau_prot_state,
                                                            tautomer.dg)
                to_write += '{0:1d}  0.000000e+00\n'.format(ref_prot_state)
            if Config.debug:
                with open('interactions.dat', 'w') as f_new:
                    f_new.write(to_write)

        Config.loadParams(self.__parameters)

        sites = self.get_all_sites(get_list=True)
        if Config.debug:
            writeDatHeader(sites)
        counter = 0
        site_interactions = []
        for site1 in sites[:-1]:
            counter += 1
            for site2 in sites[counter:]:
                site_interactions.append((site1, site2))

        npossible_states = [len(site.tautomers) + 1 for site in sites]
        Config.parallel_params.npossible_states = npossible_states

        interactions = []
        for nstate1 in range(sum(npossible_states)):
            interactions.append([])
            for _ in range(sum(npossible_states)):
                interactions[nstate1].append(-999999)

        interactions_look = []
        aux = -1
        site = -1
        for nstates in npossible_states:
            site += 1
            interactions_look.append([])
            for _ in range(nstates):
                aux += 1
                interactions_look[site].append(aux)

        Config.parallel_params.interactions_look = interactions_look
        Config.parallel_params.site_interactions = site_interactions
        Config.parallel_params.all_sites = sites
        Config.parallel_params.interactions = interactions
        ncpus = min(len(site_interactions), Config.pypka_params['ncpus'])
        results = []

        if ncpus > 0:
            results = startPoolProcesses(runInteractionCalcs,
                                         site_interactions, ncpus,
                                         assign='ordered', merged_results=True)

        interactions = Config.parallel_params.interactions
        to_write = ''
        temperature = float(Config.pypka_params['temp'])
        for interaction in results:
            site1i = interaction[0]

            site2i = interaction[1]
            if interactions[site1i][site2i] != -999999:
                interactions[site1i][site2i] += interaction[2]
                interactions[site2i][site1i] += interaction[2]
                interactions[site1i][site2i] /= 2
                interactions[site2i][site1i] /= 2

                if Config.debug:
                    col1 = interaction[3][6:18]
                    col2 = interaction[3][:6]
                    col3 = interactions[site1i][site2i] * (KBOLTZ * temperature)
                    to_write += '{0}{1} {2:13.6e}\n'.format(col1, col2, col3)
            else:
                interactions[site1i][site2i] = interaction[2]
                interactions[site2i][site1i] = interaction[2]


        if Config.debug:
            with open('interactions.dat', 'a') as f_new:
                f_new.write(to_write)

    def DelPhiLaunch(self):
        Config.loadParams(self.__parameters)

        if self.get_total_sites() < 1:
            raise Exception('At least one site has to be correctly defined.')

        all_tautomers = list(self.iterAllSitesTautomers())

        Config.parallel_params.all_tautomers_order = all_tautomers

        # Runs DelPhi simulations for all tautomers
        results = startPoolProcesses(runDelPhiSims,
                                     all_tautomers,
                                     Config.pypka_params['ncpus'])

        print('\rPB Runs Ended{:>80}'.format(""))

        # Calculates the pKint of all tautomers
        self.calcpKint(results)

    def run_mc(self):
        def resize_list_of_lists(listn, maxsize, filler=None):
            for i in listn:
                diff = maxsize - len(i)
                for ii in range(diff):
                    i.append(filler)

        def calcpKhalfs(pH, nsites, avgs, pmean, pKs, mcsteps, dpH):
            totalP = 0.0
            means = []
            for site in range(nsites):
                mean = avgs[site] / float(mcsteps)
                means.append(mean)
                totalP += mean

                p = pmean[site]

                if p > 0.5 and mean <= 0.5 or p < 0.5 and mean >= 0.5:
                    pKhalf = pH - dpH * (mean - 0.5) / (mean - p)
                    for i in range(MAXNPKHALFS):
                        if pKs[site][i] == PKAPLACEHOLDER:
                            pKs[site][i] = pKhalf
                            break

                elif p > 1.5 and mean <= 1.5 or p < 1.5 and mean >= 1.5:
                    pKhalf = pH - dpH * (mean - 1.5) / (mean - p)
                    for i in range(MAXNPKHALFS):
                        if pKs[site][i] == PKAPLACEHOLDER:
                            pKs[site][i] = pKhalf
                            break

                pmean[site] = mean
            totalP /= nsites
            return totalP, pKs, pmean, means

        Config.loadParams(self.__parameters)

        print('\nStart MC', end='\r')

        sites = self.get_all_sites(get_list=True)
        nsites = len(sites)
        #possible_states     = [[] for _ in sites]
        possible_states_g   = [[] for _ in sites]
        possible_states_occ = [[] for _ in sites]

        temperature = float(Config.pypka_params['temp'])
        isite = -1
        for site in sites:
            isite += 1
            itaut = 0
            for tautomer in site.iterOrderedTautomersWithoutRef():
                dg = tautomer.dg / (KBOLTZ * temperature)
                possible_states_g[isite].append(dg)
                #possible_states[isite].append(itaut)
                if site.type == 'c':
                    prot_state = 0
                elif site.type == 'a':
                    prot_state = 1

                possible_states_occ[isite].append(prot_state)
                itaut += 1

            if site.type == 'c':
                prot_state = 1
            elif site.type == 'a':
                prot_state = 0
            possible_states_occ[isite].append(prot_state)
            possible_states_g[isite].append(0.0)

        maxstates = max(Config.parallel_params.npossible_states)
        interactions_look = Config.parallel_params.interactions_look
        #resize_list_of_lists(possible_states, maxstates)
        resize_list_of_lists(possible_states_g, maxstates)
        resize_list_of_lists(possible_states_occ, maxstates, filler=-500)
        resize_list_of_lists(interactions_look, maxstates, filler=-500)

        params = Config.mc_params
        pHmin, pHmax = params['pHmin'], params['pHmax']
        dpH = params['pHstep']
        pHsteps = int(round(1 + (pHmax - pHmin) / dpH, 0))


        Config.parallel_params.possible_states_g   = possible_states_g
        Config.parallel_params.possible_states_occ = possible_states_occ

        ncpus = min(Config.pypka_params['ncpus'], nsites)
        results = startPoolProcesses(runMCCalcs,
                                     list(range(pHsteps)), ncpus,
                                     assign='ordered', merged_results=True)

        print('\rMC Runs Ended{:>80}\n'.format(""))

        counts_all = []
        avgs_all = []
        cur_states = []

        for i in results:
            avgs_all.append(i[0])
            counts_all.append(i[1])
            cur_states.append(i[3])

        pKs = np.array([[PKAPLACEHOLDER for ii in range(MAXNPKHALFS)] for i in range(nsites)])
        mcsteps = params['mcsteps']
        pmeans = avgs_all[0] / float(mcsteps)

        tit_curve = {}
        tit_curve[pHmin] = {}
        tit_curve[pHmin]['total'] = sum(avgs_all[0]) / mcsteps / nsites
        for i, mean in enumerate(avgs_all[0]):
            site = sites[i]
            sitenumber = site.getResNumber()
            tit_curve[pHmin][sitenumber] = mean / mcsteps

        sites = Config.parallel_params.all_sites

        for pHstep in range(1, pHsteps):
            pH = pHmin + pHstep * dpH
            totalP, pKs, pmeans, means = calcpKhalfs(pH, nsites, avgs_all[pHstep],
                                                     pmeans, pKs, mcsteps, dpH)
            tit_curve[pH] = {}
            tit_curve[pH]['total'] = totalP
            for i, mean in enumerate(means):
                site = sites[i]
                sitenumber = site.getResNumber()
                tit_curve[pH][sitenumber] = mean

        pKas = pKs

        text_pks = ''
        text_prots = '#pH       total'
        c = -1
        for i in pKas:
            c += 1
            site = sites[c]
            sitename = site.getName()
            if sitename in ('NTR', 'CTR'):
                resnumb = site.getResNumber() - TERMINAL_OFFSET
                text_prots += '     {0:3}'.format(sitename)
            else:
                resnumb = site.getResNumber()
                text_prots += '{0:5d}{1:3s}'.format(resnumb, sitename)
            text_pks += '{0} {1} {2}\n'.format(resnumb, sitename, i[0])

        final_states = {}
        state_distribution = {}
        most_prob_states = {}
        self.tit_curve = {}
        for pHstep in range(pHsteps):
            pH = pHmin + pHstep * dpH
            text_prots += '\n{pH:5.2f}'.format(pH=pH)

            self.tit_curve[pH] = tit_curve[pH]['total']

            final_states[pH] = {}
            state_distribution[pH] = {}
            most_prob_states[pH] = {}
            for c, site in enumerate(sites):
                sitename = site.getName()
                sitenumb = site.getResNumber()
                mean = tit_curve[pH][sitenumb]
                final_states[pH][sitenumb] = cur_states[pHstep][c]
                state_distribution[pH][sitenumb] = list(counts_all[pHstep][c] / mcsteps)

                ntauts = site.getNTautomers()
                ref_i = ntauts
                prot_state = site.getRefProtState()
                if (mean > 0.5 and prot_state == 1) or \
                   (mean <= 0.5 and prot_state == -1):
                    state_i = ref_i
                else:
                    max_prob = max(state_distribution[pH][sitenumb][:ref_i])
                    state_i = state_distribution[pH][sitenumb].index(max_prob)

                most_prob_state = state_i + 1
                most_prob_states[pH][sitenumb] = most_prob_state

                if mean != PKAPLACEHOLDER:
                    text_prots += '\t{mean:7.4f}'.format(mean=mean)
                else:
                    text_prots += '\t-'

                site.most_prob_states[pH] = most_prob_state
                site.final_states[pH] = final_states[pH][sitenumb]
                site.tit_curve[pH] = tit_curve[pH][sitenumb]
                site.states_prob[pH] = state_distribution[pH][sitenumb]

        if Config.pypka_params['f_out']:
            with open(Config.pypka_params['f_out'], 'w') as f:
                f.write(text_pks)
        if Config.pypka_params['f_prot_out']:
            with open(Config.pypka_params['f_prot_out'], 'w') as f:
                f.write(text_prots)

        #self.tit_curve = tit_curve
        #self.pH_values = sorted(tit_curve.keys())
        #self.final_states = final_states
        #self.state_prob = state_distribution
        #self.most_prob_states = most_prob_states

        c = -1
        for i in pKas:
            c += 1
            #site = sites[c].res_number
            #chain = sites[c].molecule.chain
            #if chain not in self.pKas.keys():
            #    self.pKas[chain] = {}
            pK = i[0]
            sites[c].setpK(pK)
            #self.pKas[chain][site] = pK

    def writeOutputStructure(self):
        def getProtomerResname(pdb_content, site, pH, ff_protomers):
            resnumb = site.getResNumber()
            resname = site.getName()
            new_state, new_state_prob = site.getMostProbTaut(pH)
            new_state_i = new_state - 1
            for ff_resname, protomers in ff_protomers[resname].items():
                if new_state_i in protomers.keys():
                    new_resname = ff_resname
                    remove_hs = protomers[new_state_i]

                    state_prob, taut_prob = site.getTautProb(new_state,
                                                              pH)

                    if resnumb > TERMINAL_OFFSET:
                        resnumb -= TERMINAL_OFFSET
                    if state_prob < 0.75:
                        warn = '{0}{1} ' \
                               'protonation state probability: {2}, ' \
                               'tautomer probability: {3}'.format(resname, resnumb, state_prob, taut_prob)
                        Config.log.report_warning(warn)

                        print(warn)
                    rounded_sprob = round(state_prob, 2)
                    rounded_tprob = round(taut_prob, 2)
                    remark_line = '{0: <5}{1: <10}{2: ^7}'\
                                  '{3: >1.2f}{4: ^13}{5: >1.2f}'.format(resname, resnumb,
                                                                        "", rounded_sprob,
                                                                        "", rounded_tprob)

                    pdb_content += 'REMARK     {text}\n'.format(text=remark_line)

            #print(resnumb, new_state, new_resname, remove_hs, state_prob, taut_prob)
            return pdb_content, new_state_i, new_resname, remove_hs

        outputname = Config.pypka_params['f_structure_out']
        pH = float(Config.pypka_params['f_structure_out_pH'])
        ff_out = Config.pypka_params['ff_structure_out']

        ff_protomer = {'amber': AMBER_protomers,
                       'gromos_cph': GROMOS_protomers}[ff_out]

        pdb_content = 'REMARK     Protonation states assigned according to PypKa\n'\
                      'REMARK     Residue    Prot State Prob    Tautomer Prob\n'


        sites = self.get_all_sites(get_list=True)
        new_states = {}
        for site in sites:
            resname = site.getName()
            resnumb = site.getResNumber()
            molecule = site.molecule
            chain = molecule.chain

            (pdb_content, new_state,
             new_resname, remove_hs) = getProtomerResname(pdb_content,
                                                          site, pH,
                                                          ff_protomer)

            if resname in ('NTR', 'CTR'):
                new_resname = site.termini_resname

            if chain not in new_states:
                new_states[chain] =  {}

            new_states[resnumb] = (resname, new_state,
                                   new_resname, remove_hs)

        new_pdb = pdb_content
        counter = 0

        tit_atoms = {}
        other_atoms = {}
        for molecule in self.molecules.values():
            for atom_numb in molecule.atoms_tit_res:
                if molecule.atoms_tit_res[atom_numb]:
                    tit_atoms[atom_numb] = molecule
                else:
                    other_atoms[atom_numb] = molecule

        for line in self.delphi_input_content:
            if line.startswith('ATOM '):
                (aname, anumb, resname, chain,
                 resnumb, x, y, z) = read_pdb_line(line)

                if anumb in tit_atoms.keys():
                    molecule = tit_atoms[anumb]

                    (oldresname, new_state,
                     resname, removeHs) = new_states[resnumb]

                    if aname in removeHs:
                        continue

                    if ff_out == 'amber' and \
                       oldresname in gromos2amber and \
                       new_state in gromos2amber[oldresname] and \
                       aname in gromos2amber[oldresname][new_state]:
                        aname = gromos2amber[oldresname][new_state][aname]
                else:
                    molecule = other_atoms[anumb]

                if resnumb > TERMINAL_OFFSET:
                    termini_site = molecule.sites[resnumb]
                    resnumb -= TERMINAL_OFFSET
                    if resnumb in molecule.sites.keys():
                        ter_resname, ter_new_state, \
                        resname, ter_removeHs = new_states[resnumb]
                    else:
                        resname = termini_site.termini_resname

                    #print(new_pdb_line(anumb, aname, resname, resnumb, x, y, z).strip())
                if resnumb in molecule.getCYS_bridges():
                    resname = 'CYX'

                counter += 1
                new_pdb += new_pdb_line(counter, aname, resname, resnumb, x, y, z)
                if resnumb in mainchain_Hs:
                    while len(mainchain_Hs[resnumb]) > 0:
                        counter += 1
                        (aname, anumb, oldresname, chain, \
                         x, y, z) = mainchain_Hs[resnumb].pop()
                        new_pdb += new_pdb_line(counter, aname, resname, resnumb, x, y, z)
                    del mainchain_Hs[resnumb]
            else:
                new_pdb += line

        with open(outputname, 'w') as f_new:
            f_new.write(new_pdb)

    def getTitrationCurve(self):
        return self.tit_curve

    def __iter__(self):
        self.itersites = self.get_all_sites(get_list=True)
        self.iternumb = -1
        self.itermax = len(self.itersites)
        return self

    def __next__(self):
        self.iternumb += 1
        if self.iternumb < self.itermax:
            site = self.itersites[self.iternumb]
            return site
        else:
            raise StopIteration

    def getParameters(self):
        """Get the parameters used in the calculations
        """
        return '{}\n{}\n{}'.format(Config.pypka_params.__str__(),
                                           Config.delphi_params.__str__(),
                                           Config.mc_params.__str__())

    def __getitem__(self, chain):
        return self.molecules[chain]

    def __str__(self):
        output = 'Chain  Site   Name      pK'
        sites = self.get_all_sites()

        for chain in sites.keys():
            for site in sites[chain]:
                pk = site.pK
                if pk:
                    pk = round(pk, 2)
                else:
                    pk = 'Not In Range'
                #site = convertTermini(site)
                output += '\n{0:>4} {1:>6}    {2:3}    {3:>}'.format(chain,
                                                                     site.res_number,
                                                                     site.res_name, pk)
        return output

def CLI():
    # Read command line arguments
    sites, parameters, debug = check_cli_args()

    Titration(parameters, sites=sites, debug=debug)
    print('CLI exited successfully')

if __name__ == "__main__":
    CLI()