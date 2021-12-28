import os
from delphi4py import DelPhi4py
from pdbmender.formats import gro2pdb, get_grobox_size, get_chains_from_file
from pdbmender.utils import identify_tit_sites

from pypka.log import checkDelPhiErrors
from pypka.clean.checksites import (
    check_sites_integrity,
    make_delphi_inputfile,
    fix_fixed_sites,
)
from pypka.clean.cleaning import cleanPDB, inputPDBCheck
from pypka.clean.pdb_out import write_output_structure
from pypka.concurrency import (
    runDelPhiSims,
    runInteractionCalcs,
    startPoolProcesses,
)
from pypka.config import Config
from pypka.constants import KBOLTZ, TITRABLETAUTOMERS, TERMINAL_OFFSET
from pypka.molecule import Molecule
from pypka.mc.run_mc import MonteCarlo

import logging

logger = logging.getLogger()
logger.setLevel(logging.DEBUG)
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(levelname)s: %(message)s")
ch.setFormatter(formatter)
logger.addHandler(ch)


class Titration:
    """
    Main PypKa class

    Attributes:
        params (Config): method parameters
        molecules (dict): titrating molecules ordered by chain
    """

    def __init__(
        self, parameters, sites="all", fixed_sites=None, debug=False, run="all"
    ):
        """
        Runs the pKa prediction

        Args:
            parameters (dict): input parameters to overwrite the default configuration.
            sites (dict, optional): tuple titrable residue numbers indexed by chain.
                                    Defaults to 'all' which will titrate all titrable residues.
            debug (boolean, optional): debug mode switch. Defaults to False.
        """
        self.molecules = {}
        self.__parameters = Config.storeParams(self, debug, parameters, sites)
        self.pKas = {}
        self.isoelectric_point = None
        self.isoelectric_point_limit = None

        print("Start Preprocessing")
        self.preprocessing(sites, fixed_sites)
        self.processDelPhiParams()

        if run == "preprocess":
            return

        print("Start PB Calculations")
        self.DelPhiLaunch()

        if run == "PB":
            return

        print("Calculating site-site interactions")
        self.calcSiteInteractionsParallel()

        print("Start MC", end="\r")
        self.run_mc()

        if Config.pypka_params["f_structure_out"]:
            sites = self.get_all_sites(get_list=True)
            write_output_structure(sites, self.molecules, self.delphi_input_content)

        print("Results")
        print(self)

        print("API exited successfully")

    def preprocessing(self, sites, fixed_sites):
        def create_tit_sites(chains_res, TMPpdb=False):
            for _, molecule in self.molecules.items():
                molecule.deleteAllSites()
            check_sites_integrity(self.molecules, chains_res, useTMPpdb=TMPpdb)

        if fixed_sites:
            for chain, fixed_sites_info in fixed_sites.items():
                if chain not in sites:
                    sites[chain] = []
                for sitenumb in fixed_sites_info.keys():
                    sites[chain].append(sitenumb)

        if Config.pypka_params["f_in_extension"] == "gro":
            groname = Config.pypka_params["f_in"]
            f_in = "TMP.pdb"
            gro2pdb(groname, f_in)
            box = get_grobox_size(groname)
            Config.pypka_params.setBox(box)
            Config.pypka_params.redefine_f_in(f_in)

        if not Config.pypka_params["f_in_extension"] == "pdb":
            f_format = Config.pypka_params["f_in_extension"]
            raise Exception(
                "{0} file format is not currently supported".format(f_format)
            )
        automatic_sites = self.create_molecule(sites)

        f_in = Config.pypka_params["f_in"]
        if not automatic_sites:
            # Reading .st files
            # If the titrable residues are defined
            clean_pdb = Config.pypka_params["clean_pdb"]
            _, chains_res = inputPDBCheck(f_in, sites, clean_pdb)
            for chain, molecule in self.molecules.items():
                molecule.loadSites(chains_res[chain])

        else:
            # If the titrable residues are not defined
            chains_res = identify_tit_sites(
                Config.pypka_params["f_in"],
                self.molecules.keys(),
                nomenclature=Config.pypka_params["ffinput"],
                add_ser_thr=Config.pypka_params["ser_thr_titration"],
            )
            self.assign_residues_to_molecule(chains_res)

            if not Config.pypka_params["clean_pdb"]:
                # if the input pdb is not missing any atoms
                create_tit_sites(chains_res)

        if Config.pypka_params["clean_pdb"]:
            # Creates a .pdb input for DelPhi
            # where residues are in their standard state
            inputpqr = "clean.pqr"
            outputpqr = "cleaned_tau.pqr"
            cleanPDB(self.molecules, chains_res, inputpqr, outputpqr, automatic_sites)
            create_tit_sites(chains_res, TMPpdb=True)
            f_in = "TMP.pdb"

        f_out = "delphi_in_stmod.pdb"
        self.sequence = make_delphi_inputfile(f_in, f_out, self.molecules)

        if fixed_sites:
            fix_fixed_sites(self.molecules, fixed_sites, f_out)

        if not Config.debug and os.path.isfile("TMP.pdb"):
            os.remove("TMP.pdb")

    def create_molecule(self, sites):
        # Creating instance of TitratingMolecule
        automatic_sites = False
        if sites == "all":
            f_in = Config.pypka_params["f_in"]
            chains = get_chains_from_file(f_in)
            sites = {chain: ["all"] for chain in chains}

        for chain, site_list in sites.items():
            if site_list == ["all"]:
                site_list = "all"
            else:
                site_list = [str(site) for site in site_list]
                sites[chain] = site_list
            self.molecules[chain] = Molecule(chain, site_list)
            if site_list == "all":
                automatic_sites = True
        return automatic_sites

    def assign_residues_to_molecule(self, chains_res):
        if not chains_res:
            f_in = Config.pypka_params["f_in"]
            raise Exception("Not one titrable residue was found in {}".format(f_in))

        for chain, residues in chains_res.items():
            molecule = self.molecules[chain]
            for resnumb, resname in residues.items():
                if resname in ("NTR", "CTR"):
                    resnumb = int(resnumb)
                    if resname == "NTR":
                        molecule.NTR = resnumb
                    elif resname == "CTR":
                        molecule.CTR = resnumb
                    resnumb += TERMINAL_OFFSET

                sID = molecule.addSite(resnumb)
                ntautomers = TITRABLETAUTOMERS[resname]
                molecule.addTautomers(sID, ntautomers, resname)

        for molecule in self.molecules.values():
            # Adding the reference tautomer to each site
            molecule.addReferenceTautomers()
            # Assigning a charge set to each tautomer
            molecule.addTautomersChargeSets()

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
                # if isinstance(chain_sites, dict):
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
        logfile = "LOG_readFiles"

        delphimol = DelPhi4py(
            Config.pypka_params["f_crg"],
            Config.pypka_params["f_siz"],
            "delphi_in_stmod.pdb",
            igrid=Config.delphi_params["gsize"],
            scale=Config.delphi_params["scaleM"],
            precision=Config.delphi_params["precision"],
            epsin=Config.delphi_params["epsin"],
            epsout=Config.delphi_params["epssol"],
            conc=Config.delphi_params["ionicstr"],
            ibctyp=Config.delphi_params["bndcon"],
            res2=Config.delphi_params["maxc"],
            nlit=Config.delphi_params["nlit"],
            nonit=Config.delphi_params["nonit"],
            relfac=0.0,
            relpar=0.0,
            pbx=Config.delphi_params["pbx"],
            pby=Config.delphi_params["pby"],
            isurftype=Config.delphi_params["nanoshaper"],
            debug=Config.debug,
            outputfile=logfile,
        )

        if Config.pypka_params["f_structure_out"]:
            with open("delphi_in_stmod.pdb") as f:
                self.delphi_input_content = f.readlines()
        if Config.pypka_params["save_pdb"]:
            with open("delphi_in_stmod.pdb") as f, open(
                Config.pypka_params["save_pdb"], "w"
            ) as f_new:
                content = f.read()
                f_new.write(content)

        if not Config.debug:
            os.remove("delphi_in_stmod.pdb")

        checkDelPhiErrors(logfile, "readFiles")

        lookup_atoms = {}
        for molecule in Config.titration.molecules.values():
            for atom_name, atom_id, atom_position in molecule.iterAtoms():
                lookup_atoms[atom_position] = (atom_id, atom_name)

        # Stores DelPhi4py object and attributes
        Config.delphi_params.store_run_params(delphimol)
        Config.delphi_params.lookup_atoms = lookup_atoms

        if Config.debug:
            for chain, molecule in self.molecules.items():
                print("### CHAIN {chain} ###".format(chain=chain))
                molecule.printAllSites()
                molecule.printAllTautomers()
            print(delphimol)

    def iterAllSitesTautomers(self):
        """
        Generator that iterates through all Tautomer instances.

        The iteration is sorted by site and within each site the first
        to be yielded is the reference tautomer and then the rest of
        the tautomers by order
        """
        for molecule in self.molecules.values():
            for site in molecule.sites_order:
                yield site.ref_tautomer
                for tautomer in sorted(site.tautomers.values()):
                    yield tautomer

    def calcpKint(self, unpacked_results):
        """Calculation the pKint of all tautomers."""
        i = -1
        if Config.debug:
            print("############ results ############")
        pkints = ""
        contributions = ""
        for tautomer in self.iterAllSitesTautomers():
            i += 1
            core_index = i % Config.pypka_params["ncpus"]
            job_index = int(i / Config.pypka_params["ncpus"])
            result = unpacked_results[core_index][job_index]

            tautname = result[0]
            # tautresnumb = result[1]

            esolvationM = result[2]
            sitpotM = result[3]
            esolvationS = result[4]
            sitpotS = result[5]

            tautomer = Config.parallel_params.all_tautomers_order[i]

            tautomer.saveDelPhiResults(esolvationS, sitpotS, esolvationM, sitpotM)
            if Config.debug:
                print("### new tautomer ###")
                print(
                    (
                        i,
                        core_index,
                        job_index,
                        tautname,
                        tautomer.name,
                        esolvationM,
                        esolvationS,
                    )
                )
                print(
                    (
                        tautomer.name,
                        tautomer.esolvationS,
                        len(tautomer.sitpotS),
                        tautomer.esolvationM,
                        len(tautomer.sitpotM),
                    )
                )

            tautomer.calcBackEnergy()
            if not tautomer.isRefTautomer():
                tautomer.calcpKint()
                if Config.debug:
                    print(("pkint", tautomer.name, tautomer.dg, id(tautomer)))
                    pkints += "{} {} {}\n".format(
                        tautomer.name, tautomer.dg, tautomer.pKint
                    )
                    contributions += "{}{} {} {} {} {}\n".format(
                        tautomer.site.res_number,
                        tautomer.name,
                        tautomer.dG_solvationM,
                        tautomer.dG_solvationS,
                        tautomer.dG_solvationM - tautomer.dG_solvationS,
                        tautomer.dG_back,
                    )
        if Config.debug:
            with open("pkint", "w") as f_new1, open("contributions", "w") as f_new2:
                f_new1.write(pkints)
                f_new2.write(contributions)

    def calcSiteInteractionsParallel(self):
        """
        Calculates the pairwise interaction energies

        Interactions are calculated using a pool of processes
        and written in a formatted .dat file

        Args:
          ncpus (int): number of cpus to be used
        """

        def writeDatHeader(sites):
            """Writes pKint energies in .dat file header."""
            to_write = "{0}\n".format(len(sites))
            for site in sites:
                to_write += "{0:3s}-{1:<7}{2:>2}  P  *\n".format(
                    site.res_name, site.res_number, len(site.tautomers) + 1
                )
                if site.type == "c":
                    tau_prot_state = 0
                    ref_prot_state = 1
                elif site.type == "a":
                    tau_prot_state = 1
                    ref_prot_state = 0

                for tautomer in site.iterOrderedTautomersWithoutRef():
                    to_write += "{0:1d} {1:13.6e}\n".format(tau_prot_state, tautomer.dg)
                to_write += "{0:1d}  0.000000e+00\n".format(ref_prot_state)
            if Config.debug:
                with open("interactions.dat", "w") as f_new:
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
        ncpus = min(len(site_interactions), Config.pypka_params["ncpus"])
        results = []

        if ncpus > 0:
            results = startPoolProcesses(
                runInteractionCalcs,
                site_interactions,
                ncpus,
                assign="ordered",
                merged_results=True,
            )

        interactions = Config.parallel_params.interactions
        to_write = ""
        temperature = float(Config.pypka_params["temp"])
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
                    to_write += "{0}{1} {2:13.6e}\n".format(col1, col2, col3)
            else:
                interactions[site1i][site2i] = interaction[2]
                interactions[site2i][site1i] = interaction[2]

        if Config.debug:
            with open("interactions.dat", "a") as f_new:
                f_new.write(to_write)

    def DelPhiLaunch(self):
        Config.loadParams(self.__parameters)

        if self.get_total_sites() < 1:
            raise Exception("At least one site has to be correctly defined.")

        all_tautomers = list(self.iterAllSitesTautomers())

        Config.parallel_params.all_tautomers_order = all_tautomers

        # Runs DelPhi simulations for all tautomers
        results = startPoolProcesses(
            runDelPhiSims, all_tautomers, Config.pypka_params["ncpus"]
        )

        print("\rPB Runs Ended{:>80}".format(""))

        # Calculates the pKint of all tautomers
        print("Calculating intrinsic pK values")
        self.calcpKint(results)

    def run_mc(self):
        Config.loadParams(self.__parameters)

        sites = self.get_all_sites(get_list=True)

        mc = MonteCarlo(sites)
        text_pks, text_prots, self.tit_curve = (
            mc.text_pks,
            mc.text_prots,
            mc.total_tit_curve,
        )
        print("\rMC Runs Ended{:>80}\n".format(""))

        if Config.pypka_params["isoelectric_point"]:
            self.getIsoelectricPoint()
            if self.isoelectric_point:
                text_pks += "\nPredicted Isoelectric Point: {0}\n".format(
                    self.isoelectric_point
                )
            elif self.isoelectric_point_limit:
                pH, limit, charge = self.isoelectric_point_limit
                text_pks += "\nPredicted Isoelectric Point: {0} Limit: {1} \nCharge at pH {0}: {2}\n".format(
                    pH, limit, charge
                )

        if Config.pypka_params["f_out"]:
            with open(Config.pypka_params["f_out"], "w") as f:
                f.write(text_pks)
        if Config.pypka_params["f_prot_out"]:
            with open(Config.pypka_params["f_prot_out"], "w") as f:
                f.write(text_prots)

    def getTotalProtonationCurve(self):
        return self.tit_curve

    def getTitrationCurve(self):
        if Config.pypka_params["sites"] != "all":
            warn_message = "The titration curve weighted by charges can only be calculated when titrating all sites."
            raise Exception(warn_message)

        sites = self.get_all_sites(get_list=True)
        nsites = len(sites)
        total_anionic_titres = 0
        for site in sites:
            if site.type == "a":
                total_anionic_titres += 1

        total_arg_charge = 0
        for _, aas in self.sequence.items():
            for _, resname in aas.items():
                if resname == "ARG":
                    total_arg_charge += 1

        titration_curve = {}
        for pH, prot in self.tit_curve.items():
            prot = prot * nsites
            charge = prot - total_anionic_titres + total_arg_charge
            titration_curve[pH] = charge
        return titration_curve

    def getIsoelectricPoint(self):
        def calc_isoelectric_point():
            # anionic_aas = ("ASP", "GLU", "CYS", "TYR", "CTR")
            # cationic_aas = ("LYS", "ARG", "HIS", "NTR")

            titration_curve = [
                (pH, prot) for pH, prot in self.getTitrationCurve().items()
            ]
            i = 0
            pH, prot = titration_curve[i]
            while prot > 0.0 and i < len(titration_curve) - 1:
                i += 1
                pH, prot = titration_curve[i]

            if i == 0 or i == len(titration_curve) - 1:
                charge = round(prot, 2)
                warn_message = (
                    "Not enough points to interpolate the Isoelectric Point\n"
                    "Lowest pH value: {0}\t"
                    "Charge: {1}".format(pH, charge)
                )
                print(warn_message)
                limit = "-" if i == 0 else "+"
                return pH, limit, charge
            else:
                last_pH, last_prot = titration_curve[i - 1]

                isoelectric_point = pH - (prot * (last_pH - pH)) / (last_prot - prot)
                return isoelectric_point

        if not self.isoelectric_point:
            result = calc_isoelectric_point()
            if isinstance(result, tuple):
                self.isoelectric_point_limit = result
            else:
                self.isoelectric_point = result
        else:
            result = self.isoelectric_point

        return result

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

    @staticmethod
    def getParameters():
        """Get the parameters used in the calculations."""
        return "{}\n{}\n{}".format(
            Config.pypka_params.__str__(),
            Config.delphi_params.__str__(),
            Config.mc_params.__str__(),
        )

    @staticmethod
    def getParametersDict():
        return (
            Config.pypka_params.get_clean_params(),
            Config.delphi_params.get_clean_params(),
            Config.mc_params.get_clean_params(),
        )

    def getSiteInteractions(self):
        Config.loadParams(self.__parameters)

        return (
            Config.parallel_params.all_sites,
            Config.parallel_params.npossible_states,
            Config.parallel_params.interactions_look,
            Config.parallel_params.interactions,
        )

    def __getitem__(self, chain):
        return self.molecules[chain]

    def __str__(self, str_format="print"):
        if str_format == "print":
            output = "Chain  Site   Name      pK"
        elif str_format == "file":
            output = ""
        sites = self.get_all_sites()

        for chain in sites.keys():
            for site in sites[chain]:
                pk = site.pK
                if pk and str_format == "print":
                    pk = "{:.2f}".format(round(pk, 2))
                elif not pk:
                    pk = "Not In Range"
                output += "\n{0:>4} {1:>6}    {2:3}    {3:>5}".format(
                    chain, site.getResNumber(), site.res_name, pk
                )

        if self.isoelectric_point:
            output += "\n\nPredicted Isoelectric Point: {0}\n".format(
                self.isoelectric_point
            )

        return output


def getTitrableSites(pdb, ser_thr_titration=True, debug=False):
    """
    Gets the all titrable sites from the pdb

    Returns an input structure to the Titration class containing all
    titrable sites found in the pdb file.

    Args:
        pdb (str): The filename a PDB file

    Returns:
        A dict mapping all titrable sites found in the pdb file
    """
    parameters = {
        "structure": pdb,
        "ncpus": 1,
        "epsin": 15,
        "ser_thr_titration": ser_thr_titration,
    }
    Config.storeParams("", debug, parameters, "all")

    chains = get_chains_from_file(pdb)
    sites = {chain: "all" for chain in chains}

    molecules = {}
    for chain, site_list in sites.items():
        molecules[chain] = Molecule(chain, site_list)

    chains_res = identify_tit_sites(pdb, chains)

    out_sites = {chain: [] for chain in chains_res.keys()}

    for chain, sites in chains_res.items():
        for resnumb, resname in sites.items():
            out_site = resnumb
            if resname == "NTR":
                out_site += "N"
            elif resname == "CTR":
                out_site += "C"

            out_sites[chain].append(out_site)

    return out_sites, chains_res
