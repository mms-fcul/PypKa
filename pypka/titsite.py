from pypka.config import Config
from pypka.constants import KBOLTZ, TERMINAL_OFFSET
from pypka.tautomer import Tautomer


class Titsite:
    """Titrable Site with more than one Tautomer"""

    def __init__(self, res_number, molecule):
        """
        Args:
          res_number (int): number of the site in the .pdb
                            and .sites file
          molecule (TitratingMolecule): molecule to
                                        which the Site belongs

        # Site Details
        self.res_name (str): residue name
        self.type (str): site type ('c' or 'a')
                          'c' if site is cationic
                          'a' if site is anionic
        self.center (list) : Site focus center
            if dimension is 0, then center is defined in cent file
            elif dimension is 2, then center (z) is defined in cent file
                                 and center (x, y) are defined in boxsize variable
        self.center_original (list): Site geometric center on
                                      cent file when dimension is 2

        # Site Tautomers
                                key is tautomer name
                                value is tautomer instance
        self.ref_tautomer (Tautomer): reference tautomer instance

        # Site Atoms
        self.atoms (dict): existing atoms in site
                            key is atom id name
                            value is atom id number
        """
        self.molecule = molecule
        self.res_number = res_number
        self.res_name = ""
        self.termini_resname = ""
        self.atoms = {}
        self.tautomers = {}
        self.ref_tautomer = ""
        self.center = []
        self.type = ""
        self.pK = None

        self.most_prob_states = {}
        self.final_states = {}
        self.tit_curve = {}
        self.states_prob = {}

    # Set Methods
    def setTautomers(self, ntauts, resname):
        """Adds Tautomers from res_tauts to self.tautomers

        Args:
            res_tauts (list): tautomers from sites file
        """
        for i in range(ntauts):
            self.res_name = resname
            correct_number = str(i)
            tautomer = resname[:2] + correct_number

            tID = Tautomer(tautomer, self, self.molecule)
            self.tautomers[tautomer] = tID

    def setpK(self, pK):
        self.pK = pK

    def setTerminiResname(self, resname):
        self.termini_resname = resname

    def addReferenceTautomer(self):
        """Gets last tautomer from .sites file adds one and saves it as the
        reference tautomer"""
        last_tautomer = max(sorted(self.tautomers.keys()))
        correct_number = str(int(last_tautomer[-1]) + 1)
        ref_tautomer = last_tautomer[:2] + correct_number
        tID = Tautomer(ref_tautomer, self, self.molecule)
        self.ref_tautomer = tID

    def addChargeSets(self):
        """Stores the charge set of each existing tautomer
        for the present Site"""
        for tautomer in list(self.tautomers.values()):
            tautomer.loadChargeSet(self.res_name, self.ref_tautomer)

    def addAtom(self, aname, anumb):
        """
        Args:
            aname (str): atom name
            anumb (int): atom id number
        """
        self.atoms[aname] = anumb

    def addCenter(self, center, boxsize=False, box_z=False):
        x = center[0]
        y = center[1]
        z = center[2]
        if boxsize:
            if Config.delphi_params["pbc_dim"] != 2:
                raise Exception(
                    "ERROR: The original center is only "
                    "needed for 2-dimensional calculation"
                )
            self.center_original = [x, y, z]
            x = boxsize / 2
            y = boxsize / 2
            z += box_z
        self.center = [x, y, z]

    def addCenterH(self, center):
        self.centerH = center

    def setType(self, stype):
        self.type = stype

    # Get Methods
    def getName(self):
        return self.res_name

    def getTautomers(self):
        """Returns list of all tautomers instances
        except the tautomers of reference"""
        return list(self.tautomers.values())

    def getNTautomers(self):
        return len(self.getTautomers())

    def getAtomNumbersList(self):
        return list(self.atoms.values())

    def getAtomNamesList(self):
        return list(self.atoms.keys())

    def getRefTautomerName(self):
        return self.ref_tautomer.getName()

    def getRefTautomer(self):
        return self.ref_tautomer

    def getMolecule(self):
        return self.molecule

    def getType(self):
        return self.type

    def getRefProtState(self):
        if self.type == "c":
            reftau_prot_state = 1
        elif self.type == "a":
            reftau_prot_state = -1
        return reftau_prot_state

    def getCenterOriginal(self):
        if Config.delphi_params["pbc_dim"] != 2:
            raise Exception(
                "ERROR: The original center is only "
                "needed for 2-dimensional calculation"
            )
        return self.center_original

    def getCenterH(self):
        return self.centerH

    def getCenter(self):
        return self.center

    def getResNumber(self, correct_icode=False):
        if self.res_name in ("NTR", "CTR"):
            return self.res_number - TERMINAL_OFFSET
        elif correct_icode and self.res_number in self.molecule.icodes:
            origin_resnumb, icode = self.molecule.icodes[self.res_number]
            return "{0}{1}".format(origin_resnumb, icode)
        else:
            return self.res_number

    def getpK(self):
        return self.pK

    def getTautProb(self, taut, pH):
        taut_i = taut - 1
        taut_prob = self.getTautsProb(pH)[taut_i]
        ntauts = self.getNTautomers() + 1
        state_prob = 0.0
        if taut_i == ntauts - 1:
            state_prob = taut_prob
        else:
            i = 0
            while i < ntauts - 1:
                prob = self.getTautsProb(pH)[i]
                state_prob += prob
                i += 1

        return state_prob, taut_prob

    def getAverageProt(self, pH):
        """Calculates the average protonation of a site at a given pH value

        Args:
            site (str): Residue number of the site with the suffix 'N'
            or 'C' to indicate a N-ter or C-ter, respectively.

            pH (float): pH value

        Returns:
            A float of the average protonation of the site at the
            selected pH value
        """
        if not self.pK:
            return "pk Not In Range"
        average_prot = 10 ** (self.pK - pH) / (1 + 10 ** (self.pK - pH))
        return average_prot

    def getProtState(self, pH):
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
        state = "undefined"
        average_prot = self.getAverageProt(pH)

        if isinstance(average_prot, str):
            return state, average_prot

        if average_prot > 0.9:
            state = "protonated"
        elif average_prot < 0.1:
            state = "deprotonated"

        return state, average_prot

    def getTitrationCurve(self):
        return self.tit_curve

    def getFinalState(self, pH):
        return self.final_states[pH]

    def getMostProbTaut(self, pH):
        most_prob_taut = self.most_prob_states[pH]
        return most_prob_taut, self.getTautProb(most_prob_taut, pH)

    def getTautsProb(self, pH):
        ntauts = len(self.tautomers) + 1
        return self.states_prob[pH][:ntauts]

    # Iter Methods
    def iterTautomers(self):
        for i in list(self.tautomers.values()):
            yield i

    def iterOrderedTautomersWithoutRef(self):
        for i in sorted(self.tautomers.keys()):
            yield self.tautomers[i]

    def getOrderedTautomersList(self):
        tmp = []
        for i in sorted(self.tautomers.keys()):
            tmp.append(self.tautomers[i])
        tmp.append(self.ref_tautomer)
        return tmp

    # Calculations
    def calc_interaction_between(self, site2):
        """Calculates the pairwise interaction energies
        related to the two sites

        Args:
          site1 (Titsite)
          site2 (Titsite)

        Ensures:
          to_write (str): site interaction energies formatted
          to be written in .dat file
        """

        def convertIntoDatFormat(tau1, tau2, interaction):
            """Returns a .dat format interaction line

            Args:
              tau1 (Tautomer): first tautomer of the pair
              tau2 (Tautomer): second tautomer of the pair
              interaction (float): interaction energy of the tautomers pair

            Ensures:
              line (str): .dat formatted interaction line
            """
            site1_index = all_sites.index(tau1.site)
            tau1_index = tau1.name[-1]
            tau1_dat = "{0:>3} {1:>2}".format(site1_index, tau1_index)
            site2_index = all_sites.index(tau2.site)
            tau2_index = tau2.name[-1]
            tau2_dat = "{0:>3} {1:>2}".format(site2_index, tau2_index)
            line = "{0}   {1}   {2:13.6e}\n".format(tau1_dat, tau2_dat, interaction)
            return line

        site1 = self
        pairs = [(site1, site2), (site2, site1)]

        all_sites = Config.parallel_params.all_sites

        interactions_look = Config.parallel_params.interactions_look

        temperature = Config.pypka_params["temp"]
        interactions = []
        for pair in pairs:
            site1 = pair[0]
            site2 = pair[1]
            ordered_tautomers_site1 = site1.getOrderedTautomersList()
            ordered_tautomers_site2 = site2.getOrderedTautomersList()
            for taut1 in ordered_tautomers_site1:
                site_atom_list = taut1.site.getAtomNumbersList()
                for taut2 in ordered_tautomers_site2:
                    interaction = taut1.calcInteractionWith(taut2, site_atom_list)

                    gg = interaction / (KBOLTZ * temperature)

                    nsite1 = all_sites.index(taut1.site)
                    nsite2 = all_sites.index(taut2.site)

                    state1 = int(taut1.name[-1])
                    state2 = int(taut2.name[-1])

                    if Config.debug:
                        print(
                            (nsite1, nsite2, state1, state2, taut1.name, taut2.name, gg)
                        )
                    site1i = interactions_look[nsite1][state1]
                    site2i = interactions_look[nsite2][state2]

                    datf = convertIntoDatFormat(taut1, taut2, interaction)
                    interactions.append((site1i, site2i, gg, datf))

        return interactions

    # Print Methods
    def __str__(self):
        tautomers = " ".join(sorted(self.tautomers.keys()))
        natoms = len(self.atoms)
        center = [round(i, 2) for i in self.center]
        out = (
            "Site Number -> {0:5}   "
            "Tautomers -> {1:30}  "
            "Reference -> {2:5} "
            "NAtoms -> {3:6} "
            "Center -> {4}".format(
                self.res_number, tautomers, self.ref_tautomer.getName(), natoms, center
            )
        )
        return out
