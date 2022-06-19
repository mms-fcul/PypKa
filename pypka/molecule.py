from pypka.concurrency import runInteractionCalcs
from pypka.constants import TERMINAL_OFFSET, TITRABLETAUTOMERS
from pypka.tautomer import Tautomer
from pypka.titsite import Titsite as Site
from pdbmender.formats import read_gro_line


class Molecule:
    """Molecule with more than one titrable sites"""

    def __init__(self, chain, sites, icodes):
        """
        #TODO: update docstring

        # DelPhi Parameters
        self.delphi_refparams (list): input DelPhi parameters

        # Sites
        self.sites (dict): key is residue number
                           value is site object instance

        self.sites_order (list): values are site object instances
                                  order respects .sites

        # Atoms
        self.natoms (int): number of atoms in the molecule
        self.atoms (dict): key is atom id number
                            value is atom name
        self.atoms_array_position (dict):
            key is atom number
            value is atom position in pqr
            Example: if the first atom in pqr is number 541
                     then self.atoms_array_position[541] = 1

        # Offset of the first site
        This is only used because the membrane part is not finished
        and uses a pqr input of the first site
        self.InputPQRoffSet (list):
            x, y, z coordinates of the offset in the input pqr

        # Site Interactions
        self.site_interactions (list):
            contains tuples with pairs of site objects
        """
        self.input_sites = sites
        self.chain = chain
        self.icodes = icodes

        self.delphi_refparams = ""
        self.sites = {}
        self.sites_order = []

        self.natoms = 0
        self.atoms = {}
        self.atoms_array_position = {}
        self.atoms_tit_res = {}

        self.InputPQRoffSet = []
        self.box = []
        self.dat_to_write = ""

        self.correct_names = {}
        self.correct_atoms = {}

        self.CYS_bridges = {}

        self.NTR = []
        self.CTR = []

    # Set Methods
    def addAtom(self, aname, anumb, position, titrable_res):
        self.atoms[anumb] = aname
        self.atoms_array_position[anumb] = position
        self.natoms += 1
        self.atoms_tit_res[position] = titrable_res

    def loadParams(self, params):
        self.delphi_refparams = params

    def saveCYSBridges(self, CYS_bridges):
        self.CYS_bridges = CYS_bridges

    def saveHISStates(self, HIS_states):
        self.HIS_states = HIS_states

    # Get Methods
    def __getitem__(self, numb):
        numb = self.correct_site_numb(numb)
        return self.sites[numb]

    def getSites(self):
        return self.sites

    def getSitesOrdered(self):
        return self.sites_order

    def getArrayPosition(self, atom_id):
        """Returns the index of the atom in DelPhi position array

        Disclamer: used only in testing
        """
        return self.atoms_array_position[atom_id]

    def getAtomsList(self):
        """Returns a list of atoms details

        Atoms details = (id, instance, atom index in DelPhi data structures)
        sorted by atom id number
        """
        atom_list = []
        for atom in sorted(self.atoms.keys()):
            atom_list.append((self.atoms[atom], atom, self.atoms_array_position[atom]))
        return atom_list

    def getTautomerInstance(self, tautname, site_resnum):
        """Return the tautomer instance named tautname

        Tautomer instance must exist in the site of the residue number site_resnum
        """
        site = self.sites[site_resnum]
        if tautname in list(site.tautomers.keys()):
            return site.tautomers[tautname]
        elif tautname == site.ref_tautomer.name:
            return site.ref_tautomer
        raise Exception("Something is very wrong!!!")

    def getNAtoms(self):
        """Return number of atoms of the Site."""
        return self.natoms

    def getTautNAtoms(self, taut_name):
        """Return number of atoms of the tautomer named taut_name

        Disclamer: it was used for testing
        """
        for site in list(self.sites.values()):
            for tautomer in site.getTautomers():
                if tautomer.name == taut_name:
                    return tautomer.natoms

    def getPQROffset(self):
        return self.InputPQRoffSet

    def getCYS_bridges(self):
        return self.CYS_bridges

    # Iter Methods
    def iterAtoms(self):
        """Generator that iterates through all atoms details (name, id, position) in the Site."""
        for atom in sorted(self.atoms.keys()):
            yield (self.atoms[atom], atom, self.atoms_array_position[atom])

    def iterNonRefSitesTautomers(self):
        """Generator that iterates through all Tautomer instances except the reference ones."""
        for site in self.sites_order:
            for tautomer in sorted(site.tautomers.values()):
                yield tautomer

    def iterAllSites(self):
        for site in self.sites_order:
            yield site

    # Printing Methods
    def printAllSites(self):
        """Prints all Site names."""
        for site in self.sites_order:
            print((site.getName()))

    def printAllTautomers(self):
        """Prints all Tautomer details"""
        for site in self.sites_order:
            print((site.getName()))
            for tautomer in site.iterTautomers():
                print(tautomer)
            print((site.ref_tautomer))

    def addSite(self, resnum):
        sID = Site(resnum, self)
        self.sites[resnum] = sID
        self.sites_order.append(sID)

        return sID

    @staticmethod
    def addTautomers(sID, ntautomers, resname, termini_resname=None):
        rootname = resname[0:2]
        sID.res_name = resname
        if termini_resname:
            sID.setTerminiResname(termini_resname)
        for itautomer in range(ntautomers):
            tautomer = rootname + str(itautomer)
            tID = Tautomer(tautomer, sID, sID.molecule)
            sID.tautomers[tautomer] = tID
        return sID.tautomers

    def reorder_sites(self):
        new_list = []
        ctr = None
        nsites = len(self.sites_order)
        for i, site in enumerate(self.sites_order):
            if site.res_name == "CTR" and i != nsites - 1:
                ctr = site
                continue
            new_list.append(site)
        if ctr:
            new_list.append(ctr)
        self.sites_order = new_list

    def loadSites(self, sites):
        def siteHasNTautomers(resnum):
            if resnum in sites:
                resname = sites[resnum]
                if resname in TITRABLETAUTOMERS:
                    ntauts = TITRABLETAUTOMERS[resname]
                    return ntauts, resname
                else:
                    for res in TITRABLETAUTOMERS:
                        if res[:-1] == resname[:-1]:
                            ntauts = TITRABLETAUTOMERS[res]
                            return ntauts, res
            raise Exception("{0}_{1} is not a titrable residue".format(resname, resnum))

        for site_number, site_name in sites.items():
            if site_name == "NTR":
                resnum = int(site_number)
                self.NTR += [resnum]
                resnum += TERMINAL_OFFSET
                res_ntauts = 3
                res_name = "NTR"
            elif site_name == "CTR":
                resnum = int(site_number)
                self.CTR += [resnum]
                resnum += TERMINAL_OFFSET
                res_ntauts = 4
                res_name = "CTR"
            else:
                resnum = int(site_number)
                res_ntauts, res_name = siteHasNTautomers(resnum)

            sID = self.addSite(resnum)
            sID.setTautomers(res_ntauts, res_name)

        self.reorder_sites()

        self.addReferenceTautomers()
        self.addTautomersChargeSets()

    def addReferenceTautomers(self):
        for site in list(self.sites.values()):
            site.addReferenceTautomer()

    def addTautomersChargeSets(self):
        for site in list(self.sites.values()):
            # reads .st files
            site.addChargeSets()

    @staticmethod
    def readIndexFile(f_ndx):
        protein_trigger = False
        protein_atoms = []
        with open(f_ndx) as f:
            for line in f:
                if "[ protein ]" in line.lower():
                    protein_trigger = True
                elif protein_trigger:
                    for atom in line.split():
                        protein_atoms.append(int(atom))
        return protein_atoms

    def deleteAllSites(self):
        for site in list(self.sites.values()):
            site.tautomers = {}
            site.ref_tautomer = ""
        self.sites = {}
        self.sites_order = []

    def correct_site_numb(self, numb):
        if numb == "NTR":
            numb = self.NTR + TERMINAL_OFFSET
        elif numb == "CTR":
            numb = self.CTR + TERMINAL_OFFSET
        if isinstance(numb, str):
            try:
                numb = int(numb)
            except ValueError:
                raise Exception("Unknown site")
        return numb
