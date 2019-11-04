from tautomer import Tautomer
import config


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
        self._res_name (str): residue name
        self._type (str): site type ('c' or 'a')
                          'c' if site is cationic
                          'a' if site is anionic
        self._center (list) : Site focus center
            if dimension is 0, then center is defined in cent file
            elif dimension is 2, then center (z) is defined in cent file
                                 and center (x, y) are defined in boxsize variable
        self._center_original (list): Site geometric center on
                                      cent file when dimension is 2

        # Site Tautomers
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
        self._atoms = {}
        self._tautomers = {}
        self._ref_tautomer = ''
        self._center = []
        self._type = ''
        self._pK = None

    # Set Methods
    def setTautomers(self, ntauts, resname):
        """Adds Tautomers from res_tauts to self._tautomers

        Args:
            res_tauts (list): tautomers from sites file
        """
        for i in range(ntauts):
            self._res_name = resname
            correct_number = str(i)
            tautomer = resname[:2] + correct_number

            tID = Tautomer(tautomer, self, self._molecule)
            self._tautomers[tautomer] = tID

    def setpK(self, pK):
        self._pK = pK

    def setTerminiResname(self, resname):
        self._termini_resname = resname
            
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
        for tautomer in list(self._tautomers.values()):
            tautomer.loadChargeSet(self._res_name,
                                   self._ref_tautomer)

    def addAtom(self, aname, anumb):
        """
        Args:
            aname (str): atom name
            anumb (int): atom id number
        """
        self._atoms[aname] = anumb

    def addCenter(self, center, boxsize=False, box_z=False):
        x = center[0]
        y = center[1]
        z = center[2]
        if boxsize:
            if config.params['pbc_dim'] != 2:
                raise Exception('ERROR: The original center is only '
                                'needed for 2-dimensional calculation')
            self._center_original = [x, y, z]
            x = boxsize * 10 / 2
            y = boxsize * 10 / 2
            z += box_z * 10
        self._center = [x, y, z]

    def addCenterH(self, center):
        self._centerH = center

    def setType(self, stype):
        self._type = stype

    # Get Methods
    def getName(self):
        return self._res_name

    def getTautomers(self):
        """Returns list of all tautomers instances
        except the tautomers of reference"""
        return list(self._tautomers.values())

    def getNTautomers(self):
        return len(self.getTautomers())

    def getAtomNumbersList(self):
        return list(self._atoms.values())

    def getAtomNamesList(self):
        return list(self._atoms.keys())

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
        if config.params['pbc_dim'] != 2:
            raise Exception('ERROR: The original center is only '
                            'needed for 2-dimensional calculation')
        return self._center_original

    def getCenterH(self):
        return self._centerH

    def getCenter(self):
        return self._center

    def getResNumber(self):
        return self._res_number

    def getpK(self):
        return self._pK

    # Iter Methods
    def iterTautomers(self):
        for i in list(self._tautomers.values()):
            yield i

    def iterOrderedTautomersWithoutRef(self):
        for i in sorted(self._tautomers.keys()):
            yield self._tautomers[i]

    def getOrderedTautomersList(self):
        tmp = []
        for i in sorted(self._tautomers.keys()):
            tmp.append(self._tautomers[i])
        tmp.append(self._ref_tautomer)
        return tmp

    # Print Methods
    def __str__(self):
        tautomers = ' '.join(sorted(self._tautomers.keys()))
        natoms = len(self._atoms)
        center = [round(i, 2) for i in self._center]
        out = 'Site Number -> {0:5}   '\
              'Tautomers -> {1:30}  '\
              'Reference -> {2:5} '\
              'NAtoms -> {3:6} '\
              'Center -> {4}'.format(self._res_number, tautomers,
                                     self._ref_tautomer.getName(),
                                     natoms, center)
        return out
