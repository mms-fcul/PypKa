import config
import numpy as np
from formats import read_pdb_line, read_gro_line, new_pdb_line
from copy import copy
from log import reportWarning
from titsite import Titsite as Site
from tautomer import Tautomer
from concurrency import startPoolProcesses, runInteractionCalcs, runMCCalcs
from mc import MCrun

MAXNPKHALFS = 5

class Molecule(object):
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
        self.box = []
        self._dat_to_write = ''
        self._correct_names = {}
        self._correct_atoms = {}

        self._NTR = None
        self._CTR = None

        self.readTermini()

    def readTermini(self):
        NTR_atoms = []
        with open("{0}/{1}/sts/NTRtau1.st".format(config.script_dir,
                                                  config.params['ffID'])) as f:
            for line in f:
                parts = line.split()
                if len(parts) > 1:
                    NTR_atoms.append(parts[1].strip())

        CTR_atoms = []
        with open("{0}/{1}/sts/CTRtau1.st".format(config.script_dir,
                                                  config.params['ffID'])) as f:
            for line in f:
                parts = line.split()
                if len(parts) > 1:
                    CTR_atoms.append(parts[1].strip())

        self._NTR_atoms = NTR_atoms
        self._CTR_atoms = CTR_atoms

    # Set Methods
    def addAtom(self, aname, anumb, position):
        self._atoms[anumb] = aname
        self._atoms_array_position[anumb] = position

    def loadParams(self, params):
        self._delphi_refparams = params

    def loadDelPhiParams(self, delphimol):
        if delphimol != 'reload':
            self._delphimol = delphimol
        else:
            delphimol = self._delphimol
        self.p_atpos  = copy(delphimol.get_atpos())
        self.p_rad3   = copy(delphimol.get_rad3())
        self.p_chrgv4 = copy(delphimol.get_chrgv4())
        self.atinf    = copy(delphimol.get_atinf())
        self.p_iatmed = copy(delphimol.get_iatmed())

    # Get Methods
    def getSites(self):
        return self._sites

    def getSitesOrdered(self):
        return self._sites_order

    def getDelPhi(self):
        return self._delphimol

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
        raise Exception('Something is wrong')

    def getTautomerInstance(self, tautname, site_resnum):
        """Return the tautomer instance named tautname
        existent in the site of the residue number site_resnum"""
        site = self._sites[site_resnum]
        if tautname in list(site._tautomers.keys()):
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
        for site in list(self._sites.values()):
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

    def iterAllSites(self):
        for site in self._sites_order:
            yield site

    # Printing Methods
    def printAllSites(self):
        """Prints all Site names"""
        for site in self._sites_order:
            print((site.getName()))

    def printAllTautomers(self):
        """Prints all Tautomer details"""
        for site in self._sites_order:
            print((site.getName()))
            for tautomer in site.iterTautomers():
                print(tautomer)
            print((site._ref_tautomer))

    # Input Files Manipulation Methods
    def makeSimpleSites(self):
        def SitesFileLine(resnumb, resname):
            for res in config.REGULARTITRATINGRES:
                if res[0:2] == resname[0:2]:
                    resname = res

            if resname in config.TITRABLETAUTOMERS:
                ntautomers = config.TITRABLETAUTOMERS[resname]
            else:
                for res in config.REGULARTITRATINGRES:
                    if res[0:2] == resname[0:2]:
                        ntautomers = config.TITRABLETAUTOMERS[res]
            text = '{0} '.format(resnumb)
            for i in range(ntautomers):
                text += '{0}tau{1} '.format(resname, i + 1)

            if resname in ('NTR', 'CTR'):
                resnumb += config.terminal_offset

            sID = self.addSite(resnumb)
            self.addTautomers(sID, ntautomers, resname)

            return text + '\n'

        def add2chain(chain, chain_res, resnumb, resname):
            if chain in chain_res:
                chain_res[chain][resnumb] = resname
            else:
                chain_res[chain] = {}
                chain_res[chain][resnumb] = resname

        sites_file = ''

        sites = []
        chain_res = {}
        with open(config.f_in) as f:
            nline = 0
            resnumb = None
            resname = None
            for line in f:
                nline += 1
                if 'ATOM ' == line[0:5]:
                    (aname, anumb, resname, chain,
                     resnumb, x, y, z) = read_pdb_line(line)
                    if chain == ' ':
                        chain = 'A'
                    if resnumb not in sites:
                        if len(sites) == 0 and 'NTR' not in sites and chain == 'A':
                            sites.append('NTR')
                            sites_file += SitesFileLine(resnumb, 'NTR')
                            self._NTR = resnumb
                        if resname in config.TITRABLERESIDUES and \
                           resname != 'NTR' and resname != 'CTR':
                            sites.append(resnumb)
                            sites_file += SitesFileLine(resnumb, resname)
                            add2chain(chain, chain_res, resnumb, resname)

                    if 'CTR' not in sites and \
                       aname in ('CT', 'OT', 'OT1', 'OT2', 'O1', 'O2', 'OXT') and chain == 'A':
                        sites.append('CTR')
                        sites_file += SitesFileLine(resnumb, 'CTR')
                        self._CTR = resnumb
                        add2chain(chain, chain_res, resnumb, 'CTR')
                        
        # Adding the reference tautomer to each site
        self.addReferenceTautomers()
        # Assigning a charge set to each tautomer
        self.addTautomersChargeSets()

        if config.debug:
            with open('tmp.sites', 'w') as f_new:
                f_new.write(sites_file)

        return chain_res, sites

    def makeSites(self, useTMPgro=None, sites=None):
        """Identifies titrable residues and checks integrity of the residue blocks
        (excluding Hydrogens)
        """
        def correctResName(resname):
            for res in config.REGULARTITRATINGRES:
                if res[0:2] == resname[0:2]:
                    return res
            return resname

        def makeSite(resnumb, resname):
            if resname in config.TITRABLETAUTOMERS:
                ntautomers = config.TITRABLETAUTOMERS[resname]
            else:
                for res in config.REGULARTITRATINGRES:
                    if res[0:2] == resname[0:2]:
                        ntautomers = config.TITRABLETAUTOMERS[res]
            sID = self.addSite(resnumb)
            self.addTautomers(sID, ntautomers, resname)

        def warning(resnumb, resname, res_atoms, mode=None):
            if mode == 'CYS' or resname == 'CYS':
                CYS_atoms = ['N', 'CA', 'CB', 'SG', 'C', 'O', 'H']
                if set(res_atoms).issubset(CYS_atoms) and \
                   set(CYS_atoms).issubset(res_atoms):
                    # no need to correct residue name
                    warn = '{0} {1} is assumed to be participating '\
                           'in a SS-bond'.format(resnumb, resname)
                    reportWarning(warn)
                    return

                CY0_atoms = ['N', 'CA', 'CB', 'SG', 'C', 'O', 'H', 'HG1']
                if set(res_atoms).issubset(CY0_atoms) and \
                   set(CY0_atoms).issubset(res_atoms):
                    self._correct_names[resnumb] = 'CY0'
                    return

                CY0_atoms = ['N', 'CA', 'CB', 'SG', 'C', 'O', 'H', 'HG']
                if set(res_atoms).issubset(CY0_atoms) and \
                   set(CY0_atoms).issubset(res_atoms):
                    self._correct_names[resnumb] = 'CY0'
                    self._correct_atoms[resnumb] = {'HG': 'HG1'}
                    return
                else:
                    warn = '{0} {1} failed integrity check'.format(resnumb,
                                                                   resname)
                    reportWarning(warn)

            elif resname not in config.TITRABLERESIDUES:
                return
            else:
                warn = '{0} {1} failed integrity check'.format(resnumb,
                                                               resname)
                reportWarning(warn)

        if config.f_in and not useTMPgro:
            filename = config.f_in
            filetype = 'pdb'
        else:
            filename = "TMP.gro"
            filetype = 'gro'

        resnumb = None
        cur_atoms = []
        prev_resnumb = None
        prev_resname = None
        with open(filename) as f:
            nline = 0
            maxnlines = 0
            for line in f:
                resname = None
                nline += 1
                if 'ATOM ' == line[0:5]:
                    (aname, anumb, resname, chain,
                     resnumb, x, y, z) = read_pdb_line(line)
                elif filetype == 'gro':
                    if nline > 2 and nline < maxnlines:
                        (aname, anumb, resname,
                         resnumb, x, y, z) = read_gro_line(line)
                    elif nline == 2:
                        natoms = int(line.strip())
                        maxnlines = natoms + 3

                if line == 'TER\n':
                    resnumb += 1
                if (prev_resnumb != resnumb or nline == maxnlines) and \
                   prev_resnumb:
                    if nline == maxnlines:
                        prev_resnumb = copy(resnumb)
                        resnumb = 'None'

                    if prev_resname in config.TITRABLERESIDUES or \
                       (prev_resnumb == self._NTR or resnumb == self._NTR) or \
                       (prev_resnumb == self._CTR or resnumb == self._CTR):
                        if prev_resnumb == self._NTR and resnumb != self._NTR:
                            prev_resname = correctResName(prev_resname)

                            if prev_resname in config.TITRABLERESIDUES:
                                res_tits = True
                            else:
                                res_tits = False

                            (integrity_nter,
                             integrity_site) = self.check_integrity(prev_resname,
                                                                    cur_atoms,
                                                                    nter=True,
                                                                    site=res_tits)
                            if integrity_nter:
                                nter_resnumb = prev_resnumb + config.terminal_offset
                                makeSite(nter_resnumb, 'NTR')
                                self._NTR = prev_resnumb
                            else:
                                warning(prev_resnumb, 'NTR', '')

                            if sites and prev_resnumb in sites:
                                if integrity_site:
                                    makeSite(prev_resnumb, prev_resname)
                                else:
                                    warning(prev_resnumb,
                                            prev_resname, cur_atoms)

                            prev_resnumb = None

                        # Dealing with the last residue and CTR
                        elif prev_resnumb == self._CTR and resnumb != self._CTR:
                            prev_resname = correctResName(prev_resname)

                            if prev_resname in config.TITRABLERESIDUES:
                                res_tits = True
                            else:
                                res_tits = False

                            integrity_cter, integrity_site = self.check_integrity(prev_resname,
                                                                                  cur_atoms,
                                                                                  cter=True,
                                                                                  site=res_tits)
                            if integrity_cter:
                                cter_resnumb = prev_resnumb + config.terminal_offset
                                makeSite(cter_resnumb, 'CTR')
                                self._CTR = prev_resnumb
                            else:
                                warning(prev_resnumb, 'CTR', '')

                            if sites and prev_resnumb in sites:
                                if integrity_site:
                                    makeSite(prev_resnumb, prev_resname)
                                else:
                                    warning(prev_resnumb, prev_resname, cur_atoms)

                        # Dealing with the previous residue
                        elif prev_resnumb:
                            if sites and prev_resnumb in sites:
                                prev_resname = correctResName(prev_resname)
                                res_atoms = copy(cur_atoms)
                                integrity_site = self.check_integrity(prev_resname,
                                                                      cur_atoms)
                                
                                if integrity_site:
                                    makeSite(prev_resnumb, prev_resname)
                                else:
                                    prev_resname = warning(prev_resnumb, prev_resname, res_atoms)
                            elif prev_resname == 'CYS':
                                # dealing with a CYS that is not in sites
                                res_atoms = copy(cur_atoms)
                                integrity_site = self.check_integrity(prev_resname,
                                                                      cur_atoms)

                                if not integrity_site:
                                    prev_resname = warning(prev_resnumb, prev_resname,
                                                           res_atoms, mode='CYS')
                    elif prev_resname == 'ALA':
                        #TODO: check residue block integrity for other non titrating residues
                        pass

                # Dealing with the new residue
                if prev_resnumb != resnumb:                    
                    cur_atoms = [aname]
                    prev_resnumb = resnumb
                    prev_resname = resname

                elif resnumb:
                    cur_atoms.append(aname)
                    if prev_resname in ('NTR', 'CTR') and \
                       prev_resname != resname:
                        prev_resname = resname

        # Adding the reference tautomer to each site
        self.addReferenceTautomers()
        # Assigning a charge set to each tautomer
        self.addTautomersChargeSets()

        # TODO: report blocks that failed the check (in .log file with
        # numbering reference to stepwise scheme)

        # TODO: add lipid residues

        if config.debug:
            print('exiting makeSites')

    def check_integrity(self, resname, res_atoms,
                        nter=False, cter=False, site=True):
        def read_anames(filename):
            res_atoms_st = []
            with open(filename) as f:
                for line in f.readlines()[1:]:
                    cols = line.strip().split()
                    aname = cols[1]
                    res_atoms_st.append(aname)
            return res_atoms_st

        def pop_atoms(anames1, anames2):
            trigger = False
            for aname in anames1:
                if aname not in anames2:
                    trigger = True
                    if config.debug:
                        print((aname, 'not in', resname))
                else:
                    anames2.remove(aname)
            if trigger:
                return anames2, False
            else:
                return anames2, True

        if config.debug:
            print('###### INTEGRITY CHECK ######')
            print(res_atoms)
        integrity = None

        if nter:
            st = '{0}/{1}/sts/NTRtau1.st'.format(config.script_dir,
                                                 config.params['ffID'])
            res_atoms_st = read_anames(st)
            res_atoms, integrity_nter = pop_atoms(res_atoms_st, res_atoms)
        elif cter:
            st = '{0}/{1}/sts/CTRtau1.st'.format(config.script_dir,
                                                 config.params['ffID'])
            res_atoms_st = read_anames(st)
            res_atoms, integrity_cter = pop_atoms(res_atoms_st, res_atoms)

        if site:
            main_chain = ('N', 'H', 'CA', 'C', 'O')
            if nter:
                main_chain = ()
            if cter:
                main_chain = ('N', 'H', 'CA')

            for aname in main_chain:
                if aname in res_atoms:
                    res_atoms.remove(aname)
                elif not nter and not cter:
                    integrity = False
            
            if config.debug:
                print(('i', integrity))

            if integrity is not False:
                filename = '{0}/{1}/sts/{2}tau1.st'.format(config.script_dir,
                                                           config.params['ffID'],
                                                           resname)
                res_atoms_st = read_anames(filename)
                res_atoms, integrity = pop_atoms(res_atoms_st, res_atoms)
                if config.debug:
                    print((res_atoms_st, res_atoms))
                    print(('i', integrity))

        if len(res_atoms) != 0 and integrity:
            raise Exception('Something is wrong')

        if nter:
            return integrity_nter, integrity
        elif cter:
            return integrity_cter, integrity
        else:
            return integrity

    def addSite(self, resnum):
        sID = Site(resnum, self)
        self._sites[resnum] = sID
        self._sites_order.append(sID)

        return sID

    def addTautomers(self, sID, ntautomers, resname):
        rootname = resname[0:2]
        sID._res_name = resname
        for itautomer in range(ntautomers):
            tautomer = rootname + str(itautomer)
            tID = Tautomer(tautomer, sID, sID._molecule)
            sID._tautomers[tautomer] = tID
        return sID._tautomers

    def loadSites(self, chains_length, chains_res):
        def siteHasNTautomers(chain, resnum):
            if resnum in chains_res[chain]:
                resname = chains_res[chain][resnum]
                if resname in config.TITRABLETAUTOMERS:
                    ntauts = config.TITRABLETAUTOMERS[resname]
                    return ntauts, resname
                else:
                    for res in config.TITRABLETAUTOMERS:
                        if res[:-1] == resname[:-1]:
                            ntauts = config.TITRABLETAUTOMERS[res]
                            return ntauts, res
            print(chains_length, chains_res, chain, resnum)
            print(type(resnum))
            print(resnum in chains_res[chain])
            print(resname in config.TITRABLETAUTOMERS)
            raise Exception('Something is wrong.')

        for chain in config.sites:
            sites = config.sites[chain]
            for site in sites:
                if 'N' in site:
                    resnum = int(site[:-1])
                    self._NTR = resnum
                    resnum += config.terminal_offset
                    res_ntauts = 3
                    res_name = 'NTR'
                elif 'C' in site:
                    resnum = int(site[:-1])
                    self._CTR = resnum
                    resnum += config.terminal_offset
                    res_ntauts = 4
                    res_name = 'CTR'
                else:
                    resnum = int(site)
                    res_ntauts, res_name = siteHasNTautomers(chain, resnum)

                sID = self.addSite(resnum)
                sID.setTautomers(res_ntauts, res_name)

        self.addReferenceTautomers()
        self.addTautomersChargeSets()


    def addReferenceTautomers(self):
        for site in list(self._sites.values()):
            site.addReferenceTautomer()

    def addTautomersChargeSets(self):
        for site in list(self._sites.values()):
            # reads .st files
            site.addChargeSets()

    def readIndexFile(self, f_ndx):
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
        for site in list(self._sites.values()):
            site._tautomers = {}
            site._ref_tautomer = ''
        self._sites = {}
        self._sites_order = []

    def readGROFile(self, groname):
        # TODO: For CpHMD read the index file:
        # self.readIndexFile()

        # Getting the "protein" with trjconv PBC mol
        # os.system('echo "Protein" | {GroDIR}/trjconv '
        #          '-f TMP_{config.sysname}.gro -s TMP_CpHMD.tpr '
        #          '-n TMP_CpHMD.ndx -o TMP_auxD1.gro '
        #          '-pbc mol -quiet'.format(GroDIR=paths['GroDIR'], sysname=config.sysname))

        new_pdb_content = """TITLE     Protein in water
REMARK    THIS IS A SIMULATION BOX
CRYST1  000.000  000.000  000.000  60.00  60.00  90.00 P 1           1
MODEL        1
"""
        new_pdb_content = ""
        site_positions = {}
        site_Hs = {}
        with open(groname) as f:
            line_counter = 0
            natoms_left = 0
            natoms = 0

            lines = f.read().splitlines()
            penul_line = lines[-2]
            last_line  = lines[-1]
            #final_res_numb = int(penul_line[:5].strip())
            #final_atom_numb = int(penul_line[15:20].strip())
            self.box = [float(i) for i in last_line.split()[:3]]

            if config.params['pbc_dim'] == 2:
                scaleP = (config.params['gsize'] - 1) / (self.box[0] * 10)
                scaleM = int(4 / scaleP + 0.5) * scaleP

                config.params['scaleP'] = scaleP
                config.params['scaleM'] = scaleM

            for line in lines:
                line_counter += 1
                if natoms_left > 0:

                    natoms_left -= 1
                    natoms += 1
                    aposition = natoms - 1

                    (aname, anumb, resname,
                     resnumb, x, y, z) = read_gro_line(line)
                    x, y, z = x * 10, y * 10, z * 10

                    #print resnumb, resname, aname, anumb, aposition, x, y, z

                    if resnumb == self._NTR and aname in self._NTR_atoms:
                        resname = 'NTR'
                        resnumb += config.terminal_offset
                    elif resnumb == self._CTR and aname in self._CTR_atoms:
                        resname = 'CTR'
                        resnumb += config.terminal_offset
                        if aname == 'C':
                            aname = 'CT'
                        #elif aname in ('O1', 'O2'):
                        #    aname = aname[0] + 'T' + aname[1]

                    if resnumb in list(self._correct_names.keys()):
                        resname = self._correct_names[resnumb]
                    if resnumb in list(self._correct_atoms.keys()) and \
                       aname in self._correct_atoms[resnumb]:
                        aname = self._correct_atoms[resnumb][aname]

                    self.addAtom(aname, anumb, aposition)
                    ref_tau_name = resname

                    if resnumb in list(self._sites.keys()) and \
                       aname in list(self._sites[resnumb]._ref_tautomer._charge_set.keys()):
                        #( aname not in ('N', 'H', 'C', 'O', 'CA') or 
                        #(aname in ('N', 'H', 'C', 'O', 'CA') and resname == 'NTR')):
                        # change res name to reference tautomer
                        ref_tau_name = self._sites[resnumb].getRefTautomerName()

                        # add atom to corresponding site
                        self._sites[resnumb].addAtom(aname, anumb)
                        if resnumb in site_positions:
                            site_positions[resnumb].append((x, y, z))
                            if aname[0] == 'H':
                                site_Hs[resnumb].append((x, y, z))
                        else:
                            site_positions[resnumb] = [(x, y, z)]
                            if aname[0] == 'H':
                                site_Hs[resnumb] = [(x, y, z)]
                            else:
                                site_Hs[resnumb] = []

                    new_pdb_content += new_pdb_line(aposition, aname,
                                                    ref_tau_name, resnumb,
                                                    x, y, z)

                elif line_counter == 2:
                    natoms_left = int(line.strip())
                    self._natoms = natoms_left

        new_pdb_content += 'TER\nENDMDL\n'
        with open('delphi_in_stmod.pdb', 'w') as f_new:
            f_new.write(new_pdb_content)

        # TODO: check Terminal_offset has to be bigger than the total number of residues
        # TODO: delete terminal_offset and use another approach to distinguish between N- and C-ter
        # TODO: check size xy > config.cutoff * 2
        # if so, raise Exception, and ask to change cutoff value

        # TODO: check if pbc_dim -> set gsizes from pdb size xy and ignore perfil

        for site in site_positions:
            if site in list(self._sites.keys()):
                pos_max = [-9999990, -999999, -999999]
                pos_min = [999999, 999999, 999999]
                focus_center = [0, 0, 0]
                for atom in site_positions[site]:
                    for i in range(3):
                        if pos_max[i] < atom[i]:
                            pos_max[i] = atom[i]
                        if pos_min[i] > atom[i]:
                            pos_min[i] = atom[i]
                focus_center[0] = (pos_max[0] + pos_min[0]) / 2
                focus_center[1] = (pos_max[1] + pos_min[1]) / 2
                focus_center[2] = (pos_max[2] + pos_min[2]) / 2
                if config.params['pbc_dim'] == 2:
                    self._sites[site].addCenter(focus_center,
                                                boxsize=self.box[0],
                                                box_z=self.box[2])
                else:
                    self._sites[site].addCenter(focus_center)
                hx, hy, hz = 0, 0, 0
                nHs = len(site_Hs[site])
                if nHs == 0:
                    sitename = self._sites[site].getName()
                    raise Exception('Site {1}{0} appears '
                                    'to have no Hydrogen atoms'.format(site,
                                                                       sitename))
                for h in site_Hs[site]:
                    hx += h[0]
                    hy += h[1]
                    hz += h[2]
                hx /= nHs
                hy /= nHs
                hz /= nHs
                Hcenter = (round(hx, 2), round(hy, 2), round(hz, 2))
                self._sites[site].addCenterH(Hcenter)

        if config.debug:
            with open('cent', 'w') as f_new:
                text = ''
                for site in site_positions:
                    if site in list(self._sites.keys()):
                        text += str(self._sites[site]._center) + '\n'
                f_new.write(text)
            if config.params['pbc_dim'] == 2:
                with open('centHs', 'w') as f_new:
                    text = ''
                    for site in site_positions:
                        if site in list(self._sites.keys()):
                            text += '{0} {1}\n'.format(site,
                                                       self._sites[site]._centerH)
                    f_new.write(text)
                with open('cent_original', 'w') as f_new:
                    text = ''
                    for site in site_positions:
                        if site in list(self._sites.keys()):
                            text += '{0} {1}\n'.format(site,
                                                       self._sites[site]._center_original)
                    f_new.write(text)

    def calcSiteInteractionsParallel(self, ncpus):
        """Calculates the pairwise interaction energies
        and writes them in a formatted .dat file
        Interactions are calculated using a pool of processes

        Args:
          ncpus (int): number of cpus to be used
        """
        if config.debug:
            self.writeDatHeader()
        counter = 0
        for site1 in self.getSitesOrdered()[:-1]:
            counter += 1
            for site2 in self.getSitesOrdered()[counter:]:
                self._site_interactions.append((site1, site2))

        sites = self.getSitesOrdered()
        self._npossible_states = [len(site._tautomers) + 1 for site in sites]

        self._interactions = []
        for nstate1 in range(sum(self._npossible_states)):
            self._interactions.append([])
            for nstate2 in range(sum(self._npossible_states)):
                self._interactions[nstate1].append(-999999)

        self._interactions_look = []
        aux = -1
        site = -1
        for nstates in self._npossible_states:
            site += 1
            self._interactions_look.append([])
            for state in range(nstates):
                aux += 1
                self._interactions_look[site].append(aux)

        ncpus = min(len(self._site_interactions), ncpus)
        if ncpus > 0:
            results = startPoolProcesses(runInteractionCalcs,
                                         self._site_interactions, ncpus,
                                         assign='ordered', merged_results=True)
        else:
            results = []
        to_write = ''
        T = float(config.params['temp'])
        for interaction in results:
            site1i = interaction[0]
            site2i = interaction[1]
            if self._interactions[site1i][site2i] != -999999:
                self._interactions[site1i][site2i] += interaction[2]
                self._interactions[site2i][site1i] += interaction[2]
                self._interactions[site1i][site2i] /= 2
                self._interactions[site2i][site1i] /= 2

                if config.debug:
                    col1 = interaction[3][6:18]
                    col2 = interaction[3][:6]
                    col3 = self._interactions[site1i][site2i] * (config.kBoltz * T)
                    to_write += '{0}{1} {2:13.6e}\n'.format(col1, col2, col3)
            else:
                self._interactions[site1i][site2i] = interaction[2]
                self._interactions[site2i][site1i] = interaction[2]

        if config.debug:
            with open('interactions.dat', 'a') as f_new:
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
        sites = self._site_interactions[inumber]
        site1 = sites[0]
        site2 = sites[1]
        pairs = ((site1, site2), (site2, site1))

        iterAtomsList = self.getAtomsList()

        T = float(config.params['temp'])
        interactions = []
        for pair in pairs:
            site1 = pair[0]
            site2 = pair[1]
            ordered_tautomers_site1 = site1.getOrderedTautomersList()
            ordered_tautomers_site2 = site2.getOrderedTautomersList()
            for taut1 in ordered_tautomers_site1:
                site_atom_list = taut1._site.getAtomNumbersList()
                for taut2 in ordered_tautomers_site2:
                    interaction = taut1.calcInteractionWith(taut2,
                                                            site_atom_list,
                                                            iterAtomsList)

                    gg = interaction / (config.kBoltz * T)

                    nsite1 = self.getSitesOrdered().index(taut1._site)
                    nsite2 = self.getSitesOrdered().index(taut2._site)

                    state1 = int(taut1._name[-1])
                    state2 = int(taut2._name[-1])

                    if config.debug:
                        print((nsite1, nsite2, state1, state2,
                               taut1._name, taut2._name, gg))
                    site1i = self._interactions_look[nsite1][state1]
                    site2i = self._interactions_look[nsite2][state2]

                    datf = self.convertIntoDatFormat(taut1, taut2, interaction)
                    interactions.append((site1i, site2i, gg, datf))

        return interactions

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
                to_write += '{0:1d} {1:13.6e}\n'.format(tau_prot_state,
                                                        tautomer._dg)
            to_write += '{0:1d}  0.000000e+00\n'.format(ref_prot_state)
        if config.debug:
            with open('interactions.dat', 'w') as f_new:
                f_new.write(to_write)

    def calcpKint(self, unpacked_results):
        """Calculation the pKint of all tautomers
        """
        i = -1
        if config.debug:
            print('############ results ############')
        pkints = ''
        contributions = ''
        for tautomer in config.tit_mole.iterAllSitesTautomers():
            i += 1
            core_index = i % config.params['ncpus']
            job_index = int(i / config.params['ncpus'])
            result = unpacked_results[core_index][job_index]

            tautname    = result[0]
            tautresnumb = result[1]

            esolvationM = result[2]
            sitpotM     = result[3]
            esolvationS = result[4]
            sitpotS     = result[5]

            tautomer = config.tit_mole.getTautomerInstance(tautname,
                                                           tautresnumb)

            tautomer.saveDelPhiResults(esolvationS, sitpotS, esolvationM,
                                       sitpotM)
            if config.debug:
                print('### new tautomer ###')
                print((i, core_index, job_index, tautname,
                       tautomer._name, esolvationM, esolvationS))
                print((tautomer._name, tautomer._esolvationS,
                       len(tautomer._sitpotS), tautomer._esolvationM,
                       len(tautomer._sitpotM)))

            tautomer.calcBackEnergy()
            if not tautomer.isRefTautomer():
                tautomer.calcpKint()
                if config.debug:
                    print(('pkint', tautomer._name, tautomer._dg,
                           id(tautomer)))
                    pkints += '{} {} {}\n'.format(tautomer._name,
                                                  tautomer._dg, tautomer.pKint)
                    contributions += '{}{} {} {} {} {}\n'.format(tautomer._site._res_number,
                                                                 tautomer._name,
                                                                 tautomer.dG_solvationM,
                                                                 tautomer.dG_solvationS,
                                                                 tautomer.dG_solvationM - tautomer.dG_solvationS,
                                                                 tautomer.dG_back)
        if config.debug:
            with open('pkint', 'w') as f_new1, \
                 open('contributions', 'w') as f_new2:
                f_new1.write(pkints)
                f_new2.write(contributions)

    def calcpKhalfs(self, pH, nsites, avgs, pmean, pKs, count, mcsteps, pHmin, dpH):
        totalP = 0.0
        for site in range(nsites):
            mean = avgs[site] / float(mcsteps)
            totalP += mean
        
            p = pmean[site]
        
            if p > 0.5 and mean <= 0.5 or p < 0.5 and mean >= 0.5:
                pKhalf = pH - dpH * (mean - 0.5) / (mean - p)
                for i in range(MAXNPKHALFS): 
                    if pKs[site][i] == 100.0:
                        pKs[site][i] = pKhalf
                        break
        
            elif p > 1.5 and mean <= 1.5 or p < 1.5 and mean >= 1.5:
                pKhalf = pH - dpH * (mean - 1.5) / (mean - p)
                for i in range(MAXNPKHALFS): 
                    if pKs[site][i] == 100.0:
                        pKs[site][i] = pKhalf
                        break
        
            pmean[site] = mean
        return totalP, pKs, pmean

                
    def parallelMCrun(self, pHstep):
        sites = self.getSitesOrdered()
        nsites = len(sites)

        pHmin, pHmax = config.params['pHmin'], config.params['pHmax']
        dpH = config.params['pHstep']
        couple_min = config.params['couple_min']
        mcsteps = config.params['mcsteps']
        eqsteps = config.params['eqsteps']
        seed = config.params['seed']

        pHsteps = int(round(1 + (pHmax - pHmin) / dpH, 0))

        pH = pHmin + pHstep * dpH
        avgs, pmean, count = MCrun(nsites, self._npossible_states,
                                   self._possible_states_g,
                                   self._possible_states_occ,
                                   self._interactions,
                                   self._interactions_look,
                                   pHsteps, mcsteps, eqsteps, seed,
                                   pHmin, dpH, couple_min, pH)

        return (avgs, count)


    def runMC(self):
        def resize_list_of_lists(listn, maxsize, filler=None):
            for i in listn:
                diff = maxsize - len(i)
                for ii in range(diff):
                    i.append(filler)

        print('\nStart MC')

        sites = self.getSitesOrdered()
        nsites = len(sites)
        self._possible_states     = [[] for site in sites]
        self._possible_states_g   = [[] for site in sites]
        self._possible_states_occ = [[] for site in sites]

        T = float(config.params['temp'])
        isite = -1
        for site in sites:
            isite += 1
            itaut = 0
            for tautomer in site.iterOrderedTautomersWithoutRef():
                dg = tautomer._dg / (config.kBoltz * T)
                self._possible_states_g[isite].append(dg)
                self._possible_states[isite].append(itaut)
                if site._type == 'c':
                    prot_state = 0
                elif site._type == 'a':
                    prot_state = 1

                self._possible_states_occ[isite].append(prot_state)
                itaut += 1

            if site._type == 'c':
                prot_state = 1
            elif site._type == 'a':
                prot_state = 0
            self._possible_states_occ[isite].append(prot_state)
            self._possible_states_g[isite].append(0.0)
        
        maxstates = max(self._npossible_states)
        resize_list_of_lists(self._possible_states, maxstates)
        resize_list_of_lists(self._possible_states_g, maxstates)
        resize_list_of_lists(self._possible_states_occ, maxstates, filler=-500)
        resize_list_of_lists(self._interactions_look, maxstates, filler=-500)

        pHmin, pHmax = config.params['pHmin'], config.params['pHmax']
        dpH = config.params['pHstep']
        couple_min = config.params['couple_min']
        mcsteps = config.params['mcsteps']
        eqsteps = config.params['eqsteps']
        seed = config.params['seed']

        pHsteps = int(round(1 + (pHmax - pHmin) / dpH, 0))

        pmeans_all = []
        counts_all = []
        avgs_all = []

        ncpus = min(config.params['ncpus'], nsites)
        results = startPoolProcesses(runMCCalcs,
                                     list(range(pHsteps)), ncpus,
                                     assign='ordered', merged_results=True)
        counter = 0
        for i in results:
            avgs_all.append(i[0])
            counts_all.append(i[1])
        
        pKs = np.array([[100.0 for ii in range(MAXNPKHALFS)] for i in range(nsites)])
        pmeans = avgs_all[0] / float(mcsteps)
        pmeans_raw = [100.0]
        for pHstep in range(1, pHsteps):
            pH = pHmin + pHstep * dpH
            totalP, pKs, pmeans = self.calcpKhalfs(pH, nsites, avgs_all[pHstep], pmeans, pKs, counts_all[pHstep], mcsteps, pHmin, dpH)
            pmeans_raw.append(totalP)

        pKas = pKs
        
        #pKas, pmeans_raw = mc.MCrun(nsites, self._npossible_states,
        #                            self._possible_states_g,
        #                            self._possible_states_occ,
        #                            self._interactions,
        #                            self._interactions_look,
        #                            pHsteps, mcsteps, eqsteps, seed,
        #                            pHmin, dpH, couple_min)

        print('\n\nResults')
        text_pks = ''
        text_prots = '# pH      total'
        c = -1
        for i in pKas:
            c += 1
            site = self.getSitesOrdered()[c]
            sitename = site.getName()
            if sitename in ('NTR', 'CTR'):
                resnumb = site.getResNumber() - config.terminal_offset
                text_prots += '     {0:3}'.format(sitename)
            else:
                resnumb = site.getResNumber()
                text_prots += '{0:5d}{1:3s}'.format(resnumb, sitename)
            text_pks += '{0} {1} {2}\n'.format(resnumb, sitename, i[0])
            if i[0] != 100.0:
                print(('{0} {1} {2}'.format(resnumb, sitename, round(i[0], 2))))
            else:
                print(('{0} {1} Not in range'.format(resnumb, sitename)))
        print('')

        if config.f_out:
            with open(config.f_out, 'w') as f:
                f.write(text_pks)

        pmeans = {}
        for pHstep in range(pHsteps):
            pH = pHmin + pHstep * dpH
            c = -1
            pmeans[pH] = {}
            text_prots += '\n{0:5.2f}'.format(pH)

            # total protonation
            #pmean = pmeans_raw[pHstep][c]
            pmean = pmeans_raw[pHstep]
            pmeans[pH]['total'] = pmean
            text_prots += '\t{0:7.4f}'.format(pmean)

            for site in self.getSitesOrdered():
                c += 1
                sitename = site.getName()
                sitenumb = site.getResNumber()
                #pmean = pmeans_raw[pHstep][c]
                pmean = avgs_all[pHstep][c]
                pmeans[pH][sitenumb] = pmean
                if pmean != 100.0:
                    text_prots += '\t{0:7.4f}'.format(pmean)
                else:
                    text_prots += '\t-'

        if config.f_prot_out:
            with open(config.f_prot_out, 'w') as f:
                f.write(text_prots)

        return pKas, pmeans
