import config
import time
import log
from copy import copy
from formats import new_pdb_line


class Tautomer(object):
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

    def __lt__(self, other):
        return self._name < other._name

    # Set Methods
    def loadChargeSet(self, res_name, ref_tautomer):
        """Reads .st file related to Tautomer with residue name res_name
        Stores charge set related to both the Tautomer and the
        Reference Tautomer

        .st file named TYRtau1.st has charge set of TY0 and reference TY2
                       TYRtau2.st has charge set of TY1 and reference TY2
        """
        tau_number = int(self._name[-1]) + 1
        fname = '{0}/{1}/sts/{2}tau{3}.st'.format(config.script_dir,
                                                  config.params['ffID'],
                                                  res_name, tau_number)
        with open(fname) as f:
            nline = 0
            charge_set1 = {}
            charge_set2 = {}
            for line in f:
                nline += 1
                cols = line.split()
                if nline > 1:
                    atom_name = cols[1]
                    if atom_name == 'C' and res_name == 'CTR':
                        atom_name = 'CT'
                    #elif atom_name in ('O1', 'O2') and res_name == 'CTR':
                    #    atom_name = atom_name[0] + 'T' + atom_name[1]
                    charge_set1[atom_name] = float(cols[-2])
                    charge_set2[atom_name] = float(cols[-1])
                else:
                    self._pKmod = float(line)

        if sum(charge_set1.values()) < 0.001 and \
           sum(charge_set2.values()) > -1.001:
            if config.debug:
                print((fname, 'anionic'))
            self._charge_set = charge_set1
            ref_tautomer._charge_set = charge_set2
            self._site.setType('a')
        elif (sum(charge_set1.values()) > 0.99 and
              sum(charge_set2.values()) < 0.001):
            if config.debug:
                print((fname, 'cationic'))
            self._charge_set = charge_set2
            ref_tautomer._charge_set = charge_set1
            self._site.setType('c')
        self._natoms = len(self._charge_set)
        ref_tautomer._natoms = len(self._charge_set)
        if config.debug:
            print((self._charge_set))
            print((ref_tautomer._charge_set))
            print(('finished reading', fname))

    def saveDelPhiResults(self, esolvationS, sitpotS, esolvationM, sitpotM):
        self._esolvationS  =  esolvationS
        self._sitpotS      =  sitpotS
        self._esolvationM  =  esolvationM
        self._sitpotM      =  sitpotM

    def getSiteAcent(self):
        if config.params['pbc_dim'] == 2:
            x = self._molecule.box[0] * 10 / 2
            y = x
            z = self._site._center[2]
            acent = [x, y, z]
        else:
            acent = self._site._center
        acent = [round(acent[0], 4), round(acent[1], 4), round(acent[2], 3)]
        return acent

    def getMoleculeDetails(self):
        molecule = self._molecule
        delphimol = molecule.getDelPhi()

        # Could have used empty structures
        # TODO: quantify time impact of this quick fix
        p_atpos  = copy(molecule.p_atpos)
        p_rad3   = copy(molecule.p_rad3)
        p_chrgv4 = copy(molecule.p_chrgv4)
        atinf    = copy(molecule.atinf)
        p_iatmed = copy(molecule.p_iatmed)

        return molecule, delphimol, p_atpos, p_rad3, p_chrgv4, atinf, p_iatmed

    def getCenteringDetails(self):
        box = self._molecule.box
        half_box_xy = self._site.getCenter()[0]
        site_center = self._site.getCenterOriginal()

        offset_x = half_box_xy - site_center[0]
        offset_y = half_box_xy - site_center[1]
        offset_z = box[2] * 10  # irrelevant, only for debug

        box_x = box[0] * 10
        box_y = box[1] * 10

        return (box, half_box_xy, site_center,
                offset_x, offset_y, offset_z, box_x, box_y)

    def putInsideBox(self, point, box_x, box_y):
        if point[0] < 0:
            point[0] += box_x
        elif point[0] >= box_x:
            point[0] -= box_x
        if point[1] < 0:
            point[1] += box_y
        elif point[1] >= box_y:
            point[1] -= box_y
        return point

    def setDetailsFromWholeMolecule(self):
        """Set DelPhi parameters to run a calculation of a whole molecule
        (all sites neutral, except one)"""
        (molecule, delphimol, p_atpos,
         p_rad3, p_chrgv4, atinf, p_iatmed) = self.getMoleculeDetails()
        if config.params['pbc_dim'] == 2:
            (box, half_box_xy, site_center,
             offset_x, offset_y, offset_z,
             box_x, box_y) = self.getCenteringDetails()
        pdb_text = ''
        crg = '!crg file created by gen_files.awk\n'\
              'atom__resnumbc_charge_\n'
        new_atoms = []
        for atom_name, atom_id, atom_position in molecule.iterAtoms():
            aID = atom_position + 1
            aname = str(atinf[atom_position].value.split()[0])
            resname = str(atinf[atom_position].value.split()[1])
            resnumb = int(atinf[atom_position].value.split()[2])
            x = float(p_atpos[atom_position][0])
            y = float(p_atpos[atom_position][1])
            z = float(p_atpos[atom_position][2])
            pdb_text += new_pdb_line(aID, aname, resname, resnumb, x, y, z)

            # TODO: quick fix, should be done only once per site
            if config.params['pbc_dim'] == 2:
                p_atpos[atom_position][0] += offset_x
                p_atpos[atom_position][1] += offset_y
                p_atpos[atom_position][2] += offset_z

                p_atpos[atom_position] = self.putInsideBox(p_atpos[atom_position],
                                                           box_x, box_y)

                # Comment 15 May 2019
                # No rounding is better, however, pypka calculations
                # are no longer comparable to delphiT ones due to this #noregrets
                #rounding_precision = decimal.Decimal('0.01')
                #x_dec = decimal.Decimal(p_atpos[atom_position][0])
                #y_dec = decimal.Decimal(p_atpos[atom_position][1])
                #z_dec = decimal.Decimal(p_atpos[atom_position][2])
                #p_atpos[atom_position][0] = float(x_dec.quantize(rounding_precision))
                #p_atpos[atom_position][1] = float(y_dec.quantize(rounding_precision))
                #p_atpos[atom_position][2] = float(z_dec.quantize(rounding_precision))

                pbc_atoms = self.add_pbc(p_atpos[atom_position][0],
                                         p_atpos[atom_position][1],
                                         p_atpos[atom_position][2], box[0],
                                         p_rad3[atom_position],
                                         p_chrgv4[atom_position],
                                         atinf[atom_position].value)
                #pbc_atoms = []
                new_atoms += pbc_atoms

            if atom_id in self._site.getAtomNumbersList():
                p_chrgv4[atom_position] = self.getCharge(atom_name)
            else:
                p_chrgv4[atom_position] = 0.0

            if float(p_chrgv4[atom_position]) < 0.0:
                signal = '-'
            else:
                signal = ' '

            crg += '{0:<6}{1:<7}  {3}{2:0<5} \n'.format(aname, resname,
                                                        round(abs(p_chrgv4[atom_position]), 3),
                                                        signal)

        natoms = delphimol.changeStructureSize(molecule._natoms,
                                               p_atpos, p_rad3,
                                               p_chrgv4, atinf,
                                               p_iatmed,
                                               extra_atoms=new_atoms)
        acent = self.getSiteAcent()

        pdb_text = ''
        for i in range(natoms):
            aID = i + 1
            aname   = str(molecule._delphimol.atinf[i].value.split()[0])
            resname = str(molecule._delphimol.atinf[i].value.split()[1])
            resnumb = int(molecule._delphimol.atinf[i].value.split()[2])
            if '-' in aname[0]:
                aname = aname[1:]
                resname = 'PBC'
                resnumb = -666
            x = float(molecule._delphimol.p_atpos[i][0])
            y = float(molecule._delphimol.p_atpos[i][1])
            z = float(molecule._delphimol.p_atpos[i][2])
            pdb_text += new_pdb_line(aID, aname, resname, resnumb, x, y, z)

        box = molecule.box
        if config.debug:
            filename = 'P_{1}-{0}.pdb'.format(self._name, self._site._res_number)
            with open(filename, 'w') as f_new:
                x, y, z = acent
                pdb_text += new_pdb_line(-1, 'P', 'CNT', -1, x, y, z)

                pdb_text += new_pdb_line(-2, 'P', 'PDB', -2, 0, 0, z)
                pdb_text += new_pdb_line(-2, 'P', 'PDB', -2, box[0] * 10, 0, z)
                pdb_text += new_pdb_line(-2, 'P', 'PDB', -2, 0, box[1] * 10, z)
                pdb_text += new_pdb_line(-2, 'P', 'PDB', -2,
                                         box[0] * 10, box[1] * 10, z)

                cutoff = box[0] * 10 * config.params['slice']
                pdb_text += new_pdb_line(-3, 'P', 'PDB', -3, cutoff, cutoff, z)
                pdb_text += new_pdb_line(-3, 'P', 'PDB', -3,
                                         box[0] * 10 - cutoff, cutoff, z)
                pdb_text += new_pdb_line(-3, 'P', 'PDB', -3,
                                         cutoff, box[1] * 10 - cutoff, z)
                pdb_text += new_pdb_line(-3, 'P', 'PDB', -3, box[0] * 10 - cutoff,
                                         box[1] * 10 - cutoff, z)
                pdb_text += new_pdb_line(-4, 'P', 'PDB', -4, -cutoff, -cutoff, z)
                pdb_text += new_pdb_line(-4, 'P', 'PDB', -4,
                                         box[0] * 10 + cutoff, - cutoff, z)
                pdb_text += new_pdb_line(-4, 'P', 'PDB', -4,
                                         -cutoff, box[1] * 10 + cutoff, z)
                pdb_text += new_pdb_line(-4, 'P', 'PDB', -4,
                                         box[0] * 10 + cutoff,
                                         box[1] * 10 + cutoff, z)

                f_new.write(pdb_text)
            with open('P_{1}-{0}.crg'.format(self._name, self._site._res_number), 'w') as f_new:
                f_new.write(crg)

        return acent

    def setDetailsFromTautomer(self):
        """Set DelPhi parameters to run a calculation
        of a single site tautomer"""
        (molecule, delphimol, p_atpos,
         p_rad3, p_chrgv4, atinf, p_iatmed) = self.getMoleculeDetails()
        if config.params['pbc_dim'] == 2:
            (box, half_box_xy, site_center,
             offset_x, offset_y, offset_z,
             box_x, box_y) = self.getCenteringDetails()

        pdb_text = ''
        site_atom_position = -1
        for atom_name, atom_id, atom_position in molecule.iterAtoms():
            if atom_id in self._site.getAtomNumbersList():
                site_atom_position += 1
                p_atpos[site_atom_position]  = molecule.p_atpos[atom_position]
                p_rad3[site_atom_position]   = molecule.p_rad3[atom_position]
                p_chrgv4[site_atom_position] = self.getCharge(atom_name)
                atinf[site_atom_position].value = molecule.atinf[atom_position].value
                # quick fix, should be done only once per site
                # TODO: fix ^
                if config.params['pbc_dim'] == 2:
                    p_atpos[site_atom_position][0] += offset_x
                    p_atpos[site_atom_position][1] += offset_y
                    p_atpos[site_atom_position][2] += offset_z

                    p_atpos[atom_position] = self.putInsideBox(p_atpos[atom_position],
                                                               box_x, box_y)

                p_chrgv4[site_atom_position] = round(p_chrgv4[site_atom_position], 3)
                p_rad3[site_atom_position] = round(p_rad3[site_atom_position], 3)

                aID = int(site_atom_position)
                aname   = str(atinf[site_atom_position].value.split()[0])
                resname = str(atinf[site_atom_position].value.split()[1])
                resnumb = int(atinf[site_atom_position].value.split()[2])
                x = round(p_atpos[site_atom_position][0], 3)
                y = round(p_atpos[site_atom_position][1], 3)
                z = round(p_atpos[site_atom_position][2], 3)
                pdb_text += new_pdb_line(aID, aname, resname, resnumb, x, y, z)

        delphimol.changeStructureSize(self._natoms, p_atpos, p_rad3,
                                      p_chrgv4, atinf, p_iatmed)
        acent = self.getSiteAcent()

        if config.debug:
            filename = '{1}-{0}.pdb'.format(self._name, self._site._res_number)
            with open(filename, 'w') as f_new:
                x, y, z = acent
                pdb_text += new_pdb_line(-1, 'P', 'CNT', -1, x, y, z)
                f_new.write(pdb_text)

        return acent

    def add_pbc(self, x, y, z, box, radius, charge, inf):
        def pdb_y(y, cutoff_y, box, new_atoms,
                  x_new, z_new, radius, charge, inf):
            if (y < cutoff_y):
                y_new = box + y
                new_atoms.append((x_new, y_new, z_new, radius, charge, inf))
            elif (y > box - cutoff_y):
                y_new = y - box
                new_atoms.append((x_new, y_new, z_new, radius, charge, inf))
            return new_atoms
            
        box = box * 10
        x_new = x
        y_new = y
        z_new = z
        inf = inf[1:]
        charge = 0.0
        cutoff_x = config.params['slice'] * box
        cutoff_y = config.params['slice'] * box
        new_atoms = []
        if (x < cutoff_x):
            x_new = box + x
            new_atoms.append((x_new, y_new, z_new, radius, charge, inf))
            new_atoms = pdb_y(y, cutoff_y, box, new_atoms,
                              x_new, z_new, radius, charge, inf)
        elif (x > box - cutoff_y):
            x_new = x - box
            new_atoms.append((x_new, y_new, z_new, radius, charge, inf))
            new_atoms = pdb_y(y, cutoff_y, box, new_atoms,
                              x_new, z_new, radius, charge, inf)
        x_new = x
        new_atoms = pdb_y(y, cutoff_y, box, new_atoms,
                          x_new, z_new, radius, charge, inf)

        return new_atoms

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
        for i in list(self._charge_set.keys()):
            out += '{0:>7.3f} {1}\n'.format(self._charge_set[i], i)
        return out

    # Assertion Methods
    def isRefTautomer(self):
        if self == self._site._ref_tautomer:
            return True
        else:
            return False

    # Calculation Methods
    def CalcPotentialTautomer(self):
        """Run DelPhi simulation of single site tautomer

        Ensures:
            self._esolvation (float): tautomer solvation energy
            self._p_sitpot (list): potential on site atoms
        """
        if config.debug:
            t0 = time.clock()
            print((self._name))
            print((self._charge_set))
        molecule = self._molecule
        delphimol = molecule.getDelPhi()
        acent = self.setDetailsFromTautomer()


        if config.params['bndcon'] == 3:
            ibctyp = 4
        else:
            ibctyp = config.params['bndcon']

        filename = '{0}_{1}.prm'.format(self._name, self._site._res_number)
        logfile = 'LOG_runDelPhi_{0}_{1}_modelcompound'.format(self._name,
                                                               self._site._res_number)
        if config.debug:
            print(('started', self._name, self._site._res_number, 'modelcompound'))

        delphimol.runDelPhi(scale=config.params['scaleM'],
                            nonit=0, nlit=500, relpar=0, relfac=0,
                            acent=acent, pbx=False, pby=False, debug=config.debug,
                            filename=filename, ibctyp=ibctyp,
                            outputfile=logfile)
        if config.debug:
            print(('ended', self._name, self._site._res_number, 'modelcompound'))
        
        log.checkDelPhiErrors(logfile, 'runDelPhi')
        
        self._esolvation = delphimol.getSolvation()
        self._p_sitpot   = delphimol.getSitePotential()
        if config.debug:
            t1 = time.clock() - t0
            filename = '{0}_{1}.profl'.format(self._name, self.getSiteResNumber())
            with open(filename, 'a') as f_new:
                f_new.write('time -> {0:10} {1:10}\n'.format(t0, t1))
        
        return self._esolvation, self._p_sitpot[:]

    def CalcPotentialTitratingMolecule(self):
        """Run DelPhi simulation of the site tautomer
        within the whole molecule

        Ensures:
            self._esolvation (float): tautomer solvation energy
            self._p_sitpot (list): potential on site atoms
        """
        
        if config.debug:
            start = time.clock()
        molecule = self._molecule
        delphimol = molecule.getDelPhi()

        acent = self.setDetailsFromWholeMolecule()

        if config.debug:
            t0 = time.clock() - start

            print((self._name, 'starting'))
            p_atpos  = delphimol.get_atpos()
            p_rad3   = delphimol.get_rad3()
            p_chrgv4 = delphimol.get_chrgv4()
            atinf    = delphimol.get_atinf()             
            #for atom_name, atom_id, atom_position in molecule.iterAtoms():
            #    print (atinf[atom_position].value, p_chrgv4[atom_position], p_rad3[atom_position],
            #           p_atpos[atom_position][0], p_atpos[atom_position][1], p_atpos[atom_position][2])

        filename = '{0}_{1}.prm'.format(self._name, self._site._res_number)
        logfile = 'LOG_runDelPhi_{0}_{1}_wholeprotein'.format(self._name,
                                                              self._site._res_number)

        if config.debug:
            print(('started', self._name, self._site._res_number, 'wholeprotein'))
        if config.params['pbc_dim'] == 2:
            delphimol.runDelPhi(scale_prefocus=config.params['scaleP'],
                                scale=config.params['scaleM'],
                                nlit_prefocus=config.params['nlit'],
                                nonit=config.params['nonit'],
                                nlit=500, acent=acent, nonit_focus=0,
                                relfac_focus=0.0, relpar_focus=0.0,
                                relpar=config.params['relpar'],
                                relfac=config.params['relfac'],
                                pbx=True,
                                pby=True, pbx_focus=False,
                                pby_focus=False, ibctyp=4,
                                debug=config.debug,
                                filename=filename,
                                outputfile=logfile)
        elif config.params['pbc_dim'] == 0 and config.params['bndcon'] == 3:
            delphimol.runDelPhi(scale_prefocus=config.params['scaleP'],
                                scale=config.params['scaleM'],
                                nlit_prefocus=config.params['nlit'],
                                nonit=0,
                                nlit=500, acent=acent, nonit_focus=0,
                                relfac_focus=0.0, relpar_focus=0.0,
                                relpar=config.params['relpar'],
                                relfac=config.params['relfac'],
                                pbx=False,
                                pby=False, pbx_focus=False,
                                pby_focus=False, ibctyp=4,
                                debug=config.debug,
                                filename=filename,
                                outputfile=logfile)
        else:
            delphimol.runDelPhi(scale=config.params['scaleM'],
                                nonit=0, nlit=500, relpar=0, relfac=0,
                                acent=acent, pbx=False, pby=False,
                                debug=config.debug,
                                filename=filename,
                                outputfile=logfile)
            
        if config.debug:
            print(('ended', self._name, self._site._res_number, 'wholeprotein'))

        log.checkDelPhiErrors(logfile, 'runDelPhi')

        self._esolvation = delphimol.getSolvation()
        self._p_sitpot   = delphimol.getSitePotential()

        if config.debug:
            filename = '{0}_{1}.frc'.format(self._name,
                                            self._site._res_number)
            with open(filename, 'w') as f:
                text = ''
                for atom_name, atom_id, atom_position in molecule.iterAtoms():
                    text += '{0} {1} {2} '\
                            '{3} {4} {5} {6}\n'.format(atinf[atom_position].value,
                                                       round(p_chrgv4[atom_position], 3),
                                                       round(p_rad3[atom_position], 4),
                                                       round(p_atpos[atom_position][0], 3),
                                                       round(p_atpos[atom_position][1], 3),
                                                       round(p_atpos[atom_position][2], 3),
                                                       self._p_sitpot[atom_position])
                f.write(text)

            t1 = time.clock() - start
            filename = '{0}_{1}.profl'.format(self._name,
                                              self.getSiteResNumber())
            with open(filename, 'a') as f_new:
                f_new.write('time -> {0:10}     {1:10}\n'.format(t0, t1))

            print((self._esolvation, self._name))

        return self._esolvation, self._p_sitpot[:]

    def calcBackEnergy(self):
        """Calculates background energy contribution"""
        if config.debug:
            print((self._name, 'background energy start'))
        molecule = self._molecule
        text = ''
        distance = -999999
        cutoff = copy(config.params['cutoff'])
        cutoff2 = (cutoff * 10) ** 2
        point_energy = -1
        for atom_name, atom_id, atom_position in molecule.iterAtoms():
            if atom_id not in self._site.getAtomNumbersList():
                if cutoff != -1:
                    distance = self.distance_to(molecule.p_atpos[atom_position])
                if cutoff != -1 or distance <= cutoff2:
                    point_energy = round(molecule.p_chrgv4[atom_position], 3) * round(self._sitpotM[atom_position], 4)
                    self._e_back += point_energy
                    if config.debug:
                        point_energy = round(molecule.p_chrgv4[atom_position], 3) * round(self._sitpotM[atom_position], 4)
                        text += '{} {} {} '\
                                '{} {} {} {} {}\n'.format(atom_name,
                                                          atom_id,
                                                          point_energy,
                                                          molecule.p_chrgv4[atom_position],
                                                          self._sitpotM[atom_position],
                                                          molecule.p_atpos[atom_position][:],
                                                          distance,
                                                          cutoff2)

        if config.debug:
            filename = '{0}_{1}_eback.xvg'.format(self._name,
                                                  self._site._res_number)
            with open(filename, 'w') as f_new:
                text += str(self._e_back / config.log10)
                f_new.write(text)

            print('e_back finished')

    def calcpKint(self):
        """Calculates the pKint of the tautomer"""
        ref_taut = self._site._ref_tautomer

        dG_solvationM = ref_taut._esolvationM - self._esolvationM
        dG_solvationS = ref_taut._esolvationS - self._esolvationS
        dG_back       = ref_taut._e_back - self._e_back

        dG_solvationM /= config.log10
        dG_solvationS /= config.log10
        dG_back       /= config.log10

        if self._site.getType() == 'a':
            pKint = self._pKmod + (dG_solvationM - dG_solvationS + dG_back)
            chargediff = -1
        elif self._site.getType() == 'c':
            pKint = self._pKmod - (dG_solvationM - dG_solvationS + dG_back)
            chargediff = 1
        else:
            raise Exception('Site files were poorly interpreted')

        dg = pKint * config.log10 * config.kBoltz * float(config.params['temp']) * chargediff
        if config.debug:
            print(('pKint ->', self._name, pKint, dg))
        self.pKint = pKint
        self.dG_solvationM = dG_solvationM
        self.dG_solvationS = dG_solvationS
        self.dG_back = dG_back
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
                charge_ref = molecule.p_chrgv4[atom_position]
                charge_tau = self.getCharge(atom_name)
                charge = round(charge_ref, 3) - round(charge_tau, 3)

                potential_ref = tau2_ref._sitpotM[atom_position]
                potential_tau2 = tautomer2._sitpotM[atom_position]
                potential = round(potential_ref, 4) - round(potential_tau2, 4)

                interaction += charge * potential

        site1_chrgtyp = self._site.getRefProtState()
        site2_chrgtyp = tautomer2._site.getRefProtState()

        dG_interaction = site1_chrgtyp * site2_chrgtyp
        dG_interaction *= abs(interaction * config.kBoltz * float(config.params['temp']))

        return dG_interaction

    def distance_to(self, atom_position):
        center = self._site.getCenterH()

        dx = abs(center[0] - atom_position[0])
        dy = abs(center[1] - atom_position[1])
        dz = abs(center[2] - atom_position[2])

        box = [i * 10 for i in self._molecule.box]

        if dx > box[0] / 2:
            dx = abs(dx - box[0])
        if dy > box[1] / 2:
            dy = abs(dy - box[1])

        return dx ** 2 + dy ** 2 + dz ** 2
