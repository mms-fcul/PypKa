import time
from copy import copy

from pypka.log import checkDelPhiErrors
from pypka.config import Config
from pypka.constants import KBOLTZ, LOG10

from pdbmender.formats import new_pdb_line


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
        self.molecule
        self.site

        # Tautomer Details
        self.name
        self.charge_set (dict): charge set of the tautomer
                                 key is atom name str
                                 value is charge float
        self.natoms (int): number of atoms of the Tautomer
        Redundant to self.site.natoms as it must be the same
        between all Tautomers belonging to the same Site.

        # Tautomer Energy Details
        self.esolvation (float): solvation energy
        self.e_back (float): background interaction energy
        self.dg (float): pKint energy
        """
        self.molecule = molecule
        self.site = site
        self.name = name
        self.esolvation = ""
        self.charge_set = {}
        self.natoms = 0
        self.e_back = 0.0
        self.dg = 0.0

    def __lt__(self, other):
        return self.name < other.name

    # Set Methods
    def loadChargeSet(self, res_name, ref_tautomer):
        """Reads .st file related to Tautomer with residue name res_name

        Stores charge set related to both the Tautomer and the
        Reference Tautomer

        .st file named TYRtau1.st has charge set of TY0 and reference TY2
                       TYRtau2.st has charge set of TY1 and reference TY2
        """
        tau_number = int(self.name[-1]) + 1
        fname = "{0}/{1}tau{2}.st".format(
            Config.pypka_params["sts_path"],
            res_name,
            tau_number,
        )
        with open(fname) as f:
            nline = 0
            charge_set1 = {}
            charge_set2 = {}
            for line in f:
                nline += 1
                cols = line.split()
                if nline > 1:
                    atom_name = cols[1]
                    # if atom_name == "C" and res_name == "CTR":
                    #    atom_name = "CT"
                    # elif atom_name in ('O1', 'O2') and res_name == 'CTR':
                    #    atom_name = atom_name[0] + 'T' + atom_name[1]
                    charge_set1[atom_name] = float(cols[-2])
                    charge_set2[atom_name] = float(cols[-1])
                else:
                    self.pKmod = float(line)

        if sum(charge_set1.values()) < 0.001 and sum(charge_set2.values()) > -1.001:
            if Config.debug:
                print((fname, "anionic"))
            self.charge_set = charge_set1
            ref_tautomer.charge_set = charge_set2
            self.site.setType("a")
        elif sum(charge_set1.values()) > 0.99 and sum(charge_set2.values()) < 0.001:
            if Config.debug:
                print((fname, "cationic"))
            self.charge_set = charge_set2
            ref_tautomer.charge_set = charge_set1
            self.site.setType("c")
        self.natoms = len(self.charge_set)
        ref_tautomer.natoms = len(self.charge_set)
        if Config.debug:
            print((self.charge_set))
            print((ref_tautomer.charge_set))
            print(("finished reading", fname))

    def saveDelPhiResults(self, esolvationS, sitpotS, esolvationM, sitpotM):
        self.esolvationS = esolvationS
        self.sitpotS = sitpotS
        self.esolvationM = esolvationM
        self.sitpotM = sitpotM

    def getSiteAcent(self):
        if Config.delphi_params["pbc_dim"] == 2:
            x = Config.pypka_params["box"][0] / 2
            y = x
            z = self.site.center[2]
            acent = [x, y, z]
        else:
            acent = self.site.center

        return acent

    def getMoleculeDetails(self):
        molecule = self.molecule
        delphimol = copy(Config.delphi_params["delphimol"])

        # Could have used empty structures
        # TODO: quantify time impact of this quick fix
        p_atpos = copy(delphimol.p_atpos)
        p_rad3 = copy(delphimol.p_rad3)
        p_chrgv4 = copy(delphimol.p_chrgv4)
        atinf = copy(delphimol.atinf)
        p_iatmed = copy(delphimol.p_iatmed)

        return molecule, delphimol, p_atpos, p_rad3, p_chrgv4, atinf, p_iatmed

    def getCenteringDetails(self):
        box = Config.pypka_params["box"]
        half_box_xy = self.site.getCenter()[0]
        site_center = self.site.getCenterOriginal()

        offset_x = half_box_xy - site_center[0]
        offset_y = half_box_xy - site_center[1]
        offset_z = box[2]  # irrelevant, only for debug

        box_x = box[0]
        box_y = box[1]

        return (
            box,
            half_box_xy,
            site_center,
            offset_x,
            offset_y,
            offset_z,
            box_x,
            box_y,
        )

    @staticmethod
    def putInsideBox(point, box_x, box_y):
        if point[0] < 0:
            point[0] += box_x
        elif point[0] >= box_x:
            point[0] -= box_x
        if point[1] < 0:
            point[1] += box_y
        elif point[1] >= box_y:
            point[1] -= box_y
        return point

    @staticmethod
    def read_atom_details(atinf, atom_position):
        atinf_str = atinf[atom_position].value.decode("ascii")
        aname = atinf_str[:5].strip()
        resname = atinf_str[6:10].strip()
        chain = atinf_str[10:11].strip()
        resnumb = int(atinf_str[11:].strip())
        return aname, resname, chain, resnumb

    @staticmethod
    def read_atom_positions(p_atpos, atom_position):
        x = float(p_atpos[atom_position][0])
        y = float(p_atpos[atom_position][1])
        z = float(p_atpos[atom_position][2])
        return x, y, z

    @staticmethod
    def regular_delphi_run(delphimol, acent, filename, logfile, ibctyp=None):
        delphimol.runDelPhi(
            scale=Config.delphi_params["scaleM"],
            nonit=0,
            nlit=500,
            relpar=0,
            relfac=0,
            acent=acent,
            pbx=False,
            pby=False,
            ibctyp=ibctyp,
            debug=Config.debug,
            filename=filename,
            outputfile=logfile,
        )

    def setDetailsFromWholeMolecule(self):
        """Set DelPhi parameters to run a calculation of a whole molecule (all sites neutral, except one)."""
        (
            _,
            delphimol,
            p_atpos,
            p_rad3,
            p_chrgv4,
            atinf,
            p_iatmed,
        ) = self.getMoleculeDetails()
        if Config.delphi_params["pbc_dim"] == 2:
            (
                box,
                _,
                _,
                offset_x,
                offset_y,
                offset_z,
                box_x,
                box_y,
            ) = self.getCenteringDetails()
        pdb_text = ""
        crg = "!crg file created by gen_files.awk\n" "atom__resnumbc_charge_\n"
        siz = "!siz file created by gen_files.awk\n" "atom__res_radius_\n"
        new_atoms = []

        lookup_atoms_keys = Config.delphi_params.lookup_atoms.keys()

        for atom_position in range(delphimol.natoms):
            aID = atom_position + 1

            aname, resname, chain, resnumb = self.read_atom_details(
                atinf, atom_position
            )
            x, y, z = self.read_atom_positions(p_atpos, atom_position)

            pdb_text += new_pdb_line(aID, aname, resname, resnumb, x, y, z, chain=chain)

            if atom_position in lookup_atoms_keys:
                atom_id, atom_name = Config.delphi_params.lookup_atoms[atom_position]

                if atom_id in self.site.getAtomNumbersList():
                    p_chrgv4[atom_position] = self.getCharge(atom_name)
                else:
                    p_chrgv4[atom_position] = 0.0
            else:
                p_chrgv4[atom_position] = 0.0

            if float(p_chrgv4[atom_position]) < 0.0:
                signal = "-"
            else:
                signal = " "

            crg += "{0:<6}{1:<7}  {3}{2:0<5} \n".format(
                aname, resname, round(abs(p_chrgv4[atom_position]), 3), signal
            )

            siz += "{0:<6}{1:<4} {2:0<5} \n".format(
                aname, resname, round(abs(p_rad3[atom_position]), 3)
            )

            # TODO: quick fix, should be done only once per site
            if Config.delphi_params["pbc_dim"] == 2:
                p_atpos[atom_position][0] += offset_x
                p_atpos[atom_position][1] += offset_y
                p_atpos[atom_position][2] += offset_z

                p_atpos[atom_position] = self.putInsideBox(
                    p_atpos[atom_position], box_x, box_y
                )

                # Comment 15 May 2019
                # No rounding is better, however, pypka calculations
                # are no longer comparable to delphiT ones due to this #noregrets
                # rounding_precision = decimal.Decimal('0.01')
                # x_dec = decimal.Decimal(p_atpos[atom_position][0])
                # y_dec = decimal.Decimal(p_atpos[atom_position][1])
                # z_dec = decimal.Decimal(p_atpos[atom_position][2])
                # p_atpos[atom_position][0] = float(x_dec.quantize(rounding_precision))
                # p_atpos[atom_position][1] = float(y_dec.quantize(rounding_precision))
                # p_atpos[atom_position][2] = float(z_dec.quantize(rounding_precision))

                pbc_atoms = self.add_pbc(
                    p_atpos[atom_position][0],
                    p_atpos[atom_position][1],
                    p_atpos[atom_position][2],
                    box[0],
                    p_rad3[atom_position],
                    p_chrgv4[atom_position],
                    atinf[atom_position].value,
                )
                # pbc_atoms = []
                new_atoms += pbc_atoms

        natoms = delphimol.changeStructureSize(
            p_atpos,
            p_rad3,
            p_chrgv4,
            atinf,
            p_iatmed,
            extra_atoms=new_atoms,
            natoms=delphimol.natoms,
        )
        acent = self.getSiteAcent()

        pdb_text = ""
        for i in range(natoms):
            aID = i + 1

            aname, resname, chain, resnumb = self.read_atom_details(delphimol.atinf, i)
            x, y, z = self.read_atom_positions(delphimol.p_atpos, i)

            if "-" in aname[0]:
                aname = aname[1:]
                resname = "PBC"
                resnumb = -666
            pdb_text += new_pdb_line(aID, aname, resname, resnumb, x, y, z, chain=chain)

        box = Config.pypka_params.box
        if Config.debug:
            filename = "P_{1}-{0}.pdb".format(self.name, self.site.res_number)
            with open(filename, "w") as f_new:
                x, y, z = acent
                pdb_text += new_pdb_line(-1, "P", "CNT", -1, x, y, z)

                pdb_text += new_pdb_line(-2, "P", "PDB", -2, 0, 0, z)
                pdb_text += new_pdb_line(-2, "P", "PDB", -2, box[0], 0, z)
                pdb_text += new_pdb_line(-2, "P", "PDB", -2, 0, box[1], z)
                pdb_text += new_pdb_line(-2, "P", "PDB", -2, box[0], box[1], z)

                cutoff = box[0] * Config.pypka_params["slice"]
                pdb_text += new_pdb_line(-3, "P", "PDB", -3, cutoff, cutoff, z)
                pdb_text += new_pdb_line(-3, "P", "PDB", -3, box[0] - cutoff, cutoff, z)
                pdb_text += new_pdb_line(-3, "P", "PDB", -3, cutoff, box[1] - cutoff, z)
                pdb_text += new_pdb_line(
                    -3, "P", "PDB", -3, box[0] - cutoff, box[1] - cutoff, z
                )
                pdb_text += new_pdb_line(-4, "P", "PDB", -4, -cutoff, -cutoff, z)
                pdb_text += new_pdb_line(
                    -4, "P", "PDB", -4, box[0] + cutoff, -cutoff, z
                )
                pdb_text += new_pdb_line(
                    -4, "P", "PDB", -4, -cutoff, box[1] + cutoff, z
                )
                pdb_text += new_pdb_line(
                    -4, "P", "PDB", -4, box[0] + cutoff, box[1] + cutoff, z
                )

                f_new.write(pdb_text)
            with open(
                "P_{1}-{0}.crg".format(self.name, self.site.res_number), "w"
            ) as f_new:
                f_new.write(crg)
            with open(
                "P_{1}-{0}.siz".format(self.name, self.site.res_number), "w"
            ) as f_new:
                f_new.write(siz)

        return delphimol, acent

    def setDetailsFromTautomer(self):
        """Set DelPhi parameters to run a calculation of a single site tautomer."""
        (
            molecule,
            delphimol,
            p_atpos,
            p_rad3,
            p_chrgv4,
            atinf,
            p_iatmed,
        ) = self.getMoleculeDetails()

        if Config.delphi_params["pbc_dim"] == 2:
            (
                _,
                half_box_xy,
                site_center,
                offset_x,
                offset_y,
                offset_z,
                box_x,
                box_y,
            ) = self.getCenteringDetails()

        pdb_text = ""
        site_atom_position = -1
        for atom_name, atom_id, atom_position in molecule.iterAtoms():
            if atom_id in self.site.getAtomNumbersList():
                site_atom_position += 1
                p_atpos[site_atom_position] = delphimol.p_atpos[atom_position]
                p_rad3[site_atom_position] = delphimol.p_rad3[atom_position]
                p_chrgv4[site_atom_position] = self.getCharge(atom_name)
                atinf[site_atom_position].value = delphimol.atinf[atom_position].value
                # TODO: should be done only once per site
                if Config.delphi_params["pbc_dim"] == 2:
                    p_atpos[site_atom_position][0] += offset_x
                    p_atpos[site_atom_position][1] += offset_y
                    p_atpos[site_atom_position][2] += offset_z

                    p_atpos[atom_position] = self.putInsideBox(
                        p_atpos[atom_position], box_x, box_y
                    )

                # p_chrgv4[site_atom_position] = round(p_chrgv4[site_atom_position], 3)
                # p_rad3[site_atom_position] = round(p_rad3[site_atom_position], 3)

                aID = int(site_atom_position)

                aname, resname, chain, resnumb = self.read_atom_details(
                    atinf, site_atom_position
                )
                x, y, z = self.read_atom_positions(p_atpos, site_atom_position)

                pdb_text += new_pdb_line(
                    aID, aname, resname, resnumb, x, y, z, chain=chain
                )

        delphimol.changeStructureSize(
            p_atpos, p_rad3, p_chrgv4, atinf, p_iatmed, natoms=self.natoms
        )
        acent = self.getSiteAcent()

        if Config.debug:
            filename = "{1}-{0}.pdb".format(self.name, self.site.res_number)
            with open(filename, "w") as f_new:
                x, y, z = acent
                pdb_text += new_pdb_line(-1, "P", "CNT", -1, x, y, z)
                f_new.write(pdb_text)

        return delphimol, acent

    @staticmethod
    def add_pbc(x, y, z, box, radius, charge, inf):
        def pdb_y(y, cutoff_y, box, new_atoms, x_new, z_new, radius, charge, inf):
            if y < cutoff_y:
                y_new = box + y
                new_atoms.append((x_new, y_new, z_new, radius, charge, inf))
            elif y > box - cutoff_y:
                y_new = y - box
                new_atoms.append((x_new, y_new, z_new, radius, charge, inf))
            return new_atoms

        x_new = x
        y_new = y
        z_new = z
        inf = inf[1:]
        charge = 0.0
        cutoff_x = Config.pypka_params["slice"] * box
        cutoff_y = Config.pypka_params["slice"] * box
        new_atoms = []
        if x < cutoff_x:
            x_new = box + x
            new_atoms.append((x_new, y_new, z_new, radius, charge, inf))
            new_atoms = pdb_y(
                y, cutoff_y, box, new_atoms, x_new, z_new, radius, charge, inf
            )
        elif x > box - cutoff_y:
            x_new = x - box
            new_atoms.append((x_new, y_new, z_new, radius, charge, inf))
            new_atoms = pdb_y(
                y, cutoff_y, box, new_atoms, x_new, z_new, radius, charge, inf
            )
        x_new = x
        new_atoms = pdb_y(
            y, cutoff_y, box, new_atoms, x_new, z_new, radius, charge, inf
        )

        return new_atoms

    # Get Methods
    def getSiteResNumber(self):
        return self.site.getResNumber()

    def getCharge(self, atom_name):
        return self.charge_set[atom_name]

    def getName(self):
        return self.name

    # Print Methods
    def __str__(self):
        out = self.name + "\n"
        for i in list(self.charge_set.keys()):
            out += "{0:>7.3f} {1}\n".format(self.charge_set[i], i)
        return out

    # Assertion Methods
    def isRefTautomer(self):
        if self == self.site.ref_tautomer:
            return True
        else:
            return False

    # Calculation Methods
    def CalcPotentialTautomer(self):
        """
        Run DelPhi simulation of single site tautomer

        Ensures:
            self.esolvation (float): tautomer solvation energy
            self.p_sitpot (list): potential on site atoms
        """
        if Config.debug:
            t0 = time.process_time()
            print((self.name))
            print((self.charge_set))
        _ = self.molecule

        delphimol, acent = self.setDetailsFromTautomer()

        if Config.delphi_params["bndcon"] == 3:
            ibctyp = 4
        else:
            ibctyp = Config.delphi_params["bndcon"]

        filename = "{0}_{1}.prm".format(self.name, self.site.res_number)
        logfile = "LOG_runDelPhi_{0}_{1}_{2}_modelcompound".format(
            self.name, self.site.res_number, self.molecule.chain
        )
        if Config.debug:
            print(("started", self.name, self.site.res_number, "modelcompound"))

        self.regular_delphi_run(delphimol, acent, filename, logfile, ibctyp=ibctyp)

        if Config.debug:
            print(("ended", self.name, self.site.res_number, "modelcompound"))

        checkDelPhiErrors(logfile, "runDelPhi")

        self.esolvation = delphimol.getSolvation()
        self.p_sitpot = delphimol.getSitePotential()
        if Config.debug:
            t1 = time.process_time() - t0
            filename = "{0}_{1}.profl".format(self.name, self.getSiteResNumber())
            with open(filename, "a") as f_new:
                f_new.write("time -> {0:10} {1:10}\n".format(t0, t1))

        return self.esolvation, self.p_sitpot[:]

    def CalcPotentialTitratingMolecule(self):
        """
        Run DelPhi simulation of the site tautomer
        within the whole molecule

        Ensures:
            self.esolvation (float): tautomer solvation energy
            self.p_sitpot (list): potential on site atoms
        """

        if Config.debug:
            start = time.process_time()
        molecule = self.molecule

        delphimol, acent = self.setDetailsFromWholeMolecule()

        if Config.debug:
            t0 = time.process_time() - start

            print((self.name, "starting"))
            p_atpos = delphimol.get_atpos()
            p_rad3 = delphimol.get_rad3()
            p_chrgv4 = delphimol.get_chrgv4()
            atinf = delphimol.get_atinf()
            # for atom_name, atom_id, atom_position in molecule.iterAtoms():
            #    print (atinf[atom_position].value, p_chrgv4[atom_position], p_rad3[atom_position],
            #           p_atpos[atom_position][0], p_atpos[atom_position][1], p_atpos[atom_position][2])

        filename = "{0}_{1}.prm".format(self.name, self.site.res_number)
        logfile = "LOG_runDelPhi_{0}_{1}_{2}_wholeprotein".format(
            self.name, self.site.res_number, self.molecule.chain
        )

        if Config.debug:
            print(("started", self.name, self.site.res_number, "wholeprotein"))
        if Config.delphi_params["pbc_dim"] == 2:
            delphimol.runDelPhi(
                scale_prefocus=Config.delphi_params["scaleP"],
                scale=Config.delphi_params["scaleM"],
                nlit_prefocus=Config.delphi_params["nlit"],
                nonit=Config.delphi_params["nonit"],
                nlit=500,
                acent=acent,
                nonit_focus=0,
                relfac_focus=0.0,
                relpar_focus=0.0,
                relpar=Config.delphi_params["relpar"],
                relfac=Config.delphi_params["relfac"],
                pbx=True,
                pby=True,
                pbx_focus=False,
                pby_focus=False,
                ibctyp=4,
                debug=Config.debug,
                filename=filename,
                outputfile=logfile,
            )
        elif (
            Config.delphi_params["pbc_dim"] == 0 and Config.delphi_params["bndcon"] == 3
        ):
            delphimol.runDelPhi(
                scale_prefocus=Config.delphi_params["scaleP"],
                scale=Config.delphi_params["scaleM"],
                nlit_prefocus=Config.delphi_params["nlit"],
                nonit=Config.delphi_params["nonit"],
                nlit=500,
                acent=acent,
                nonit_focus=0,
                relfac_focus=0.0,
                relpar_focus=0.0,
                relpar=Config.delphi_params["relpar"],
                relfac=Config.delphi_params["relfac"],
                pbx=False,
                pby=False,
                pbx_focus=False,
                pby_focus=False,
                ibctyp=4,
                debug=Config.debug,
                filename=filename,
                outputfile=logfile,
            )
        else:
            self.regular_delphi_run(delphimol, acent, filename, logfile)

        if Config.debug:
            print(("ended", self.name, self.site.res_number, "wholeprotein"))

        checkDelPhiErrors(logfile, "runDelPhi")

        self.esolvation = delphimol.getSolvation()
        self.p_sitpot = delphimol.getSitePotential()

        if Config.debug:
            filename = "{0}_{1}.frc".format(self.name, self.site.res_number)
            with open(filename, "w") as f:
                text = ""
                for _, atom_id, atom_position in molecule.iterAtoms():
                    text += "{0} {1} {2} " "{3} {4} {5} {6}\n".format(
                        atinf[atom_position].value,
                        round(p_chrgv4[atom_position], 3),
                        round(p_rad3[atom_position], 4),
                        round(p_atpos[atom_position][0], 3),
                        round(p_atpos[atom_position][1], 3),
                        round(p_atpos[atom_position][2], 3),
                        self.p_sitpot[atom_position],
                    )
                f.write(text)

            t1 = time.process_time() - start
            filename = "{0}_{1}.profl".format(self.name, self.getSiteResNumber())
            with open(filename, "a") as f_new:
                f_new.write("time -> {0:10}     {1:10}\n".format(t0, t1))

            print((self.esolvation, self.name))

        return self.esolvation, self.p_sitpot[:]

    def calcBackEnergy(self):
        """Calculates background energy contribution."""
        if Config.debug:
            print((self.name, "background energy start"))
        _ = self.molecule
        delphimol = Config.delphi_params["delphimol"]
        text = ""
        distance = -999999
        cutoff = Config.pypka_params["cutoff"] * 10
        cutoff_sq = (cutoff) ** 2
        point_energy = -1

        lookup_atoms_keys = Config.delphi_params.lookup_atoms.keys()
        # for atom_name, atom_id, atom_position in molecule.iterAtoms():
        for atom_position in range(delphimol.natoms):
            atom_id = None
            if atom_position in lookup_atoms_keys:
                atom_id, atom_name = Config.delphi_params.lookup_atoms[atom_position]

            if atom_id == None or atom_id not in self.site.getAtomNumbersList():
                if cutoff != -10:  # -1 is off times 10 from the conversion
                    distance = self.distance_to(delphimol.p_atpos[atom_position])
                    if distance > cutoff_sq:
                        continue

                point_energy = (
                    delphimol.p_chrgv4[atom_position] * self.sitpotM[atom_position]
                )
                self.e_back += point_energy
                if Config.debug:
                    text += "{} {} {} {} " "{} {} {} {} {}\n".format(
                        atom_position,
                        atom_name,
                        atom_id,
                        point_energy,
                        delphimol.p_chrgv4[atom_position],
                        self.sitpotM[atom_position],
                        delphimol.p_atpos[atom_position][:],
                        distance,
                        cutoff_sq,
                    )

        if Config.debug:
            filename = "{0}_{1}_eback.xvg".format(self.name, self.site.res_number)
            with open(filename, "w") as f_new:
                text += str(self.e_back / LOG10)
                f_new.write(text)

            print("e_back finished")

    def calcpKint(self):
        """Calculates the pKint of the tautomer."""
        ref_taut = self.site.ref_tautomer

        dG_solvationM = ref_taut.esolvationM - self.esolvationM
        dG_solvationS = ref_taut.esolvationS - self.esolvationS
        dG_back = ref_taut.e_back - self.e_back

        dG_solvationM /= LOG10
        dG_solvationS /= LOG10
        dG_back /= LOG10

        if self.site.getType() == "a":
            pKint = self.pKmod + (dG_solvationM - dG_solvationS + dG_back)
            chargediff = -1
        elif self.site.getType() == "c":
            pKint = self.pKmod - (dG_solvationM - dG_solvationS + dG_back)
            chargediff = 1
        else:
            raise Exception("Site files were poorly interpreted")

        dg = pKint * LOG10 * KBOLTZ * float(Config.pypka_params["temp"]) * chargediff
        if Config.debug:
            print(("pKint ->", self.name, pKint, dg))
        self.pKint = pKint
        self.dG_solvationM = dG_solvationM
        self.dG_solvationS = dG_solvationS
        self.dG_back = dG_back
        self.dg = dg

    def calcInteractionWith(self, tautomer2, site_atom_list):
        """
        Calculates the interaction energy between self tautomer and tautomer2

        Args:
            tautomer2 (Tautomer)
            site_atom_list (list): atom numbers belonging to the site
        """
        delphimol = Config.delphi_params["delphimol"]
        molecule = self.molecule
        interaction = 0.0
        tau2_ref = tautomer2.site.ref_tautomer
        for atom_name, atom_id, atom_position in molecule.getAtomsList():
            if atom_id in site_atom_list:
                charge_ref = delphimol.p_chrgv4[atom_position]
                charge_tau = self.getCharge(atom_name)
                charge = charge_ref - charge_tau

                potential_ref = tau2_ref.sitpotM[atom_position]
                potential_tau2 = tautomer2.sitpotM[atom_position]
                potential = potential_ref - potential_tau2

                interaction += charge * potential

        site1_chrgtyp = self.site.getRefProtState()
        site2_chrgtyp = tautomer2.site.getRefProtState()

        dG_interaction = site1_chrgtyp * site2_chrgtyp
        dG_interaction *= abs(interaction * KBOLTZ * float(Config.pypka_params["temp"]))

        return dG_interaction

    def distance_to(self, atom_position):
        center = self.site.getCenterH()

        dx = abs(center[0] - atom_position[0])
        dy = abs(center[1] - atom_position[1])
        dz = abs(center[2] - atom_position[2])

        box = Config.pypka_params["box"]

        if dx > box[0] / 2:
            dx = abs(dx - box[0])
        if dy > box[1] / 2:
            dy = abs(dy - box[1])

        return dx ** 2 + dy ** 2 + dz ** 2
