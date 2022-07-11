from copy import copy
from pypka.config import Config
from pypka.constants import (
    TERMINAL_OFFSET,
    PROTEIN_RESIDUES,
    TITRABLERESIDUES,
    TITRABLETAUTOMERS,
    REGULARTITRATINGRES,
)
from pypka.clean.ffconverter import main_chains
from pdbmender.formats import new_pdb_line, read_pdb_line
import logging

logger = logging.getLogger(__name__)


def warn_integrity_check_failed(resnumb, resname):
    warn = "{0} {1} failed integrity check".format(resnumb, resname)
    logger.warning(warn)


def addSiteAndTautomers(molecule, resnumb, resname, termini_resname=None):
    if resname in TITRABLETAUTOMERS:
        ntautomers = TITRABLETAUTOMERS[resname]
    else:
        for res in REGULARTITRATINGRES:
            if res[0:2] == resname[0:2]:
                ntautomers = TITRABLETAUTOMERS[res]
    sID = molecule.addSite(resnumb)
    molecule.addTautomers(sID, ntautomers, resname, termini_resname=termini_resname)


def check_CYS_bridge(res_atoms, resnumb):
    if Config.pypka_params["ff_family"] == "GROMOS":
        CYS_atoms = ["N", "CA", "CB", "SG", "C", "O", "H"]
    elif Config.pypka_params["ff_family"] == "CHARMM":
        CYS_atoms = ["N", "CA", "CB", "SG", "C", "O", "HN", "HA", "HB1", "HB2"]

    if set(res_atoms).issubset(CYS_atoms) and set(CYS_atoms).issubset(res_atoms):
        warn = "CYS-{0} is assumed to be participating in a SS-bond".format(resnumb)
        logger.warning(warn)
        return True
    else:
        return False


def check_site(molecule, resnumb, sites, resname, cur_atoms, ter=None):
    # Corrects weird tautomer residue names like GL4
    if resname not in PROTEIN_RESIDUES:
        for res in REGULARTITRATINGRES:
            if res[0:2] == resname[0:2]:
                resname = res

    is_titrating = True
    if ter:
        if not Config.pypka_params["ser_thr_titration"] and resname in (
            "SER",
            "THR",
        ):
            is_titrating = False
        else:
            is_titrating = bool(resname in TITRABLERESIDUES)

    elif not Config.pypka_params["ser_thr_titration"] and resname in (
        "SER",
        "THR",
    ):
        return

    res_atoms = copy(cur_atoms)
    (integrity_terminal, integrity_site) = check_integrity(
        resname, res_atoms, is_titrating, ter=ter
    )

    if integrity_terminal:
        ter_resnumb = resnumb + TERMINAL_OFFSET
        addSiteAndTautomers(molecule, ter_resnumb, ter, termini_resname=resname)
        if ter == "NTR":
            molecule.NTR += [resnumb]
        elif ter == "CTR":
            molecule.CTR += [resnumb]
    elif ter:
        warn_integrity_check_failed(resnumb, ter)

    # print(resname, resnumb, integrity_site, resnumb in molecule.CYS_bridges)

    if resnumb in sites:
        if integrity_site:
            addSiteAndTautomers(molecule, resnumb, resname)
            if resname == "CYS":
                molecule.correct_names[resnumb] = "CY0"
        else:
            warn_integrity_check_failed(resnumb, resname)
    elif resnumb in molecule.CYS_bridges or resname == "CYS":
        if ter == "NTR":
            cur_atoms = set(cur_atoms) - set(["H1", "H2", "H3"])
            cur_atoms.add("H")
        elif ter == "CTR":
            cur_atoms = set(cur_atoms) - set(
                ["HO11", "HO12", "HO21", "HO22", "O1", "O2"]
            )
            cur_atoms.update(("H", "O"))

        is_cys_bridge = check_CYS_bridge(cur_atoms, resnumb)
        if not is_cys_bridge and resnumb in molecule.CYS_bridges:
            warn_integrity_check_failed(resnumb, resname)
        elif is_cys_bridge and Config.pypka_params["ff_family"] == "CHARMM":
            molecule.correct_names[resnumb] = "SSB"


def check_sites_integrity(filename, molecules, chains_res):
    """Identifies titrable residues and checks integrity of the residue blocks
    (excluding Hydrogens)
    """
    resnumb = None
    cur_atoms = []
    prev_resnumb = None
    prev_resname = None
    last_chain = None
    chain = None
    with open(filename) as f:
        nline = 0
        f_lines = f.readlines()
        maxnlines = len(f_lines)
        for line in f_lines:
            resname = None
            nline += 1

            if "ATOM " == line[0:5]:
                (aname, anumb, resname, chain, resnumb, x, y, z) = read_pdb_line(line)

                if last_chain in molecules or (chain in molecules and not last_chain):
                    if not last_chain:
                        last_chain = chain
                    molecule = molecules[last_chain]
                    sites = chains_res[last_chain]
                    last_molecule = molecule

                if nline == maxnlines:
                    cur_atoms.append(aname)

            if line == "TER\n":
                resnumb += 1

            if (
                prev_resnumb != resnumb or nline == maxnlines or chain != last_chain
            ) and prev_resnumb is not None:

                if nline == maxnlines:
                    prev_resnumb = copy(resnumb)
                    resnumb = "None"

                if last_chain in molecules:
                    if (
                        prev_resname in TITRABLERESIDUES
                        or (prev_resnumb in molecule.NTR or resnumb in molecule.NTR)
                        or (prev_resnumb in molecule.CTR or resnumb in molecule.CTR)
                    ):
                        if prev_resnumb in molecule.NTR and resnumb not in molecule.NTR:
                            check_site(
                                molecule,
                                prev_resnumb,
                                sites,
                                prev_resname,
                                cur_atoms,
                                ter="NTR",
                            )
                            prev_resnumb = None
                        # Dealing with the last residue and CTR
                        elif (
                            prev_resnumb in molecule.CTR and resnumb not in molecule.CTR
                        ) or (prev_resnumb in molecule.CTR and chain != last_chain):
                            check_site(
                                molecule,
                                prev_resnumb,
                                sites,
                                prev_resname,
                                cur_atoms,
                                ter="CTR",
                            )
                            prev_resnumb = None
                            last_chain = None
                        # Dealing with the previous residue
                        elif (
                            prev_resnumb is not None
                            and prev_resname in TITRABLERESIDUES
                        ):
                            check_site(
                                molecule,
                                prev_resnumb,
                                sites,
                                prev_resname,
                                cur_atoms,
                            )

                    elif prev_resname == "ALA":
                        # TODO: check residue block integrity for other non titrating residues
                        pass
                elif (
                    last_molecule
                    and prev_resnumb in last_molecule.CTR
                    and resnumb not in last_molecule.CTR
                ):
                    check_site(
                        molecule,
                        prev_resnumb,
                        sites,
                        prev_resname,
                        cur_atoms,
                        ter="CTR",
                    )
                elif prev_resname in TITRABLERESIDUES and prev_resnumb is not None:
                    check_site(molecule, prev_resnumb, sites, prev_resname, cur_atoms)

            # Dealing with the new residue
            if prev_resnumb != resnumb:
                cur_atoms = [aname]
                prev_resnumb = resnumb
                prev_resname = resname
                last_chain = chain
            elif resnumb is not None:
                cur_atoms.append(aname)
                if prev_resname in ("NTR", "CTR") and prev_resname != resname:
                    prev_resname = resname

    for molecule in molecules.values():
        # Adding the reference tautomer to each site
        molecule.addReferenceTautomers()
        # Assigning a charge set to each tautomer
        molecule.addTautomersChargeSets()

    # TODO: report blocks that failed the check (in .log file with
    # numbering reference to stepwise scheme)
    # TODO: add lipid residues
    if Config.debug:
        print("exiting check_sites_integrity")


def check_integrity(resname, res_atoms, is_titrating, ter=None):
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
                if Config.debug:
                    print((aname, "not in", resname))
            else:
                anames2.remove(aname)
        if trigger:
            return anames2, False
        else:
            return anames2, True

    if Config.debug:
        print("###### INTEGRITY CHECK ######")
        print(resname)
        print(res_atoms)
    integrity = None
    integrity_ter = None

    if ter:
        st = "{}/{}tau1.st".format(Config.pypka_params["sts_path"], ter)
        res_atoms_st = read_anames(st)
        res_atoms, integrity_ter = pop_atoms(res_atoms_st, res_atoms)

    if is_titrating:
        main_chain_atoms = main_chains[Config.pypka_params["ff_family"]]
        if ter == "NTR":
            main_chain = main_chain_atoms["NTR"]
        elif ter == "CTR":
            main_chain = main_chain_atoms["CTR"]
        else:
            main_chain = main_chain_atoms["main"]

        for aname in main_chain:
            if aname in res_atoms:
                res_atoms.remove(aname)
            elif not ter:
                integrity = False

        if Config.debug:
            print(("i", integrity))

        if integrity is not False:
            filename = "{0}/{1}/sts/{2}tau1.st".format(
                Config.pypka_params["ffs_dir"], Config.pypka_params["ffID"], resname
            )
            res_atoms_st = read_anames(filename)
            res_atoms, integrity = pop_atoms(res_atoms_st, res_atoms)
            if Config.debug:
                print((res_atoms_st, res_atoms))
                print(("i", integrity))

    if len(res_atoms) != 0 and integrity:
        raise Exception("Something is wrong")

    return integrity_ter, integrity


def make_delphi_inputfile(f_in, f_out, molecules, fixed_sites):
    def getMaxCoords(coords, max_coords):
        x, y, z = coords
        max_x, max_y, max_z = max_coords
        if max_x < x:
            max_x = x
        if max_y < y:
            max_y = y
        if max_z < z:
            max_z = z
        return max_x, max_y, max_z

    def correct_termini(resnumb, resname, aname, ntr_res, ctr_res):
        if resnumb in ntr_res and aname in Config.pypka_params["NTR_atoms"]:
            resname = "NTR"
            resnumb += TERMINAL_OFFSET
        elif resnumb in ctr_res and aname in Config.pypka_params["CTR_atoms"]:
            resname = "CTR"
            resnumb += TERMINAL_OFFSET
            # if aname == "C":
            #    aname = "CT"
        return resnumb, resname, aname

    def correct_res_names(molecule, chain, resnumb, resname, aname):
        if resnumb in list(molecule.correct_names.keys()):
            resname = molecule.correct_names[resnumb]
        if (
            resnumb in list(molecule.correct_atoms.keys())
            and aname in molecule.correct_atoms[resnumb]
        ):
            aname = molecule.correct_atoms[resnumb][aname]

        if Config.pypka_params["ff_family"] == "CHARMM" and resname == "HIS":
            if (
                not fixed_sites
                or (
                    fixed_sites
                    and chain in fixed_sites
                    and resnumb not in fixed_sites[chain]
                )
                or (fixed_sites and chain not in fixed_sites)
            ):
                resname = "HSP"

        return resnumb, resname, aname

    def assign_atoms(sites, resnumb, aname, site_Hs, site_positions):
        ref_tau_name = resname
        # print(aname in list(sites[resnumb].getRefTautomer().charge_set.keys())
        if resnumb in list(sites.keys()) and aname in list(
            sites[resnumb].getRefTautomer().charge_set.keys()
        ):
            # ( aname not in ('N', 'H', 'C', 'O', 'CA') or
            # (aname in ('N', 'H', 'C', 'O', 'CA') and resname == 'NTR')):
            # change res name to reference tautomer
            ref_tau_name = sites[resnumb].getRefTautomerName()

            # add atom to corresponding site
            sites[resnumb].addAtom(aname, anumb)

            if chain not in site_positions:
                site_positions[chain] = {}
                site_Hs[chain] = {}
            if resnumb not in site_positions[chain]:
                site_positions[chain][resnumb] = []
                site_Hs[chain][resnumb] = []

            if resnumb in site_positions[chain]:
                site_positions[chain][resnumb].append((x, y, z))
                if aname[0] == "H":
                    site_Hs[chain][resnumb].append((x, y, z))

        return site_Hs, ref_tau_name, site_positions

    new_pdb_content = ""
    site_positions = {}
    site_Hs = {}
    max_box = [0.0, 0.0, 0.0]
    aposition = -1
    sequence = {}
    with open(f_in) as f:
        for line in f:
            if line.startswith("ATOM"):
                aposition += 1
                (aname, anumb, resname, chain, resnumb, x, y, z) = read_pdb_line(line)

                max_box = getMaxCoords([x, y, z], max_box)

                if chain not in sequence:
                    sequence[chain] = {}
                if resnumb not in sequence[chain]:
                    sequence[chain][resnumb] = resname

                if chain in molecules:
                    molecule = molecules[chain]
                    ntr_res = molecule.NTR
                    ctr_res = molecule.CTR
                    sites = molecule.sites

                    if (
                        resname == "HIS"
                        and aname == "HD1"
                        and resnumb not in sites.keys()
                    ):
                        aposition -= 1
                        continue

                    resnumb, resname, aname = correct_termini(
                        resnumb, resname, aname, ntr_res, ctr_res
                    )

                    resnumb, resname, aname = correct_res_names(
                        molecule, chain, resnumb, resname, aname
                    )

                    titrable_res = False
                    if resnumb in sites.keys():
                        titrable_res = True

                    molecule.addAtom(aname, anumb, aposition, titrable_res)

                    (site_Hs, resname, site_positions) = assign_atoms(
                        sites, resnumb, aname, site_Hs, site_positions
                    )

                else:
                    if resname == "HIS" and aname == "HD1":
                        aposition -= 1
                        continue
                    resnumb, resname, aname = correct_res_names(
                        molecule, chain, resnumb, resname, aname
                    )

                new_pdb_content += new_pdb_line(
                    aposition, aname, resname, resnumb, x, y, z, chain=chain
                )

            elif line.startswith("CRYST1"):
                parts = line.split()
                box = [float(i) for i in parts[1:4]]

    if box == [1.0, 1.0, 1.0]:
        box = max_box

    if Config.pypka_params["box"]:
        box = Config.pypka_params["box"]
    else:
        Config.pypka_params.setBox(box)

    if Config.delphi_params["pbc_dim"] == 2:
        Config.delphi_params.redefineScale()

    new_pdb_content += "TER\nENDMDL\n"
    with open(f_out, "w") as f_new:
        f_new.write(new_pdb_content)

    # TODO: check terminal_offset has to be bigger than the total number of residues
    # TODO: delete terminal_offset and use another approach to distinguish between N- and C-ter
    # TODO: check size xy > config.cutoff * 2
    # if so, raise Exception, and ask to change cutoff value

    # TODO: check if pbc_dim -> set gsizes from pdb size xy and ignore perfil

    for chain in site_positions.keys():
        molecule = molecules[chain]
        for site in site_positions[chain]:
            if site in list(molecule.sites.keys()):
                pos_max = [-9999990, -999999, -999999]
                pos_min = [999999, 999999, 999999]
                focus_center = [0, 0, 0]
                for atom in site_positions[chain][site]:
                    for i in range(3):
                        if pos_max[i] < atom[i]:
                            pos_max[i] = atom[i]
                        if pos_min[i] > atom[i]:
                            pos_min[i] = atom[i]
                focus_center[0] = (pos_max[0] + pos_min[0]) / 2
                focus_center[1] = (pos_max[1] + pos_min[1]) / 2
                focus_center[2] = (pos_max[2] + pos_min[2]) / 2

                if Config.delphi_params["pbc_dim"] == 2:
                    molecule.sites[site].addCenter(
                        focus_center, boxsize=box[0], box_z=box[2]
                    )
                else:
                    molecule.sites[site].addCenter(focus_center)
                hx, hy, hz = 0, 0, 0
                nHs = len(site_Hs[chain][site])
                if nHs == 0:
                    sitename = molecule.sites[site].getName()
                    raise Exception(
                        "Site {1}{0} appears "
                        "to have no Hydrogen atoms".format(site, sitename)
                    )
                for h in site_Hs[chain][site]:
                    hx += h[0]
                    hy += h[1]
                    hz += h[2]
                hx /= nHs
                hy /= nHs
                hz /= nHs
                Hcenter = (hx, hy, hz)
                molecule.sites[site].addCenterH(Hcenter)

    return sequence


def fix_fixed_sites(molecules, fixed_sites, fname):
    for chain in fixed_sites:
        for site, state in list(fixed_sites[chain].items()):
            if isinstance(site, str) and site[-1] in "NC":
                del fixed_sites[chain][site]
                site = int(site[:-1]) + TERMINAL_OFFSET
                fixed_sites[chain][site] = state

    for molecule in molecules.values():
        chain = molecule.chain
        for sitenumb, site in list(molecule.sites.items()):
            if sitenumb in fixed_sites[chain]:
                del molecule.sites[sitenumb]
                site_i = molecule.sites_order.index(site)
                del molecule.sites_order[site_i]

    new_pdb_content = ""
    with open(fname) as f:
        for line in f:
            if line.startswith("ATOM "):
                (aname, anumb, resname, chain, resnumb, x, y, z) = read_pdb_line(line)
                if resnumb in fixed_sites[chain]:
                    resname = "{}{}".format(
                        resname[:-1], str(fixed_sites[chain][resnumb])
                    )
                new_pdb_content += new_pdb_line(
                    anumb, aname, resname, resnumb, x, y, z, chain=chain
                )
            else:
                new_pdb_content += line
    with open(fname, "w") as f:
        f.write(new_pdb_content)
