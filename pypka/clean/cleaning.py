import os
from copy import copy

from pypka.config import Config
from pypka.constants import *

from pypka.clean.ffconverter import *
from pypka.clean.formats import (
    correct_names,
    new_gro_line,
    new_pdb_line,
    new_pqr_line,
    read_gro_line,
    read_pdb_line,
    read_pqr_line,
)
from pypka.clean.utils import add_tautomers, rinse_pdb


def inputPDBCheck(filename, sites, clean_pdb):
    """
    Returns: chains_length, chains_res
    """
    if filename[-3:] in ("pdb", "pqr"):
        filetype = "pdb"
    elif filename[-3:] == "gro":
        filetype = "gro"
    else:
        raise Exception("Input file must be either a pdb or a gro.")

    chains_length = {}
    chains_res = {chain: {} for chain in sites.keys()}
    chains = []
    if filetype == "pdb" and not clean_pdb:
        new_gro_header = "CREATED within PyPka\n"
        new_gro_body = ""
    with open(filename) as f:
        last_chain = ""
        chain_length = 0

        nline = 0
        maxnlines = 0
        atom_number = 0
        for line in f:
            nline += 1
            atom_line = False
            if filetype == "pdb":
                if line.startswith("ATOM "):
                    atom_line = True
                    chain_length += 1
                    (aname, anumb, resname, chain, resnumb, x, y, z) = read_pdb_line(
                        line
                    )
                    if chain not in chains:
                        chains.append(chain)
                    atom_number += 1
                    if not clean_pdb:
                        if (
                            len(aname) > 2
                            and aname[1] == "H"
                            and aname[0] in ("1", "2")
                        ):
                            aname = aname[1:] + aname[0]
                        new_gro_body += new_gro_line(
                            anumb, aname, resname, resnumb, x / 10.0, y / 10, z / 10
                        )
                elif line.startswith("CRYST1"):
                    tmp = line.split()[1:4]
                    box = (float(tmp[0]), float(tmp[1]), float(tmp[2]))
                    new_gro_footer = "{0:10.5f}{1:10.5f}{2:10.5f}\n".format(
                        box[0] / 10.0, box[1] / 10.0, box[2] / 10.0
                    )

            elif filetype == "gro":
                if nline > 2 and nline < maxnlines:
                    (aname, anumb, resname, resnumb, x, y, z) = read_gro_line(line)
                    chain = "A"
                    atom_line = True
                elif nline == 2:
                    natoms = int(line.strip())
                    maxnlines = natoms + 3

            if atom_line:
                if chain_length == 1:
                    last_chain = chain

                if chain != last_chain and chain_length != 1:
                    chains_length[last_chain] = chain_length
                    # chains_res[chain] = done[chain]
                    chain_length = 0
                    last_chain = chain

                if (
                    chain in sites
                    and resnumb not in chains_res[chain]
                    and str(resnumb) in sites[chain]
                ):
                    if Config.pypka_params["ffinput"] == "CHARMM":
                        if resname in ("HSD", "HSE", "HSP"):
                            resname = "HIS"
                    chains_res[chain][resnumb] = resname

    # if filetype == 'pdb' and not clean_pdb:
    #    new_gro_header += '{0}\n'.format(atom_number)
    #    with open('TMP.gro', 'w') as f:
    #        f.write(new_gro_header + new_gro_body + new_gro_footer)

    chains_length[last_chain] = chain_length
    # chains_res[chain] = done[chain]

    # tmp_chains_res is an ugly hack so that test cases hold
    # TODO: remove tmp_chains_res and update tests
    tmp_chains_res = {}
    for chain in sites.keys():
        if chain not in chains:
            continue

        tmp_chains_res[chain] = {}
        for site in sites[chain]:
            if site[-1] == "C":
                resnumb = site[:-1]
                tmp_chains_res[chain][resnumb] = "CTR"
            elif site[-1] == "N":
                resnumb = site[:-1]
                tmp_chains_res[chain][resnumb] = "NTR"
        chains_res[chain] = {**tmp_chains_res[chain], **chains_res[chain]}

    skipped_sites = {}
    for chain, resnumbs in sites.items():
        if chain not in chains_res:
            skipped_sites[chain] = resnumbs
            continue
        for resnumb in resnumbs:
            skipped = False
            if resnumb[-1] in "NC":
                termini = "NTR" if resnumb[-1] == "N" else "CTR"
                resnumb = resnumb[:-1]
                if (
                    resnumb not in chains_res[chain]
                    or chains_res[chain][resnumb] != termini
                ):
                    print(
                        "{1} in chain '{0}' not found or not titratable.".format(
                            chain, termini
                        )
                    )
            else:
                resnumb = int(resnumb)
                if resnumb not in chains_res[chain].keys():
                    print(
                        "Residue #{1} in chain '{0}' not found or not titratable.".format(
                            chain, resnumb
                        )
                    )

    nsites = sum([len(res) for res in chains_res.values()])
    if not nsites:
        raise Exception(
            "No titrable residues found. Please check the residue number and chain."
        )

    return chains_length, chains_res


def cleanPDB(molecules, chains_res, inputpqr, outputpqr, automatic_sites):

    pdb_filename = Config.pypka_params["f_in"]

    remove_membrane_n_rna(pdb_filename, Config.pypka_params["pdb2pqr_inputfile"])
    if " " in molecules.keys():
        molecules["_"] = molecules[" "]
        del molecules[" "]
        chains_res["_"] = chains_res[" "]
        del chains_res[" "]

    logfile = "LOG_pdb2pqr"

    # CTR O1/O2 will be deleted and a O/OXT will be added
    # CYS will be turned into CY0
    # All HIS will become fully protonated
    rinse_pdb(
        Config.pypka_params["pdb2pqr_inputfile"],
        inputpqr,
        Config.pypka_params["ffinput"],
        Config.pypka_params["ff_family"],
        logfile=logfile,
        hopt=Config.pypka_params["pdb2pqr_h_opt"],
    )

    CYS_bridges = get_cys_bridges(logfile, molecules)

    if automatic_sites:
        chains_res = identify_cter(inputpqr, molecules, chains_res)

    if Config.pypka_params["f_structure_out"]:
        ff_out = Config.pypka_params["ff_structure_out"]
        if ff_out == "gromos_cph":
            ff_out = "gromos"
        rinse_pdb(
            Config.pypka_params["pdb2pqr_inputfile"],
            "Hs.pqr",
            ff_out,
            ff_out,
            logfile=logfile,
        )

    if Config.pypka_params["f_structure_out"]:
        with open("Hs.pqr") as f:
            for line in f:
                if line.startswith("ATOM "):
                    (
                        aname,
                        anumb,
                        resname,
                        chain,
                        resnumb,
                        x,
                        y,
                        z,
                        charge,
                        radius,
                    ) = read_pqr_line(line)

                    if chain not in chains_res:
                        continue
                    elif (
                        chain in chains_res
                        and resnumb in chains_res[chain]
                        and resname in AMBER_Hs
                        and aname in AMBER_Hs[resname]
                        and aname not in AMBER_mainchain_Hs
                    ):
                        continue
                    elif aname[0] == "H" and aname not in ("H1", "H2", "H3"):
                        if not chain in mainchain_Hs:
                            mainchain_Hs[chain] = {}
                        if not resnumb in mainchain_Hs[chain]:
                            mainchain_Hs[chain][resnumb] = []
                        mainchain_Hs[chain][resnumb].append(
                            (aname, anumb, resname, chain, x, y, z)
                        )

    new_pdb_text = ""
    removed_pdb_text = ""
    removed_pdb_lines = []
    # ntr_trigger = True
    resnumb_max = 0
    chains = list(molecules.keys())

    with open(inputpqr) as f:
        NUCLEIC_ACIDS = list(DNA_RESIDUES.values()) + list(RNA_RESIDUES.values())
        for line in f:
            if "ATOM" in line[:4]:
                (
                    aname,
                    anumb,
                    resname,
                    chain,
                    resnumb,
                    x,
                    y,
                    z,
                    charge,
                    radius,
                ) = read_pqr_line(line)

                termini_trigger = False

                if chain in chains:
                    molecule = molecules[chain]
                    NTR_numb = molecule.NTR
                    CTR_numb = molecule.CTR

                    sites_numbs = molecule.sites.keys()

                    aname, resname = correct_names(
                        resnumb, resname, aname, sites_numbs, NTR_numb, CTR_numb
                    )
                    resnumb_max = resnumb

                    # if resname == "CY0":
                    #    resname = "CYS"
                    #    if aname == "HG1":
                    #        continue

                    if resnumb in (NTR_numb, CTR_numb):
                        termini_trigger = True

                if resname == "CYS":
                    change_atoms = {"1CB": "CB", "1SG": "SG"}
                    if aname in change_atoms.keys():
                        aname = change_atoms[aname]
                if (
                    aname in ("O1", "O2", "OT1", "OT2", "H1", "H2", "H3")
                    and not termini_trigger
                    and resname in PROTEIN_RESIDUES
                ):
                    if aname == "O1":
                        aname = "O"
                    elif aname == "H1":
                        aname = "H"
                    else:
                        continue

                if line[26] != " ":
                    resnumb += TERMINAL_OFFSET

                new_line = new_pqr_line(
                    anumb,
                    aname,
                    resname,
                    resnumb,
                    x,
                    y,
                    z,
                    charge,
                    radius,
                    chain=chain,
                )
                if chain in molecules:
                    new_pdb_text += new_line
                elif (
                    aname not in ("O1", "O2", "OT1", "OT2", "H1", "H2", "H3")
                    or resname in NUCLEIC_ACIDS
                ):
                    if resname in ("SER", "THR") and aname == "HG1":
                        aname = "HG"
                    removed_pdb_lines.append(new_line)

    resnumb_old = resnumb_max + 1
    for line in removed_pdb_lines:
        (
            aname,
            anumb_old,
            resname,
            chain,
            resnumb,
            x,
            y,
            z,
            charge,
            radius,
        ) = read_pqr_line(line)
        anumb += 1
        resnumb += resnumb_old
        while resnumb < resnumb_max:
            resnumb += resnumb_old
        removed_pdb_text += new_pqr_line(
            anumb, aname, resname, resnumb, x, y, z, charge, radius, chain=chain
        )
        resnumb_max = resnumb
    #            if ntr_trigger and :
    #                config.tit_mole._NTR = resnumb

    # chains = copy(chains_res['A'])
    # for res in chains:
    #    if str(res) > config.terminal_offset:
    #        print((res, chains_res['A'][res]))
    #        del chains_res['A'][res]
    # chains_res['A'][config.tit_mole._NTR] = 'NTR'
    # chains_res['A'][config.tit_mole._CTR] = 'CTR'

    with open("cleaned.pqr", "w") as f_new:
        f_new.write(new_pdb_text)

    with open("removed.pqr", "w") as f_new:
        f_new.write(removed_pdb_text)

    sites_addHtaut = ""
    for chain in molecules:
        for res in chains_res[chain]:
            if chains_res[chain][res] == "NTR" or (
                chain in CYS_bridges and res in CYS_bridges[chain]
            ):
                continue
            sites_addHtaut += "{0}--{1}--{2},".format(
                res, chains_res[chain][res], chain.replace(" ", "_")
            )

    if len(sites_addHtaut) > 0 and sites_addHtaut[-1] == ",":
        sites_addHtaut = sites_addHtaut[:-1]

    if len(sites_addHtaut.strip()) == 0:
        os.system("cp cleaned.pqr {}".format(outputpqr))
    else:
        add_tautomers(sites_addHtaut, Config.pypka_params["ff_family"], outputpqr)

    with open(outputpqr) as f:
        content = f.read()
    with open(outputpqr, "w") as f_new:
        f_new.write(content + removed_pdb_text)

    rna_pqr = None
    # if rna_file:
    #    rna_pqr = rna_inputpqr
    final_pdb = "TMP.pdb"
    write_final_pdb(pdb_filename, outputpqr, final_pdb, rna_pqr)

    keep_membrane = False
    keep_ions = False
    if Config.delphi_params["pbc_dim"] == 2:
        keep_membrane = True

    if Config.pypka_params["keep_ions"]:
        keep_ions = True
    if keep_ions or keep_membrane:
        add_non_protein(
            pdb_filename, final_pdb, keep_membrane=keep_membrane, keep_ions=keep_ions
        )

    tmpfiles = (
        "LOG_pdb2pqr",
        "LOG_pdb2pqr_err",
        "clean.pqr",
        "cleaned.pqr",
        "cleaned_tau.pqr",
        "removed.pqr",
        "Hs.pqr",
    )
    for filename in tmpfiles:
        if not Config.debug and os.path.isfile(filename):
            os.remove(filename)
    if (
        not Config.debug
        and os.path.isfile(Config.pypka_params["pdb2pqr_inputfile"])
        and not Config.pypka_params.structure_output
    ):
        os.remove(Config.pypka_params["pdb2pqr_inputfile"])


def get_cys_bridges(f_pdb2pqr_log, molecules):
    CYS_bridges = {}
    with open(f_pdb2pqr_log) as f:
        trigger_error = ""
        for line in f:
            if "patched with CYX" in line:
                parts = line.split("patched")[0].replace("PATCH INFO: ", "").split()
                resname, chain, resnumb = parts
                # chain = chain.replace("_", " ")
                if not chain in CYS_bridges:
                    CYS_bridges[chain] = []
                cys_res_numb = int(parts[-1])
                if cys_res_numb not in CYS_bridges[chain]:
                    CYS_bridges[chain].append(cys_res_numb)
            if (
                "error" in line.lower()
                and "Error parsing line: invalid literal for int() with base 10"
                not in line
            ):
                trigger_error += line
        if trigger_error:
            raise Exception(
                "Found errors while parsing the input structure to PDB2PQR:\n"
                + trigger_error
            )

    for chain, molecule in molecules.items():
        if chain in CYS_bridges:
            molecule.saveCYSBridges(CYS_bridges[chain])
    return CYS_bridges


def identify_cter(inputpqr, molecules, chains_res):
    with open(inputpqr) as f:
        for line in f:
            if line.startswith("ATOM "):
                (aname, anumb, resname, chain, resnumb, x, y, z) = read_pdb_line(line)
                if (
                    "CTR" not in chains_res
                    and aname
                    in (
                        "CT",
                        "OT",
                        "OT1",
                        "OT2",
                        "O1",
                        "O2",
                        "OXT",
                    )
                    and chain in chains_res
                    and str(resnumb) not in chains_res[chain]
                    and resname in [*TITRABLERESIDUES, *PROTEIN_RESIDUES]
                ):
                    chains_res[chain][str(resnumb)] = "CTR"
                    molecules[chain].CTR = resnumb
                    sID = molecules[chain].addSite(resnumb + TERMINAL_OFFSET)
                    molecules[chain].addTautomers(sID, TITRABLETAUTOMERS["CTR"], "CTR")
    return chains_res


def remove_membrane_n_rna(pdbfile, outfile):
    protein_lines = ""
    to_remove = LIPID_RESIDUES + list(Config.pypka_params.LIPIDS.values())

    with open(pdbfile) as f:
        for line in f:
            if line.startswith("ATOM"):
                (aname, anumb, resname, chain, resnumb, x, y, z) = read_pdb_line(line)

                if chain == " ":
                    chain = "_"  # workaround to deal with pdb2pqr
                insertion_code = line[26].strip()
                if insertion_code:
                    continue
                if aname[0] == "H" and Config.pypka_params["remove_hs"]:
                    continue

                if Config.pypka_params["ffinput"] == "CHARMM":
                    if resname in ("HSD", "HSE"):
                        resname = "HSP"

                if resname not in to_remove:
                    if resname in PDB_RNA_RESIDUES:
                        resname = PDB_RNA_RESIDUES[resname]
                    protein_lines += new_pdb_line(
                        anumb, aname, resname, resnumb, x, y, z, chain=chain
                    )

            elif line.startswith("ENDMDL"):
                break

    with open(outfile, "w") as f_new:
        f_new.write(protein_lines)

    # rna_fname = None
    # if rna_lines:
    #    rna_fname = "tmp_rna.pdb"
    #    with open(rna_fname, "w") as f_new:
    #        f_new.write(rna_lines)

    return  # rna_fname


def write_final_pdb(pdb_filename, outputpqr, final_pdb, rna_inputpqr):
    def pqr2pdb(line, counter):
        counter += 1
        (
            aname,
            anumb_old,
            resname,
            chain,
            resnumb,
            x,
            y,
            z,
            charge,
            radius,
        ) = read_pqr_line(line)

        # if chain == "_":
        #    chain = " "

        if resname in RNA_RESIDUES:
            resname = RNA_RESIDUES[resname]

            incorrect_rna_ops = {"OP1": "O1P", "OP2": "O2P"}

            if aname in incorrect_rna_ops:
                aname = incorrect_rna_ops[aname]

        return (
            new_pdb_line(counter, aname, resname, resnumb, x, y, z, chain=chain),
            counter,
        )

    with open(pdb_filename) as f:
        box = None
        for line in f:
            if line.startswith("CRYST1"):
                tmp = line.split()[1:4]
                box = (float(tmp[0]), float(tmp[1]), float(tmp[2]))
        if not box:
            box = (1.0, 1.0, 1.0)

    with open(outputpqr) as f, open(final_pdb, "w") as fnew:
        counter = 0
        new_pdb_text = (
            "CRYST1  {:3.3f}  {:3.3f}  {:3.3f}  "
            "60.00  60.00  90.00 P 1           1\n".format(box[0], box[1], box[2])
        )
        for line in f:
            next_pdb_line, counter = pqr2pdb(line, counter)
            new_pdb_text += next_pdb_line

        if rna_inputpqr:
            with open(rna_inputpqr) as f_rna:
                for line in f_rna:
                    if line.startswith("ATOM "):
                        next_pdb_line, counter = pqr2pdb(line, counter)
                        new_pdb_text += next_pdb_line
        fnew.write(new_pdb_text)


def add_non_protein(pdbfile_origin, add_to_pdb, keep_membrane=False, keep_ions=False):
    new_file_body = ""

    with open(add_to_pdb) as f:
        for line in f:
            if line.startswith("ATOM "):
                (aname, anumb, resname, chain, resnumb, x, y, z) = read_pdb_line(line)
                last_anumb = anumb
                last_resnumb = resnumb

    # Read the original pdb with the membrane
    with open(pdbfile_origin) as f:
        for line in f:
            if "ATOM " == line[0:5]:
                (aname, anumb, resname, chain, resnumb, x, y, z) = read_pdb_line(line)
                if keep_membrane:
                    if resname in LIPID_RESIDUES:
                        last_anumb += 1
                        new_file_body += new_pdb_line(
                            last_anumb, aname, resname, resnumb, x, y, z, chain=" "
                        )

                    if resname in list(Config.pypka_params.LIPIDS.values()):
                        aname, resname, to_include = convert_FF_atomnames(
                            aname, resname
                        )
                        if to_include:
                            last_anumb += 1
                            resnumb += last_resnumb
                            new_file_body += new_pdb_line(
                                last_anumb,
                                aname,
                                resname,
                                resnumb,
                                x,
                                y,
                                z,
                                chain=" ",
                            )
                if keep_ions and aname in IONS and resname == aname:
                    last_anumb += 1
                    resnumb += last_resnumb
                    new_file_body += new_pdb_line(
                        last_anumb, aname, resname, resnumb, x, y, z, chain=chain
                    )

    with open(add_to_pdb, "a") as f_new:
        f_new.write(new_file_body)


def convert_FF_atomnames(aname, resname):
    popc_resname = Config.pypka_params.LIPIDS["POPC"]
    chol_resname = Config.pypka_params.LIPIDS["cholesterol"]

    lookup = [
        ((" N  ", popc_resname), ("NTM", "CHL")),
        ((" C12", popc_resname), ("CB ", "CHL")),
        ((" C13", popc_resname), ("CN1", "CHL")),
        ((" C14", popc_resname), ("CN2", "CHL")),
        ((" C15", popc_resname), ("CN3", "CHL")),
        ((" C11", popc_resname), ("CA ", "PJ2")),
        ((" P  ", popc_resname), ("P  ", "PJ2")),
        ((" O13", popc_resname), ("OG ", "PJ2")),
        ((" O14", popc_resname), ("OB ", "PJ2")),
        ((" O12", popc_resname), ("OA ", "PJ2")),
        ((" O11", popc_resname), ("OD ", "PJ2")),
        ((" C1 ", popc_resname), ("CD ", "POX")),
        ((" C2 ", popc_resname), ("CE ", "POX")),
        ((" O21", popc_resname), ("OE ", "POX")),
        ((" C21", popc_resname), ("C1A", "POX")),
        ((" O22", popc_resname), ("O1A", "POX")),
        ((" C22", popc_resname), ("C1B", "POX")),
        ((" C3 ", popc_resname), ("CZ ", "POX")),
        ((" O31", popc_resname), ("OZ ", "POX")),
        ((" C31", popc_resname), ("C2A", "POX")),
        ((" O32", popc_resname), ("O2A", "POX")),
        ((" C32", popc_resname), ("C2B", "POX")),
        ((" C23", popc_resname), ("C1C", "POX")),
        ((" C24", popc_resname), ("C1D", "POX")),
        ((" C25", popc_resname), ("C1E", "POX")),
        ((" C26", popc_resname), ("C1F", "POX")),
        ((" C27", popc_resname), ("C1G", "POX")),
        ((" C28", popc_resname), ("C1H", "POX")),
        ((" C29", popc_resname), ("C1I", "POX")),
        (("0C21", popc_resname), ("C1J", "POX")),
        (("1C21", popc_resname), ("C1K", "POX")),
        (("2C21", popc_resname), ("C1L", "POX")),
        (("3C21", popc_resname), ("C1M", "POX")),
        (("4C21", popc_resname), ("C1N", "POX")),
        (("5C21", popc_resname), ("C1O", "POX")),
        (("6C21", popc_resname), ("C1P", "POX")),
        (("7C21", popc_resname), ("C1Q", "POX")),
        (("8C21", popc_resname), ("C1R", "POX")),
        ((" C33", popc_resname), ("C2C", "POX")),
        ((" C34", popc_resname), ("C2D", "POX")),
        ((" C35", popc_resname), ("C2E", "POX")),
        ((" C36", popc_resname), ("C2F", "POX")),
        ((" C37", popc_resname), ("C2G", "POX")),
        ((" C38", popc_resname), ("C2H", "POX")),
        ((" C39", popc_resname), ("C2I", "POX")),
        (("0C31", popc_resname), ("C2J", "POX")),
        (("1C31", popc_resname), ("C2K", "POX")),
        (("2C31", popc_resname), ("C2L", "POX")),
        (("3C31", popc_resname), ("C2M", "POX")),
        (("4C31", popc_resname), ("C2N", "POX")),
        (("5C31", popc_resname), ("C2O", "POX")),
        (("6C31", popc_resname), ("C2P", "POX")),
        ((" H3", chol_resname), (" H3", "CHO")),
        ((" C3", chol_resname), (" C3", "CHO")),
        ((" C2", chol_resname), (" C2", "CHO")),
        ((" C1", chol_resname), (" C1", "CHO")),
        (("C10", chol_resname), ("C10", "CHO")),
        (("C19", chol_resname), ("C19", "CHO")),
        ((" C9", chol_resname), (" C9", "CHO")),
        (("C11", chol_resname), ("C11", "CHO")),
        (("C12", chol_resname), ("C12", "CHO")),
        (("C13", chol_resname), ("C13", "CHO")),
        (("C18", chol_resname), ("C18", "CHO")),
        (("C17", chol_resname), ("C17", "CHO")),
        (("C20", chol_resname), ("C20", "CHO")),
        (("C22", chol_resname), ("C22", "CHO")),
        (("C23", chol_resname), ("C23", "CHO")),
        (("C24", chol_resname), ("C24", "CHO")),
        (("C25", chol_resname), ("C25", "CHO")),
        (("C27", chol_resname), ("C27", "CHO")),
        (("C26", chol_resname), ("C26", "CHO")),
        (("C21", chol_resname), ("C21", "CHO")),
        (("C16", chol_resname), ("C16", "CHO")),
        (("C15", chol_resname), ("C15", "CHO")),
        (("C14", chol_resname), ("C14", "CHO")),
        ((" C8", chol_resname), (" C8", "CHO")),
        ((" C7", chol_resname), (" C7", "CHO")),
        ((" C6", chol_resname), (" C6", "CHO")),
        ((" H6", chol_resname), (" H6", "CHO")),
        ((" C5", chol_resname), (" C5", "CHO")),
        ((" C4", chol_resname), (" C4", "CHO")),
        ((" O3", chol_resname), (" O3", "CHO")),
        (("H3'", chol_resname), ("H3'", "CHO")),
    ]

    for swap in lookup:
        old_aname = swap[0][0].strip()
        old_resname = swap[0][1]
        new_aname = swap[1][0].strip()
        new_resname = swap[1][1]
        if old_aname == aname and old_resname == resname:
            aname = new_aname
            resname = new_resname
            return aname, resname, True

    return aname, resname, False
