import os
import logging
from pypka.clean.ffconverter import (
    AMBER_protomers,
    GROMOS_protomers,
    gromos2amber,
    mainchain_Hs,
)
from pypka.config import Config
from pypka.constants import *
from pdbmender.formats import new_pdb_line, read_pdb_line, read_pqr_line
from pdbmender.utils import mend_pdb

logger = logging.getLogger(__name__)


def write_output_structure(sites, molecules, delphi_input_content):
    def getProtomerResname(pdb_content, site, pH, ff_protomers):
        resnumb = site.getResNumber()
        resname = site.getName()
        new_state, (state_prob, taut_prob) = site.getMostProbTaut(pH)
        new_state_i = new_state - 1
        for ff_resname, protomers in ff_protomers[resname].items():
            if new_state_i in protomers.keys():
                new_resname = ff_resname
                remove_hs = protomers[new_state_i]
                average_prot = site.getTitrationCurve()[pH]
                if state_prob < 0.75:
                    warn = (
                        "{0}{1} "
                        "protonation state probability: {2}, "
                        "tautomer probability: {3}".format(
                            resname, resnumb, state_prob, taut_prob
                        )
                    )
                    logger.warn(warn)

                rounded_sprob = round(state_prob, 2)
                rounded_tprob = round(taut_prob, 2)
                rounded_avgprot = round(average_prot, 2)
                remark_line = (
                    "{0: <5}{1: <6}    {2: >1.2f}         "
                    "{3: >1.2f}         {4: >1.2f}".format(
                        resname,
                        resnumb,
                        rounded_avgprot,
                        rounded_sprob,
                        rounded_tprob,
                    )
                )
                pdb_content += "REMARK     {text}\n".format(text=remark_line)
        # print(resnumb, new_state, new_resname, remove_hs, state_prob, taut_prob)
        return pdb_content, new_state_i, new_resname, remove_hs

    outputname = Config.pypka_params["f_structure_out"]
    pH = float(Config.pypka_params["f_structure_out_pH"])
    ff_out = Config.pypka_params["ff_structure_out"]
    ff_protomer = {"amber": AMBER_protomers, "gromos_cph": GROMOS_protomers}[ff_out]
    pdb_content = (
        "REMARK     Protonation states assigned according to PypKa\n"
        "REMARK     Residue    Avg Prot   State Prob    Taut Prob\n"
    )
    new_states = {}
    for site in sites:
        resname = site.getName()
        resnumb = site.res_number
        molecule = site.molecule
        chain = molecule.chain
        (pdb_content, new_state, new_resname, remove_hs) = getProtomerResname(
            pdb_content, site, pH, ff_protomer
        )
        if resname in ("NTR", "CTR"):
            new_resname = site.termini_resname
        if chain not in new_states:
            new_states[chain] = {}
        new_states[chain][resnumb] = (resname, new_state, new_resname, remove_hs)
    new_pdb = pdb_content
    counter = 0
    tit_atoms = {}
    other_atoms = {}
    for molecule in molecules.values():
        for atom_numb in molecule.atoms_tit_res:
            if molecule.atoms_tit_res[atom_numb]:
                tit_atoms[atom_numb] = molecule
            else:
                other_atoms[atom_numb] = molecule

    in_delphi_pdb = {}
    for line in delphi_input_content:
        if line.startswith("ATOM "):
            (aname, anumb, resname, chain, resnumb, x, y, z) = read_pdb_line(line)
            if chain not in in_delphi_pdb:
                in_delphi_pdb[chain] = {}
            if resnumb not in in_delphi_pdb[chain]:
                in_delphi_pdb[chain][resnumb] = []
            in_delphi_pdb[chain][resnumb].append(aname)

    for line in delphi_input_content:
        if line.startswith("ATOM "):
            (aname, anumb, resname, chain, resnumb, x, y, z) = read_pdb_line(line)

            if anumb in tit_atoms.keys():
                molecule = tit_atoms[anumb]
                (oldresname, new_state, resname, removeHs) = new_states[chain][resnumb]
                if aname in removeHs:
                    continue
                if (
                    ff_out == "amber"
                    and oldresname in gromos2amber
                    and new_state in gromos2amber[oldresname]
                    and aname in gromos2amber[oldresname][new_state]
                ):
                    aname = gromos2amber[oldresname][new_state][aname]
            elif anumb in other_atoms:
                molecule = other_atoms[anumb]
            else:
                continue
            if resnumb > TERMINAL_OFFSET:
                termini_site = molecule.sites[resnumb]
                resnumb -= TERMINAL_OFFSET
                if resnumb in molecule.sites.keys():
                    _, ter_new_state, resname, ter_removeHs = new_states[chain][resnumb]
                else:
                    resname = termini_site.termini_resname
                # print(new_pdb_line(anumb, aname, resname, resnumb, x, y, z).strip())
            if resnumb in molecule.getCYS_bridges():
                resname = "CYX"
            counter += 1
            new_pdb += new_pdb_line(
                counter, aname, resname, resnumb, x, y, z, chain=chain
            )
            if chain in mainchain_Hs and resnumb in mainchain_Hs[chain]:
                while len(mainchain_Hs[chain][resnumb]) > 0:
                    counter += 1
                    (aname, anumb, oldresname, chain, x, y, z) = mainchain_Hs[chain][
                        resnumb
                    ].pop()
                    if (
                        resnumb not in in_delphi_pdb[chain]
                        or aname not in in_delphi_pdb[chain][resnumb]
                    ):
                        new_pdb += new_pdb_line(
                            counter, aname, resname, resnumb, x, y, z, chain=chain
                        )
                del mainchain_Hs[chain][resnumb]
        elif not line.startswith("ENDMDL"):
            new_pdb += line

    outputpqr = "leftovers.pqr"
    logfile = "LOG_pdb2pqr_nontitrating"
    if ff_out == "gromos_cph":
        ff_out = "GROMOS"
    mend_pdb(
        Config.pypka_params["pdb2pqr_inputfile"],
        outputpqr,
        ff_out,
        ff_out,
        logfile=logfile,
    )
    os.system("rm -f input_clean_fixed.pdb")

    with open(outputpqr) as f:
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
                if chain not in mainchain_Hs:
                    counter += 1
                    new_pdb += new_pdb_line(
                        counter, aname, resname, resnumb, x, y, z, chain=chain
                    )
    to_remove = (logfile, outputpqr, Config.pypka_params["pdb2pqr_inputfile"])
    for f in to_remove:
        os.remove(f)
    with open(outputname, "w") as f_new:
        f_new.write(new_pdb)
