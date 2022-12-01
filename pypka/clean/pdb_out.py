from pypka.config import Config
from pypka.constants import *
from pdbmender.postprocess import fix_structure_states


def write_output_structure(sites, molecules, delphi_input_content):
    outputname = Config.pypka_params["f_structure_out"]
    pH = float(Config.pypka_params["f_structure_out_pH"])
    ff_out = Config.pypka_params["ff_structure_out"]

    sites_ = {}
    for site in sites:
        resname = site.getName()
        resnumb = site.res_number
        molecule = site.molecule
        chain = molecule.chain
        termini_resname = site.termini_resname
        tit_curve = site.getTitrationCurve()
        most_prob_taut, (state_prob, taut_prob) = site.getMostProbTaut(pH)
        if chain not in sites_:
            sites_[chain] = []
        sites_[chain].append(
            (
                resname,
                resnumb,
                termini_resname,
                most_prob_taut,
                state_prob,
                taut_prob,
                tit_curve,
            )
        )

    tit_atoms = []
    other_atoms = []
    cys_bridges = {}
    renumbered = {}
    for molecule in molecules.values():
        cys_numbs = molecule.getCYS_bridges()
        cys_bridges[chain] = cys_numbs
        renumbered[chain] = molecule.icodes

        for atom_numb in molecule.atoms_tit_res:
            if molecule.atoms_tit_res[atom_numb]:
                tit_atoms.append(atom_numb)
            else:
                other_atoms.append(atom_numb)

    fix_structure_states(
        outputname,
        pH,
        ff_out,
        sites_,
        tit_atoms,
        other_atoms,
        cys_bridges,
        delphi_input_content,
        renumbered,
        TERMINAL_OFFSET,
        Config.pypka_params["pdb2pqr_inputfile"]
    )
