import os
import log
from formats import (read_pdb_line, read_pqr_line, read_gro_line,
                     correct_names, new_pqr_line, new_pdb_line, new_gro_line, pdb2gro)

from config import Config
from constants import *
from ffconverter import *
from copy import copy

def inputPDBCheck(filename, sites, clean_pdb):
    """
    Returns: chains_length, chains_res
    """
    if filename[-3:] in ('pdb', 'pqr'):
        filetype = 'pdb'
    elif filename[-3:] == 'gro':
        filetype = 'gro'
    else:
        raise Exception('Input file must be either a pdb or a gro.')

    chains_length = {}
    chains_res = {}

    for chain in sites.keys():
        chains_res[chain] = {}
        for site in sites[chain]:
            if site[-1] == 'C':
                resnumb = site[:-1]
                chains_res[chain][resnumb] = 'CTR'
            elif site[-1] == 'N':
                resnumb = site[:-1]
                chains_res[chain][resnumb] = 'NTR'

    if filetype == 'pdb' and not clean_pdb:
        new_gro_header = 'CREATED within PyPka\n'
        new_gro_body = ''
    with open(filename) as f:
        last_chain = ''
        chain_length = 0

        nline = 0
        maxnlines = 0
        atom_number = 0
        for line in f:
            nline += 1
            atom_line = False
            if filetype == 'pdb':
                if 'ATOM ' == line[0:5]:
                    atom_line = True
                    chain_length += 1
                    (aname, anumb, resname,
                     chain, resnumb, x, y, z) = read_pdb_line(line)
                    atom_number += 1
                    if not clean_pdb:
                        if len(aname) > 2 and \
                           aname[1] == 'H' and \
                           aname[0] in ('1', '2'):
                            aname = aname[1:] + aname[0]
                        new_gro_body += new_gro_line(anumb, aname,
                                                     resname, resnumb,
                                                     x / 10.0, y / 10, z / 10)
                elif 'CRYST1' in line:
                    tmp = line.split()[1:4]
                    box = (float(tmp[0]), float(tmp[1]), float(tmp[2]))
                    new_gro_footer = '{0:10.5f}{1:10.5f}{2:10.5f}\n'.format(box[0] / 10.0,
                                                                            box[1] / 10.0,
                                                                            box[2] / 10.0)

            elif filetype == 'gro':
                if nline > 2 and nline < maxnlines:
                    (aname, anumb, resname, resnumb, x, y, z) = read_gro_line(line)
                    chain = 'A'
                    atom_line = True
                elif nline == 2:
                    natoms = int(line.strip())
                    maxnlines = natoms + 3

            if atom_line:
                if chain_length == 1:
                    last_chain = chain

                if chain != last_chain and chain_length != 1:
                    chains_length[last_chain] = chain_length
                    #chains_res[chain] = done[chain]
                    chain_length = 0
                    last_chain = chain

                if chain in sites and \
                   resnumb not in chains_res[chain] and \
                   str(resnumb) in sites[chain]:
                    chains_res[chain][resnumb] = resname

    #if filetype == 'pdb' and not clean_pdb:
    #    new_gro_header += '{0}\n'.format(atom_number)
    #    with open('TMP.gro', 'w') as f:
    #        f.write(new_gro_header + new_gro_body + new_gro_footer)

    chains_length[last_chain] = chain_length
    #chains_res[chain] = done[chain]

    return chains_length, chains_res


def cleanPDB(molecules, chains_res, inputpqr, outputpqr):
    """
    """
    pdb_filename = Config.pypka_params['f_in']
    pdb2pqr_path = Config.pypka_params['pdb2pqr']
    inputpdbfile = removeMembrane(pdb_filename)

    logfile = 'LOG_pdb2pqr'
    #errfile = 'LOG_pdb2pqr_err'

    #log.redirectOutput("start", logfile)
    #log.redirectErr("start", errfile)

    # CTR O1/O2 will be deleted and a O/OXT will be added
    os.system('python2 {0} {1} {2} --ff {4} --ffout {4} '
              '--drop-water -v --chain > {5} 2>&1 '.format(pdb2pqr_path,
                                                           inputpdbfile, inputpqr,
                                                           Config.pypka_params['ffinput'],
                                                           Config.pypka_params['ff_family'],
                                                           logfile))
    if Config.pypka_params['f_structure_out']:
        ff_out = Config.pypka_params['ff_structure_out']
        if ff_out == 'gromos_cph':
            ff_out = 'gromos'
        os.system('python2 {0} {1} Hs.pqr --ff {3} --ffout {3} '
                  '--drop-water -v --chain >> {2} 2>&1 '.format(pdb2pqr_path,
                                                                inputpdbfile,
                                                                logfile,
                                                                ff_out))

    #log.redirectOutput("stop", logfile)
    #log.redirectErr("stop", errfile)

    if Config.pypka_params['f_structure_out']:
        with open('Hs.pqr') as f:
            for line in f:
                if line.startswith('ATOM '):
                    (aname, anumb, resname, chain, resnumb, x, y,
                     z, charge, radius) = read_pqr_line(line)

                    if resnumb in chains_res['A'] and \
                       resname in AMBER_Hs and \
                       aname in AMBER_Hs[resname] and \
                       aname not in AMBER_mainchain_Hs:
                        continue
                    elif aname[0] == 'H' and aname not in ('H1', 'H2', 'H3'):
                        if not resnumb in mainchain_Hs:
                            mainchain_Hs[resnumb] = []
                        mainchain_Hs[resnumb].append((aname, anumb, resname, chain,
                                                      x, y, z))

    CYS_bridges = {}
    with open('LOG_pdb2pqr') as f:
        for line in f:
            if 'patched with CYX' in line:
                parts = line.split('patched')[0].replace('PATCH INFO: ', '').split()
                resname, chain, resnumb = parts
                chain = chain.replace('_', ' ')
                if not chain in CYS_bridges:
                    CYS_bridges[chain] = []
                cys_res_numb = int(parts[-1])
                if cys_res_numb not in CYS_bridges[chain]:
                    CYS_bridges[chain].append(cys_res_numb)

    for chain, molecule in molecules.items():
        if chain in CYS_bridges:
            molecule.saveCYSBridges(CYS_bridges[chain])

    new_pdb_text = ''
    removed_pdb_text = ''
    removed_pdb_lines = []
    #ntr_trigger = True
    resnumb_max = 0
    chains = list(molecules.keys())

    with open(inputpqr) as f:
        for line in f:
            if "ATOM" in line[:4]:
                (aname, anumb, resname, original_chain, resnumb, x, y,
                 z, charge, radius) = read_pqr_line(line)

                if original_chain == '_':
                    chain = ' '
                else:
                    chain = original_chain

                termini_trigger = False
                if chain in chains:
                    molecule = molecules[chain]
                    NTR_numb = molecule.NTR
                    CTR_numb = molecule.CTR

                    sites_numbs = molecule.sites.keys()

                    aname, resname = correct_names(resnumb,
                                                   resname, aname,
                                                   sites_numbs, NTR_numb, CTR_numb)
                    resnumb_max = resnumb

                    if resnumb in (NTR_numb, CTR_numb):
                        termini_trigger = True

                if aname in ('O1', 'O2',
                             'OT1', 'OT2',
                             'H1', 'H2', 'H3') and \
                   not termini_trigger and resname in PROTEIN_RESIDUES:
                    if aname == 'O1':
                        aname = 'O'
                    elif aname == 'H1':
                        aname = 'H'
                    else:
                        continue

                if line[26] != ' ':
                    resnumb += TERMINAL_OFFSET

                new_line = new_pqr_line(anumb, aname, resname,
                                        resnumb, x, y, z, charge, radius, chain=original_chain)
                if chain in molecules:
                    new_pdb_text += new_line
                elif aname not in ('O1', 'O2',
                                   'OT1', 'OT2',
                                   'H1', 'H2', 'H3'):
                    if resname in ('SER', 'THR') and aname == 'HG1':
                        aname = 'HG'
                    removed_pdb_lines.append(new_line)

    resnumb_old = resnumb_max + 1
    for line in removed_pdb_lines:
        (aname, anumb_old, resname, chain, resnumb, x, y,
         z, charge, radius) = read_pqr_line(line)
        anumb += 1
        resnumb += resnumb_old
        while resnumb < resnumb_max:
            resnumb += resnumb_old
        removed_pdb_text += new_pqr_line(anumb, aname, resname,
                                         resnumb, x, y, z, charge, radius)
        resnumb_max = resnumb
    #            if ntr_trigger and :
    #                config.tit_mole._NTR = resnumb

    #chains = copy(chains_res['A'])
    #for res in chains:
    #    if str(res) > config.terminal_offset:
    #        print((res, chains_res['A'][res]))
    #        del chains_res['A'][res]
    #chains_res['A'][config.tit_mole._NTR] = 'NTR'
    #chains_res['A'][config.tit_mole._CTR] = 'CTR'

    with open('cleaned.pqr', 'w') as f_new:
        f_new.write(new_pdb_text)

    with open('removed.pqr', 'w') as f_new:
        f_new.write(removed_pdb_text)

    sites_addHtaut = ''
    for chain in molecules:
        for res in chains_res[chain]:
            if chains_res[chain][res] == 'NTR' or \
               (chain in CYS_bridges and res in CYS_bridges[chain]):
               continue
            sites_addHtaut += '{0}-{1}-{2},'.format(res, chains_res[chain][res],
                                                    chain.replace(' ', '_'))

    if len(sites_addHtaut) > 0 and sites_addHtaut[-1] == ',':
        sites_addHtaut = sites_addHtaut[:-1]

    logfile = 'LOG_addHtaut'
    #log.redirectErr("start", logfile)

    script_dir = Config.pypka_params['script_dir']
    os.system('{}/addHtaut_{} cleaned.pqr {} > {} 2> {}'.format(script_dir,
                                                                Config.pypka_params['ff_family'],
                                                                sites_addHtaut,
                                                                outputpqr, logfile))

    #log.redirectErr("stop", logfile)
    log.checkDelPhiErrors(logfile, 'addHtaut')

    with open(outputpqr) as f:
        content = f.read()
    with open(outputpqr, 'w') as f_new:
        f_new.write(content + removed_pdb_text)

    with open(pdb_filename) as f:
        box = None
        for line in f:
            if line.startswith('CRYST1'):
                tmp = line.split()[1:4]
                box = (float(tmp[0]), float(tmp[1]), float(tmp[2]))
        if not box:
            box = (1.0, 1.0, 1.0)

    #pdb2gro(outputpqr, "TMP.gro", chains_res, box, pqr=True, fix_termini=False)

    with open(outputpqr) as f,\
         open('TMP.pdb', 'w') as fnew:
        counter = 0
        new_pdb_text = "CRYST1  {:3.3f}  {:3.3f}  {:3.3f}  "\
                       "60.00  60.00  90.00 P 1           1\n".format(box[0],
                                                                      box[1], box[2])
        for line in f:
            counter += 1
            (aname, anumb_old, resname, chain, resnumb, x, y,
             z, charge, radius) = read_pqr_line(line)

            if chain == '_':
                chain = ' '

            new_pdb_text += new_pdb_line(counter, aname, resname, resnumb, x, y, z, chain=chain)
        fnew.write(new_pdb_text)

    keep_membrane = False
    keep_ions = False
    if Config.delphi_params['pbc_dim'] == 2:
        keep_membrane = True

    if Config.pypka_params['keep_ions']:
        keep_ions = True
    if keep_ions or keep_membrane:
        add_non_protein(pdb_filename, 'TMP.pdb',
                        keep_membrane=keep_membrane, keep_ions=keep_ions)

    tmpfiles = ('LOG_pdb2pqr', 'LOG_pdb2pqr_err', 'clean.pqr', 'cleaned.pqr',
                'cleaned_tau.pqr', 'input_clean.pdb', 'removed.pqr', 'Hs.pqr')
    for filename in tmpfiles:
        if not Config.debug and os.path.isfile(filename):
            os.remove(filename)

def removeMembrane(pdbfile):
    nomembrane_text = ''
    with open(pdbfile) as f:
        for line in f:
            if 'ATOM ' == line[0:5]:
                (aname, anumb, resname,
                 chain, resnumb, x, y, z) = read_pdb_line(line)

                if chain == ' ':
                    chain = '_' # workaround to deal with pdb2pqr

                if resname not in LIPID_RESIDUES:
                    nomembrane_text += new_pdb_line(anumb, aname, resname, resnumb,
                                                    x, y, z, chain=chain)

            else:
                nomembrane_text += line
    with open('tmp.tmp', 'w') as f_new:
        f_new.write(nomembrane_text)
    os.rename('tmp.tmp', 'input_clean.pdb')
    return 'input_clean.pdb'


def add_non_protein(pdbfile_origin, add_to_pdb, keep_membrane=False, keep_ions=False):
    new_file_body = ''

    with open(add_to_pdb) as f:
        for line in f:
            if line.startswith('ATOM '):
                (aname, anumb, resname,
                 chain, resnumb, x, y, z) = read_pdb_line(line)
                last_anumb = anumb
                last_resnumb = resnumb

    # Read the original pdb with the membrane
    with open(pdbfile_origin) as f:
        for line in f:
            if 'ATOM ' == line[0:5]:
                (aname, anumb, resname,
                 chain, resnumb, x, y, z) = read_pdb_line(line)
                if keep_membrane:
                    if resname in LIPID_RESIDUES:
                        last_anumb += 1
                        new_file_body += new_pdb_line(last_anumb, aname, resname,
                                                      resnumb, x, y, z, chain=chain)

                    if resname in list(Config.pypka_params.LIPIDS.values()):
                        aname, resname, to_include = convert_FF_atomnames(aname, resname)
                        if to_include:
                            last_anumb += 1
                            resnumb += last_resnumb
                            new_file_body += new_pdb_line(last_anumb, aname, resname,
                                                          resnumb, x, y, z, chain=chain)
                if keep_ions and aname in IONS and resname == aname:
                    last_anumb += 1
                    resnumb += last_resnumb
                    new_file_body += new_pdb_line(last_anumb, aname, resname,
                                                  resnumb, x, y, z, chain=chain)

    with open(add_to_pdb, 'a') as f_new:
        f_new.write(new_file_body)

def convert_FF_atomnames(aname, resname):
    popc_resname = Config.pypka_params.LIPIDS['POPC']
    chol_resname = Config.pypka_params.LIPIDS['cholesterol']

    lookup = [((' N  ', popc_resname), ('NTM', 'CHL')),
              ((' C12', popc_resname), ('CB ', 'CHL')),
              ((' C13', popc_resname), ('CN1', 'CHL')),
              ((' C14', popc_resname), ('CN2', 'CHL')),
              ((' C15', popc_resname), ('CN3', 'CHL')),
              ((' C11', popc_resname), ('CA ', 'PJ2')),
              ((' P  ', popc_resname), ('P  ', 'PJ2')),
              ((' O13', popc_resname), ('OG ', 'PJ2')),
              ((' O14', popc_resname), ('OB ', 'PJ2')),
              ((' O12', popc_resname), ('OA ', 'PJ2')),
              ((' O11', popc_resname), ('OD ', 'PJ2')),
              ((' C1 ', popc_resname), ('CD ', 'POX')),
              ((' C2 ', popc_resname), ('CE ', 'POX')),
              ((' O21', popc_resname), ('OE ', 'POX')),
              ((' C21', popc_resname), ('C1A', 'POX')),
              ((' O22', popc_resname), ('O1A', 'POX')),
              ((' C22', popc_resname), ('C1B', 'POX')),
              ((' C3 ', popc_resname), ('CZ ', 'POX')),
              ((' O31', popc_resname), ('OZ ', 'POX')),
              ((' C31', popc_resname), ('C2A', 'POX')),
              ((' O32', popc_resname), ('O2A', 'POX')),
              ((' C32', popc_resname), ('C2B', 'POX')),
              ((' C23', popc_resname), ('C1C', 'POX')),
              ((' C24', popc_resname), ('C1D', 'POX')),
              ((' C25', popc_resname), ('C1E', 'POX')),
              ((' C26', popc_resname), ('C1F', 'POX')),
              ((' C27', popc_resname), ('C1G', 'POX')),
              ((' C28', popc_resname), ('C1H', 'POX')),
              ((' C29', popc_resname), ('C1I', 'POX')),
              (('0C21', popc_resname), ('C1J', 'POX')),
              (('1C21', popc_resname), ('C1K', 'POX')),
              (('2C21', popc_resname), ('C1L', 'POX')),
              (('3C21', popc_resname), ('C1M', 'POX')),
              (('4C21', popc_resname), ('C1N', 'POX')),
              (('5C21', popc_resname), ('C1O', 'POX')),
              (('6C21', popc_resname), ('C1P', 'POX')),
              (('7C21', popc_resname), ('C1Q', 'POX')),
              (('8C21', popc_resname), ('C1R', 'POX')),
              ((' C33', popc_resname), ('C2C', 'POX')),
              ((' C34', popc_resname), ('C2D', 'POX')),
              ((' C35', popc_resname), ('C2E', 'POX')),
              ((' C36', popc_resname), ('C2F', 'POX')),
              ((' C37', popc_resname), ('C2G', 'POX')),
              ((' C38', popc_resname), ('C2H', 'POX')),
              ((' C39', popc_resname), ('C2I', 'POX')),
              (('0C31', popc_resname), ('C2J', 'POX')),
              (('1C31', popc_resname), ('C2K', 'POX')),
              (('2C31', popc_resname), ('C2L', 'POX')),
              (('3C31', popc_resname), ('C2M', 'POX')),
              (('4C31', popc_resname), ('C2N', 'POX')),
              (('5C31', popc_resname), ('C2O', 'POX')),
              (('6C31', popc_resname), ('C2P', 'POX')),
              ((" H3", chol_resname), (" H3", 'CHO')),
              ((" C3", chol_resname), (" C3", 'CHO')),
              ((" C2", chol_resname), (" C2", 'CHO')),
              ((" C1", chol_resname), (" C1", 'CHO')),
              (("C10", chol_resname), ("C10", 'CHO')),
              (("C19", chol_resname), ("C19", 'CHO')),
              ((" C9", chol_resname), (" C9", 'CHO')),
              (("C11", chol_resname), ("C11", 'CHO')),
              (("C12", chol_resname), ("C12", 'CHO')),
              (("C13", chol_resname), ("C13", 'CHO')),
              (("C18", chol_resname), ("C18", 'CHO')),
              (("C17", chol_resname), ("C17", 'CHO')),
              (("C20", chol_resname), ("C20", 'CHO')),
              (("C22", chol_resname), ("C22", 'CHO')),
              (("C23", chol_resname), ("C23", 'CHO')),
              (("C24", chol_resname), ("C24", 'CHO')),
              (("C25", chol_resname), ("C25", 'CHO')),
              (("C27", chol_resname), ("C27", 'CHO')),
              (("C26", chol_resname), ("C26", 'CHO')),
              (("C21", chol_resname), ("C21", 'CHO')),
              (("C16", chol_resname), ("C16", 'CHO')),
              (("C15", chol_resname), ("C15", 'CHO')),
              (("C14", chol_resname), ("C14", 'CHO')),
              ((" C8", chol_resname), (" C8", 'CHO')),
              ((" C7", chol_resname), (" C7", 'CHO')),
              ((" C6", chol_resname), (" C6", 'CHO')),
              ((" H6", chol_resname), (" H6", 'CHO')),
              ((" C5", chol_resname), (" C5", 'CHO')),
              ((" C4", chol_resname), (" C4", 'CHO')),
              ((" O3", chol_resname), (" O3", 'CHO')),
              (("H3'", chol_resname), ("H3'", 'CHO'))]

    for swap in lookup:
        old_aname = swap[0][0].strip()
        old_resname = swap[0][1]
        new_aname = swap[1][0].strip()
        new_resname = swap[1][1]
        if old_aname == aname and \
           old_resname == resname:
            aname = new_aname
            resname = new_resname
            return aname, resname, True

    return aname, resname, False
