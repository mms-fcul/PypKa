import os
import log
from formats import (read_pdb_line, read_pqr_line, read_gro_line,
                     correct_names, new_pqr_line, new_gro_line, pdb2gro)
import config
from copy import copy

def inputPDBCheck(filename, sites):
    """
    chains_length, chains_res
    """
    if filename[-3:] in ('pdb', 'pqr'):
        filetype = 'pdb'
    elif filename[-3:] == 'gro':
        filetype = 'gro'
    else:
        raise Exception('Input file must be either a pdb or a gro.')

    chains_length = {}
    chains_res = {}
    done = {}

    for site in sites['A']:
        if site[-1] == 'C':
            resnumb = site[:-1]
            done[resnumb] = 'CTR'
        if site[-1] == 'N':
            resnumb = site[:-1]
            done[resnumb] = 'NTR'

    if filetype == 'pdb' and not config.params['clean_pdb']:
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
                    if not config.params['clean_pdb']:
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
                if chain == ' ':
                    chain = 'A'

                if chain_length == 1:
                    last_chain = chain

                if chain != last_chain and chain_length != 1:
                    chains_length[last_chain] = chain_length
                    chains_res[chain] = done
                    chain_length = 0
                    last_chain = chain

                if chain in sites and \
                   resnumb not in done and \
                   str(resnumb) in sites[chain]:
                    done[resnumb] = resname

    if filetype == 'pdb' and not config.params['clean_pdb']:
        new_gro_header += '{0}\n'.format(atom_number)
        with open('TMP.gro', 'w') as f:
            f.write(new_gro_header + new_gro_body + new_gro_footer)

    chains_length[last_chain] = chain_length
    chains_res[chain] = done

    return chains_length, chains_res


def cleanPDB(pdb_filename, pdb2pqr_path, chains_res,
             inputpqr, outputpqr, sites):
    """
    """

    inputpdbfile = removeMembrane(pdb_filename)

    logfile = 'LOG_pdb2pqr'
    errfile = 'LOG_pdb2pqr_err'

    log.redirectOutput("start", logfile)
    log.redirectErr("start", errfile)

    sites_numbs = list(sites.keys())

    # CTR O1/O2 will be deleted and a O/OXT will be added
    os.system('python2 {0} {1} {2} --ff {3} --ffout GROMOS '
              '--drop-water -v --chain'.format(pdb2pqr_path,
                                               inputpdbfile, inputpqr,
                                               config.params['ffinput']))

    if config.f_structure_out:
        os.system('python2 {0} {1} Hs.pqr --ff AMBER --ffout AMBER '
                  '--drop-water -v --chain'.format(pdb2pqr_path,
                                                   inputpdbfile))

    log.redirectOutput("stop", logfile)
    log.redirectErr("stop", errfile)

    if config.f_structure_out:
        with open('Hs.pqr') as f:
            for line in f:
                if line.startswith('ATOM '):
                    (aname, anumb, resname, chain, resnumb, x, y,
                     z, charge, radius) = read_pqr_line(line)

                    if resnumb in chains_res['A'] and \
                       resname in config.AMBER_Hs and \
                       aname in config.AMBER_Hs[resname] and \
                       aname not in config.AMBER_mainchain_Hs:
                        continue
                    elif aname[0] == 'H' and aname not in ('H1', 'H2', 'H3'):
                        if not resnumb in config.mainchain_Hs:
                            config.mainchain_Hs[resnumb] = []
                        config.mainchain_Hs[resnumb].append((aname, anumb, resname, chain, 
                                                             x, y, z))
    
    CYS_bridges = {'A': []}
    with open('LOG_pdb2pqr') as f:
        for line in f:
            if 'patched with CYX' in line:
                parts = line.split('patched')[0].replace('PATCH INFO: ', '').split()
                resname, chain, resnumb = parts
                if not chain in CYS_bridges:
                    CYS_bridges[chain] = []
                CYS_bridges[chain].append(int(parts[-1]))
    config.tit_mole.saveCYSBridges(CYS_bridges)
    
    new_pdb_text = ''
    removed_pdb_text = ''
    removed_pdb_lines = []
    #ntr_trigger = True
    with open(inputpqr) as f:
        for line in f:
            if "ATOM" in line[:4]:
                (aname, anumb, resname, chain, resnumb, x, y,
                 z, charge, radius) = read_pqr_line(line)
                
                if chain in (' ', 'A'):
                    aname, resname = correct_names(sites_numbs, resnumb,
                                                   resname, aname, sites_numbs)
                    resnumb_max = resnumb
                if line[26] != ' ':
                    resnumb += config.terminal_offset
                new_line = new_pqr_line(anumb, aname, resname,
                                        resnumb, x, y, z, charge, radius)
                if chain in (' ', 'A'):
                    new_pdb_text += new_line
                elif aname not in ('H1', 'H2', 'H3', 'O1', 'O2'):
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
    
    chains = copy(chains_res['A'])
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
    for res in chains_res['A']:
        if chains_res['A'][res] != 'NTR' and res not in CYS_bridges['A']:
            sites_addHtaut += '{0}_{1},'.format(res, chains_res['A'][res])

    if len(sites_addHtaut) > 0 and sites_addHtaut[-1] == ',':
        sites_addHtaut = sites_addHtaut[:-1]

    logfile = 'LOG_addHtaut'
    log.redirectErr("start", logfile)

    os.system(f'{config.script_dir}/addHtaut cleaned.pqr {sites_addHtaut} > {outputpqr}')

    log.redirectErr("stop", logfile)
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

    pdb2gro(outputpqr, "TMP.gro", box, sites, pqr=True, fix_termini=False)

    if config.params['pbc_dim'] == 2:
        addMembrane("TMP.gro", pdb_filename)

    if config.params['keep_ions']:
        addIons("TMP.gro", pdb_filename)

    tmpfiles = ('LOG_pdb2pqr', 'LOG_pdb2pqr_err', 'clean.pqr', 'cleaned.pqr',
                'cleaned_tau.pqr', 'input_clean.pdb', 'removed.pqr', 'Hs.pqr')
    for filename in tmpfiles:
        if not config.debug and os.path.isfile(filename):
            os.remove(filename)

def removeMembrane(pdbfile):
    nomembrane_text = ''
    with open(pdbfile) as f:
        for line in f:
            if 'ATOM ' == line[0:5]:
                (aname, anumb, resname,
                 chain, resnumb, x, y, z) = read_pdb_line(line)
                if resname not in config.lipid_residues:
                    nomembrane_text += line
            else:
                nomembrane_text += line
    with open('tmp.tmp', 'w') as f_new:
        f_new.write(nomembrane_text)
    os.rename('tmp.tmp', 'input_clean.pdb')
    return 'input_clean.pdb'


def addMembrane(grofile, pdbfile):
    atom_number = 0
    new_file_header = ''
    new_file_body = ''
    new_file_footer = ''
    # Read the grofile with only the protein
    with open(grofile) as f:
        nline = 0
        maxnlines = 0
        for line in f:
            nline += 1
            if nline > 2 and nline < maxnlines:
                atom_number += 1
                (aname, anumb, resname, resnumb, x, y, z) = read_gro_line(line)
                new_file_body += new_gro_line(atom_number, aname, resname,
                                              resnumb, x, y, z)
            elif nline == 2:
                natoms = int(line.strip())
                maxnlines = natoms + 3
            elif nline == 1:
                new_file_header += line
            else:
                new_file_footer += line

    # Read the original pdb with the membrane
    with open(pdbfile) as f:
        for line in f:
            if 'ATOM ' == line[0:5]:
                (aname, anumb, resname,
                 chain, resnumb, x, y, z) = read_pdb_line(line)
                if resname in config.lipid_residues:
                    atom_number += 1
                    x, y, z = x/10, y/10, z/10
                    new_file_body += new_gro_line(atom_number, aname, resname,
                                                  resnumb, x, y, z)

                if resname in list(config.lipids.values()):
                    aname, resname, to_include = convert_FF_atomnames(aname, resname)
                    if to_include:
                        atom_number += 1
                        x, y, z = x/10, y/10, z/10
                        new_file_body += new_gro_line(atom_number, aname, resname,
                                                      resnumb, x, y, z)

    with open('tmp.tmp', 'w') as f_new:
        new_file_header += str(atom_number) + '\n'
        f_new.write(new_file_header + new_file_body + new_file_footer)
    os.rename('tmp.tmp', grofile)


def addIons(grofile, pdbfile):
    with open(grofile) as f:
        content = f.readlines()
        header = content[0]
        natoms = int(content[1].strip())
        coordinates = content[2:-1]
        box = content[-1]

    (aname, last_anumb, resname, last_resnumb, x, y, z) = read_gro_line(coordinates[-1])
    extra_lines = ''
    with open(pdbfile) as f:
        for line in f:
            if line.startswith('ATOM'):
                (aname, anumb, resname,
                 chain, resnumb, x, y, z) = read_pdb_line(line)
                if aname in config.ions and resname == aname:
                    last_anumb += 1
                    last_resnumb +=1
                    natoms += 1
                    x, y, z = x/10, y/10, z/10
                    extra_lines += new_gro_line(last_anumb, aname, resname,
                                                last_resnumb, x, y, z)

    with open(grofile, 'w') as f_new:
        natoms = f'{natoms}\n'
        coordinates = ''.join(coordinates)
        content = header + natoms + coordinates + extra_lines + box
        f_new.write(content)

def convert_FF_atomnames(aname, resname):
    popc_resname = config.lipids['POPC']
    chol_resname = config.lipids['cholesterol']

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
