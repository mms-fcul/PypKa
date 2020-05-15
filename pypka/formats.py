from config import Config
from constants import TERMINAL_OFFSET, TITRABLETAUTOMERS, PROTEIN_RESIDUES

def new_pdb_line(aID, aname, resname, resnumb, x, y, z, chain=' '):
    pdb_format = "ATOM  {:5d} {:4s} {:4s}{:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}\n"
    return pdb_format.format(aID, aname, resname, chain, resnumb, x, y, z)

def new_pqr_line(aID, aname, resname, resnumb, x, y, z, charge, radius, chain=' '):
    pdb_format = "ATOM  {:5d} {:4s} {:4s}{:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}{:8.3f}{:8.3f}\n"
    return pdb_format.format(aID, aname, resname, chain, resnumb, x, y, z, charge, radius)

def new_gro_line(aID, aname, resname, resnumb, x, y, z):
    gro_format = '{:5}{:5}{:>5}{:5}{:8.3f}{:8.3f}{:8.3f}\n'
    return gro_format.format(resnumb, resname, aname, aID, x, y, z)


def correct_longHs(aname):
    if aname[1] == 'H' and aname[0] in '12' \
       and aname[3] in '12':
        return aname[1:3] + aname[3] + aname[0]
    return aname


def read_pqr_line(line):
    aname   = line[12:16].strip()
    anumb   = int(line[7:11].strip())
    resname = line[17:21].strip()
    chain   = line[21]
    resnumb = int(line[22:26])
    x       = float(line[30:38])
    y       = float(line[38:46])
    z       = float(line[46:54])
    charge  = float(line[54:62])
    radius  = float(line[62:70])
    if len(aname) == 4:
        aname = correct_longHs(aname)
    return (aname, anumb, resname, chain, resnumb, x, y, z, charge,
            radius)

def read_pdb_line(line):
    aname   = line[12:16].strip()
    anumb   = int(line[5:11].strip())
    resname = line[17:21].strip()
    chain   = line[21]
    resnumb = int(line[22:26])
    x       = float(line[30:38])
    y       = float(line[38:46])
    z       = float(line[46:54])
    if len(aname) == 4:
        aname = correct_longHs(aname)
    return (aname, anumb, resname, chain, resnumb, x, y, z)

def read_gro_line(line):
    resnumb = int(line[:5].strip())
    resname = line[5:10].strip()
    aname   = line[10:15].strip()
    anumb   = int(line[15:20].strip())
    x = float(line[20:28].strip())
    y = float(line[28:36].strip())
    z = float(line[36:44].strip())
    return (aname, anumb, resname, resnumb, x, y, z)

def gro2pdb(f_in, f_out, save_box=False):
    new_pdb_content = "REMARK    CONVERTED from {} by PypKa\n".format(f_in)
    with open(f_in) as f:
        line_counter = 0
        natoms_left = 0
        natoms = 0
        lines = f.read().splitlines()
        penul_line = lines[-2]
        last_line  = lines[-1]
        box = [float(i) * 10 for i in last_line.split()[:3]]
        if save_box:
            Config.pypka_params.setBox(box)
        new_pdb_content += "CRYST1  {0:3.3f}  {1:3.3f}  {2:3.3f}  "\
                           "60.00  60.00  90.00 P 1           1\n".format(box[0],
                                                                          box[1], box[2])
        for line in lines:
            line_counter += 1
            if natoms_left > 0:
                natoms_left -= 1
                natoms += 1
                (aname, anumb, resname,
                 resnumb, x, y, z) = read_gro_line(line)

                x, y, z = x * 10, y * 10, z * 10

                if resname[-1] == 'X' and Config.pypka_params.CpHMD_mode:
                    resname = resname[:-1] + '0'

                new_pdb_content += new_pdb_line(anumb, aname,
                                                resname, resnumb,
                                                x, y, z, chain='A')
            elif line_counter == 2:
                natoms_left = int(line.strip())
    new_pdb_content += 'TER\nENDMDL\n'
    with open(f_out, 'w') as f_new:
        f_new.write(new_pdb_content)

def pdb2gro(filename_in, filename_out, chains_res, box=[], pqr=False,
            renumber_res=False, fix_termini=True):
    """
    Returns
      - aposition (int) of last atom id in filename_out
    """
    NTR_atoms = Config.pypka_params['NTR_atoms']
    CTR_atoms = Config.pypka_params['CTR_atoms']

    header = 'CREATED within PyPka\n'
    new_pdb_text = ''
    aposition = 0
    rposition = 1
    prev_resnumb = None

    with open(filename_in) as f:
        for line in f:
            if 'CRYST1' in line[:6]:
                parts = line.split()
                pdb_box = [float(parts[1]), float(parts[2]), float(parts[3])]
                continue
            elif 'ATOM' != line[:4]:
                continue
            elif pqr:
                (aname, anumb, resname, chain, resnumb, x, y,
                 z, charge, radius) = read_pqr_line(line)
            else:
                (aname, anumb, resname, chain,
                 resnumb, x, y, z) = read_pdb_line(line)
            aposition += 1

            chains_sites_numbs = []
            if chain in chains_res:
                chains_sites_numbs = chains_res[chain]
                for site_numb, site_name in chains_sites_numbs.items():
                    if site_name == 'NTR':
                        NTR_numb = int(site_numb)
                    elif site_name == 'CTR':
                        CTR_numb = int(site_numb)

            if fix_termini and resnumb in chains_sites_numbs:
                if (resnumb == NTR_numb and aname in NTR_atoms) or \
                   (resnumb == CTR_numb and aname in CTR_atoms):
                    site_numb = resnumb
                    resname = chains_res[chain][site_numb]

            elif (resname == 'HIS' and aname == 'HD1' and
                  resnumb not in chains_sites_numbs):
                aposition -= 1
                continue

            if prev_resnumb and prev_resnumb != resnumb:
                rposition += 1

            x /= 10.0
            y /= 10.0
            z /= 10.0
            if not renumber_res:
                rposition = resnumb
            new_pdb_text += new_gro_line(aposition, aname, resname,
                                         rposition, x, y, z)
            prev_resnumb = resnumb

    header += '{0}\n'.format(aposition)
    if box == []:
        footer = '{0:10.5f}{1:10.5f}{2:10.5f}\n'.format(pdb_box[0] / 10.0,
                                                        pdb_box[1] / 10.0,
                                                        pdb_box[2] / 10.0)
    else:
        footer = '{0:10.5f}{1:10.5f}{2:10.5f}\n'.format(box[0] / 10.0,
                                                        box[1] / 10.0,
                                                        box[2] / 10.0)

    new_pdb = header + new_pdb_text + footer
    with open(filename_out, 'w') as f_new:
        f_new.write(new_pdb)

    return aposition


def correct_names(resnumb, resname, aname,
                  titrating_sites, NTR_numb, CTR_numb):
    def change_aname(aname, restype, mode='regular'):
        if mode == 'titrating':
            not_correct_names = list(correct_atoms_sites_table[restype].keys())
        else:
            not_correct_names = list(correct_atoms_table[restype].keys())
        for not_corrected in not_correct_names:
            if aname == not_corrected:
                if mode == 'titrating':
                    aname = correct_atoms_sites_table[restype][not_corrected]
                else:
                    aname = correct_atoms_table[restype][not_corrected]

        return aname

    def change_resname(resname):
        for tit_res in TITRABLETAUTOMERS.keys():
            if tit_res[:2] == resname[:2]:
                return tit_res
        return resname


    # no longer used as it is done by pdb2pqr
    correct_atoms_table = {'CTR': {'O': 'O1',
                                   'OXT': 'O2'},
                           'NTR': {'H': 'H1'},
                           'ILE': {'CD1': 'CD'}}

    correct_atoms_sites_table = {'CYS': {'HG': 'HG1'},
                                 'SER': {'HG': 'HG1'},
                                 'TYR': {'HH': 'HH1'}}

    correct_residues_table = {'HSD': 'HI0',
                              'HSE': 'HI1',
                              'HSP': 'HI2',
                              'ARGN': 'AR0',
                              'ASPH': 'AS0',
                              'CYS2': 'CYS',
                              'CYSH': 'CY0',
                              'GLUH': 'GL0',
                              'HISD': 'HI0',
                              'HISE': 'HI1',
                              'HISH': 'HI2',
                              'HISA': 'HI0',
                              'HISB': 'HI1',
                              'HISP': 'HI2',
                              'LYSH': 'LY3',
                              'LYSN': 'LY0',
                              'ISU': 'CYS'}

    if resnumb == CTR_numb:
        restype = 'CTR'
        aname = change_aname(aname, restype)

    if resnumb == NTR_numb:
        restype = 'NTR'
        aname = change_aname(aname, restype)

    if resname in list(correct_atoms_table.keys()):
        aname = change_aname(aname, resname)

    if resname in list(correct_residues_table.keys()):
        resname = correct_residues_table[resname]

    if resnumb in titrating_sites:
        if resname not in PROTEIN_RESIDUES:
            resname = change_resname(resname)
        if resname in list(correct_atoms_sites_table.keys()):
            aname = change_aname(aname, resname, mode='titrating')


    return aname, resname


def readPBPFile(f_dat):
    with open(f_dat) as f:
        for line in f:
            if line[0] != "#" and len(line) > 1:
                parts = line.strip().split('=')
                param = parts[0].split()[-1]
                value = parts[1].replace('"', '').replace("'", '')
                config.params[param] = value.strip()


def convertTermini(site_numb):
    if site_numb >= TERMINAL_OFFSET:
        return site_numb - TERMINAL_OFFSET
    return site_numb
