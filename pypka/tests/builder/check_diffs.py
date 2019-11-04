import sys
import os
sys.path.insert(1, os.getcwd() +  '/../../')

from formats import read_pdb_line

def storeResidues(filename):
    residues = {}
    with open(filename) as f_original:
        for line in f_original:
            if line.startswith('ATOM '):
                (aname, anumb, resname,
                 chain, resnumb, x, y, z) = read_pdb_line(line)
                if resnumb not in residues:
                    residues[resnumb] = {}
                residues[resnumb][aname] = (resname, x, y, z)
    return residues

def roundXYZ(x, y, z):
    x = round(x, 2)
    y = round(y, 2)
    z = round(z, 2)
    return x, y, z

def areCloseEnough(p1, p2):
    x1, y1, z1 = p1
    x2, y2, z2 = p2
    return (abs(x1 - x2) <= 0.1 and 
            abs(y1 - y2) <= 0.1 and
            abs(z1 - z2) <= 0.1)

def dist3D(p1, p2):
    x1, y1, z1 = p1
    x2, y2, z2 = p2    
    
    return ((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2) ** 0.5


def compareFiles(f1, f2):
    original = storeResidues(f1)
    new = storeResidues(f2)

    assert(len(original) == len(new))

    problems = False
    
    for resid in original:
        res = original[resid]
        resnew = new[resid]
        len_orig = len(res)
        len_new = len(resnew)
        if len_orig != len_new:
            print(f'Error in res #{resid} {f1}: {len_orig} {f2}: {len_new}')
            problems = True
            for i in res:
                if i not in resnew:
                    print(f'{i} not in {f2}')
            for i in resnew:
                if i not in res:
                    print(f'{i} not in {f1}')

        for atom in res:
            if atom not in resnew:
                print(f'{atom} not in {res[atom][0]}{resid} {f2}')
                problems = True
                continue
            (resname1, x1, y1, z1) = res[atom]
            (resname2, x2, y2, z2) = resnew[atom]
            if (resname1 != resname2):
                print(f'Warning: resname {resname1} != {resname2}')
                problems = True
            p1 = roundXYZ(x1, y1, z1)
            p2 = (x2, y2, z2)
            if not areCloseEnough(p1, p2):
                WarningTrigger = True
                for atom2 in resnew:
                    (resname2, x2, y2, z2) = resnew[atom2]
                    p2 = (x2, y2, z2)
                    if areCloseEnough(p1, p2):
                        WarningTrigger = False
                        break
                if WarningTrigger:
                    print(f'{resnew[atom][0]}{resid} {atom}: {p1} != {resnew[atom][1:]} {dist3D(p1, p2)}')

        for atom in resnew:
            if atom not in res:
                print(f'{atom} not in {f2} {resnew[atom][0]}{resid}')
                problems = True
                continue
    return problems
