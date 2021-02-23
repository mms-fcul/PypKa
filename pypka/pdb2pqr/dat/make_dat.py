fname_crg = "../../G54A7/DataBaseT.crg"
fname_siz = "../../G54A7/DataBaseT.siz"

with open(fname_crg) as f_crg, open(fname_siz) as f_siz:
    for l_crg, l_siz in zip(f_crg, f_siz):
        if l_crg.startswith("!") or "atom_" in l_crg:
            continue
        l_crg, l_siz = l_crg.strip(), l_siz.strip()

        aname, resname, charge = l_crg.split()
        aname, resname, radii = l_siz.split()

        line = f"{resname:3}     {aname:3}      {charge:8}       {radii:6}  {aname:3}"

        print(line)
