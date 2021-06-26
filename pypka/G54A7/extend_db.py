### USER INPUT ##

rt = 2  # Options: 2 5 'zero' 'min'

out_siz = "tmp.siz"
out_crg = "tmp.crg"

###

water_atype = "OW"

f_nonbonded = "ffG54a7pHtnb.itp"
f_rtp = "ffG54a7pHt.rtp"
old_siz_db = "DataBaseT_template.siz"
old_crg_db = "DataBaseT_template.crg"  # NTR and CTR charges are copied from here

ignore_list = ()

### MAIN ###

# when adding new titratable residues please include them here
sizes = {
    "HI2": ("HI0", "HI1"),
    "AS4": ("AS0", "AS1", "AS2", "AS3"),
    "LY3": ("LY0", "LY1", "LY2"),
    "GL4": ("GL0", "GL1", "GL2", "GL3"),
    "TY2": ("TY0", "TY1"),
    "CY2": ("CY0", "CY1"),
    "GD4": ("GD0", "GD1", "GD2", "GD3"),
    "PA4": ("PA0", "PA1", "PA2", "PA3"),
    "AR4": ("AR0", "AR1", "AR2", "AR3"),
    "PJ2": ("PJ0", "PJ1"),
}

rules = {
    "NTR": {
        "N": "nl",
        "H1": "h",
        "H2": "h",
        "H3": "h",
        "CA": "ch1",
        "C": "c",
        "O": "o",
    },
    "CTR": {
        "C": "c",
        "O1": "om",
        "O2": "om",
        "HO11": "h",
        "HO12": "h",
        "HO21": "h",
        "HO22": "h",
    },
}

rules["NT0"] = rules["NTR"]
rules["NT1"] = rules["NTR"]
rules["NT2"] = rules["NTR"]
rules["NT3"] = rules["NTR"]
rules["CT0"] = rules["CTR"]
rules["CT1"] = rules["CTR"]
rules["CT2"] = rules["CTR"]
rules["CT3"] = rules["CTR"]
rules["CT4"] = rules["CTR"]


def radius_from_lennard(c6, c12, rt="min", water=False):
    def calc_radius(rmin, c6, c12, rt):
        if isinstance(rt, int):
            kT = 2.49432
            radius = (rmin ** (-6) + (rt * kT / c12) ** 0.5) ** (-1 / 6)
        elif rt == "zero":
            radius = rmin * 2 ** (-1 / 6)
        else:
            radius = rmin
        return radius

    if not c6:
        return 0.001

    rmin = (2 * c12 / c6) ** (1.0 / 6)

    radius = calc_radius(rmin, c6, c12, rt)

    if water:
        return 10 * radius / 2
    else:
        return 10 * radius - water_radius


db_lj = {}
db_lj_atomic = {}
with open(f_nonbonded) as f:
    trigger = False
    for line in f:
        if ";" in line:
            comment_i = line.index(";")
            line = line[:comment_i]
        cols = line.strip().split()
        if not line.startswith(";") and "nonbond_params" in line:
            trigger = True
        elif "[ " in line:
            trigger = False

        if trigger and len(cols) == 5:
            if cols[0] == water_atype:
                atype = cols[1].lower()
            elif cols[1] == water_atype:
                atype = cols[0].lower()
            else:
                continue

            c6, c12 = float(cols[-2]), float(cols[-1])
            db_lj[atype] = (c6, c12)

db_res = {}
db_charges = {}
trigger = False
with open(f_rtp) as f:
    for line in f:
        if line.startswith(";"):
            continue

        if "[ atoms ]" in line:
            trigger = True
            continue

        elif "[" in line:
            trigger = False

            resname = line.split(";")[0].strip("[ ]\n")

            if resname not in (
                "angles",
                "exclusions",
                "dihedrals",
                "impropers",
                "bonds",
            ):
                db_res[resname] = {}
                db_charges[resname] = {}

        if trigger and line.strip():
            cols = line.strip().split()
            aname, atype, charge = cols[0], cols[1].lower(), cols[2]
            db_res[resname][aname] = atype
            db_charges[resname][aname] = float(charge)

rules_crg = {}

rules["HIS"] = db_res["HI2"]
rules["ASP"] = db_res["AS4"]
rules["LYS"] = db_res["LY3"]
rules["GLU"] = db_res["GL4"]
rules["TYR"] = db_res["TY2"]
rules["CYS"] = db_res["CYS2"]
rules["ARG"] = db_res["AR4"]

for ref, taus in sizes.items():
    for tau in taus:
        db_res[tau] = db_res[ref]

# rules_crg["NTR"] = db_charges["NT3"]
# rules_crg["CTR"] = db_charges["CT4"]
rules_crg["HIS"] = db_charges["HI2"]
rules_crg["ASP"] = db_charges["AS4"]
rules_crg["LYS"] = db_charges["LY3"]
rules_crg["GLU"] = db_charges["GL4"]
rules_crg["TYR"] = db_charges["TY2"]
rules_crg["CYS"] = db_charges["CYS2"]
rules_crg["ARG"] = db_charges["AR4"]

# c6, c12 = db_lj[water_atype.lower()]
# water_radius = radius_from_lennard(c6, c12, rt=rt, water=True)
water_radius = 1.4241233005024403
new_siz = f"""! aminoacids from 54a7
! created by extend_db.py
! water_atype = {water_atype}; rt = {rt}
atom__res_radius_\n"""
with open(old_siz_db) as f:
    for line in f:
        if line.startswith("!") or "atom__res_radius_" in line:
            continue
        cols = line.strip().split()
        aname, resname = cols[0:2]

        if resname in ignore_list:
            continue

        if resname in db_res and aname in db_res[resname]:
            atype = db_res[resname][aname]
        elif resname in rules and aname in rules[resname]:
            atype = rules[resname][aname]
        else:
            print(aname, resname, "ignored on SIZ")
            continue

        sigma, epsilon = db_lj[atype]

        radius = radius_from_lennard(sigma, epsilon, rt=rt)
        out = f"{aname:<4}  {resname:3} {radius:6.3f}\n"

        new_siz += out

with open(out_siz, "w") as f:
    f.write(new_siz)

new_crg = """! aminoacids from 54a7
! created by extend_db.py
!
atom__resnumbc_charge_\n"""
with open(old_crg_db) as f:
    for line in f:
        if line.startswith("!") or "atom__resnumbc_charge_" in line:
            continue

        cols = line.strip().split()
        aname, resname = cols[0:2]

        if resname in ignore_list:
            continue

        if resname in db_charges and aname in db_charges[resname]:
            charge = db_charges[resname][aname]
            if cols[2] != "X.XXX" and float(cols[2]) != charge:
                print("ERROR CRG", aname, resname, charge, cols)
        elif resname[:2] in ("CT", "NT"):
            charge = float(cols[2])
        elif resname in rules_crg and aname in rules_crg[resname]:
            charge = rules_crg[resname][aname]
        else:
            print(aname, resname, "ignored on CRG")
            continue

        out = f"{aname:<4}  {resname:8} {charge:6.3f}\n"

        new_crg += out

with open(out_crg, "w") as f:
    f.write(new_crg)
