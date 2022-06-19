### USER INPUT ##

rt = 2  # Options: 2 5 'zero' 'min'

out_siz = "tmp.siz"
out_crg = "tmp.crg"

###

atomistic = False
water_atype = "OW"

f_nonbonded = "ffnonbonded.itp"
f_rtp = "aa.rtp"
template_siz_db = "DataBaseT_template.siz"
template_crg_db = "DataBaseT_template.crg"  # NTR and CTR charges are copied from here

ignore_list = ()

### MAIN ###
rules = {
    "NTR": {
        "N": "nh3",
        "H1": "hc",
        "H2": "hc",
        "H3": "hc",
        "CA": "ct1",
        "C": "c",
        "O": "o",
        "HA": "ha",
    },
    "CTR": {
        "C": "ca",
        "OT1": "o",
        "OT2": "o",
        "HT11": "h",
        "HT12": "h",
        "HT21": "h",
        "HT22": "h",
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


def radius_from_lennard(sigma, epsilon, rt="min", water=False, atomistic=False):
    def calc_radius(rmin, c6, c12, rt):
        if isinstance(rt, int):
            kT = 2.49432
            radius = (rmin ** (-6) + (rt * kT / c12) ** 0.5) ** (-1 / 6)
        elif rt == "zero":
            radius = rmin * 2 ** (-1 / 6)
        else:
            radius = rmin
        return radius

    def water_interaction(sigma, epsilon):
        sigma_w, epsilon_w = db_lj[water_atype.lower()]
        sigma_i = (sigma_w + sigma) / 2
        epsilon_i = (epsilon * epsilon_w) ** 0.5

        return sigma_i, epsilon_i

    if not epsilon or not sigma:
        return 0.001

    if not atomistic:
        sigma, epsilon = water_interaction(sigma, epsilon)

    c6 = 4 * epsilon * (sigma ** 6)
    c12 = 4 * epsilon * (sigma ** 12)

    rmin = (2 * c12 / c6) ** (1.0 / 6)

    radius = calc_radius(rmin, c6, c12, rt)

    if water or atomistic:
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
        if not line.startswith(";") and "atomtypes" in line:
            trigger = True
        elif "[ " in line:
            trigger = False

        if trigger and len(cols) == 7:
            atype = cols[0].lower()
            sigma, epsilon = float(cols[-2]), float(cols[-1])
            db_lj[atype] = (sigma, epsilon)

db_res = {}
db_charges = {}
trigger = False
with open(f_rtp) as f:
    for line in f:
        if line.startswith("["):  # begin residue block
            resname = line.split(";")[0].strip("[ ]\n")
            db_res[resname] = {}
            db_charges[resname] = {}
        elif "[ atoms ]" in line:
            trigger = True
        elif " [ " in line:
            trigger = False
        elif trigger and line.strip():
            cols = line.strip().split()
            aname, atype, charge = cols[0], cols[1].lower(), cols[2]
            db_res[resname][aname] = atype
            db_charges[resname][aname] = float(charge)

for i in range(4):
    rules[f"SE{i}"] = db_res["SER"]
    for ii in range(4):
        rules[f"SE{i}"][f"HG{ii}"] = db_res["SER"]["HG1"]
for i in range(4):
    rules[f"TH{i}"] = db_res["THR"]
    for ii in range(4):
        rules[f"TH{i}"][f"HG{ii}"] = db_res["THR"]["HG1"]

rules["HIS"] = db_res["HI2"]
rules["ASP"] = db_res["AS4"]
rules["LYS"] = db_res["LY3"]
rules["GLU"] = db_res["GL4"]
rules["TYR"] = db_res["TY2"]
rules["CYS"] = db_res["CY3"]
rules["SSB"] = db_res["CYS2"]

rules_crg = {}
rules_crg["HIS"] = db_charges["HI2"]
rules_crg["ASP"] = db_charges["AS4"]
rules_crg["LYS"] = db_charges["LY3"]
rules_crg["GLU"] = db_charges["GL4"]
rules_crg["TYR"] = db_charges["TY2"]
rules_crg["CYS"] = db_charges["CY3"]
rules_crg["SSB"] = db_charges["CYS2"]

sigma, epsilon = db_lj[water_atype.lower()]
water_radius = radius_from_lennard(sigma, epsilon, rt=rt, water=True)

new_siz = """!siz file
!Hydrogen atoms have 0.001 radii because
! DelPhi complains if H has charge and no
! radii. It does not change final results
atom__res_radius_\n"""
with open(template_siz_db) as f:
    for line in f.readlines()[5:]:
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

        radius = radius_from_lennard(sigma, epsilon, rt=rt, atomistic=atomistic)
        out = f"{aname:>4}  {resname:3} {radius:6.3f}\n"

        new_siz += out

with open(out_siz, "w") as f:
    f.write(new_siz)

new_crg = """!crg file
!aminoacids from charmm36m
atom__resnumbc_charge_\n"""
with open(template_crg_db) as f:
    for line in f.readlines()[3:]:
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
        out = f"{aname:>4}  {resname:8} {charge:6.3f}\n"

        new_crg += out

with open(out_crg, "w") as f:
    f.write(new_crg)
