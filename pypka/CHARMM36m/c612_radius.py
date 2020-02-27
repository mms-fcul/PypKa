rt = 'min' # 2 5 'zero' 'min'

water_atype = 'OW'

ignore_list = ()

rules = {
    'NTR': {
        'N': 'nh3',
        'HT1': 'hc',
        'HT2': 'hc',
        'HT3': 'hc',
        'CA': 'ct1',
        'C': 'c',
        'O': 'o',
        'HA': 'ha'
    },
    'CTR': {
        'CT': 'ca',
        'OT1': 'o',
        'OT2': 'o',
        'HT11': 'h',
        'HT12': 'h',
        'HT21': 'h',
        'HT22': 'h'
    }
}

rules['NT0'] = rules['NTR']
rules['NT1'] = rules['NTR']
rules['NT2'] = rules['NTR']
rules['NT3'] = rules['NTR']
rules['CT0'] = rules['CTR']
rules['CT1'] = rules['CTR']
rules['CT2'] = rules['CTR']
rules['CT3'] = rules['CTR']
rules['CT4'] = rules['CTR']


def radius_from_lennard(sigma, epsilon, rt='min', water=False):
    def calc_radius(rmin, c6, c12, rt):
        if isinstance(rt, int):
            kT = 2.49432
            radius = (rmin ** (-6) + (rt * kT / c12) ** 0.5) ** (-1/6)
        elif rt == 'zero':
            radius = rmin * 2 ** (-1/6)
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

    #sigma, epsilon = water_interaction(sigma, epsilon)

    c6  = 4 * epsilon * (sigma ** 6)
    c12 = 4 * epsilon * (sigma ** 12)

    rmin = (2 * c12 / c6) ** (1.0 / 6)

    radius = calc_radius(rmin, c6, c12, rt)

    if water:
        return 10 * radius / 2
    else:
        return 10 * radius / 2 #10 * radius - water_radius

db_lj = {}
db_lj_atomic = {}
with open('ffnonbonded.itp') as f:
    trigger = False
    for line in f:
        if ';' in line:
            comment_i = line.index(';')
            line = line[:comment_i]
        cols = line.strip().split()
        if not line.startswith(';') and "atomtypes" in line:
            trigger = True
        elif "[ " in line:
            trigger = False

        if trigger and len(cols) == 7:
            atype = cols[0].lower()
            sigma, epsilon = float(cols[-2]), float(cols[-1])
            db_lj[atype] = (sigma, epsilon)

db_res = {}
trigger = False
with open('merged.rtp') as f:
    for line in f:
        if line.startswith('['):      # begin residue block
            resname = line.split(';')[0].strip('[ ]\n')
            db_res[resname] = {}
        elif '[ atoms ]' in line:
            trigger = True
        elif ' [ ' in line:
            trigger = False
        elif trigger and line.strip():
            cols = line.strip().split()
            aname, atype = cols[0], cols[1].lower()
            db_res[resname][aname] = atype

for i in range(4):
    rules[f'SE{i}'] = db_res['SER']
    for ii in range(4):
        rules[f'SE{i}'][f'HG{ii}'] = db_res['SER']['HG1']
for i in range(4):
    rules[f'TH{i}'] = db_res['THR']
    for ii in range(4):
        rules[f'TH{i}'][f'HG{ii}'] = db_res['THR']['HG1']
rules['HI2'] = db_res['HSP']
rules['HIS'] = db_res['HI2']
rules['HI0'] = db_res['HI2']
rules['HI1'] = db_res['HI2']
rules['LYS'] = db_res['LY3']

db_res['LYS'] = db_res['LY3']
db_res['LY0'] = db_res['LY3']
db_res['LY1'] = db_res['LY3']
db_res['LY2'] = db_res['LY3']

db_res['GLU'] = db_res['GL4']
db_res['GL0'] = db_res['GL4']
db_res['GL1'] = db_res['GL4']
db_res['GL2'] = db_res['GL4']
db_res['GL3'] = db_res['GL4']

db_res['ASP'] = db_res['AS4']
db_res['AS0'] = db_res['AS4']
db_res['AS1'] = db_res['AS4']
db_res['AS2'] = db_res['AS4']
db_res['AS3'] = db_res['AS4']

sigma, epsilon = db_lj[water_atype.lower()]
water_radius = radius_from_lennard(sigma, epsilon, rt=rt, water=True)

##sigma_OW = 3.16557e-01
##epsilon_OW = 6.50194e-01
##c6_OW = 0.002617064809528995
##c12_OW = 2.633455636808028e-06
##radius = radius_from_lennard(c6_OW, c12_OW, rt='min', water=False)
#print(water_radius)
#exit()

#water_radius = 1.4241233005024403
print("""!siz file modified by VT 2013/03/15
!Hydrogen atoms have 0.001 radii because
! DelPhi complains if H has charge and no
! radii. It does not change final results
atom__res_radius_""")
with open('DataBaseT_old.siz') as f:
    for line in f.readlines()[5:]:
        cols = line.strip().split()
        aname, resname = cols[0:2]

        if resname in ignore_list:
            continue

        if resname in db_res and \
           aname in db_res[resname]:
            atype = db_res[resname][aname]
        elif resname in rules and \
             aname in rules[resname]:
            atype = rules[resname][aname]
        else:
            print(aname, resname, 'ignored')
            exit()
            continue

        sigma, epsilon = db_lj[atype]

        radius = radius_from_lennard(sigma, epsilon, rt=rt)

        out = f'{aname:4}  {resname:3}   {radius:1.3f}'

        print(out)
