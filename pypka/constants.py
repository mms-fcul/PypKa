""""Constants file
"""

# Constants
KBOLTZ = 5.98435e-6  # e^2/(Angstrom*K)
LOG10 = 2.302585092994046

# Lipids
LIPIDS = {"cholesterol": "CHO", "POPC": "POP"}  # to edit
LIPID_RESIDUES = ["DMX", "POX", "PJ2", "PJ1", "PJ0", "CHL"]  # allowed residue names

DNA_RESIDUES = {
    "DA": "DAD",
    "DC": "DCY",
    "DT": "DTH",
    "DG": "DGU",
}

PDB_RNA_RESIDUES = {
    "A": "RA",
    "C": "RC",
    "U": "RU",
    "G": "RG",
}

RNA_RESIDUES = {
    "RA": "RAD",
    "RC": "RCY",
    "RU": "RUR",
    "RG": "RGU",
}

NUCLEIC_ACIDS = list(DNA_RESIDUES.values()) + list(RNA_RESIDUES.values())

# Ions
IONS = ["CL-", "NA+"]

# Tautomers and residues
TERMINAL_OFFSET = 5000

TITRABLETAUTOMERS = {
    "LYS": 3,
    "HIS": 2,
    "ASP": 4,
    "GLU": 4,
    "SER": 3,
    "THR": 3,
    "CYS": 3,
    "CTR": 4,
    "NTR": 3,
    "TYR": 2,
}

PROTEIN_RESIDUES = (
    "ALA",
    "CYS",
    "ASP",
    "GLU",
    "PHE",
    "GLY",
    "HIS",
    "ILE",
    "LYS",
    "LEU",
    "MET",
    "ASN",
    "PRO",
    "GLN",
    "ARG",
    "SER",
    "THR",
    "VAL",
    "TRP",
    "TYR",
    "NTR",
    "CTR",
)


TITRABLERESIDUES = list(TITRABLETAUTOMERS.keys())
REGULARTITRATINGRES = list(TITRABLETAUTOMERS.keys())

for res in REGULARTITRATINGRES:
    ntautomers = TITRABLETAUTOMERS[res]
    for i in range(ntautomers + 1):
        TITRABLERESIDUES.append(res[0:2] + str(i))

MANDATORY_PARAMS = [
    "structure",
    "epsin",
    "ionicstr",
    "pbc_dimensions",
    "ncpus",
    "sites",
]

IGNORED_PARAMS = [
    "tmpsites",
    "pid",
    "file_dir",
    "pdb2pqr",
    "userff",
    "usernames",
    "delphi_params",
    "mc_params",
    "pH",
    "lipid_definition",
    "f_dat",
]

MAXNPKHALFS = 3
PKAPLACEHOLDER = None
