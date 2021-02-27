import sys

pKs = {}
with open(sys.argv[1]) as f:
    for line in f:
        line = line.strip()
        cols = line.split()

        resnumb, resname, pK = int(cols[0]), cols[1], float(cols[2])

        pKs[resname] = pK

diff1 = pKs['NTR'] - pKs['TYR']
diff2 = pKs['CTR'] - pKs['TYR']
diff3 = pKs['NTR'] - pKs['CTR']

print(diff1, diff2, diff3)
