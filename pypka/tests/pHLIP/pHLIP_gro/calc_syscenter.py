import sys

x_min = 99999
x_max = -1
y_min = 99999
y_max = -1
z_min = 99999
z_max = -1

with open(sys.argv[1]) as f:
    for line in f:
        if line.startswith('ATOM'):
            cols = line.split()
            if cols[3] not in ('PDB', 'CNT'):
                x = float(cols[5])
                y = float(cols[6])
                z = float(cols[7])

                if x < x_min:
                    x_min = x
                if x > x_max:
                    x_max = x

                if y < y_min:
                    y_min = y
                if y > y_max:
                    y_max = y

                if z < z_min:
                    z_min = z
                if z > z_max:
                    z_max = z

print x_min, x_max
print (x_min + x_max) / 2.0, (y_min + y_max) / 2.0, (z_min + z_max) / 2.0

