import sys

c6, c12 = sys.argv[1:3]
sigma = ((c12) / (c6)) ** (1.0 / 6.0)
radius = sigma * 10 / 2.0
print(radius)
