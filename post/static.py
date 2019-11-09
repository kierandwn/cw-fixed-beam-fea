#!/usr/bin/python
from sys import argv

import numpy as np
import matplotlib.pyplot as plt
# import math

L  = False
Nx = False
A  = False
I  = False
E  = False

for i in range(1, len(argv)):
	if (argv[i] == '--length'):
		L  = float(argv[i + 1])
	elif (argv[i] == '--Nx'):
		Nx = float(argv[i + 1])
	elif (argv[i] == '--area'):
		A  = float(argv[i + 1])
	elif (argv[i] == '--I'):
		I  = float(argv[i + 1])
	elif (argv[i] == '--E'):
		E  = float(argv[i + 1])

if (not L):
	L = 10.0
if (not Nx):
	Nx = 24
if (not A):
	A = 0.012
if (not I):
	I = 0.0000144
if (not E):
	E = 210e+09;

# Define loading conditions
qy = -   1.0
Fy = -1000.0

fResult = open("./output/static.txt")

X = []
u = []

i = 0;

for line in fResult: 
	vals = line.split(' ')
	X.append(float(vals[0]))
	u.append(float(vals[2]))

# Enforce boundary conditions
X = [0] + X + [L]
u = [0] + u + [0]

ua = []

# Mesh resolution for analytical solution.
na = 1000
xa = np.linspace(0, L, na)

# Solve the analytical solution
for i in range(0, na):
	# Contribution from distributed load
	ua.append( ((qy * xa[i] * xa[i]) / (24 * E * I)) * (L - xa[i]) * (L - xa[i]) )
	
	# Contribution from concentrated load
	if (xa[i] <= L / 2):
		ua[i] += ((Fy * xa[i] * xa[i]) / (48 * E * I)) * (3 * L - 4 * xa[i])
	else:
		ua[i] += ((Fy * (L - xa[i]) * (L - xa[i])) / (48 * E * I)) * (3 * L - 4 * (L - xa[i]))

plt.plot(X, u, 'r', xa, ua, 'k--')
plt.grid('on') 

plt.xlabel('Axial Co-ordinate, (mm)')
plt.ylabel('Transverse Deflection (mm)')

# plt.show()

plt.savefig("./output/task1.png")