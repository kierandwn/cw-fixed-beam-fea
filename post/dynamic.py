#!/usr/bin/python
from sys import argv
from time import sleep

# import numpy as np
import matplotlib.pyplot as plt
import math

fResult = open("./output/dynamic.txt")

for i in range(1, len(argv)):
	if (argv[i] == '--length'):
		L  = float(argv[i + 1])
	elif (argv[i] == '--Nx'):
		Nx = float(argv[i + 1])


u = []
t = []

for line in fResult:
	vals = line.split(' ')
	
	t.append(float(vals[0]))
	u.append(float(vals[1]))

fResult.close()

plt.plot(t, u, "-")

plt.grid('on') 

plt.xlabel('Time (s)')
plt.ylabel('Maximum Deflection (mm)')

# plt.show()

plt.savefig("./output/task_.png")