#!/usr/bin/python
import matplotlib.pyplot as plt
import numpy as np
import os

def max_value (array, r):
	maximum = -1 * float("inf");

	for i in range(0, r):
		if (array[i] > maximum):
			maximum = array[i]

	return maximum

def min_value (array, r):
	minimum = float("inf");

	for i in range(0, r):
		if (array[i] < minimum):
			minimum = array[i]

	return minimum

os.system("make")

sample_size = 100

Tl  = np.linspace(1, 50, sample_size)

dt = 1e-04

amplitudes = []

for tl in Tl:

	T = 1.2 * tl
	Nt = tl / dt

	print "Iteration:", len(amplitudes)
	print "...Loading Time:", tl
	print "...Number of Timesteps:", Nt

	CMD = "./main --dynamic --implicit --Nt " + str(int(Nt)) + " --Tl " + str(tl) + " --T " + str(T)
	os.system(CMD)

	fResult = open("./output/dynamic.txt")

	u = []
	t = []

	for line in fResult:
		vals = line.split(' ')
		
		if float(vals[0]) >= tl:
			t.append(float(vals[0]))
			u.append(float(vals[1]))

	fResult.close()

	os.system("rm ./.build/dynamic.txt")

	amp = max_value(u, len(u)) - min_value(u, len(u))
	amplitudes.append(amp)

	print "...Complete."


plt.loglog(Tl, amplitudes, "-")


plt.grid('on') 

plt.xlabel('Loading Time (s)')
plt.ylabel('Amplitude of Oscillation (mm)')

plt.savefig("./output/amplitudes3c.png")