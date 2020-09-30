#!/usr/bin/env python3

import numpy as np
from math import sin, cos, pi

seq_N = [10, 25, 50, 100, 200, 350, 750]
seq_rho = [i/2 for i in range(8, 15)]

exact = [3, 7, 11, 15]

M = np.zeros((len(seq_N), len(seq_rho)))
for i, N in enumerate(seq_N):
	for j, rho in enumerate(seq_rho):
		with open('lambda_%u_%f.csv' % (N, rho)) as f: csv = f.read()
		sq = np.array([float(x) for x in csv.strip().split('\n')]) - exact
		sq = np.sqrt(np.sum(sq ** 2))
		M[i, j] = sq

np.set_printoptions(precision=6, suppress=True)
print (M)
