#!/usr/bin/env python3

import matplotlib.pyplot as plt, numpy as np
from matplotlib import rc
from math import sin, cos, pi

# Set opp LaTeX
rc('font', family='sans-serif', weight='bold', size=12)
rc('text.latex', preamble=[r'\usepackage{tgheros}', r'\usepackage{sansmathfonts}'])
rc('text', usetex=True)
rc('legend', fontsize=12)

def lj(j, N):
	return 2 -2 * cos(j * pi / N)

def uj(j, N):
	u = np.array([sin(i * j * pi / N) for i in range(1, N)])

	# Normaliser eigenvektorane
	return u / np.sqrt((u**2).sum())

# Plott sann løysing mot numerisk
ev = np.genfromtxt('evec_256.csv', delimiter=',')
vec = ev[:,0]
N = len(vec) + 1
exact = uj(1, N)
x=np.linspace(0, 1, N+1)[1:-1]
fig = plt.figure()
plt.plot(x, vec, linestyle='-')
plt.plot(x, exact, linestyle='--')
plt.legend(['jacobi', 'analytisk'])
plt.xlabel('$u_j$')
plt.ylabel('')
fig.tight_layout(pad=1.5)
w, h = fig.get_size_inches()
fig.set_size_inches([w, w/1.8])
plt.savefig('jacobi-plot.pdf', orientation='landscape')

C = np.genfromtxt('gstate_350_5.500000.csv', delimiter=',')
a,b,c,d=C[:,0],C[:,1],C[:,2],C[:,3]
f=np.linspace(0, 5.5, 351)[1:-1]

# Plott u(ρ) for ulike ω
fig = plt.figure()
plt.plot(f,a,f,b,f,c,f,d)
plt.legend(['$\omega_r=0.05$', '$\omega_r=0.5$', '$\omega_r=1$', '$\omega_r=5$'], loc='upper right')
plt.xlabel(r'$\rho$')
plt.ylabel(r'$u(\rho)$');
fig.tight_layout(pad=1.5)
w, h = fig.get_size_inches()
fig.set_size_inches([w, w/1.8])
plt.savefig('schrod2-plot.pdf', orientation='landscape')
