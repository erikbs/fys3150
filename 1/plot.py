#!/usr/bin/env python3

import matplotlib.pyplot as plt, numpy as np
from matplotlib import rc

# Analytisk løysing
x = np.linspace(0, 1, 1000)
u = 1 - (1 - np.exp(-10))*x - np.exp(-10*x)

# Set opp LaTeX
rc('font', family='sans-serif', weight='bold', size=12)
rc('text.latex', preamble=[r'\usepackage{tgheros}', r'\usepackage{sansmathfonts}'])
rc('text', usetex=True)
rc('legend', fontsize=12)

# Lagar eit tri-i-eitt-plott til u og v for n = 10, 100, 1000
fig = plt.figure()
E = (1, 2, 3)
for e in E:
	n = 10**e
	h = 1 / (n + 1.)
	xi = [i * h for i in range(n + 2)]
	with open('v-gn_%u.csv' % n) as f: v = [0,] + [float(i) for i in f.read().split()] + [0,]

	plt.subplot(1, len(E), e)
	plt.xlim(0, 1)
	plt.xlabel('$x$')
	plt.ylabel('$y$')
	plt.title('$n=%u$' % n)
	plt.plot(x, u, linewidth=2, linestyle='-')
	plt.plot(xi, v, linewidth=1, linestyle='none', marker='.')

fig.legend(['$y=u(x)$', '$y=v(x)$'], loc='lower center', ncol=2)
fig.tight_layout(pad=0.1)
fig.subplots_adjust(top=0.9, bottom=0.3)
w, h = fig.get_size_inches()
fig.set_size_inches([w, w/2])
plt.savefig('uv-plot.pdf', orientation='landscape')

# Lagar plott til køyretidene
n    = [   10**1,    10**2,    10**3,     10**4,    10**5,    10**6,    10**7,    10**8]
t_gn = [0.000005, 0.000012, 0.000092,  0.000852, 0.008176, 0.084792, 0.222927, 2.861431]
t_sp = [0.000002, 0.000003, 0.000064,  0.000514, 0.007103, 0.064417, 0.197545, 2.263835]
t_lu = [0.000196, 0.001877, 0.127761, 23.819581]

fig = plt.figure()
plt.plot(n, t_gn, n, t_sp, n[:len(t_lu)], t_lu, marker='.')
plt.legend(['generell algoritme', 'spesialisert algoritme', 'lu-faktorisering'], loc='lower right')
plt.xlabel('$n$')
plt.ylabel('Køyretid [s]')
plt.xscale('log', basex=10, subsx=[])
plt.yscale('log')
fig.tight_layout(pad=1)
w, h = fig.get_size_inches()
fig.set_size_inches([w, w/1.8])
plt.savefig('t-plot.pdf', orientation='landscape')

# Lagar plott til relativt fråvik
h = [-np.log10(i + 1) for i in n]
eps = [-1.179698, -3.087196, -4.442024, -4.159437, -2.527016, -1.600981, -0.632154, -0.229417]

fig = plt.figure()
plt.plot(h, eps)
plt.xlabel(r'$\log_{10}(h)$')
plt.ylabel(r'$\epsilon=\log_{10}\left(\displaystyle\max_{i=1}^n\left|\frac{v_i-u_i}{u_i}\right|\right)$')
fig.tight_layout(pad=1.5)
w, h = fig.get_size_inches()
fig.set_size_inches([w, w/2])
plt.savefig('eps-plot.pdf', orientation='landscape')
