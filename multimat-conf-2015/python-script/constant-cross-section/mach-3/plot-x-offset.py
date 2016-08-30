#### import python modules ####
import sys, os, math, csv
import numpy as np
import itertools as it
from compute_mass_difference import compute_mass_diff
from compute_error_norms import compute_error_norms
from scipy.optimize import fmin
import matplotlib.pyplot as plt
from Lagrange_test_function import b
from decimal import *
#### import python modules ####

nb_cells = []
x_offset_mass, x_offset_energy = [], []

file_mass = sys.argv[1]
file_energy = sys.argv[2]
out_file = sys.argv[3]

# open file and read first line
file_data_mass=open(file_mass, 'r')
file_data_energy=open(file_energy, 'r')
line_head = file_data_mass.readline()
print 'Variables in file', file_mass,':', line_head[:-1]
line_head = file_data_energy.readline()
print 'Variables in file', file_energy,':', line_head[:-1]

# read remaining of the file and store data
for line_mass in file_data_mass:
  row = line_mass.split()
  nb_cells.append(float(row[0]))
  x_offset_mass.append(row[1])

for line_energy in file_data_energy:
  row = line_energy.split()
  x_offset_energy.append(row[1])

plt.plot(nb_cells, x_offset_mass, '-o', markersize=8, linewidth=2, label=r'$\Delta \rho$')
plt.plot(nb_cells, x_offset_energy, '-+', markersize=8, linewidth=2, label=r'$\Delta (\rho E_{tot})$')
plt.legend(loc='best', fontsize=20, frameon=False)
plt.xlabel(r'$cells$', fontsize=20)
plt.ylabel(r'$x_{offset}$', fontsize=20)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.margins(0.1, 0.1)
fig_name=out_file+'-x-offset.eps'
print 'Saving plot using Matplotlib:', fig_name
plt.savefig(fig_name)
plt.clf()