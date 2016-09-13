# This script aims at reading csv files and store the data for postprocessing in 1-D.

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

#### define function ####
def plot_error_norms(nb_cells, l1_norm_mass, l1_norm_energy, variable):
  nb_cells_log = [math.log(float(x)) for x in nb_cells]
  l1_norm_mass_log = [math.log(float(x)) for x in l1_norm_mass]
  l1_norm_energy_log = [math.log(float(x)) for x in l1_norm_energy]
  plt.plot(nb_cells_log, l1_norm_mass_log, '+-', label=r'$L_1^{error} norm \ \Delta \rho$', linewidth=2, markersize=8)
  plt.plot(nb_cells_log, l1_norm_energy_log, 'o-', label=r'$L_1^{error} norm \ \Delta (\rho E)_{tot}$', linewidth=2, markersize=8)
  x1 = [math.log(nb_cells[0]), math.log(nb_cells[-1])]
  a = 0.5*(float(l1_norm_mass[-1])+float(l1_norm_energy[-1]))
  y1 = [-math.log(nb_cells[0])+math.log(a)+math.log(nb_cells[-1]), math.log(a)]
  plt.plot(x1, y1, '-', label=r'$line \ of  \ slope \ 1$', linewidth=2)
  #  plt.plot(nb_cells_log, L2_norm_log, 'o-', label=r'$L_2^{error} norm$')
  #  y2 = [-math.log(nb_cells[0])+math.log(L2_norm[-1])+math.log(nb_cells[-1]), math.log(L2_norm[-1])]
  #  plt.plot(x1, y2, '-', label=r'$line \ of  \ slope \ 2$')
  plt.legend(loc='best', fontsize=20, frameon=False)
  plt.xlabel(r'$\log (cells)$', fontsize=20)
  if variable=='density':
    y_label=r'$\rho$'
  elif variable=='mat-temp':
    y_label=r'$T$'
  elif variable=='radiation':
    y_label=r'$\epsilon$'
  elif variable=='mach-number':
    y_label=r'$Mach$'
  elif variable=='total':
    y_label=r'$Total$'
  else:
    print 'ERROR: unvalid variable name'
    sys.exit()
  plt.ylabel(r'$\log (L_{1}^{error}($'+y_label+'$))$', fontsize=20)
  fig_name=out_file+'-'+variable+'-convergence.eps'
  print fig_name
  plt.savefig(fig_name)
  plt.clf()

file_mass = sys.argv[1]
file_energy = sys.argv[2]
variable = sys.argv[3]
out_file = variable

nb_cells = []
l1_norm_mass_density = []
l1_norm_mass_mach = []
l1_norm_mass_radiation = []
l1_norm_mass_temp = []
l1_norm_energy_density = []
l1_norm_energy_mach = []
l1_norm_energy_radiation = []
l1_norm_energy_temp = []

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
  l1_norm_mass_density.append(row[1])
  l1_norm_mass_temp.append(row[2])
  l1_norm_mass_mach.append(row[3])
  l1_norm_mass_radiation.append(row[4])
for line_energy in file_data_energy:
  row = line_energy.split()
  l1_norm_energy_density.append(row[1])
  l1_norm_energy_temp.append(row[2])
  l1_norm_energy_mach.append(row[3])
  l1_norm_energy_radiation.append(row[4])

plot_error_norms(nb_cells, l1_norm_mass_density, l1_norm_energy_density, 'density')
plot_error_norms(nb_cells, l1_norm_mass_temp, l1_norm_energy_temp, 'mat-temp')
plot_error_norms(nb_cells, l1_norm_mass_temp, l1_norm_energy_temp, 'mach-number')
plot_error_norms(nb_cells, l1_norm_mass_radiation, l1_norm_energy_radiation, 'radiation')