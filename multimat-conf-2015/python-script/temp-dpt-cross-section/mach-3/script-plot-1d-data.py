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
def plot_error_norms(nb_cells, L1_norm, L2_norm, variable):
  nb_cells_log = [math.log(x) for x in nb_cells]
  L1_norm_log = [math.log(x) for x in L1_norm]
#  L2_norm_log= [math.log(x) for x in L2_norm]
  plt.plot(nb_cells_log, L1_norm_log, '+-', label=r'$L_1^{error} norm$')
  x1 = [math.log(nb_cells[0]), math.log(nb_cells[-1])]
  y1 = [-math.log(nb_cells[0])+math.log(L1_norm[-1])+math.log(nb_cells[-1]), math.log(L1_norm[-1])]
  plt.plot(x1, y1, '-', label=r'$line \ of  \ slope \ 1$')
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
    y_label=r'$T_r$'
  elif variable=='mach-number':
    y_label=r'$Mach$'
  else:
    print 'ERROR: unvalid variable name'
    sys.exit()
  plt.ylabel(r'$\log (L_{1,2}^{error}($'+y_label+'$))$', fontsize=20)
  fig_name=out_file+'-'+variable+'-convergence.eps'
  plt.savefig(fig_name)
  plt.clf()

def plot_solution(x_num, y_num, x_anal, y_anal, x_offset, variable, nb_cells, text_on=False, mass_diff=0):
  x_anal_offset = [float(x)+float(x_offset) for x in x_anal]
  x_num = [float(x) for x in x_num]
  plt.plot(x_num, y_num, '+-', markevery=30, markersize=8, label=r'$numerical \ solution$', linewidth=2)
  plt.plot(x_anal_offset, y_anal, 'o-', markevery=5000, markersize=8, label=r'$exact \ solution$', linewidth=1.2)
  plt.legend(loc='best', fontsize=20, frameon=False)
  plt.xlabel(r'$x$', fontsize=20)
  if variable=='density':
    y_label=r'$\rho$'
  elif variable=='mat-temp':
    y_label=r'$T$'
  elif variable=='radiation':
    y_label=r'$T_r$'
  elif variable=='mach-number':
    y_label=r'$Mach$'
  else:
    print 'ERROR: unvalid variable name'
    sys.exit()
  plt.ylabel(y_label, fontsize=20)
#  plt.ylim(min(y_num), max(y_num))
  plt.xlim(float(min(x_num)), float(max(x_num)))
#  plt.xlim(-0.015, 0.015)
  plt.margins(0.1, 0.1)
  if text_on:
    plt.figtext(0.2,0.5,r'$\Delta \rho$='+str('{:.4e}'.format(mass_diff)), fontsize=14)
  fig_name=out_file+'-'+variable+'-nel-'+str(nb_cells)+'-plot.eps'
  plt.savefig(fig_name)
  print 'Saving plot using Matplotlib:', fig_name
  plt.clf()

def plot_visc_coeff(x_num, file, var_index, nb_cells):
  # open file and read first line
  visc_e = []
  visc_max = []
  file_data=open(file, 'r')
  line_head = file_data.readline()
  for line in file_data:
    row = line.split(',')
    visc_e.append(row[var_index[5]])
    visc_max.append(row[var_index[6]])
  x_num = [float(x) for x in x_num[:-1]]
  plt.plot(x_num, visc_max, '+-', markevery=30, markersize=8, label=r'$\kappa_{max}$', linewidth=2)
  plt.plot(x_num, visc_e, 'o-', markevery=30, markersize=8, label=r'$\kappa_e$', linewidth=2)
  plt.legend(loc='best', fontsize=20, frameon=False)
  plt.xlabel(r'$x$', fontsize=20)
  plt.ylabel(r'$\kappa$', fontsize=20)
  plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
  plt.xlim(float(min(x_num)), float(max(x_num)))
  plt.margins(0.1, 0.1)
  fig_name=out_file+'-visc-nel-'+str(nb_cells)+'-plot.eps'
  plt.savefig(fig_name)
  print 'Saving plot using Matplotlib:', fig_name
  plt.clf()
#### define function ####

# READ EXACT SOLUTION
file_exact_list = []
# x-coordinates
#file_exact_list.append('data_x.dat')
file_exact_list.append('x_data.txt')
x_coord_exact = []
file_data_exact=open(file_exact_list[-1], 'r')
x_coord_exact[:] = [ line[:-1] for line in file_data_exact]
# material density
#file_exact_list.append('data_Density.dat')
file_exact_list.append('Density_data.txt')
mat_density_exact = []
file_data_exact=open(file_exact_list[-1], 'r')
mat_density_exact[:] = [ line[:-1] for line in file_data_exact]
# radiation energy density
#file_exact_list.append('data_RED.dat')
file_exact_list.append('Er_data.txt')
radiation_exact = []
file_data_exact=open(file_exact_list[-1], 'r')
radiation_exact[:] = [ line[:-1] for line in file_data_exact]
# mach number or fluid velocity
#file_exact_list.append('data_Mach.dat')
file_exact_list.append('Mach_Data.txt')
mach_nb_exact = []
file_data_exact=open(file_exact_list[-1], 'r')
mach_nb_exact[:] = [ line[:-1] for line in file_data_exact]
# material temperature
#file_exact_list.append('data_Temp.dat')
file_exact_list.append('Tm_data.txt')
mat_temp_exact = []
file_data_exact=open(file_exact_list[-1], 'r')
mat_temp_exact[:] = [ line[:-1] for line in file_data_exact]
nb_nodes_exact = len(x_coord_exact)
nb_exact_files = len(file_exact_list)
# remove duplicate values in x_coord_exact and consistently in other lists
x_coord_exact_temp, idx = np.unique(np.asarray(x_coord_exact), return_index=True)
idx = sorted(idx)
x_coord_exact = [x_coord_exact[i] for i in idx]
mat_density_exact = [mat_density_exact[i] for i in idx]
radiation_exact = [radiation_exact[i] for i in idx]
mach_nb_exact = [mach_nb_exact[i] for i in idx]
mat_temp_exact = [mat_temp_exact[i] for i in idx]
# normalize mach number
mach_nb_exact = [ float(i)/float(mach_nb_exact[0]) for i in mach_nb_exact]

# SET INPUT FILES
file_list = []

## works well with cfl=0.1
file_list.append('mach-3-nel-200-points0.csv')
#file_list.append('mach-3-nel-300-points0.csv')
#file_list.append('mach-3-nel-400-points0.csv')
#file_list.append('mach-3-nel-500-points0.csv')
#file_list.append('mach-3-nel-600-points0.csv')
#file_list.append('mach-3-nel-700-points0.csv')
#file_list.append('mach-3-nel-800-points0.csv')
#file_list.append('mach-3-nel-900-points0.csv')
#file_list.append('mach-3-nel-1000-points0.csv')
#file_list.append('mach-3-nel-1100-points0.csv')
#file_list.append('mach-3-nel-1200-points0.csv')
#file_list.append('mach-3-nel-1300-points0.csv')

# SET SOME VARIABLES
dir_path = os.getcwd()
quad_order = 10
interp_kind = 'linear'
nb_files = len(file_list)
var_index = [11, 5, 1, 2, 8, 4, 3] # [x, rho, radiation, mach, mat temp]
var_index[:] = [i -1 for i in var_index] # convert to python index

# OUTPUT SOME INFORMATION
print '------------------------------'
print 'Number of files with exact solution to read:', nb_exact_files
print 'Number of nodes in exact solution:', nb_nodes_exact
print 'Number of files to read:', nb_files
print 'Index of variables in each file:', var_index

# VARIABLES TO STORE DATA
L1_norm_density = []
L2_norm_density = []
L1_norm_radiation = []
L2_norm_radiation = []
L1_norm_mach = []
L2_norm_mach = []
L1_norm_mat_temp = []
L2_norm_mat_temp = []
nb_cells = []
x_offset = []

mat_temp = []
mach_nb = []
radiation = []
mat_density = []
x_coord = []

# LOOP OVER FILES TO COMPUTE L2 and L1 norms
for file in file_list:
  print '------------------------------'
  print 'Name of the input file:', file
  # set/reset data
  out_file = file[:6]
  mat_temp[:] = []
  mach_nb[:] = []
  radiation[:] = []
  mat_density[:] = []
  x_coord[:] = []
  # open file and read first line
  file_data=open(file, 'r')
  line_head = file_data.readline()
  print 'Variables in file', file,':', line_head[:-1]
  # read remaining of the file and store data
  for line in file_data:
    row = line.split(',')
    x_coord.append(row[var_index[0]])
    mat_density.append(row[var_index[1]])
    radiation.append(row[var_index[2]])
    mach_nb.append(row[var_index[3]])
    mat_temp.append(row[var_index[4]])
  # normalize each physical variable
  mat_density = [ float(i)/float(mat_density[0]) for i in mat_density]
  radiation = [ float(i)/float(radiation[0]) for i in radiation]
  mach_nb = [ float(i)/float(mach_nb[0]) for i in mach_nb]
  mat_temp = [ float(i)/float(mat_temp[0]) for i in mat_temp]
  # output number of nodes for numerical mesh
  nb_cells.append(len(x_coord)-1)
  print'Number of cells in file', file, ':', nb_cells[-1]

#  mass_diff = compute_mass_diff(0., x_coord, mat_density, x_coord_exact, mat_density_exact, quad_order, interp_kind)
#  print mass_diff
#  sys.exit()
#  res = [0,mass_diff]
  # minimize the mass difference between the exact and numerical solutions to get 'x_offset'
##  res = fmin(compute_mass_diff, 0., args=(x_coord, mat_density, x_coord_exact, mat_density_exact, quad_order, interp_kind,), xtol=1., ftol=1., full_output=True, disp=True, retall=True)[0:2]
  res = fmin(compute_mass_diff, 0., args=(x_coord, mat_density, x_coord_exact, mat_density_exact, quad_order, interp_kind,), xtol=1.e-20, ftol=1e-10, full_output=True, disp=True, retall=True)[0:2]
##  res = minimize(compute_mass_diff, 2.e-4, args=(x_coord, mat_density, x_coord_exact, mat_density_exact, quad_order, interp_kind,), method='nelder-mead', options={'xtol': 1e-4, 'disp': True, 'maxiter' : 10000})
  x_offset.append(float(res[0]))
  mass_diff = res[1]
  print 'x offset for', file, 'is', x_offset[-1]
  print 'mass difference for', file, 'is', mass_diff
  # save density plot
  plot_solution(x_coord, mat_density, x_coord_exact, mat_density_exact, x_offset[-1], 'density', nb_cells[-1], True, mass_diff)
  # compute the error norms for density
  l1_norm, l2_norm = compute_error_norms(x_offset[-1], x_coord, mat_density, x_coord_exact, mat_density_exact, quad_order, interp_kind)
  print l1_norm
  L1_norm_density.append(l1_norm)
  L2_norm_density.append(l2_norm)
  # compute the error norms for radiation
  l1_norm, l2_norm = compute_error_norms(x_offset[-1], x_coord, radiation, x_coord_exact, radiation_exact, quad_order, interp_kind)
  print l1_norm
  L1_norm_radiation.append(l1_norm)
  L2_norm_radiation.append(l2_norm)
  # compute the error norms for mach number
  l1_norm, l2_norm = compute_error_norms(x_offset[-1], x_coord, mach_nb, x_coord_exact, mach_nb_exact, quad_order, interp_kind)
  print l1_norm
  L1_norm_mach.append(l1_norm)
  L2_norm_mach.append(l2_norm)
  # compute the error norms for material temperature
  l1_norm, l2_norm = compute_error_norms(x_offset[-1], x_coord, mat_temp, x_coord_exact, mat_temp_exact, quad_order, interp_kind)
  print l1_norm
  L1_norm_mat_temp.append(l1_norm)
  L2_norm_mat_temp.append(l2_norm)

  file_data.close()
#  del mat_temp, mach_nb, radiation, mat_density, x_coord

# PLOT L1 AND L2 NORMS
plot_error_norms(nb_cells, L1_norm_density, L2_norm_density, 'density')
plot_error_norms(nb_cells, L1_norm_radiation, L2_norm_radiation, 'radiation')
plot_error_norms(nb_cells, L1_norm_mach, L2_norm_mach, 'mach-number')
plot_error_norms(nb_cells, L1_norm_mat_temp, L2_norm_mat_temp, 'mat-temp')

# PLOT MATERIAL DENSITY AND TEMPERATURE, RADIATION TEMPERATURE AND MACH NUMBER
#plot_solution(x_coord, mat_density, x_coord_exact, mat_density_exact, x_offset[-1], 'density')
#plot_solution(x_coord, radiation, x_coord_exact, radiation_exact, x_offset[-1], 'radiation', nb_cells[-1])
#plot_solution(x_coord, mach_nb, x_coord_exact, mach_nb_exact, x_offset[-1], 'mach-number', nb_cells[-1])
#plot_solution(x_coord, mat_temp, x_coord_exact, mat_temp_exact, x_offset[-1], 'mat-temp', nb_cells[-1])

plot_solution(x_coord, mat_temp, x_coord_exact, mat_temp_exact, x_offset[-1], 'density', nb_cells[-1])
plot_solution(x_coord, mat_temp, x_coord, mat_temp, 0.0, 'radiation', 200)
plot_solution(x_coord_exact, mat_temp_exact, x_coord_exact, mat_temp_exact, 0.0, 'mat-temp', 200)
exit()

#plot_visc_coeff(x_coord, 'visc-coeff-nel-1000-points0.csv', var_index, nb_cells[-1])
#plot_visc_coeff(x_coord, 'visc-coeff-nel-1600-points0.csv', var_index, nb_cells[-1])

# PLOT X_OFFSET AND SAVE VALUES IN FILE
# plot
plt.plot(nb_cells, x_offset, '-o', markersize=8, linewidth=2)
plt.xlabel(r'$cells$', fontsize=20)
plt.ylabel(r'$x_{offset}$', fontsize=20)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
plt.margins(0.1, 0.1)
fig_name=out_file+'-x-offset.eps'
plt.savefig(fig_name)
plt.clf()
# save
datafile_id = open('x_offset.txt', 'w+')
datafile_id.write('nb_cells '+' x_offset \n')
data = np.asarray([nb_cells, x_offset])
data = data.T
np.savetxt(datafile_id, data, fmt=['%d','%.4e'])
datafile_id.close()