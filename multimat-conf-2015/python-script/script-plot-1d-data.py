# This script aims at reading csv files and store the data for postprocessing in 1-D.

#### import python modules ####
import sys, os, math, csv
import numpy as np
import itertools as it
from compute_mass_difference import compute_mass_diff
from compute_error_norms import compute_error_norms
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from Lagrange_test_function import b
from decimal import *
#### import python modules ####

#### define function ####
def plot_error_norms(nb_cells, L1_norm, L2_norm, variable):
  nb_cells_log = [math.log(x) for x in nb_cells]
  L1_norm_log = [math.log(x) for x in L1_norm]
  L2_norm_log= [math.log(x) for x in L2_norm]
  plt.plot(nb_cells_log, L1_norm_log, '+-', label=r'$L_1^{error} norm$')
  x1 = [math.log(nb_cells[0]), math.log(nb_cells[-1])]
  y1 = [-2*math.log(nb_cells[0])+math.log(L1_norm[-1])+2*math.log(nb_cells[-1]), math.log(L1_norm[-1])]
  plt.plot(x1, y1, '-', label='line of slope 2')
  plt.plot(nb_cells_log, L2_norm_log, 'o-', label=r'$L_2^{error} norm$')
  y2 = [-2*math.log(nb_cells[0])+math.log(L2_norm[-1])+2*math.log(nb_cells[-1]), math.log(L2_norm[-1])]
  plt.plot(x1, y2, '-', label='line of slope 2')
  plt.legend(loc='upper right')
  plt.xlabel(r'$\log (cells)$')
  plt.ylabel(r'$\log (L_{1,2}^{error})$ for '+variable)
  fig_name=out_file+'-'+variable+'-convergence.eps'
  plt.savefig(fig_name)
  plt.clf()
#### define function ####

# READ EXACT SOLUTION
file_exact_list = []
path_to_exact_files = 'test'
# x-coordinates
file_exact_list.append('x_1p05.txt')
x_coord_exact = []
file_data_exact=open(file_exact_list[-1], 'r')
x_coord_exact[:] = [ line[:-1] for line in file_data_exact]
# material density
file_exact_list.append('rho_1p05.txt')
mat_density_exact = []
file_data_exact=open(file_exact_list[-1], 'r')
mat_density_exact[:] = [ line[:-1] for line in file_data_exact]
# radiation energy density
file_exact_list.append('Er_1p05.txt')
radiation_exact = []
file_data_exact=open(file_exact_list[-1], 'r')
radiation_exact[:] = [ line[:-1] for line in file_data_exact]
# mach number or fluid velocity
file_exact_list.append('mach_1p05.txt')
mach_nb_exact = []
file_data_exact=open(file_exact_list[-1], 'r')
mach_nb_exact[:] = [ line[:-1] for line in file_data_exact]
# material temperature
file_exact_list.append('T_1p05.txt')
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
#file_list.append('mach-1p05-nel-20-points0.csv')
#file_list.append('mach-1p05-nel-30-points0.csv')
#file_list.append('mach-1p05-nel-40-points0.csv')
#file_list.append('mach-1p05-nel-50-points0.csv')
#file_list.append('mach-1p05-nel-60-points0.csv')
#file_list.append('mach-1p05-nel-70-points0.csv')
#file_list.append('mach-1p05-nel-80-points0.csv')
#file_list.append('mach-1p05-nel-90-points0.csv')
file_list.append('mach-1p05-nel-100-points0.csv')
#file_list.append('mach-1p05-nel-110-points0.csv')
#file_list.append('mach-1p05-nel-120-points0.csv')
#file_list.append('mach-1p05-nel-130-points0.csv')
#file_list.append('mach-1p05-nel-140-points0.csv')
#file_list.append('mach-1p05-nel-150-points0.csv')
#file_list.append('mach-1p05-nel-160-points0.csv')
#file_list.append('mach-1p05-nel-170-points0.csv')
#file_list.append('mach-1p05-nel-180-points0.csv')
#file_list.append('mach-1p05-nel-190-points0.csv')
#file_list.append('mach-1p05-nel-200-points0.csv')
#file_list.append('mach-1p05-nel-210-points0.csv')
#file_list.append('mach-1p05-nel-220-points0.csv')
#file_list.append('mach-1p05-nel-230-points0.csv')
#file_list.append('mach-1p05-nel-240-points0.csv')
#file_list.append('mach-1p05-nel-250-points0.csv')

# SET SOME VARIABLES
dir_path = os.getcwd()
quad_order = 20
nb_files = len(file_list)
var_index = [11, 5, 1, 2, 8] # [x, rho, radiation, mach, mat temp]
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

# LOOP OVER FILES TO COMPUTE L2 and L1 norms
for file in file_list:
  print '------------------------------'
  print 'Name of the input file:', file
  # set/reset data
  out_file = file[:9]
  mat_temp = []
  mach_nb = []
  radiation = []
  mat_density = []
  x_coord = []
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

#  mass_diff = compute_mass_diff(2.e-4, x_coord, mat_density, x_coord_exact, mat_density_exact, quad_order)
  # minimize the mass difference between the exact and numerical solutions to get 'x_offset'
  res = minimize(compute_mass_diff, 2.e-4, args=(x_coord, mat_density, x_coord_exact, mat_density_exact, quad_order,), method='nelder-mead',
                 options={'xtol': 1e-20, 'disp': True, 'maxiter' : 1000})
  x_offset.append(res.x[0])
  mass_diff = res.fun
  print 'x offset for', file, 'is', x_offset[-1]
  print 'mass difference for', file, 'is', mass_diff
  # save density plot
  x_coord_exact_offset = [float(x)+float(x_offset[-1]) for x in x_coord_exact]
  plt.plot(x_coord_exact_offset, mat_density_exact, label='exact solution')
  plt.plot(x_coord, mat_density, label = 'numerical solution')
  plt.xlabel('x (m)')
  plt.ylabel(r'$\rho \ (g/cm^3)$')
  fig_name=out_file+'-density-nel-'+str(nb_cells[-1])+'.eps'
  plt.legend(loc='upper left')
  plt.figtext(0.2,0.7,r'$x_{offset}$='+str('{:.4e}'.format(x_offset[-1])), fontsize=14)
  plt.figtext(0.2,0.5,r'$\Delta \rho$='+str('{:.4e}'.format(mass_diff)), fontsize=14)
  plt.savefig(fig_name)
  plt.clf()
  print 'Saving plot using Matplotlib:', fig_name
  # compute the error norms for density
  l1_norm, l2_norm = compute_error_norms(x_offset[-1], x_coord, mat_density, x_coord_exact, mat_density_exact, quad_order)
  L1_norm_density.append(l1_norm)
  L2_norm_density.append(l2_norm)
  # compute the error norms for radiation
  l1_norm, l2_norm = compute_error_norms(x_offset[-1], x_coord, radiation, x_coord_exact, radiation_exact, quad_order)
  L1_norm_radiation.append(l1_norm)
  L2_norm_radiation.append(l2_norm)
  # compute the error norms for mach number
  l1_norm, l2_norm = compute_error_norms(x_offset[-1], x_coord, mach_nb, x_coord_exact, mach_nb_exact, quad_order)
  L1_norm_mach.append(l1_norm)
  L2_norm_mach.append(l2_norm)
  # compute the error norms for material temperature
  l1_norm, l2_norm = compute_error_norms(x_offset[-1], x_coord, mat_temp, x_coord_exact, mat_temp_exact, quad_order)
  L1_norm_mat_temp.append(l1_norm)
  L2_norm_mat_temp.append(l2_norm)

  file_data.close()
  del mat_temp, mach_nb, radiation, mat_density, x_coord

# PLOT L1 AND L2 NORMS
plot_error_norms(nb_cells, L1_norm_density, L2_norm_density, 'density')
plot_error_norms(nb_cells, L1_norm_radiation, L2_norm_radiation, 'radiation')
plot_error_norms(nb_cells, L1_norm_mach, L2_norm_mach, 'mach-number')
plot_error_norms(nb_cells, L1_norm_mat_temp, L2_norm_mat_temp, 'mat-temp')

# PLOT X_OFFSET AND SAVE VALUES IN FILE
# plot
plt.plot(nb_cells, x_offset, '-o')
plt.xlabel('number of degree of freedom')
plt.ylabel(r'$x_{offset}$ '+r'$cm$')
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