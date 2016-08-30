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
  elif variable=='total':
    y_label=r'$Total$'
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
  plt.plot(x_anal_offset, y_anal, 'o-', markevery=500, markersize=8, label=r'$exact \ solution$', linewidth=1.2)
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
path_to_exact_files = 'test'
# x-coordinates
file_exact_list.append('x_data.txt')
#file_exact_list.append('data_x.dat')
x_coord_exact = []
file_data_exact=open(file_exact_list[-1], 'r')
x_coord_exact[:] = [ line[:-1] for line in file_data_exact]
# material density
file_exact_list.append('Density_data.txt')
#file_exact_list.append('data_Density.dat')
mat_density_exact = []
file_data_exact=open(file_exact_list[-1], 'r')
mat_density_exact[:] = [ line[:-1] for line in file_data_exact]
# radiation energy density
file_exact_list.append('Er_data.txt')
#file_exact_list.append('data_RED.dat')
radiation_exact = []
file_data_exact=open(file_exact_list[-1], 'r')
radiation_exact[:] = [ line[:-1] for line in file_data_exact]
# mach number or fluid velocity
file_exact_list.append('Mach_Data.txt')
#file_exact_list.append('data_Mach.dat')
mach_nb_exact = []
file_data_exact=open(file_exact_list[-1], 'r')
mach_nb_exact[:] = [ line[:-1] for line in file_data_exact]
# material temperature
file_exact_list.append('Tm_data.txt')
#file_exact_list.append('data_Temp.dat')
mat_temp_exact = []
file_data_exact=open(file_exact_list[-1], 'r')
mat_temp_exact[:] = [ line[:-1] for line in file_data_exact]
nb_nodes_exact = len(x_coord_exact)
nb_exact_files = len(file_exact_list)
# remove duplicate values in x_coord_exact and consistently in other lists
x_coord_exact_temp, idx = np.unique(np.asarray(x_coord_exact), return_index=True)
idx = sorted(idx)
x_coord_exact = [x_coord_exact[i] for i in idx]
mat_density_exact = [float(mat_density_exact[i])/float(mat_density_exact[0]) for i in idx]
radiation_exact = [float(radiation_exact[i])/float(radiation_exact[0]) for i in idx]
mach_nb_exact = [float(mach_nb_exact[i])/float(mach_nb_exact[0]) for i in idx]
mat_temp_exact = [float(mat_temp_exact[i])/float(mat_temp_exact[0]) for i in idx]
gamma = 5/3
P_inf = 8.5319737603665362e-05 # 1.e-5
total_nrg_exact = []
total_nrg_exact = [float(mat_temp_exact[i])*float(mat_density_exact[i])*(1+0.5*gamma*(gamma-1)*float(mach_nb_exact[i])**2)+P_inf*float(radiation_exact[i]) for i in range(len(mat_temp_exact))]
# normalize mach number
#mach_nb_exact = [ float(i)/float(mach_nb_exact[0]) for i in mach_nb_exact]

# SET INPUT FILES
file_list = []

#file_list.append('mach-3-nel-250-points0.csv')
#file_list.append('mach-3-nel-300-points0.csv')

# above 1e-10 mass difference
#file_list.append('mach-3-nel-350-points0.csv')
#file_list.append('mach-3-nel-400-points0.csv')
#file_list.append('mach-3-nel-500-points0.csv')
##file_list.append('mach-3-nel-600-points0.csv')
#file_list.append('mach-3-nel-650-points0.csv')
#file_list.append('mach-3-nel-700-points0.csv')
#file_list.append('mach-3-nel-850-points0.csv')
#file_list.append('mach-3-nel-900-points0.csv')
##file_list.append('mach-3-nel-1000-points0.csv')
##file_list.append('mach-3-nel-1500-points0.csv')
#file_list.append('mach-3-nel-2200-points0.csv')

# below 1e-10 mass difference
#file_list.append('mach-3-nel-450-points0.csv')
#file_list.append('mach-3-nel-550-points0.csv')
#file_list.append('mach-3-nel-750-points0.csv')
#file_list.append('mach-3-nel-800-points0.csv')
#file_list.append('mach-3-nel-950-points0.csv')
#file_list.append('mach-3-nel-1200-points0.csv')
#file_list.append('mach-3-nel-2000-points0.csv')
#file_list.append('mach-3-nel-2500-points0.csv')

#file_list.append('mach-3-nel-300-points0.csv')
#file_list.append('mach-3-nel-400-points0.csv')
#file_list.append('mach-3-nel-500-points0.csv')
##file_list.append('mach-3-nel-600-points0.csv')
#file_list.append('mach-3-nel-700-points0.csv')
#file_list.append('mach-3-nel-800-points0.csv')
#file_list.append('mach-3-nel-900-points0.csv')
#file_list.append('mach-3-nel-1000-points0.csv')
#file_list.append('mach-3-nel-2000-points0.csv')

## works well with cfl=0.1
#file_list.append('mach-3-nel-100-points0.csv')
#file_list.append('mach-3-nel-200-points0.csv') # mass diff = 1.e-8
#file_list.append('mach-3-nel-300-points0.csv')
#file_list.append('mach-3-nel-400-points0.csv')
#file_list.append('mach-3-nel-500-points0.csv') # mass diff = 1.e-8 energy diff = 1.e-6
#file_list.append('mach-3-nel-600-points0.csv')
#file_list.append('mach-3-nel-700-points0.csv') # energy diff = 1.e-6
##file_list.append('mach-3-nel-800-points0.csv')
#file_list.append('mach-3-nel-900-points0.csv')
##file_list.append('mach-3-nel-1000-points0.csv')
#file_list.append('mach-3-nel-1100-points0.csv')
##file_list.append('mach-3-nel-1200-points0.csv')
#file_list.append('mach-3-nel-1300-points0.csv')
#file_list.append('mach-3-nel-1400-points0.csv') # mass diff = 1.e-7
#
#file_list.append('mach-3-nel-1500-points0.csv')
##file_list.append('mach-3-nel-1600-points0.csv')
##file_list.append('mach-3-nel-1700-points0.csv')

# for energy and mass conservation with 50k nodes in semi-analytical solutions
file_list.append('mach-3-nel-300-points0.csv')
file_list.append('mach-3-nel-500-points0.csv')
file_list.append('mach-3-nel-700-points0.csv')
file_list.append('mach-3-nel-900-points0.csv')
file_list.append('mach-3-nel-1100-points0.csv')
file_list.append('mach-3-nel-1300-points0.csv')
file_list.append('mach-3-nel-1500-points0.csv')
file_list.append('mach-3-nel-1700-points0.csv')

# SET SOME VARIABLES
dir_path = os.getcwd()
quad_order = 10 # 20 # 70 # 100
interp_kind = 'linear'
nb_files = len(file_list)
var_index = [11, 5, 1, 2, 8, 4, 3] # [x, rho, radiation, mach, mat temp]
#var_index = [14, 8, 1, 2, 11, 7, 6] # [x, rho, radiation, mach, mat temp, rad temp, ]
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
L1_norm_total = []
L2_norm_total = []
nb_cells = []
x_offset = []

mat_temp = []
mach_nb = []
radiation = []
mat_density = []
x_coord = []

out_file_base = sys.argv[1]

# LOOP OVER FILES TO COMPUTE L2 and L1 norms
for file in file_list:
  print '------------------------------'
  print 'Name of the input file:', file
  # set/reset data
  out_file = out_file_base+'-'+file[:6]
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
  total_nrg = [float(mat_temp[i])*float(mat_density[i])*(1+0.5*gamma*(gamma-1)*float(mach_nb[i])**2)+P_inf*float(radiation[i]) for i in range(len(mat_temp))]

#  print mat_temp[-1], mat_density[-1], mach_nb[-1], radiation[-1]
  # output number of nodes for numerical mesh
  nb_cells.append(len(x_coord)-1)
  print'Number of cells in file', file, ':', nb_cells[-1]

#  plot_solution(x_coord, total_nrg, x_coord_exact, total_nrg_exact, 0., 'density', nb_cells[-1])
#  print 'done plotting'

# minimize the energy difference between the exact and numerical solutions to get 'x_offset'
#  res = fmin(compute_mass_diff, 0., args=(x_coord, total_nrg, x_coord_exact, total_nrg_exact, quad_order, interp_kind,), xtol=1.e-20, ftol=1e-10, full_output=True, disp=True, retall=True, maxiter=10000000, maxfun=1000)[0:2]

  # minimize the mass difference between the exact and numerical solutions to get 'x_offset'
  res = fmin(compute_mass_diff, 0., args=(x_coord, mat_density, x_coord_exact, mat_density_exact, quad_order, interp_kind,), xtol=1.e-20, ftol=1e-10, full_output=True, disp=True, retall=True, maxiter=10000000, maxfun=1000)[0:2]

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

  L1_norm_total.append(L1_norm_density[-1]+L1_norm_radiation[-1]+L1_norm_mach[-1]+L1_norm_mat_temp[-1])
  print L1_norm_total[-1]
  L2_norm_total.append(L2_norm_density[-1]+L2_norm_radiation[-1]+L2_norm_mach[-1]+L2_norm_mat_temp[-1])

  file_data.close()
#  del mat_temp, mach_nb, radiation, mat_density, x_coord

# PLOT L1 AND L2 NORMS
out_file = out_file_base
plot_error_norms(nb_cells, L1_norm_density, L2_norm_density, 'density')
plot_error_norms(nb_cells, L1_norm_radiation, L2_norm_radiation, 'radiation')
plot_error_norms(nb_cells, L1_norm_mach, L2_norm_mach, 'mach-number')
plot_error_norms(nb_cells, L1_norm_mat_temp, L2_norm_mat_temp, 'mat-temp')
plot_error_norms(nb_cells, L1_norm_total, L2_norm_total, 'total')

# PLOT MATERIAL DENSITY AND TEMPERATURE, RADIATION TEMPERATURE AND MACH NUMBER
#plot_solution(x_coord, mat_density, x_coord_exact, mat_density_exact, x_offset[-1], 'density')
plot_solution(x_coord, radiation, x_coord_exact, radiation_exact, x_offset[-1], 'radiation', nb_cells[-1])
plot_solution(x_coord, mach_nb, x_coord_exact, mach_nb_exact, x_offset[-1], 'mach-number', nb_cells[-1])
plot_solution(x_coord, mat_temp, x_coord_exact, mat_temp_exact, x_offset[-1], 'mat-temp', nb_cells[-1])
#plot_visc_coeff(x_coord, 'visc-coeff-nel-1000-points0.csv', var_index, nb_cells[-1])
#plot_visc_coeff(x_coord, 'visc-coeff-nel-1600-points0.csv', var_index, nb_cells[-1])

# save L1 norm values
file_name = out_file+'-l1-norm.txt'
datafile_id = open(file_name, 'w+')
datafile_id.write('nb_cells '+' l1_rho '+' l1_temp '+' l1_mach '+' l1_rad \n')
data = np.asarray([nb_cells, L1_norm_density, L1_norm_mat_temp, L1_norm_mach, L1_norm_radiation])
data = data.T
np.savetxt(datafile_id, data, fmt=['%d','%.10e','%.10e','%.10e','%.10e'])
datafile_id.close()

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
# save x_offset value
file_name = out_file+'-x-offset.txt'
datafile_id = open(file_name, 'w+')
datafile_id.write('nb_cells '+' x_offset \n')
data = np.asarray([nb_cells, x_offset])
data = data.T
np.savetxt(datafile_id, data, fmt=['%d','%.4e'])
datafile_id.close()