# This script aims at reading csv files and store the data for postprocessing in 1-D.

#### import python modules ####
import sys, os, math, csv
import numpy as np
import itertools as it
import matplotlib.pyplot as plt
from decimal import *
#### import python modules ####

#### define function ####
def plot_solution(x_num, y_num, x_anal, y_anal, variable, text_on=False, mass_diff=0):
  x_anal_offset = [float(x) for x in x_anal]
  x_num = [float(x) for x in x_num]
  plt.plot(x_num, y_num, '+-', markevery=5000, markersize=8, label=r'$new \ solution$', linewidth=2)
  plt.plot(x_anal_offset, y_anal, 'o-', markevery=5000, markersize=8, label=r'$old \ solution$', linewidth=1.2)
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
  fig_name=variable+'-plot.eps'
  plt.savefig(fig_name)
  print 'Saving plot using Matplotlib:', fig_name
  plt.clf()

#### define function ####

# READ EXACT SOLUTION
file_exact_list_0 = []
# x-coordinates
file_exact_list_0.append('new-semi-analytical-sol/x_data.txt')
#file_exact_list.append('data_x.dat')
x_coord_exact_0 = []
file_data_exact=open(file_exact_list_0[-1], 'r')
x_coord_exact_0[:] = [ line[:-1] for line in file_data_exact]
# material density
file_exact_list_0.append('new-semi-analytical-sol/Density_data.txt')
#file_exact_list_0.append('data_Density.dat')
mat_density_exact_0 = []
file_data_exact=open(file_exact_list_0[-1], 'r')
mat_density_exact_0[:] = [ line[:-1] for line in file_data_exact]
# radiation energy density
file_exact_list_0.append('new-semi-analytical-sol/Er_data.txt')
#file_exact_list_0.append('data_RED.dat')
radiation_exact_0 = []
file_data_exact=open(file_exact_list_0[-1], 'r')
radiation_exact_0[:] = [ line[:-1] for line in file_data_exact]
## mach number or fluid velocity
#file_exact_list_0.append('Mach_Data.txt')
##file_exact_list_0.append('data_Mach.dat')
#mach_nb_exact = []
#file_data_exact=open(file_exact_list_0[-1], 'r')
#mach_nb_exact[:] = [ line[:-1] for line in file_data_exact]
# material temperature
file_exact_list_0.append('new-semi-analytical-sol/Tm_data.txt')
#file_exact_list_0.append('data_Temp.dat')
mat_temp_exact_0 = []
file_data_exact=open(file_exact_list_0[-1], 'r')
mat_temp_exact_0[:] = [ line[:-1] for line in file_data_exact]
nb_nodes_exact = len(x_coord_exact_0)
nb_exact_files = len(file_exact_list_0)
# remove duplicate values in x_coord_exact and consistently in other lists
x_coord_exact_temp, idx = np.unique(np.asarray(x_coord_exact_0), return_index=True)
idx = sorted(idx)
x_coord_exact_0 = [x_coord_exact_0[i] for i in idx]
mat_density_exact_0 = [mat_density_exact_0[i] for i in idx]
radiation_exact_0 = [radiation_exact_0[i] for i in idx]
#mach_nb_exact = [mach_nb_exact[i] for i in idx]
mat_temp_exact_0 = [mat_temp_exact_0[i] for i in idx]
## normalize mach number
#mach_nb_exact = [ float(i)/float(mach_nb_exact[0]) for i in mach_nb_exact]

# READ EXACT SOLUTION
file_exact_list_1 = []
# x-coordinates
file_exact_list_1.append('old-semi-analytical-sol/data_x.dat')
x_coord_exact_1 = []
file_data_exact=open(file_exact_list_1[-1], 'r')
x_coord_exact_1[:] = [ line[:-1] for line in file_data_exact]
# material density
file_exact_list_1.append('old-semi-analytical-sol/data_Density.dat')
mat_density_exact_1 = []
file_data_exact=open(file_exact_list_1[-1], 'r')
mat_density_exact_1[:] = [ line[:-1] for line in file_data_exact]
# radiation energy density
file_exact_list_1.append('old-semi-analytical-sol/data_RED.dat')
radiation_exact_1 = []
file_data_exact=open(file_exact_list_1[-1], 'r')
radiation_exact_1[:] = [ line[:-1] for line in file_data_exact]
# material temperature
file_exact_list_1.append('old-semi-analytical-sol/data_Temp.dat')
mat_temp_exact_1 = []
file_data_exact=open(file_exact_list_1[-1], 'r')
mat_temp_exact_1[:] = [ line[:-1] for line in file_data_exact]

# PLOT MATERIAL DENSITY AND TEMPERATURE, RADIATION TEMPERATURE AND MACH NUMBER
plot_solution(x_coord_exact_0, mat_density_exact_0, x_coord_exact_1, mat_density_exact_1, 'density')
plot_solution(x_coord_exact_0, radiation_exact_0, x_coord_exact_1, radiation_exact_1, 'radiation')
#plot_solution(x_coord, mach_nb, x_coord_exact, mach_nb_exact, x_offset[-1], 'mach-number', nb_cells[-1])
plot_solution(x_coord_exact_0, mat_temp_exact_0, x_coord_exact_1, mat_temp_exact_1, 'mat-temp')
#plot_visc_coeff(x_coord, 'visc-coeff-nel-1000-points0.csv', var_index, nb_cells[-1])
#plot_visc_coeff(x_coord, 'visc-coeff-nel-1600-points0.csv', var_index, nb_cells[-1])

## PLOT X_OFFSET AND SAVE VALUES IN FILE
## plot
#plt.plot(nb_cells, x_offset, '-o', markersize=8, linewidth=2)
#plt.xlabel(r'$cells$', fontsize=20)
#plt.ylabel(r'$x_{offset}$', fontsize=20)
#plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
#plt.margins(0.1, 0.1)
#fig_name=out_file+'-x-offset.eps'
#plt.savefig(fig_name)
#plt.clf()
## save
#datafile_id = open('x_offset.txt', 'w+')
#datafile_id.write('nb_cells '+' x_offset \n')
#data = np.asarray([nb_cells, x_offset])
#data = data.T
#np.savetxt(datafile_id, data, fmt=['%d','%.4e'])
#datafile_id.close()