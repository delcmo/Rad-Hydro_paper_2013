#!python

import radshock
from utils import AnalyticShockProfiles
import matplotlib.pyplot as plt

def plot_solution(x_anal, y_anal, variable):
  plt.plot(x_anal, y_anal, 'o-', markevery=500, markersize=8, label=r'$exact \ solution$', linewidth=1.2)
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
#  plt.xlim(-0.015, 0.015)
  plt.margins(0.1, 0.1)
  fig_name=variable+'-plot.eps'
  plt.savefig(fig_name)
  print 'Saving plot using Matplotlib:', fig_name
  plt.clf()

## for mach 1.05 (the one simulated with Rhea)
#prob = radshock.greyNED_RadShock(M0 = 1.05, sigA = 577.3502692, sigS = 0., expDensity_abs = 0., expDensity_scat = 0., expTemp_abs = 0., expTemp_scat = 0., rho0 = 1., Tref = 100.)

## for mach 3 (the one simulated with Rhea)
#prob = radshock.greyNED_RadShock(M0 = 3.0, sigA = 577.3502692, sigS = 0., expDensity_abs = 0., expDensity_scat = 0., expTemp_abs = 0., expTemp_scat = 0., rho0 = 1., Tref = 100.)

## for mach 3 with temp-dep-opacity (the one simulated with Rhea)
#prob = radshock.greyNED_RadShock(M0 = 3.0, sigA = 7216.875, sigS = 0., expDensity_abs = 1., expDensity_scat = 0., expTemp_abs = -3.5, expTemp_scat = 0., rho0 = 1., Tref = 100., dumpnum = 100)

prob = radshock.greyNED_RadShock(M0 = 3.0, sigA = 577.3502692, sigS = 0., expDensity_abs = 1., expDensity_scat = 0., expTemp_abs = -3.5, expTemp_scat = 0., rho0 = 1., Tref = 100., dumpnum = 100)

# others
#prob = radshock.greyNED_RadShock(M0 = 3.0, sigA = 577.3502692, sigS = 0., expDensity_abs = 0., expDensity_scat = 0., expTemp_abs = 0., expTemp_scat = 0., rho0 = 5.516636721, Tref = 186.298821)

#prob = radshock.greyNED_RadShock(M0 = 3.0, sigA = 577.3502692, sigS = 0., expDensity_abs = 0., expDensity_scat = 0., expTemp_abs = 0., expTemp_scat = 0., rho0 = 1., Tref = 179.9539188)

#prob = radshock.greyNED_RadShock(M0 = 3.0, sigA = 7216.875, sigS = 0., expDensity_abs = 1., expDensity_scat = 0., expTemp_abs = -3.5, expTemp_scat = 0., rho0 = 5.516636721, Tref = 186.298821, dumpnum = 100)

#test = AnalyticShockProfiles(prob)
#test.write_data()

file_exact_list=[]
# x-coordinates
file_exact_list.append('x_data.txt')
x_coord_exact = []
file_data_exact=open(file_exact_list[-1], 'r')
x_coord_exact[:] = [ line[:-1] for line in file_data_exact]
# material density
file_exact_list.append('Density_data.txt')
mat_density_exact = []
file_data_exact=open(file_exact_list[-1], 'r')
mat_density_exact[:] = [ line[:-1] for line in file_data_exact]
# radiation energy density
file_exact_list.append('Er_data.txt')
radiation_exact = []
file_data_exact=open(file_exact_list[-1], 'r')
radiation_exact[:] = [ line[:-1] for line in file_data_exact]
# mach number or fluid velocity
file_exact_list.append('Mach_Data.txt')
mach_nb_exact = []
file_data_exact=open(file_exact_list[-1], 'r')
mach_nb_exact[:] = [ line[:-1] for line in file_data_exact]
# material temperature
file_exact_list.append('Tm_data.txt')
mat_temp_exact = []
file_data_exact=open(file_exact_list[-1], 'r')
mat_temp_exact[:] = [ line[:-1] for line in file_data_exact]

plot_solution(x_coord_exact, mat_density_exact, 'density')
plot_solution(x_coord_exact, radiation_exact, 'radiation')
plot_solution(x_coord_exact, mach_nb_exact, 'mach-number')
plot_solution(x_coord_exact, mat_temp_exact, 'mat-temp')