# import numpy
# import scipy.optimize
# import scipy.integrate
# import math
# import matplotlib.pyplot

# execfile('input_deck.py')

# some constants are redefined for brevity
Cp = 1. / (gamma - 1.)
M0 = initialMach
P0 = initialRadratio

beta_eq = M0 / C0
sigma_a_eq = sigma_a(M0, 1)
sigma_t_eq = sigma_t(M0, 1)

execfile('final_equilibrium_values.py')
# execfile('define_functions.py')
Mach = numpy.append(Mach_left, Mach_right[::-1])
fRED = Mach
RED = 3. * fRED
Mach_RT = numpy.append(Mach_left, Mach_right[::-1])
f = Mach * 0 + 1./3.
x = 0 * f + 1.
f_old = f
f_iters = 0
max_iters = 100
f_err = 1.
RED_err = 1.
fRED_err = 1.
L2_f = 1.
L2_RED = 1.
L2_fRED = 1.
errs = numpy.zeros((max_iters + 1,9))
errs[0,:] = [f_iters, f_err, L2_f, RED_err, L2_RED, fRED_err, L2_fRED, 0, 0]

print 'going into RH'
log_file = open('collect_data.log', 'a+')
print >> log_file, 'going into RH'
log_file.close()
execfile('nonED_driver.py')
x_left2 = x_left1
x_fill_left = numpy.linspace(x_left2[0], x_left2[9], points, endpoint = True)
x_left3 = scipy.append(x_fill_left, x_left2[10:])
x_right2 = x_right1[::-1]
x_fill_right = numpy.linspace(x_right2[0], x_right2[9], points, endpoint = True)
x_right3 = scipy.append(x_fill_right, x_right2[10:])
x = scipy.append(x_left2, x_right2[::-1])
x_RT = scipy.append(x_left3, x_right3[::-1])

f = x_RT * 0 + 1./3.
Mach_RT = M_interp(x_RT)
Mach_RT[1991] = Mach[1000]
RED_RT = RED_interp(Mach_RT)
radTemp_RT = RED_RT**(1./4.)
fRED_RT = fRED_interp(Mach_RT)
Density_RT = fnctnDensity(Mach_RT, fRED_RT)
Speed_RT = fnctnSpeed(Mach_RT, fRED_RT)
radFlux_RT = fnctnradFlux(Mach_RT, fRED_RT)
Temp_RT = fnctnTemp(Mach_RT, fRED_RT)
RED_RT[len(x_RT) / 2. - 1] = (RED_RT[len(x_RT) / 2. - 2] + RED_RT[len(x_RT) / 2. + 1]) / 2.
RED_RT[len(x_RT) / 2.] = (RED_RT[len(x_RT) / 2. - 2] + RED_RT[len(x_RT) / 2. + 1]) / 2.
fRED_RT[len(x_RT) / 2. - 1] = (fRED_RT[len(x_RT) / 2. - 2] + fRED_RT[len(x_RT) / 2. + 1]) / 2.
fRED_RT[len(x_RT) / 2.] = (fRED_RT[len(x_RT) / 2. - 2] + fRED_RT[len(x_RT) / 2. + 1]) / 2.
radTemp_RT[len(x_RT) / 2. - 1] = (radTemp_RT[len(x_RT) / 2. - 2] + radTemp_RT[len(x_RT) / 2. + 1]) / 2.
radTemp_RT[len(x_RT) / 2.] = (radTemp_RT[len(x_RT) / 2. - 2] + radTemp_RT[len(x_RT) / 2. + 1]) / 2.
radFlux_RT[len(x_RT) / 2. - 1] = (radFlux_RT[len(x_RT) / 2. - 2] + radFlux_RT[len(x_RT) / 2. + 1]) / 2.
radFlux_RT[len(x_RT) / 2.] = (radFlux_RT[len(x_RT) / 2. - 2] + radFlux_RT[len(x_RT) / 2. + 1]) / 2.

x_diff = x_RT
Mach_diff = Mach_RT
RED_diff = RED_RT
radTemp_diff = radTemp_RT
fRED_diff = fRED_RT
Density_diff = Density_RT
Speed_diff = Speed_RT
radFlux_diff = radFlux_RT
Temp_diff = Temp_RT

f = 0 * x_RT + 1./3.
f_new = f
print 'going into RT'
log_file = open('collect_data.log', 'a+')
print >> log_file, 'going into RT'
log_file.close()
execfile('Sn_driver.py')
print 'f[:10] =', f[:10]
print 'f[-10:] =', f[-10:]
log_file = open('collect_data.log', 'a+')
print >> log_file, 'f[:10] =', f[:10]
print >> log_file, 'f[-10:] =', f[-10:]
log_file.close()
f_old = f_new
f_new = f
f_err = max(abs(f_new - f_old) / f)
L2_f = numpy.sqrt(sum(((f_new[2:-1] + f_new[1:-2]) / 2. - (f_old[2:-1] + f_old[1:-2]) / 2.)**2 * (x_RT[2:-1] - x_RT[1:-2])) / (x_RT[-2] - x_RT[1]))
f_err_low = f_err
L2_f_low = L2_f
RED_err = max(abs(RED_RT - 2. * numpy.pi * f_den) / RED_RT)
L2_RED = numpy.sqrt(sum(((RED_RT[2:-1] + RED_RT[1:-2]) / 2. - (2. * numpy.pi * f_den[2:-1] + 2. * numpy.pi * f_den[1:-2]) / 2.)**2 * (x_RT[2:-1] - x_RT[1:-2])) / (x_RT[-2] - x_RT[1]))
L2_RED_low = L2_RED
RED_err_low = RED_err
fRED_err = max(abs(fRED_RT - 2. * numpy.pi * f_num) / fRED_RT)
L2_fRED = numpy.sqrt(sum(((fRED_RT[2:-1] + fRED_RT[1:-2]) / 2. - (2. * numpy.pi * f_num[2:-1] + 2. * numpy.pi * f_num[1:-2]) / 2.)**2 * (x_RT[2:-1] - x_RT[1:-2])) / (x_RT[-2] - x_RT[1]))
L2_fRED_low = L2_fRED
fRED_err_low = fRED_err
hold_f = 1
errs[f_iters,:] = [f_iters, f_err, L2_f, RED_err, L2_RED, fRED_err, L2_fRED, int_deltaM_left, int_deltaM_right]
print errs[:f_iters + 1,:]
log_file = open('collect_data.log', 'a+')
print >> log_file, 'errs =', errs[:f_iters + 1,:]
log_file.close()
keep_running = 1

d_x = numpy.zeros((max_iters + 1, len(x_RT)))
d_M = numpy.zeros((max_iters + 1, len(x_RT)))
d_f = numpy.zeros((max_iters + 1, len(x_RT)))
d_Erh = numpy.zeros((max_iters + 1, len(x_RT)))
d_fErh = numpy.zeros((max_iters + 1, len(x_RT)))
d_Ert = numpy.zeros((max_iters + 1, len(x_RT)))
d_fErt = numpy.zeros((max_iters + 1, len(x_RT)))

d_x[f_iters,:] = x_RT
d_M[f_iters,:] = Mach_RT
d_f[f_iters,:] = f
d_Erh[f_iters,:] = RED_RT
d_fErh[f_iters,:] = fRED_RT
d_Ert[f_iters,:] = 2. * numpy.pi * f_den
d_fErt[f_iters,:] = 2. * numpy.pi * f_num

while (f_iters < max_iters) & (int(f_err > f_tol) | int(RED_err > RED_tol) | int(fRED_err > fRED_tol)) & (keep_running == 1):
  f_iters += 1
  print 'going into RH'
  log_file = open('collect_data.log', 'a+')
  print >> log_file, 'going into RH'
  log_file.close()
  execfile('nonED_driver_1.py')
  if (sum(numpy.diff(x_left1) < 0) == 0):
    x_left2 = x_left1
    x_fill_left = numpy.linspace(x_left2[0], x_left2[9], points, endpoint = True)
    x_left3 = scipy.append(x_fill_left, x_left2[10:])
  if (sum(numpy.diff(x_right1) < 0) == 0):
    x_right2 = x_right1[::-1]
    x_fill_right = numpy.linspace(x_right2[0], x_right2[9], points, endpoint = True)
    x_right3 = scipy.append(x_fill_right, x_right2[10:])
  x = scipy.append(x_left2, x_right2[::-1])
  x_old = x_RT
  x_RT = scipy.append(x_left3, x_right3[::-1])
  print 'x[1] =', x[0]
  print 'x[-2] =', x[-1]
  log_file = open('collect_data.log', 'a+')
  print >> log_file, 'x[1] =', x[1]
  print >> log_file, 'x[-2] =', x[-2]
  log_file.close()
  if ((sum(numpy.diff(x_left1) < 0) == len(x_left1) - 1) | (sum(numpy.diff(x_right1) < 0) == len(x_right1) - 1)):
    log_file = open('collect_data.log', 'a+')
    print >> log_file, 'Are the x-directions bad?'
    log_file.close()
    print 'Are the x-directions bad?'
    x_left1_bad = x_left1
    x_right1_bad = x_right1
    break
  Mach_old = Mach_RT
  Mach_RT = M_interp(x_RT)
  Mach_RT[1991] = Mach[1000]
  RED_RT = RED_interp(Mach_RT)
  fRED_RT = fRED_interp(Mach_RT)
  radTemp_RT = RED_RT**(1./4.)
  Density_RT = fnctnDensity(Mach_RT, fRED_RT)
  Speed_RT = fnctnSpeed(Mach_RT, fRED_RT)
  radFlux_RT = fnctnradFlux(Mach_RT, fRED_RT)
  Temp_RT = fnctnTemp(Mach_RT, fRED_RT)
  RED_RT[len(x_RT) / 2. - 1] = (RED_RT[len(x_RT) / 2. - 2] + RED_RT[len(x_RT) / 2. + 1]) / 2.
  RED_RT[len(x_RT) / 2.] = (RED_RT[len(x_RT) / 2. - 2] + RED_RT[len(x_RT) / 2. + 1]) / 2.
  fRED_RT[len(x_RT) / 2. - 1] = (fRED_RT[len(x_RT) / 2. - 2] + fRED_RT[len(x_RT) / 2. + 1]) / 2.
  fRED_RT[len(x_RT) / 2.] = (fRED_RT[len(x_RT) / 2. - 2] + fRED_RT[len(x_RT) / 2. + 1]) / 2.
  radTemp_RT[len(x_RT) / 2. - 1] = (radTemp_RT[len(x_RT) / 2. - 2] + radTemp_RT[len(x_RT) / 2. + 1]) / 2.
  radTemp_RT[len(x_RT) / 2.] = (radTemp_RT[len(x_RT) / 2. - 2] + radTemp_RT[len(x_RT) / 2. + 1]) / 2.
  radFlux_RT[len(x_RT) / 2. - 1] = (radFlux_RT[len(x_RT) / 2. - 2] + radFlux_RT[len(x_RT) / 2. + 1]) / 2.
  radFlux_RT[len(x_RT) / 2.] = (radFlux_RT[len(x_RT) / 2. - 2] + radFlux_RT[len(x_RT) / 2. + 1]) / 2.
  f_old = f_new
  print 'going into RT'
  log_file = open('collect_data.log', 'a+')
  print >> log_file, 'going into RT'
  log_file.close()
  execfile('Sn_driver.py')
  print 'f[:10] =', f[:10]
  print 'f[-10:] =', f[-10:]
  log_file = open('collect_data.log', 'a+')
  print >> log_file, 'f[:10] =', f[:10]
  print >> log_file, 'f[-10:] =', f[-10:]
  log_file.close()
  f_old = f_new
  f_new = f
  f_push = f_old_interp(M_old_interp(x_RT))
  f_err = max(abs((f - f_push) / f))
  L2_f = numpy.sqrt(sum(((f_push[2:-1] + f_push[1:-2]) / 2. - (f[2:-1] + f[1:-2]) / 2.)**2 * (x_RT[2:-1] - x_RT[1:-2])) / (x_RT[-2] - x_RT[1]))
  RED_err = max(abs(RED_RT - 2. * numpy.pi * f_den) / RED_RT)
  L2_RED = numpy.sqrt(sum(((RED_RT[2:-1] + RED_RT[1:-2]) / 2. - (2. * numpy.pi * f_den[2:-1] + 2. * numpy.pi * f_den[1:-2]) / 2.)**2 * (x_RT[2:-1] - x_RT[1:-2])) / (x_RT[-2] - x_RT[1]))
  fRED_err = max(abs(fRED_RT - 2. * numpy.pi * f_num) / fRED_RT)
  L2_fRED = numpy.sqrt(sum(((fRED_RT[2:-1] + fRED_RT[1:-2]) / 2. - (2. * numpy.pi * f_num[2:-1] + 2. * numpy.pi * f_num[1:-2]) / 2.)**2 * (x_RT[2:-1] - x_RT[1:-2])) / (x_RT[-2] - x_RT[1]))
  if (int(f_err < f_err_low)):
    f_err_low = f_err
  if (int(L2_f < L2_f_low) & int(L2_f > 0)):
    L2_f_low = L2_f
  if (int(RED_err < RED_err_low) & int(RED_err > 0)):
    RED_err_low = RED_err
  if (int(L2_RED < L2_RED_low) & int(L2_RED > 0)):
    L2_RED_low = L2_RED
  if (int(fRED_err < fRED_err_low) & int(fRED_err > 0)):
    fRED_err_low = fRED_err
  if (int(L2_fRED < L2_fRED_low) & int(L2_fRED > 0)):
    L2_fRED_low = L2_fRED
  errs[f_iters,:] = [f_iters, f_err, L2_f, RED_err, L2_RED, fRED_err, L2_fRED, int_deltaM_left, int_deltaM_right]
  d_x[f_iters,:] = x_RT
  d_M[f_iters,:] = Mach_RT
  d_f[f_iters,:] = f
  d_Erh[f_iters,:] = RED_RT
  d_fErh[f_iters,:] = fRED_RT
  d_Ert[f_iters,:] = 2. * numpy.pi * f_den
  d_fErt[f_iters,:] = 2. * numpy.pi * f_num
  if f_iters < 5:
    print errs[:f_iters + 1,:]
    log_file = open('collect_data.log', 'a+')
    print >> log_file, 'errs =', errs[:f_iters + 1,:]
    log_file.close()
  else:
    print errs[f_iters - 4:f_iters + 1,:]
    print f_err_low, L2_f_low, RED_err_low, L2_RED_low, fRED_err_low, L2_fRED_low
    log_file = open('collect_data.log', 'a+')
    print >> log_file, 'errs =', errs[f_iters - 4:f_iters + 1,:]
    print >> log_file, 'lows =', f_err_low, L2_f_low, RED_err_low, L2_RED_low, fRED_err_low, L2_fRED_low
    log_file.close()
  if int(f_err > f_err_low):
    keep_running = 0
  if (int(L2_f == L2_f_low) | int(RED_err == RED_err_low) | int(L2_RED == L2_RED_low) | int(fRED_err == fRED_err_low) | int(L2_fRED == L2_fRED_low)):
    keep_running = 1
  if (int(L2_f > L2_f_low) & (int(RED_err > RED_err_low) | int(RED_err < RED_tol)) & int(L2_RED > L2_RED_low) & (int(fRED_err > fRED_err_low) | int(fRED_err < fRED_tol)) & int(L2_fRED > L2_fRED_low)):
    keep_running = 0
  if int(f_err < f_tol):
    keep_running = 0

print 'RT complete'
log_file = open('collect_data.log', 'a+')
print >> log_file, 'RT complete'
log_file.close()

i = 0
name_list = ['x', 'M', 'f', 'REDrh', 'fREDrh', 'REDrt', 'fREDrt', 'errs']
for write_var in [d_x, d_M, d_f, d_Erh, d_fErh, d_Ert, d_fErt, errs]:
  open_file = open('data_' + name_list[i] + '.dat', 'w')
  numpy.savetxt(open_file, write_var[:f_iters + 1,:])
  open_file.close()
  i += 1

open_file = open('data_Im_fwd.dat', 'w')
numpy.savetxt(open_file, Im_fwd)
open_file.close()

open_file = open('data_Im_rev.dat', 'w')
numpy.savetxt(open_file, Im_rev)
open_file.close()

# Mach_RT_hold = Mach_RT
# 
# #################################################
# Mach_left = numpy.linspace(M0, 1, points, endpoint = False)
# fRED_left = numpy.zeros(2)
# fRED_left[0] = initialTemperature**4 / 3.
# Mach_right = numpy.linspace(finalMach, 1, points, endpoint = False)
# fRED_right = numpy.zeros(2)
# fRED_right[0] = finalTemperature**4 / 3.
# Mach = numpy.append(Mach_left, Mach_right[::-1])
# fRED = Mach
# RED = 3. * fRED
# Mach_RT = numpy.append(Mach_left, Mach_right[::-1])
# f = Mach * 0 + 1./3.
# x = 0 * f + 1.
# print 'going into RH'
# execfile('nonED_driver.py')
# # x_left1 = scipy.append(x0, x_left1)
# # x_right1 = scipy.append(x1, x_right1)
# # x = scipy.append(x_left1, x_right1[::-1])
# x = scipy.append(x_left1, x_right1)
# Mach_RT = Mach_RT_hold
# f = f_new
# #################################################

fig = matplotlib.pyplot.figure()
ax1 = fig.add_subplot(111)
matplotlib.pyplot.plot(x_RT, Temp_RT, label = r'$T$')
matplotlib.pyplot.plot(x_diff, Temp_diff, '--', label = r'$T_{diff}$')
matplotlib.pyplot.plot(x_RT, radTemp_RT, label = r'$\theta_{RH}$')
matplotlib.pyplot.plot(x_RT, (2. * numpy.pi * f_den)**(1./4.), '--', label = r'$\theta_{RT}$')
matplotlib.pyplot.plot(x_diff, radTemp_diff, '--', label = r'$\theta_{diff}$')
params = {'legend.fontsize':16}
matplotlib.pyplot.rcParams.update(params)
matplotlib.pyplot.legend(loc = 'upper left', title = 'Temperatures')
ax2 = ax1.twinx()
matplotlib.pyplot.plot(x_RT, f, label = '$VEF$')
matplotlib.pyplot.plot(x_RT, f * 0 + 1./3., label = r'$f = \frac{1}{3}$')
matplotlib.pyplot.legend(title = 'Eddington factors')
# matplotlib.pyplot.show()
var_name = 'M' + str(M0) + '.png'
matplotlib.pyplot.savefig(var_name)
