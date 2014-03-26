import numpy
import scipy.optimize
import scipy.integrate
import math
import matplotlib.pyplot

execfile('define_functions.py')
execfile('input_deck.py')

# some constants are redefined for brevity
Cp = 1. / (gamma - 1.)
M0 = initialMach
P0 = initialRadratio

beta_eq = M0 / C0
sigma_a_eq = sigma_a(M0, 1)
sigma_t_eq = sigma_t(M0, 1)

execfile('final_equilibrium_values.py')
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
