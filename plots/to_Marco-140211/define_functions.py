sigma_t = lambda M, fRED: sigt * C0
sigma_a = lambda M, fRED: siga * C0
sigma_s = lambda M, fRED: sigma_t(M, fRED) - sigma_a(M, fRED)


# define density functions
fnctnDensity_num = lambda M: 3. * M0**2 * (gamma * M**2 + 1.)
fnctnDensity_den = lambda M, fRED:  M**2 * (3. * (gamma * M0**2 + 1.) + gamma * P0 * (1. - 3. * fRED))
fnctnDensity = lambda M, fRED: fnctnDensity_num(M) / fnctnDensity_den(M, fRED)

rho_plus = lambda fE, T: ((M0**2 + 1. / gamma + P0 / 3. - P0 * fE) + numpy.sqrt((M0**2 + 1. / gamma + P0 / 3. - P0 * fE)**2 - 4. * M0**2 * T / gamma)) * gamma / 2. / T
rho_minus = lambda fE, T: ((M0**2 + 1. / gamma + P0 / 3. - P0 * fE) - numpy.sqrt((M0**2 + 1. / gamma + P0 / 3. - P0 * fE)**2 - 4. * M0**2 * T / gamma)) * gamma / 2. / T

# define material temperature function
fnctnTemp = lambda M, fRED: M0**2 / M**2 / fnctnDensity(M, fRED)**2

# define speed function
fnctnSpeed = lambda M, fRED: M0 / fnctnDensity(M, fRED)

# define the relativistic measure, beta
fnctnBeta = lambda M, fRED: fnctnSpeed(M, fRED) / C0

# define radiative flux function
fnctnradFlux = lambda M, fRED: - kappa(M, fRED) * slope_fRED(M, fRED) + fnctnBeta(M, fRED) * (sigma_a(M, fRED) / sigma_t(M, fRED) * (fRED + fnctnTemp(M, fRED)**4) + sigma_s(M, fRED) / sigma_t(M, fRED) * (1. / f_interp(M) + 1.) * fRED) / (1. - 2. * fnctnBeta(M, fRED)**2 * sigma_a(M, fRED) / sigma_t(M, fRED)) / fnctnDensity(M, fRED)

# in numpy.interp the Mach_left term must be an increasing list, thus the reversal of terms
fRED_interp_left = lambda M: numpy.interp(M, Mach_left1[::-1], fRED_left1[::-1])
fRED_interp_right = lambda M: numpy.interp(M, Mach_right1, fRED_right1)
fRED_interp = lambda M: numpy.interp(M, Mach[::-1], fRED[::-1])
RED_interp = lambda M: numpy.interp(M, Mach[::-1], RED[::-1])
M_interp = lambda x_in: numpy.interp(x_in, x, Mach)
M_old_interp = lambda x_in: numpy.interp(x_in, x_old, Mach_old)
rF_interp = lambda M: numpy.interp(M, Mach[::-1], radFlux[::-1])
Speed_interp = lambda M: numpy.interp(M, Mach[::-1], Speed[::-1])
f_interp = lambda M: numpy.interp(M, Mach_RT[::-1], f[::-1])
f_old_interp = lambda M: numpy.interp(M, Mach_old[::-1], f_old[::-1])
fnum_interp = lambda M: numpy.interp(M, Mach_RT[::-1], f_num[::-1])
fden_interp = lambda M: numpy.interp(M, Mach_RT[::-1], f_den[::-1])

kappa = lambda M, fRED: C0**2 / sigma_t(M, fRED)

# beta_eq = M0 / C0
# sigma_a_eq = sigma_a(M0, 1)
# sigma_t_eq = sigma_t(M0, 1)

# define radTemp slope function
slope_fRED1 = lambda M, fRED: Cp * (fnctnTemp(M, fRED) - 1)
slope_fRED2 = lambda M, fRED: M0**2 * (1. - fnctnDensity(M, fRED)**2) / fnctnDensity(M, fRED)**2 / 2.
slope_fRED3 = lambda M, fRED: P0 * ((sigma_a(M, fRED) / sigma_t(M, fRED) * (fRED + fnctnTemp(M, fRED)**4) + sigma_s(M, fRED) / sigma_t(M, fRED) * (1. / f_interp(M) + 1.) * fRED) / fnctnDensity(M, fRED) / (1. - 2. * fnctnBeta(M, fRED)**2 * sigma_a(M, fRED) / sigma_t(M, fRED)) - 4. / 3. / (1. - 2. * beta_eq**2 * sigma_a_eq / sigma_t_eq))
slope_fRED = lambda M, fRED: M0 * (slope_fRED1(M, fRED) + slope_fRED2(M, fRED) + slope_fRED3(M, fRED)) / P0 / kappa(M, fRED)

# define Mach slope function
slope_Mach = lambda M, fRED:  M * P0 * (gamma + 1.) * (M0 * slope_fRED(M, fRED) + (gamma - 1.) / (gamma + 1.) * (gamma * M**2 + 1.) * fnctnDensity(M, fRED) * sigma_a(M, fRED) * (fRED / f_interp(M) - fnctnTemp(M, fRED)**4)) / (- 2. * fnctnDensity(M, fRED) * fnctnTemp(M, fRED) * M0 * (M**2 - 1.))

# define the radiative flux
rad_flux = lambda left, right: fnctnradFlux(Mach_right[- 1 - right], fRED_right[- 1 - right]) - fnctnradFlux(Mach_left[- 1 - left], fRED_left[- 1 - left])

DfRED_DM = lambda fRED, M: slope_fRED(M, fRED) / slope_Mach(M, fRED)

dx_dM_left = lambda x, M: - 2. * fnctnDensity(M, fRED_interp_left(M)) * fnctnTemp(M, fRED_interp_left(M)) * M0 * (M**2 - 1.) / (M * P0 * (gamma + 1.) * (M0 * slope_fRED(M, fRED_interp_left(M)) + (gamma - 1.) / (gamma + 1.) * (gamma * M**2 + 1.) * fnctnDensity(M, fRED_interp_left(M)) * sigma_a(M, fRED_interp_left(M)) * (fRED_interp_left(M) / f_interp(M) - fnctnTemp(M, fRED_interp_left(M))**4)))
dx_dM_right = lambda x, M: - 2. * fnctnDensity(M, fRED_interp_right(M)) * fnctnTemp(M, fRED_interp_right(M)) * M0 * (M**2 - 1.) /(M * P0 * (gamma + 1.) * (M0 * slope_fRED(M, fRED_interp_right(M)) + (gamma - 1.) / (gamma + 1.) * (gamma * M**2 + 1.) * fnctnDensity(M, fRED_interp_right(M)) * sigma_a(M, fRED_interp_right(M)) * (fRED_interp_right(M) / f_interp(M) - fnctnTemp(M, fRED_interp_right(M))**4)))
