Mach_fwd = Mach_RT
Mach_rev = Mach_RT[::-1]

x_fwd = x_RT
x_rev = x_RT[::-1]

RED_fwd = RED_RT
RED_rev = RED_RT[::-1]

Temp_fwd = Temp_RT
Temp_rev = Temp_RT[::-1]

Speed_fwd = Speed_RT
Speed_rev = Speed_RT[::-1]

radFlux_fwd = radFlux_RT
radFlux_rev = radFlux_RT[::-1]

Im = numpy.zeros((Sn, len(x_RT)))
Im_fwd = numpy.zeros((Sn / 2, len(x_RT)))
Im_rev = numpy.zeros((Sn / 2, len(x_RT)))
Im_fwd[:,0] = RED_fwd[0] / 4. / numpy.pi * (1 + 3 * mus[0:Sn/2] * Speed_fwd[0] / C0)
Im_rev[:,0] = RED_rev[0] / 4. / numpy.pi * (1 + 3 * mus[Sn/2:Sn] * Speed_rev[0] / C0)
Im[:,0] = RED_RT[-1] / 4. / numpy.pi * (1 + 3 * mus * Speed_rev[0] / C0)

# break the transport equation into parts
# first, no material motion corrections
# the zeroth order terms
# (these terms are angle independent)
trans1_0 = lambda Im, M: - sigma_t(M, fRED_interp(M)) * Im
trans2_0 = lambda M: sigma_s(M, fRED_interp(M)) * RED_interp(M) / 4. / numpy.pi
trans3_0 = lambda M: sigma_a(M, fRED_interp(M)) * fnctnTemp(M, fRED_interp(M))**4 / 4. / numpy.pi
transport_0 = lambda Im, M: trans1_0(Im, M) + trans2_0(M) + trans3_0(M)
# second, the material motion corrections
# the first order terms
# (the first term is angle independent, the other three are angle dependent)
transport1_A = lambda M: - 2. * sigma_s(M, fRED_interp(M)) * rF_interp(M) / 4. / numpy.pi
trans2_1 = lambda Im, M: sigma_t(M, fRED_interp(M)) * Im
trans3_1 = lambda M: 3. * sigma_s(M, fRED_interp(M)) * RED_interp(M) / 4. / numpy.pi
trans4_1 = lambda M: 3. * sigma_a(M, fRED_interp(M)) * fnctnTemp(M, fRED_interp(M))**4 / 4. / numpy.pi
transport1_B = lambda Im, M: trans2_1(Im, M) + trans3_1(M) + trans4_1(M)
# third, the energy and momentum conservation corrections
# (the first term is angle independent, the second term is angle dependent)
transport2_A = lambda M: Speed_interp(M)**2 / C0**2 * sigma_t(M, fRED_interp(M)) * (RED_interp(M) + fRED_interp(M)) / 4. / numpy.pi
transport2_B = lambda M: 3. * (2. * Speed_interp(M)**2 / C0**2 * sigma_a(M, fRED_interp(M)) * rF_interp(M)) / 4. / numpy.pi

k_m = 0

# traverse positive mu directions
for mu in mus[:len(mus)/2]:
  print 'k_m = ', k_m
  log_file = open('collect_data.log', 'a+')
  print >> log_file, 'k_m =', k_m
  log_file.close()
  # the transport equation to be solved with mu_m on the RHS
  DIm_Dx = lambda Im, x_in: transport_0(Im, M_interp(x_in)) / mu / C0 + Speed_interp(M_interp(x_in)) / C0 * (transport1_A(M_interp(x_in)) / mu + transport1_B(Im, M_interp(x_in))) / C0 + transport2_A(M_interp(x_in)) / mu / C0 + transport2_B(M_interp(x_in)) / C0
  a = scipy.integrate.odeint(DIm_Dx, Im_fwd[k_m, 0], x_fwd, mxstep = 5000, atol = int_atol, rtol = int_rtol)
  Im_fwd[k_m, :] = numpy.transpose(a)
  k_m += 1

k_m = 0

# traverse negative mu directions
for mu in mus[len(mus)/2:]:
  print 'k_m = ', k_m + Sn / 2 
  log_file = open('collect_data.log', 'a+')
  print >> log_file, 'k_m =', k_m + Sn / 2
  log_file.close()
  # the transport equation to be solved with mu_m on the RHS
  DIm_Dx = lambda Im, x_in: transport_0(Im, M_interp(x_in)) / mu / C0 + Speed_interp(M_interp(x_in)) / C0 * (transport1_A(M_interp(x_in)) / mu + transport1_B(Im, M_interp(x_in))) / C0 + transport2_A(M_interp(x_in)) / mu / C0 + transport2_B(M_interp(x_in)) / C0
  a = scipy.integrate.odeint(DIm_Dx, Im_rev[k_m, 0], x_rev, mxstep = 5000, atol = int_atol, rtol = int_rtol)
  Im_rev[k_m, :] = numpy.transpose(a)
  k_m += 1

f_num = 0 * Im_fwd[0,:]
f_den = 0 * Im_fwd[0,:]

for i in range(Sn/2):
  f_num += mus[i]**2 * weights[i] * Im_fwd[i,:]
  f_den += weights[i] * Im_fwd[i,:]

for i in range(Sn/2, Sn):
  f_num += mus[i]**2 * weights[i] * Im_rev[i - Sn / 2,::-1]
  f_den += weights[i] * Im_rev[i - Sn / 2,::-1]

f = f_num / f_den

Im_fwd = numpy.zeros((Sn / 2, len(x_RT)))
Im_rev = numpy.zeros((Sn / 2, len(x_RT)))
Im_fwd_init = RED_fwd[0] / 4. / numpy.pi * (1 + 3 * mus[0:Sn/2] * Speed_fwd[0] / C0)
Im_rev_init = RED_rev[0] / 4. / numpy.pi * (1 + 3 * mus[Sn/2:Sn] * Speed_rev[0] / C0)

def deriv_fwd(Im, x):
  d = 4. * numpy.pi
  M_i = M_interp(x)
  fRED_i = fRED_interp(M_i)
  sum_Im = 2. * numpy.pi * (weights[0] * Im[0] + weights[1] * Im[1]) + RED_interp(M_i) / 2.
  F = rF_interp(M_i)
  T = fnctnTemp(M_i, fRED_i)
  t = sigma_t(M_i, fRED_i)
  s = sigma_s(M_i, fRED_i)
  a = sigma_a(M_i, fRED_i)
  beta = Speed_interp(M_i) / C0
  return numpy.array([
    ( - t * Im[0] + s * sum_Im / d + a * T**4 / d - 2. * s * beta * F / d + mus[0] * beta * (t * Im[0] + 3. * s * sum_Im / d + 3. * a * T**4 / d)) / mus[0],
    ( - t * Im[1] + s * sum_Im / d + a * T**4 / d - 2. * s * beta * F / d + mus[1] * beta * (t * Im[1] + 3. * s * sum_Im / d + 3. * a * T**4 / d)) / mus[1]
  ])

def deriv_rev(Im, x):
  d = 4. * numpy.pi
  M_i = M_interp(x)
  fRED_i = fRED_interp(M_i)
  sum_Im = 2. * numpy.pi * (weights[2] * Im[0] + weights[3] * Im[1]) + RED_interp(M_i) / 2.
  F = rF_interp(M_i)
  T = fnctnTemp(M_i, fRED_i)
  t = sigma_t(M_i, fRED_i)
  s = sigma_s(M_i, fRED_i)
  a = sigma_a(M_i, fRED_i)
  beta = Speed_interp(M_i) / C0
  return numpy.array([
    ( - t * Im[0] + s * sum_Im / d + a * T**4 / d - 2. * s * beta * F / d + mus[2] * beta * (t * Im[0] + 3. * s * sum_Im / d + 3. * a * T**4 / d)) / mus[2],
    ( - t * Im[1] + s * sum_Im / d + a * T**4 / d - 2. * s * beta * F / d + mus[3] * beta * (t * Im[1] + 3. * s * sum_Im / d + 3. * a * T**4 / d)) / mus[3]
  ])

