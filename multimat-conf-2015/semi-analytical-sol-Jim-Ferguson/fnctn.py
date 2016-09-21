import numpy
import scipy.integrate

def sigma_a(P, M, self):
    val = self.sigA * mat_density(P, M, self)**self.expDensity_abs \
        * mat_temp(P, M, self)**self.expTemp_abs
    return val

def sigma_s(P, M, self):
    val = self.sigS * mat_density(P, M, self)**self.expDensity_scat \
        * mat_temp(P, M, self)**self.expTemp_scat
    return val

def sigma_t(P, M, self):
    val = sigma_a(P, M, self) + sigma_s(P, M, self)
    return val

# checked 160321_1208
def C_0(P, M, self):
  if (self.problem == 'JMF'):
    val = mat_beta(P, M, self)**2 \
        * ((sigma_s(P, M, self) - sigma_a(P, M, self)) \
           * (rad_energy_density(P, M, self) + P) \
           + sigma_a(P, M, self) \
           * (mat_temp(P, M, self)**4 - rad_energy_density(P, M, self)))
  elif (self.problem == 'LE'):
    val = 0.
  return val

# checked 160614_1324
def Qeq_Sre(P, M, self):
  if (self.problem == 'JMF'):
    if (M > 1.):
        P = self.Pr0
        M = self.M0
    else:
        P = self.Pr1
        M = self.M1
    val  = 4. * mat_beta(P, M, self)**2 * P
    val *= (sigma_s(P, M, self) - sigma_a(P, M, self))
  elif (self.problem == 'LE'):
    val = 0.
  return val

# checked 160523_1504
def Sre(P, M, self):
    if (self.problem == 'JMF'):
      val = sigma_a(P, M, self) \
          * (mat_temp(P, M, self)**4 - rad_energy_density(P, M, self)) \
          + mat_beta(P, M, self) \
          * (sigma_a(P, M, self) - sigma_s(P, M, self)) * rad_flux(P, M, self) \
          + Qeq_Sre(P, M, self)
#           + C_0(P, M, self)
    elif (self.problem == 'LE'):
      val = sigma_a(P, M, self) \
          * (mat_temp(P, M, self)**4 - rad_energy_density(P, M, self)) \
          + mat_beta(P, M, self) * dPdx(P, M, self)
    return val

# checked 160321_1208
def Srp(P, M, self):
    if (self.problem == 'JMF'):
      val = - sigma_t(P, M, self) * rad_flux(P, M, self) \
          + mat_beta(P, M, self) \
          * (sigma_t(P, M, self) * P \
             + sigma_s(P, M, self) * rad_energy_density(P, M, self) \
             + sigma_a(P, M, self) * mat_temp(P, M, self)**4)
    elif (self.problem == 'LE'):
      val = - sigma_t(P, M, self) * rad_flux(P, M, self)
    return val

# checked 160321_1208
def Srie(P, M, self):
    if (self.problem == 'JMF'):
      val = Sre(P, M, self) - mat_beta(P, M, self) * Srp(P, M, self)
    elif (self.problem == 'LE'):
      val = sigma_a(P, M, self) \
          * (mat_temp(P, M, self)**4 - rad_energy_density(P, M, self))
    return val

# checked 160321_1238
def mat_density(P, M, self):
    val = self.M0**2 * (self.gamma * M**2 + 1.) / M**2 \
        / (self.gamma * self.M0**2 + 1. + self.gamma * self.P0 * (self.Pr0 - P))
    return val

# checked 160321_1238
def mat_temp(P, M, self):
    val = self.M0**2 / M**2 / mat_density(P, M, self)**2
    return val

# checked 160321_1238
def mat_speed(P, M, self):
    val = self.M0 / mat_density(P, M, self)
    return val

# checked 160321_1238
def mat_beta(P, M, self):
    val = mat_speed(P, M, self) / self.C0
    return val

# checked 160321_1238
def mat_pres(P, M, self):
    val = mat_density(P, M, self) * mat_temp(P, M, self) / self.gamma
    return val

# checked 160321_1238
def mat_internal_energy(P, M, self):
    val = mat_temp(P, M, self) / self.gamma / (self.gamma - 1.)
    return val

# checked 160321_1238
def mat_total_energy(P, M, self):
    val = 0.5 * mat_density(P, M, self) * mat_speed(P, M, self)**2 \
        + mat_density(P, M, self) * mat_internal_energy(P, M, self) \
        + mat_pres(P, M, self)
    return val    

def momentum_flux(P, M, self):
    val  = mat_density(P, M, self) * mat_speed(P, M, self)**2
    val += mat_pres(P, M, self)
    return val

def energy_flux(P, M, self):
    val = mat_beta(P, M, self) * mat_total_energy(P, M, self)
    return val

# checked 160321_1238
def rad_energy_density(P, M, self):
    val = P / f_interp(P, M, self)
    return val

def rad_temp(P, M, self):
    val = rad_energy_density(P, M, self)**(1./4.)
    return val

# checked 160321_1158
def rad_flux2(P, M, self):
    if (self.problem == 'JMF'):
      val = (sigma_t(P, M, self) * P \
             + sigma_s(P, M, self) * rad_energy_density(P, M, self) \
             + sigma_a(P, M, self) * mat_temp(P, M, self)**4) \
          / sigma_t(P, M, self)
    elif (self.problem == 'LE'):
      val = 4. * P
    return val

# checked 160321_1158
def rad_flux(P, M, self):
    val = - dPdx(P, M, self) / sigma_t(P, M, self) \
        + mat_beta(P, M, self) * rad_flux2(P, M, self)
    return val

# checked 160321_1158
def dPdx(P, M, self):
    P0 = self.P0
    if (numpy.size(M) > 1):
        Pr0 = self.Pr0
        M0 = self.M0
    elif (M > 1):
        Pr0 = self.Pr0
        M0 = self.M0
    else:
        Pr0 = self.Pr1
        M0 = self.M1
    val = sigma_t(P, M, self) / P0 \
        * (mat_beta(P, M, self) \
           * (mat_total_energy(P, M, self) + P0 * rad_flux2(P, M, self)) \
         - mat_beta(Pr0, M0, self)
           * (mat_total_energy(Pr0, M0, self) + P0 * rad_flux2(Pr0, M0, self)))
    return val

# checked 160321_1319
def drhodx(P, M, self):
    val = self.P0 * (self.M0 * dPdx(P, M, self) \
                     - (self.gamma - 1.) * mat_density(P, M, self) \
                       * self.C0 * Srie(P, M, self)) \
        / self.M0 / mat_temp(P, M, self) / (M**2 - 1.)
    return val

# checked 160321_1323
def dTdx(P, M, self):
    val = (mat_temp(P, M, self) * (self.gamma * M**2 - 1) * drhodx(P, M, self) \
        - self.gamma * self.P0 * dPdx(P, M, self)) / mat_density(P, M, self)
    return val

def dMdx(P, M, self):
    val = - M * (drhodx(P, M, self) / mat_density(P, M, self) \
                 + dTdx(P, M, self) / 2. / mat_temp(P, M, self))
    return val

def dPdM(P, M, self):
    val = dPdx(P, M, self) / dMdx(P, M, self)
    return val

def dxdM(x, M, self):
    if (M > 1.):
        val = numpy.interp(M, \
                           self.Mach_precursor[::-1], \
                           self.Pr_precursor[::-1])
    else:
        val = numpy.interp(M, \
                           self.Mach_relaxation, \
                           self.Pr_relaxation)
    return (1. / dMdx(val, M, self))

def ddM_ode(M, vals, self):
    P, x = vals
    val_dPdM = dPdM(P, M, self)
    val_dxdM = 1. / dMdx(P, M, self)
    return [val_dPdM, val_dxdM]

# break the transport equation into parts
# first, no material motion corrections: transx_0
# the zeroth order terms are independent of angle
def trans1_0(Im, M, self):
    P = fRED_interp_RT(M, self)
    return - sigma_t(P, M, self) * Im

def trans2_0(M, self):
    P = fRED_interp_RT(M, self)
    return sigma_s(P, M, self) * rad_energy_density(P, M, self) / 4. / numpy.pi

def trans3_0(M, self):
    P = fRED_interp_RT(M, self)
    return sigma_a(P, M, self) * mat_temp(P, M, self)**4 / 4. / numpy.pi

def transport_0(Im, M, self):
    return trans1_0(Im, M, self) + trans2_0(M, self) + trans3_0(M, self)

# second, the material motion corrections
# the first order terms
# (the first term is angle independent, the other three are angle dependent)
def transport1_A(M, self):
    P = fRED_interp_RT(M, self)
    return - 2. * sigma_s(P, M, self) * rad_flux(P, M, self) / 4. / numpy.pi

def trans2_1(Im, M, self):
    P = fRED_interp_RT(M, self)
    return sigma_t(P, M, self) * Im

def trans3_1(M, self):
    P = fRED_interp_RT(M, self)
    val  = 3. * sigma_s(P, M, self) * rad_energy_density(P, M, self)
    val /= 4. * numpy.pi
    return val

def trans4_1(M, self):
    P = fRED_interp_RT(M, self)
    val  = 3. * sigma_a(P, M, self) * mat_temp(P, M, self)**4
    val /= 4. * numpy.pi
    return val

def transport1_B(Im, M, self):
    return trans2_1(Im, M, self) + trans3_1(M, self) + trans4_1(M, self)

# third, the second-order correction term
# which imposes that the radiation transport equation should be 0 in equilibrium
# there are two terms:
# the first is angle independent and the second is quadratic in angle
def transport2_0_EQ(M, self):
    if (M > 1.):
        P = self.Pr0
        M = self.M0
    else:
        P = self.Pr1
        M = self.M1
# because we are subtracting this term,
# the term below is positive instead of negative
    val = 2. * mat_beta(P, M, self)**2 * P / numpy.pi * sigma_s(P, M, self)
    return val

def transport2_2_EQ(M, self):
    if (M > 1.):
        P = self.Pr0
        M = self.M0
    else:
        P = self.Pr1
        M = self.M1
# because we are subtracting this term,
# the term below is negative instead of positive
    val = - 3. * mat_beta(P, M, self)**2 * P / numpy.pi * sigma_t(P, M, self)
    return val

def DIm_Dx_ode(x_in, Im, self):
    M = M_interp_RT(x_in, self)
    P = fRED_interp_RT(M, self)
    val = transport_0(Im, M, self) / self.mu + mat_beta(P, M, self) * (transport1_A(M, self) / self.mu + transport1_B(Im, M, self)) + transport2_0_EQ(M, self) / self.mu + transport2_2_EQ(M, self) * self.mu
    return val

def fRED_interp_RT(M, self):
    return numpy.interp(M, self.Mach_RT[::-1], self.Pr_RT[::-1])

def P_RT_interp(M, self):
    return numpy.interp(M, self.Mach_RT[::-1], self.P_RT[::-1])

def M_interp_RT(x_in, self):
    return numpy.interp(x_in, self.x_RT, self.Mach_RT)

def f_interp(P, M, self):
    if hasattr(self, 'f'):
#         P = numpy.interp(M, self.d_M[-1][::-1], self.d_Prh[-1][::-1])
#         return numpy.interp(P, self.d_Prt[-1], self.d_f[-1])
        return numpy.interp(M, self.d_M[-1][::-1], self.d_f[-1][::-1])
    else:
        return 1. / 3.

def Q_mu_x(mu, x, self): 
    return trans2_0(M_interp_RT(x, self), self) + trans3_0(M_interp_RT(x, self), self) + mat_beta(fRED_interp_RT(M_interp_RT(x, self), self), M_interp_RT(x, self), self) * transport1_A(M_interp_RT(x, self), self) + mu * mat_beta(M_interp_RT(x, self), fRED_interp_RT(M_interp_RT(x, self), self), self) * (trans3_1(M_interp_RT(x, self), self) + trans4_1(M_interp_RT(x, self), self)) + transport2(M_interp_RT(x, self), self)

def Ieq_x(mu, self):
    if (mu >= 0.):
        val  = self.Ert[-1][0] / 4. / numpy.pi
        val *= (1. + 4. * mu * mat_beta(self.Pr_RT[0], self.Mach_RT[0], self))
    else:
        val  = self.Ert[-1][0] / 4. / numpy.pi
        val *= (1. + 4. * mu * mat_beta(self.Pr_RT[0], self.Mach_RT[0], self))
    return val 

def exp_factor(mu, x, self): 
    val  = 1. - mat_beta(fRED_interp_RT(M_interp_RT(x, self), self), M_interp_RT(x, self), self) * mu
    val *= sigma_t(fRED_interp_RT(M_interp_RT(x, self), self), M_interp_RT(x, self), self) * x / mu
    return val

def I_mu_x(mus, x, self): 
    vals = numpy.zeros(numpy.size(mus))
    for i in range(numpy.size(vals)):
        print 'i = ', i
        if (mus[i] >= 0.):
            vals[i] = Ieq_x(mus[i], self) * numpy.exp(exp_factor(mus[i], self.x_RT[0], self)) * numpy.exp(- exp_factor(mus[i], x, self)) + numpy.exp(- exp_factor(mus[i], x, self)) * scipy.integrate.quad(lambda x_var: Q_mu_x(mus[i], x_var, self) / mus[i] * numpy.exp(exp_factor(mus[i], x_var, self)), self.x_RT[0], x)[0]
        else:
            vals[i] = Ieq_x(mus[i], self) * numpy.exp(exp_factor(mus[i], self.x_RT[-1], self)) * numpy.exp(- exp_factor(mus[i], x, self)) + numpy.exp(- exp_factor(mus[i], x, self)) * scipy.integrate.quad(lambda x_var: Q_mu_x(mus[i], x_var, self) / mus[i] * numpy.exp(exp_factor(mus[i], x_var, self)), self.x_RT[-1], x)[0]
    return vals
