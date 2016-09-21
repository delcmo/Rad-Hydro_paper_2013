"""
Because this is the first item in the file, this is a docstring.
The (!) difference is that 'help' on this file will grab this docstring, and return it to the terminal (useful!).  The comments below go into overhead, but otherwise are unaccessible (useless!).  Use docstrings to explain how to use file/function/class/method...'  Use comments for longterm maintenance to describe internal methods of code.
"""

import os
# import AMHC
import copy
import numpy
import fnctn
import pickle
# import rage_utils
import scipy.optimize
import scipy.integrate
import matplotlib.pyplot

# Hi, I'm a comment.

class AnalyticShockProfiles:
    '''
    '''
    def __init__(self, incoming):
        print 'entered AnalyticShockProfiles __init__ for M0 = ', incoming.M0
        print '\n'
        self.M0 = incoming.M0
        self.P0 = incoming.P0
        self.C0 = incoming.C0
        self.sigA = incoming.sigA
        self.sigS = incoming.sigS
        self.expDensity_abs = incoming.expDensity_abs
        self.expTemp_abs = incoming.expTemp_abs
        self.expDensity_scat = incoming.expDensity_scat
        self.expTemp_scat = incoming.expTemp_scat
        self.gamma = incoming.gamma
        self.problem = incoming.problem
        self.T0 = incoming.T0
        self.Pr0 = self.T0**4 / 3.
        self.num_pts = 1e3 # 1e4 # 5e3 # 1e3
        self.left_pts = self.num_pts * self.M0 + 1
        self.right_pts = 2. * self.num_pts + 1
        self.max_steps = 1e4
        self.int_atol = incoming.int_tol
        self.int_rtol = incoming.int_tol
        self.full_output = 1
        self.eps_precursor_equil  = incoming.eps
        self.eps_precursor_ASP    = incoming.eps
        self.eps_relaxation_equil = incoming.eps
        self.eps_relaxation_ASP   = incoming.eps
        self.integrator = incoming.integrator
        self.Pr_init = numpy.zeros(2)
        self.Mach_precursor_baseline = 0.
        self.Mach_relaxation_baseline = 0.
        self.left1 = 0
        self.right1 = 0
        self.x0_mult = 1.
        print 'leaving AnalyticShockProfiles __init__ for M0 = ', incoming.M0
        print '\n'

    def downstream_equilibrium(self):
        print 'entered downstream_equilibrium'
        print '\n'
        M0 = self.M0
        P0 = self.P0
        gamma = self.gamma
        def a2(T):
            return T / gamma
        def a1(T):
            return P0 * (T**4 - 1.) / 3. - M0**2 - 1. / gamma
        def a0(T):
            return M0**2
        def b2(T):
            return (T - 1.) / (gamma - 1.) - M0**2 / 2. - 4. * P0 / 3. 
        def b1(T):
            return 4. * P0 * T**4 / 3. 
        def b0(T):
            return M0**2 / 2. 
        def discriminant_a(T):
            return a1(T)**2 - 4. * a2(T) * a0(T)
        def discriminant_b(T):
            return b1(T)**2 - 4. * b2(T) * b0(T)
        zero_discriminant_a = scipy.optimize.fsolve(discriminant_a, 1., \
                                                    xtol = 1.e-13)
        zero_discriminant_b = scipy.optimize.fsolve(discriminant_b, 1., \
                                                    xtol = 1.e-13)
        def momentum_and_energy(x):
            rho, T   = x
            momentum = M0**2 + rho**2 * T / gamma + P0 * rho * T**4 / 3. \
                     - rho * (M0**2 + 1. / gamma + P0 / 3.)
            energy   = M0**2 / 2. + rho**2 * T / (gamma - 1.) \
                     + 4. * P0 * rho * T**4 / 3. \
                     - rho**2 * (M0**2 / 2. + 1. / (gamma - 1.) + 4. * P0 / 3.)
            return [momentum, energy]
        self.rho1, self.T1 = scipy.optimize.fsolve(momentum_and_energy, \
                                                  [zero_discriminant_b, \
                                                   zero_discriminant_a], \
                                                   xtol = 1.e-13)
        self.M1 = M0 / self.rho1 / numpy.sqrt(self.T1)
        self.Pr1 = self.T1**4 / 3.
        self.speed1 = M0 / self.rho1
        print 'leaving downstream_equilibrium'
        print '\n'
        
    def print_RH_sources(self):
        print 'entered print_sources'
        print "\n"
        val  = fnctn.mat_temp(self.Pr0, self.M0, self)**4
        val -= fnctn.rad_energy_density(self.Pr0, self.M0, self)
        print_str  = 'fnctn.mat_temp(self.Pr0, self.M0, self)**4 - '
        print_str += 'fnctn.rad_energy_density(self.Pr0, self.M0, self) = '
        print_str += str(val)
        print print_str
        print '\n'
        val  = fnctn.sigma_a(self.Pr0, self.M0, self)
        val -= fnctn.sigma_s(self.Pr0, self.M0, self)
        val *= fnctn.rad_flux(self.Pr0, self.M0, self)
        val *= fnctn.mat_beta(self.Pr0, self.M0, self)
        val += fnctn.Qeq_Sre(self.Pr0, self.M0, self)
        print_str  = 'fnctn.rad_flux(self.Pr0, self.M0, self) + '
        print_str += 'fnctn.Qeq_Sre(self.Pr0, self.M0, self) = '
        print_str += str(val)
        print print_str
        print '\n'
        print "fnctn.Sre(self.Pr0, self.M0, self) =", \
               fnctn.Sre(self.Pr0, self.M0, self)
        print "\n"
        print "fnctn.Srp(self.Pr0, self.M0, self) =", \
               fnctn.Srp(self.Pr0, self.M0, self)
        print "\n"
        print "fnctn.Srie(self.Pr0, self.M0, self) =", \
               fnctn.Srie(self.Pr0, self.M0, self)
        print "\n"
        print "fnctn.dPdx(self.Pr0, self.M0, self) =", \
               fnctn.dPdx(self.Pr0, self.M0, self)
        print "\n"
        print "fnctn.dMdx(self.Pr0, self.M0, self) =", \
               fnctn.dMdx(self.Pr0, self.M0, self)
        print "\n"
        print "fnctn.dPdM(self.Pr0, self.M0, self) =", \
               fnctn.dPdM(self.Pr0, self.M0, self)
        print "\n"
        print 'self.T1**4                                        = ', self.T1**4
        val = fnctn.mat_temp(self.Pr1, self.M1, self)**4 
        print 'fnctn.mat_temp(self.Pr1, self.M1, self)**4        = ', val
        val = fnctn.rad_energy_density(self.Pr1, self.M1, self)
        print 'fnctn.rad_energy_density(self.Pr1, self.M1, self) = ', val
        val  = fnctn.mat_temp(self.Pr1, self.M1, self)**4
        val -= fnctn.rad_energy_density(self.Pr1, self.M1, self)
        print '\n'
        print_str  = 'fnctn.mat_temp(self.Pr1, self.M1, self)**4 - '
        print_str += 'fnctn.rad_energy_density(self.Pr1, self.M1, self) = '
        print_str += str(val)
        print print_str
        print '\n'
        val  = fnctn.sigma_a(self.Pr1, self.M1, self)
        val -= fnctn.sigma_s(self.Pr1, self.M1, self)
        val *= fnctn.rad_flux(self.Pr1, self.M1, self)
        val *= fnctn.mat_beta(self.Pr1, self.M1, self)
        val += fnctn.Qeq_Sre(self.Pr1, self.M1, self)
        print_str  = 'fnctn.rad_flux(self.Pr1, self.M1, self) + '
        print_str += 'fnctn.Qeq_Sre(self.Pr1, self.M1, self) = '
        print_str += str(val)
        print print_str
        print '\n'
        print "fnctn.Sre(self.Pr1, self.M1, self) =", \
               fnctn.Sre(self.Pr1, self.M1, self)
        print "\n"
        print "fnctn.Srp(self.Pr1, self.M1, self) =", \
               fnctn.Srp(self.Pr1, self.M1, self)
        print "\n"
        print "fnctn.Srie(self.Pr1, self.M1, self) =", \
               fnctn.Srie(self.Pr1, self.M1, self)
        print "\n"
        print "fnctn.dPdx(self.Pr1, self.M1, self) =", \
               fnctn.dPdx(self.Pr1, self.M1, self)
        print "\n"
        print "fnctn.dMdx(self.Pr1, self.M1, self) =", \
               fnctn.dMdx(self.Pr1, self.M1, self)
        print "\n"
        print "fnctn.dPdM(self.Pr1, self.M1, self) =", \
               fnctn.dPdM(self.Pr1, self.M1, self)
        print '\n'
        print 'leaving print_sources'
        print '\n'

    def make_mach_arrays(self):
        print 'entered make_mach_arrays'
        print '\n'
        self.precursor_M_ends = numpy.array([self.M0 \
                                             - self.eps_precursor_equil, \
                                             1. + self.eps_precursor_ASP])
        self.Mach_precursor = numpy.linspace(self.precursor_M_ends[0], \
                              self.precursor_M_ends[1], \
                              self.left_pts, \
                              endpoint = True)
        self.relaxation_M_ends = numpy.array([self.M1 \
                                              + self.eps_relaxation_equil, \
                                              1. - self.eps_relaxation_ASP])
        self.Mach_relaxation = numpy.linspace(self.relaxation_M_ends[0], \
                               self.relaxation_M_ends[1], \
                               self.right_pts, \
                               endpoint = True)
        if (hasattr(self, 'f')):
            self.add_to_mach_arrays()
        print 'leaving make_mach_arrays'
        print '\n'

    def add_to_mach_arrays(self):
        print 'entered add_to_mach_arrays'
        print '\n'
        iP = numpy.argmin(numpy.abs(self.d_M[-1] - 1))
        if (self.d_M[-1][iP] < 1.):
            iP -= 1
        iS = iP + 1
        for M in self.Machs_to_add:
            print 'inserting points around M = ', M
            print '\n'
            iM = numpy.argmin(numpy.abs(self.d_M[-1] - M))
            print 'iP = ', iP
            print 'iS = ', iS
            print 'iM = ', iM
            if (iM <= iP):
                iMin = int(iM - self.left_pts / 10)
                iMax = min(int(iM + self.left_pts / 10), iP)
            elif (iM >= iS):
                iMin = max(int(iM - self.right_pts / 10), iS)
                iMax = int(iM + self.right_pts / 10)
            machs_to_add = numpy.linspace(self.d_M[-1][iMin], \
                                          self.d_M[-1][iMax], \
                                          max(int(self.num_pts), 1001))
            print 'machs_to_add = ', machs_to_add
            print 'numpy.min(machs_to_add) = ', numpy.min(machs_to_add)
            print 'numpy.max(machs_to_add) = ', numpy.max(machs_to_add)
            if (M > 1):
                print 'self.Mach_precursor = ', self.Mach_precursor
                machs_to_add = scipy.append(self.Mach_precursor, machs_to_add)
                machs_to_add = numpy.sort(machs_to_add)
                self.Mach_precursor = machs_to_add[::-1]
                print 'self.Mach_precursor = ', self.Mach_precursor
            else:
                print 'self.Mach_relaxation = ', self.Mach_relaxation
                machs_to_add = scipy.append(self.Mach_relaxation, machs_to_add)
                print 'machs_to_add = ', machs_to_add
                machs_to_add = numpy.sort(machs_to_add)
                print 'machs_to_add = ', machs_to_add
                self.Mach_relaxation = machs_to_add
                print 'self.Mach_relaxation = ', self.Mach_relaxation
        print 'leaving add_to_mach_arrays'
        print '\n'

    def linearize_away_from_equilibrium(self):
        print 'entered linearize_away_from_equilibrium'
        print '\n'
        if hasattr(self, 'f'):
            self.Pr_init[0] = fnctn.P_RT_interp(self.Mach_precursor[0], self)
            self.Pr_init[1] = fnctn.P_RT_interp(self.Mach_relaxation[0], self)
        else:
            eps = self.eps_precursor_equil
            init_val = 1. / 3. + eps
            print 'init_val = ', init_val
            M = self.precursor_M_ends[0]
            val = scipy.optimize.fsolve( \
                  lambda P: P - 1./3. + eps * fnctn.dPdM(P, M, self), \
                  init_val, xtol = 1.e-13)[0]
            print 'val = ', val
            if (val < 1./3.):
                print "if (val < 1./3.):"
                val += 2. * numpy.abs(1./3. - val)
            self.Pr_init[0] = val
            eps = self.eps_relaxation_equil
            init_val = self.T1**4 / 3. - eps
            M = self.relaxation_M_ends[0]
            val = scipy.optimize.fsolve(lambda P: P - self.T1**4 / 3. \
                                        - eps * fnctn.dPdM(P, M, self), \
                                        init_val, xtol = 1.e-13)[0]
            print 'val = ', val
            print 'self.T1**4 / 3. = ', self.T1**4 / 3.
            if (val > self.T1**4 / 3.):
                print "if (val > self.T1**4 / 3.):"
                val -= 2. * numpy.abs(self.T1**4 / 3. - val)
                print 'val = ', val
            self.Pr_init[1] = val
            print '\n'
        print 'leaving linearize_away_from_equilibrium'
        print '\n'

    def integrate_ddM(self):
        print 'entered integrate_ddM'
        print '\n'
        mult_exp = 0
        mult_val = 1.
        for side in ('precursor', 'relaxation'):
            if (side == 'precursor'):
                init_vals = [self.Pr_init[0], 0]
                integration_array = self.Mach_precursor
                print 'self.Mach_precursor = ', self.Mach_precursor
                print 'self.Pr_init[0] = ', self.Pr_init[0]
                print '\n'
            elif (side == 'relaxation'):
                init_vals = [self.Pr_init[1], 0]
                integration_array = self.Mach_relaxation
                print 'self.Mach_relaxation = ', self.Mach_relaxation
                print 'self.Pr_init[1] = ', self.Pr_init[1]
                print '\n'
            while (mult_exp < 5):
                print 'integrating ' + side + ' region'
                print '\n'
                result_arrays = numpy.zeros((numpy.size(integration_array), 2))
                i = scipy.integrate.ode(fnctn.ddM_ode)
                i.set_integrator('vode', atol = self.int_atol, \
                                 rtol = self.int_rtol, method = 'bdf', \
                                 nsteps = self.max_steps)
                i.set_initial_value(init_vals, integration_array[0])
                result_arrays[0,:] = init_vals
                i.set_f_params(self)
                k = 1
                while (i.successful() & (k < numpy.size(integration_array))):
                    dt = numpy.diff(integration_array[k - 1:k+1])
                    i.integrate(i.t + dt)
                    result_arrays[k, 0] = i.y[0]
                    result_arrays[k, 1] = i.y[1]
                    k += 1
                P = result_arrays[:, 0]
                x = result_arrays[:, 1]
                got_zeros = numpy.any(P == 0.)
                print 'got_zeros = ', got_zeros
                print '\n'
                got_spurious_solution = False
                if ((side == 'relaxation') & (hasattr(self, 'f_iters'))):
                    bool1 = numpy.any(P > \
                                      fnctn.mat_temp(P, integration_array, self)**4 * numpy.interp(integration_array, self.d_M[-1][::-1], self.d_f[-1][::-1])) 
                    print 'bool1 = ', bool1
                    bool2 = numpy.any(numpy.diff(P) > 0)
                    print 'bool2 = ', bool2
                    got_spurious_solution = (bool1 | bool2)
                print 'got_spurious_solution = ', got_spurious_solution
                print '\n'
                if (got_zeros & (not got_spurious_solution)):
                    while got_zeros:
                        P = result_arrays[:, 0]
                        x = result_arrays[:, 1]
                        if (side == 'precursor'):
                            self.Mach_precursor  = numpy.delete(self.Mach_precursor,  numpy.argmin(P))
                        else:
                            self.Mach_relaxation  = numpy.delete(self.Mach_relaxation,  numpy.argmin(P))
                        x = numpy.delete(x, numpy.argmin(P))
                        P = numpy.delete(P, numpy.argmin(P))
                        M = numpy.delete(integration_array, numpy.argmin(P))[-1]
                        result_arrays = numpy.array([P,x]).transpose()
                        got_zeros = numpy.any(result_arrays[:,0] == 0.)
                    if (side == 'precursor'):
                        self.Mach_precursor_baseline = M - 1.
                    elif (side == 'relaxation'):
                        self.Mach_relaxation_baseline = 1. - M
                if (got_spurious_solution):
                    print 'got a spurious solution'
                    print 'spurious solution = ', result_arrays[:, 0]
                    print '\n'
                    if (mult_exp == 0):
                        mult_exp += 1
                        mult_val -= 10**(- mult_exp)
                        print 'mult_exp = ', mult_exp
                        print 'mult_val = ', mult_val
                        init_vals[0] = mult_val * self.Pr_init[1]
                    else:
                        mult_val -= 10**(- mult_exp)
                        print 'mult_val = ', mult_val
                        init_vals[0] = mult_val * self.Pr_init[1]
                    print '\n'
                elif (mult_exp != 0):
                    mult_val += 10**(- mult_exp)
                    mult_exp += 1
                    mult_val -= 10**(- mult_exp)
                    print 'mult_exp = ', mult_exp
                    print 'mult_val = ', mult_val
                    print '\n'
                    init_vals[0] = mult_val * self.Pr_init[1]
                else:
                    break
            if (side == 'precursor'):
                self.Pr_precursor = result_arrays[:,0]
                self.x_precursor = result_arrays[:,1] - result_arrays[-1,1]
                if (self.integrator == 'odeint'):
                    self.precursor_report = result_report
            else:
                self.Pr_relaxation = result_arrays[:,0]
                self.x_relaxation = result_arrays[:,1] - result_arrays[-1,1]
                if (self.integrator == 'odeint'):
                    self.relaxation_report = result_report
        mult_val = 1. - (1. - mult_val) / 4. / numpy.pi
        self.mult_val = mult_val
        print 'self.Mach_precursor = ', self.Mach_precursor
        print 'self.Mach_relaxation = ', self.Mach_relaxation
        print '\n'
        print 'self.Pr_precursor  = ', self.Pr_precursor
        print 'self.Pr_relaxation = ', self.Pr_relaxation
        print '\n'
        print 'leaving integrate_ddM'
        print '\n'

    def check_Pr_overlap(self):
        print 'entered check_Pr_overlap'
        print '\n'
        self.continuous_shock = 0
        if (self.Pr_precursor[-1] > self.Pr_relaxation[-1] + self.int_atol):
            Pcont = scipy.optimize.bisect( \
                    lambda P: \
                    fnctn.rad_flux(P, \
                                   numpy.interp(P, \
                                                self.Pr_precursor, \
                                                self.Mach_precursor), \
                                   self) \
                  - fnctn.rad_flux(P, \
                                   numpy.interp(P, \
                                                self.Pr_relaxation[::-1], \
                                                self.Mach_relaxation[::-1]), \
                                   self), \
                    self.Pr_relaxation[-1], \
                    self.Pr_precursor[-1])
        else:
            self.continuous_shock = 1
            print 'seeing a continuous shock'
            print '\n'
            self.integrate_ddM_again = 0
            Pcont = numpy.interp(1., \
                                      numpy.array([self.Mach_relaxation[-2], \
                                                   self.Mach_precursor[-2]]), \
                                      numpy.array([self.Pr_relaxation[-1], \
                                                   self.Pr_precursor[-1]]))
            self.Pr_precursor[-1] = Pcont
            self.Pr_relaxation[-1] = Pcont
        self.left1 = sum(self.Pr_precursor <= Pcont)
        self.right1 = sum(self.Pr_relaxation >= Pcont)
        print 'self.left1 = ', self.left1
        print 'self.right1 = ', self.right1
        if (self.left1 == self.left_pts):
            self.integrate_ddM_again = 0
        if (self.right1 == self.right_pts):
            self.integrate_ddM_again = 0
        if (self.continuous_shock == 1):
            self.integrate_ddM_again = 0
            print 'self.left1 = ', self.left1
            print 'self.right1 = ', self.right1
        self.Pcont = Pcont
        print 'self.Pcont = ', self.Pcont
        self.eps_precursor_ASP = self.Mach_precursor[:self.left1 + 2][-1] - 1.
        print 'self.eps_precursor_ASP = ', self.eps_precursor_ASP
        self.eps_relaxation_ASP = 1.- self.Mach_relaxation[:self.right1 + 2][-1]
        print 'self.eps_relaxation_ASP = ', self.eps_relaxation_ASP
        print '\n'
        print 'leaving check_Pr_overlap'
        print '\n'

    def splice_precursor_and_relaxation(self):
        print 'entered splice_precursor_and_relaxation'
        print '\n'
        Mach_precursor_last = numpy.interp(self.Pcont, \
                                           self.Pr_precursor, \
                                           self.Mach_precursor)
        Mach_relaxation_last = numpy.interp(self.Pcont, \
                                            self.Pr_relaxation[::-1], \
                                            self.Mach_relaxation[::-1])
        x_precursor_last = numpy.interp(self.Pcont, \
                                        self.Pr_precursor, \
                                        self.x_precursor)
        x_relaxation_last = numpy.interp(self.Pcont, \
                                         self.Pr_relaxation[::-1], \
                                         self.x_relaxation[::-1])
        self.Pr_precursor = self.Pr_precursor[:self.left1]
        self.Pr_relaxation = self.Pr_relaxation[:self.right1]
        self.Pr_precursor[-1] = self.Pcont
        self.Pr_relaxation[-1] = self.Pcont
        self.Mach_precursor = self.Mach_precursor[:self.left1]
        self.Mach_relaxation = self.Mach_relaxation[:self.right1]
        self.Mach_precursor[-1] = Mach_precursor_last
        self.Mach_relaxation[-1] = Mach_relaxation_last
        self.x_precursor = self.x_precursor[:self.left1] - x_precursor_last
        self.x_relaxation = self.x_relaxation[:self.right1] - x_relaxation_last
        self.x_precursor[-1] = 0.
        self.x_relaxation[-1] = 0.
        self.Mach = scipy.append(scipy.append(self.M0, self.Mach_precursor), \
                                 scipy.append(self.Mach_relaxation[::-1], \
                                              self.M1))
        self.Pr = scipy.append(scipy.append(self.Pr0, self.Pr_precursor), \
                               scipy.append(self.Pr_relaxation[::-1], \
                                            self.Pr1))
        x0  = 5. * self.x_precursor[0]
        x1  = max(numpy.abs(x0), 5. * self.x_relaxation[0])
        if (hasattr(self, 'f')):
            x0 = min(self.x0, x0)
            x1 = max(self.x1, x1)
        self.x = scipy.append(scipy.append(x0, self.x_precursor), \
                              scipy.append(self.x_relaxation[::-1], x1))
        self.Tm = fnctn.mat_temp(self.Pr, self.Mach, self)
        self.Er = fnctn.rad_energy_density(self.Pr, self.Mach, self)
        self.Tr = fnctn.rad_temp(self.Pr, self.Mach, self)
        self.Speed = fnctn.mat_speed(self.Pr, self.Mach, self)
        self.Density = fnctn.mat_density(self.Pr, self.Mach, self)
        self.Pressure = fnctn.mat_pres(self.Pr, self.Mach, self)
        self.Fr = fnctn.rad_flux(self.Pr, self.Mach, self)
        print 'self.x = ', self.x
        print '\n'
        multiplier = 1. / 2.
        if (hasattr(self, 'f_iters')):
            multiplier = (self.f_iters + 1.) / (self.f_iters + 2.)
        self.eps_precursor_ASP  = Mach_precursor_last - 1.
        self.eps_precursor_ASP -= self.Mach_precursor_baseline
        self.eps_precursor_ASP *= multiplier
        self.eps_precursor_ASP += self.Mach_precursor_baseline
        self.eps_relaxation_ASP  = 1. - Mach_relaxation_last
        self.eps_relaxation_ASP -= self.Mach_relaxation_baseline
        self.eps_relaxation_ASP *= multiplier
        self.eps_relaxation_ASP += self.Mach_relaxation_baseline
        print 'leaving splice_precursor_and_relaxation'
        print '\n'

    def fill_in_xs(self):
        print 'entered fill_in_xs'
        print '\n'
        x0 = self.x[0]
        x1 = self.x[-1]
        x0  = self.x[0] * self.M0**3
        x1  = self.M0**3 * self.x[-1]
        x0 *= self.x0_mult
        self.x0 = x0
        x1 *= self.x0_mult
        self.x1 = x1
        x_linfill_left = numpy.linspace(x0, self.x_precursor[9], \
                                        self.left_pts, endpoint = True)
        x_logfill_left = numpy.logspace(numpy.log10(numpy.abs(x0)), \
                                        numpy.log10(self.x_precursor[9] \
                                                    + 100. * numpy.abs(x0)), \
                                        self.left_pts, endpoint = True)
        x_logfill_left_shock = ((x_logfill_left - x_logfill_left[0]) \
                              / (x_logfill_left - x_logfill_left[0])[-1] \
                              * (self.x_precursor[-10] - self.x_precursor[-1]) \
                              + self.x_precursor[-1])[::-1]
        x_logfill_left_equil = ((x_logfill_left - x_logfill_left[0]) \
                              / (x_logfill_left - x_logfill_left[0])[-1] \
                              * (- x0 + self.x_precursor[9]) + x0)
        x_left_RT = scipy.append(\
                    scipy.append(\
                    scipy.append(x_logfill_left_equil, x_logfill_left_shock), \
                                 x_linfill_left), self.x_precursor)
        x_linfill_right = numpy.linspace(x1, self.x_relaxation[9], \
                                         5 * self.right_pts, endpoint = True)
        x_logfill_right = numpy.logspace(numpy.log10(x1), \
                                         numpy.log10(100. * x1), \
                                         5 * self.right_pts, \
                                         endpoint = True)[::-1]
        x_logfill_right_shock = ((x_logfill_right[-1] - x_logfill_right) \
                               / (x_logfill_right[-1] - x_logfill_right)[0] \
                               * (self.x_relaxation[-10] \
                                  - self.x_relaxation[-1]) \
                               + self.x_relaxation[-1])
        x_logfill_right_equil = ((x_logfill_right[-1] - x_logfill_right) \
                               / (x_logfill_right[-1] - x_logfill_right)[0] \
                               * (self.x_relaxation[9] - x1) + x1)[::-1]
        x_right_RT = scipy.append(\
                     scipy.append(\
                     scipy.append(x_logfill_right_equil, \
                                  x_logfill_right_shock), \
                                  x_linfill_right), self.x_relaxation)
        x_left_RT.sort()
        x_right_RT.sort()
        x_left_RT = x_left_RT[1:-1]
        x_right_RT = x_right_RT[1:-1]
        x_RT = scipy.append(x_left_RT, x_right_RT)
        x_RT.sort()
        self.x_RT = x_RT
        x_RT_zero = numpy.argmin(numpy.abs(x_RT))
        Mach_RT = numpy.interp(x_RT, self.x, self.Mach)
        Mach_RT[x_RT_zero] = self.Mach_precursor[-1]
        Mach_RT[x_RT_zero + 1] = self.Mach_relaxation[-1]
        self.Mach_RT = Mach_RT
        self.Pr_RT = numpy.interp(self.Mach_RT, self.Mach[::-1], self.Pr[::-1])
        print 'leaving fill_in_xs'
        print '\n'

    def make_Ims(self):
        print 'entered make_Ims'
        print '\n'
        if (self.f_iters == 0):
            self.mus = numpy.loadtxt("mus_" + str(self.Sn) + ".txt")
            self.weights = numpy.loadtxt("weights_" + str(self.Sn) + ".txt")
        self.Im_fwd = numpy.zeros((self.Sn / 2, numpy.size(self.x_RT)))
        self.Im_rev = numpy.zeros((self.Sn / 2, numpy.size(self.x_RT)))
        val  = self.mus[:self.Sn / 2] * fnctn.mat_beta(self.Pr0, self.M0, self)
        val *= 4.
        val += 1.
        val *= fnctn.mat_temp(self.Pr0, self.M0, self)**4 / 4. / numpy.pi
        self.Im_fwd[:,0] = val
        val  = self.mus[self.Sn / 2:] * fnctn.mat_beta(self.Pr1, self.M1, self)
        val *= 4.
        val += 1.
        val *= fnctn.mat_temp(self.Pr1, self.M1, self)**4 / 4. / numpy.pi
        self.Im_rev[:,-1] = val
        print 'leaving make_Ims'
        print '\n'

    def linearize_transport_equation(self):
        print 'entered linearize_transport_equation'
        print '\n'
        x = self.x_RT
        for i, self.mu in enumerate(self.mus[:self.Sn / 2]):
            print 'self.mu = ', self.mu
            eps = 1.e-4
            init_val = self.Im_fwd[i, 0]
            print 'init_val = ', init_val
            x1 = x[1]
            diff = x1 - x[0]
            val = scipy.optimize.fsolve( \
                  lambda I: I - init_val - diff*fnctn.DIm_Dx_ode(x1, I, self), \
                  init_val + eps, xtol = 1.e-13)[0]
            print 'val =      ', val
            self.Im_fwd[i, 1] = val
            print '\n'
        for i, self.mu in enumerate(self.mus[self.Sn / 2:]):
            print 'self.mu = ', self.mu
            eps = 1.e-4
            init_val = self.Im_rev[i, -1]
            print 'init_val = ', init_val
            x1 = x[-1]
            diff = x1 - x[-1]
            val = scipy.optimize.fsolve( \
                  lambda I: I - init_val - diff*fnctn.DIm_Dx_ode(x1, I, self), \
                  init_val - eps, xtol = 1.e-13)[0]
            print 'val =      ', val
            self.Im_rev[i, -2] = val
            print '\n'
        print 'leaving linearize_transport_equation'
        print '\n'

    def integrate_Sn(self):
        print 'entered integrate_Sn'
        print '\n'
        print 'Sn = ', self.Sn
        for direction in ['fwd', 'rev']:
            T0 = self.T0
            u0 = self.M0
            T1 = self.T1
            u1 = self.M0 / self.rho1
            if (direction == 'fwd'):
                mus = self.mus[:self.Sn / 2]
                init_vals = self.Im_fwd[:, 1]
                integration_array = self.x_RT[1:]
            elif (direction == 'rev'):
                mus = self.mus[self.Sn / 2:]
                init_vals = self.Im_rev[:, -2]
                integration_array = self.x_RT[::-1][1:]
            else:
                print 'bad direction passed into integrate_Sn'
            result_array = numpy.zeros(numpy.size(integration_array))
            for j, mu in enumerate(mus):
                self.mu = mu
                mu_i = (j + self.Sn / 2) if (direction == 'rev') else j
                print 'i = ' + str(mu_i) + ' and mu = ' + str(mu)
                i = scipy.integrate.ode(fnctn.DIm_Dx_ode)
                i.set_integrator('vode', atol = self.int_atol, rtol = self.int_rtol, \
                                 method = 'bdf', nsteps = self.max_steps)
                i.set_initial_value(init_vals[j], integration_array[0])
                result_array[0] = init_vals[j]
                i.set_f_params(self)
                k = 1
                while (i.successful() & (k < numpy.size(integration_array))):
                    dt = numpy.diff(integration_array)[k - 1]
                    i.integrate(i.t + dt)
                    result_array[k] = i.y[0]
                    k += 1
                val0  = (1. + 4. * self.mu * u0 / self.C0)
                val0 *= T0**4 / 4. / numpy.pi
                val1  = (1. + 4. * self.mu * u1 / self.C0)
                val1 *= T1**4 / 4. / numpy.pi
                if (direction == 'fwd'):
                    self.Im_fwd[j,1:] = result_array
                    print 'self.Im_fwd[j, :]  = ', self.Im_fwd[j, :]
                    print 'self.Im_fwd[j, 0]  = ', self.Im_fwd[j, 0]
                    print 'first elem should  = ', val0
                    absval0 = numpy.abs(val0 - self.Im_fwd[j, 0])
                    print 'first elem |diff|  = ', absval0
                    print 'self.Im_fwd[j, -2] = ', self.Im_fwd[j, -2]
                    print 'self.Im_fwd[j, -1] = ', self.Im_fwd[j, -1]
                    print 'last elem should   = ', val1
                    absval1 = numpy.abs(val1 - self.Im_fwd[j, -1])
                    print 'last elem |diff|   = ', absval1
                    if (absval1 > 1.e-9):
                        self.x0_mult = numpy.sqrt(absval1 / 1.e-9)
                    else:
                        self.x0_mult = 1.
                else:
                    self.Im_rev[j,:-1] = result_array[::-1]
                    print 'self.Im_rev[j, :]  = ', self.Im_rev[j, :]
                    print 'self.Im_rev[j, 1]  = ', self.Im_rev[j, 1]
                    print 'self.Im_rev[j, 0]  = ', self.Im_rev[j, 0]
                    print 'first elem should  = ', val0
                    absval0 = numpy.abs(val0 - self.Im_rev[j, 0])
                    print 'first elem |diff|  = ', absval0
                    print 'self.Im_rev[j, -1] = ', self.Im_rev[j, -1]
                    print 'last elem should   = ', val1
                    absval1 = numpy.abs(val1 - self.Im_rev[j, -1])
                    print 'last elem |diff|   = ', absval1
                print '\n'
        print 'leaving integrate_Sn'
        print '\n'

    def Sn_angular_moments(self):
        print 'entered Sn_angular_moments'
        print '\n'
        self.E_RT = numpy.zeros(numpy.size(self.Im_fwd[0,:]))
        self.F_RT = numpy.zeros(numpy.size(self.Im_fwd[0,:]))
        self.P_RT = numpy.zeros(numpy.size(self.Im_fwd[0,:]))
        for i in range(self.Sn / 2):
            self.E_RT += 2. * numpy.pi * self.weights[i] * self.Im_fwd[i,:]
            self.F_RT += 2. * numpy.pi * self.weights[i] * self.Im_fwd[i,:] \
                       * self.mus[i]
            self.P_RT += 2. * numpy.pi * self.weights[i] * self.Im_fwd[i,:] \
                       * self.mus[i]**2
        for i in range(self.Sn / 2, self.Sn):
            self.E_RT += 2. * numpy.pi * self.weights[i] \
                       * self.Im_rev[i - self.Sn / 2,:]
            self.F_RT += 2. * numpy.pi * self.weights[i] \
                       * self.Im_rev[i - self.Sn / 2,:] * self.mus[i]
            self.P_RT += 2. * numpy.pi * self.weights[i] \
                       * self.Im_rev[i - self.Sn / 2,:] * self.mus[i]**2
        self.f = self.P_RT / self.E_RT
        print 'upstream equilibrium radiation values:\n'
        print 'Er        = ', self.E_RT[0]
        print 'should    = ', self.T0**4
        val  = numpy.abs(self.E_RT[0] - self.T0**4)
        val /= min(self.E_RT[0], self.T0**4)
        print '|rel err| = ', val
        print 'Fr        = ', self.F_RT[0]
        Fval = 4./3. * self.M0 / self.C0 * self.T0**4
        print 'should    = ', Fval
        val  = numpy.abs(self.F_RT[0] - Fval) / min(self.F_RT[0], Fval)
        print '|rel err| = ', numpy.abs(val)
        print 'Pr        = ', self.P_RT[0]
        print 'should    = ', 1./3.
        val  = numpy.abs(self.P_RT[0] - self.Pr0) / min(self.P_RT[0], self.Pr0)
        print '|rel err| = ', val
        print 'f         = ', self.f[0]
        print 'should    = ', 1./3.
        val  = numpy.abs(self.f[0] - 1./3.) / min(self.f[0], 1./3.)
        print '|rel err| = ', val
        print '\n'
        print 'downstream equilibrium radiation values:\n'
        print 'Er        = ', self.E_RT[-1]
        print 'should    = ', self.T1**4 
        val  = numpy.abs(self.E_RT[-1] - self.T1**4)
        val /= min(self.E_RT[-1], self.T1**4)
        print '|rel err| = ', val
        print 'Fr        = ', self.F_RT[-1]
        Fval = 4./3. * self.M0 / self.rho1 / self.C0 * self.T1**4 
        print 'should    = ', Fval
        val  = numpy.abs(self.F_RT[-1] - Fval) / min(self.F_RT[-1], Fval)
        print '|rel err| = ', numpy.abs(val)
        print 'Pr        = ', self.P_RT[-1]
        print 'should    = ', self.T1**4 / 3.
        val = numpy.abs(self.P_RT[-1] - self.Pr1) / min(self.P_RT[-1], self.Pr1)
        print '|rel err| = ', val
        print 'f         = ', self.f[-1]
        print 'should    = ', 1./3.
        val  = numpy.abs(self.f[-1] - 1./3.) / min(self.f[-1], 1./3.)
        print '|rel err| = ', val
        print 'leaving Sn_angular_moments'
        print '\n'

    def make_RH_solution(self):
        print 'entered make_RH_solution'
        print '\n'
        self.integrate_ddM_again = 2
        if (self.print_the_sources):
#         Ensure equations are written correctly such that
#         Sre, Srp, Srie, dPdx and dMdx all equal 0
#         when evaluated at the equilibrium points
            self.print_RH_sources()
        while self.integrate_ddM_again:
            self.integrate_ddM_again -= 1
            self.make_mach_arrays()
            if (self.integrate_ddM_again == 1):
                self.linearize_away_from_equilibrium()
            self.integrate_ddM()
            self.check_Pr_overlap()
            if (self.continuous_shock == 1.):
                print 'in make_RH_solution, seeing a continuous shock'
                print 're-integrating with extended eps_ASPs'
                print '\n'
                self.eps_precursor_ASP /= 2.
                self.eps_relaxation_ASP /= 2.
                self.make_mach_arrays()
                self.integrate_ddM()
                self.check_Pr_overlap()
            if (self.integrate_ddM_again == 0):
                self.splice_precursor_and_relaxation()
        if (self.continuous_shock == 1.):
            if (not hasattr(self, 'f_iters')):
                os.system('touch continuous_shock_fiter0')
            else:
                os.system('touch continuous_shock_fiter' + str(self.f_iters))
        print 'leaving make_RH_solution'
        print '\n'

    def make_RT_solution(self):
        print 'entered make_RT_solution'
        print '\n'
        self.fill_in_xs()
        self.make_Ims()
        self.linearize_transport_equation()
        if (self.print_the_sources):
#         Ensure equations are written correctly such that
#         Sre, Srp, Srie, dPdx and dMdx all equal 0
#         when evaluated at the equilibrium points
            self.print_RT_sources()
        self.integrate_Sn()
        self.Sn_angular_moments()
        self.update_dictionaries()
        self.compute_errors()
        self.pickle_dictionary_variables()
        self.pickle_prob()
        print 'leaving make_RT_solution'
        print '\n'

    def print_RT_sources(self):
        print 'entered print_RT_sources'
        print '\n'
        if (bool(self.print_the_sources)):
            for i, self.mu in enumerate(self.mus[:self.Sn / 2]):
                print 'for self.mu = ', self.mu
                print 'fnctn.DIm_Dx_ode(self.x_RT[0], self.Im_fwd[i, 0], self) =', \
                       fnctn.DIm_Dx_ode(self.x_RT[0], self.Im_fwd[i, 0], self)
                print '\n'
            for i, mu in enumerate(self.mus[self.Sn / 2:]):
                self.mu = mu
                print 'for self.mu = ', self.mu
                print 'fnctn.DIm_Dx_ode(self.x_RT[-1], self.Im_rev[i, -1], self) =', \
                       fnctn.DIm_Dx_ode(self.x_RT[-1], self.Im_rev[i, -1], self)
                print '\n'
        print 'leaving print_RT_sources'
        print '\n'

    def make_dictionaries(self):
        print 'entered make_dictionaries'
        print '\n'
        self.f_iters = 0
        self.d_f = []
        self.d_x = []
        self.d_M = []
        self.d_Erh = []
        self.d_Frh = []
        self.d_Prh = []
        self.d_Ert = []
        self.d_Frt = []
        self.d_Prt = []
        self.d_Im_fwd = []
        self.d_Im_rev = []
        self.d_errs = []
        self.f_err = []
        self.f_L2 = []
        self.Er_err = []
        self.Er_L2 = []
        self.Fr_err = []
        self.Fr_L2 = []
        self.Pr_err = []
        self.Pr_L2 = []
        self.Im_err = []
        self.Machs_to_add = []
        print 'leaving make_dictionaries'
        print '\n'

    def update_dictionaries(self):
        print 'entered update_dictionaries'
        print '\n'
        self.d_f.append(self.f)
        self.d_x.append(self.x_RT)
        self.d_M.append(self.Mach_RT)
        self.d_Prt.append(self.P_RT)
        self.d_Prh.append(self.Pr_RT)
        Er = fnctn.rad_energy_density(self.Pr_RT, self.Mach_RT, self)
        self.d_Erh.append(Er)
        Fr = fnctn.rad_flux(self.Pr_RT, self.Mach_RT, self)
        self.d_Frh.append(Fr)
        self.d_Ert.append(self.E_RT)
        self.d_Frt.append(self.F_RT)
        self.d_Im_fwd.append(self.Im_fwd)
        self.d_Im_rev.append(self.Im_rev)
        self.d_errs.append([self.f_err, self.f_L2, self.Er_err, self.Er_L2, \
                            self.Fr_err, self.Fr_L2, self.Pr_err, self.Pr_L2])
        print 'leaving update_dictionaries'
        print '\n'

    def compute_errors(self):
    # the method 'update_dictionaries' has already been called
        print 'entered compute_errors'
        print '\n'
        f = self.f
        x = self.x_RT
        Sn = self.Sn
        Erh = fnctn.rad_energy_density(self.Pr_RT, self.Mach_RT, self)
        Frh = fnctn.rad_flux(self.Pr_RT, self.Mach_RT, self)
        Prh = self.Pr_RT
        Ert = self.E_RT
        Frt = self.F_RT
        Prt = self.P_RT
        Im_fwd = self.Im_fwd
        Im_rev = self.Im_rev
        f_iters = self.f_iters
        Im_fwd_err = []
        Im_rev_err = []
        bool_Im_err = 0
        if (f_iters == 0):
            f_OldOnNewx = numpy.ones(numpy.shape(f)) / 3.
            for i in range(Sn / 2):
                Im_fwd_err.append(numpy.max(Im_fwd[i]))
                Im_rev_err.append(numpy.max(Im_rev[i]))
        else:
            f_OldOnNewx = numpy.interp(x, self.d_x[f_iters - 1], \
                                          self.d_f[f_iters - 1])
            for i in range(Sn / 2):
                Im_fwd_OldOnNewx = numpy.interp(x, \
                                                self.d_x[f_iters - 1], \
                                                self.d_Im_fwd[f_iters - 1][i])
                Im_rev_OldOnNewx = numpy.interp(x, \
                                                self.d_x[f_iters - 1], \
                                                self.d_Im_rev[f_iters - 1][i])
                val = numpy.abs(Im_fwd[i] - Im_fwd_OldOnNewx) / Im_fwd[i]
                Im_fwd_err.append(numpy.max(val))
                val = numpy.abs(Im_rev[i] - Im_rev_OldOnNewx) / Im_rev[i]
                Im_rev_err.append(numpy.max(val))
            print 'Im_fwd_err = ', Im_fwd_err
            print '\n'
            print 'Im_rev_err = ', Im_rev_err
            print '\n'
        Im_err = scipy.append(Im_fwd_err, Im_rev_err)
        print 'Im_err = ', Im_err
        print '\n'
        f_err = numpy.max(numpy.abs(f - f_OldOnNewx) / f)
        val  = (f[2:-1] + f[1:-2]) / 2.
        val -= (f_OldOnNewx[2:-1] + f_OldOnNewx[1:-2]) / 2.
        val *= val * (x[2:-1] - x[1:-2])
        val  = numpy.sqrt(numpy.sum(val) / (x[-1] - x[0]))
        f_L2 = val
        Er_err = numpy.max(numpy.abs(Erh - Ert) / Erh)
        val  = (Erh[2:-1] + Erh[1:-2]) / 2. - (Ert[2:-1] + Ert[1:-2]) / 2.
        val *= val * (x[2:-1] - x[1:-2])
        val  = numpy.sqrt(numpy.sum(val) / (x[-1] - x[0]))
        Er_L2 = val
        Fr_err = numpy.max(numpy.abs((Frh - Frt) / Frh))
        val  = (Frh[2:-1] + Frh[1:-2]) / 2. - (Frt[2:-1] + Frt[1:-2]) / 2.
        val *= val * (x[2:-1] - x[1:-2])
        val  = numpy.sqrt(numpy.sum(val) / (x[-1] - x[0]))
        Fr_L2 = val
        Pr_err = numpy.max(numpy.abs(Prh - Prt) / Prh)
        val  = (Prh[2:-1] + Prh[1:-2]) / 2. - (Prt[2:-1] + Prt[1:-2]) / 2.
        val *= val * (x[2:-1] - x[1:-2])
        val  = numpy.sqrt(numpy.sum(val) / (x[-1] - x[0]))
        Pr_L2 = val
        print_bool = "The following errors are new lows:\n"
        if (self.f_iters == 0):
            self.f_err_low = f_err
            self.f_L2_low = f_L2
            self.Er_err_low = Er_err
            self.Er_L2_low = Er_L2
            self.Fr_err_low = Fr_err
            self.Fr_L2_low = Fr_L2
            self.Pr_err_low = Pr_err
            self.Pr_L2_low = Pr_L2
            self.Im_err_low = numpy.append(Im_fwd_err, Im_rev_err)
            bool_f_err = 1
            bool_f_L2 = 1
            bool_Er_err = 1
            bool_Er_L2 = 1
            bool_Fr_err = 1
            bool_Fr_L2 = 1
            bool_Pr_err = 1
            bool_Pr_L2 = 1
            bool_Im_err = 1
        else:
            for i in range(Sn/2):
               if (int(Im_fwd_err[i] < self.Im_err_low[i])):
                   self.Im_err_low[i] = Im_fwd_err[i]
                   bool_Im_err = 1
                   print_bool += "Im_err_fwd[" + str(i) + "]\n"
               if (int(Im_rev_err[i] < self.Im_err_low[i + Sn/2])):
                   self.Im_err_low[i + Sn/2] = Im_rev_err[i]
                   bool_Im_err = 1
                   print_bool += "Im_err_rev[" + str(i) + "]\n"
            if (int(f_err < self.f_err_low)):
                self.f_err_low = f_err
                bool_f_err = 1
                print_bool += "f_err\n"
            else:
                bool_f_err = 0
            if (int(f_L2 < self.f_L2_low) & int(f_L2 > 0)):
                self.f_L2_low = f_L2
                bool_f_L2 = 1
                print_bool += "f_L2\n"
            else:
                bool_f_L2 = 0
            if (int(Er_err < self.Er_err_low) & int(Er_err > 0)):
                self.Er_err_low = Er_err
                bool_Er_err = 1
                print_bool += "Er_err\n"
            else:
                bool_Er_err = 0
            if (int(Er_L2 < self.Er_L2_low) & int(Er_L2 > 0)):
                self.Er_L2_low = Er_L2
                bool_Er_L2 = 1
                print_bool += "Er_L2\n"
            else:
                bool_Er_L2 = 0
            if (int(Fr_err < self.Fr_err_low) & int(Fr_err > 0)):
                self.Fr_err_low = Fr_err
                bool_Fr_err = 1
                print_bool += "Fr_err\n"
            else:
                bool_Fr_err = 0
            if (int(Fr_L2 < self.Fr_L2_low) & int(Fr_L2 > 0)):
                self.Fr_L2_low = Fr_L2
                bool_Fr_L2 = 1
                print_bool += "Fr_L2\n"
            else:
                bool_Fr_L2 = 0
            if (int(Pr_err < self.Pr_err_low) & int(Pr_err > 0)):
                self.Pr_err_low = Pr_err
                bool_Pr_err = 1
                print_bool += "Pr_err\n"
            else:
                bool_Pr_err = 0
            if (int(Pr_L2 < self.Pr_L2_low) & int(Pr_L2 > 0)):
                self.Pr_L2_low = Pr_L2
                bool_Pr_L2 = 1
                print_bool += "Pr_L2\n"
            else:
                bool_Pr_L2 = 0
        self.make_RT_solution_bool = 0
        print 'self.make_RT_solution_bool = ', self.make_RT_solution_bool
        print '\n'
        if ((f_err > self.f_tol) & \
            (bool_f_err | bool_f_L2 | bool_Er_err | bool_Er_L2 | bool_Fr_err | \
             bool_Fr_L2 | bool_Pr_err | bool_Pr_L2 | bool_Im_err)):
            self.make_RT_solution_bool = 1
            print 'self.make_RT_solution_bool = ', self.make_RT_solution_bool
        if (not self.make_RT_solution_override):
            self.make_RT_solution_bool = 0
        self.f_err.append(f_err)
        self.f_L2.append(f_L2)
        self.Er_err.append(Er_err)
        self.Er_L2.append(Er_L2)
        self.Fr_err.append(Fr_err)
        self.Fr_L2.append(Fr_L2)
        self.Pr_err.append(Pr_err)
        self.Pr_L2.append(Pr_L2)
        if (self.f_iters != 0):
            self.Im_err.append(Im_err)
        if self.f_iters < 5:
            for i_iters in range(self.f_iters + 1):
                self.print_errors(i_iters)
        else:
            for i_iters in range(self.f_iters - 5, self.f_iters + 1):
                self.print_errors(i_iters)
        print print_bool
        print 'leaving compute_errors'
        print '\n'

    def print_errors(self, i_iters):
        print 'entered print_errors'
        print '\n'
        print 'For iteration ', i_iters
        print '\n'
        print 'The VEF relative error              is ', self.f_err[i_iters]
        if (i_iters == 0):
            ith  = numpy.abs(self.d_f[i_iters] - 1. / 3.)
            ith /= self.d_f[i_iters]
            ith  = numpy.argmax(ith)
        else:
            ith = numpy.argmax(numpy.abs(self.d_f[i_iters] - numpy.interp(self.d_x[i_iters], self.d_x[i_iters - 1], self.d_f[i_iters - 1])) / self.d_f[i_iters])
        print "This ocurrs at the " + str(ith) + " value of f"
        len_f = len(self.d_f[i_iters])
        print "and f contains " + str(len(self.d_f[i_iters])) + " elements"
        if (i_iters == 0):
            f_old_ith = 1./3.
        else:
            f_old_ith = numpy.interp(self.d_x[i_iters][ith], \
                                   self.d_x[i_iters - 1], self.d_f[i_iters - 1])
        print 'The value of f_old at that position is ', f_old_ith
        print 'The value of f     at that position is ', self.d_f[i_iters][ith]
        print 'and the value of x at that position is ', self.d_x[i_iters][ith]
        print 'The VEF L2 error (x normalized)     is ', self.f_L2[i_iters]
        print '\n'
        if (i_iters == self.f_iters):
            self.Machs_to_add.append(self.d_M[-1][ith])
        print 'The Er relative error               is ', self.Er_err[i_iters]
        ith  = abs(self.d_Erh[i_iters] - self.d_Ert[i_iters])
        ith /= self.d_Erh[i_iters]
        ith  = numpy.argmax(ith)
        print "This is the " + str(ith) + " value of Er,"
        len_RED_RT = len(self.d_Erh[i_iters])
        print "and Er contains " + str(len_RED_RT) + " elements."
        print 'The value of Erh at that position   is ', self.d_Erh[i_iters][ith]
        print 'The value of Ert at that position   is ', self.d_Ert[i_iters][ith]
        print 'and the value of x at that position is ', self.d_x[i_iters][ith]
        print 'The Er L2 error (x normalized)      is ', self.Er_L2[i_iters]
        print '\n'
        print 'The Fr relative error               is ', self.Fr_err[i_iters]
        ith  = abs(self.d_Frh[i_iters] - self.d_Frt[i_iters]) 
        ith /= self.d_Frh[i_iters]
        ith  = numpy.argmax(ith)
        print "This is the " + str(ith) + " value of Fr,"
        len_radFlux_RT = len(self.d_Frh[i_iters])
        print "and Fr contains " + str(len_radFlux_RT) + " elements."
        print 'The value of Frh at that position   is ', self.d_Frh[i_iters][ith]
        print 'The value of Frt at that position   is ', self.d_Frt[i_iters][ith]
        print 'and the value of x at that position is ', self.d_x[i_iters][ith]
        print 'The Fr L2 error (x normalized)      is ', self.Fr_L2[i_iters]
        print '\n'
        print 'The P_RT relative error             is ', self.Pr_err[i_iters]
        ith  = abs(self.d_Prh[i_iters] - self.d_Prt[i_iters])
        ith /= self.d_Prh[i_iters]
        ith  = numpy.argmax(ith)
        print "This is the " + str(ith) + " value of Pr,"
        len_fRED_RT = len(self.d_Prh[i_iters])
        print "and Pr contains " + str(len_fRED_RT) + " elements."
        print 'The value of Prh at that position   is ', self.d_Prh[i_iters][ith]
        print 'The value of Prt at that position   is ', self.d_Prt[i_iters][ith]
        print 'and the value of x at that position is ', self.d_x[i_iters][ith]
        print 'The Pr L2 error (x normalized)      is ', self.Pr_L2[i_iters]
        print '\n'
        print 'leaving print_errors'
        print '\n'

    def pickle_dictionary_variables(self):
        print 'entered pickle_dictionary_variables'
        print '\n'
        data = {'d_x':self.d_x, 'd_M':self.d_M, 'd_Erh':self.d_Erh, \
                'd_Frh':self.d_Frh, 'd_Prh':self.d_Prh, 'd_Ert':self.d_Ert, \
                'd_Frt':self.d_Frt, 'd_Prt':self.d_Prt, 'd_f':self.d_f, \
                'd_Im_fwd':self.d_Im_fwd, 'd_Im_rev':self.d_Im_rev, \
                'd_errs':self.d_errs}
        file_name = "data_dictionary_fiter" + str(self.f_iters) + ".pickle" 
        open_file = open(file_name, 'w')
        pickle.dump(data, open_file)
        open_file.close()
        print 'leaving pickle_dictionary_variables'
        print '\n'

    def pickle_prob(self):
        print 'entered pickle_prob'
        print '\n'
        if (self.f_iters != 0):
            os.system('mv prob.pickle removing_prob.pickle')
        open_file = open('prob.pickle', 'w')
        pickle.dump(self, open_file)
        open_file.close()
        if (self.f_iters != 0):
            os.system('rm -fr removing_prob.pickle')
        print 'leaving pickle_prob'
        print '\n'

    def append_metadata(self):
        print 'entered append_metadata'
        print '\n'
        index_x_is_0 = numpy.argmin(numpy.abs(self.x))
        self.index_x_is_0 = index_x_is_0
        self.Tp_val = self.Tm[index_x_is_0]
        self.Ts_val = self.Tm[index_x_is_0 + 1]
        self.Tf_val = self.T1
        self.Tmax_val = numpy.max(self.Tm)
        self.theta_ps_val = self.Tr[index_x_is_0]
        self.theta_max_val = numpy.max(self.Tr)
        self.Mp_val = self.Mach[index_x_is_0]
        self.Ms_val = self.Mach[index_x_is_0 + 1]
        self.Mf_val = self.Mach[-1]
        self.Mmax_val = numpy.max(self.Mach)
        self.x_FWHM_val = 0.
        if (self.Tmax_val > self.Tf_val + self.int_atol):
            Thalf = (self.Tmax_val + self.Tf_val) / 2.
            index_Tmax = numpy.argmax(\
                         fnctn.mat_temp(self.Pr_relaxation[::-1], \
                                        self.Mach_relaxation[::-1], self))
            self.index_Tmax = index_Tmax
            if (Thalf < self.Tp_val):
                x_left = numpy.interp(Thalf, \
                                      fnctn.mat_temp(self.Pr_precursor, \
                                                     self.Mach_precursor, \
                                                     self), \
                                      self.x_precursor)
            elif (self.Tp_val < Thalf < self.Ts_val):
                x_left = 0.
            else:
                x_left = numpy.interp(Thalf, \
                                      fnctn.mat_temp( \
                                      self.Pr_relaxation[::-1][:index_Tmax], \
                                      self.Mach_relaxation[::-1][:index_Tmax], \
                                      self), \
                                      self.x_relaxation[::-1][:index_Tmax])
            x_right = numpy.interp(Thalf, \
                                   fnctn.mat_temp( \
                                   self.Pr_relaxation[index_Tmax:], \
                                   self.Mach_relaxation[index_Tmax:], \
                                   self), \
                                   self.x_relaxation[index_Tmax:])
            self.x_left = x_left
            self.x_right = x_right
            self.x_FWHM_val = x_right - x_left
        print 'leaving append_metadata'
        print '\n'

    def plot_polar(self):
        print 'entered create polar plots'
        mus = numpy.linspace(-1., -0.005, 1.e2)
        self.mus = scipy.append(mus, -mus[::-1])
        self.x_slice = 0.
        index_val = numpy.argmin(numpy.abs(self.x_RT - self.x_slice)) + 1
        I_vals = fnctn.I_mu_x(self.mus, self.x_RT[index_val], self)
        val  = fnctn.mat_temp(self.Pr_RT[index_val], \
                              self.Mach_RT[index_val], \
                              self)
        val *= val**3 / 4. / numpy.pi
        self.I_vals_frac = I_vals / val
        matplotlib.pyplot.polar(numpy.arccos(self.mus), self.I_vals_frac, \
                                'k', linewidth = 3.)
        matplotlib.pyplot.polar(\
           - numpy.arccos(self.mus), self.I_vals_frac, 'k', linewidth = 3.,\
           label = r'$\frac{I \left( \mu \right)}{\left( \frac{T^4}{4 \pi} \right)}$')
        matplotlib.pyplot.xticks(scipy.append(numpy.arccos(self.mus_16), -numpy.arccos(self.mus_16)), [' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' '], fontsize = 20.)
        matplotlib.pyplot.yticks([0.5, 1.0], [0.5, 1.0], fontsize = 20.)
        matplotlib.pyplot.legend(loc = 'upper right', frameon = False, \
                                 bbox_to_anchor = (1., 1.25))
        matplotlib.pyplot.show()

    def plot_compare_to_LE(self):
        open_file = open('dictionary.pickle', 'r')
        self.dic_data = pickle.load(open_file)
        open_file.close()
        self.f_iters = numpy.shape(self.dic_data['d_f'])[0]
        self.x_RT = self.dic_data['d_x'][-1]
        self.xdiff = self.dic_data['d_x'][0]
        self.Mach_RT = self.dic_data['d_M'][-1]
        self.Mdiff = self.dic_data['d_M'][0]
        self.Erh = self.dic_data['d_Erh'][-1]
        self.Erhdiff = self.dic_data['d_Erh'][0]
        self.Ert = self.dic_data['d_Ert'][-1]
        self.Frh = self.dic_data['d_Frh'][-1]
        self.Frhdiff = self.dic_data['d_Frh'][0]
        self.Frt = self.dic_data['d_Frt'][-1]
        self.Prh = self.dic_data['d_Prh'][-1]
        self.Prhdiff = self.dic_data['d_Prh'][0]
        self.Prt = self.dic_data['d_Prt'][-1]
        self.f = self.dic_data['d_f'][-1]
        os.chdir('../')

        self.Tm = fnctn.mat_temp(self.Pr, self.Mach, self)
        self.Er = fnctn.rad_energy_density(self.Pr, self.Mach, self)
        self.Tr = fnctn.rad_temp(self.Pr, self.Mach, self)
        self.Speed = fnctn.mat_speed(self.Pr, self.Mach, self)
        self.Density = fnctn.mat_density(self.Pr, self.Mach, self)
        self.Pressure = fnctn.mat_pres(self.Pr, self.Mach, self)
        self.Fr = fnctn.rad_flux(self.Pr, self.Mach, self)

    def write_data(self):
        print 'entered write_data'
        print '\n'
        open_file = open('x_data.txt', 'wb')
        numpy.savetxt(open_file, self.x)
        open_file.close()
        open_file = open('Mach_data.txt', 'wb')
        numpy.savetxt(open_file, self.Mach)
        open_file.close()
        open_file = open('Tm_data.txt', 'wb')
        numpy.savetxt(open_file, self.Tm)
        open_file.close()
        open_file = open('Tr_data.txt', 'wb')
        numpy.savetxt(open_file, self.Tr)
        open_file.close()
        open_file = open('Er_data.txt', 'wb')
        numpy.savetxt(open_file, self.Er)
        open_file.close()
        open_file = open('Density_data.txt', 'wb')
        numpy.savetxt(open_file, self.Density)
        open_file.close()
        open_file = open('Speed_data.txt', 'wb')
        numpy.savetxt(open_file, self.Speed)
        open_file.close()
        open_file = open('Pressure_data.txt', 'wb')
        numpy.savetxt(open_file, self.Pressure)
        open_file.close()
        open_file = open('Fr_data.txt', 'wb')
        numpy.savetxt(open_file, self.Fr)
        open_file.close()
        print 'leaving write_data'
        print '\n'

class RageShockProfiles:
    def __init__(self, incoming, nED_profile):
        print 'entered RageShockProfiles __init__ for M0 = ', incoming.M0
        print '\n'
        self.M0 = copy.deepcopy(incoming.M0)
        self.sigA = copy.deepcopy(incoming.sigA)
        self.sigS = copy.deepcopy(incoming.sigS)
        self.expDensity_abs = copy.deepcopy(incoming.expDensity_abs)
        self.expTemp_abs = copy.deepcopy(incoming.expTemp_abs)
        self.expDensity_scat = copy.deepcopy(incoming.expDensity_scat)
        self.expTemp_scat = copy.deepcopy(incoming.expTemp_scat)
        self.gamma = copy.deepcopy(incoming.gamma)
        self.sound = copy.deepcopy(incoming.sound)
        self.Tref = copy.deepcopy(incoming.Tref)
        self.rho0 = copy.deepcopy(incoming.rho0)
        self.Cv = copy.deepcopy(incoming.Cv)
        self.T0 = copy.deepcopy(incoming.T0)
        self.ar = copy.deepcopy(incoming.ar)
        self.P0 = copy.deepcopy(incoming.P0)
        self.dxset = copy.deepcopy(incoming.dxset)
        self.runtime = copy.deepcopy(incoming.runtime)
        self.freezeWidth = copy.deepcopy(incoming.freezeWidth)
        self.dumpnum = copy.deepcopy(incoming.dumpnum)
        self.rho1 = copy.deepcopy(nED_profile.rho1)
        self.T1 = copy.deepcopy(nED_profile.T1)
        self.speed1 = copy.deepcopy(nED_profile.speed1)
        self.numlev = copy.deepcopy(incoming.numlev)
        self.x = copy.deepcopy(nED_profile.x)
        self.freezeWidth = 2. * self.x[-1]
        self.shockTubeLength = self.sound * (self.M0 - self.speed1) \
                               * self.runtime * 2. + 4. * self.freezeWidth \
                               + (self.x[-1] - self.x[0])
        self.imxset = int(numpy.ceil(self.shockTubeLength / self.dxset))
        if (self.imxset % 2 == 1):
            self.imxset += 1
        self.shockTubeLength = self.imxset * self.dxset
        self.x = scipy.append(self.x[0] - self.freezeWidth, self.x)
        self.x = scipy.append(self.x, self.x[-1] + self.freezeWidth)
        self.Mach = copy.deepcopy(nED_profile.Mach) 
        self.Mach = scipy.append(self.Mach[0], self.Mach)
        self.Mach = scipy.append(self.Mach, self.Mach[-1])
        self.Tm = copy.deepcopy(nED_profile.Tm)
        self.Tm = scipy.append(self.Tm[0], self.Tm)
        self.Tm = scipy.append(self.Tm, self.Tm[-1])
        self.Tr = copy.deepcopy(nED_profile.Tr)
        self.Tr = scipy.append(self.Tr[0], self.Tr)
        self.Tr = scipy.append(self.Tr, self.Tr[-1])
        self.Density = copy.deepcopy(nED_profile.Density)
        self.Density = scipy.append(self.Density[0], self.Density)
        self.Density = scipy.append(self.Density, self.Density[-1])
        self.Pressure = copy.deepcopy(nED_profile.Pressure)
        self.Pressure = scipy.append(self.Pressure[0], self.Pressure)
        self.Pressure = scipy.append(self.Pressure, self.Pressure[-1])
        self.Fr = copy.deepcopy(nED_profile.Fr)
        self.Fr = scipy.append(self.Fr[0], self.Fr)
        self.Fr = scipy.append(self.Fr, self.Fr[-1])
        self.Pr = copy.deepcopy(nED_profile.Pr)
        self.Pr = scipy.append(self.Pr[0], self.Pr)
        self.Pr = scipy.append(self.Pr, self.Pr[-1])
        self.Speed = copy.deepcopy(nED_profile.Speed)
        self.Speed = scipy.append(self.Speed[0], self.Speed)
        self.Speed = scipy.append(self.Speed, self.Speed[-1])
        print 'leaving RageShockProfiles __init__ for M0 = ', incoming.M0
        print '\n'

    def write_fndx2x(self):
        print 'entered write_fndx2x'
        print '\n'
        x = self.x + self.shockTubeLength - 3. * self.freezeWidth
        v = self.sound * (self.Speed - self.Speed[0])
        rho = self.Density
        e = self.Cv * self.Tm * self.Tref
        Er = self.ar * (self.Tref * self.Tr)**4
        for i in range(numpy.size(self.x)):
            if (i == 0):
                print_statement = str(numpy.size(self.x)) + ' 1 0\n'
            print_statement += str(i + 1) + ' ' + str(x[i]) + ' ' + str(v[i])
            print_statement += ' ' + str(rho[i]) + ' ' + str(e[i]) + ' ' + str(Er[i]) + '\n'
        open_file = open('fndx2x.dat', 'wb')
        print >> open_file, print_statement
        open_file.close()

    def write_RAGE_runjob_script(self):
        print 'entered write_RAGE_runjob_script'
        print '\n'
        numpe = 2**int(numpy.floor(numpy.log(self.imxset) / numpy.log(10.)) + 1)
        numpe = min(numpe, 512)
        kappa0  = self.sigA * (self.Tref)**(- self.expTemp_abs)
        kappa0 *= (1. / self.rho0)**self.expDensity_abs
        open_file = open('run_job.csh', 'wb')
        print_statement  = "#!/bin/tcsh -f\n"
        print_statement += "#RJ BATCH      = yes\n"
        print_statement += "#RJ CHAIN      = yes\n"
        print_statement += "#RJ EXEC       = /yellow/usr/projects/eap/"
        print_statement += "releases/1603.01/xrage/xrage1603.01_debug_MP.x\n"
        print_statement += "#RJ GROUP      = dacodes\n"
        print_statement += "#RJ NUMPE      = " + str(numpe) + "\n"
        if (self.dxset > 5.e-3):
          print_statement += "#RJ TIME       = 01:00:00\n"
        else:
          print_statement += "#RJ TIME       = 10:00:00\n"
        print_statement += "\n"
        print_statement += "# source the eap .cshrc to get modules and stuff\n"
        print_statement += "source /yellow/usr/projects/eap/tools/../dotfiles/"
        print_statement += ".cshrc\n"
        print_statement += "\n"
        print_statement += "# cleanup from previous run\n"
        print_statement += "run_job_cleanup.pl\n"
        print_statement += "\n"
        print_statement += "# generate teos_out file\n"
        print_statement += "if ( -e  'teos.in' && ! -e 'deflect.teos' ) then\n"
        print_statement += "  &&&RJ_CMD_PRUN&&& teos.in\n"
        print_statement += "endif\n"
        print_statement += "\n"
        print_statement += "# run\n"
        print_statement += "&&&RJ_CMD_PRUN&&& radshock_RAGE.input"
        problem_name = "radshock_m"
        if (str(self.M0).split('.')[-1] == '0'):
            problem_name += str(self.M0).split('.')[0]
        else:
            problem_name += str(self.M0)
        problem_name += "_dxset" + '{:.10f}'.format(self.dxset)
        print_statement += ' -v pname=' + str(problem_name)
        print_statement += ' -v finalDensity=' + str(self.rho1)
        finalLabFrameSpeed = self.sound \
                             * (fnctn.mat_speed(self.Pr[-1], \
                                                self.Mach[-1], \
                                                self) \
                              - fnctn.mat_speed(self.Pr[0], \
                                                self.Mach[0], \
                                                self))
        self.finalLabFrameSpeed = finalLabFrameSpeed
        print_statement += ' -v finalLabFrameSpeed=' + str(finalLabFrameSpeed)
        print_statement += ' -v finalTemperature=' + str(self.T1 * self.Tref)
        print_statement += ' -v dxset=' + str(self.dxset)
        print_statement += ' -v imxset=' + str(self.imxset)
        print_statement += ' -v numlev=' + str(self.numlev)
        print_statement += ' -v gamma=' + str(self.gamma)
        print_statement += ' -v Cv=' + str(self.Cv)
        print_statement += ' -v rho0=' + str(self.rho0)
        print_statement += ' -v Tref=' + str(self.Tref)
        print_statement += ' -v sigma_scat=' + str(self.sigS)
        print_statement += ' -v kappa0=' + str(kappa0)
        print_statement += ' -v powd=' + str(self.expDensity_abs - 1)
        print_statement += ' -v powt=' + str(numpy.abs(self.expTemp_abs))
        print_statement += ' -v runtime=' + str(self.runtime)
        print_statement += ' -v freezeWidth=' + str(self.freezeWidth)
        print_statement += ' -v dumpnum=' + str(self.dumpnum)
        print_statement += ' -v shockTubeLength=' + str(self.shockTubeLength)
        print_statement += ' -v secmax=$RJ_TIME_REMAINING_B\n'
        print_statement += '\n'
        print_statement += '# cleanup and check for continue\n'
        print_statement += 'run_job_cleanup.pl --check\n'
        print >> open_file, print_statement
        open_file.close()

    def run_RAGE(self):
        print 'entered run_RAGE'
        print '\n'
        os.system("run_job.pl go")

    def read_RAGE_AMHC(self):
        print 'entered read_RAGE_AMHC'
        print '\n'
        os.system('cp *-dmp* dumps/')
        os.chdir('dumps')
        dump_files = os.listdir(os.getcwd())
        dump_files.sort()
        self.data = []
        self.Tm_RAGE = []
        self.Tr_RAGE = []
        self.Speed_RAGE = []
        self.Mach_RAGE = []
        self.Density_RAGE = []
        self.Pressure_RAGE = []
        self.hist_time = []
        for file in dump_files:
            self.data.append(rage_utils.rage_dump1D(file))
            tmp = self.data[-1]['tev'] / self.Tref
            self.Tm_RAGE.append(tmp)
            tmp = (self.data[-1]['rade'] / self.ar)**(1./4.) / self.Tref
            self.Tr_RAGE.append(tmp)
            tmp = self.data[-1]['cell_momentum'] / self.data[-1]['mass']
            tmp = tmp / self.sound + self.M0
            self.Speed_RAGE.append(tmp)
            tmp = self.Speed_RAGE[-1] / numpy.sqrt(self.Tm_RAGE[-1][0])
            self.Mach_RAGE.append(tmp)
            self.Density_RAGE.append(self.data[-1]['dens'])
            self.Pressure_RAGE.append(self.Density_RAGE[-1] \
                                      * self.Tm_RAGE[-1] / self.gamma)
            self.hist_time.append(self.data[-1]['time'])
        os.chdir('../')

    def shift_data(self):
        print 'entered shift_data'
        print '\n'
        metric = lambda delta, var_data_x, var_data_y, var_model_x, var_model_y: numpy.sqrt(numpy.average((numpy.interp(var_data_x, var_model_x + delta, var_model_y) - var_data_y)**2))
        self.shift_Tm_midpoint = numpy.zeros(numpy.size(self.hist_time))
        self.shift_Tm_galilean = numpy.zeros(numpy.size(self.hist_time))
        self.shift_Tm_slope = numpy.zeros(numpy.size(self.hist_time))
        self.shift_Tm_peak = numpy.zeros(numpy.size(self.hist_time))
        self.shift_Tm_background = numpy.zeros(numpy.size(self.hist_time))
        self.shift_Tr = numpy.zeros(numpy.size(self.hist_time))
        self.shift_Speed = numpy.zeros(numpy.size(self.hist_time))
        self.shift_Mach = numpy.zeros(numpy.size(self.hist_time))
        self.shift_Density = numpy.zeros(numpy.size(self.hist_time))
        self.L2error_Tm_midpoint = numpy.zeros(numpy.size(self.hist_time))
        self.L2error_Tm_galilean = numpy.zeros(numpy.size(self.hist_time))
        self.L2error_Tm_slope = numpy.zeros(numpy.size(self.hist_time))
        self.L2error_Tm_peak = numpy.zeros(numpy.size(self.hist_time))
        self.L2error_Tm_background = numpy.zeros(numpy.size(self.hist_time))
        self.L2error_Tr = numpy.zeros(numpy.size(self.hist_time))
        self.L2error_Speed = numpy.zeros(numpy.size(self.hist_time))
        self.L2error_Mach = numpy.zeros(numpy.size(self.hist_time))
        self.L2error_Density = numpy.zeros(numpy.size(self.hist_time))
        self.L2_avg_Tm = 0.
        self.L2_stddev_Tm = 0.
        self.L2_avg_Tr = 0.
        self.L2_stddev_Tr = 0.
        self.L2_avg_Speed = 0.
        self.L2_stddev_Speed = 0.
        self.L2_avg_Mach = 0.
        self.L2_stddev_Mach = 0.
        self.L2_avg_Density = 0.
        self.L2_stddev_Density = 0.
        self.Tm_galilean_L2 = numpy.zeros(numpy.size(self.hist_time))
        self.galilean_guess = numpy.zeros(numpy.size(self.hist_time))
        self.peak_max_guess = numpy.zeros(numpy.size(self.hist_time))
        self.slope_max_guess = numpy.zeros(numpy.size(self.hist_time))
        self.background_guess = numpy.zeros(numpy.size(self.hist_time))
        self.initial_guess = numpy.zeros(numpy.size(self.hist_time))
#         self.shift_slidex = numpy.zeros(numpy.size(self.hist_time))
#         self.L2_slidex = numpy.zeros(numpy.size(self.hist_time))
#         self.opt_slidex = numpy.zeros(numpy.size(self.hist_time))
#         self.L2opt_slidex = numpy.zeros(numpy.size(self.hist_time))
        for elem in range(numpy.size(self.hist_time)):
            metric = lambda delta, var_data_x, var_data_y, var_model_x, var_model_y: numpy.sqrt(numpy.average((numpy.interp(var_data_x, var_model_x + delta, var_model_y) - var_data_y)**2))
#             slidex = numpy.linspace(self.data[elem]['x'][0], self.data[elem]['x'][-1], 4 * self.imxset)
#             L2s_slidex = numpy.zeros(len(slidex))
#             for i in range(len(slidex)):
#                 L2s_slidex[i] = metric(slidex[i], self.data[elem]['x'], self.Tm_RAGE[elem], self.x, self.Tm)
#             self.shift_slidex[elem] = slidex[numpy.argmin(L2s_slidex)]
#             self.L2_slidex[elem] = L2s_slidex[numpy.argmin(L2s_slidex)]
#             print "self.hist_time[elem]             = ", self.hist_time[elem]
#             print "L2_slidex[elem, 1]          = ", L2s_slidex[numpy.argmin(L2s_slidex)]
#             print "brute_slide[elem, 1]          = ", slidex[numpy.argmin(L2s_slidex)]
#             L2slidex_background = 0.9 * numpy.min(L2s_slidex) + 0.1 * numpy.max(L2s_slidex)
#             slidex_index = numpy.argmin(L2s_slidex)
#             left_slidex_index = numpy.argmin(numpy.abs(L2slidex_background - L2s_slidex[:slidex_index - 1]))
#             right_slidex_index = numpy.argmin(numpy.abs(L2s_slidex[slidex_index + 1:] - L2slidex_background)) + slidex_index + 1
#             slidex = numpy.arange(slidex[left_slidex_index], slidex[right_slidex_index] + 1.e-4, 1.e-4)
#             L2s_slidex = numpy.zeros(len(slidex))
#             for i in range(len(slidex)):
#                 L2s_slidex[i] = metric(slidex[i], self.data[elem]['x'], self.Tm_RAGE[elem], self.x, self.Tm)
#             print "numpy.min(L2s_slidex)            = ", numpy.min(L2s_slidex)
#             print "brute_slide[elem, 1]             = ", slidex[numpy.argmin(L2s_slidex)]
#             self.opt_slidex[elem] = slidex[numpy.argmin(L2s_slidex)]
#             self.L2opt_slidex[elem] = L2s_slidex[numpy.argmin(L2s_slidex)]
#             print "self.L2opt_slidex[elem]       = ", L2s_slidex[numpy.argmin(L2s_slidex)]
            guess1  = self.shockTubeLength - 2.5 * self.freezeWidth
            guess1 -= self.hist_time[elem] * self.sound * (self.M0 - self.speed1)
            self.galilean_guess[elem] = guess1
            guess2 = numpy.argmax(numpy.diff(self.Tm_RAGE[elem]))
            guess2 = (self.data[elem]['x'][guess2] + self.data[elem]['x'][guess2 + 1]) / 2.
            self.slope_max_guess[elem] = guess2
            guess3 = self.data[elem]['x'][numpy.argmax(self.Tm_RAGE[elem])]
            guess3 -= self.x[numpy.argmax(self.Tm)]
            self.peak_max_guess[elem] = guess3
            guess4 = 0.9 * self.T0 + 0.1 * self.T1
            guess4a = numpy.argmin(numpy.abs(self.Tm_RAGE[elem] - guess4))
            guess4a = self.data[elem]['x'][guess4a]
            guess4b = self.x[numpy.argmin(numpy.abs(self.Tm - guess4))]
            guess4 = guess4a - guess4b
            self.background_guess[elem] = guess4
            initial_guess = min(min(guess1, guess2), min(guess3, guess4))
            initial_guess += max(max(guess1, guess2), min(guess3, guess4))
            initial_guess /= 2.
            self.initial_guess[elem] = initial_guess
            self.shift_Tm_midpoint[elem] = scipy.optimize.fsolve(lambda tau: metric(tau, self.data[elem]['x'], self.Tm_RAGE[elem], self.x, self.Tm), initial_guess, xtol = 1e-14)
            self.shift_Tm_galilean[elem] = scipy.optimize.fsolve(lambda tau: metric(tau, self.data[elem]['x'], self.Tm_RAGE[elem], self.x, self.Tm), guess1, xtol = 1e-14)
            self.shift_Tm_slope[elem] = scipy.optimize.fsolve(lambda tau: metric(tau, self.data[elem]['x'], self.Tm_RAGE[elem], self.x, self.Tm), guess2, xtol = 1e-14)
            self.shift_Tm_peak[elem] = scipy.optimize.fsolve(lambda tau: metric(tau, self.data[elem]['x'], self.Tm_RAGE[elem], self.x, self.Tm), guess3, xtol = 1e-14)
            self.shift_Tm_background[elem] = scipy.optimize.fsolve(lambda tau: metric(tau, self.data[elem]['x'], self.Tm_RAGE[elem], self.x, self.Tm), guess4, xtol = 1e-14)
            self.shift_Tr[elem] = scipy.optimize.fsolve(lambda tau: metric(tau, self.data[elem]['x'], self.Tr_RAGE[elem], self.x, self.Tr), initial_guess)
            self.shift_Speed[elem] = scipy.optimize.fsolve(lambda tau: metric(tau, self.data[elem]['x'], self.Speed_RAGE[elem], self.x, self.Speed), initial_guess)
            self.shift_Mach[elem] = scipy.optimize.fsolve(lambda tau: metric(tau, self.data[elem]['x'], self.Mach_RAGE[elem], self.x, self.Mach), initial_guess)
            self.shift_Density[elem] = scipy.optimize.fsolve(lambda tau: metric(tau, self.data[elem]['x'], self.Density_RAGE[elem], self.x, self.Density), initial_guess)
            self.L2error_Tm_midpoint[elem] = metric(self.shift_Tm_midpoint[elem], self.data[elem]['x'], self.Tm_RAGE[elem], self.x, self.Tm)
            self.L2error_Tm_galilean[elem] = metric(self.shift_Tm_galilean[elem], self.data[elem]['x'], self.Tm_RAGE[elem], self.x, self.Tm)
            self.L2error_Tm_slope[elem] = metric(self.shift_Tm_slope[elem], self.data[elem]['x'], self.Tm_RAGE[elem], self.x, self.Tm)
            self.L2error_Tm_peak[elem] = metric(self.shift_Tm_peak[elem], self.data[elem]['x'], self.Tm_RAGE[elem], self.x, self.Tm)
            self.L2error_Tm_background[elem] = metric(self.shift_Tm_background[elem], self.data[elem]['x'], self.Tm_RAGE[elem], self.x, self.Tm)
            self.Tm_galilean_L2[elem] = metric(self.galilean_guess[elem], self.data[elem]['x'], self.Tm_RAGE[elem], self.x, self.Tm)
            print '\n'
            print "self.initial_guess[elem]   = ", self.initial_guess[elem]
            print "self.galilean_guess[elem]   = ", self.galilean_guess[elem]
            print "self.slope_max_guess[elem]      = ", self.slope_max_guess[elem]
            print "self.peak_max_guess[elem]       = ", self.peak_max_guess[elem]
            print "self.background_guess[elem] = ", self.background_guess[elem]
            print '\n'
            print "self.shift_Tm_midpoint[elem]   = ", self.shift_Tm_midpoint[elem]
            print "self.shift_Tm_galilean[elem]   = ", self.shift_Tm_galilean[elem]
            print "self.shift_Tm_slope[elem]      = ", self.shift_Tm_slope[elem]
            print "self.shift_Tm_peak[elem]       = ", self.shift_Tm_peak[elem]
            print "self.shift_Tm_background[elem] = ", self.shift_Tm_background[elem]
            print '\n'
            print "self.L2error_Tm_midpoint[elem]   = ", self.L2error_Tm_midpoint[elem]
            print "self.L2error_Tm_galilean[elem]   = ", self.L2error_Tm_galilean[elem]
            print "self.L2error_Tm_slope[elem]      = ", self.L2error_Tm_slope[elem]
            print "self.L2error_Tm_peak[elem]       = ", self.L2error_Tm_peak[elem]
            print "self.L2error_Tm_background[elem] = ", self.L2error_Tm_background[elem]
            print '\n'
            self.L2error_Tr[elem] = metric(self.shift_Tr[elem], self.data[elem]['x'], self.Tr_RAGE[elem], self.x, self.Tr)
            self.L2error_Speed[elem] = metric(self.shift_Speed[elem], self.data[elem]['x'], self.Speed_RAGE[elem], self.x, self.Speed)
            self.L2error_Mach[elem] = metric(self.shift_Mach[elem], self.data[elem]['x'], self.Mach_RAGE[elem], self.x, self.Mach)
            self.L2error_Density[elem] = metric(self.shift_Density[elem], self.data[elem]['x'], self.Density_RAGE[elem], self.x, self.Density)
        self.L2_avg_Tm = numpy.average(self.L2error_Tm_midpoint)
        self.L2_stddev_Tm = numpy.sqrt(numpy.average((self.L2_avg_Tm - self.L2error_Tm_midpoint)**2))
        self.L2_avg_Tr = numpy.average(self.L2error_Tr)
        self.L2_stddev_Tr = numpy.sqrt(numpy.average((self.L2_avg_Tr - self.L2error_Tr)**2))
        self.L2_avg_Speed = numpy.average(self.L2error_Speed)
        self.L2_stddev_Speed = numpy.sqrt(numpy.average((self.L2_avg_Speed - self.L2error_Speed)**2))
        self.L2_avg_Mach = numpy.average(self.L2error_Mach)
        self.L2_stddev_Mach = numpy.sqrt(numpy.average((self.L2_avg_Mach - self.L2error_Mach)**2))
        self.L2_avg_Density = numpy.average(self.L2error_Density)
        self.L2_stddev_Density = numpy.sqrt(numpy.average((self.L2_avg_Density - self.L2error_Density)**2))

    def write_RAGE_firstmidlast_data(self):
        print 'entered write_RAGE_data'
        print '\n'
        first_time = 0
        mid_time = int(numpy.size(self.hist_time) / 2.)
        last_time = -1
        name_list = ['a','b','c']
        for num, val in enumerate([first_time, mid_time, last_time]):
            file_name = 'RAGE_data_x_' + name_list[num] + '.txt'
            open_file = open(file_name, 'wb')
            numpy.savetxt(open_file, self.data[num]['x'])
            open_file.close()
            file_name = 'RAGE_data_Tm_' + name_list[num] + '.txt'
            open_file = open(file_name, 'wb')
            numpy.savetxt(open_file, self.Tm_RAGE[num])
            open_file.close()
            file_name = 'RAGE_data_Tr_' + name_list[num] + '.txt'
            open_file = open(file_name, 'wb')
            numpy.savetxt(open_file, self.Tr_RAGE[num])
            open_file.close()
            file_name = 'RAGE_data_Speed_' + name_list[num] + '.txt'
            open_file = open(file_name, 'wb')
            numpy.savetxt(open_file, self.Speed_RAGE[num])
            open_file.close()
            file_name = 'RAGE_data_Mach_' + name_list[num] + '.txt'
            open_file = open(file_name, 'wb')
            numpy.savetxt(open_file, self.Mach_RAGE[num])
            open_file.close()
            file_name = 'RAGE_data_Density_' + name_list[num] + '.txt'
            open_file = open(file_name, 'wb')
            numpy.savetxt(open_file, self.Density_RAGE[num])
            open_file.close()
            file_name = 'RAGE_data_Pressure_' + name_list[num] + '.txt'
            open_file = open(file_name, 'wb')
            numpy.savetxt(open_file, self.Pressure_RAGE[num])
            open_file.close()
        fig = matplotlib.pyplot.figure()
        for i in range(numpy.size(self.hist_time)):
            if (i % 7 == 0):
                matplotlib.pyplot.plot(self.data[i]['x'], self.Tm_RAGE[i], 'b')
                matplotlib.pyplot.plot(self.x + self.shift_Tm_peak[i], self.Tm, 'g')
                matplotlib.pyplot.plot(self.data[i]['x'], numpy.interp(self.data[i]['x'], self.x + self.shift_Tm_peak[i], self.Tm), '--g')
        matplotlib.pyplot.xlim((self.x[1] + self.shift_Tm_peak[-1], self.x[-2] + self.shift_Tm_peak[0]))
        matplotlib.pyplot.ylim((0.95, 1.1 * numpy.max(self.Tm)))
        matplotlib.pyplot.title(r'dx = ' + str(self.dxset))
        save_str = "dxset_" + '{:.10f}'.format(self.dxset) + ".pdf"
        fig.savefig(save_str)
        mkdir_str  = 'RAGE_and_analytic_data_M' + str(self.M0) + '_dxset'
        mkdir_str += '{:.10f}'.format(self.dxset)
        if (mkdir_str not in os.listdir(os.getcwd())):
            mkdir_str += '/'
            os.mkdir(mkdir_str)
        os.system('cp RAGE_data* *_data.txt ' + save_str + mkdir_str)
