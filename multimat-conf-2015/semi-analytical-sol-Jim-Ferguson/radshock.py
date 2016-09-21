import os
import copy
import fnctn
import utils
import numpy
import pickle
import scipy
import scipy.integrate
import matplotlib.pyplot

class init_RadShock():
    def __init__(self, M0 = 1.2, rho0 = 1., Sn = 16, \
                 Tref = 100., Cv = 1.4472799784454e12, gamma = 5. / 3., \
                 sigA = 44.87940546534333, sigS = 0.4006204354152789, \
                 expDensity_abs = 2., expTemp_abs = -3.5, \
                 expDensity_scat = 1., expTemp_scat = 0., \
                 dxset = 5, numlev = 1, runtime = 5.e-7, \
                 freezeWidth = 10., dumpnum = 100, \
                 problem = 'JMF', \
#                  make_RH_solution_bool = True, \
                 make_RT_solution_override = True, \
                 eps = 1.e-6, int_tol = 1.e-13, integrator = 'ode'):
#: Some comment describing M0.  Sphinx will pick up this comment and the documentation for RadShock will present this comment for M0 nicely.  It's the #: that's important.
#: The initial Mach number.
        self.M0 = M0
#: The ambient, upstream equilibrium material density
        self.rho0 = rho0
#: The ambient, upstream equilibrium temperature
        self.Tref = Tref
        self.Cv = Cv
        self.gamma = gamma
        self.sigA = sigA
        self.sigS = sigS
        self.expDensity_abs = expDensity_abs
        self.expTemp_abs = expTemp_abs
        self.expDensity_scat = expDensity_scat
        self.expTemp_scat = expTemp_scat
        self.T0 = 1.
        self.Pr0 = self.T0**4 / 3.
        self.c = 2.99792458e10# [cm / s]
        self.ar = 137.20172# [erg / cm^3 - eV^4]
        self.sound = numpy.sqrt(self.gamma * (self.gamma - 1.) * self.Cv * \
                                self.Tref)
        self.C0 = self.c / self.sound
        self.P0 = self.ar * self.Tref**4 / (self.rho0 * self.sound**2)
        self.Tref = Tref
        self.dxset = dxset
        self.numlev = numlev
        self.runtime = runtime
        self.freezeWidth = freezeWidth
        self.dumpnum = dumpnum
        self.problem = problem
        self.eps = eps
        self.int_tol = int_tol
        self.integrator = integrator

class greyNED_RadShock(init_RadShock):
    '''
    Define the radiative-shock problem, and drive the solution
    '''
# prob = radshock.greyNED_RadShock(M0 = 3., sigA = 44.93983839817290, sigS = 0.4006, expDensity_abs = 1, expTemp_abs = -3.5)
    def __init__(self, M0 = 1.2, rho0 = 1., Sn = 16, \
                 Tref = 100., Cv = 1.4472799784454e12, gamma = 5. / 3., \
                 sigA = 44.87940546534333, sigS = 0.4006204354152789, \
                 expDensity_abs = 2., expTemp_abs = -3.5, \
                 expDensity_scat = 1., expTemp_scat = 0., \
                 dxset = 5, numlev = 1, runtime = 5.e-7, \
                 freezeWidth = 10., dumpnum = 100, \
                 problem = 'JMF', \
#                  make_RH_solution_bool = True, \
                 make_RT_solution_override = True, \
                 eps = 1.e-6, int_tol = 1.e-10, integrator = 'ode', \
                 print_the_sources = True):
        init_RadShock.__init__(self, M0 = M0, rho0 = rho0, Sn = Sn, \
                               Tref = Tref, Cv = Cv, gamma = gamma, \
                               eps = eps, sigA = sigA, sigS = sigS, \
                               expDensity_abs = expDensity_abs, \
                               expTemp_abs = expTemp_abs, \
                               expDensity_scat = expDensity_scat, \
                               expTemp_scat = expTemp_scat, \
                               dxset = dxset, numlev = numlev,
                               runtime = runtime, freezeWidth = freezeWidth,\
                               dumpnum = dumpnum, problem = problem, \
#                                make_RH_solution_bool = make_RH_solution_bool, \
                               make_RT_solution_override = make_RT_solution_override, \
                               int_tol = int_tol, integrator = integrator)
        self.nED_profile = utils.AnalyticShockProfiles(self)
#         self.nED_profile.make_RH_solution_bool = make_RH_solution_bool
        self.nED_profile.print_the_sources = print_the_sources
        self.nED_profile.downstream_equilibrium()
        self.nED_profile.make_RH_solution()
        self.nED_profile.write_data()

class greySn_RadShock(greyNED_RadShock):
    '''
    Define the grey-Sn radiative-shock problem, and drive the solution
    '''
    def __init__(self, M0 = 1.2, rho0 = 1., Sn = 16, \
                 Tref = 100., Cv = 1.4472799784454e12, gamma = 5. / 3., \
#                  sigA = 44.87940546534333, sigS = 0.4006204354152789, \
#                  expDensity_abs = 2., expTemp_abs = -3.5, \
#                  expDensity_scat = 1., expTemp_scat = 0., \
                 sigA = 577.35, sigS = 0., \
                 expDensity_abs = 0., expTemp_abs = 0., \
                 expDensity_scat = 0., expTemp_scat = 0., \
                 dxset = 5, numlev = 1, runtime = 5.e-7, \
                 freezeWidth = 10., dumpnum = 100, \
                 make_RT_solution_override = True, \
                 problem = 'JMF', f_tol = 1.e-6, \
                 eps = 1.e-6, int_tol = 1.e-10, integrator = 'ode', \
                 print_the_sources = True):
        self.Sn = Sn
        greyNED_RadShock.__init__(self, M0 = M0, rho0 = rho0, Sn = Sn, \
                                  Tref = Tref, Cv = Cv, gamma = gamma, \
                                  eps = eps, sigA = sigA, sigS = sigS, \
                                  expDensity_abs = expDensity_abs, \
                                  expTemp_abs = expTemp_abs, \
                                  expDensity_scat = expDensity_scat, \
                                  expTemp_scat = expTemp_scat, \
                                  dxset = dxset, numlev = numlev,
                                  runtime = runtime, freezeWidth = freezeWidth,\
                                  dumpnum = dumpnum, problem = problem, \
                                  make_RT_solution_override = make_RT_solution_override, \
                                  int_tol = int_tol, integrator = integrator)
        self.Sn_profile = copy.deepcopy(self.nED_profile)
        self.Sn_profile.Sn = Sn
        self.Sn_profile.f_tol = f_tol
        self.Sn_profile.make_RT_solution_override = make_RT_solution_override
        self.Sn_profile.print_the_sources = print_the_sources
        self.Sn_profile.make_dictionaries()
        self.Sn_profile.make_RT_solution()
        while self.Sn_profile.make_RT_solution_bool:
            self.Sn_profile.f_iters += 1
            self.Sn_profile.make_RH_solution()
            self.Sn_profile.make_RT_solution()
        self.Sn_profile.append_metadata()

class polar_plot(init_RadShock):
    def __init__(self, M0 = 1.2, rho0 = 1., Sn = 16, \
                 Tref = 100., Cv = 1.4472799784454e12, gamma = 5. / 3., \
#                  sigA = 577.35, sigS = 0., \
#                  expDensity_abs = 0., expTemp_abs = 0., \
#                  expDensity_scat = 0., expTemp_scat = 0., \
                 sigA = 44.87940546534333, sigS = 0.4006204354152789, \
                 expDensity_abs = 2., expTemp_abs = -3.5, \
                 expDensity_scat = 1., expTemp_scat = 0., \
                 dxset = 5, numlev = 1, runtime = 5.e-7, \
                 freezeWidth = 10., dumpnum = 100, \
                 make_new_solution = True, problem = 'JMF', f_tol = 1.e-6, \
                 eps = 1.e-6, int_tol = 1.e-10, integrator = 'ode', \
                 plot = 'polar'):
        init_RadShock.__init__(self, M0 = M0, rho0 = rho0, Sn = Sn, \
                                  Tref = Tref, Cv = Cv, gamma = gamma, \
                                  eps = eps, sigA = sigA, sigS = sigS, \
                                  expDensity_abs = expDensity_abs, \
                                  expTemp_abs = expTemp_abs, \
                                  expDensity_scat = expDensity_scat, \
                                  expTemp_scat = expTemp_scat, \
                                  dxset = dxset, numlev = numlev,
                                  runtime = runtime, freezeWidth = freezeWidth,\
                                  dumpnum = dumpnum, problem = problem, \
                                  make_new_solution = make_new_solution, \
                                  int_tol = int_tol, integrator = integrator)
        def latexify(fig_width = None, fig_height = None, columns = 1):
            assert(columns in [1,2])
            if fig_width is None:
                fig_width = 3.39 if columns == 1 else 6.9
            if fig_height is None:
                golden_mean = (numpy.sqrt(5.) - 1.) / 2.
                fig_height = fig_width * golden_mean
            MAX_HEIGHT_INCHES = 8.
            if fig_height > MAX_HEIGHT_INCHES:
                print_stmnt  = "WARNING: fig_height too large: " + str(fig_height)
                print_stmnt += "so will reduce to " + str(MAX_HEIGHT_INCHES) + "inches."
                print print_stmnt
                fig_height = MAX_HEIGHT_INCHES
            params = {'backend': 'ps',
                      'text.latex.preamble':['\usepackage{gensymb}'],
                      'axes.labelsize':8,
                      'axes.titlesize':8,
                      'text.fontsize':8,
                      'legend.fontsize':8,
                      'xtick.labelsize':8,
                      'ytick.labelsize':8,
                      'text.usetex':True,
                      'figure.figsize':[fig_width, fig_height],
                      'font.family':'serif'}
            matplotlib.rcParams.update(params)
        latexify()
        file_name = '/home/jmferguson/Documents/code-verification/2Trad_shocks/'
        file_name += 'grey_Sn_transport/OOP/TempDensity_crosssections/M3.0/'
        file_name += 'data_dictionary_fiter7.pickle'
        open_file = open(file_name, 'r')
        self.data = pickle.load(open_file)
        open_file.close()
        self.M0 = self.data['d_M'][0][0]
        self.gamma = 5./3.
        self.Sn = 2 * numpy.shape(self.data['d_Im_fwd'][0])[0]
        self.c = 2.99792458e10
        self.ar = 137.20172
        self.Cv = 1.4472799784454e12
        self.rho0 = 1.
        self.Tref = 100.
        val = numpy.sqrt(self.gamma * (self.gamma - 1.) * self.Cv * self.Tref)
        self.sound_speed = val
        self.P0 = self.ar * self.Tref**4 / self.rho0 / self.sound_speed**2
        self.C0 = self.c / self.sound_speed
        self.mus_16 = numpy.loadtxt('mus_' + str(self.Sn) + '.txt')
        self.weights_16 = numpy.loadtxt('weights_' + str(self.Sn) + '.txt')
        self.Ert = self.data['d_Ert']
        self.Frt = self.data['d_Frt']
        self.Prt = self.data['d_Prt']
        self.f = self.data['d_f'][-1]
        self.M = self.data['d_M']
        self.x = self.data['d_x']
        self.Fr_RT = self.Frt[-1]
        self.Pr_RT = self.Prt[-1]
        self.Mach_RT = self.M[-1]
        self.x_RT = self.x[-1]
        self.Temp_RT = fnctn.mat_temp(self.Pr_RT, self.Mach_RT, self)
        self.radTemp_RT = fnctn.rad_temp(self.Pr_RT, self.Mach_RT, self)
        self.Density_RT = fnctn.mat_density(self.Pr_RT, self.Mach_RT, self)
        self.Pressure_RT = fnctn.mat_pres(self.Pr_RT, self.Mach_RT, self)
        self.Im_fwd = self.data['d_Im_fwd'][-1]
        self.Im_rev = self.data['d_Im_rev'][-1]
        self.Pr0 = self.Pr_RT[0]
        if (plot == 'polar'):
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
        elif(plot == 'compare_LE'):
            print 'entered compare_LE plots'
        elif(plot == 'Sns'):
            print 'entered create Sn plots'
