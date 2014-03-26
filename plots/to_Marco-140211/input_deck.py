# initialMach = 1.05
x0 = -.02
x1 = .02
sigt = 577.35
siga = 577.35

"""  T O L E R A N C E S  """
f_tol = 1.e-6
RED_tol = 1.e-5
fRED_tol = 1.e-5
int_atol = 1.e-10
int_rtol = 1.e-10

"""  I N P U T S  """""
Sn = 16
weights = numpy.loadtxt('weights_16.txt')
mus = numpy.loadtxt('mus_16.txt')
# execfile('Gauss-Legendre_weights_nodes.py')
# [weights, mus, err] = GaussLegendreWeights(Sn)
# mus = mus[::-1]

"""  D I M E N S I O N A L         C O N S T A N T S    """
gamma = 5. / 3
R = 8.3145 #[Joules /  mole / Kelvin]

"""  D I M E N S I O N A L         Q U A N T I T I E S  """
L_ref = 1 
T_inf = 300 #[Kelvin]
M_air = 0.0289645 #[kg / mole], the mean molar mass of air
# for Lowrie Edwards u_inf = 173205 m/s
u_inf = 173205
# u_inf = numpy.sqrt(gamma * R * T_inf / M_air) #[meters / sec], the speed of sound used as a reference velocity
                                              #                the speed of sound in air at 293 Kelvin is 343.2 [meters / sec]

"""  D I M E N S I O N L E S S     Q U A N T I T I E S  """

points = 1e3 + 1
# C0 = 3.e8 / u_inf
C0 = numpy.sqrt(3.e6)
initialRadratio = 1.e-4
initialTemperature = 1.
initialDensity = 1.
initialSpeed = initialMach
