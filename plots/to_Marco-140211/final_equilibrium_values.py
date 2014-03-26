# DETERMINE FINAL EQUILIBRIUM VALUES
guessTemperature = ((18. * M0 ** 2. / P0 - 1.) / 7.) ** (.25)
eq1 = lambda T: 3. * (gamma + 1.) * (T - 1.) - P0 * gamma * (gamma - 1.) * (T ** 4 + 7.)
eq2 = lambda T: 12. * (gamma - 1.) ** 2 * T * (3. + gamma * P0 * (1. + 7. * T ** 4))
# eq3 is the  equilibrium density
eq3 = lambda T: (eq1(T) + (eq1(T) ** 2 + eq2(T)) ** (.5)) / (6. * (gamma - 1.) * T)
# eqRoot determines the final equilibrium temperature
eqRoot = lambda T: 3. * eq3(T) * (eq3(T) * T - 1.) + gamma * P0 * eq3(T) * (T ** 4 - 1.) - 3. * gamma * (eq3(T) - 1.) * M0 ** 2
finalTemperature = scipy.optimize.bisect(eqRoot, initialTemperature + 1e-15, guessTemperature)
finalDensity = eq3(finalTemperature)
finalMach = M0 / finalDensity / (finalTemperature ** (.5))
finalSpeed = finalMach * (finalTemperature ** (.5))

Mach_left = numpy.linspace(M0, 1, points, endpoint = False)
fRED_left = numpy.zeros(2)
fRED_left[0] = initialTemperature**4 / 3.
Mach_right = numpy.linspace(finalMach, 1, points, endpoint = False)
fRED_right = numpy.zeros(2)
fRED_right[0] = finalTemperature**4 / 3.
