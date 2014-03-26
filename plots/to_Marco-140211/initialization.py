# print 'left initialization'
M = M0
P = P0 / 3.
vel = initialSpeed
fred = initialTemperature / 3.
D1 = M**2 - 1
Dg = gamma * M**2 - 1

a = Cp * M0 * D1 / P
b = 12. * Dg * sigma_a(M, fred)
c = M0
d = Dg * sigma_a(M, fred) * (RED_interp(Mach_left[1]) - RED_interp(Mach_left[0])) / (fRED_interp(Mach_left[1]) - fRED_interp(Mach_left[0]))
e = Cp - (M**2 + 4. * P / (1. - 2. * vel**2 / C0**2 * sigma_a(M, fred) / sigma_t(M, fred))) / Dg + 12. * P * sigma_a(M, fred) / (sigma_t(M, fred) - 2. * vel**2 / C0**2 * sigma_a(M, fred))
g = P * C0**2 / M0 / (sigma_t(M, fred) - 2. * vel**2 / C0**2 * sigma_a(M, fred))
h = gamma * P * (M**2 + 4. * P / (1. - 2. * vel**2 / C0**2 * sigma_a(M, fred) / sigma_t(M, fred))) / Dg - P * (sigma_t(M, fred) + sigma_s(M, fred) * (RED_interp(Mach_left[1]) - RED_interp(Mach_left[0])) / (fRED_interp(Mach_left[1]) - fRED_interp(Mach_left[0]))) / (sigma_t(M, fred) - 2. * vel**2 / C0**2 * sigma_a(M, fred))

# alpha_plus_left = - ((a * h + b * g - e * c) + numpy.sqrt((e * c - a * h - b * g)**2 + 4. * a * g * (e * d - b * h))) / 2. / a / g

alpha_minus_left = - ((a * h + b * g - e * c) - numpy.sqrt(abs((e * c - a * h - b * g)**2 + 4. * a * g * (e * d - b * h)))) / 2. / a / g
delta_fRED_left = - (Mach_left[1] - Mach_left[0]) / M0 * fred * Dg / ((gamma * M**2 + 1.) * (c * alpha_minus_left + d) / (a * alpha_minus_left + b) / 2. + gamma * P)

# print 'a =', a
# print 'b =', b
# print 'c =', c
# print 'd =', d
# print 'e =', e
# print 'g =', g
# print 'h =', h
# print 'alpha_minus_left =', alpha_minus_left
# print 'delta_fRED_left =', delta_fRED_left

# print 'right initialization'
M = finalMach
P = P0 * finalTemperature**3 / finalDensity / 3.
vel = finalSpeed
fred = finalTemperature**4 / 3.
D1 = M**2 - 1
Dg = gamma * M**2 - 1

a = Cp * M0 * D1 / P
b = 12. * Dg * sigma_a(M, fred)
c = M0
d = Dg * sigma_a(M, fred) * (RED_interp(Mach_right[1]) - RED_interp(Mach_right[0])) / (fRED_interp(Mach_right[1]) - fRED_interp(Mach_right[0]))
e = Cp - (M**2 + 4. * P / (1. - 2. * vel**2 / C0**2 * sigma_a(M, fred) / sigma_t(M, fred))) / Dg + 12. * P * sigma_a(M, fred) / (sigma_t(M, fred) - 2. * vel**2 / C0**2 * sigma_a(M, fred))
g = P * C0**2 / M0 / (sigma_t(M, fred) - 2. * vel**2 / C0**2 * sigma_a(M, fred))
h = gamma * P * (M**2 + 4. * P / (1. - 2. * vel**2 / C0**2 * sigma_a(M, fred) / sigma_t(M, fred))) / Dg - P * (sigma_t(M, fred) + sigma_s(M, fred) * (RED_interp(Mach_right[1]) - RED_interp(Mach_right[0])) / (fRED_interp(Mach_right[1]) - fRED_interp(Mach_right[0]))) / (sigma_t(M, fred) - 2. * vel**2 / C0**2 * sigma_a(M, fred))

# alpha_plus_right = - ((a * h + b * g - e * c) + numpy.sqrt((e * c - a * h - b * g)**2 + 4. * a * g * (e * d - b * h))) / 2. / a / g

alpha_minus_right = - ((a * h + b * g - e * c) - numpy.sqrt(abs((e * c - a * h - b * g)**2 + 4. * a * g * (e * d - b * h)))) / 2. / a / g
delta_fRED_right = - (Mach_right[1] - Mach_right[0]) / M0 * fred * Dg / ((gamma * M**2 + 1.) * (c * alpha_minus_right + d) / (a * alpha_minus_right + b) / 2. + gamma * P)

# print 'a =', a
# print 'b =', b
# print 'c =', c
# print 'd =', d
# print 'e =', e
# print 'g =', g
# print 'h =', h
# print 'alpha_minus_right =', alpha_minus_right
# print 'delta_fRED_right =', delta_fRED_right
