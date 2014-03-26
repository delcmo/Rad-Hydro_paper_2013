left1 = 1
right1 = 1
hold_left_right = 0
Mach_left = numpy.linspace(M0, 1, points, endpoint = False)
Mach_right = numpy.linspace(finalMach, 1, points, endpoint = False)
int_deltaM_left = 1
int_deltaM_right = 1

while (left1 != 0) & (right1 != 0):
  if hold_left_right == 1:
    Mach_left = numpy.linspace(M0, Mach_left1[-1], points, endpoint = True)
    Mach_right = numpy.linspace(finalMach, Mach_right1[-1], points, endpoint = True) 
# delta_M_left = Mach_left[1] - Mach_left[0]
# delta_M_right = Mach_right[1] - Mach_right[0]
# execfile('initialization.py')
# check_left = scipy.integrate.odeint(DfRED_DM, fRED_left[0] + delta_fRED_left, Mach_left[1:50], mxstep = 5000, atol = int_atol, rtol = int_rtol)
# check_left = scipy.append(fRED_left[0], check_left)
# while sum(numpy.diff(check_left) < 0) != 0:
#   int_deltaM_left += 1
#   if hold_left_right == 0:
#     Mach_left = numpy.linspace(M0 + delta_M_left * int_deltaM_left, 1, points - 1, endpoint = False)
#   else:
#     Mach_left = numpy.linspace(M0 + delta_M_left * int_deltaM_left, Mach_left1[-1], points - 1, endpoint = True)
#   Mach_left = scipy.append(M0, Mach_left)
#   execfile('initialization.py')
#   check_left = scipy.integrate.odeint(DfRED_DM, fRED_left[0] + delta_fRED_left, Mach_left[1:50], mxstep = 5000, atol = int_atol, rtol = int_rtol)
#   check_left = scipy.append(fRED_left[0], check_left)
  a = scipy.integrate.odeint(DfRED_DM, 2. * numpy.pi * fnum_interp(Mach_left[1]), Mach_left[1:], mxstep = 5000, atol = int_atol, rtol = int_rtol)
  fRED_left = scipy.append(fRED_left[0], a)
  print 'LEFT INTEGRATION COMPLETE'
#  check_right = scipy.integrate.odeint(DfRED_DM, fRED_right[0] + delta_fRED_right, Mach_right[1:50], mxstep = 5000, atol = int_atol, rtol = int_rtol)
#  check_right = scipy.append(fRED_right[0], check_right)
#  while sum(numpy.diff(check_right) > 0) != 0:
#    int_deltaM_right += 1
#    if hold_left_right == 0:
#      Mach_right = numpy.linspace(finalMach + delta_M_right * int_deltaM_right, 1, points - 1, endpoint = False)
#    else:
#      Mach_right = numpy.linspace(finalMach + delta_M_right * int_deltaM_right, Mach_right1[-1], points - 1, endpoint = True)
#    Mach_right = scipy.append(finalMach, Mach_right)
#    execfile('initialization.py')
#    check_right = scipy.integrate.odeint(DfRED_DM, fRED_right[0] + delta_fRED_right, Mach_right[1:50], mxstep = 5000, atol = int_atol, rtol = int_rtol)
#    check_right = scipy.append(fRED_right[0], check_right)
  print 'RIGHT INTEGRATION COMPLETE'
  b = scipy.integrate.odeint(DfRED_DM, 2. * numpy.pi * fnum_interp(Mach_right[1]), Mach_right[1:], mxstep = 5000, atol = int_atol, rtol = int_rtol)
  fRED_right = scipy.append(fRED_right[0], b)
  diff_fRED_left = numpy.diff(fRED_left)
  val = 0
  k_val = 0
  if int(sum(diff_fRED_left < 0)) != 0:
    val = 1
  if val > 0:
    for val_iter in diff_fRED_left:
      val = val_iter
      k_val += 1
      if val < 0:
        break
  if k_val > 0:
    k_val -= 1
    Mach_left_ext = numpy.linspace(Mach_left[k_val], 1, 10 * (len(fRED_left) - k_val), endpoint = False)
    fRED_left_ext = scipy.integrate.odeint(DfRED_DM, fRED_left[k_val], Mach_left_ext, mxstep = 5000, atol = int_atol, rtol = int_atol)
    fRED_left = scipy.append(fRED_left[:k_val], fRED_left_ext)
    Mach_left = scipy.append(Mach_left[:k_val], Mach_left_ext)
RED_left = fRED_left / f_interp(Mach_left)
RED_left[0] = initialTemperature**4
RED_right = fRED_right / f_interp(Mach_right)
RED_right[0] = finalTemperature**4
  right = 0
  if min(RED_right) < 1:
    while RED_right[right] > 1:
      right = right + 1
  if right > 0:
    RED_right = RED_right[:right]
    Mach_right = Mach_right[:right]
left1 = 0
right1 = 0
if RED_left[-1]  > RED_right[-1]:
   while RED_left[-1 - left1] > RED_right[-1 - right1]:
      left1 = left1 + 1
#  end while
   flux_sign = rad_flux(left1, right1)
   while math.copysign(1, rad_flux(left1, right1)) == math.copysign(1, flux_sign):
      right1 = right1 + 1
      while int(RED_left[-1 - left1] < RED_right[-1 - right1]) & int(left1 > 0):
            left1 = left1 - 1
#     end while
      left1 = left1 + 1
      if right1 == len(RED_right) - 1:
         right1 = 0
         break
      if left1 == points - 1:
         left1 = left1 - 2
         break
  #     end if
  #  end while
if left1 == 0:
  Mach_left1 = Mach_left
  RED_left1 = RED_left
  fRED_left1 = fRED_left
else:
  Mach_left1 = Mach_left[:- left1]
  RED_left1 = RED_left[: -left1]
  fRED_left1 = fRED_left[: -left1]
if right1 == 0:
  Mach_right1 = Mach_right
  RED_right1 = RED_right
  fRED_right1 = fRED_right
else:
  Mach_right1 = Mach_right[:- right1]
  RED_right1 = RED_right[:- right1]
  fRED_right1 = fRED_right[:- right1]
hold_left_right = 1 
# end while

print 'int_deltaM_left =', int_deltaM_left
print 'int_deltaM_right =', int_deltaM_right

x_left = scipy.integrate.odeint(dx_dM_left, 0, Mach_left1[1:], mxstep = 5000, atol = int_atol, rtol = int_rtol)
x_left1 = x_left - x_left[-1]

x_right = scipy.integrate.odeint(dx_dM_right, 0, Mach_right1[1:], mxstep = 5000, atol = int_atol, rtol = int_rtol)
x_right1 = x_right - x_right[-1]

Mach = scipy.append(Mach_left1, Mach_right1[::-1])

RED = scipy.append(RED_left1, RED_right1[::-1])
RED[0] = initialTemperature**4
RED[-1] = finalTemperature**4

radTemp = RED**(1./4.)
radTemp[0] = initialTemperature
radTemp[-1] = finalTemperature

fRED = scipy.append(fRED_left1, fRED_right1[::-1])
fRED[0] = RED[0] / 3.
fRED[-1] = RED[-1] / 3.

Density = fnctnDensity(Mach, fRED)
Density[0] = initialDensity
Density[-1] = finalDensity

Speed = fnctnSpeed(Mach, fRED)
Speed[0] = initialSpeed
Speed[-1] = finalSpeed

radFlux = fnctnradFlux(Mach, fRED)

Temp = fnctnTemp(Mach, fRED)
Temp[0] = initialTemperature
Temp[-1] = finalTemperature
