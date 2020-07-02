import hepak
import he3pak
import math
import scipy.integrate
import scipy.optimize
import scipy.special
import numpy
import matplotlib.pyplot as plt
import CoolProp.CoolProp

def CuThermalConductivity(T):
  # NIST cryogenic material properties, RRR = 50 [W/m K]
  if not 4. <= T <= 301.:
    print('Evaluating thermal conductivity of Cu outside valid range ({0} K)!'.format(T))
  a = 1.8743
  b = -0.41538
  c = -0.6018
  d = 0.13294
  e = 0.26426
  f = -0.0219
  g = -0.051276
  h = 0.0014871
  i = 0.003723
  logk = (a + c * T**0.5 + e * T + g * T**1.5 + i * T**2) / (1. + b * T**0.5 + d * T + f * T**1.5 + h * T**2)
  return 10**logk
  
def CuSpecificHeat(T):
  # NIST cryogenic material properties [J/kg K]
  if not 4. <= T <= 300.:
    print('Evaluating specific heat of Cu outside valid range!')
  a = -1.91844
  b = -0.15973
  c = 8.61013
  d = -18.996
  e = 21.9661
  f = -12.7328
  g = 3.54322
  h = -0.3797
  i = 0.
  logT = math.log10(T)
  logCP = a + b * logT + c * logT**2 + d * logT**3 + e * logT**4 + f * logT**5 + g * logT**6 + e * logT**7 + i * logT**8
  return 10**logCP

# calculate velocity of He4, He3, or N2 gas flow through annulus
def velocity(massFlow, ID, OD, T, P, gas):
  if gas == 'He4':
    density = hepak.HeCalc('D', 0, 'T', T, 'P', P, 1)
  elif gas == 'He3':
    density = he3pak.He3Density(P, T)
  elif gas == 'N2':
    density = CoolProp.CoolProp.PropsSI('D', 'T', T, 'P', P, 'Nitrogen')
  elif gas == 'He':
    density = CoolProp.CoolProp.PropsSI('D', 'T', T, 'P', P, 'Helium')
  else:
    raise NameError('Unknow gas ' + gas + '!')
  crossSection = (OD**2 - ID**2)*math.pi/4
  return massFlow / density / crossSection
  
# Reynolds number of He4, He3, or N2 gas flowing through annulus
def reynoldsNumber(massFlow, ID, OD, T, P, gas):
  crossSection = (OD**2 - ID**2) * math.pi / 4.
  hydraulicDiameter = OD - ID
  if gas == 'He4':
    viscosity = hepak.HeCalc(25, 0, 'T', T, 'P', P, 1)
  elif gas == 'He3':
    density = he3pak.He3Density(P, T)
    viscosity = he3pak.He3Prop(29, density, T)
  elif gas == 'N2':
    viscosity = CoolProp.CoolProp.PropsSI('viscosity', 'T', T, 'P', P, 'Nitrogen')
  elif gas == 'He':
    viscosity = CoolProp.CoolProp.PropsSI('viscosity', 'T', T, 'P', P, 'Helium')
  else:
    raise NameError('Unknow gas ' + gas + '!')
  return massFlow * hydraulicDiameter / viscosity / crossSection  

# Darcy friction formula for gas flowing through annulus with wall roughness
def frictionFactor(ID, OD, roughness, reynoldsNumber):
#  if reynoldsNumber < 1000:
#    return 64./reynoldsNumber
#  a = 2.51 / reynoldsNumber
#  b = roughness/(OD - ID)/3.7
#  f = (2. * scipy.special.lambertw(math.log(10.)/2/a * 10.**(b/2*a)).real / math.log(10.) - b/a)**-2 # Colebrook-White equation
#  return f
  a = 1./(1. + (reynoldsNumber/2712.)**8.4)
  b = 1./(1. + (reynoldsNumber/150/(OD - ID)*roughness)**1.8)
  return (64./reynoldsNumber)**a * (0.75*math.log(reynoldsNumber/5.37))**(2*(a-1)*b) * (0.88*math.log(3.41*(OD-ID)/roughness)**(2*(a-1)*(1-b))) # Bellos, Nalbantis, Tsakiris (Wikipedia)


# Nusselt number for He4, He3, or N2 gas flowing through annulus with wall temperature T_wall
def nusseltNumber(massFlow, ID, OD, T, T_wall, P, roughness, gas):
  Re = reynoldsNumber(massFlow, ID, OD, T, P, gas)
  
  if Re < 1000: # laminar flow
    return 4.36
  
  if gas == 'He4':
    prandtlNumber = hepak.HeCalc(27, 0, 'T', T, 'P', P, 1)
    viscosity = hepak.HeCalc(25, 0, 'T', T, 'P', P, 1)
    viscosity_wall = hepak.HeCalc(25, 0, 'T', T_wall, 'P', P, 1)
  elif gas == 'He3':
    density = he3pak.He3Density(P, T)
    prandtlNumber = he3pak.He3Prop(32, density, T)
    viscosity = he3pak.He3Prop(29, density, T)
    viscosity_wall = he3pak.He3Prop(29, density, T_wall)
  elif gas == 'N2':
    prandtlNumber = CoolProp.CoolProp.PropsSI('Prandtl', 'T', T, 'P', P, 'Nitrogen')
    viscosity = CoolProp.CoolProp.PropsSI('viscosity', 'T', T, 'P', P, 'Nitrogen')
    viscosity_wall = CoolProp.CoolProp.PropsSI('viscosity', 'T', T_wall, 'P', P, 'Nitrogen')
  elif gas == 'He':
    prandtlNumber = CoolProp.CoolProp.PropsSI('Prandtl', 'T', T, 'P', P, 'Helium')
    viscosity = CoolProp.CoolProp.PropsSI('viscosity', 'T', T, 'P', P, 'Helium')
    viscosity_wall = CoolProp.CoolProp.PropsSI('viscosity', 'T', T_wall, 'P', P, 'Helium')
  else:
    raise NameError('Unknow gas ' + gas + '!')
  
#  if Re > 10000.:
#  if not 0.7 <= prandtlNumber <= 16700.:
#    print('Prandtl number {0} outside valid range for Sieder-Tate correlation'.format(prandtlNumber))
  return 0.023 * Re**0.8 * prandtlNumber**(1./3.) * (viscosity/viscosity_wall)**0.14 # Sieder-Tate correlation
  
#  if not 0.5 <= prandtlNumber <= 2000.:
#    print('Prandtl number {0} outside valid range for Gnielinski correlation'.format(prandtlNumber))
#  f = frictionFactor(ID, OD, roughness, Re)
#  return (f/8.)*(Re - 1000.)*prandtlNumber / (1 + 12.7*(f/8.)**0.5 * (prandtlNumber**(2./3.) - 1)) # Gnielinski correlation
  

# convective heat transfer coefficient for He4, He3, or N2 gas flowing through annulus with wall temperature T_wall
def heatTransferCoeff(massFlow, ID, OD, T, T_wall, P, roughness, gas):
  Nu = nusseltNumber(massFlow, ID, OD, T, T_wall, P, roughness, gas)
  thermalConductivity = 0.
  if gas == 'He4':
    thermalConductivity = hepak.HeCalc(26, 0, 'T', T, 'P', P, 1)
  elif gas == 'He3':
    density = he3pak.He3Density(P, T)
    thermalConductivity = he3pak.He3Prop(28, density, T)
  elif gas == 'N2':
    thermalConductivity = CoolProp.CoolProp.PropsSI('conductivity', 'T', T, 'P', P, 'Nitrogen')
  elif gas == 'He':
    thermalConductivity = CoolProp.CoolProp.PropsSI('conductivity', 'T', T, 'P', P, 'Helium')
  else:
    raise NameError('Unknow gas ' + gas + '!')
  return Nu * thermalConductivity / (OD - ID)

# heat per length transferred into He4, He3, or N2 gas flowing through annulus
# interfaceLength is the wetted perimeter with temperature T_wall
def dQdx(massFlow, ID, OD, interfaceLength, T, T_wall, P, roughness, gas):
  h = heatTransferCoeff(massFlow, ID, OD, T, T_wall, P, roughness, gas)
  return h * interfaceLength * (T_wall - T)

# heat per length transferred into He4, He3, or N2 gas flowing through annulus with different temperatures for inner and outer wall
def dQdx2(massFlow, ID, OD, T, T_innerWall, T_outerWall, P, roughness, gas):
  hInner = heatTransferCoeff(massFlow, ID, OD, T, T_innerWall, P, roughness, gas)
  hOuter = heatTransferCoeff(massFlow, ID, OD, T, T_outerWall, P, roughness, gas)
  return hInner * ID*math.pi * (T_innerWall - T) + hOuter * OD*math.pi * (T_outerWall - T)
  
# heat per length transferred through copper pipe wall
def dQCudx(T1, T2, diameter, thickness):
  return CuThermalConductivity((T1 + T2)/2.) * (T2 - T1) / thickness * diameter*math.pi
  
# flow of N2 evaporated by heat load from a flow of warm He4 or He3 gas
def LN2evaporationFlow(vaporPressure, massFlow, T_incoming, P_incoming, incomingGas):
  T = CoolProp.CoolProp.PropsSI('T', 'P', vaporPressure, 'Q', 0, 'Nitrogen')
  latentHeat = CoolProp.CoolProp.PropsSI('H', 'P', vaporPressure, 'Q', 1, 'Nitrogen') - CoolProp.CoolProp.PropsSI('H', 'P', vaporPressure, 'Q', 0, 'Nitrogen')
  if incomingGas == 'He4':
    enthalpyDifference = hepak.HeCalc('H', 0, 'T', T_incoming, 'P', P_incoming, 1) - hepak.HeCalc('H', 0, 'T', T, 'P', P_incoming, 1)
  elif incomingGas == 'He3':
    incomingDensity = he3pak.He3Density(P_incoming, T_incoming)
    coldDensity = he3pak.He3Density(P_incoming, T)
    enthalpyDifference = he3pak.He3Prop(6, incomingDensity, T_incoming) - he3pak.He3Prop(6, coldDensity, T)
  elif incomingGas == 'He':
    enthalpyDifference = CoolProp.CoolProp.PropsSI('enthalpy', 'T', T_incoming, 'P', P_incoming, 'Helium') - CoolProp.CoolProp.PropsSI('enthalpy', 'T', T, 'P', P_incoming, 'Helium')
  else:
    raise NameError('Unknow gas ' + gas + '!')
  heatLoad = massFlow*enthalpyDifference
  return heatLoad / latentHeat
  

# pressure drop per length of He4 or He3 gas flowing through annulus with wall roughness
def dPdx(massFlow, ID, OD, roughness, T, P, gas):
  Re = reynoldsNumber(massFlow, ID, OD, T, P, gas)
  f = frictionFactor(ID, OD, roughness, Re)
  v = velocity(massFlow, ID, OD, T, P, gas)
  crossSection = (OD**2 - ID**2)*math.pi/4
  dPdx = f * 2.*massFlow/crossSection/math.pi * v/(OD - ID)
  return dPdx

# temperature gradient of He4, He3, or N2 gas flowing through annulus
# interfaceLength is the wetted perimeter with temperature T_wall
def dTdx(massFlow, ID, OD, interfaceLength, T, T_wall, P, roughness, gas):
  if gas == 'He4':
    cP = hepak.HeCalc(14, 0, 'T', T, 'P', P, 1)
  elif gas == 'He3':
    density = he3pak.He3Density(P, T)
    cP = he3pak.He3Prop(9, density, T)
  elif gas == 'N2':
    cP = CoolProp.CoolProp.PropsSI('Cpmass', 'T', T, 'P', P, 'Nitrogen')
  elif gas == 'He':
    cP = CoolProp.CoolProp.PropsSI('Cpmass', 'T', T, 'P', P, 'Helium')
  else:
    raise NameError('Unknow gas ' + gas + '!')
  return dQdx(massFlow, ID, OD, interfaceLength, T, T_wall, P, roughness, gas) / cP / massFlow
  
# temperature gradient of He4, He3, or N2 gas flowing through annulus with different temperatures for inner and outer wall
def dTdx2(massFlow, ID, OD, T, T_innerWall, T_outerWall, P, roughness, gas):
  if gas == 'He4':
    cP = hepak.HeCalc(14, 0, 'T', T, 'P', P, 1)
  elif gas == 'He3':
    density = he3pak.He3Density(P, T)
    cP = he3pak.He3Prop(9, density, T)
  elif gas == 'N2':
    cP = CoolProp.CoolProp.PropsSI('Cpmass', 'T', T, 'P', P, 'Nitrogen')
  elif gas == 'He':
    cP = CoolProp.CoolProp.PropsSI('Cpmass', 'T', T, 'P', P, 'Helium')
  else:
    raise NameError('Unknow gas ' + gas + '!')
  return dQdx2(massFlow, ID, OD, T, T_innerWall, T_outerWall, P, roughness, gas) / cP / massFlow  
  
# temperature of thin wall between annular He4, He3, or N2 gas flows (set wall temperature so heat transfer on both sides is equal)
def Twall(innerMassFlow, innerID, innerOD, T_inner, P_inner, innerGas, outerMassFlow, outerID, outerOD, T_outer, P_outer, outerGas, roughness):
  diameter = (innerOD + outerID)/2.
  thickness = outerID - innerOD
  fun = lambda T: dQdx(innerMassFlow, innerID, innerOD, innerOD*math.pi, T_inner, T, P_inner, roughness, innerGas) \
                + dQdx(outerMassFlow, outerID, outerOD, outerID*math.pi, T_outer, T, P_outer, roughness, outerGas)
  T_wall = scipy.optimize.root_scalar(fun, bracket = [T_inner, T_outer])
  if not T_wall.converged:
    print('Could not calculate wall temperature!')
  return T_wall.root
  
  
# plot various gas properties along the length of the heat exchanger
def plotTubeInTubeHEX(sol, innerMassFlow, innerID, innerOD, innerGas, outerMassFlow, outerID, outerOD, outerGas, N2ID, N2OD, roughness):
  fig, axes = plt.subplots(3, 2, figsize = (9.6, 10.8))
  fig.set_tight_layout(True)
  axes[0][0].plot(sol.x, sol.y[0], color = 'tab:red')
  axes[0][0].plot(sol.x, sol.y[1], color = 'tab:blue')
  y = numpy.array([Twall(innerMassFlow, innerID, innerOD, T_He, P_He, innerGas, outerMassFlow, outerID, outerOD, T_He3, P_He3, outerGas, roughness) \
                   for T_He, T_He3, P_He, P_He3 in zip(sol.y[0], sol.y[1], sol.y[2], sol.y[3])])
  axes[0][0].plot(sol.x, y, color = 'tab:brown')
  if len(sol.y) > 4:
    N2massFlow = LN2evaporationFlow(sol.y[5][0], He3MassFlow, sol.y[1][0], sol.y[3][0], outerGas)
    axes[0][0].plot(sol.x, sol.y[4], color = 'tab:green')    
    y = numpy.array([Twall(innerMassFlow, innerID, innerOD, T_He, P_He, outerGas, N2massFlow, 0., N2OD, T_N2, P_N2, 'N2', roughness) \
                     for T_N2, T_He, P_N2, P_He in zip(sol.y[0], sol.y[4], sol.y[2], sol.y[5])])
    axes[0][0].plot(sol.x, y, color = 'tab:purple')
  axes[0][0].set_yscale('log')
  axes[0][0].set_ylabel('Temperature (K)')
  axes[0][0].grid(True, 'both')
  
  y = numpy.array([reynoldsNumber(innerMassFlow, innerID, innerOD, T, P, innerGas) for T, P in zip(sol.y[0], sol.y[2])])
  axes[0][1].plot(sol.x, y, color = 'tab:red')
  y = numpy.array([reynoldsNumber(outerMassFlow, outerID, outerOD, T, P, outerGas) for T, P in zip(sol.y[1], sol.y[3])])
  axes[0][1].plot(sol.x, y, color = 'tab:blue')
  if len(sol.y) > 4:
    y = numpy.array([reynoldsNumber(N2massFlow, 0., N2OD, T, P, 'N2') for T, P in zip(sol.y[4], sol.y[5])])
    axes[0][1].plot(sol.x, y, color = 'tab:green')
  y = numpy.array([nusseltNumber(outerMassFlow, outerID, outerOD, T_outer, Twall(innerMassFlow, innerID, innerOD, T_inner, P_inner, innerGas, outerMassFlow, outerID, outerOD, T_outer, P_outer, outerGas, roughness), P_outer, roughness, outerGas) \
                   for T_inner, P_inner, T_outer, P_outer in zip(sol.y[0], sol.y[2], sol.y[1], sol.y[3])])
  axes[0][1].plot(sol.x, y, color = 'tab:cyan')
  axes[0][1].set_yscale('log')
  axes[0][1].set_ylabel('Reynolds number')
  axes[0][1].grid(True, 'both')
  
  #y = numpy.array([hepak.HeCalc('D', 0, 'P', P_He, 'T', T, 1) for T, P_He in zip(sol.y[0], sol.y[2])])
  #axes[1][0].plot(sol.x, y, color = 'tab:red')
  y = numpy.array([he3pak.He3Density(P_He3, T) for T, P_He3 in zip(sol.y[1], sol.y[3])])
  axes[1][0].plot(sol.x, y, color = 'tab:blue')
  axes[1][0].set_ylabel('Density (kg/m3)')
  axes[1][0].grid(True, 'both')
  
  y = numpy.array([velocity(innerMassFlow, innerID, innerOD, T, P, innerGas) for T, P in zip(sol.y[0], sol.y[2])])
  axes[1][1].plot(sol.x, y, color = 'tab:red')
  y = numpy.array([velocity(outerMassFlow, outerID, outerOD, T, P, innerGas) for T, P in zip(sol.y[1], sol.y[3])])
  axes[1][1].plot(sol.x, y, color = 'tab:blue')
  if len(sol.y) > 4:
    y = numpy.array([velocity(N2massFlow, N2ID, N2OD, T, P, 'N2') for T, P in zip(sol.y[4], sol.y[5])])
    axes[1][1].plot(sol.x, y, color = 'tab:green')
  axes[1][1].set_yscale('log')
  axes[1][1].set_ylabel('Velocity (m/s)')
  axes[1][1].grid(True, 'both')
  
  axes[2][0].plot(sol.x, sol.y[2], color = 'tab:red')
  axes[2][0].plot(sol.x, sol.y[3], color = 'tab:blue')
  if len(sol.y) > 4:
    axes[2][0].plot(sol.x, sol.y[5], color = 'tab:green')
  axes[2][0].set_ylabel('Pressure (Pa)')
  axes[2][0].grid(True, 'both')  
  
  while True:
    try:
      plt.savefig('HEX7_{0:.4g}K_{1:.4g}K_{2:.4g}gs_{3:.4g}gs.pdf'.format(sol.y[0][0], sol.y[1][-1], innerMassFlow, outerMassFlow))
    except PermissionError:
      input('Could not save plots. If the target file is still opened in another application please close it. Press ENTER to try again.')
    else:
      break
  plt.close()


# calculate tube-in-tube counter-flow heat exchanger with cold He4 flowing through the inner tube and warm He3 flowing through the outer tube
def HEX7(HeMassFlow, He3MassFlow, T_He_in, T_He3_in):
  innerID = 0.
  innerOD = 7.747e-3
  outerID = 9.525e-3
  outerOD = 16.561e-3
  P_He3_in = 51120
  P_He_in = 110300
  meshSize = 10
  HEXlength = 19.
  roughness = 0.002e-3
  
  # differential-equation set
  # x = [position of grid points along heat exchanger]
  # y = [[temperatures of gas in inner tube at grid points], [temperatures of outer gas], [pressures of inner gas], [pressures of outer gas]]
  def dydx(x, y):
    dydx = [[],[],[],[]]
    for T_He, T_He3, P_He, P_He3 in zip(y[0], y[1], y[2], y[3]):
      T_He = max(T_He, T_He_in)
      T_He3 = max(T_He3, T_He_in)
      P_He = max(P_He, 1000.)
      P_He3 = max(P_He3, 1000.)
      T_wall = Twall(HeMassFlow, innerID, innerOD, T_He, P_He, 'He', He3MassFlow, outerID, outerOD, T_He3, P_He3, 'He3', roughness)
      dydx[0].append(dTdx(HeMassFlow, innerID, innerOD, innerOD*math.pi, T_He, T_wall, P_He, roughness, 'He'))
      dydx[1].append(-dTdx(He3MassFlow, outerID, outerOD, outerID*math.pi, T_He3, T_wall, P_He3, roughness, 'He3'))
      dydx[2].append(-dPdx(HeMassFlow, innerID, innerOD, roughness, T_He, P_He, 'He'))
      dydx[3].append(dPdx(He3MassFlow, outerID, outerOD, roughness, T_He3, P_He3, 'He3'))
    print('4He out: {0:.4g} K, 3He out: {1:.4g} K, 4He dP: {2:.4g} Pa, 3He dP: {3:.4g} Pa'.format(y[0][-1], y[1][0], y[2][-1] - y[2][0], y[3][0] - y[3][-1]))
    return dydx
  
  # boundary conditions: fixed temperatures and pressures at gas inlets
  def bc(ya, yb):
    return [ya[0] - T_He_in, yb[1] - T_He3_in, ya[2] - P_He_in, yb[3] - P_He3_in]
    
  sol = scipy.integrate.solve_bvp(dydx, bc, numpy.linspace(0., HEXlength, meshSize), \
                                  [numpy.linspace(T_He_in, T_He3_in, meshSize), \
  								   numpy.linspace(T_He_in, T_He3_in, meshSize), \
  								   numpy.linspace(P_He_in, P_He_in - 5000., meshSize), \
  								   numpy.linspace(P_He3_in - 1000., P_He3_in, meshSize)],
								  max_nodes = 100)
  print(sol.message)
  plotTubeInTubeHEX(sol, HeMassFlow, innerID, innerOD, 'He', He3MassFlow, outerID, outerOD, 'He3', 0., 0., roughness)


# calculate tube-in-tube counter-flow heat exchanger with cold He4 flowing through the inner tube and warm He4 flowing through the outer tube
def HEX7test(innerMassFlow, outerMassFlow, T_inner_in, T_outer_in):
  innerID = 0.
  innerOD = 7.747e-3
  outerID = 9.525e-3
  outerOD = 16.561e-3
  P_inner_in = 110e3
  P_outer_out = 110e3
  meshSize = 10
  HEXlength = 18.
  roughness = 0.004e-3

  # differential-equation set
  # x = [position of grid points along heat exchanger]
  # y = [[temperatures of gas in inner tube at grid points], [temperatures of outer gas], [pressures of inner gas], [pressures of outer gas]]
  def dydx(x, y):
    dydx = [[],[],[],[]]
    for T_inner, T_outer, P_inner, P_outer in zip(y[0], y[1], y[2], y[3]):
      T_inner = max(T_inner, T_inner_in)
      T_outer = max(T_outer, T_inner_in)
      P_inner = max(P_inner, 1000.)
      P_outer = max(P_outer, 1000.)
      T_wall = Twall(innerMassFlow, innerID, innerOD, T_inner, P_inner, 'He', outerMassFlow, outerID, outerOD, T_outer, P_outer, 'He', roughness)
      dydx[0].append(dTdx(innerMassFlow, innerID, innerOD, innerOD*math.pi, T_inner, T_wall, P_inner, roughness, 'He'))
      dydx[1].append(-dTdx(outerMassFlow, outerID, outerOD, outerID*math.pi, T_outer, T_wall, P_outer, roughness, 'He'))
      dydx[2].append(-dPdx(innerMassFlow, innerID, innerOD, roughness, T_inner, P_inner, 'He'))
      dydx[3].append(dPdx(outerMassFlow, outerID, outerOD, roughness, T_outer, P_outer, 'He'))
    print('Inner out: {0:.4g} K, Outer out: {1:.4g} K, 4He dP: {2:.4g} Pa, 3He dP: {3:.4g} Pa'.format(y[0][-1], y[1][0], y[2][-1] - y[2][0], y[3][0] - y[3][-1]))
    return dydx
  
  # boundary conditions: fixed temperatures and pressures at gas inlets
  def bc(ya, yb):
    return [ya[0] - T_inner_in, yb[1] - T_outer_in, ya[2] - P_inner_in, ya[3] - P_outer_out]
    
  sol = scipy.integrate.solve_bvp(dydx, bc, numpy.linspace(0., HEXlength, meshSize), \
                                  [numpy.linspace(T_inner_in, T_outer_in, meshSize), \
  								   numpy.linspace(T_inner_in, T_outer_in, meshSize), \
  								   numpy.linspace(P_inner_in, P_inner_in - 5000., meshSize), \
  								   numpy.linspace(P_outer_out, P_outer_out - 1000., meshSize)], \
								  max_nodes = 100)
  print(sol.message)
  plotTubeInTubeHEX(sol, innerMassFlow, innerID, innerOD, 'He', outerMassFlow, outerID, outerOD, 'He', 0., 0., roughness)


# calculate tube-in-tube-in-tube heat exchanger with cold He3 flowing out through the inner tube, warm He3 flowing in through second tube, and cold N2 flowing out through outer tube
def purifier(He3MassFlow, T_warm_in):
  innerID = 0.
  innerOD = 13.843e-3
  outerID = 15.875e-3
  outerOD = 26.035e-3
  N2ID = 28.575e-3
  N2OD = 32.131e-3
  P_in = 60e3
  P_N2_out = 101e3
  meshSize = 10
  HEXlength = 6.
  roughness = 0.002e-3
  
  # differential-equation set
  # x = [position of grid points along heat exchanger]
  # y = [[temperatures of outflowing He3], [temperatures of inflowing He3], [pressures of outflowing He3], [pressures of inflowing He3], [temperatures of ouflowing N2], [pressures of outflowing N2]]
  def dydx(x, y):
    dydx = [[],[],[],[],[],[]]
    N2massFlow = LN2evaporationFlow(y[5][0], He3MassFlow, y[1][0], y[3][0], 'He3')
    for T_inner, T_outer, P_inner, P_outer, T_N2, P_N2 in zip(y[0], y[1], y[2], y[3], y[4], y[5]):
      T_inner = max(T_inner, 77.)
      T_outer = max(T_outer, 77.)
      P_inner = max(P_inner, 1000.)
      P_outer = max(P_outer, 1000.)
      T_N2 = max(T_N2, 77.)
      P_N2 = max(P_N2, 1000.)
      T_wall = Twall(He3MassFlow, innerID, innerOD, T_inner, P_inner, 'He3', He3MassFlow, outerID, outerOD, T_outer, P_outer, 'He3', roughness)
      T_wall_N2 = Twall(He3MassFlow, outerID, outerOD, T_outer, P_outer, 'He3', N2massFlow, N2ID, N2OD, T_N2, P_N2, 'N2', roughness)
      dydx[0].append(dTdx(He3MassFlow, innerID, innerOD, innerOD*math.pi, T_inner, T_wall, P_inner, roughness, 'He3'))
      dydx[1].append(-dTdx2(He3MassFlow, outerID, outerOD, T_outer, T_wall, T_wall_N2, P_outer, roughness, 'He3'))
      dydx[2].append(-dPdx(He3MassFlow, innerID, innerOD, roughness, T_inner, P_inner, 'He3'))
      dydx[3].append(dPdx(He3MassFlow, outerID, outerOD, roughness, T_outer, P_outer, 'He3'))
      dydx[4].append(dTdx(N2massFlow, N2ID, N2OD, N2ID*math.pi, T_N2, T_wall_N2, P_N2, roughness, 'N2'))
      dydx[5].append(-dPdx(N2massFlow, N2ID, N2OD, roughness, T_N2, P_N2, 'N2'))
    print('Inner out: {0:.4g} K, Outer out: {1:.4g} K, Inner dP: {2:.4g} Pa, Outer dP: {3:.4g} Pa, N2 flow: {4:.4g} kg/s, N2 out: {5:.4g} K, N2 dP: {6:.4g} Pa'.format(y[0][-1], y[1][0], y[2][-1] - y[2][0], y[3][-1] - y[3][0], N2massFlow, y[4][-1], y[5][-1] - y[5][0]))
    return dydx
  
  # boundary conditions: fixed temperatures at gas inlets (cold-inlet temperature is given by N2 vapor pressure),
  # fixed pressure for warm-He3 inlet, cold-He3 inlet pressure = warm-He3 outlet pressure, cold-N2 inlet pressure = vapor pressure
  def bc(ya, yb):
    T_N2gas = CoolProp.CoolProp.PropsSI('T', 'P', ya[5], 'Q', 1, 'Nitrogen') + 2.
    return numpy.array([ya[0] - T_N2gas, yb[1] - T_warm_in, ya[2] - ya[3], yb[3] - P_in, ya[4] - T_N2gas, yb[5] - P_N2_out])
    
  sol = scipy.integrate.solve_bvp(dydx, bc, numpy.linspace(0., HEXlength, meshSize), \
                                  [numpy.linspace(80., T_warm_in, meshSize), \
  								   numpy.linspace(80., T_warm_in, meshSize), \
  								   numpy.linspace(P_in - 5000., P_in - 10000., meshSize), \
  								   numpy.linspace(P_in - 5000., P_in, meshSize), \
								   numpy.linspace(80., T_warm_in, meshSize), \
								   numpy.linspace(P_N2_out + 5000, P_N2_out, meshSize)], \
								  max_nodes = 1000)
  print(sol.message)
  N2massFlow = LN2evaporationFlow(sol.y[5][0], He3MassFlow, sol.y[1][0], sol.y[3][0], 'He3')
  print('Inner out: {0:.4g} K, Outer out: {1:.4g} K, Inner dP: {2:.4g} Pa, Outer dP: {3:.4g} Pa, N2 flow: {4:.4g} kg/s, N2 out: {5:.4g} K, N2 dP: {6:.4g} Pa'.format(sol.y[0][-1], sol.y[1][0], sol.y[2][-1] - sol.y[2][0], sol.y[3][-1] - sol.y[3][0], N2massFlow, sol.y[4][-1], sol.y[5][-1] - sol.y[5][0]))
  plotTubeInTubeHEX(sol, He3MassFlow, innerID, innerOD, 'He3', He3MassFlow, outerID, outerOD, 'He3', N2ID, N2OD, roughness)
  

# calculate tube-in-tube counter-flow heat exchanger with cold He3 flowing through the inner tube and warm He3 flowing through the outer tube
def purifierWithoutN2(He3MassFlow, T_warm_in):
  innerID = 0.
  innerOD = 10.922e-3
  outerID = 12.7e-3
  outerOD = 19.939e-3
  P_in = 80e3
  meshSize = 10
  HEXlength = 5.
  roughness = 0.002e-3
  
  # differential-equation set
  # x = [position of grid points along heat exchanger]
  # y = [[temperatures of gas in inner tube at grid points], [temperatures of outer gas], [pressures of inner gas], [pressures of outer gas]]
  def dydx(x, y):
    dydx = [[],[],[],[]]
    for T_inner, T_outer, P_inner, P_outer in zip(y[0], y[1], y[2], y[3]):
      T_inner = max(T_inner, 77.)
      T_outer = max(T_outer, 77.)
      P_inner = max(P_inner, 1000.)
      P_outer = max(P_outer, 1000.)
      T_wall = Twall(He3MassFlow, innerID, innerOD, T_inner, P_inner, 'He3', He3MassFlow, outerID, outerOD, T_outer, P_outer, 'He3', roughness)
      dydx[0].append(dTdx(He3MassFlow, innerID, innerOD, innerOD*math.pi, T_inner, T_wall, P_inner, roughness, 'He3'))
      dydx[1].append(-dTdx(He3MassFlow, outerID, outerOD, outerID*math.pi, T_outer, T_wall, P_outer, roughness, 'He3'))
      dydx[2].append(-dPdx(He3MassFlow, innerID, innerOD, roughness, T_inner, P_inner, 'He3'))
      dydx[3].append(dPdx(He3MassFlow, outerID, outerOD, roughness, T_outer, P_outer, 'He3'))
    N2massFlow = LN2evaporationFlow(110e3, He3MassFlow, y[1][0], y[3][0], 'He3')
    print('Inner out: {0:.4g} K, Outer out: {1:.4g} K, Inner dP: {2:.4g} Pa, Outer dP: {3:.4g} Pa, N2 flow: {4:.4g} kg/s'.format(y[0][-1], y[1][0], y[2][-1] - y[2][0], y[3][-1] - y[3][0], N2massFlow))
    return dydx
  
  # boundary conditions: fixed temperatures at gas inlets,
  # fixed pressure for warm-He3 inlet, cold-He3 inlet pressure = warm-He3 outlet pressure
  def bc(ya, yb):
    return numpy.array([ya[0] - 80., yb[1] - T_warm_in, ya[2] - ya[3], yb[3] - P_in])
    
  sol = scipy.integrate.solve_bvp(dydx, bc, numpy.linspace(0., HEXlength, meshSize), \
                                  [numpy.linspace(80., T_warm_in, meshSize), \
  								   numpy.linspace(80., T_warm_in, meshSize), \
  								   numpy.linspace(P_in - 5000., P_in - 10000., meshSize), \
  								   numpy.linspace(P_in - 5000., P_in, meshSize)], \
								  max_nodes = 1000)
  print(sol.message)
  plotTubeInTubeHEX(sol, He3MassFlow, innerID, innerOD, 'He3', He3MassFlow, outerID, outerOD, 'He3', 0., 0., roughness)
  


# nominal flows and temperatures of HEX7
He3MassFlow = 0.275e-3
HeMassFlow = 0.478e-3
T_He3_in = 300.
T_He_in = 5.
HEX7(HeMassFlow, He3MassFlow, T_He_in, T_He3_in)
#HEX7(HeMassFlow/2., He3MassFlow/2., T_He_in, T_He3_in)
#HEX7(HeMassFlow/4., He3MassFlow/4., T_He_in, T_He3_in)
#HEX7(HeMassFlow/13., He3MassFlow/10., T_He_in, T_He3_in)

# flows and temperatures for test at KEK
#HEX7test(0.4e-3, 0.35e-3, 5., 295.)

# nominal flow for He3 purifier
#purifier(1.e-3, 300.)
#purifierWithoutN2(1.e-3, 300.)