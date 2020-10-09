import hepak
import scipy.integrate
import scipy.optimize
import matplotlib.pyplot as plt
import numpy
import math

length = 50.
verticalDrop = 4.
heatLoad = 0.05 # W/m
pressureDrop = 3000./50. # Pa/m
massFlow = 3.5e-3 # kg/s
crossSection = 0.005**2*math.pi # m^2
reservoirCrossSection = 0.2**2*math.pi #m^2

inletPressure = 1.225e5
valvePressureDrop = 60

def dydx(x, y):
  H = y[0]
  P = max(100., y[1])
  density = hepak.HeCalc('D', 0, 'P', P, 'H', H, 1)
  dHdx = heatLoad/massFlow + pressureDrop/density
  dPdx = -pressureDrop
  if x < verticalDrop:
    dHdx = dHdx + 9.81
    dPdx = dPdx + density*9.81
  return [dHdx, dPdx]
  
result = scipy.integrate.solve_ivp(dydx, (0, length), [hepak.HeCalc('H', 0, 'P', inletPressure, 'SL', 0, 1), inletPressure - valvePressureDrop], max_step = 1.)

exitEnthalpy = result.y[0][-1]
exitPressure = result.y[1][-1]
exitDensity = hepak.HeCalc('D', 0, 'P', exitPressure, 'H', exitEnthalpy, 1)
exitLiquidDensity = hepak.HeCalc('D', 0, 'P', exitPressure, 'SL', 0, 1)
exitVaporDensity = hepak.HeCalc('D', 0, 'P', exitPressure, 'SV', 0, 1)
exitVaporMassFraction = hepak.HeCalc('X', 0, 'H', exitEnthalpy, 'P', exitPressure, 1)
exitVaporVolumeFraction = exitVaporMassFraction/exitVaporDensity*exitDensity

fullLiquidVelocity = massFlow/exitLiquidDensity/crossSection

exitLiquidVelocity = fullLiquidVelocity*(1 - exitVaporMassFraction)/(1 - exitVaporVolumeFraction)
exitVaporVelocity = fullLiquidVelocity*exitVaporMassFraction/exitVaporVolumeFraction*exitLiquidDensity/exitVaporDensity

print(exitPressure)


plotRows = 2
plotCols = 2
fig, axes = plt.subplots(plotRows, plotCols, figsize=(plotCols*5,plotRows*5), sharex='all')
fig.set_tight_layout(True)
axes[0][0].plot(result.t, result.y[1], color = 'tab:blue', label = 'Pressure')
axes[0][0].set_ylabel('Pressure (Pa)')
axes[0][1].plot(result.t, result.y[0], color = 'tab:blue', label = 'Enthalpy')
axes[0][1].set_ylabel('Enthalpy (J/kg)')
y = [max(0., hepak.HeCalc('X', 0, 'H', H, 'P', P, 1)) for H, P in zip(result.y[0], result.y[1])]
axes[1][0].plot(result.t, y, color = 'tab:blue', label = 'Vapor quality')
axes[1][0].set_ylabel('Vapor quality')
axes[1][0].set_xlabel('Length (m)')
y = [massFlow/hepak.HeCalc('D', 0, 'H', H, 'P', P, 1)/crossSection for H, P in zip(result.y[0], result.y[1])]
axes[1][1].plot(result.t, y, color = 'tab:blue', label = 'Velocity')
axes[1][1].set_ylabel('Velocity (m/s)')
axes[1][1].set_xlabel('Length (m)')



while True:
  try:
    plt.savefig('HeTransfer.pdf')
  except PermissionError:
    input('Could not save plots. If the target file is still opened in another application please close it. Press ENTER to try again.')
  else:
    break
plt.close()