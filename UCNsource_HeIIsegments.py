import UCNsource
import hepak
import math
import scipy.integrate
import numpy
import matplotlib.pyplot as plt
import csv

hbar = 1.054571817e-34 # J s
m_n = 1.67492749804e-27 # kg
kB = 1.380649e-23 # J/K
massHe = 4.002602 *  1.66053906660e-27 # kg
thermalHeCrossSection = 0.76e-28 # m2
HeScatteringLength = 3.26e-15 # m
JpereV = 1.602176634e-19 # J/eV

def liquidLifetime(T, YoshikiParameter):
  return 1./(YoshikiParameter*T**7)

def vaporLifetime(temperature, pressure):
  return math.sqrt(math.pi/8. * kB * massHe * temperature) / pressure / thermalHeCrossSection

def realFermiPotential(density):
  return 2.*math.pi*hbar**2/m_n*density/massHe*HeScatteringLength/JpereV*1e9
  
def imaginaryFermiPotential(lifetime):
  return 0.5*hbar/lifetime/JpereV*1e9

def segmentedVaporTemperature(T_low, T_high, lengths):
  dTdx = (T_high - T_low)/sum(lengths)
  pressure = hepak.HeCalc('P', 0, 'T', T_low, 'SV', 0, 1)
  x = 0.
  fermiPotentials = []
  for length in lengths:
    x2 = x + length
    meanDensity = scipy.integrate.quad(lambda y: hepak.HeCalc('D', 0, 'T', T_low + y*dTdx, 'P', pressure, 1), x, x2)[0]/length
    meanTemperature = T_low + (dTdx*x2 + dTdx*x)/2.
    meanLifetime = scipy.integrate.quad(lambda y: vaporLifetime(T_low + dTdx*y, pressure), x, x2)[0]/length
    fermiImag = imaginaryFermiPotential(meanLifetime)
    fermiPotentials.append( (realFermiPotential(meanDensity), fermiImag) )
#    print('{0:.2f}m: {1:.3f}K, {2:.1f}s'.format(x + length/2., meanTemperature, meanLifetime))
    x = x2
  
  return fermiPotentials
  

def simpleHeatFluxModel(T_low, pressure, heatLoad, channelDiameter, channelLength, conductivityModel):
  heatFluxModel = lambda x: -heatLoad/channelDiameter**2*4./math.pi

  fig, ax = plt.subplots()
  xdata = numpy.linspace(0., channelLength, 100)
  ax.plot(xdata, [heatFluxModel(x) for x in xdata])
  ax.set_xlabel('Distance from HEX1 center (m)')
  ax.set_ylabel('Heat flux in He-II (W/m2)')
  ax.set_title(r'$T_\mathrm{{low}}$ = {0:.3f}K, heat load = {1:.1f}W, {2}'.format(T_low, heatLoad, conductivityModel))
  ax.grid()
  filename = '{0:.3f}K_{1:.2f}W_{2}_Q'.format(T_low, heatLoad, conductivityModel)
  plt.savefig(filename + '.pdf')

  return UCNsource.GorterMellink(T_low, pressure, heatFluxModel, 0., channelLength, conductivityModel)

def refinedHeatFluxModel(T_low, pressure, beamHeating, staticHeat, funneledHeat, bulbLength, bulbDiameter, channelLength, channelDiameter, HEX1length, funnelLength, funnelDiameter, conductivityModel):
  x0 = UCNsource.HeIIlowestTemperaturePosition(beamHeating, bulbLength, channelLength, HEX1length, staticHeat, funneledHeat, funnelLength)
  heatFluxModel = lambda x: UCNsource.HeIIheatFlux(x, beamHeating, bulbDiameter/2., bulbLength, channelDiameter**2/4.*math.pi, channelLength, HEX1length, staticHeat, funneledHeat, funnelDiameter**2/4.*math.pi, funnelLength)

  fig, ax = plt.subplots()
  xdata = numpy.linspace(-funnelLength - HEX1length/2., HEX1length/2. + channelLength + bulbLength, 100)
  ax.plot(xdata, [heatFluxModel(x) for x in xdata])
  ax.set_xlabel('Distance from HEX1 center (m)')
  ax.set_ylabel('Heat flux in He-II (W/m2)')
  ax.set_title(r'$T_\mathrm{{low}}$ = {0:.3f}K, heat load = {1:.2f} + {2:.2f} + {3:.2f}W'.format(T_low, beamHeating, staticHeat, funneledHeat))
  ax.grid()
  filename = '{0:.3f}K_{1:.2f}+{2:.2f}+{3:.2f}W_Q'.format(T_low, beamHeating, staticHeat, funneledHeat)
  plt.savefig(filename + '.pdf')

  profileUpstream = UCNsource.GorterMellink(T_low, pressure, heatFluxModel, x0, HEX1length/2. + channelLength + bulbLength, conductivityModel)
  profileDownstream = UCNsource.GorterMellink(T_low, pressure, heatFluxModel, x0, -HEX1length/2. - funnelLength, conductivityModel)
  return profileUpstream, profileDownstream

def segmentedLiquidTemperature(temperatureProfile, pressure, segmentation, heatLoad, conductivityModel, YoshikiParameter):
  with open('{0:.3f}K_{1}_{2}_{3:.3f}.csv'.format(T_low, heatLoad, conductivityModel, YoshikiParameter), 'w', newline = '') as csvfile:
    csvWriter = csv.writer(csvfile)
    csvWriter.writerow(['Position (m)', 'Local temperature (K)', 'Mean temperature (K)', 'Mean UCN lifetime (s)', 'Mean real Fermi potential (neV)', 'Mean imaginary Fermi potential (neV)'])
    fermiPotentials = []
    for x1, x2 in zip(segmentation[:-1], segmentation[1:]):
      length = x2 - x1
      meanDensity = scipy.integrate.quad(lambda x: hepak.HeCalc('D', 0, 'T', temperatureProfile(x)[0], 'P', pressure, 1), x1, x2)[0]/length
      meanTemperature = scipy.integrate.quad(lambda x: temperatureProfile(x)[0], x1, x2)[0]/length
      meanLifetime = scipy.integrate.quad(lambda x: liquidLifetime(temperatureProfile(x)[0], YoshikiParameter), x1, x2)[0]/length
      fermiPotentials.append( (realFermiPotential(meanDensity), imaginaryFermiPotential(meanLifetime)) )
      xmid = (x1 + x2)/2.
#      print('{0:.2f}m: {1:.3f}K, {2:.1f}s'.format(xmid,meanTemperature, meanLifetime))
      csvWriter.writerow([xmid, temperatureProfile(xmid)[0], meanTemperature, meanLifetime, fermiPotentials[-1][0], fermiPotentials[-1][1]])
      csvWriter.writerow([x2, temperatureProfile(x2)[0], '', '', '', ''])
      x = x2
      
  fig, ax = plt.subplots()
  xdata = numpy.linspace(segmentation[0], segmentation[-1], 100, endpoint = True)
  ax.plot(xdata, [temperatureProfile(x)[0] for x in xdata])
  ax.set_xlabel('Distance from HEX1 center (m)')
  ax.set_ylabel('Temperature (K)')
  ax.set_title(r'$T_\mathrm{{low}}$ = {0:.3f}K, heat load = {1}, {2}'.format(T_low, heatLoad, conductivityModel))
  ax.grid()
  plt.savefig('{0:.3f}K_{1}_{2}.pdf'.format(T_low, heatLoad, conductivityModel))
  


  return fermiPotentials


#liquidSegmentation = numpy.concatenate([numpy.linspace(-0.55, -0.3, 6, True), numpy.linspace(-0.2077, 0.3, 5, False), numpy.linspace(0.3, 2.65, 20, False), numpy.linspace(2.65, 3.05, 6, True)])
liquidSegmentation = numpy.concatenate([numpy.linspace(-0.3985 - 0.283, -0.3985, 10, False), numpy.linspace(-0.3985, 3.0515, 69, True)])
roomTemperature = 300.
#vaporLengths = [0.251, 0.058, 0.1384, 0.1384, 0.1384, 0.1384]
vaporLengths = [0.0283]*10 + [0.05]*22

parameters = {
'Beam heating':                   8.1, # Beam on
#'Beam heating':                   0., # Beam off
'He-II static heat':              0.25,
'He-II funnel heat':              1.25,
'He reservoir static load':       0.6,
'1K pot static load':             0.05,
#'3He pumping speed':              4300., # 1.14g/s @ 5.85 torr (Busch proposal)
#'He pumping speed':               2000., # 1.2g/s @ 9.9 torr (Busch proposal)
#'Pump inlet temperature':         290.,
'3He pressure drop':              100.,
#'3He pressure drop':              20000., # Standby mode
'He pressure drop':               100.,
#'He pressure drop':               10000., # Standby mode
'He reservoir inlet temperature': 10.,
'HEX4 exit temperature':          2.8,
'HEX5 exit temperature':          2.8,
#'1K pot inlet temperature':       2.8,
'3He inlet pressure':             50000.,
'Channel diameter':               0.148,
'Channel length':                 2.356,
'HEX1 length':                    0.6,
'HEX1 surface':                   1.67,
'HeII overfill':                  0.05,
'He reservoir pressure':          1.2e5,
'Isopure He flow':                0.,
#'Isopure He flow':                0.00014/3, # isopure condensation
'Isopure He pressure':            100e2,
'20K shield temperature':         20.,
'100K shield temperature':        100.,
}

result = UCNsource.calcUCNSource(parameters)

print('3He flow: {0:.4g} g/s'.format(result['3He flow']*1000))
print('He reservoir evaporation: {0:.4g} g/s'.format(result['He reservoir flow']*1000))
print('1K pot: {0:.4} K, {1:.4g} g/s'.format(result['1K pot temperature'], result['1K pot flow']*1000))
print('He consumption: {0:.3g} L/h'.format((max(result['He reservoir flow'], result['20K shield flow'], result['100K shield flow']) + result['1K pot flow'])/hepak.HeCalc('D', 0, 'P', 1013e2, 'SL', 0., 1)*1000*3600))
print('Required shield flows: {0:.3g} g/s, {1:.3g} g/s'.format(result['20K shield flow']*1000, result['100K shield flow']*1000))
print('He-II pressure: {0:.4g} + {1:.4g} Pa'.format(result['HeII vapor pressure'], result['HeII pressure head']))
print('T: {0:.4g} -> {1:.4g} -> {2:.4g} -> {3:.4g} -> {4[0]:.4g} to {4[1]:.4g} K'.format(result['T_3He'], result['T_HEX1_low'], result['T_HEX1_high'], result['T_HeII_low'], result['T_HeII_high']))

### Steve's parameters: T_low = 0.982 and 0.876, each for VanSciver and HEPAK, UCNproduction 1.62e7

#for T_low, heatLoad, conductivityModel, yoshikiParameter in zip([1.023, 0.9, 1.1, 0.8, 1.0, 1.023, 1.023, 1.023],\
#                                                            [9.6, 9.6, 9.6, 0.2, 0.2, 9.6, 9.6, 9.6],\
#                                                            ['VanSciver', 'VanSciver', 'VanSciver', 'VanSciver', 'VanSciver', 'HEPAK', 'VanSciver', 'VanSciver'],\
#                                                            [0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.008, 0.024]):
for T_low, beamHeating, conductivityModel, yoshikiParameter in zip([0.982, 0.982, 0.876, 0.876], [8.1, 8.1, 8.1, 8.1], ['VanSciver', 'HEPAK', 'VanSciver', 'HEPAK'], [0.016, 0.016, 0.016, 0.016]):
#for T_low, beamHeating, conductivityModel, yoshikiParameter in zip([result['T_HeII_low']], [parameters['Beam heating']], ['VanSciver'], [0.016]):
  vaporPressure = hepak.HeCalc('P', 0, 'T', T_low, 'SV', 0, 1)
  liquidPressure = vaporPressure + hepak.HeCalc('D', 0, 'T', T_low, 'SL', 0., 1)*9.81*parameters['HeII overfill']
  
#  simpleProfile = simpleHeatFluxModel(T_low, liquidPressure, beamHeating, parameters['Channel diameter'], parameters['Channel length'], conductivityModel)
#  liquidPotentials = segmentedLiquidTemperature(lambda x: simpleProfile(x) if x > 0. else simpleProfile(simpleProfile.t_min), liquidPressure, liquidSegmentation, '{0:.2f}W'.format(heatLoad), conductivityModel, yoshikiParameter)
  refinedProfileUpstream, refinedProfileDownstream = refinedHeatFluxModel(T_low, liquidPressure, beamHeating, parameters['He-II static heat'], parameters['He-II funnel heat'], 0.386, 0.36, parameters['Channel length'], parameters['Channel diameter'], parameters['HEX1 length'], 0.283, 0.0955, conductivityModel)
  liquidPotentials = segmentedLiquidTemperature(lambda x: refinedProfileUpstream(x) if x >= refinedProfileUpstream.t_min else refinedProfileDownstream(x), \
                                                liquidPressure, liquidSegmentation, '{0:.2f}+{1:.2f}+{2:.2f}W'.format(beamHeating, parameters['He-II static heat'], parameters['He-II funnel heat']), conductivityModel, yoshikiParameter)

  vaporPotentials = segmentedVaporTemperature(T_low, roomTemperature, vaporLengths)

  print('# T_low {0:.3f}K, {1}, heat load {2:.2f}+{3:.2f}+{4:.2f}W, B {5:.3f}/sK^7'.format(T_low, conductivityModel, beamHeating, parameters['He-II static heat'], parameters['He-II funnel heat'], yoshikiParameter))
  T_high = refinedProfileUpstream(refinedProfileUpstream.t_max)[0]
  print('# resulting T_high {0:.3f}K'.format(T_high))
  for i, W in enumerate(liquidPotentials):
    print('LHe{0}   {1[0]:.3g}     {1[1]:.3g}     0 0 0 0'.format(i, W))

  for i, W in enumerate(vaporPotentials):
    print('He{0}  {1[0]:.3g}     {1[1]:.3g}   0 0 0 0'.format(i, W))
  print('HeRT    {0:.3g}      {1:.3g}      0 0 0 0'.format(realFermiPotential(hepak.HeCalc('D', 0, 'P', vaporPressure, 'T', roomTemperature, 1)), imaginaryFermiPotential(vaporLifetime(roomTemperature, vaporPressure))))
  print('\n\n')