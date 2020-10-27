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
    print('{0:.2f}m: {1:.3f}K, {2:.1f}s'.format(x + length/2., meanTemperature, meanLifetime))
    x = x2
  
  return fermiPotentials
  

def printSegmentedTemperature(T_low, pressure, heatLoad, channelDiameter, channelLength, conductivityModel, segmentLength, YoshikiParameter):
  profile = UCNsource.GorterMellink(T_low, pressure, heatLoad, channelDiameter, channelLength, conductivityModel)
  
  A = math.pi/4.*channelDiameter**2
  bottleVolume = 0.034
  temperatureVolumeIntegral = scipy.integrate.quad(lambda x: profile(x)[0], 0., channelLength)[0]*A + profile(channelLength)[0]*bottleVolume
  meanTemperature = temperatureVolumeIntegral/(channelLength*A + bottleVolume)
  lifetimeVolumeIntegral = scipy.integrate.quad(lambda x: liquidLifetime(profile(x)[0], YoshikiParameter), 0., channelLength)[0]*A + liquidLifetime(profile(channelLength)[0], YoshikiParameter)*bottleVolume
  meanLifetime = lifetimeVolumeIntegral/(channelLength*A + bottleVolume)
  print('Mean: {0:.3f} K, {1:.1f} s'.format(meanTemperature, meanLifetime))

  with open('{0:.3f}K_{1:.2f}W_{2}_{3:.3f}.csv'.format(T_low, heatLoad, conductivityModel, YoshikiParameter), 'w', newline = '') as csvfile:
    csvWriter = csv.writer(csvfile)
    csvWriter.writerow(['Position (m)', 'Mean temperature (K)'])
    x = 0.
    fermiPotentials = []
    while x < channelLength:
      meanDensity = scipy.integrate.quad(lambda y: hepak.HeCalc('D', 0, 'T', profile(y)[0], 'P', pressure, 1), x, x + segmentLength)[0]/segmentLength
      meanTemperature = scipy.integrate.quad(lambda y: profile(y)[0], x, x + segmentLength)[0]/segmentLength
      meanLifetime = scipy.integrate.quad(lambda y: liquidLifetime(profile(y)[0], YoshikiParameter), x, x + segmentLength)[0]/segmentLength
      fermiPotentials.append( (realFermiPotential(meanDensity), imaginaryFermiPotential(meanLifetime)) )
      print('{0:.2f}m: {1:.3f}K, {2:.1f}s'.format(x + segmentLength/2.,meanTemperature, meanLifetime))
      csvWriter.writerow([x + segmentLength/2., meanTemperature])
      x = x + segmentLength
  
  return profile, fermiPotentials


segmentLength = 0.1341
numberOfSegments = 17
roomTemperature = 300.
vaporLengths = [0.251, 0.058, 0.1384, 0.1384, 0.1384, 0.1384]

parameters = {
'Beam heating':                   8.1, # Beam on
#'Beam heating':                   0., # Beam off
'Isopure static heat':            1.,
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
'Channel length':                 numberOfSegments*segmentLength,
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

it = 0
for T_low, heatLoad, GorterMellink, yoshikiParameter in zip([1.023, 0.9, 1.1, 0.8, 1.0, 1.023, 1.023, 1.023],\
                                                            [9.6, 9.6, 9.6, 0.2, 0.2, 9.6, 9.6, 9.6],\
                                                            ['VanSciver', 'VanSciver', 'VanSciver', 'VanSciver', 'VanSciver', 'HEPAK', 'VanSciver', 'VanSciver'],\
                                                            [0.016, 0.016, 0.016, 0.016, 0.016, 0.016, 0.008, 0.024]):
  vaporPressure = hepak.HeCalc('P', 0, 'T', T_low, 'SV', 0, 1)
  liquidPressure = vaporPressure + hepak.HeCalc('D', 0, 'T', T_low, 'SL', 0., 1)*9.81*parameters['HeII overfill']
  profile, liquidPotentials = printSegmentedTemperature(T_low, liquidPressure, heatLoad, parameters['Channel diameter'], parameters['Channel length'], GorterMellink, segmentLength, yoshikiParameter)

  fig, ax = plt.subplots()
  xdata = numpy.linspace(0., parameters['Channel length'], int(parameters['Channel length']/segmentLength), endpoint = True)
  ax.plot(xdata, [profile(x)[0] for x in xdata])
  ax.set_xlabel('Distance from HEX1 (m)')
  ax.set_ylabel('Temperature (K)')
  ax.set_title(r'$T_\mathrm{{low}}$ = {0:.3f}K, heat load = {1:.1f}W, {2}, B = {3:.3f}/(s K$^7$)'.format(T_low, heatLoad, GorterMellink, yoshikiParameter))
  ax.grid()
  plt.savefig('temperatureProfile_{0}.pdf'.format(it))
  it = it + 1

  vaporPotentials = segmentedVaporTemperature(T_low, roomTemperature, vaporLengths)

  print('# T_low {0:.3f}K, {1}, heat load {2:.1f}W, B {3:.3f}/sK^7'.format(T_low, GorterMellink, heatLoad, yoshikiParameter))
  T_high = profile(parameters['Channel length'])[0]
  print('# resulting T_high {0:.3f}K'.format(T_high))
  print('LHeHEX1   {0:.3g}     {1:.3g}     0 0 0 0'.format(realFermiPotential(hepak.HeCalc('D', 0, 'T', T_low, 'SL', 0, 1)), imaginaryFermiPotential(liquidLifetime(T_low, yoshikiParameter))))
  for i, W in enumerate(liquidPotentials):
    print('LHe{0}   {1[0]:.3g}     {1[1]:.3g}     0 0 0 0'.format(i, W))
  print('LHeBottle   {0:.3g}     {1:.3g}     0 0 0 0'.format(realFermiPotential(hepak.HeCalc('D', 0, 'T', T_high, 'P', liquidPressure, 1)), imaginaryFermiPotential(liquidLifetime(T_high, yoshikiParameter))))

  print('HeSurface    {0:.3g}      {1:.3g}      0 0 0 0'.format(realFermiPotential(hepak.HeCalc('D', 0, 'T', T_low, 'SV', 0, 1)), imaginaryFermiPotential(vaporLifetime(T_low, vaporPressure))))
  for i, W in enumerate(vaporPotentials):
    print('He{0}  {1[0]:.3g}     {1[1]:.3g}   0 0 0 0'.format(i, W))
  print('HeRT    {0:.3g}      {1:.3g}      0 0 0 0'.format(realFermiPotential(hepak.HeCalc('D', 0, 'P', vaporPressure, 'T', roomTemperature, 1)), imaginaryFermiPotential(vaporLifetime(roomTemperature, vaporPressure))))
  print('\n\n')