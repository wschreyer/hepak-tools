import UCNsource
import hepak
import math
import scipy.integrate

def printSegmentedTemperature(T_low, pressure, heatLoad, channelDiameter, channelLength, conductivityModel, segmentLength, YoshikiParameter):
  profile = UCNsource.GorterMellink(T_low, pressure, heatLoad, channelDiameter, channelLength, conductivityModel)

  A = math.pi/4.*channelDiameter**2
  bottleVolume = 0.034
  temperatureVolumeIntegral = scipy.integrate.quad(lambda x: profile(x)[0], 0., channelLength)[0]*A + profile(channelLength)[0]*bottleVolume
  meanTemperature = temperatureVolumeIntegral/(channelLength*A + bottleVolume)
  lifetimeVolumeIntegral = scipy.integrate.quad(lambda x: 1./(YoshikiParameter*profile(x)[0]**7), 0., channelLength)[0]*A + 1./(YoshikiParameter*profile(channelLength)[0]**7)*bottleVolume
  meanLifetime = lifetimeVolumeIntegral/(channelLength*A + bottleVolume)
  print('Mean: {0:.3f} K, {1:.1f} s'.format(meanTemperature, meanLifetime))

  x = 0.
  while x < channelLength:
    meanLifetime = scipy.integrate.quad(lambda y: 1./(YoshikiParameter*profile(y)[0]**7), x, x + segmentLength)[0]/segmentLength
    print('{0:.2f}m: {1:.3f}K, {2:.1f}s, {3:.3g}neV'.format(x + segmentLength/2., \
                                                scipy.integrate.quad(lambda y: profile(y)[0], x, x + segmentLength)[0]/segmentLength, \
                                                meanLifetime, 0.5*6.582e-16/meanLifetime*1e9))
    x = x + segmentLength


parameters = {
'Beam heating':                   8.1, # Beam on
#'Beam heating':                   0., # Beam off
'Static heat':                    1.,
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
'Channel length':                 2.5,
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
#print('He-II pressure: {0:.4g} + {1:.4g} Pa'.format(result['HeII vapor pressure'], result['HeII pressure head']))
print('T: {0:.4g} -> {1:.4g} -> {2:.4g} -> {3:.4g} -> {4[0]:.4g} to {4[1]:.4g} K'.format(result['T_3He'], result['T_HEX1_low'], result['T_HEX1_high'], result['T_HeII_low'], result['T_HeII_high']))
printSegmentedTemperature(result['T_HeII_low'], result['HeII vapor pressure'] + result['HeII pressure head'], parameters['Beam heating'], parameters['Channel diameter'], parameters['Channel length'], 'HEPAK', 0.1, 0.016)


