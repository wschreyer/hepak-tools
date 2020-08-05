import hepak
import he3pak
import math
import numpy
import scipy.integrate
import scipy.optimize
import scipy.interpolate
import matplotlib.pyplot as plt
import pumpdata
import HEXdata

# cooling power required to cool He flow from inlet temperature to outlet temperature
def HeCoolingLoad(flow, pressure, inletTemperature, outletTemperature):
  outletTemperature = max(hepak.HeConst(3), outletTemperature)
  return flow*(hepak.HeCalc('H', 0, 'P', pressure, 'T', inletTemperature, 1) - hepak.HeCalc('H', 0, 'P', pressure, 'T', outletTemperature, 1))

# flow required to remove heat load
def HeCoolingFlow(heatLoad, pressure, inletTemperature, outletTemperature):
  return heatLoad/(hepak.HeCalc('H', 0, 'P', pressure, 'T', outletTemperature, 1) - hepak.HeCalc('H', 0, 'P', pressure, 'T', inletTemperature, 1))

# cooling power required to cool He3 flow from inlet temperature to outlet temperature
def He3CoolingLoad(flow, pressure, inletTemperature, outletTemperature):
  inletDensity = he3pak.He3Density(pressure, inletTemperature)
  inletEnthalpy = he3pak.He3Prop(6, inletDensity, inletTemperature)
  outletDensity = he3pak.He3Density(pressure, outletTemperature)
  outletEnthalpy = he3pak.He3Prop(6, outletDensity, outletTemperature)
  return flow*(inletEnthalpy - outletEnthalpy) # heat load on 1K pot = He3 gas flow * He3 enthalpy difference

# calculate evaporation rate of helium reservoir when flowing He3 through it (assuming He3 is cooled to reservoir temperature)
def HeReservoirEvaporation(He3Flow, He3Pressure, He3InletTemperature, IPFlow, IPPressure, IPInletTemperature, reservoirPressure):
  reservoirTemperature = hepak.HeCalc('T', 0, 'P', reservoirPressure, 'SV', 0., 1)
  He3HeatLoad = He3CoolingLoad(He3Flow, He3Pressure, He3InletTemperature, reservoirTemperature)
  IPHeatLoad = HeCoolingLoad(IPFlow, IPPressure, IPInletTemperature, reservoirTemperature)
  return (He3HeatLoad + IPHeatLoad)/hepak.HeCalc(7, 0, 'P', reservoirPressure, 'SL', 0., 1) # evaporation flow = heatLoad/latent heat
  
# calculate required flow to cool 20K and 100K shields, with additional heat load by isopure He flowing over shields
# static heat loads on shields are assumed fixed
def shieldFlow(IPFlow, IPPressure, inletPressure):
  IPHeatLoad_100 = HeCoolingLoad(IPFlow, IPPressure, 300., 100.)
  IPHeatLoad_20 = HeCoolingLoad(IPFlow, IPPressure, 100., 20.)
  flow_20 = HeCoolingFlow(4.3 + 0.7 + 1.1 + 2. + 6.6 + IPHeatLoad_20, inletPressure, 10., 20.)
  flow_100 = HeCoolingFlow(49. + 8.9 + 5.7 + 14. + IPHeatLoad_100, inletPressure, 30., 100.)
  return flow_20, flow_100

pumpData = pumpdata.loadPumpData('pumpdata/Busch_2stage_He3_torqueControl.csv')
	  
# calculate 1K pot temperature when pumping away a certain gas flow
def OneKPotTemperature(HeFlow, pumpingPressureDrop):
  HeFlow = max(0., min(HeFlow, 0.1))
#  inletPressure = scipy.optimize.root_scalar(lambda P: hepak.HeCalc(5, 0, 'P', P, 'T', pumpInletTemperature, 0)*gasFlow*hepak.HeConst(4)*inletTemperature/(pumpingSpeed/3600) - P, bracket = (1., 1e5)) # He pressure = Z(P, T)*dm/dt*R_s*T/(dV/dt)
#  if not inletPressure.converged:
#    print('Could not determine pump inlet pressure')
#  pumpInletPressure = HeFlow*hepak.HeConst(4)*pumpInletTemperature/(pumpingSpeed/3600) # P = dm/dt R_s T / (dV/dt)

  pumpInletPressure = scipy.optimize.root_scalar(lambda P: pumpData(P)/0.75 - HeFlow, bracket = (15., 12000.)) # find inlet pressure for given flow in pump data for He3 (divide by 0.75 to scale to He4)
  if not pumpInletPressure.converged:
    print("Could not calculate He4 pump inlet pressure!")
  T = hepak.HeCalc('T', 0, 'P', pumpInletPressure.root + pumpingPressureDrop, 'SL', 0., 1) # T_sat(P + dP)
  return T

# calculate evaporation rate from 1K pot at certain temperature and He3 flow
def OneKPotEvaporation(He3Flow, He3Pressure, He3InletTemperature, IPFlow, IPPressure, IPInletTemperature, HeInletPressure, HeInletTemperature, OneKPotTemperature):
  He3HeatLoad = He3CoolingLoad(He3Flow, He3Pressure, He3InletTemperature, OneKPotTemperature)
  IPHeatLoad = HeCoolingLoad(IPFlow, IPPressure, IPInletTemperature, OneKPotTemperature)

  HeInletEnthalpy = hepak.HeCalc('H', 0, 'T', HeInletTemperature, 'P', HeInletPressure, 1)
  HeVaporEnthalpy = hepak.HeCalc('H', 0, 'T', OneKPotTemperature, 'SV', 0., 1)
  HeLiquidEnthalpy = hepak.HeCalc('H', 0, 'T', OneKPotTemperature, 'SL', 0., 1)
  liquidFraction = (HeInletEnthalpy - HeVaporEnthalpy)/(HeLiquidEnthalpy - HeVaporEnthalpy)
  return (He3HeatLoad + IPHeatLoad)/hepak.HeCalc(7, 0, 'T', OneKPotTemperature, 'SL', 0., 1)/liquidFraction # He mass flow = heat load / latent heat(T) / (1 - vapor mass fraction(H,T))

# calculate He3 flow through JT valve required to remove heatLoad
def He3Flow(inletTemperature, inletPressure, temperature, heatLoad):
  inletDensity = he3pak.He3Density(inletPressure, inletTemperature)
  inletEnthalpy = he3pak.He3Prop(6, inletDensity, inletTemperature)
  vaporEnthalpy = he3pak.He3SaturatedVaporProp(6, temperature)
  liquidEnthalpy = he3pak.He3SaturatedLiquidProp(6, temperature)
  liquidFraction = (inletEnthalpy - vaporEnthalpy)/(liquidEnthalpy - vaporEnthalpy)
  latentHeat = he3pak.He3SaturatedLiquidProp(13, temperature)
  liquidFlow = heatLoad/latentHeat
  return liquidFlow/liquidFraction

# calculate He3 temperature when pumping away a certain gas flow
# pumpingSpeed dV/dt is assumed to be in m3/h
def He3Temperature(He3Flow, pumpingPressureDrop):
#  specificGasConstant = 2756.7579 # gas constant/molar mass (J/kg/K)
#  inletPressure = scipy.optimize.root_scalar(lambda P: he3pak.He3Prop(5, he3pak.He3Density(P, inletTemperature), inletTemperature)*gasFlow*specificGasConstant*inletTemperature/(pumpingSpeed/3600) - P, bracket = (0., 1e5)) # solve Z(P, T)*dm/dt*R_s*T/(dV/dt) = P for P
#  if not inletPressure.converged:
#    print('Could not determine pump inlet pressure')
#  pumpInletPressure = He3Flow*specificGasConstant*pumpInletTemperature/(pumpingSpeed/3600) # P = dm/dt R_s T / (dV/dt)
  pumpInletPressure = scipy.optimize.root_scalar(lambda P: pumpData(P)*2. - He3Flow, bracket = (15., 12000.)) # find inlet presure for given flow in pump data (using two parallel pumps with performance scaled by 75% to He3)
  if not pumpInletPressure.converged:
    print("Could not calculate He3 pump inlet pressure!")

  return he3pak.He3SaturatedTemperature(pumpInletPressure.root + pumpingPressureDrop) # T(P + dP)


# calculate temperature of HEX1 by interpolating measured He3 boiling curve
He3boilingData = HEXdata.loadHe3boilingData()
def HEX1Temperature(T_He3, HEX1length, HEX1diameter, HEX1surface, heatLoad):
  if heatLoad <= 0.:
    return T_He3
#  numberFins = HEX1length/finPitch
#  area = HEX1diameter*math.pi*HEX1length/2 + (HEX1diameter + 2*finHeight)*math.pi*HEX1length/2 + ((HEX1diameter/2 + finHeight)**2 - (HEX1diameter/2)**2)*math.pi*numberFins*2
  area = HEX1surface * HEX1length
  q = heatLoad/area
  
  dT = He3boilingData[0](T_He3, q)
  if math.isnan(dT):
    print('He3 temperature outside valid boiling temperature range!')
    dT = He3boilingData[1](T_He3, q)
  return T_He3 + dT

  
# calculate temperature of HeII at HEX1 from empirical Kapitza conductance
def He4Temperature(T_Cu, HEX1length, HEX1diameter, heatLoad):
  A = HEX1diameter*math.pi*HEX1length
  h = 900.*T_Cu**3
  return T_Cu + heatLoad/A/h
#  return (T_Cu**3.46 + heatLoad/A/460.)**(1./3.46) # for polished/oxidized Cu, see van Sciver table 7.4


# calculate temperature profile of He-II in conduction channel from Gorter-Mellink equation
def HeIItemperature(x, T_low, pressure, heatLoad, channelDiameter):
  if heatLoad == 0.:
    return T_low, T_low
	
  A = math.pi*channelDiameter**2/4.
  
  def dTdx(x, T):
    if T[0] < hepak.HeConst(3):
      print('T_4He {0:.3g} outside HEPAK range!'.format(T[0]))
      return 0.
    k = hepak.HeCalc(38, 0, 'P', pressure, 'T', T[0], 1) # heat conductivity from HEPAK
    if k > 0:
      return (heatLoad/A)**3 / k
    else:
      print('Conductivity <= 0!')
      return 0.
	  
  def dTdx2(x, T):
    if T[0] < hepak.HeConst(3):
      print('T_4He {0:.3g} outside HEPAK range!'.format(T[0]))
      return 0.
    g_lambda = hepak.HeCalc('D', 0, 'P', pressure, 'T', T[0], 1)**2 * 1559.**4 * hepak.HeConst(9)**3 / 1450.
    f_inverse = g_lambda * ((T[0]/hepak.HeConst(9))**5.7 * (1 - (T[0]/hepak.HeConst(9))**5.7))**3 # heat conductivity from vanSciver
    return (heatLoad/A)**3 / f_inverse
	  
  result = scipy.integrate.solve_ivp(dTdx, (0., x), [T_low]) # integrate ODE dT/dx = (q/A)^3 / k(T)
  result2 = scipy.integrate.solve_ivp(dTdx2, (0., x), [T_low])
  
  channelVolume = A*x
  bottleVolume = 0.034
  meanTemperature = (numpy.mean(result.y[0])*channelVolume + result.y[0][-1]*bottleVolume)/(channelVolume + bottleVolume)
  meanLifetime = (numpy.mean(numpy.reciprocal(numpy.power(result.y[0], 7)))/0.016 * channelVolume + 1./result.y[0][-1]**7/0.016 * bottleVolume)/(channelVolume + bottleVolume)
  return result.y[0][-1], result2.y[0][-1], meanTemperature, meanLifetime


# x: [3He flow, 1K pot temperature, 4He temperature at HEX1]
def equationSet(x, parameters):
  heatLoad = parameters['Beam heating']['value'] + parameters['Static heat']['value'] + HeCoolingLoad(parameters['Isopure He flow']['value'], parameters['Isopure He pressure']['value'], x[1], x[2])
  T_He3 = He3Temperature(x[0], parameters['3He pressure drop']['value'])

  flow = He3Flow(x[1], parameters['3He inlet pressure']['value'], T_He3, heatLoad)

  reservoirTemperature = hepak.HeCalc('T', 0, 'P', parameters['He reservoir pressure']['value'], 'SV', 0., 1)
  OneKPotFlow = OneKPotEvaporation(x[0], parameters['3He inlet pressure']['value'], parameters['HEX4 exit temperature']['value'], \
                                   parameters['Isopure He flow']['value'], parameters['Isopure He pressure']['value'], reservoirTemperature, \
								   parameters['He reservoir pressure']['value'], parameters['HEX5 exit temperature']['value'], x[1])
  T_1Kpot = OneKPotTemperature(OneKPotFlow, parameters['He pressure drop']['value'])

  T_HEX1 = HEX1Temperature(T_He3, parameters['HEX1 length']['value'], parameters['Channel diameter']['value'], parameters['HEX1 surface']['value'], heatLoad)
  T_He4 = He4Temperature(T_HEX1, parameters['HEX1 length']['value'], parameters['Channel diameter']['value'], heatLoad)
  
  return [flow - x[0], T_1Kpot - x[1], T_He4 - x[2]]
  
# solve equation set for UCN source and return them with some other calculated properties
def calcUCNSource(parameters):
  sol = scipy.optimize.root(equationSet, x0 = [0.0008, 1.8, 1.], args = parameters, method = 'hybr')
  if not sol.success:
    print('Solution did not converge!')
    return
  result = {}
  result['3He flow'] = sol.x[0]
  result['1K pot temperature'] = sol.x[1]
  result['T_4He'] = sol.x[2]
  result['He reservoir temperature'] = hepak.HeCalc('T', 0, 'P', parameters['He reservoir pressure']['value'], 'SV', 0., 1)
  result['1K pot flow'] = OneKPotEvaporation(result['3He flow'], parameters['3He inlet pressure']['value'], parameters['HEX4 exit temperature']['value'], \
                                             parameters['Isopure He flow']['value'], parameters['Isopure He pressure']['value'], result['He reservoir temperature'], \
											 parameters['He reservoir pressure']['value'], parameters['HEX5 exit temperature']['value'], result['1K pot temperature'])
  result['He reservoir flow'] = HeReservoirEvaporation(result['3He flow'], parameters['3He inlet pressure']['value'], parameters['He reservoir inlet temperature']['value'], \
                                                       parameters['Isopure He flow']['value'], parameters['Isopure He pressure']['value'], parameters['Isopure He inlet temperature']['value'], \
													   parameters['He reservoir pressure']['value'])
  result['He consumption'] = (result['He reservoir flow'] + result['1K pot flow'])/hepak.HeCalc('D', 0, 'P', 1013e2, 'SL', 0., 1)*1000*3600
  result['20K shield flow'] = shieldFlow(parameters['Isopure He flow']['value'], parameters['Isopure He pressure']['value'], parameters['He reservoir pressure']['value'])[0]
  result['100K shield flow'] = shieldFlow(parameters['Isopure He flow']['value'], parameters['Isopure He pressure']['value'], parameters['He reservoir pressure']['value'])[1]
  result['Shield consumption'] = max(result['20K shield flow'], result['100K shield flow'])/hepak.HeCalc('D', 0, 'P', 1013e2, 'SL', 0., 1)*1000*3600
 
  result['T_3He'] = He3Temperature(result['3He flow'], parameters['3He pressure drop']['value'])
  heatLoad = parameters['Beam heating']['value'] + parameters['Static heat']['value'] + HeCoolingLoad(parameters['Isopure He flow']['value'], parameters['Isopure He pressure']['value'], result['1K pot temperature'], result['T_4He'])
  result['T_Cu'] = HEX1Temperature(result['T_3He'], parameters['HEX1 length']['value'], parameters['Channel diameter']['value'], parameters['HEX1 surface']['value'], heatLoad)

  result['HeII vapor pressure'] = 0.
  result['HeII pressure head'] = 0.
  result['T_HeII'] = result['T_4He'], result['T_4He'], result['T_4He'], 1./result['T_4He']**7/0.016
  if parameters['Beam heating']['value'] > 0. and result['T_4He'] > hepak.HeConst(3):
    result['HeII vapor pressure'] = hepak.HeCalc('P', 0, 'T', result['T_4He'], 'SL', 0., 1)
    result['HeII pressure head'] = hepak.HeCalc('D', 0, 'T', result['T_4He'], 'SL', 0., 1)*9.81*parameters['HeII overfill']['value']
    result['T_HeII'] = HeIItemperature(parameters['Channel length']['value'], result['T_4He'], result['HeII vapor pressure'] + result['HeII pressure head'], parameters['Beam heating']['value'], parameters['Channel diameter']['value'])
  return result


parameters = {
'Beam heating': 					{'value': 8.1,		'range': (0., 10.), 		'unit': 'W'}, # Beam on
#'Beam heating': 					{'value': 0.,		'range': (0., 1.), 			'unit': 'W'}, # Beam off
'Static heat': 						{'value': 1., 		'range': (0.1, 2.),         'unit': 'W'},
#'3He pumping speed':				{'value': 4300.,	'range': (300., 10000.), 	'unit': r'm$^{3}$/h'}, # 1.14g/s @ 5.85 torr (Busch proposal)
#'He pumping speed':					{'value': 2000.,	'range': (500., 5000.), 	'unit': r'm$^{3}$/h'}, # 1.2g/s @ 9.9 torr (Busch proposal)
#'Pump inlet temperature':			{'value': 290.,		'range': (100., 300.), 		'unit': 'K'},
'3He pressure drop': 				{'value': 100., 	'range': (0., 1000.), 		'unit': 'Pa'},
#'3He pressure drop': 				{'value': 20000., 	'range': (0., 50000.), 		'unit': 'Pa'}, # Standby mode
'He pressure drop': 				{'value': 100., 	'range': (0., 1000.), 		'unit': 'Pa'},
#'He pressure drop': 				{'value': 10000., 	'range': (0., 10000.), 		'unit': 'Pa'}, # Standby mode
'He reservoir inlet temperature': 	{'value': 10., 		'range': (6., 12.,), 		'unit': 'K'},
'HEX4 exit temperature':			{'value': 2.8,		'range': (2., 3.5),			'unit': 'K'},
'HEX5 exit temperature':			{'value': 2.8,		'range': (2., 3.5),			'unit': 'K'},
#'1K pot inlet temperature': 		{'value': 2.8, 		'range': (2., 3.5),			'unit': 'K'},
'3He inlet pressure': 				{'value': 50000.,	'range': (20000., 100000.), 'unit': 'Pa'},
'Channel diameter': 				{'value': 0.148, 	'range': (0.12, 0.2), 		'unit': 'm'},
'Channel length': 					{'value': 2.5, 		'range': (1., 4.), 			'unit': 'm'},
'HEX1 length': 						{'value': 0.6, 		'range': (0.1, 0.8), 		'unit': 'm'},
'HEX1 surface': 					{'value': 1.67,	 	'range': (0.21, 2.),		'unit': r'm$^{2}$/m'},
'HeII overfill': 					{'value': 0.05, 	'range': (0.02, 0.2), 		'unit': 'm'},
'He reservoir pressure': 			{'value': 1.2e5, 	'range': (900e2, 1500e2), 	'unit': 'Pa'},
'Isopure He flow':					{'value': 0.,		'range': (0., 0.),			'unit': 'kg/s'},
#'Isopure He flow':                  {'value': 0.00014,  'range': (0., 0.0005),      'unit': 'kg/s'}, # isopure condensation
'Isopure He inlet temperature':     {'value': 20.,      'range': (10., 40.),        'unit': 'K'},
'Isopure He pressure':              {'value': 100e2,    'range': (10e2, 500e2),     'unit': 'Pa'}
}

result = calcUCNSource(parameters)

print('3He flow: {0:.4g} g/s'.format(result['3He flow']*1000))
print('He reservoir evaporation: {0:.4g} g/s'.format(result['He reservoir flow']*1000))
print('1K pot: {0:.4} K, {1:.4g} g/s'.format(result['1K pot temperature'], result['1K pot flow']*1000))
print('He consumption: {0:.3g} L/h'.format(result['He consumption']))
print('Required shield flows: {0:.3g} g/s, {1:.3g} g/s'.format(result['20K shield flow']*1000, result['100K shield flow']*1000))
#print('He-II pressure: {0:.4g} + {1:.4g} Pa'.format(result['HeII vapor pressure'], result['HeII pressure head']))
print('T: {0:.4g} -> {1:.4g} -> {2:.4g} -> {3[0]:.4g} to {3[1]:.4g} K'.format(result['T_3He'], result['T_Cu'], result['T_4He'], result['T_HeII']))
print('UCN lifetime: {0:.4g}'.format(result['T_HeII'][3]))

# scan parameter ranges and plot temperatures and flows
plotRows = 3
plotCols = 6
fig, axes = plt.subplots(plotRows, plotCols, figsize=(plotCols*5,plotRows*5))
axes2 = axes.copy()
axes3 = axes.copy()
fig.set_tight_layout(True)
for i, p in enumerate(parameters):
  params = {}
  for pp in parameters:
    params[pp] = dict(parameters[pp]) # make copy of original parameter set that we can modify
  print(p)
  x = []
  y = []
  y2 = []
  y2_2 = []
  y3 = []
  y4 = []
  y5 = []
  y6 = []
  y7 = []
  y8 = []
  y9 = []
  y10 = []
  for j in range(11):
    value = params[p]['range'][0] + float(j)/10.*(params[p]['range'][1] - params[p]['range'][0])
    print(value, params[p]['unit'])
    params[p]['value'] = value
    res = calcUCNSource(params)
    if res:
      print((res['He reservoir flow'] + res['1K pot flow'])/hepak.HeCalc('D', 0, 'P', 1000e2, 'SL', 0, 1)*1000.*3600., 'L/h')
      x.append(params[p]['value'])
      y.append(res['He reservoir flow']*1000)
      y2.append(res['T_HeII'][0])
      y2_2.append(res['T_HeII'][1])
      y3.append(res['T_3He'])
      y4.append(res['T_4He'])
      y5.append(res['3He flow']*1000)
      y6.append(res['1K pot flow']*1000)
      y7.append(max(res['20K shield flow'], res['100K shield flow'])*1000)
      y8.append(res['1K pot temperature'])
      y9.append(res['T_Cu'])
      y10.append(res['T_HeII'][3])
  ax = axes[i % plotRows][int(i/plotRows)]
  axes[0][0].get_shared_y_axes().join(axes[0][0], ax)
  ax.plot(x, y, color = 'tab:green', label = 'He reservoir')
  ax.plot(x, y6, color = 'tab:olive', label = '1K pot')
  ax.plot(x, y7, color = 'tab:cyan', label = 'Req. for shields')
  ax.plot(x, y5, color = 'tab:orange', label = '3He')
  ax.set_xlabel(p + ' (' + params[p]['unit'] + ')')
  ax.set_ylabel('Gas flow (g/s)')

  axes2[i % plotRows][int(i/plotRows)] = ax.twinx()
  ax2 = axes2[i % plotRows][int(i/plotRows)]
  axes2[0][0].get_shared_y_axes().join(axes2[0][0], ax2)
  ax2.plot(x, y8, color = 'tab:pink', label = '1K pot')
  ax2.plot(x, y3, color = 'tab:red', label = '3He')
  ax2.plot(x, y9, color = 'tab:gray', label = 'HEX1')
  ax2.plot(x, y4, color = 'tab:purple', label = '4He at HEX1')
  ax2.fill_between(x, y2, y2_2, color = 'tab:blue', label = '4He in bottle', alpha = 0.2)
  ax2.set_ylabel('Temperature (K)')
  
#  axes3[i % plotRows][int(i/plotRows)] = ax.twinx()
#  ax3 = axes3[i % plotRows][int(i/plotRows)]
#  axes3[0][0].get_shared_y_axes().join(axes3[0][0], ax3)
#  ax3.set_ylim(10., 60.)
#  ax3.plot(x, y10, color = 'tab:brown', label = 'UCN lifetime')
  
  plt.axvline(parameters[p]['value'], 0, 1, color = 'k', dashes = (4, 2))
  if i == 9:
    ax2.legend(title = 'Temperatures', loc = 'upper right')
    ax.legend(title = 'Flows', loc = 'upper left')
  
while True:
  try:
    plt.savefig('UCNsource.pdf')
  except PermissionError:
    input('Could not save plots. If the target file is still opened in another application please close it. Press ENTER to try again.')
  else:
    break
plt.close()