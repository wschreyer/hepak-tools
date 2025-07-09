import hepak
import he3pak
import math
import numpy
import scipy.integrate
import scipy.optimize
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
def HeReservoirEvaporation(He3Flow, He3Pressure, He3InletTemperature, IPFlow, IPPressure, IPInletTemperature, staticLoad, reservoirPressure):
  reservoirTemperature = hepak.HeCalc('T', 0, 'P', reservoirPressure, 'SV', 0., 1)
  He3HeatLoad = He3CoolingLoad(He3Flow, He3Pressure, He3InletTemperature, reservoirTemperature)
  IPHeatLoad = HeCoolingLoad(IPFlow, IPPressure, IPInletTemperature, reservoirTemperature)
  return (He3HeatLoad + IPHeatLoad + staticLoad)/hepak.HeCalc(7, 0, 'P', reservoirPressure, 'SL', 0., 1) # evaporation flow = heatLoad/latent heat
  
# calculate required flow to cool 20K and 100K shields, with additional heat load by isopure He flowing over shields
# static heat loads on shields are assumed fixed
def shieldFlow(temp20K, temp100K, inlet20K, inlet100K, isopureFlow, isopurePressure):
  IPHeatLoad_100 = HeCoolingLoad(isopureFlow, isopurePressure, 300., temp100K)
  IPHeatLoad_20 = HeCoolingLoad(isopureFlow, isopurePressure, temp100K, temp20K)
  flow_20 = HeCoolingFlow(4.3 + 0.7 + 1.1 + 2. + 6.6 + IPHeatLoad_20, isopurePressure, inlet20K, temp20K)
  flow_100 = HeCoolingFlow(49. + 8.9 + 5.7 + 14. + IPHeatLoad_100, isopurePressure, inlet100K, temp100K)
  return flow_20, flow_100

pumpDataHe4 = pumpdata.loadPumpData('pumpdata/Busch_2stage_He4_torqueControl_70Hz.csv')
	  
# calculate 1K pot temperature when pumping away a certain gas flow
def OneKPotTemperature(HeFlow, pumpingPressureDrop):
  HeFlow = max(0., min(HeFlow, 0.1))
#  inletPressure = scipy.optimize.root_scalar(lambda P: hepak.HeCalc(5, 0, 'P', P, 'T', pumpInletTemperature, 0)*gasFlow*hepak.HeConst(4)*inletTemperature/(pumpingSpeed/3600) - P, bracket = (1., 1e5)) # He pressure = Z(P, T)*dm/dt*R_s*T/(dV/dt)
#  if not inletPressure.converged:
#    print('Could not determine pump inlet pressure')
#  pumpInletPressure = HeFlow*hepak.HeConst(4)*pumpInletTemperature/(pumpingSpeed/3600) # P = dm/dt R_s T / (dV/dt)

  pumpInletPressure = pumpDataHe3(HeFlow) #scipy.optimize.root_scalar(lambda P: pumpDataHe4(P) - HeFlow, bracket = (1., 20000.)) # find inlet pressure for given flow in pump data
#  if not pumpInletPressure.converged:
#    print("Could not calculate He4 pump inlet pressure!")
  T = hepak.HeCalc('T', 0, 'P', pumpInletPressure + pumpingPressureDrop, 'SL', 0., 1) # T_sat(P + dP)
  return T

# calculate evaporation rate from 1K pot at certain temperature and He3 flow
def OneKPotEvaporation(He3Flow, He3Pressure, He3InletTemperature, IPFlow, IPPressure, IPInletTemperature, HeInletPressure, HeInletTemperature, staticLoad, OneKPotTemperature):
  He3HeatLoad = He3CoolingLoad(He3Flow, He3Pressure, He3InletTemperature, OneKPotTemperature)
  IPHeatLoad = HeCoolingLoad(IPFlow, IPPressure, IPInletTemperature, OneKPotTemperature)

  HeInletEnthalpy = hepak.HeCalc('H', 0, 'T', HeInletTemperature, 'P', HeInletPressure, 1)
  HeVaporEnthalpy = hepak.HeCalc('H', 0, 'T', OneKPotTemperature, 'SV', 0., 1)
  HeLiquidEnthalpy = hepak.HeCalc('H', 0, 'T', OneKPotTemperature, 'SL', 0., 1)
  liquidFraction = (HeInletEnthalpy - HeVaporEnthalpy)/(HeLiquidEnthalpy - HeVaporEnthalpy)
  return (He3HeatLoad + IPHeatLoad + staticLoad)/hepak.HeCalc(7, 0, 'T', OneKPotTemperature, 'SL', 0., 1)/liquidFraction # He mass flow = heat load / latent heat(T) / (1 - vapor mass fraction(H,T))

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

pumpDataHe3 = pumpdata.loadPumpData('pumpdata/Busch_2stage_He3_torqueControl_90Hz.csv')
	  
# calculate He3 temperature when pumping away a certain gas flow
# pumpingSpeed dV/dt is assumed to be in m3/h
def He3Temperature(He3Flow, pumpingPressureDrop):
#  specificGasConstant = 2756.7579 # gas constant/molar mass (J/kg/K)
#  inletPressure = scipy.optimize.root_scalar(lambda P: he3pak.He3Prop(5, he3pak.He3Density(P, inletTemperature), inletTemperature)*gasFlow*specificGasConstant*inletTemperature/(pumpingSpeed/3600) - P, bracket = (0., 1e5)) # solve Z(P, T)*dm/dt*R_s*T/(dV/dt) = P for P
#  if not inletPressure.converged:
#    print('Could not determine pump inlet pressure')
#  pumpInletPressure = He3Flow*specificGasConstant*pumpInletTemperature/(pumpingSpeed/3600) # P = dm/dt R_s T / (dV/dt)
  pumpInletPressure = pumpDataHe3(He3Flow) #scipy.optimize.root_scalar(lambda P: pumpDataHe3(P)*2. - He3Flow, bracket = (1., 20000.)) # find inlet presure for given flow in pump data (using two parallel pumps)
#  if not pumpInletPressure.converged:
#    print("Could not calculate He3 pump inlet pressure!")

  return he3pak.He3SaturatedTemperature(pumpInletPressure + pumpingPressureDrop) # T(P + dP)


# calculate temperature of HEX1 by interpolating measured He3 boiling curve
He3boilingData = HEXdata.loadHe3boilingData()
def HEX1TemperatureLow(T_He3, HEX1length, HEX1diameter, HEX1surface, heatLoad):
  if heatLoad <= 0.:
    return T_He3
    
#  numberFins = HEX1length/finPitch
#  area = HEX1diameter*math.pi*HEX1length/2 + (HEX1diameter + 2*finHeight)*math.pi*HEX1length/2 + ((HEX1diameter/2 + finHeight)**2 - (HEX1diameter/2)**2)*math.pi*numberFins*2
  area = HEX1surface * HEX1length
  q = heatLoad/area
  
#  kG = 35
#  h = kG/2.6 * 20 * T_He3**3 # Simple Kapitza model (1.2 to 2.6 times higher Kapitza resistance in He3 compared to He-II)
#  return T_He3 + q/h
  
  dT = He3boilingData[0](T_He3, q)
  if math.isnan(dT):
    print('He3 temperature outside valid boiling temperature range!')
    dT = He3boilingData[1](T_He3, q)
  return T_He3 + dT

# Calculate temperature gradient across copper of HEX1 (based on 2D simulation with FEMM)
def HEX1TemperatureHigh(T_HEX1_low, HEX1length, HEX1conductivity, heatLoad):
  return T_HEX1_low + heatLoad/HEX1conductivity/HEX1length  # temperature difference is proportional to heat load, indirectly proportional to length
  
# calculate temperature of HeII at HEX1 from empirical Kapitza conductance
def HeIITemperatureLow(T_Cu, HEX1length, HEX1diameter, heatLoad, KapitzaKG, HeIIfillLevel):
  HEX1fillLevel = HeIIfillLevel - (0.18 + 0.085 - HEX1diameter/2.)
  if HEX1fillLevel >= HEX1diameter:
    filledPortion = 1.
  else:
    filledPortion = math.acos(1. - HEX1fillLevel/(HEX1diameter/2.))/math.pi
  A = HEX1diameter*math.pi*HEX1length*filledPortion
  h = KapitzaKG * 20 * T_Cu**3 # use kG def (0.61 factor for Ni plating)
  return T_Cu + heatLoad/A/h
#  return (T_Cu**3.46 + heatLoad/A/460.)**(1./3.46) # for polished/oxidized Cu, see van Sciver table 7.4


# calculate position where heat flux becomes zero (assuming constant heat flux removed along HEX1)
def HeIIlowestTemperaturePosition(HEX1length, beamHeating, staticHeat, funneledHeat):
  upstreamHeating = beamHeating + staticHeat # total heat introduced upstream of HEX1 (beam heating + upstream fraction of static heat)
  downstreamHeating = funneledHeat # total heat introduced downstream of HEX1 (funneled heating + upstream fraction of static heat)
  return HEX1length/2.*(downstreamHeating - upstreamHeating)/(beamHeating + staticHeat + funneledHeat) # position where heat flux becomes zero (assuming constant heat flux removed along HEX1)


# calculate cross section of a circular segment (partially filled circular cross section)
def partialCrossSection(radius, fillLevel):
  if fillLevel > radius*2.:
    return radius**2*math.pi
  else:
    return radius**2*math.acos(1. - fillLevel/radius) - (radius - fillLevel)*math.sqrt(radius**2 - (radius - fillLevel)**2)

# calculate heat flux in He-II, x is distance from center of HEX1
def HeIIheatFlux(x, beamHeating, bulbRadius, bulbLength, bulbChannelOffset, channelDiameter, channelLength, HEX1diameter, HEX1length, staticHeat, funneledHeat, funnelDiameter, HeIIfillLevel):
  channelFillLevel = HeIIfillLevel - (bulbRadius + bulbChannelOffset - channelDiameter/2.) # HeII fill level is measured from bottom of bulb
  if channelFillLevel <= 0:
    print('HeII fill level insufficient to reach channel')
    return 0.
  channelCrossSection = partialCrossSection(channelDiameter/2., channelFillLevel) # cross section of channel filled with HeII
  HEX1fillLevel = HeIIfillLevel - (bulbRadius + bulbChannelOffset - HEX1diameter/2.)
  HEX1crossSection = partialCrossSection(HEX1diameter/2., HEX1fillLevel) # cross section of HEX1 filled with HeII
  funnelFillLevel = HeIIfillLevel - (bulbRadius + bulbChannelOffset - funnelDiameter/2.)
  if funnelFillLevel > 0.3: # we're assuming the S-bend has a curvature radius of 15cm
    print('HeII overfilled')
    return 0.
  elif funnelFillLevel <= 0:
    funnelLength = 0.
    funnelCrossSection = 0.
  else:
    funnelLength = math.acos(1. - funnelFillLevel/0.15)*0.15 # approximate center line length of S-bend filled with HeII
    funnelCrossSection = partialCrossSection(funnelDiameter/2., funnelFillLevel) # roughly approximate cross section of S-bend filled with HeII
  x0 = HeIIlowestTemperaturePosition(HEX1length, beamHeating, staticHeat, funneledHeat) 
  staticHeatLength = bulbLength + channelLength # total length that receives static heating
  if x < -HEX1length/2. - funnelLength:
    return 0.
  elif x < -HEX1length/2.:
    return funneledHeat/funnelCrossSection # downstream of HEX1 heat flux consists of funneled heat plus static heat added along funnel
  elif x < x0:
    return funneledHeat*(x0 - x)/(x0 + HEX1length/2.)/HEX1crossSection # along HEX1 downstream heat flux is removed at a constant rate
  elif x < HEX1length/2.:
    return -(beamHeating + staticHeat)*(x - x0)/(HEX1length/2. - x0)/HEX1crossSection # along HEX1 upstream heat flux is removed at a constant rate
  elif x < HEX1length/2. + channelLength:
    return -(beamHeating + staticHeat*(HEX1length/2. + channelLength + bulbLength - x)/staticHeatLength)/channelCrossSection # in channel heat flux consists of beam heating plus static heat load added in channel and bulb upstream of x
  elif x < HEX1length/2. + channelLength + bulbLength:
    distanceFromBulbCenter = x - HEX1length/2. - channelLength - bulbLength/2. # model cross section of bulb as hemispheres with cylindrical piece in between
    if abs(distanceFromBulbCenter) < bulbLength/2. - bulbRadius:
      bulbCrossSection = partialCrossSection(bulbRadius, HeIIfillLevel) # cylindrical piece in center
    elif distanceFromBulbCenter < 0:
      distanceFromBulbCenter = abs(distanceFromBulbCenter) - (bulbLength/2. - bulbRadius)
      radius = math.sqrt(bulbRadius**2 - distanceFromBulbCenter**2)
      bulbCrossSection = max(channelCrossSection, partialCrossSection(radius, HeIIfillLevel - (bulbRadius - radius))) # cross section of sphere at given distance from center, but use channel cross section when downstream hemisphere shrinks below that
    else:
      distanceFromBulbCenter = abs(distanceFromBulbCenter) - (bulbLength/2. - bulbRadius)
      radius = math.sqrt(bulbRadius**2 - distanceFromBulbCenter**2)
      bulbCrossSection = partialCrossSection(radius, HeIIfillLevel - (bulbRadius - radius)) # cross section of sphere at given distance from center
    return -(beamHeating*(HEX1length/2. + channelLength + bulbLength - x)/bulbLength + staticHeat*(HEX1length/2. + channelLength + bulbLength - x)/staticHeatLength)/bulbCrossSection # in bulb heat flux consists of beam and static heating added in bulb upstream of x
  else:
    return 0.


# calculate temperature profile according to Gorter-Mellink equation, starting at temperature T_low and position x0
# along channel with given heat flux and length filled with He-II at given pressure, and transporting given heat load
# conductivity is either calculated using 'HEPAK' or 'VanSciver'
# returns interpolation function for temperature profile
def GorterMellink(T_low, pressure, heatFluxModel, x0, length, conductivityModel):
  def dTdx(x, T):
    heatFlux = heatFluxModel(x)
    if heatFlux == 0.:
      return 0.
    if T[0] < hepak.HeConst(3):
      print('T_4He {0:.3g} outside HEPAK range!'.format(T[0]))
      return 0.
    
    if conductivityModel == 'HEPAK':
      k = hepak.HeCalc(38, 0, 'P', pressure, 'T', T[0], 1) # heat conductivity from HEPAK
      if k > 0:
        return -heatFlux**3 / k
      else:
        print('HEPAK conductivity <= 0!')
        return 0.
    elif conductivityModel == 'VanSciver':	  
      g_lambda = hepak.HeCalc('D', 0, 'P', pressure, 'T', T[0], 1)**2 * 1559.**4 * hepak.HeConst(9)**3 / 1450.
      f_inverse = g_lambda * ((T[0]/hepak.HeConst(9))**5.7 * (1 - (T[0]/hepak.HeConst(9))**5.7))**3 # heat conductivity from vanSciver
      return -heatFlux**3 / f_inverse
    else:
      raise 'Invalid Gorter-Mellink model!'
  
  result = scipy.integrate.solve_ivp(dTdx, (x0, length), [T_low], max_step = 0.05, dense_output = True, rtol = 1e-5, atol = 1e-8) # integrate ODE dT/dx = (q/A)^3 / k(T)
  if not result.success:
    print(result.message)
  return result.sol
 

# calculate temperature profile of He-II in conduction channel from Gorter-Mellink equation
def HeIItemperatureHigh(channelLength, HEX1length, T_low, pressure, beamHeating, staticHeat, funneledHeat, channelDiameter, HEX1diameter, HeIIfillLevel):
  heatFluxModel = lambda x: HeIIheatFlux(x, beamHeating, 0.18, 0.386, 0.085, channelDiameter, channelLength, HEX1diameter, HEX1length, staticHeat, funneledHeat, 0.0955, HeIIfillLevel)
  x0 = HeIIlowestTemperaturePosition(HEX1length, beamHeating, staticHeat, funneledHeat)
  profileHEPAK = GorterMellink(T_low, pressure, heatFluxModel, x0, channelLength, 'HEPAK')
  profileVanSciver = GorterMellink(T_low, pressure, heatFluxModel, x0, channelLength, 'VanSciver')
  
  return profileVanSciver(channelLength)[0], profileHEPAK(channelLength)[0]


# x: [3He flow, 1K pot temperature, 4He temperature at HEX1]
def equationSet(x, parameters):
  heatLoad = parameters['Beam heating'] + parameters['He-II static heat'] + parameters['He-II funnel heat'] + HeCoolingLoad(parameters['Isopure He flow'], parameters['Isopure He pressure'], x[1], x[2])
  T_He3 = He3Temperature(x[0], parameters['3He pressure drop'])

  flow = He3Flow(x[1], parameters['3He inlet pressure'], T_He3, heatLoad)

  reservoirTemperature = hepak.HeCalc('T', 0, 'P', parameters['He reservoir pressure'], 'SV', 0., 1)
  OneKPotFlow = OneKPotEvaporation(x[0], parameters['3He inlet pressure'], parameters['HEX4 exit temperature'], \
                                   parameters['Isopure He flow'], parameters['Isopure He pressure'], reservoirTemperature, \
								                   parameters['He reservoir pressure'], parameters['HEX5 exit temperature'], parameters['1K pot static load'], x[1])
  T_1Kpot = OneKPotTemperature(OneKPotFlow, parameters['He pressure drop'])

  T_HEX1_low = HEX1TemperatureLow(T_He3, parameters['HEX1 length'], parameters['HEX1 diameter'], parameters['HEX1 3He surface'], heatLoad)
  T_HEX1_high = HEX1TemperatureHigh(T_HEX1_low, parameters['HEX1 length'], parameters['HEX1 conductivity'], heatLoad)
  T_HeII_low = HeIITemperatureLow(T_HEX1_high, parameters['HEX1 length'], parameters['Channel diameter'], heatLoad, parameters['HEX1 Kapitza kG'], parameters['HeII fill level'])
  
  return [flow - x[0], T_1Kpot - x[1], T_HeII_low - x[2]]
  
# solve equation set for UCN source and return them with some other calculated properties
def calcUCNSource(parameters):
  sol = scipy.optimize.root(equationSet, x0 = [0.0008, 1.8, 1.], args = parameters, method = 'hybr')
  if not sol.success:
    print('Solution did not converge!')
    return
  result = {}
  result['3He flow'] = sol.x[0]
  result['1K pot temperature'] = sol.x[1]
  result['T_HeII_low'] = sol.x[2]
  result['He reservoir temperature'] = hepak.HeCalc('T', 0, 'P', parameters['He reservoir pressure'], 'SV', 0., 1)
  result['1K pot flow'] = OneKPotEvaporation(result['3He flow'], parameters['3He inlet pressure'], parameters['HEX4 exit temperature'], \
                                             parameters['Isopure He flow'], parameters['Isopure He pressure'], result['He reservoir temperature'], \
											                       parameters['He reservoir pressure'], parameters['HEX5 exit temperature'], parameters['1K pot static load'], result['1K pot temperature'])
  result['He reservoir flow'] = HeReservoirEvaporation(result['3He flow'], parameters['3He inlet pressure'], parameters['He reservoir inlet temperature'], \
                                                       parameters['Isopure He flow'], parameters['Isopure He pressure'], parameters['20K shield temperature'], \
													                             parameters['He reservoir static load'], parameters['He reservoir pressure'])
  result['He consumption'] = (result['He reservoir flow'] + result['1K pot flow'])/hepak.HeCalc('D', 0, 'P', 1013e2, 'SL', 0., 1)*1000*3600
  result['20K shield flow'], result['100K shield flow'] = shieldFlow(parameters['20K shield temperature'], parameters['100K shield temperature'], result['He reservoir temperature'] + 0.1, \
                                                                     parameters['20K shield temperature'], parameters['Isopure He flow'], parameters['Isopure He pressure'])
  result['Shield consumption'] = max(result['20K shield flow'], result['100K shield flow'])/hepak.HeCalc('D', 0, 'P', 1013e2, 'SL', 0., 1)*1000*3600
 
  result['T_3He'] = He3Temperature(result['3He flow'], parameters['3He pressure drop'])
  result['3He vapor pressure'] = he3pak.He3SaturatedVaporProp(2, result['T_3He'])
  heatLoad = parameters['Beam heating'] + parameters['He-II static heat'] + parameters['He-II funnel heat'] + HeCoolingLoad(parameters['Isopure He flow'], parameters['Isopure He pressure'], result['1K pot temperature'], result['T_HeII_low'])
  result['He3 liquid fraction'] = heatLoad/result['3He flow']/he3pak.He3SaturatedLiquidProp(13, result['T_3He'])
  result['T_HEX1_low'] = HEX1TemperatureLow(result['T_3He'], parameters['HEX1 length'], parameters['HEX1 diameter'], parameters['HEX1 3He surface'], heatLoad)
  result['T_HEX1_high'] = HEX1TemperatureHigh(result['T_HEX1_low'], parameters['HEX1 length'], parameters['HEX1 conductivity'], heatLoad)

  result['HeII vapor pressure'] = 0.
  result['HeII pressure head'] = 0.
  result['T_HeII_high'] = result['T_HeII_low'], result['T_HeII_low']
  if parameters['Beam heating'] > 0. and result['T_HeII_low'] > hepak.HeConst(3):
    result['HeII vapor pressure'] = hepak.HeCalc('P', 0, 'T', result['T_HeII_low'], 'SL', 0., 1)
    result['HeII pressure head'] = hepak.HeCalc('D', 0, 'T', result['T_HeII_low'], 'SL', 0., 1)*9.81*parameters['HeII fill level']/2 # hydrostatic pressure at half the fill level (we use this to calculate the average HeII pressure)
    result['T_HeII_high'] = HeIItemperatureHigh(parameters['Channel length'], parameters['HEX1 length'], result['T_HeII_low'], result['HeII vapor pressure'] + result['HeII pressure head'], \
                                                parameters['Beam heating'], parameters['He-II static heat'], parameters['He-II funnel heat'], parameters['Channel diameter'], parameters['HEX1 diameter'], parameters['HeII fill level'])

  return result

