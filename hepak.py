import win32com.client
import sys
import os
import atexit

# THERMODYNAMIC PROPERTIES AVAILABLE IN HEPAK, see HEPAK user guide
# 0 ('X'): Quality = vapor mass fraction
# 1 ('P'): Pressure
# 2 ('T'): Temperature
# 3 ('D'): Density
# 4 ('V'): Specific Volume
# 5: Z = PV/RT
# 6: dP/dT (saturation line)
# 7: Latent Heat
# 8 ('S'): Entropy
# 9 ('H'): Enthalpy
# 10 ('A'): Helmholtz
# 11 ('U'): Internal Energy
# 12 ('G'): Gibbs Energy
# 13: Fugacity
# 14: Cp
# 15: Cv
# 16: Gamma = Cp/Cv
# 17: Expansivity = (T/V)(dV/dT)
# 18: Gruneisen parameter = (V/Cv)(dP/dT) at constant V
# 19: Isothermal compressibility = (1/D) dD/dP
# 20: Sound Velocity
# 21: Joule-Thomson Coefficient = dT/dP at constant H
# 22: dP/dD at constant T
# 23: dP/dT at constant D
# 24: V*dH/dV at constant P
# 25: Viscosity
# 26: Conductivity
# 27: Prandtl number
# 28: Thermal Diffusivity
# 29: Surface Tension
# 30: Dielectric constant - 1
# 31: Refractive index - 1
# 32: Isochoric dT to the lambda line dT(V)
# 33: Isobaric dT to the lambda line dT(P)
# 34: Superfluid density fraction RhoS/Rho
# 35: 2nd sound velocity
# 36: 4th sound velocity
# 37: Gorter-Mellink mutual friction constant
# 38: Superfluid thermal conductivity function
# 39: Lambda line temperature (isochoric to the state point) (T-Tlambda)

pwd = os.path.abspath(os.path.dirname(sys.argv[0]))

xl = win32com.client.gencache.EnsureDispatch('Excel.Application')
print('Loading HEPAK Excel Add-In. You may have to acknowledge a popup dialog. Launch Excel with HEPAK Add-In already installed to suppress dialog.')
hepak = xl.Workbooks.Open(os.path.join(pwd, 'hepak.xla'))

def close():
  if hepak:
    hepak.Close(True)
  
atexit.register(close)

# Return thermodynamic property of 4He in specified phase based on two input parameters
# property: integer -1..39, see above/HEPAK user guide
# phase: 0 (auto), 1 (single- or mixed-phase, non-saturated), 2 (single-phase only), 3 (mixed-phase only), 4 (saturated liquid), 5 (saturated vapor)
# input1: first input parameter, one of 1/'P', 2/'T', 3/'D', 4/'V', 8/'S', 9/'H', 11/'U', 12/'G', 0/'X', 'dT', 'M', 'SL', 'SV', 'L'
# value1: value of first input parameter
# input2: second input parameter, see input1
# value2: value of second input parameter
# units: units of input parameters and returned property, 1 (SI units), 2 (mixed SI-cgs), 3 (mixed SI-molar), 4 (imperial)
def HeCalc(property, phase, input1, value1, input2, value2, units):
  return xl.Application.Run('HeCalc', property, phase, input1, value1, input2, value2, units)

# supposed to return index of refraction for specified wavelength, but parameters are the same as for HeCalc???
# parameters see HeCalc
def HeRefrac(property, phase, input1, value1, input2, value2, units):
  return xl.Application.Run('HeRefrac', property, phase, input1, value1, input2, value2, units)

# return text description of fluid state
# messageID: ID returned by HeCalc when called with property -1
def HeMsg(messageID):
  return xl.Application.Run('HeMsg', messageID)

# return label for unit of selected property
# property: integer 0..40, see above/HEPAK user guide
# units: integer 1..4, see HeCalc
def HeUnit(parameter, units):
  return xl.Application.Run('HeUnit', parameter, units)

# returns a range of constants used in HEPAK calculations
# index: 1 (max valid pressure), 2 (max valid temperature), 3 (min valid temperature), 4 (gas constant),
#        5 (critical pressure), 6 (critical temperature), 7 (critical density),
#        8 (lambda point pressure), 9 (lambda point temperature), 10 (vapor density at lambda point), 11 (liquid density at lambda point),
#        12 (molecular weight), 13 (reference pressure for entropy scale), 14 (reserved), 15 (HEPAK version)
def HeConst(index):
  return xl.Application.Run('HeConst', index)

# return name of selected property
# property: integer 0..40, see HEPAK user guide
def HeProperty(property):
  return xl.Application.Run('HeProperty', property)

# return message indicating if input pair is valid for use in HeCalc
# input1: integer 0..40
# input2: integer 0..40
def HeValidate(input1, input2):
  return xl.Application.Run('HeValidate', input1, input2)
