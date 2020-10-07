import sys
import os

if any([sys.platform.startswith(os_name) for os_name in ['linux', 'darwin', 'freebsd']]):
  # Linux version:
  # import zugbruecke as ctypes
  import zugbruecke as ctypes
elif sys.platform.startswith('win'):
  import ctypes
else:
  # Handle unsupported platforms
  print("Unknown platform")
  exit(1)
  
# THERMODYNAMIC PROPERTIES AVAILABLE IN HE3PAK:
# 1: temperature
# 2: pressure
# 3: density
# 4: volm
# 5: compfactor
# 6: enthalpy
# 7: entropy
# 8: C_v
# 9: C_p
# 10: Helmholtz free energy
# 11: Gibbs free energy
# 12: sound
# 13: latent heat of evaporation
# 14: Joule-Thompson coefficient
# 15: inte
# 16: adiacomp
# 17: isocomp
# 18: volexp
# 19: isenexp
# 20: sndvirial
# 21: dBdT
# 22: trdvirial
# 23: dpdd
# 24: dpdt
# 25: dddt
# 26: dpdds
# 27: grun
# 28: thcon
# 29: viscosity
# 30: viscosity/density
# 31: surface tension
# 32: viscosity*C_p/thcon
# 33: thcon/density/C_p
# 34: dielectric constant
# 35: refractive index
# 36: 2nd sound velocity
# 37: 4th sound velocity
# 38: superfluid density fraction
# 39: Gorter Mellink mutual friction parameter
# 40: superfluid thermal conductivity

pwd = os.path.abspath(os.path.dirname(sys.argv[0]))

he3pak = ctypes.WinDLL(os.path.join(pwd, 'he3eos.dll'))

he3pak.DFPTdll.argtypes = (ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_int),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double))
he3pak.Fundtdll.argtypes = (ctypes.POINTER(ctypes.c_double*40),ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double))
he3pak.SatXFunTdll.argtypes = (ctypes.POINTER(ctypes.c_double*40),ctypes.POINTER(ctypes.c_double*40),ctypes.POINTER(ctypes.c_double))
he3pak.TIPsatdll.argtypes = (ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double))

# return density of He3 at specified pressure and temperature
def He3Density(pressure, temperature):
  density = ctypes.c_double(0.0)
  idid = ctypes.c_int(0)
  P = ctypes.c_double(pressure)
  T = ctypes.c_double(temperature)
  he3pak.DFPTdll(density, idid, P, T)
  return density.value
  
# returns thermodynamic property of He3 at specified density and temperature
def He3Prop(property, density, temperature):
  D = ctypes.c_double(density)
  T = ctypes.c_double(temperature)
  Xprop = (ctypes.c_double*40)()
  he3pak.Fundtdll(Xprop, D, T)
  return Xprop[property - 1]
  
# returns thermodynamic property of saturated He3 liquid and vapor at specified temperature
def He3SaturatedProp(property, temperature):
  T = ctypes.c_double(temperature)
  XpropL = (ctypes.c_double*40)()
  XpropV = (ctypes.c_double*40)()
  he3pak.SatXFunTdll(XpropL, XpropV, T)
  return XpropL[property - 1], XpropV[property - 1]
 
# returns thermodynamic property of saturated He3 liquid at specified temperature
def He3SaturatedLiquidProp(property, temperature):
  return He3SaturatedProp(property, temperature)[0]
  
# returns thermodynamic property of saturated He3 vapor at specified temperature
def He3SaturatedVaporProp(property, temperature):
  return He3SaturatedProp(property, temperature)[1]
  
# returns temperature of He3 at specified saturated vapor pressure
def He3SaturatedTemperature(pressure):
  P = ctypes.c_double(pressure)
  T = ctypes.c_double(0.)
  he3pak.TIPsatdll(T, P)
  return T.value