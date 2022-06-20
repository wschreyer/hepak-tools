import numpy
import CoolProp.CoolProp
import scipy.integrate


hbar = 1.054571817e-34 # J s
boltzmannConstant = 1.380649e-23 # J/K
massNeutron = 1.67492749804e-27 # kg
massHe = 4.002602 *  1.66053906660e-27 # kg
HeScatteringLength = 3.26e-15 # m (https://www.ncnr.nist.gov/resources/n-lengths/elements/he.html)
HeBoundCrossSection = 1.34e-28 # m^2 (https://www.ncnr.nist.gov/resources/n-lengths/elements/he.html)
JpereV = 1.602176634e-19 # J/eV


# calculate He vapor pressure (Pa) from temperature (K)
def HeVaporPressure(T):
  # from Clement, Logan, Gaffney, Phys. Rev. 100, 743
  # https://doi.org/10.1103/PhysRev.100.743
  if not numpy.all(0.66 <= T <= 5.2):
    raise Exception('Tried to evaluate vapor pressure at T = {0}. Formula only valid between 0.66 and 5.2K!'.format(T))
  I = 4.6202
  A = 6.399
  B = 2.541
  C = 0.00612
  D = 0.5197
  a = 7.
  b = 14.14
  lnP = I - A/T + B*numpy.log(T) + C/2*T**2 - D*(a*b/(b**2 + 1) - 1./T)*numpy.arctan(a*T - b) - a*D/2/(b**2 + 1)*numpy.log(T**2/(1 + (a*T - b)**2))
  return 133.322387415*numpy.exp(lnP)


# calculate He density (kg/m3) at temperature T (K) and pressure P (Pa)
def HeDensity(T, P):
  return CoolProp.CoolProp.PropsSI('D', 'T', T, 'P', P, 'Helium')


# calculate real Fermi potential (neV) of He gas at temperature T (K) and pressure P (Pa)
def realFermiPotential(T, P):
  density = HeDensity(T, P)
  return 2. * numpy.pi * hbar**2 / massNeutron * density/massHe * HeScatteringLength / JpereV*1e9


# calculate imaginary Fermi potential (neV) of He gas at temperature T (K) and pressure P (Pa)
def imaginaryFermiPotential(T, P):
  density = HeDensity(T, P)
  
  freeCrossSection = HeBoundCrossSection/(1 + massNeutron/massHe)**2
  # upscattering cross section = "free" cross section * average He atom velocity / neutron velocity
  # see https://doi.org/10.1103/PhysRevC.92.065501
  averageHeVelocity = 2 * numpy.sqrt(2 * boltzmannConstant * T / numpy.pi / massHe)
  
  # imaginary Fermi potential according to Golub's book
  return 0.5 * hbar * density/massHe * freeCrossSection * averageHeVelocity / JpereV*1e9


def averageRealFermiPotential(Tmin, Tmax, P):
  integral = scipy.integrate.quad(lambda T: realFermiPotential(T, P), Tmin, Tmax)
  return integral[0] / numpy.abs(Tmax - Tmin)


def averageImaginaryFermiPotential(Tmin, Tmax, P):
  integral = scipy.integrate.quad(lambda T: imaginaryFermiPotential(T, P), Tmin, Tmax)
  return integral[0] / numpy.abs(Tmax - Tmin)


if __name__ == '__main__':
  T = 1.
  P = HeVaporPressure(T)
  Tsegments = [2.18, 15., 20., 30., 45., 60., 75., 90., 105., 120., 130., 135., 140., 145., 150., 155., 160., 165., 170., 175., 180., 185., 190., 195., 200., 205., 210., 215., 220., 230., 240., 250., 260., 290.]
  for Tmin, Tmax in zip(Tsegments[:-1], Tsegments[1:]):
    print(averageRealFermiPotential(Tmin, Tmax, P), averageImaginaryFermiPotential(Tmin, Tmax, P))