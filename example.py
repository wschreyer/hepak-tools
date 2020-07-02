import hepak
import he3pak

print('He4 Standard density:', hepak.HeCalc('D', 0, 'P', 101300, 'T', 273.15, 1))
print(hepak.HeValidate(1, 2), hepak.HeMsg(hepak.HeCalc(-1, 0, 'P', 101300, 'T', 273.15, 1)))

print('\nList of HEPAK properties and their units:')
for p in range(40):
  print(p, hepak.HeProperty(p), [hepak.HeUnit(p, u) for u in range(1, 5)])
  
print('\nList of HEPAK constants:')
for p in range(1, 16):
  print(p, hepak.HeConst(p))
print('\n')

print('He3 standard density:', he3pak.He3Density(101300, 273.15), 'kg/m3')
print('He3 saturated density at 3K: {0} kg/m3 (liquid), {1} kg/m3 (vapor)'.format(he3pak.He3SaturatedLiquidProp(3, 0.8), he3pak.He3SaturatedVaporProp(3, 0.8)))
print('He3 temperature at saturated vapor pressure of 1 atm:', he3pak.He3SaturatedTemperature(101300), 'K')
