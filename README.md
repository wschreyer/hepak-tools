# HEPAK Python modules

## HEPAK interface

hepak.py provides a Python interface to the [HEPAK Excel Add-In](http://www.htess.com/hepak.htm). It uses the Microsoft Excel COM API, so can only be used on Windows. It requires the pywin32 package for Python.

The module provides the same functions that are available in Excel when using the Add-In. See the HEPAK user guide and some basic information in the module itself.

### Thermodynamic properties available in HEPAK, see HEPAK user guide:

- 0 ('X'): Quality = vapor mass fraction
- 1 ('P'): Pressure
- 2 ('T'): Temperature
- 3 ('D'): Density
- 4 ('V'): Specific Volume
- 5: Compressibility factor Z = PV/RT
- 6: dP/dT (saturation line)
- 7: Latent Heat
- 8 ('S'): Entropy
- 9 ('H'): Enthalpy
- 10 ('A'): Helmholtz
- 11 ('U'): Internal Energy
- 12 ('G'): Gibbs Energy
- 13: Fugacity
- 14: Cp
- 15: Cv
- 16: Gamma = Cp/Cv
- 17: Expansivity = (T/V)(dV/dT)
- 18: Gruneisen parameter = (V/Cv)(dP/dT) at constant V
- 19: Isothermal compressibility = (1/D) dD/dP
- 20: Sound Velocity
- 21: Joule-Thomson Coefficient = dT/dP at constant H
- 22: dP/dD at constant T
- 23: dP/dT at constant D
- 24: V*dH/dV at constant P
- 25: Viscosity
- 26: Conductivity
- 27: Prandtl number
- 28: Thermal Diffusivity
- 29: Surface Tension
- 30: Dielectric constant - 1
- 31: Refractive index - 1
- 32: Isochoric dT to the lambda line dT(V)
- 33: Isobaric dT to the lambda line dT(P)
- 34: Superfluid density fraction RhoS/Rho
- 35: 2nd sound velocity
- 36: 4th sound velocity
- 37: Gorter-Mellink mutual friction constant
- 38: Superfluid thermal conductivity function
- 39: Lambda line temperature (isochoric to the state point) (T-Tlambda)

## HE3PAK interface

he3pak.py provides a Python interface to the [HE3PAK DLL library](http://www.htess.com/he3pak.htm). It works out-of-the-box with 32bit Python on Windows. It will NOT WORK with 64bit Python (unless you have a 64bit version of HE3PAK I guess?). It should in principle work on Linux as well with the zugbruecke package (not tested).

The module provides five functions calculating various properties of He3, see the module for all available properties:

- He3Density(pressure, temperature)
- He3Prop(property, density, temperature)
- He3SaturatedLiquidProp(property, temperature)
- He3SaturatedVaporProp(property, temperature)
- He3SaturatedTemperature(pressure)

### Thermodynamic properties available in HE3PAK:

- 1: Temperature
- 2: Pressure
- 3: Density
- 4: Specific volume (1/density)
- 5: Compressibility factor Z = PV/RT
- 6: Enthalpy
- 7: Entropy
- 8: C_v
- 9: C_p
- 10: Helmholtz free energy
- 11: Gibbs free energy
- 12: Sound velocity
- 13: Latent heat of evaporation
- 14: Joule-Thompson coefficient
- 15: inte
- 16: adiacomp
- 17: isocomp
- 18: volexp
- 19: isenexp
- 20: sndvirial
- 21: dBdT
- 22: trdvirial
- 23: dP/dD
- 24: dP/dT
- 25: dD/dT
- 26: dpdds
- 27: Gruneisen parameter
- 28: thcon
- 29: Viscosity
- 30: Viscosity/density
- 31: Surface tension
- 32: Viscosity*C_p/thcon
- 33: thcon/density/C_p
- 34: Dielectric constant
- 35: Refractive index
- 36: 2nd sound velocity
- 37: 4th sound velocity
- 38: Superfluid density fraction
- 39: Gorter Mellink mutual friction parameter
- 40: Superfluid thermal conductivity


## Installation on Windows

1. Download and install python (32bit!) for Windows.
2. Download the python modules, HEPAK Excel Add-In, and HE3PAK library all into the same folder.
3. Start a command line or PowerShell (shift+right-click in folder) and install pywin32 with `python -m pip install pywin32`.
4. If you've downloaded 64bit Python install msl-loadlib with `python -m pip install msl-loadlib`.
5. Run `python example.py` to execute the example script. You may have to acknowledge the dialog showing HEPAK License information (Starting Excel with the Add-In already installed beforehand will suppress this dialog). Then the script will print a few calculated He4 and He3 properties and list some information about the data available in HEPAK.


# Python scripts using the HEPAK modules

## Example script

example.py demonstrates the use of all the functions provided by both interfaces.

## UCN-source model

UCNsource.py provides a model of the He3 fridge for the new UCN source at TRIUMF.

UCNsource_parameterSweep.py does a scan of all experimental parameters and plots the resulting temperatures and flows in the fridge. Requires SciPy to solve the model equations and matplotlib to plot the results.

UCNsource_HeIIsegments.py calculates the temperature profile in the HeII channel and prints out a list of temperatures, UCN storage lifetimes, and imaginary Fermi potential averaged over segments along the channel.

### Supplemental data

The folder HEXdata contains He3 boiling curves from [Maeda, et al; Cryogenics 40 (2000) 713-719](https://doi.org/10.1016/S0011-2275(01)00002-9) and [Tanaka and Kodama; Cryogenics 29 (1989) 203-205](https://doi.org/10.1016/0011-2275(89)90085-4). This data is interpolated to model He3 boiling in the UCN fridge.

It also contains results of simulations of several heat exchangers in the fridge provided by [T. Okamura](https://kds.kek.jp/indico/event/31409/contributions/117125/attachments/91366/108628/report_TOkamura_KEK_20190624.pdf) for reference.

## Model of counter-flow heat exchanger HEX7

HEX7.py provides a model for a tube-in-tube counter-flow heat exchanger, similar to 'HEX7' in the fridge for the new UCN source at TRIUMF. It also includes an extension to a tube-in-tube-in-tube heat exchanger that can be used for a liquid-nitrogen-bath-cooled gas purifier. It cools the incoming flow to liquid nitrogen temperature using the return flow and the evaporated liquid nitrogen.

Requires SciPy to numerically solve the boundary-value problem and matplotlib to plot the results. Also uses [CoolProp](http://www.coolprop.org/) to calculate nitrogen-gas properties and, optionally, as an interface to helium-gas properties faster than the HEPAK Excel Add-In.