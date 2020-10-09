import csv
import numpy
import scipy.interpolate
import matplotlib.pyplot as plt
import matplotlib.ticker
import hepak
import he3pak

def plotPumpData(Pdata, flowData):
  fig, ax = plt.subplots(constrained_layout = True)
  ax.plot(Pdata, numpy.multiply(flowData, 1000.*2.), label = '3He')
  ax.plot(Pdata, numpy.multiply(flowData, 1000./0.75), label = '4He')
  ax.grid(which = 'both')
  ax.set_xscale('log')
  #ax.set_yscale('log')
  ax.set_xlabel('Pressure (Pa)')
  ax.set_ylabel('Mass flow (g/s)')
  ax.legend(loc = 'upper left')
  ax.set_xlim(100, 10000)
  ax.set_ylim(0, 2)
  ax2 = ax.secondary_xaxis('top', functions = (lambda P: numpy.vectorize(he3pak.He3SaturatedTemperature)(P + 100.), \
                                               lambda T: numpy.vectorize(lambda x: he3pak.He3SaturatedLiquidProp(2, x))(T) - 100.))
  ax2.set_xlabel('Equivalent 3He temperature (K), with 100 Pa pressure drop added')
  ax2.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))
  ax2.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1f'))
  ax3 = ax.secondary_xaxis(1.15, functions = (lambda P: numpy.vectorize(hepak.HeSaturatedTemperature)(P + 100.), \
                                              lambda T: numpy.vectorize(hepak.HeSaturatedPressure)(T) - 100.))
  ax3.set_xlabel('Equivalent 1K pot temperature (K), with 100 Pa pressure drop added')
  ax3.xaxis.set_major_locator(matplotlib.ticker.MultipleLocator(0.1))
  ax3.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('%.1f'))
  plt.savefig('pumpingCurve.pdf')


# load pumping curves of Busch 2-stage system
def loadPumpData(filename):
  Pdata = []
  flowData = []
  with open(filename) as csvfile:
    csvreader = csv.reader(csvfile, delimiter = ',')
    for line in csvreader:
      Pdata.append(float(line[0])*133.322) # pressure data is in torr
      flowData.append(float(line[1])/1000.) # flow data is in g/s
#  plotPumpData(Pdata, flowData)
  return scipy.interpolate.interp1d(Pdata, flowData)
  
