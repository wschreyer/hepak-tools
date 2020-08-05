import csv
import numpy
import scipy.interpolate
import matplotlib.pyplot as plt
import matplotlib.colors

# interpolate measured He3 boiling curves
# from Maeda, Beppu, Fujii, Shigi, Cryogenics 40 (2000) 713-719)
# and Tanaka, Kodama, Cryogenics 29 (1989) 203
def loadHe3boilingData():
  TQdata = []
  dTdata = []
  for T, filename in zip([0.55, 0.6, 0.7, 0.8, 1.0, 1.48, 2.], ['0.55K.csv', '0.6K.csv', '0.7K.csv', '0.8K.csv', '1K.csv', '1.48K.csv', '2K.csv']):
    with open('HEXdata/He3boiling_' + filename) as csvfile:
      csvreader = csv.reader(csvfile, delimiter = ',')
      for line in csvreader:
        TQdata.append([T, float(line[1])*10000.]) # data for heat flux is in W/cm2, convert to W/m2		
        dTdata.append(float(line[0]))
  spline = scipy.interpolate.LinearNDInterpolator(TQdata, dTdata, rescale = True)
  nearest = scipy.interpolate.NearestNDInterpolator(TQdata, dTdata, rescale = True)
  T, q = numpy.meshgrid(numpy.linspace(0.5, 2.1), numpy.logspace(-2, 4))
  h = numpy.divide(q, numpy.vectorize(spline)(T, q))
  fig, ax = plt.subplots(1, 1)
  im = ax.pcolormesh(T, q, h, cmap = 'hsv', norm = matplotlib.colors.LogNorm())
  ax.set_yscale('log')
  ax.set_xlabel('Temperature (K)')
  ax.set_ylabel(r'Heat flux (W/m$^{2}$)')
  zax = fig.colorbar(im, ax = ax)
  zax.set_label(r'Heat transfer coefficient (W/m$^{2}$ K)')
  plt.savefig('He3boiling.pdf')
  return spline, nearest
