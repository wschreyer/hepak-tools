import UCNsource
import hepak
import matplotlib.pyplot as plt

parameters = {
'Beam heating':                   {'value': 8.1,         'range': (0., 10.),        'unit': 'W'}, # Beam on
#'Beam heating':                   {'value': 0.05,          'range': (0.05, 0.05),         'unit': 'W'}, # Beam off
'He-II static heat':   	          {'value': 0.25,        'range': (0., 1.),         'unit': 'W'},
'He-II funnel heat':              {'value': 1.,          'range': (0.5, 2.),        'unit': 'W'},
#'3He pumping speed':              {'value': 4300.,       'range': (300., 10000.),   'unit': r'm$^{3}$/h'}, # 1.14g/s @ 5.85 torr (Busch proposal)
#'He pumping speed':	             {'value': 2000.,       'range': (500., 5000.),    'unit': r'm$^{3}$/h'}, # 1.2g/s @ 9.9 torr (Busch proposal)
#'Pump inlet temperature':         {'value': 290.,      	'range': (100., 300.), 	   'unit': 'K'},
'3He pressure drop':              {'value': 100.,        'range': (0., 1000.),      'unit': 'Pa'},
#'3He pressure drop':              {'value': 20000.,      'range': (0., 50000.),     'unit': 'Pa'}, # Standby mode
'He pressure drop':               {'value': 100.,        'range': (0., 1000.),      'unit': 'Pa'},
#'He pressure drop':               {'value': 10000.,      'range': (0., 10000.),     'unit': 'Pa'}, # Standby mode
'He reservoir inlet temperature': {'value': 10.,         'range': (6., 12.,),       'unit': 'K'},
'He reservoir static load':       {'value': 0.6,         'range': (0., 2.),         'unit': 'W'},
'1K pot static load':             {'value': 0.05,        'range': (0., 1.),         'unit': 'W'},
'HEX4 exit temperature':          {'value': 2.8,         'range': (2.5, 4.),	      'unit': 'K'},
'HEX5 exit temperature':          {'value': 2.8,         'range': (2.5, 4.),        'unit': 'K'},
#'1K pot inlet temperature':       {'value': 2.8,         'range': (2., 3.5),        'unit': 'K'},
'3He inlet pressure':             {'value': 50000.,      'range': (30000., 80000.), 'unit': 'Pa'},
'Channel diameter':               {'value': 0.148,       'range': (0.12, 0.2),      'unit': 'm'},
'Channel length':                 {'value': 2.5,         'range': (1., 4.),         'unit': 'm'},
'HEX1 diameter':                  {'value': 0.148,       'range': (0.12, 0.2),      'unit': 'm'},
'HEX1 length': 	                  {'value': 0.6,         'range': (0.1, 0.8),       'unit': 'm'},
'HEX1 surface':                   {'value': 1.67,        'range': (0.21, 2.),       'unit': r'm$^{2}$/m'},
'HeII overfill':                  {'value': 0.05,        'range': (0.02, 0.2),      'unit': 'm'},
'He reservoir pressure':          {'value': 1.2e5,       'range': (700e2, 1300e2),  'unit': 'Pa'},
'Isopure He flow':                {'value': 0.,          'range': (0., 0.),         'unit': 'kg/s'},
#'Isopure He flow':                {'value': 0.00014/2,   'range': (0., 0.00014),     'unit': 'kg/s'}, # isopure condensation
'Isopure He pressure':            {'value': 100e2,       'range': (50e2, 500e2),    'unit': 'Pa'},
'20K shield temperature':         {'value': 20.,         'range': (15., 40.),       'unit': 'K'},
'100K shield temperature':        {'value': 100.,        'range': (80., 150.),      'unit': 'K'},
}

# scan parameter ranges and plot temperatures and flows
plotRows = 4
plotCols = 6
datapoints = 10
fig, axes = plt.subplots(plotRows, plotCols, figsize=(plotCols*6,plotRows*5))
axes2 = axes.copy()
axes3 = axes.copy()
fig.set_tight_layout(True)
for i, p in enumerate(parameters):
  tempParameters = {i: parameters[i]['value'] for i in parameters}
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
  y72 = []
  y8 = []
  y9 = []
  y9_2 = []
  for j in range(datapoints):
    value = parameters[p]['range'][0] + float(j)/(datapoints - 1)*(parameters[p]['range'][1] - parameters[p]['range'][0])
    tempParameters[p] = value
    res = UCNsource.calcUCNSource(tempParameters)
    if res:
      print('{0:.3g} {1}: {2:.1f} L/h, T_low {3:.3f}'.format(value, 
                                                             parameters[p]['unit'], 
                                                             (max(res['He reservoir flow'], res['20K shield flow'], res['100K shield flow']) + res['1K pot flow'])/hepak.HeCalc('D', 0, 'P', 1013e2, 'SL', 0., 1)*1000*3600, 
                                                             res['T_HeII_low'])
                                                            )
      x.append(value)
      y.append(res['He reservoir flow']*1000)
      y2.append(res['T_HeII_high'][0])
      y2_2.append(res['T_HeII_high'][1])
      y3.append(res['T_3He'])
      y4.append(res['T_HeII_low'])
      y5.append(res['3He flow']*1000)
      y6.append(res['1K pot flow']*1000)
      y7.append(res['20K shield flow']*1000)
      y72.append(res['100K shield flow'] * 1000)
      y8.append(res['1K pot temperature'])
      y9.append(res['T_HEX1_low'])
      y9_2.append(res['T_HEX1_high'])
  ax = axes[i % plotRows][int(i/plotRows)]
  axes[0][0].get_shared_y_axes().join(axes[0][0], ax)
  ax.plot(x, y, color = 'tab:green', label = 'He reservoir')
  ax.plot(x, y6, color = 'tab:olive', label = '1K pot')
  ax.plot(x, y7, color = 'tab:cyan', label = '20K shields')
  ax.plot(x, y72, color = 'darkcyan', label = '100K shields')
  ax.plot(x, y5, color = 'tab:orange', label = '3He')
  ax.set_xlabel(p + ' (' + parameters[p]['unit'] + ')')
  ax.set_ylabel('Gas flow (g/s)')
  ax.minorticks_on()
  ax.grid(True, 'both', 'both', alpha = 0.2)

  axes2[i % plotRows][int(i/plotRows)] = ax.twinx()
  ax2 = axes2[i % plotRows][int(i/plotRows)]
  axes2[0][0].get_shared_y_axes().join(axes2[0][0], ax2)
  ax2.plot(x, y8, color = 'tab:pink', label = '1K pot')
  ax2.plot(x, y3, color = 'tab:red', label = '3He')
  ax2.fill_between(x, y9, y9_2, color = 'tab:gray', label = 'HEX1', alpha = 0.3)
  ax2.plot(x, y4, color = 'tab:purple', label = '4He at HEX1')
  ax2.fill_between(x, y2, y2_2, color = 'tab:blue', label = '4He in bottle', alpha = 0.3)
  ax2.set_ylabel('Temperature (K)')
  ax2.minorticks_on()
  
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