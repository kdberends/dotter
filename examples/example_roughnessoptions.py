# Dotter model prototype
#
#
# Author: Koen Berends
# Contact: k.d.berends@utwente.nl & koen.berends@utwente.nl
# Copyright (c) 2017 Deltares
# Copyright (c) 2017 University of Twente


# =============================================================================
# Imports & Function definitions
# =============================================================================
from dotter import dotter
from dotter import notebook

case = 'example_03'

stream = dotter.build_model_from_config('cases/{}/config.ini'.format(case))

measurements = np.loadtxt('cases/{}/measurements.csv'.format(case), skiprows=1, delimiter=',')

Q = measurements.T[0]
h_up = measurements.T[1]
h_down = measurements.T[2]

frictionvalues = []
for q, hup, hdown in zip(Q, h_down, h_up):
    stream.parameters['h'] = hdown - 5.13
    stream.parameters['Q'] = q
    stream.generate_grid()
    friction, res = dotter.estimate_roughness(stream, hup, 'waterlevel')
    print ("Q: {}, hup: {}, n:{}\n".format(q, hup, friction[0]))
    frictionvalues.append(friction[0])

np.savetxt('frictionvalues.csv', frictionvalues, delimiter=',')

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(list(range(len(frictionvalues))), frictionvalues, '.')

plt.show()