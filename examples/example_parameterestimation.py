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


case = 'example_03'

stream = dotter.build_model_from_config('cases/{}/config.ini'.format(case))

#measurements = book = xlrd.open_workbook('cases/{}/measurements.xlsx'.format(case))

stream.parameters['h'] = 6.65 - 5.13
stream.parameters['Q'] = 1.27
stream.generate_grid()
friction, res = dotter.estimate_roughness(stream, 6.694, 'waterlevel')

print (friction)
print (res)