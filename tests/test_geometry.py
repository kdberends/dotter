#!/usr/bin/env python
# -*- coding: utf-8 -*-

# =============================================================================
# Imports & Function definitions
# =============================================================================
from dotter.models import DotterModel
import numpy as np 
import matplotlib.pyplot as plt

# =============================================================================
# Parameters
# =============================================================================

# Equilibrium
#depth = 1.5681
#h = deltabeek.grid.bedlevel[0] + depth
#A = deltabeek.grid.wet_area[0](h)
#R = deltabeek.grid.hydraulic_radius[0](h)
#i = deltabeek.grid.bedslope[0]
#C = R ** (1 / 6.) / 0.04
#Q = A * C * np.sqrt(R * i)
#print (Q)
#ax = fig.add_subplot(212)
#ax.plot(deltabeek.grid.chainage, deltabeek.output.waterlevel.iloc[0] - deltabeek.grid.bedlevel)
#ax.plot(deltabeek.grid.chainage, deltabeek.grid.chainage ** 0 + depth -1 , '--r')
#




#print (deltabeek.grid.wet_area[-1](3.298))
#ax.plot(deltabeek.grid.chainage, deltabeek.output.waterlevel.iloc[0] - deltabeek.grid.bedlevel)
#ax.plot(deltabeek.grid.chainage, deltabeek.output.waterlevel.iloc[0])
#ax.plot(deltabeek.grid.chainage, deltabeek.grid.bedlevel)



"""

from dotter import dotter

case = 'example_02'

stream = dotter.build_model_from_config('cases/{}/config.ini'.format(case))
results, stream = dotter.run_model(stream, stoptime=1)

ax.plot(results['x'].iloc[0], results['waterdepth'].iloc[0], '-r')
"""


plt.show()

