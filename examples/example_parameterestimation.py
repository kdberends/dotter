# Dotter model prototype
#
# This file: high-level run file
#
# Author: Koen Berends
# Contact: k.d.berends@utwente.nl & koen.berends@utwente.nl
# Copyright (c) 2017 Deltares
# Copyright (c) 2017 University of Twente


# =============================================================================
# Imports & Function definitions
# =============================================================================
from dotter.models import DotterModel
import matplotlib.pyplot as plt
# =============================================================================
# Parameters
# =============================================================================

deltabeek = DotterModel('cases/vegetation/config.ini')


fig, ax = plt.subplots(1)
#ax.plot(deltabeek.output.blockage.iloc[:, 0])
#ax.pcolor(deltabeek.grid.friction)
ax.pcolor(deltabeek.grid.time, deltabeek.grid.chainage, deltabeek.output.blockage.T)
plt.show()

#deltabeek.run()
#deltabeek.dash(dashtype=2, show=True)
