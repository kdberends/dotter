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
from dotter import dotter
import matplotlib.pyplot as plt
import pandas as pd
from dotter import settings
import numpy as np
# =============================================================================
# Parameters
# =============================================================================

interactive = False
RunCaseNumber = 2  # (python has zero-based indexing, so '0' means first case)

discharge_levels = [3., 4., 10.]

# =============================================================================
# Main code
# =============================================================================



stream = dotter.build_model_from_config('cases/example_01/config.ini')


fig, ax = plt.subplots(1)

ax.set_title('MaaiBOS')
ax.set_xlabel('Dagen')
ax.set_ylabel('Verhang [m]')
# Gradient ()
normal_slope = 1.5
min_slope = 0
max_slope = 4
ns = 100

zs = np.linspace(min_slope, max_slope, ns)
ts = [0, 365]
vs = np.interp(zs, [min_slope, normal_slope, max_slope], [+1, 0, -1])

vs_matrix = np.array([vs] * 2)

ax.pcolor(ts, zs, vs_matrix.T, cmap="RdYlGn")

# Run the model
for q in discharge_levels:
    stream.parameters['Q'] = q
    stream.parameters['n'] = 0.1
    stream.generate_grid()
    results, stream = dotter.run_model(stream)

    time = pd.to_datetime(results['waterdepth'].index.values).values
    time = np.linspace(0, 365, len(time))
    waterlevel_up = results['waterlevel'].T.values[0]
    waterlevel_down = results['waterlevel'].T.values[-1]
    slope = waterlevel_up - waterlevel_down

    ax.plot(time, slope, color='k',
                         linestyle='--')

    ax.text(time[-1], slope[-1], '{:.1f} $m^3/s$'.format(q))

settings.two_axes_style(ax)
plt.show()
