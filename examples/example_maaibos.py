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

# =============================================================================
# Parameters
# =============================================================================

interactive = False
RunCaseNumber = 2  # (python has zero-based indexing, so '0' means first case)

discharge_levels = [3., 4.0]

# =============================================================================
# Main code
# =============================================================================



stream = dotter.build_model_from_config('cases/example_01/config.ini')


fig, ax = plt.subplots(1)

# Run the model
for q in discharge_levels:
    stream.parameters['Q'] = q
    stream.parameters['n'] = 0.1
    stream.generate_grid()
    results, stream = dotter.run_model(stream)

    time = pd.to_datetime(results['waterdepth'].index.values)

    ax.plot(time, (results['waterlevel'].T.values[0]))

plt.show()
