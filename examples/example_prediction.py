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
from dotter import Functions as F
import matplotlib.pyplot as plt
import os
import shutil
import sys

# =============================================================================
# Parameters
# =============================================================================

interactive = False
RunCaseNumber = 2  # (python has zero-based indexing, so '0' means first case)

# =============================================================================
# Main code
# =============================================================================


def main(case):

    # Set plotting style
    F.set_style()

    # =============================================================================
    # Load and solve the model
    # =============================================================================

    stream = dotter.build_model_from_config('cases/{}/config.ini'.format(case))

    # Run the model
    stream.parameters['Q'] = 1.27
    stream.parameters['n'] = 0.1
    stream.generate_grid()
    results, stream = dotter.run_model(stream)


    # =============================================================================
    # Plot results
    # =============================================================================

    prefs = F.load_plotting_config(
        'cases/{case}/plotting.ini'.format(case=case))
    if prefs['1d']['make_plot']:
        plot_moments = prefs['1d']['plot_moments']
        plot_locations = prefs['1d']['plot_locations']
        dotter.plot(stream, results, plot_moments, plot_locations,
                    plottype=int(prefs['1d']['type']))

    if prefs['2d']['make_plot']:
        data = F.read_results_file(
            'cases/{}/output/waterdepth.csv'.format(case))
        fig, ax = F.plot_2d(
            data, colormap=prefs['2d']['colormap'], threshold=None)
        ax.set_title('Waterdiepte')

    if prefs['1d']['make_plot'] or prefs['2d']['make_plot']:
        plt.show()


if __name__ == '__main__':
    sys.stdout.write('=========================\n')
    sys.stdout.write('== Dotter model (2017) ==\n')
    sys.stdout.write('=========================\n\n')
    sys.stdout.write('Available cases: \n')

    cases = os.listdir('cases')
    runcase = cases[RunCaseNumber]
    sys.exit(main(runcase))

# ~end script
