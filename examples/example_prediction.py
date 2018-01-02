# Dotter model prototype
#
# This file: high-level run file
# 
# Author: Koen Berends
# Contact: k.d.berends@utwente.nl & koen.berends@utwente.nl
# Copyright (c) 2017 University of Twente & Deltares
#
# Permission is hereby granted, free of charge, to any person obtaining a 
# copy of this software and associated documentation files (the "Software"), 
# to deal in the Software without restriction, including without limitation 
# the rights to use, copy, modify, merge, publish, distribute, sublicense, 
# and/or sell copies of the Software, and to permit persons to whom the Software 
# is furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all 
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, 
# INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR 
# PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE 
# FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR 
# OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
# DEALINGS IN THE SOFTWARE.

# =============================================================================
# Imports & Function definitions
# =============================================================================
from Dotter import dotter 
from Dotter import Functions as F
import matplotlib.pyplot as plt
import os
import shutil
import sys 

# =============================================================================
# Parameters
# =============================================================================

interactive = False
RunCaseNumber = 1  # (python has zero-based indexing, so '0' means first case)

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
    results, stream = dotter.run_model(stream)

    # Move results to case_folder
    try:
        os.mkdir('cases/{}/output'.format(case))
    except OSError:
        # Dir already exists
        pass

    for filename in ['waterdepth', 'roughness', 'waterlevel', 'max_allowed_waterlevel', 'percentagebegroeiing', 'discharge']:
        shutil.move('{}.csv'.format(filename), 'cases/{case}/output/{fname}.csv'.format(case=case, fname=filename ))

    # =============================================================================
    # Plot results
    # =============================================================================
    
    prefs = F.load_plotting_config('cases/{case}/plotting.ini'.format(case=case))
    if prefs['1d']['make_plot']:
        plot_moments = prefs['1d']['plot_moments']
        plot_locations = prefs['1d']['plot_locations']
        dotter.plot(stream, results, plot_moments, plot_locations, plottype=int(prefs['1d']['type']))

    if prefs['2d']['make_plot']:
        data = F.read_results_file('cases/{}/output/waterdepth.csv'.format(case))
        fig, ax = F.plot_2d(data, colormap=prefs['2d']['colormap'], threshold=None)
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
