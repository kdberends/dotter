# Dotter model prototype
#
# This file: provides high-level functions to run the model
# 
# Author: Koen Berends
# Contact: k.d.berends@utwente.nl or koen.berends@utwente.nl
# Copyright (c) 2016 University of Twente & Deltares
#

# =============================================================================
# Imports
# =============================================================================

from configparser import ConfigParser
from .Classes import Stream
from . import Functions as F
import sys
import os
import numpy as np
from datetime import datetime
from scipy import optimize

# =============================================================================
# I Function definitions
# =============================================================================


def build_model_from_config(filename, verbose=True):
    """
    Returns 
    """

    (files, parameters, config) = _load_config_file(filename)

    # initialise geometry
    if verbose: sys.stdout.write('Initialising geometry...\n')
    stream = Stream(**parameters)
    stream.load_from_excel(files['geometry'])

    if verbose: sys.stdout.write('Loading lateral sources...\n')
    stream.load_laterals(files['laterals'])

    if verbose: sys.stdout.write('Generating computational grid...\n')
    stream.generate_grid()  

    if verbose: sys.stdout.write('Planting vegetation...\n')
    stream.set_vegetation(files['vegetation'])

    if verbose: sys.stdout.write('Loading events...\n')
    stream.load_events(files['events'])

    if verbose: sys.stdout.write('Model loaded\n')
    return stream


def run_model(stream, restart=None, starttime=0, stoptime=365, write_output=True, verbose=True):
    results, stream = F.qs_solver(stream,
                                  restart=restart,
                                  starttime=starttime,
                                  stoptime=stoptime,
                                  write_to_file=write_output,
                                  verbose=verbose)
    return results, stream


def plot(geometry, results=None, plot_moments=[], plot_locations=[], plottype=1, xlim=[]):
    if plottype == 1: 
        F.plot_1d_type01(geometry, results, plot_moments, plot_locations)
    if plottype == 2: 
        fig, ax = F.plot_1d_type02(geometry)
    if plottype == 3: 
        fig, ax = F.plot_1d_type03(geometry, xlim)


def estimate_roughness(stream, measurement, variable):
    """
    args:
        measurement: float
        variable: e.g. 'water depth'
    """

    res = optimize.minimize(_gof, 0.03, args=[stream, measurement, variable],
                                              tol=0.001,
                                              method="Nelder-Mead",
                                              options=dict(maxiter=50))

    return res.x, res


def _gof(n, *args):
    """
    goodness of fit helper function for optimization
    """
    n = n[0]
    stream = args[0][0]
    measurement = args[0][1]
    variable = args[0][2]

    # Set roughness and run model
    stream.parameters['constant_friction'] = n 

    results, stream = run_model(stream, stoptime=1, verbose=False)

    #sys.stdout.write("n: {}\n measurement: {}\n model: {}\n".format(n, measurement, results[variable].values[0][0]))
    gof = np.abs(results[variable].values[0][0] - measurement)
    return gof


def _load_config_file(input_file):
    [pathname, filename] = os.path.split(input_file)

    config = ConfigParser()
    config.read(input_file)

    files = dict()
    for key in config['files']:
        files[key] = os.path.join(pathname, config['files'][key])

    # Parse constant friction (=friction everywhere the same)
    if config['parameters']['use_constant_friction'].strip().lower() == 'true':
        homofriction = True
    else:
        homofriction = False
    parameters = {
                  'path': pathname,
                  'Q': float(config['boundary']['Q']),             # Discharge per unit width
                  'g': float(config['parameters']['g']),            # Gravitational acceleration
                  'h': float(config['boundary']['h']),        # downstream boundary
                  't': datetime.strptime(config['parameters']['t'], '%d/%m/%Y'), # current time
                  'dt': float(config['parameters']['dt']), # timestep (days)
                  'dx': float(config['parameters']['dx']), # spatial step (meters)
                  'friction_model': config['parameters']['friction_model'],
                  'use_constant_friction': homofriction,
                  'constant_friction': float(config['parameters']['constant_friction'])
                  }

    return files, parameters, config