# Dotter model prototype
#
# This file: provides high-level functions to run the model
# 
# Author: Koen Berends
# Contact: k.d.berends@utwente.nl or koen.berends@utwente.nl
# Copyright (c) 2016 University of Twente & Deltares
#

# region // imports
from configparser import ConfigParser
from .Classes import Stream
from . import Functions as F
import sys
import os
from datetime import datetime
# endregion

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
    results, stream = F.qs_solver(stream, restart=restart, starttime=starttime, stoptime=stoptime, write_to_file=write_output, verbose=verbose)
    return results, stream

def _load_config_file(input_file):
    [pathname, filename] = os.path.split(input_file)

    config = ConfigParser()
    config.read(input_file)
    files = {'events': os.path.join(pathname, config['files']['events']),
             'vegetation': os.path.join(pathname, config['files']['vegetation']),
             'geometry': os.path.join(pathname, config['files']['geometry']),
             'laterals': os.path.join(pathname, config['files']['laterals'])
             }

    # Parse homogeneous friction (=friction everywhere the same)
    if config['parameters']['use_homogeneous_friction'].strip().lower() == 'true':
        homofriction = True
    else:
        homofriction = False
    parameters = {
                  'x': [],
                  'z': [],
                  'Q': float(config['boundary']['Q']),             # Discharge per unit width
                  'ib': 0,           # Bed slope
                  'g': float(config['parameters']['g']),            # Gravitational acceleration
                  'n': [],        # Manning coefficient
                  'R': 0,               # Hydraulic radius
                  'A': 0,                # Hydraulic area
                  'h': float(config['boundary']['h']),        # downstream boundary
                  't': datetime.strptime(config['parameters']['t'], '%d/%m/%Y'), # current time
                  'dt': float(config['parameters']['dt']), # timestep (days)
                  'dx': float(config['parameters']['dx']), # spatial step (meters)
                  'friction_model': config['parameters']['friction_model'],
                  'use_homogeneous_friction': homofriction,
                  'constant_friction': float(config['parameters']['constant_friction'])
                  }

    return files, parameters, config

def plot(geometry, results=None, plot_moments=[], plot_locations=[], plottype=1, xlim=[]):
    if plottype == 1: 
        F.plot_1d_type01(geometry, results, plot_moments, plot_locations)
    if plottype == 2: 
        fig, ax = F.plot_1d_type02(geometry)
    if plottype == 3: 
        fig, ax = F.plot_1d_type03(geometry, xlim)
