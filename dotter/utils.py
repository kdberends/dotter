"""
Utilities, convenience functions
"""

# =============================================================================
# Imports
# =============================================================================

import os
import sys
import logging
import configparser
import seaborn as sns
import calendar
from datetime import datetime
if (sys.version_info > (3, 0)):
    from io import StringIO
else:
    from cStringIO import StringIO

flatui = ["#9b59b6", "#3498db", "#95a5a6", "#e74c3c", "#34495e", "#2ecc71"]
flatui_r = list(reversed(flatui))

# Configuration file parsing types
configtypes = """[parameters]
geometrytype=int
g=float
tstart=datetime
tstop=datetime
dt=int
dx=int
frictionmodel=str
growthmodel=str
blockagemodel=str
"""

eventypes= """[event]
eventtype=str
tstart=datetime
minchainage=float
maxchainage=float
reduce_to=float
maximum_blockage=float
triggered=bool
name=str
"""

datetimeformat = '%d/%m/%Y'
# =============================================================================
# Definitions
# =============================================================================

def parse_dict(input_dict, typedict=configtypes):
    """
    In: configparser object, out: dict with floats/int/booleans
    """
    P = dict()
    ptypes = StringIO(typedict)
    config = configparser.ConfigParser()
    config.readfp(ptypes)
    for section in config:
        P[section] = dict()
        for key in config[section]:
            if config[section][key] == 'bool':
                P[section][key] = input_dict[section].getboolean(key)
            elif config[section][key] == 'int':
                P[section][key] = input_dict[section].getint(key)
            elif config[section][key] == 'float':
                P[section][key] = input_dict[section].getfloat(key)
            elif config[section][key] == 'datetime':
                P[section][key] = datetime.strptime(input_dict[section][key], datetimeformat),
                P[section][key] = P[section][key][0]
            else:
                P[section][key] = input_dict[section][key]

    return P

def get_logger(outputpath=os.getcwd(), logfile='dotter.log', overwrite=False):
    """
    Returns a logger object which:
    -
    """
    filepath = os.path.join(outputpath, logfile)
    if overwrite and os.path.isfile(filepath):
        mode = 'w'
    else:
        mode = 'a'

    logger = logging.getLogger(__name__)
    if not logger.handlers:
        logger.setLevel(logging.DEBUG)

        # the filehandler prints every level to file
        filehandler = logging.FileHandler(filepath, mode=mode)
        filehandler.setLevel(logging.DEBUG)

        # the stream prints info and above (error, warning, critical) to console
        streamhandler = logging.StreamHandler()
        streamhandler.setLevel(logging.DEBUG)

        # create a logging format
        formatter = logging.Formatter('%(asctime)s - %(filename)s - %(levelname)s - %(message)s',
                                      "%Y-%m-%d %H:%M:%S")

        filehandler.setFormatter(formatter)
        streamhandler.setFormatter(formatter)

        # add the handlers to the logger
        logger.addHandler(filehandler)
        logger.addHandler(streamhandler)

        logger.info('Start logging to {}'.format(filepath))

    return logger

def set_plotstyle(style='', palette=flatui_r, scale=1.2):
    """
    Arguments:

    style:
          'paper': for use in journal
          'digital': for use on-screen

    palette:
          name of colorpalette, e.g.:
          - flatui
          - paired
          - deep, muted, pastel, bright, dark, colorblind
    """
    sns.set_context("paper")

    if style.lower() == 'paper':
        # Set the font to be serif, rather than sans
        sns.set(font='serif',
                palette=palette,
                font_scale=2,
                style='whitegrid')
    elif style.lower() == 'digital':
        sns.set(font='serif',
                palette=palette,
                font_scale=1.2,
                style='whitegrid')
    else:
        sns.set(font='serif',
                palette=palette,
                font_scale=scale,
                style='whitegrid')

def two_axes_style(ax):
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_position(('outward', 10))
    ax.spines['left'].set_linewidth(2)
    ax.spines['left'].set_edgecolor('k')
    ax.spines['bottom'].set_position(('outward', 10))
    ax.spines['bottom'].set_edgecolor('k')
    ax.spines['bottom'].set_linewidth(2)
    ax.grid(False)

def DatetimeToTimestamp(dates):
    for d in dates:
        yield calendar.timegm(d.timetuple())

def TimestampToDatetime(dates):
    for d in dates:
        yield datetime.fromtimestamp(d)
