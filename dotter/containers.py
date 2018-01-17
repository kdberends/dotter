"""
Custom container classes to store data
"""

# =============================================================================
# Imports
# =============================================================================
from configparser import ConfigParser
from datetime import datetime, timedelta
from collections import namedtuple
import pandas as pd
import numpy as np
from scipy.interpolate import UnivariateSpline, griddata, interp1d
from . import utils
import matplotlib.pyplot as plt
# =============================================================================
# Classes
# =============================================================================

ParameterContainer = namedtuple('Parameters', ['g', 'dx', 'dt',
                                               'tstart', 'tstop',
                                               'geometrytype',
                                               'frictionmodel', 'growthmodel',
                                               'blockagemodel'])

FileContainer = namedtuple('Files', ['geometry', 'vegetation',
                                     'events', 'laterals', 'measurements'])

SamplesContainer = namedtuple('Samples', ['X', 'Y', 'Z',
                                          'friction'])

ResultsContainer = namedtuple('Results', ['waterlevel', 'waterdepth', 'friction', 'blockage'])

BoundaryContainer = namedtuple('Boundary', ['type', 'q', 'h'])

VegetationContainer = namedtuple('VegetationContainer', ['name', 'tstart', 'tstop',
                                                         'growthspeed', 'decayspeed',
                                                         'maximumcover', 'density',
                                                         'stemheight', 'stemdiameter'])

EventContainer = namedtuple('Event', ['eventtype', 'tstart', 'minchainage',
                                              'maxchainage', 'reduce_to', 'maximum_blockage',
                                              'triggered', 'name'])

for container in [ParameterContainer, FileContainer, SamplesContainer, ResultsContainer]:
    container.__new__.__defaults__ = (None,) * len(container._fields)


class GeometryGrid:
    """

    """
    def __init__(self, parameters, samples, logger=None):
        """
        inputkwargs: chainage, bedlevel, discharge, friction, hydraulics, blockage
        """
        if not logger:
            self.logger = utils.get_logger()
        else:
            self.logger = logger

        self.parameters = parameters
        # Samples
        self.samples = SamplesContainer(**samples)

        # Grid variables
        self.time = self.get_timevector()
        self.chainage = self.get_chainagevector()
        self.h_resolution = 50
        self.max_depth = 5
        self.X = None
        self.Y = None
        self.Z = None
        self.bedlevel = None
        self.bedslope = None
        #self.discharge = pd.DataFrame(index=time, columns=chainage)
        #self.hydraulics = pd.DataFrame(index=time, columns=chainage)
        self.generate_grid()

    def generate_grid(self):
        """
        generates a grid
        """

        # Two-dimensional (X, Y)
        ygrid = self.samples.Y[0]
        xgrid = self.chainage

        self.X = np.array([self.chainage] * len(ygrid))
        self.Y = np.array([ygrid] * len(xgrid)).T
        self.Z = griddata(np.array([self.samples.X.flatten(), self.samples.Y.flatten()]).T,
                          self.samples.Z.flatten(), (self.X, self.Y), method='linear')


        # one-dimensional (time invariant)
        self.vegetation = None
        self.bedlevel = self.Z.min(axis=0)
        self.bedslope = np.abs(np.diff(self.bedlevel) / np.diff(self.chainage))
        self.bedslope = np.round(np.append([self.bedslope[0]], self.bedslope), decimals=10)
        #self.friction = [np.interp(self.chainage, self.samples.X.T[0], self.samples.friction[0])]

        (self.wet_area, self.hydraulic_radius) = self.__calculate_hydraulics()

    def get_chainagevector(self):
        """
        if there is no integer number of steps between minimum and maximum
        chainage, the dx is changed to adapt.
        """
        maxx = np.max(self.samples.X)
        minx = np.min(self.samples.X)
        self.nx = int((maxx - minx) / self.parameters.dx)
        self.dx = (maxx - minx) / self.nx
        self.logger.debug('using spatial step of {} m'.format(self.dx))

        return np.linspace(minx, maxx, self.nx + 1)

    def get_timevector(self):
        """
        timestep is dominant. end time is changed when necessary
        """
        t = self.parameters.tstart
        time = [t]
        nt = (self.parameters.tstop - self.parameters.tstart) / self.parameters.dt
        self.logger.debug('number of timesteps: {}'.format(nt.days))

        for i in range(nt.days):
            t += timedelta(days=self.parameters.dt)
            time.append(t)

        self.logger.debug('start time: {}'.format(time[0]))
        self.logger.debug('top time: {}'.format(time[-1]))

        return np.array(time)

    def set_discharge(self, discharge):
        """
        sets discharge, taking into account laterals.

        Always call after set_boundaries
        """

        self.discharge[:] = discharge
        if self.laterals is not None:
            for lateralchainage, factor in zip(self.laterals.x, self.laterals.factor):
                ind = np.where(self.chainage > lateralchainage)[0][0]
                indlen = len(self.chainage) - ind
                self.discharge.iloc[:, ind:] = np.tile(self.discharge.iloc[:, 0] * factor, (indlen, 1)).T

    def set_boundaries(self, data, laterals=None):
        """
        Interpolates boundary conditions to the grid
        """
        datatime = list(utils.DatetimeToTimestamp(pd.to_datetime(data['date'])))
        gridtime = list(utils.DatetimeToTimestamp(self.time))

        # Discharge
        # ---------------------------------------------------------------------
        Q = np.interp(gridtime, datatime, data['discharge'])
        self.discharge = pd.DataFrame(index=self.time,
                                      columns=self.chainage,
                                      data=np.array([Q] * len(self.chainage)).T,
                                      dtype='float')

        self.laterals = laterals
        if laterals is not None:
            for lateralchainage, lateralname in zip(laterals.chainage, laterals.name):
                ind = np.where(self.chainage > lateralchainage)[0][0]
                indlen = len(self.chainage) - ind
                self.logger.debug('Lateral {} inserted at chainage {}'.format(lateralname, lateralchainage))

                # interpolate and add lateral discharge to grid
                qlat = np.interp(gridtime, datatime, data[lateralname])
                self.discharge.iloc[:, ind:] = np.tile(self.discharge.iloc[:, 0] + qlat, (indlen, 1)).T
        # Downstream
        # ---------------------------------------------------------------------
        H = np.interp(gridtime, datatime, data['downstream'])
        self.downstream = pd.Series(index=self.time,
                                    data=H,
                                    dtype='float')

        # Downstream
        # ---------------------------------------------------------------------
        H = np.interp(gridtime, datatime, data['upstream'])
        self.upstream = pd.Series(index=self.time,
                                  data=H,
                                  dtype='float')

    # =============================================================================
    # private methods
    # =============================================================================

    def __calculate_hydraulics(self):
        """ calculate hydraulic radius, area from arbitrary cross-section"""

        hydraulic_radius = []
        wet_area = []
        # Loop over cross-sections
        for i, ix in enumerate(self.chainage):
            A = list()  # Wet area
            P = list()  # Wet perimeter
            waterlevels = np.linspace(np.min(self.Z.T[i]), np.min(self.Z.T[i]) + self.max_depth, self.h_resolution)
            for j, h in enumerate(waterlevels):
                if h > np.min(self.Z.T[i]):
                    # add new point in cross-section
                    (y, z) = self.__hdep(self.Y.T[i], self.Z.T[i], h)

                    yy = np.append(y, list(reversed(y)))
                    zz = np.append(z, np.ones(len(z)) * h)
                    if (i==999) and (j == 20):
                        fig, ax = plt.subplots(1)
                        ax.plot(yy, zz, '.-r')
                        ax.plot(y, np.ones(len(z)) * h, '.-b')
                        plt.show()
                    A.append(self.__polyarea(yy, zz))
                    P.append(np.sum(np.sqrt(np.diff(y)**2 + np.diff(z)**2)))
                else:
                    A.append(0)
                    P.append(1)
            A = np.array(A)
            P = np.array(P)
            R = A / P


            # Construct univariate splines for fast lookup
            wet_area.append(interp1d(waterlevels, A, bounds_error=False, fill_value=0))
            hydraulic_radius.append(interp1d(waterlevels, R, bounds_error=False, fill_value=0))

        return (wet_area, hydraulic_radius)

    def __hdep(self, y, z, h):
        """ insert point at h """
        zmask = np.where(z < h)[0]

        if zmask[0] > 0:
            yleft = np.interp(h, [z[zmask[0]], z[zmask[0] - 1]],
                                 [y[zmask[0]], y[zmask[0] - 1]])
        else:
            yleft = y[zmask[0]]
        if zmask[-1] < len(z) - 1:
            yright = np.interp(h, [z[zmask[-1]], z[zmask[-1] + 1]],
                                  [y[zmask[-1]], y[zmask[-1] + 1]])

        else:
            yright = y[zmask[-1]]
        # add point left
        z = np.concatenate([[h], z[zmask], [h]])
        y = np.concatenate([[yleft], y[zmask], [yright]])


        return (y, z)

    @staticmethod
    def __polyarea(x, y):
        return 0.5 * np.abs(np.dot(x, np.roll(y, 1)) - np.dot(y, np.roll(x, 1)))
