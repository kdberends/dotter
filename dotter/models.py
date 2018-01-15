"""
Dottermodel main class
"""

# =============================================================================
# Imports
# =============================================================================

from . import utils
from . import containers
import os
from tqdm import tqdm
import numpy as np
from configparser import ConfigParser
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import timedelta
# ===================atom-minimap==========================================================
# Classes
# =============================================================================

class DotterModel:

    def __init__(self, configfilepath):

        self.__set_outputpath_from_config(configfilepath)
        self.logger = utils.get_logger(self.outputpath, overwrite=True)
        self.maxwidth = 20   # to generate trapezoidal bathy
        self.vegetationdefs = list()
        self.blockagemodel = self.pglinear
        self.dateformat = '%d/%m/%Y'
        self.frictionmodel = self.__manning
        self.load(configfilepath)

    def load(self, configfilepath):
        """
        loads model from configfile
        """
        self.logger.info('Loading configuration file {}'.format(configfilepath))
        self.__parseconfigfile(configfilepath)

        # initialise geometry
        # ---------------------------------------------------------------------
        self.logger.info('Loading geometry from {}'.format(self.files.geometry))
        self.__parsegeometry()

        # initialize results
        # ---------------------------------------------------------------------
        rmat = pd.DataFrame(index=self.grid.time, columns=self.grid.chainage, dtype='float')
        self.output = containers.ResultsContainer(waterlevel=rmat.copy(),
                                                  waterdepth=rmat.copy(),
                                                  friction=rmat.copy(),
                                                  blockage=rmat.copy())

        # initialise boundary conditions
        # ---------------------------------------------------------------------
        self.logger.info('Loading boundary conditions')
        self.__parsemeasurements()

        # load events
        self.__parseevents()

        # build friction matrix
        # ---------------------------------------------------------------------
        self.logger.info('Planting vegetation')
        self.__parse_vegetation()
        self.__build_friction()

        self.logger.info('set output path to: {}'.format(self.outputpath))
        self.logger.info('Initialised')

    def dash(self, dashtype=1, show=False):
        """ aaa """
        if dashtype == 1:
            self.__dash_geometry(show=show)
        elif dashtype == 2:
            self.__dash_output(show=show)

    def run(self, timesteps='all', progressbar=True):
        """
        solves the ode and stores results
        """
        if progressbar: self.logger.info('start model run')
        if timesteps == 'all':
            timesteps = self.grid.time

        # Run solver
        self.__bacsfort(timesteps, progressbar=progressbar)

        # Write output
        self.__write_output()

    # =============================================================================
    # private methods
    # =============================================================================
    def __parse_vegetation(self):
        config = ConfigParser()
        config.read(self.files.vegetation)
        vegetationnumbers = []
        for name in config:
            try:
                vegetationnumbers.append(int(name))
            except ValueError:
                pass
        if vegetationnumbers:
            self.logger.debug('Found vegetation numbers: {}'.format(vegetationnumbers))
            self.vegetationdefs = list(range(np.max(vegetationnumbers)))

        for name in config:
            try:
                int(name)
                d = config[name]
                self.vegetationdefs[int(name) - 1] = containers.VegetationContainer(name=d['name'],
                                                                          tstart=float(d['tstart']),
                                                                          tstop=float(d['tstop']),
                                                                          growthspeed=float(d['growthspeed']),
                                                                          decayspeed=float(d['decayspeed']),
                                                                          maximumcover=float(d['maximumcover']),
                                                                          density=float(d['density']),
                                                                          stemheight=float(d['stemheight']),
                                                                          stemdiameter=float(d['stemdiameter']))
            except ValueError:
                pass

    def __build_friction(self):
        """
        constructs the friction matrix

        todo: make friction a function, not a value, for compatibility with 2d
        friction uav maps
        """
        if self.parameters.frictionmodel.lower() == 'stationary':
            f = np.interp(self.grid.chainage,
                          self.grid.samples.X.T[0],
                          self.grid.samples.friction[0])
            self.grid.friction = pd.DataFrame(index=self.grid.time,
                                              columns=self.grid.chainage,
                                              data=[f] * len(self.grid.time))

        elif self.parameters.frictionmodel.lower() == 'vegetationgrowth':
            self.grid.vegetationdef = np.round(np.interp(self.grid.chainage,
                                                self.grid.samples.X.T[0],
                                                self.grid.samples.friction[1]))

            initial_blockage = np.interp(self.grid.chainage,
                                         self.grid.samples.X.T[0],
                                         self.grid.samples.friction[0])
            # SET INITIAL
            self.output.blockage[:] = np.tile(initial_blockage, (183, 1))

            # Solve

            tstart = [self.parameters.tstart + timedelta(days=self.vegetationdefs[int(vdef) - 1].tstart) for vdef in self.grid.vegetationdef]
            tstop = [self.parameters.tstart + timedelta(days=self.vegetationdefs[int(vdef) - 1].tstop) for vdef in self.grid.vegetationdef]
            growthspeed = np.array([self.vegetationdefs[int(vdef) - 1].growthspeed for vdef in self.grid.vegetationdef])
            decayspeed = np.array([self.vegetationdefs[int(vdef) - 1].decayspeed for vdef in self.grid.vegetationdef])
            maximumcover = np.array([self.vegetationdefs[int(vdef) - 1].maximumcover for vdef in self.grid.vegetationdef])
            dt = self.parameters.dt

            #r = pd.DataFrame(index=self.grid.time, columns=self.grid.chainage, data=0)

            for it, t in enumerate(self.grid.time[1:]):
                tmask = np.where((np.array(tstart) < t) & (np.array(tstop) >= t))
                dmask = np.where((np.array(tstart) < t) & (np.array(tstop) < t))

                # Growth cycle
                # --------------------------------------------------------------
                for mask, speed  in zip([tmask, dmask], [growthspeed, decayspeed]):
                    mask = mask[0]
                    N = self.output.blockage.iloc[it, mask]
                    N0 = initial_blockage[mask] - 1e-3
                    K = maximumcover[mask]
                    r = speed[mask]
                    #self.output.blockage.iloc[it + 1, mask] = r
                    self.output.blockage.iloc[it + 1, mask] = N + dt * r * (N - N0) * (1 - N / K)

                # Events
                # --------------------------------------------------------------
                for ie, event in enumerate(self.events):
                    if not(event.triggered) and (t > event.tstart):
                        self.logger.info('Triggered event {} on {}'.format(event, t))
                        chainagemask = (self.grid.chainage > event.minchainage) &\
                                       (self.grid.chainage < event.maxchainage)
                        #self.logger.debug(self.grid.chainagemask)
                        self.output.blockage.loc[t, self.grid.chainage[chainagemask]] = event.reduce_to
                        # Remove events
                        # TODO: make eventclass with mutable flag
                        self.events.pop(ie)
            """
            for ix, x in enumerate(self.grid.chainage):
05
                vdef = self.vegetationdefs[int(self.grid.vegetationdef[ix]) - 1]
                tstart = self.parameters.tstart + timedelta(days=vdef.tstart)
                tstop = self.parameters.tstart + timedelta(days=vdef.tstop)
                N0 = self.output.blockage.iloc[0, ix] - 1e-2
                (r, rd, K) = (vdef.growthspeed, vdef.decayspeed, vdef.maximumcover)
                dt = self.parameters.dt

                for it, t in enumerate(self.grid.time[1:]):
                    N = self.output.blockage.iloc[it, ix]
                    self.output.blockage.loc[t, x] = N

                    # growth model
                    # --------------------------------------------------------=
                    if (t > tstart) and (t <= tstop):
                        self.output.blockage.loc[t, x] = N + dt * r * (N - N0) * (1 - N / K)
                    elif t > tstop:
                        self.output.blockage.loc[t, x] = N + dt * rd * (N - N0) * (1 - N / K)

                    # Check for events
                    # --------------------------------------------------------=
                    for ie, event in enumerate(self.events):
                        if not(event.triggered) and (t > event.tstart):
                            self.logger.info('Triggered event {} on {}'.format(event, t))
                            chainagemask = (self.grid.chainage > event.minchainage) &\
                                           (self.grid.chainage < event.maxchainage)
                            #self.logger.debug(self.grid.chainagemask)
                            self.output.blockage.loc[t, self.grid.chainage[chainagemask]] = event.reduce_to
                            # Remove events
                            # TODO: make eventclass with mutable flag
                            self.events.pop(ie)
            self.logger.debug('output block shape {}'.format(self.output.blockage.shape))

            """
            self.grid.friction = self.blockagemodel(self.output.blockage)
    def __set_outputpath_from_config(self, configfilepath):
        # Set outputpath
        # ---------------------------------------------------------------------
        path, file = os.path.split(configfilepath)
        self.outputpath = os.path.join(path, 'output')
        try:
            os.mkdir(self.outputpath)
        except FileExistsError:
            pass

    def __dash_output(self, show=False):
        """
        """

        fig = plt.figure(figsize=(10, 4))
        axs = []
        axs.append(fig.add_subplot(211))
        axs.append(fig.add_subplot(223))
        axs.append(fig.add_subplot(224))

        axs[0].pcolormesh(self.grid.chainage, self.grid.time, self.output.waterdepth)
        myFmt = mdates.DateFormatter('%m')
        axs[0].yaxis.set_major_formatter(myFmt)
        axs[0].set_ylabel('Time [months]')
        axs[0].set_title('Waterdepth')

        axs[1].plot(self.grid.time, self.output.waterlevel.iloc[:, 0] - self.output.waterlevel.iloc[:, -1],
                    label='Model')
        axs[1].plot(self.grid.time, self.grid.upstream - self.grid.downstream, '.',
                    label='Measured')
        axs[1].legend()

        axs[2].plot(self.grid.time, self.grid.friction.iloc[:, 0])
        axs[2].set_title('friction')


        for ax in axs:
            utils.two_axes_style(ax)

        if show:
            plt.show()

    def __dash_geometry(self, show=False):
        fig = plt.figure(figsize=(10, 4))
        axs = []
        axs.append(fig.add_subplot(211))
        axs.append(fig.add_subplot(223))
        axs.append(fig.add_subplot(224))
        axs[0].pcolormesh(self.grid.X,
                        self.grid.Y,
                        self.grid.Z)

        axs[0].set_title('bathymetry')

        for ax in axs:
            utils.two_axes_style(ax)

        axs[1].plot(self.grid.time, self.grid.upstream - self.grid.downstream)
        axs[1].set_ylabel('Waterlevel slope [m]')
        axs[2].plot(self.grid.time, self.grid.discharge[0])
        axs[2].set_ylabel('Discharge $\mathrm{m^3/s}$')

        if show:
            plt.show()

    def __bacsfort(self, timesteps, progressbar=False):
        """
        Euler algorithm, backward in space, EACH TIME
        """
        if progressbar:
            iterator = enumerate(tqdm(timesteps))
        else:
            iterator = enumerate(timesteps)

        for i, t in iterator:
            self.__bacs(t)

    def __bacs(self, time):
        """
        Euler algorithm, backward in space
        """

        h = self.grid.downstream[time]
        d = h - self.grid.bedlevel[-1]
        waterlevels = [h]

        for i in list(reversed(range(len(self.grid.chainage))))[:-1]:
            Q = self.grid.discharge[self.grid.chainage[i]][time]
            depth = h - self.grid.bedlevel[i]
            slope = self.grid.bedslope[i]
            u = Q / self.grid.wet_area[i](h)
            n = self.grid.friction[self.grid.chainage[i]][time]
            R = self.grid.hydraulic_radius[i](h)
            dhdx = self.__belanger(slope, u, n, R, depth)
            d += -self.grid.dx * dhdx
            h = d + self.grid.bedlevel[i - 1]
            waterlevels.append(h)

        self.output.waterlevel.loc[time] = list(reversed(waterlevels))
        self.output.waterdepth.loc[time] = list(reversed(waterlevels)) - self.grid.bedlevel
        return waterlevels

    def __belanger(self, slope, u, n, R, depth):
        ode1 = self.parameters.g * slope
        ode2 = self.frictionmodel(u, n, R)
        ode3 = self.parameters.g - u ** 2 / depth
        return (ode1 + ode2) / ode3

    def __manning(self, u, n, R):
        return -self.parameters.g * (u * n) ** 2 / (R**(4 / 3.))

    def __parsegeometry(self):
        """
        Excel file readers
        """
        if self.parameters.geometrytype == 1:
            (X, Y, Z, friction) = self.__read_exceltype_1()
            self.grid = containers.GeometryGrid(parameters=self.parameters,
                                                logger=self.logger,
                                                samples = dict(X=X, Y=Y, Z=Z, friction=friction))
        else:
            self.logger.error('unknown geometry filetype')

    def __read_exceltype_1(self):
        """
        Trapezoidal profile excel input
        """

        geom = pd.read_excel(self.files.geometry)

        X = list()
        Y = list()
        Z = list()
        for i in geom.index:
            ygpoints = np.array([-self.maxwidth / 2, 0 - geom['W'][i] / 2, 0 + geom['W'][i] / 2, self.maxwidth / 2])
            ypoints = np.append(ygpoints, np.linspace(-self.maxwidth / 2, self.maxwidth / 2, 20))
            ypoints = np.sort(ypoints, axis=0)
            Y.append(ypoints)

            # height of profile at maxwidth
            dz = (self.maxwidth - geom['W'][i]) / 2 * np.tan(np.deg2rad(geom['angle'][i]))

            zgpoints = np.array([dz + geom['z'][i], geom['z'][i], geom['z'][i], dz + geom['z'][i]])
            zpoints = np.interp(ypoints, ygpoints, zgpoints)
            Z.append(zpoints)

        X = [geom['x']] * len(ypoints)

        return np.array(X).T, np.array(Y), np.array(Z), [geom['friction_01'], geom['friction_02']]

    def __parseconfigfile(self, input_file):
        [pathname, filename] = os.path.split(input_file)

        config = ConfigParser()
        config.read(input_file)

        files = dict()
        for key in config['files']:
            files[key] = os.path.join(pathname, config['files'][key])

        parameters = utils.parse_dict({'parameters': config['parameters']},
                                      typedict=utils.configtypes)

        self.files = containers.FileContainer(**files)
        self.parameters = containers.ParameterContainer(**parameters['parameters'])

    def __parsemeasurements(self):
        """
        parses the 'measurements' file and sets boundary conditions
        """
        laterals = None
        if os.path.isfile(self.files.measurements):
            parser = lambda date: pd.datetime.strptime(date, self.dateformat)
            data = pd.read_csv(self.files.measurements, header=0,
                                                        parse_dates=[0],
                                                        date_parser=parser)
            if (self.files.laterals is not None) and (os.path.isfile(self.files.laterals)):
                laterals = pd.read_csv(self.files.laterals)

            self.grid.set_boundaries(data, laterals)

        else:
            self.logger.error('Measurement file not found')

    def __write_output(self):
        """
        Export output to csv file
        """
        for variable, data in self.output._asdict().items():
            fpath = os.path.join(self.outputpath, '{}.csv'.format(variable))
            data.to_csv(fpath)
            self.logger.info('Written {} to {}'.format(variable, fpath))

    def __parseevents(self):
        self.events = list()
        config = ConfigParser()
        config.read(self.files.events)
        for event in config:
            if not(event=="DEFAULT"):
                self.logger.debug('Parsing event {}'.format(event))
                eventdict = {'event': config[event]}
                eventdict['event']['name'] = event
                eventdict['event']['triggered'] = 'False'
                pdict = utils.parse_dict(eventdict, typedict=utils.eventypes)['event']
                self.events.append(containers.EventContainer(**pdict))
        self.logger.debug('events: {}'.format(self.events))

    @staticmethod
    def pglinear(blockage):
        return 0.0333 / (1 - blockage)
