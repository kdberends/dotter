# Dotter model prototype
#
# This file: Supporting classes
# 
# Author: Koen Berends
# Contact: k.d.berends@utwente.nl or koen.berends@utwente.nl
# Copyright (c) 2016 University of Twente & Deltares
#

# region // imports
from xlrd import open_workbook
import numpy as np
import matplotlib.pyplot as plt
from configparser import ConfigParser
from datetime import datetime
import csv
# endregion


class Event:
    """
    Container for data of events. 
    """
    def __init__(self, eventdict):
        self.name = eventdict['name']
        self.effectivity = float(eventdict['effectivity'])
        self.min = float(eventdict['range_min'])
        self.max = float(eventdict['range_max'])
        self.time = datetime.strptime(eventdict['time'], '%d/%m/%Y')
        self.triggered = False
        

class Stream:
    """
    Main 'model/geometry' class. 
    """
    def __init__(self, **parameters):
        self.path = parameters['path']
        self.result = {'x': [], 'waterdepth': [], 'A': [], 'R': []}
        self.parameters = parameters
        self.grid_rg = []
        self.grid_rd = []
        self.grid_K = []
        self.laterals = []
        self.grid_grow_season = []
        self.grid_grow_month = []
        self.events = np.array([])
        self.events_trigger = np.array([])

    def load_from_excel(self, filename):
        """
        Loads excel file with stream geometry
        """
        x = list()
        theta = list()
        w = list()
        blockage = list()
        vg = list()
        zb = list()
        sl = list()
        n = list()
        maxh = list()
        wb = open_workbook(filename)
        for s in wb.sheets():
            self.name = s.name
            for row in range(1, s.nrows):
                x.append(s.cell(row, 0).value) # afstand langs beek
                w.append(s.cell(row, 1).value) # bodembreedte
                theta.append(s.cell(row, 2).value) # taludhoek
                zb.append(s.cell(row, 3).value) # bodemhoogte
                sl.append(s.cell(row, 4).value) # slibhoogte
                n.append(s.cell(row, 5).value) # manning coefficient
                blockage.append(s.cell(row, 6).value)  # == 1-vo (dus fa = percentage begroeiing)
                vg.append(s.cell(row, 7).value) # vegetatie-id
                try:
                    maxh.append(s.cell(row, 8).value) # maximale waterstand (ter visualisatie)
                except:
                    maxh.append(0)

        self.x = np.array(x)
        self.zb = np.array(zb)
        self.sl = np.array(sl)
        self.n = np.array(n)
        self.bedwidth = np.array(w) 
        self.theta = np.array(theta)
        self.blockage = np.array(blockage)
        self.vg = np.array(vg) 
        self.maxh = np.array(maxh)
        self.length = max(self.x)

    def generate_grid(self):
        user_stepsize = self.parameters['dx']
        number_of_steps = int((np.max(self.x)-np.min(self.x))/user_stepsize)
        self.grid_x = np.linspace(min(self.x), max(self.x), number_of_steps+1)
        self.grid_z = np.interp(self.grid_x, self.x, self.zb + self.sl)
        self.grid_n = np.interp(self.grid_x, self.x, self.n)
        self.grid_i = np.abs(np.diff(self.grid_z) / np.diff(self.grid_x))
        self.grid_i = np.append([self.grid_i[0]], self.grid_i)
        self.grid_i = np.round(self.grid_i, decimals=10)
        self.grid_width = np.interp(self.grid_x, self.x, self.bedwidth)
        self.grid_theta = np.interp(self.grid_x, self.x, self.theta)
        self.grid_blockage = np.interp(self.grid_x, self.x, self.blockage)
        self.grid_vg = np.interp(self.grid_x, self.x, self.vg)
        self.grid_maxh = np.interp(self.grid_x, self.x, self.maxh)

        self.grid_Q = np.array([self.parameters['Q']] * len(self.grid_x))
        for lateral in self.laterals:
            self.grid_Q[np.where(self.grid_x >= lateral[0])] += lateral[1]

    def get_data_at_index(self, index):
        """
        Returns spatial data at point along stream
        """
        return (self.grid_x[index], 
                self.grid_z[index], 
                self.grid_i[index], 
                self.grid_n[index],
                self.grid_Q[index],
                self.grid_width[index],
                self.grid_theta[index])

    def hydraulic_parameters(self, width, theta, depth):
        """
        Hydraulic radius of trapezoidal cross=section
        theta = outer angle
        """
        A = width * depth + depth ** 2 / np.tan(theta * np.pi / 180.)
        P = width + 2 * depth / np.sin(theta * np.pi / 180.)
        R = A / P
        return A, P, R

    def get_parameters_at_index(self, index=0, waterdepth=0):
        """
        Updates the parameter dictionary
        """
        (x, z, i, n, Q, width, theta) = self.get_data_at_index(index)
        (A, P, R) = self.hydraulic_parameters(width, theta, waterdepth)
        self.parameters['x'] = x
        self.parameters['z'] = z
        self.parameters['R'] = R
        self.parameters['A'] = A
        self.parameters['ib'] = i
        self.parameters['n'] = float(n)
        self.parameters['Q'] = Q
        return self.parameters

    def set_vegetation(self, filename):
        for i, vg in enumerate(self.grid_vg):
            (rg, rd, K, grow_month, grow_season) = self.get_vegetation_parameters(str(int(vg)), filename)
            self.grid_rg.append(rg)
            self.grid_rd.append(rd)
            self.grid_K.append(K)
            self.grid_grow_month.append(grow_month)
            self.grid_grow_season.append(grow_season)

    @staticmethod
    def get_vegetation_parameters(name, filename):
        config = ConfigParser()
        config.read(filename)

        rg = float(config[name]['growth_malthusian'])
        rd = float(config[name]['death_malthusian'])
        K = float(config[name]['msp'])
        grow_month = int(config[name]['grow_month'])
        grow_season = int(config[name]['grow_season'])

        return rg, rd, K, grow_month, grow_season

    def plot(self, language='English'):
        # ---------------------------------
        # VISUALISATION
        # ----------------------------------
        labeldict = {'english': {'waterdepth': 'waterdepth [m]',
                                 'waterlevel': 'waterlevel [m]',
                                 'distance': 'Distance [m]'},
                     'dutch': {'waterdepth': 'Waterdiepte [m]',
                               'waterlevel': 'Waterstand [m]',
                               'distance': 'Afstand [m]',
                               'bed': 'Bodem',
                               'silt': 'Sliblaag',
                               'manning': 'Manning coefficient'}
                     }
        ax = [[], []]
        labels = labeldict[language.lower()]

        fig = plt.figure(self.name, figsize = (20, 10), dpi=60)
        ax[1] = fig.add_subplot(212)
        ax[0] = fig.add_subplot(211, sharex=ax[1])
        #ax[1] = fig.add_subplot(312, sharex=ax[2])
        
        ax[0].plot(self.x, 
                   self.zb,
                   color=[0.35, 0.25, 0.15],
                   linestyle='-', label='Bodemhoogte')
        ax[0].plot(self.x, 
                   self.maxh,
                   color=[0.9, 0.25, 0.15],
                   linestyle='-', label='Maximaal toelaatbare waterstand')
        ax[0].plot(self.result['x'], 
                   self.result['waterlevel'],
                   color='b', label='Waterstand')
        tax = ax[0].twinx()
        tax.plot(self.grid_x, self.result['waterlevel'] - self.grid_maxh, '--k')
        tax.set_ylim([-0.5, 0.5])
        ax[0].plot([np.max(self.result['x'])] * 2,
                   [self.zb[-1], self.result['waterlevel'][-1]],
                   color='k',
                   linewidth=5)
        ax[0].set_ylabel(labels['waterlevel'])
        ax[0].legend()

        ax[1].plot(self.grid_x, 
                   self.grid_n,
                   color=[0.2, 0.8, 0.2])
        ax[1].set_ylabel(labels['manning'])
        ax[1].set_xlabel(labels['distance'])
        ax[1].set_ylim([0, max(self.grid_n) * 2])
        
        plt.setp(ax[0].get_xticklabels(), visible=False)
        #plt.tight_layout()
        
        return fig, ax

    def event_callback(self, current_time):
        triggered_events = None
        if any(current_time > self.events_trigger):
            triggered_events_indices = np.where(current_time > self.events_trigger)
            triggered_events = self.events[triggered_events_indices]
            self.events_trigger = np.delete(self.events_trigger, triggered_events_indices)
            self.events = np.delete(self.events, triggered_events_indices)

        return triggered_events

    def load_events(self, filename):
        
        config = ConfigParser()
        config.read(filename)

        for section in config:
            try:
                float(section)
                new_event = Event(config[section])
                self.events = np.append(self.events, new_event)
                self.events_trigger = np.append(self.events_trigger, new_event.time)
            except ValueError:
                pass

    def load_laterals(self, filename):
        laterals = list()
        
        with open(filename, 'r') as f:
            reader = csv.reader(f)
            header = True
            for line in reader:
                if header:
                    header = False
                else:
                    laterals.append(list(map(float, line)))
        self.laterals = np.array(laterals)