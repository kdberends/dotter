# Dotter model prototype
#
# This file: Supporting functions
# 
# Author: Koen Berends
# Contact: k.d.berends@utwente.nl or koen.berends@utwente.nl
# Copyright (c) 2016 University of Twente & Deltares
#


# region // Imports
import numpy as np
from matplotlib import animation, image
import matplotlib.pyplot as plt
import sys
import pandas as pd
from configparser import ConfigParser
from datetime import datetime, timedelta
import seaborn as sns
from . import settings
# endregion

# Roughness Models
# ===========================================================
def chezy(waterdepth=1, velocity=1, **parameters):
    """
    Chezy friction formula. 

    Keyword arguments:
        waterdepth : float, the water depth (m)
        velocity   : float, flow velocity (m/s)
        parameters : dict, dictionary containing parameters   

    Returns:
        Friction acceleration term (m/s**2)
    """
    return -(parameters['g']*velocity**2)/(parameters['C']**2*waterdepth)

def prony(waterdepth=1, velocity=1, **parameters):
    """
    Prony friction formula. 

    Keyword arguments:
        waterdepth : float, the water depth (m)
        velocity   : float, flow velocity (m/s)
        parameters : dict, dictionary containing parameters   

    Returns:
        Friction acceleration term (m/s**2)
    """
    return -parameters['g']/waterdepth*(4.44499e-5*velocity+3.093140e-4*velocity**2)

def manning(waterdepth=1, velocity=1, **parameters):
    """
    Manning friction formula. 

    Keyword arguments:
        waterdepth : float, the water depth (m)
        velocity   : float, flow velocity (m/s)
        parameters : dict, dictionary containing parameters   

    Returns:
        Friction acceleration term (m/s**2)
    """

    #return -parameters['g']*(velocity*parameters['n'])**2/(waterdepth**(4/3.))
    return -parameters['g']*(velocity*parameters['n'])**2/(parameters["R"]**(4/3.))

def darcyweisbach(waterdepth=1, velocity=1, **parameters):
    """
    Darcy-Weisbach friction formula. 

    Keyword arguments:
        waterdepth : float, the water depth (m)
        velocity   : float, flow velocity (m/s)
        parameters : dict, dictionary containing parameters   

    Returns:
        Friction acceleration term (m/s**2)
    """
    return -parameters['fd']*velocity**2/(8*waterdepth)


# Roughness Predictors
# ===========================================================

def malthusian_growth(current_value, malthusian, capacity):
    """
    Logistic differential function or "Malthusian growth model"
    dv/dt = rv(K-v)/K
    
    Arguments:
    current_value - float
    malthusian - float, 'growth rate'
    capacity - float, maximum value
    """
    return malthusian * current_value * (capacity + 0.01 - current_value) / capacity

def pitlogriffioenlinear(geometry, time=None, timestep=None):
    """
    Pitlo-Griffioen roughness predictor based on percentage open water. 

    Roughness model - Manning
    

    Arguments:
        geometry - 'Dotter.Stream' object
        time - datetime
        timestep - float (days)

    Returns:
        roughness formula - func
        geometry - modified Dotter.Stream object
    """

    if time is not None:
        for i, fa in enumerate(geometry.grid_fa):
            # Check whether in growth season
            if (time.month > geometry.grid_grow_month[i]) & (time.month < (geometry.grid_grow_month[i] + geometry.grid_grow_season[i])):
                # In growth season, use 'rg'
                geometry.grid_fa[i] += malthusian_growth(geometry.grid_fa[i], geometry.grid_rg[i], geometry.grid_K[i]) * timestep
            else:
                # not in growth season,  use 'rd'
                geometry.grid_fa[i] += malthusian_growth(geometry.grid_fa[i], geometry.grid_rd[i], geometry.grid_K[i]) * timestep
                geometry.grid_fa[i] = np.max([geometry.grid_fa[i], 0.01])
                geometry.grid_fa[i] += malthusian_growth(geometry.grid_fa[i], geometry.grid_rd[i], geometry.grid_K[i]) * timestep
            
            # 'Pitlo-Griffioen model' (see Pitlo & Griffioen, 1991)
            # Oorspronkelijk model: vo = 3.33*km 
            # 100 * fa = 100-vo
            # n = 0.0333 / (1 - fa)
            #   
            geometry.grid_n[i] = 0.0333 / (1 - geometry.grid_fa[i])
            
    return manning, geometry

def constantmanning(geometry, time=None, timestep=None):
    """

    """
    if geometry.parameters['use_homogeneous_friction']:
        geometry.grid_n = [geometry.parameters['constant_friction']] * len(geometry.grid_n)
    else:
        pass

    return manning, geometry 

# Data functions
# ===========================================================

def dateparse(time_in):    
    return datetime.strptime(time_in.strip(), '%Y-%m-%d')

def read_results_file(filename):
    """
    Reads result file from Dottermodel

    Arguments:

    filename: str, path to file

    Returns:

    pandas dataframe
    """
    return pd.read_csv(filename, 
                       parse_dates=True, 
                       date_parser=dateparse, 
                       index_col=0)

def load_plotting_config(filename):
    
    config = ConfigParser()
    config.read(filename)
    prefs = {'1d': dict(), '2d': dict()}
    for section in config:
        for key in config[section]:
            value = config[section][key]
            if key.lower()=='make_plot':
                prefs[section.lower()][key] = bool(value)
            elif key.lower() == 'type':
                prefs[section.lower()][key] = int(value)
            elif key.lower() == 'plot_moments':
                prefs[section.lower()][key] = map(dateparse, value.split(','))
            elif key.lower() == 'plot_locations':
                prefs[section.lower()][key] = map(float, value.split(','))
            elif key.lower() == 'colormap':
                prefs[section.lower()][key] = value.strip()
                
    return prefs

# Numerical solvers
# ===========================================================

def belanger(waterdepth=0, friction=chezy, **parameters):
    """
    The Belanger equation

    Keyword arguments:
        waterdepth : float, the water depth (m) at previous step
        friction   : func, a friction function
        parameters : dict, dictionary containing parameters   

    Returns:
        dh/dx
    """
    #velocity = parameters['Q'] / float(waterdepth)
    velocity = parameters['Q'] / parameters['A']
    return (parameters['g'] * parameters['ib'] + 
            friction(waterdepth, velocity, **parameters)) / \
           (parameters['g'] - velocity ** 2 / waterdepth)

def spatial_solver(ode, frictionformula=chezy, geometry=None, verbose=False):
    """
    Solves the ode numerically using the (inverse forward) Euler method.

    Arguments:
        ode     : func, a function of the ODE

    Keyword arguments:
        initial_value  : float, the initial value 
        stepsize       : int, the stepsize. If solution is unstable, decrease the stepsize
        number_of_steps: int, maximum number of steps to take
        frictionformula: func, a function which returns the friction acceleration term
        parameters     : dict, a dictionary containing parameters


    Returns:
        xlist   : numpy array, containing the x locations
        ylist   : numpy array, containing the water depths
    """
    # Initialise the lists 
    ylist = [geometry.parameters['h']]

    # Set the initial value (y=waterdepth)
    y = ylist[0]

    # We start from zero and loop over all the steps...
    
    number_of_steps = len(geometry.grid_x)-1

    if verbose: sys.stdout.write('[')
    for i in range(number_of_steps):
        if not(i==number_of_steps-1):
            stepsize = geometry.grid_x[i] - geometry.grid_x[i + 1]
        parameters = geometry.get_parameters_at_index(index=number_of_steps - i - 1, waterdepth=y)
        
        # Calculate the value of the next point
        y += stepsize * ode(y, frictionformula, **parameters)

        # Append the values to the output lists
        ylist.append(y)
        
        if verbose & (100 * i / float(number_of_steps)).is_integer():
            sys.stdout.write('|')
    if verbose: sys.stdout.write(']')
    if verbose: print('')

    # Write output to geometry object
    geometry.result['x'] = geometry.grid_x
    geometry.result['waterdepth'] = np.flipud(np.array(ylist)   )
    geometry.result['waterlevel'] = np.flipud(np.array(ylist) + np.flipud(np.array(geometry.grid_z)) )

    return geometry

def qs_solver(geometry=None, restart=None, verbose=True, write_to_file=True, starttime=0, stoptime=365):
    """
    Quasi-stationary solver. 
    """

    # Roughness predictor
    if geometry.parameters['friction_model'].lower() == 'pitlogriffioenlinear':
        roughness_predictor = pitlogriffioenlinear
    elif geometry.parameters['friction_model'].lower() == 'manning':
        roughness_predictor = constantmanning
        
    frictionformula, forget = roughness_predictor(geometry)

    # Generate time grid 
    # ------------------------------------------------------------------------
    t = geometry.parameters['t']
    time = []
    timestep = geometry.parameters['dt']
    for i in np.arange(starttime, stoptime, timestep):
        t += timedelta(days=timestep)
        time.append(t)

    # Pre-allocate spatiotemporal result files
    # ------------------------------------------------------------------------
    if restart is None:
        results = {'x': pd.DataFrame(index=time, columns=geometry.grid_x),
                   'n': pd.DataFrame(index=time, columns=geometry.grid_x),
                   'fa': pd.DataFrame(index=time, columns=geometry.grid_x),
                   'waterdepth': pd.DataFrame(index=time, columns=geometry.grid_x),
                   'waterlevel': pd.DataFrame(index=time, columns=geometry.grid_x),
                   'max_allowed_waterlevel': pd.DataFrame(index=time, columns=geometry.grid_x),
                   'discharge': pd.DataFrame(index=time, columns=geometry.grid_x)
                   }
    else:
        sys.stdout.write('Using previous results \n')
        results = restart
        
    # Solve initial condition
    # ------------------------------------------------------------------------
    geometry = spatial_solver(belanger, 
                              frictionformula=frictionformula,
                              geometry=geometry)

    results['n'].loc[time[0]] = geometry.grid_n

    results['waterdepth'].loc[time[0]] = geometry.result['waterdepth']
    results['waterlevel'].loc[time[0]] = geometry.result['waterlevel']

    # Start time-loop
    # ------------------------------------------------------------------------
    if verbose: sys.stdout.write('[Starting calculation] \n')
    
    updnext = 0
    for it, t in enumerate(time):
        
        # Pass events
        # ---------------------------------------------------------------------
        if np.any(geometry.events):
            events = geometry.event_callback(pd.to_datetime(t))
            if np.any(events):
                for event in events:
                    if verbose: sys.stdout.write('Event %s Triggered on %s\n' % (event.name, t))
                    geometry.grid_fa[(geometry.grid_x > event.min) & (geometry.grid_x < event.max)] *= (1 - event.effectivity)

        # Vegetation growth model 
        # ---------------------------------------------------------------------         
        forget, geometry = roughness_predictor(geometry, t, timestep)

        # Hydraulic solver
        # ---------------------------------------------------------------------  
        geometry = spatial_solver(belanger, 
                                  frictionformula=frictionformula,
                                  geometry=geometry)

        # write results
        # ---------------------------------------------------------------------  
        results['x'].loc[t] = geometry.grid_x
        results['n'].loc[t] = geometry.grid_n
        results['discharge'].loc[t] = geometry.grid_Q
        results['fa'].loc[t] = geometry.grid_fa
        results['waterdepth'].loc[t] = geometry.result['waterdepth']
        results['waterlevel'].loc[t] = geometry.result['waterlevel']
        results['max_allowed_waterlevel'].loc[t] = geometry.grid_maxh
        # Update on progress
        # ---------------------------------------------------------------------  
        progress = (100 * it / float(len(time))) 
        if verbose & (progress > updnext):
            updnext += 10
            sys.stdout.write(':: %i%%\n' % updnext)

    # Write output
    # -------------------------------------------------------------------------
    if verbose: sys.stdout.write('[Computation done]\n')
    if write_to_file:
        if verbose: sys.stdout.write('Writing output to file...\n')
        results['n'].to_csv('roughness.csv', sep=',')
        if verbose: sys.stdout.write('Roughness written to roughness.csv\n')
        results['discharge'].to_csv('discharge.csv', sep=',')
        if verbose: sys.stdout.write('Discharge written to discharge.csv\n')
        results['waterdepth'].to_csv('waterdepth.csv', sep=',')
        if verbose: sys.stdout.write('Waterdepth written to waterdepth.csv\n')
        results['waterlevel'].to_csv('waterlevel.csv', sep=',')
        if verbose: sys.stdout.write('Waterlevel written to waterlevel.csv\n')
        results['max_allowed_waterlevel'].to_csv('max_allowed_waterlevel.csv', sep=',')
        if verbose: sys.stdout.write('Max. waterlevel written to max_allowed_waterlevel.csv\n')
        results['fa'].to_csv('percentagebegroeiing.csv', sep=',')
        if verbose: sys.stdout.write('percentagebegroeiing written to percentagebegroeiing.csv\n')

    return results, geometry

# Plotting functions
# ===========================================================
flatui = ["#9b59b6", "#e74c3c", "#3498db", "#95a5a6", "#34495e", "#2ecc71"]

def set_style():
    """
    'paper-ready' plotting style
    """
    sns.set(context='paper', 
            style='whitegrid', 
            palette=flatui, 
            font='serif', 
            font_scale=2, 
            color_codes=False, 
            rc=None)

def plot_1d_type01(geometry, results, plot_moments, plot_places):
    """
    Plots results from time_loop function
    """
    r = results

    time = pd.to_datetime(r['waterdepth'].index.values)

    fig = plt.figure()
    ax1 = fig.add_subplot(311)

    moments = []
    for m in plot_moments:
        moments.append(time[time > m])

    for m in moments:
        r['waterdepth'].T[m[0]].plot(ax = ax1, label=m[0])

    ylim = ax1.get_ylim()
    for x in plot_places:
        ax1.plot([x, x], ylim, '-k')

    ax1.legend()
    ax1.set_title('Profiel')
    ax1.set_ylabel('Waterdiepte [m]')

    ax2 = fig.add_subplot(312)

    for ix in plot_places:
        ix = np.where(geometry.grid_x > ix)
        ax2.plot(time, r['n'].T.iloc[ix[0][0]].values, label="X={} km".format(geometry.grid_x[ix[0][0]]))


    ylim = ax2.get_ylim()

    for m in moments:
        ax2.plot([m[0], m[0]], [ylim[0], ylim[1]], '-k')

    ax2.legend()
    ax2.set_ylabel('Ruwheid (Manning)')
    ax3 = fig.add_subplot(313)
    x = r['waterdepth'].columns.values
    d = r['waterlevel'][x[0]] - r['waterlevel'][x[-1]]
    d.plot(ax=ax3)
    ylim = ax3.get_ylim()

    for m in moments:
        ax3.plot([m[0], m[0]], [ylim[0], ylim[1]], '-k')

    ax3.set_ylabel('Verhang')

def plot_1d_type02(geometry):
    fig, ax = geometry.plot(language='Dutch')
    ax[0].set_title('Resultaten op laatste tijdstap')
    fig.tight_layout()  
    return fig, ax

def plot_1d_type03(gm, xlim):
    """
    This plot type is made for the virtual river serious game
    """
    ax = [[], [], []]
    fig = plt.figure(gm.name, figsize=(20, 10), dpi=60)
    ax[0] = fig.add_subplot(311)  # Upper figure, water level goal constrained
    ax[1] = fig.add_subplot(312)  # middle figure, roughness, constrained
    ax[2] = fig.add_subplot(313)  # bottom figure, full view, constrained

    # Top figure
    rwvl = gm.result['waterlevel'] - gm.grid_maxh
    mask = rwvl < 0
    ax[0].plot(gm.grid_x[mask], rwvl[mask] , '--', color='darkblue', linewidth=2)
    ax[0].plot(gm.grid_x[~mask], rwvl[~mask] , '--', color='darkred', linewidth=2)
    ax[0].plot(gm.x, [0] * len(gm.x), '-k', linewidth=2)
    ax[0].set_xlim(xlim)
    ax[0].set_title('Critical water level - actual water level')
    ax[0].set_ylim([-0.2, 0.2])

    # Centre figure
    ax[1].plot(gm.grid_x, gm.grid_n, '-', color='seagreen')
    ax[1].set_xlim(xlim)
    ax[1].set_ylim([0, 0.08])
    ax[1].set_ylabel('Manning coefficient')
    
    # Bottom figure
    ax[2].plot(gm.x, gm.zb)
    ax[2].plot(gm.grid_x, gm.result['waterlevel'], '-b', label='Waterlevel')
    ax[2].plot(gm.grid_x, gm.grid_maxh, '--r', label='Dike height')
    ax[2].set_ylabel('z [m + NAP]')
    
    settings.set_plotstyle(scale=2)
    for a in ax:
        settings.two_axes_style(a)


    return fig, ax

def plot_2d(data, colormap='viridis', threshold=False):
    """

    """
    x = list(map(float, data.columns))
    t = data.index
    values = np.array(data.values)

    # generate 2 2d grids for the x & y bounds

    fig = plt.figure()

    ax = fig.add_subplot(111)

    if threshold:
        im = ax.pcolor(x, t, values > threshold, cmap=colormap)
    else:
        im = ax.pcolor(x, t, values, cmap=colormap)
    ax.set_xlabel('Afstand [m]')
    fig.colorbar(im)
    fig.tight_layout()
    return fig, ax



def animate_results(wdata, rdata, bpdata, tofile=False, fname='default', locs=[1000], vformat='wmv'):
    """
    Animates results
    Use 'read_results_file' function to obtain data objects from 
    files

    Arguments:
    wdata = data object from waterlevel.csv or waterdepth.csv
    rdata = data object from roughness.csv
    bpdata= data object from percentagebegroeiing.csv

    Keyword arguments:
    tofile: boolean,  if true, write to file. If false, opens interactive window
    fname : str, if writing to file, the filename
    vformat: str, either 'wmv' or 'mp4'
    locs : list of floats. Locations to plot growcurve


    """
    locs = [1000.0, 2000.0]
    x = np.array(map(float, wdata.columns))
    t = wdata.index
    values = np.array(wdata.values)
    rvalues = np.array(rdata.values)
    bpvalues = np.array(bpdata.values)
    fig = plt.figure(figsize=(10, 6), dpi=90)
    #fig = plt.figure()
    ax = fig.add_subplot(311)
    ax2 = fig.add_subplot(312, sharex=ax)
    ax3 = fig.add_subplot(313)
    plt.xticks(rotation=30, ha='right')
    plt.gcf().subplots_adjust(bottom=0.15)
    ax3.get_xaxis().set_tick_params(direction='out')
    idlocs = list()
    for loc in locs:
        idloc = np.argwhere(np.array(x) == loc)[0][0] 
        idlocs.append(idloc)
        ax3.plot(t, bpvalues.T[idloc] * 100, '--', color=[0.3] * 3)

    ax.set_ylabel('Waterstand [m]')
    ax2.set_ylabel('Manning')
    for idloc in idlocs:
        ax2.plot(x[idloc], [0.01], 'o', color=[0.3] * 3)
    ax3.set_ylabel('Begroeiing [%]')
    bplines = list()
    bpmarker = list()
    line, = ax.plot([], [], '-', color=flatui[2])
    iline, = ax.plot([], [], '--', color=flatui[2])
    iline2, = ax.plot([], [], '.-', color='r', linewidth=1)
    rline, = ax2.plot([], [], '-', color=flatui[1])
    for i in locs:
        temp, = ax3.plot([], [], '-', color=flatui[0])
        bplines.append(temp)
        temp, = ax3.plot([], [], 'or')
        bpmarker.append(temp)
    ax.set_xlim(np.min(x), np.max(x))
    ax.set_ylim(np.min(values)*0.9, np.max(values) * 1.2)
    ax2.set_ylim(0, np.max(rvalues) * 1.2)
    ax3.set_ylim(-20, 100)
    fig.suptitle('Dottermodel (2017)')
    ax.set_title('Contact: koen.berends@deltares.nl', fontsize=10)
    tstamp = ax.text(np.mean(x), np.max(values) * 1.19, '', horizontalalignment='center', verticalalignment='top')
    inum = ax.text(x[0]+10, 0, '', horizontalalignment='left', verticalalignment='center')

    im = image.imread('deltares.jpg')
    ax.imshow(im, aspect='auto', extent=(2200, 2700, 4, 5.7), zorder=-1)

    def init():
        tstamp.set_text('')
        inum.set_text('')
        line.set_data([0], [0])
        iline.set_data([0], [0])
        iline2.set_data([0], [0])
        rline.set_data([0], [0])
        for i, bpline in enumerate(bplines):
            bpline.set_data(t[0], [0])
            bpmarker[i].set_data(t[0], [0])
        return tuple(bplines) + tuple(bpmarker) + (tstamp, )+ (inum, ) + (line, ) + (rline, ) + (iline, ) + (iline2, )

    def update_line(it, values, rvalues, bpvalues):
        y = values[it]
        ry = rvalues[it]
        
        if it > 1:
            for i, bpline in enumerate(bplines):
                bp = bpvalues.T[idlocs[i]] * 100
                bp = bp[:it]
                bplines[i].set_data(t[:it], bp)
                bpmarker[i].set_data(t[it-1], bp[-1])
        rline.set_data(x, ry)
        tstamp.set_text('{}'.format(t[it]))
        
        line.set_data(x, y)   
        iline.set_data(x, x ** 0 * y[-1])
        iline2.set_data(x[0] * 2 + 10, [y[-1], y[0]])
        inum.set_y(np.mean([y[-1], y[0]]))
        inum.set_text('{:.2f}'.format(y[0] - y[-1]))
        
        #return line, rline, bpline, bpmarker, tstamp
        return tuple(bplines) + tuple(bpmarker) + (rline, ) + (tstamp, )+ (inum, ) + (line, ) + (iline, ) + (iline2, )

    


    anim = animation.FuncAnimation(fig, update_line,
                                   fargs=(values, rvalues, bpvalues),
                                   frames=len(t),
                                   interval=50,
                                   init_func=init,
                                   blit=True,
                                   repeat=True)
    if tofile:
        if vformat == 'mp4':
            anim.save('{}.mp4'.format(fname), fps=25, extra_args=['-vcodec', 'libx264'])
            sys.stdout.write(':: animation written to {}.mp4'.format(fname))
        elif vformat == 'wmv':
            anim.save('{}.wmv'.format(fname), fps=25, extra_args=['-q:a', ' 2', '-acodec', 'wmav2'])
            sys.stdout.write(':: animation written to {}.wmv'.format(fname))
        else:
            sys.stdout.write(":: ERROR: unknown video format. Use either 'wmv' or 'mp4'")
    else:
        plt.show()