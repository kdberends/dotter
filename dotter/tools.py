"""
Dottermodel tools
"""

# =============================================================================
# Imports
# =============================================================================
from copy import copy
from . import utils
from scipy import optimize
import numpy as np 
from tqdm import tqdm
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
# =============================================================================
# Tools
# =============================================================================

def maaibos(model, discharges=None, normative_friction=0.03, show=False):
    """
    predictive use (at different Q)
    """
    
    # generate QH
    # -------------------------------------------------------------------------
    discharge_restore = copy(model.grid.discharge)
    friction_restore = copy(model.grid.friction)

    model.logger.info('Generating QH for normative friction')
    QHdischarges = np.linspace(model.grid.discharge.min().min(), 
                               np.max([model.grid.discharge.max().max(),
                                       np.max(discharges)]),
                               20)
    
    model.grid.friction[:] = normative_friction
    upstream = list()
    for discharge in tqdm(QHdischarges):
        model.grid.set_discharge(discharge)
        model.run(timesteps=[model.grid.time[0]], progressbar=False)
        upstream.append(model.output.waterlevel.iloc[0, 0])

    qh = interp1d(QHdischarges, upstream, bounds_error=False, fill_value=0)

    # Restore previous conditions
    # -------------------------------------------------------------------------
    model.grid.discharge = copy(discharge_restore)
    model.grid.friction = copy(friction_restore)

    # Get scenario results
    # -------------------------------------------------------------------------
    discharge_levels = list()
    for discharge in discharges:
        model.grid.set_discharge(discharge)
        model.logger.info('Running model for q = {}'.format(discharge))
        model.run()
        discharge_levels.append(copy(model.output.waterlevel.iloc[:, 0]))


    # Restore previous conditions
    # -------------------------------------------------------------------------
    model.grid.discharge = copy(discharge_restore)
    model.grid.friction = copy(friction_restore)

    model.run()
    modelled = model.output.waterlevel.iloc[:, 0]

    # Get measurements
    # -------------------------------------------------------------------------
    fig = plt.figure(figsize=(10,4))
    axs = list()
    axs.append(fig.add_subplot(211))
    axs.append(fig.add_subplot(212))

    scenslopes = []
    for i, discharge in enumerate(discharges):
        scenslopes.append(discharge_levels[i] - model.grid.downstream)

    measured_slope = model.grid.upstream - model.grid.downstream
    normative_slope = np.array(qh(model.grid.discharge.iloc[:, 0].values) - model.grid.downstream)
    relative_slope = measured_slope / normative_slope
    
    relscenslopes = []
    for i, discharge in enumerate(discharges):
        relscenslopes.append((discharge_levels[i] - model.grid.downstream) / normative_slope)
    relscenslopes = np.array(relscenslopes)
    
    # Plot gradient
    minslope = np.min([np.min(relative_slope), relscenslopes.min().min()])
    maxslope = np.max([np.max(relative_slope), relscenslopes.max().max()])
    slope_factors = np.linspace(minslope, maxslope, 25)
    
    vs = np.interp(slope_factors, [minslope, 0, maxslope], [+1, 0, -1])
    vs_matrix = np.array([vs] * 2)
    
    # Plot ever
    # ------------------------------------------------------------------------
    axs[1].pcolormesh([model.grid.time[0], model.grid.time[-1]], 
                      slope_factors, vs_matrix.T, cmap="RdYlGn")

    axs[0].plot(model.grid.time, modelled - model.grid.downstream, '.g')
    axs[0].plot(model.grid.time, model.grid.upstream - model.grid.downstream, '.')
    axs[0].set_ylabel('$\Delta H$')
    for scenslope in scenslopes:
        axs[0].plot(model.grid.time, scenslope, '--k', linewidth=0.5)
        axs[0].text(model.grid.time[-1], scenslope[-1], 'Q={} $m^3/s$'.format(discharge))


    for relscenslope in relscenslopes:
        axs[1].plot(model.grid.time, relscenslope, '--k', linewidth=0.5)
        axs[1].text(model.grid.time[-1], relscenslope[-1], 'Q={} $m^3/s$'.format(discharge))
    
    axs[1].plot(model.grid.time, relative_slope, '.-k', label='measurements')
    axs[1].legend()
    axs[1].set_ylabel('$\Delta H / \Delta H_n$')

    for ax in axs:
        utils.two_axes_style(ax)
    if show:
        plt.show()

def estimate_roughness(model, every=1):
    """
    > run dottermodel each nth timestep
    > optimise manning value until upstream waterlevel matches measurement
    > output calibrated friction and model that uses that friction
    """
    model.logger.info('Optimising model-wide roughness factor')
    model.grid.friction[:] = np.nan
    for i, t in enumerate(tqdm(model.grid.time)):
        if i % every == 0:
            res = optimize.minimize(__objectivefunction, 0.03, args=[model,t],
                                                               tol=0.001,
                                                               method="Nelder-Mead",
                                                               options=dict(maxiter=50))

    # Clip roughness to values higher than 0
    model.grid.friction = model.grid.friction.clip(lower=0)
    # Linearly interpolate between optimized times
    model.grid.friction = model.grid.friction.interpolate(method='time', axis=0)

def blockage_analysis(model):
    """

    """
    green = np.loadtxt('../data/green_2005.csv', skiprows=1, delimiter=',')
    fig, axs = plt.subplots(2, 2, figsize=(10, 6))
    axs = axs.flatten()
    phi = np.linspace(0, 0.9, 100)
    axs[0].plot(model.grid.time, model.grid.friction.iloc[:, 0])

    axs[1].plot(green.T[1], green.T[2], 'dk', label='Green 2005')
    axs[1].plot(phi * 100, __pglinear(phi), '--k', label='pglinear (1991)')
    axs[1].plot(phi * 100, __green(phi), '-k', label='Green (2005)')
    axs[1].plot(phi * 100, __linneman(phi), '-.k', label='Linneman (2017)')
    axs[1].plot([10, 23], [0.03, 0.18], 'dr', label='Linneman (2017)')
    axs[1].legend(loc=2)
    
    model_n = model.grid.friction.iloc[:, 0].values
    axs[3].plot(model.grid.time, __linneman(model_n, reverse=True) * 100, '-', label='linneman')
    axs[3].plot(model.grid.time, __green(model_n, reverse=True) * 100, '-', label='green')
    axs[3].plot(model.grid.time, __pglinear(model_n, reverse=True) * 100, '-', label='pglinear')
    axs[3].set_ylim([0, 100])
    axs[3].legend()

    axs[1].set_ylabel('Blockage factor (%)')
    axs[1].set_xlabel('Blockage factor (%)')
    axs[1].set_ylabel('Manning coefficient ($sm^{1/3}$)')
    axs[0].set_ylabel('Manning coefficient ($sm^{1/3}$)')
    axs[0].set_ylim(axs[1].get_ylim())
    
    for ax in axs:
        utils.two_axes_style(ax)
    axs[2].set_axis_off()
    plt.show()

def __linneman(data, reverse=False):
    if reverse:
        return (data - 0.033) / 0.312
    else: 
        return 0.312 * data + 0.033

def __green(data, reverse=False):
    if reverse:
        return np.log(data / 0.0438) / 2.
    else: 
        return 0.0438 * np.exp(2 * data)

def __pglinear(data, reverse=False):
    if reverse:
        return 1 - 0.033 / data
    else: 
        return 0.033 * (1 - data)**-1
# =============================================================================
# Helper functions
# =============================================================================

def __objectivefunction(n, *args):
    """
    goodness of fit helper function for optimization
    """
    n = n[0]
    model = args[0][0]
    time = args[0][1]

    # Set roughness and run model
    model.grid.friction.loc[time] = n
    model.run(timesteps=[time], progressbar=False)

    modelled = model.output.waterlevel.iloc[:, 0][time]
    measured = model.grid.upstream[time]

    # Goodness of fit function
    gof = np.abs(modelled - measured)
    #sys.stdout.write("n: {}\n measurement: {}\n model: {}\n".format(n, measured, modelled))
    return gof
