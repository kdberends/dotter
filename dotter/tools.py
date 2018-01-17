"""
Dottermodel tools
"""

# =============================================================================
# Imports
# =============================================================================
from copy import copy, deepcopy
from . import utils
from . import models
from scipy import optimize
import numpy as np
from tqdm import tqdm
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
# =============================================================================
# Tools
# =============================================================================

def maaibos(model, discharges=None, critical_friction=0.15, show=False, configfile=None):
    """
    predictive use (at different Q)
    """

    # generate QH at 'critical friction'
    # -------------------------------------------------------------------------
    discharge_restore = copy(model.grid.discharge)
    friction_restore = copy(model.grid.friction)

    model.logger.info('Generating QH for normative friction')
    QHdischarges = np.linspace(0, model.grid.discharge.max().max(), 20)

    model.grid.friction[:] = critical_friction

    upstream = list()
    for discharge in tqdm(QHdischarges):
        model.grid.quickset_discharge(discharge)
        model.run(timesteps=[model.grid.time[0]], progressbar=False)
        upstream.append(model.output.waterlevel.iloc[0, 0])

    qh = interp1d(QHdischarges, upstream, bounds_error=False, fill_value=0)

    # Restore previous conditions
    # -------------------------------------------------------------------------
    model.grid.discharge = copy(discharge_restore)
    model.grid.friction = copy(friction_restore)

    model.run()
    modelled_upstream = model.output.waterlevel.iloc[:, 0]

    if configfile is not None:
        calmodel = models.DotterModel(configfile)
        calibrated_roughness = estimate_roughness(calmodel, every=10)



    # Restore previous conditions
    # -------------------------------------------------------------------------
    model.grid.discharge = copy(discharge_restore)
    model.grid.friction = copy(friction_restore)

    # Get measurements and scenario results
    # -------------------------------------------------------------------------
    fig = plt.figure(figsize=(10,4))
    axs = list()
    axs.append(fig.add_subplot(211))
    axs.append(fig.add_subplot(223))
    axs.append(fig.add_subplot(224))


    slope_factors = np.linspace(-1, 1, 20)
    vs = np.interp(slope_factors, [-1, 0, 1], [-1, 0, +1])

    vs_matrix = np.array([vs] * 2)

    # Plot ever
    # ------------------------------------------------------------------------
    axs[0].pcolormesh([model.grid.time[0], model.grid.time[-1]],
                      slope_factors, vs_matrix.T, cmap="RdYlGn")
    axs[0].plot(model.grid.time, [0] * len(model.grid.time), '-k')
    axs[0].plot(model.grid.time, qh(model.grid.discharge.iloc[:, 0]) - model.grid.upstream, '--k')
    axs[0].plot(model.grid.time, qh(model.grid.discharge.iloc[:, 0]) - modelled_upstream, '-c')

    axs[0].set_ylabel('$\Delta H$')

    q = np.linspace(0, model.grid.discharge.max().max(), 100)
    axs[1].plot(q, qh(q), '-r', label='Critical QH')
    axs[1].plot(model.grid.discharge.iloc[:, 0], model.grid.upstream, '.k')

    axs[2].plot(model.grid.time, model.grid.friction.iloc[:, 0], '-c')
    if configfile is not None:
        axs[2].plot(model.grid.time, calibrated_roughness.iloc[:, 0], '--k')

    axs[2].plot(model.grid.time, [critical_friction] * len(model.grid.time), '--r')
    for ax in axs:
        utils.two_axes_style(ax)
    if show:
        plt.show()

    return True

def estimate_roughness(model, every=1):
    """
    > run dottermodel each nth timestep
    > optimise manning value until upstream waterlevel matches measurement
    > output calibrated friction and model that uses that friction
    """
    model.logger.info('Optimising model-wide roughness factor')
    model.grid.friction[:] = np.nan
    for i, t in enumerate(tqdm(model.grid.time)):
        guess = 0.03
        if i % every == 0:
            res = optimize.minimize(__objectivefunction, guess, args=[model,t],
                                                               tol=0.001,
                                                               method="Nelder-Mead",
                                                               options=dict(maxiter=50))
            guess = res.x[0]

    # Clip roughness to values higher than 0
    model.grid.friction = model.grid.friction.clip(lower=0)
    # Linearly interpolate between optimized times
    model.grid.friction = model.grid.friction.interpolate(method='time', axis=0)

    return model.grid.friction

def blockage_analysis(model):
    """

    """

    fig, axs = plt.subplots(2, 2, figsize=(10, 6))
    axs = axs.flatten()
    phi = np.linspace(0, 0.9, 100)
    axs[0].plot(model.grid.time, model.grid.friction.iloc[:, 0])


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
    model.run(timesteps=[time], progressbar=False, output_to_file=False)

    modelled = model.output.waterlevel.iloc[:, 0][time]
    measured = model.grid.upstream[time]

    # Goodness of fit function
    gof = np.abs(modelled - measured)
    #sys.stdout.write("n: {}\n measurement: {}\n model: {}\n".format(n, measured, modelled))
    return gof
