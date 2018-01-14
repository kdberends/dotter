#!/usr/bin/env python
# -*- coding: utf-8 -*-

# =============================================================================
# Imports & Function definitions
# =============================================================================
from dotter.models import DotterModel
import numpy as np 
import matplotlib.pyplot as plt
# =============================================================================
# Parameters
# =============================================================================

def test_trapezoidal_profile():
    """
    test whether the numerical approximation of a trapezoidal profile is
    good compared to analytical results of wet area and hydraulic radius
    """

    deltabeek = DotterModel('tests/testcases/trapezoidal/config.ini')
    
    # The accuracy of the numerical resolution depends on the h_resolution
    deltabeek.grid.h_resolution = 100
    deltabeek.grid.generate_grid()
    waterlevels = np.linspace(0, 4, 1000)

    # Analytical solution
    A_anal = 3 * waterlevels + waterlevels ** 2 / np.tan(30 * np.pi / 180.)
    P_anal = 3 + 2 * waterlevels / np.sin(30 * np.pi / 180.)
    R_anal = A_anal / P_anal

    error = deltabeek.grid.wet_area[-1](waterlevels) - A_anal
    assert (np.max(np.abs(error)) < 0.002)
    
    error = deltabeek.grid.hydraulic_radius[-1](waterlevels) - R_anal
    assert (np.max(np.abs(error)) < 0.0005)
    
