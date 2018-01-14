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


def test_backwater():
    """
    Tests whether the numerical approximation tends to the equilibrium
    """
    deltabeek = DotterModel('tests/testcases/backwater/config.ini')

    # Equilibrium depth
    depth = 1.568138

    # Check whether above depth is equilibrium
    h = deltabeek.grid.bedlevel[0] + depth
    A = deltabeek.grid.wet_area[0](h)
    R = deltabeek.grid.hydraulic_radius[0](h)
    i = deltabeek.grid.bedslope[0]
    C = R ** (1 / 6.) / 0.04

    assert(np.abs(3.50 - A * C * np.sqrt(R * i)) < 0.001)

    deltabeek.run(timesteps=[deltabeek.grid.time[0]])
    error = np.abs(deltabeek.output.waterdepth[0][0] - depth)

    assert (error < 0.001)

test_backwater()
