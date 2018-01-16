#!/usr/bin/env python
# -*- coding: utf-8 -*-

# =============================================================================
# Imports & Function definitions
# =============================================================================
from dotter.models import DotterModel
from dotter import tools

# =============================================================================
# Parameters
# =============================================================================


def test_maaibos():
    """
    Tests whether the numerical approximation tends to the equilibrium
    """
    deltabeek = DotterModel('tests/testcases/trapezoidal/config.ini')

    assert (tools.maaibos(model=deltabeek, discharges=[0.5, 1.0], show=False))
