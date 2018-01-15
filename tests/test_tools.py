#!/usr/bin/env python
# -*- coding: utf-8 -*-

# =============================================================================
# Imports & Function definitions
# =============================================================================
import matplotlib as mpl
if os.environ.get('DISPLAY','') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')
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
