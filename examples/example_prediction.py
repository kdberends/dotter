# Dotter model prototype
#
# This file: high-level run file
#
# Author: Koen Berends
# Contact: k.d.berends@utwente.nl & koen.berends@utwente.nl
# Copyright (c) 2017 Deltares
# Copyright (c) 2017 University of Twente


# =============================================================================
# Imports & Function definitions
# =============================================================================
from dotter.models import DotterModel
from dotter import tools
import matplotlib.pyplot as plt
# =============================================================================
# Parameters
# =============================================================================

deltabeek = DotterModel('cases/grotebeek/config.ini')
tools.estimate_roughness(deltabeek, every=2)
deltabeek.run()
deltabeek.dash(dashtype=2, show=True)
