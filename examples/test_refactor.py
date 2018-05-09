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

# =============================================================================
# Parameters
# =============================================================================

deltabeek = DotterModel('cases/grotebeek/config.ini')

print (deltabeek.output._asdict().keys())


    print (j)
#tools.estimate_roughness(deltabeek, every=15)
#deltabeek.run()
#deltabeek.dash(dashtype=2, show=True)
#tools.maaibos(deltabeek, discharges=[0.25, 0.5, 1.0, 2.0], show=True)
