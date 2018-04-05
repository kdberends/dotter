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
configfile = 'cases/grotebeek_maaibos/config.ini'


deltabeek = DotterModel(configfile)
#tools.estimate_roughness(model, every=15)

deltabeek.run()
#deltabeek.dash(dashtype=1, show=False)
deltabeek.dash(dashtype=2, show=False)

tools.maaibos(model=deltabeek, 
              critical_friction=0.10, 
              show=True, 
              every=5,
              configfile=configfile)
