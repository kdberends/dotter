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
from copy import deepcopy, copy
# =============================================================================
# Parameters
# =============================================================================
configfile = 'cases/vegetation/config.ini'


deltabeek = DotterModel(configfile)
#tools.estimate_roughness(model, every=15)

tools.maaibos(model=deltabeek, discharges=[1, 2], show=True, configfile=configfile)
