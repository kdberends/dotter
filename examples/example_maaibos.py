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

# =============================================================================
# Parameters
# =============================================================================

deltabeek = DotterModel('cases/grotebeek/config.ini')

deltabeek.run()
