###################################################################
#
# Init file for the genotyping algorithm
# 
# This tool is a sample and distillation of the real application
# hosted at: https://github.com/theMPatel/functional_genomics_tools
#
# Author: Milan Patel
# Contact: https://github.com/theMPatel
# Version 1.0
#
###################################################################

# Human friendly version for the genotyping algorithm
__version__ = '1.0.0'

# We want pdb to be available for all files if the developer is
# running this outside of a production environment.
if __debug__:
    import pdb
