###################################################################
#
# Init file for the tools
# 
# Author: Milan Patel
# Contact: mpatel5@cdc.gov
# Version 1.0
#
###################################################################

# We want pdb to be available for all files if the developer is
# running this outside of a production environment.
if __debug__:
    import pdb