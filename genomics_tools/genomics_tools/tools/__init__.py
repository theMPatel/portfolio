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

import warnings
import functools

def deprecated(func):
    """
    This is a decorator which can be used to mark
    functions as deprecated. It will result in a warning
    being emitted when the function is used.
    """

    @functools.wraps
    def new_func(*args, **kwargs):
        warnings.simplefilter('always', DeprecationWarning) # Turn off filter
        warnings.warn("Call to deprecated function {}".format(func.__name__),
                        category=DeprecationWarning,
                        stacklevel=2)
        warnings.simplefilter("default", DeprecationWarning)
        return func(*args, **kwargs)

    return new_func