###################################################################
#
# A reorganization of tools specifically for command line related
# things. Everything should eventually get migrated here.
# 
# Author: Milan Patel
# Contact: mpatel5@cdc.gov
# Version 1.0
#
###################################################################

import os
import sys
import subprocess
import threading

def add_cmdline_args(executable_name, cmd_args, args):
    """
    Modifies a list used for commandline args in place while
    doing some basic validations
    """

    if not isinstance(cmd_args, list):
        raise TypeError('Commandline args destination '
                        'must be of type list, got'
                        ' {} instead'.format(type(cmd_args))
                        )

    if not isinstance(args, (tuple, list)):
        raise TypeError('Args to add must be in list or tuple format'
                        ' got {} instead'.format(type(args))
                        )
    for arg in args:

        if not isinstance(arg, basestring) or not \
            arg.startswith('-'):

            log_warning('Invalid argument provided for {}:'
                        ' {}, skipping..'.format(executable_name, arg))
            continue

        cmd_args.append(arg)

def add_cmdline_kwargs(executable_name, cmd_args, kwargs):
    """
    Modifies a list used for commandline args in place while
    doing some basic validations
    """

    if not isinstance(cmd_args, list):
        raise TypeError('Commandline args destination '
                        'must be of type list, got'
                        ' {} instead'.format(type(cmd_args))
                        )

    if not isinstance(kwargs, dict):
        raise TypeError('Key word args to add must be in dict format'
                        ' got {} instead'.format(type(kwargs))
                        )

    for arg, value in kwargs.iteritems():
        if not isinstance(arg, basestring) or not \
            arg.startswith('-'):

            log_warning('Invalid argument provided for {}:'
                        ' {} -. {}, skipping..'.format(
                                                    executable_name, 
                                                    arg,
                                                    value
                                                )
                        )

            continue

        cmd_args.extend([str(arg), str(value)])

def popen(args, stdout=None, stderr=None, cwd=None, shell=False):

    if not isinstance(args, list):
        raise RuntimeError('Provided arguments must be of type list')

    if not stderr:
        stderr = subprocess.PIPE

    if not stdout:
        stdout = subprocess.PIPE

    child = subprocess.Popen(
                            args,
                            stdout=stdout,
                            stderr=stderr,
                            cwd=cwd,
                            shell=shell
                        )

    out, err = child.communicate()

    return child.returncode, out, err