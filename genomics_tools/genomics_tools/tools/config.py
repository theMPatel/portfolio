###################################################################
#
# Methods and classes for module configurations
# 
# Author: Milan Patel
# Contact: mpatel5@cdc.gov
# Version 1.0
#
###################################################################

import os
import sys
import json

class Config(object):

    def __init__(self, filepath):
        self.load(filepath)

    def load(self, filepath):

        if isinstance(filepath, basestring):

            if not os.path.exists(filepath):
                raise RuntimeError('Not a valid filepath to load'
                    ' configuration from: {}'.format(filepath))

            with open(filepath, 'r') as f:
                data = json.load(f)

            self._config = Settings(args=data)

        else:

            data = json.load(filepath)
            self._config = Settings(args=data)

    def organism_config(self, organism):
        if organism in self._config:
            return self._config[organism]

        else:
            return None

    def __getitem__(self, setting):

        if setting in self._config:
            return self._config[setting]
        else:
            return None

    def __contains__(self, setting):
        return setting in self._config

class Settings(object):
    # Use this to store the dictionary based
    # settings
    def __init__(self, settingspath=None, args=None):

        self._args = {}

        if args is None and settingspath is None:
            return

        if args is not None:
            self.update(args)

        else:
            self.load(settingspath)

    def load(self, settingspath):

        if not os.path.exists(settingspath):
            raise RuntimeError('Passed invalid filepath'
                ' for requested configuration')

        with open(settingspath, 'r') as f:
            args = json.load(f)

        if not isinstance(args, dict):
            raise TypeError('Loaded args is not a dictionary')

        self._args = args

        for setting_name, setting_value in self._args.iteritems():

            if isinstance(setting_value, dict):

                setattr(self, setting_name, Settings(args=setting_value))

            else:

                setattr(self, setting_name, setting_value)

    def update(self, args):
        self._args.update(args)

        for setting_name, setting_value in self._args.iteritems():

            if isinstance(setting_value, dict):

                setattr(self, setting_name, Settings(args=setting_value))

            else:

                setattr(self, setting_name, setting_value)

    def keys(self):
        return self._args.iterkeys()

    def iteritems(self):
        return self._args.iteritems()

    def __setitem__(self, setting_name, setting_value):
        self._args[setting_name] = setting_value
        setattr(self, setting_name, setting_value)

    def __getitem__(self, setting_name):

        if hasattr(self, setting_name):
            return getattr(self, setting_name)

        elif setting_name in self._args:
            return self._args[setting_name]

        else:
            return None

    def __contains__(self, setting_name):
        return hasattr(self, setting_name) or \
             setting_name in self._args
