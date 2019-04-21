#!/usr/bin/env python3

from setuptools import setup

__version__ = '0.0.1-dev'
__author__ = 'Milan Patel'
__title__ = 'genomics_tools'

packages = [
            'genomics_tools',

]

requires = [

]

setup(
        name=__title__,
        version=__version__,
        packages=packages,
        author=__author__,
        python_requires=">=3.5",
        install_requires=requires,
        entry_points={
            'console_scripts': [
                'genomics_tools=genomics_tools.__main__:_main'
                ]
            }
        )