#!/usr/bin/env python3
###################################################################
#
# Here, I'll provide some of the logic needed to make sure that
# you have the appropriate tools and that they are on the path
# to use!
# 
# This tool is a sample and distillation of the real application
# hosted at: https://github.com/theMPatel/functional_genomics_tools
#
# Author: Milan Patel
# Contact: https://github.com/theMPatel
# Version 1.0
#
###################################################################

from setuptools import setup
from setuptools.command.develop import develop
from setuptools.command.install import install

import ftplib
import functools
import hashlib
import os
import sys
import tempfile

from tools.environment import find_executable_path

def prompt_continue(message=None):

    if message is None:
        message = 'Continue (Y/N)?: '

    answer = input(message).lower()

    while answer != 'y' and answer != 'n':
        answer = input("Please enter (Y/N): ").lower()

    return answer == 'y'

def download_and_write(data, file_handle=None):

    if file_handle is None:
        raise RuntimeError("No file handle to write to!")

    file_handle.write(data)

required_tools = [
    "makeblastdb",
    "blastn"
]

NCBI_ENDPOINT = "ftp.ncbi.nlm.nih.gov"
BLAST_DIR = "blast/executables/blast+/LATEST/"

fetch_tools = False
if not all(find_executable_path(tool) for tool in required_tools):
    message = "It looks like I need to download some tools to get" \
        " things working, do I have permission to do that? (Y/N)"

    if not prompt_continue(message=message):
        print("Maybe next time, see ya!")
        sys.exit(0)

    fetch_tools = True

if fetch_tools:
    with ftplib.FTP(NCBI_ENDPOINT) as ftp_conn:

        ftp_conn.login()
        ftp_conn.cwd(BLAST_DIR)
        
        # The latests tools are available as tars which is nice since
        # I won't be able to invoke your system installer from here.
        available_versions = [file for file in ftp_conn.nlst() if \
                                file.endswith("tar.gz") and "src" not in file]

        # Determine the OS we're on. Fingers crossed that you're using something
        # relatively new!
        search_string = ''
        if "linux" in platform:
            search_string = "linux"

        elif "darwin" in platform:
            search_string = "macosx"

        elif "win" in platform:
            search_string = "win"

        package = next(file for file in available_versions \
                        if search_string in file)

        tempdir = tempfile.TemporaryDirectory()
        local_file_path = os.path.join(tempdir.name, package)
        md5_path = os.path.join(tempdir.name, package+".md5")

        tar_handle = open(local_file_path, "wb")
        md5_handle = open(md5_path, "wb")

        tar_downloader = partial(download_and_write, file_handle=tar_handle)
        md5_downloader = partial(download_and_write, file_handle=md5_handle)

        print("Beginning to download package: {} to {}".format(
            os.path.join(NCBI_ENDPOINT, BLAST_DIR, package), local_file_path))

        ftp_conn.retrbinary("RETR "+package, tar_downloader)
        ftp_conn.retrbinary("RETR "+package+".md5", md5_downloader)


class PostDevelopCommand(develop):
    """
    Post installation tasks when installed in development
    mode
    """

    def run(self):
        """
        Calls the parent class' implementation of this function
        before continuing on to install the NCBI blast tools
        """
        super().run()

__version__ = '0.0.1-dev'
__author__ = 'Milan Patel'
__title__ = 'genomics_tools'

packages = [
            'genomics_tools',
]

requires = [
    'requests>=2.21.0',
    'six>=1.12.0',
    'tqdm>=4.29.0',
    'urllib3>=1.24.1',
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