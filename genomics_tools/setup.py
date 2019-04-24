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

def install_package(package):
    """
    Infrastructure to install packages at runtime in case they
    are not present in the system libraries.
    """

    if not isinstance(package, str):
        raise TypeError('String required for package installation')

    if not package:
        raise ValueError('Recieved empty string for package installation')

    import pip

    if hasattr(pip, 'main'):
        pip.main(['install', package])
    else:
        pip._internal.main(['install', package])

try:
    # There's a chance you don't have this installed
    import setuptools
except ImportError:
    # If we fail here, not much else I can do
    install_package("setuptools")

finally:
    from setuptools import setup
    from setuptools.command.develop import develop
    from setuptools.command.install import install

from distutils import log
import ftplib
import functools
import hashlib
import io
import os
import shutil
import sys
import tarfile
import tempfile
import traceback


def prompt_continue(message=None):

    if message is None:
        message = 'Continue (Y/N)?: '

    answer = input(message).lower()

    while answer != 'y' and answer != 'n':
        answer = input("Please enter (Y/N): ").lower()

    return answer == 'y'

def chunked_file_reader(file_obj, chunk_size=io.DEFAULT_BUFFER_SIZE):
    while True:
        data = file_obj.read(chunk_size)

        if not data:
            break

        yield data

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

_tools_dirs = [
    os.path.expanduser("~/.tools"),
    os.path.expanduser("~/.bin"),
    os.path.expanduser("~/.tmp")
]

final_ncbi_foldername = "ncbi_blast"

def get_discriminator_string_by_platform():
    """
    This is specifically for the NCBI tools since
    they have the OS a particular package was intended
    for in the filename
    """

    if "linux" in sys.platform:
        return "linux"

    elif "darwin" in sys.platform:
        return "macosx"

    elif "win" in sys.platform:
        return "win"

    raise ValueError("I don't know what platform this is: {}"
        .format(sys.platform))

def verify_download_md5(downloaded_asset):
    """
    Validates the downloaded asset. Assumes that the
    downloaded file has another file with an extra
    md5 extension that contains the digest information.
    """
    log.info("Verifying: " + downloaded_asset)
    md5_file = downloaded_asset+".md5"
    
    if not os.path.exists(md5_file) or \
        not os.path.exists(downloaded_asset):
        return False

    with open(md5_file, "r") as f:
        digest = f.read().strip()
        digest = digest.split()[0]

    md5_digest = hashlib.md5()
    with open(downloaded_asset, 'rb') as f:

        for block in chunked_file_reader(f):
            md5_digest.update(block)

    return md5_digest.hexdigest() == digest

def untar_to_directory(source_tar, tarfile_final_path):
    """
    Untars the _contents_ of a tar file to a directory, in
    the case that you already have a directory made. This
    functionality is similar to the -C flag you can pass
    to tar on the command-line
    """

    if not os.path.exists(source_tar):
        raise OSError("{} does not exist!".format(tarfile_final_path))

    with tarfile.open(source_tar, mode="r:gz") as tfile:
        
        for member in tfile.getmembers():

            # Drop the root directory. The package looks like
            # ncbi-blast-version/
            #   - File
            #   - File
            #   - Dir
            #
            # I want the files in my target directory and not also
            # the original root directory. It's easier to find
            # the bin directory if I know exactly how it resolves.
            member.name = '/'.join(member.name.split('/')[1:])
            tfile.extract(member, path=tarfile_final_path)

def fetch_ncbi_tools(tarfile_final_path, retry_count=0):
    """
    Fetches the ncbi tools by logging into their FTP site
    anonymously, changing directories to the latest tools path,
    determining the correct package based on OS and downloading
    both the tar and the md5. This is all assuming you're using
    a x64 OS (sorry if that's not the case)
    """
    if retry_count >= 3:
        raise RuntimeError("Tried to download the NCBI tools three "
            "times. Now I'm giving up sorry!")

    os.makedirs(os.path.dirname(tarfile_final_path), exist_ok=True)
    tempdir_obj = tempfile.TemporaryDirectory()
    tempdir = tempdir_obj.name

    log.info("Beginning to download NCBI tools from: " + NCBI_ENDPOINT)
    with ftplib.FTP(NCBI_ENDPOINT) as ftp_conn:

        ftp_conn.login()
        ftp_conn.cwd(BLAST_DIR)

        # The latest tools are available as tars which is nice since
        # I won't be able to invoke your system installer from here.
        available_versions = [file for file in ftp_conn.nlst() if \
                                file.endswith("tar.gz") and "src" not in file]

        search_string = get_discriminator_string_by_platform()

        package = next(file for file in available_versions \
                        if search_string in file)

        tar_file_path = os.path.join(tempdir, package)
        
        mapping = {
                    package: os.path.join(tempdir, package),
                    package+".md5": os.path.join(tempdir, package+'.md5')
                }

        for file, path in mapping.items():

            with open(path, "wb") as handle:
                downloader = functools.partial(download_and_write, file_handle=handle)
                
                log.info("Beginning to download: {}"
                    .format(file))

                ftp_conn.retrbinary("RETR "+file, downloader)


    if not verify_download_md5(mapping[package]):
        log.error("Download failed md5 verification, will try a total of "
            " 3 times before quitting.")

        time.sleep(2)
        return fetch_ncbi_tools(tarfile_final_path, retry_count=retry_count+1)

    untar_to_directory(tar_file_path, tarfile_final_path)
    tempdir_obj.cleanup()

def retrieve_necessary_deps():
    """
    Retrieves the appropriate tools from the remote endpoints
    if they are needed. As of this writing, it is not the case
    that we need anything more than the NCBI blast tools.
    If that were the case, I would refactor this method
    to rather choose a tools directory, iterate over a
    hashmap of tool name to function objects passing in
    the tools directory to the function and letting it
    decide how it wants to decide whether a tool has been
    already been downloaded.
    """

    if all(shutil.which(tool) for tool in required_tools):
        return

    chosen_tools_dir = ''
    for path in _tools_dirs:

        ncbi_path = os.path.join(path, final_ncbi_foldername)

        if not os.path.exists(path):
            os.makedirs(path, exist_ok=True)
            chosen_tools_dir = ncbi_path
            break

        # Either you've run this before or wow you use NCBI!
        if os.path.exists(ncbi_path):
            bin_path = os.path.join(ncbi_path, 'bin')

            if os.path.exists(bin_path):
                log.info("Looks like you have some NCBI tools already!")
                return

            # It's quite possible that the last install was not valid,
            # This is a clean-up step
            try:
                shutil.rmtree(ncbi_path)
            
            except OSError as e:
                log.error("Failed to remove folder: {}".format(ncbi_path))
                traceback.print_exc()
                continue
            # Since we've elected to use this one before, let's go ahead
            # and pick it again.
            chosen_tools_dir = ncbi_path
            break

        chosen_tools_dir = ncbi_path
        break

    if not chosen_tools_dir:
        raise RuntimeError("I have no idea what would be a good"
            " place to put the ncbi tools since all the good ones"
            " seem to be taken:\n {}".format("\n".join(_tools_dirs)))

    fetch_ncbi_tools(chosen_tools_dir)

class PostDevelopCommand(develop):
    """
    Post installation tasks when installed in development
    mode i.e. `pip install -e .`
    """

    def run(self):
        """
        Calls the parent class' implementation of this function
        before continuing on to install any non-pip dependencies.
        """
        develop.run(self)
        retrieve_necessary_deps()

class PostInstallCommand(install):
    """
    Post installation tasks when installed in normal
    install mode i.e. `pip install {package}`
    """

    def run(self):
        """
        Calls the parent class' implementation of this function
        before continuing on to install any non-pip dependencies.
        """
        install.run(self)
        retrieve_necessary_deps()

__version__ = '0.0.1-dev'
__author__ = 'Milan Patel'
__title__ = 'genomics_tools'

packages = [
            'genomics_tools',
            'genomics_tools.tools',
            'genomics_tools.genotyping'
]

requires = [
    'requests>=2.21.0',
    'six>=1.12.0',
    'tqdm>=4.29.0',
]

setup(
        name=__title__,
        version=__version__,
        packages=packages,
        author=__author__,
        python_requires=">=3.5",
        install_requires=requires,
        cmdclass={
            "develop" : PostDevelopCommand,
            "install" : PostInstallCommand
        },
        entry_points={
            'console_scripts': [
                'genomics_tools=genomics_tools.__main__:_main'
                ]
            }
        )