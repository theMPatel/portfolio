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

import contextlib
from distutils import log
import ftplib
import functools
import gzip
import hashlib
import io
import itertools
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile
import traceback

def install_package(package):
    """
    Infrastructure to install packages at runtime in case they
    are not present in the system libraries.
    """

    if not isinstance(package, str):
        raise TypeError('String required for package installation')

    if not package:
        raise ValueError('Recieved empty string for package installation')

    # This will fail if you are on windows and do not have pip installed
    # which can happen if you skirted your sys admin and installed
    # python yourself from a portable python version
    try:
        import pip
    except ImportError as e:
        message = "I need pip to install setuptools, please " \
            "install pip for {} before continuing".format(sys.executable)

        raise ImportError(message) from e

    if hasattr(pip, 'main'):
        pip.main(['install', package])
    else:
        import pip._internal
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
    from setuptools.command.test import test as TestCommand

def prompt_continue(message=None):

    if message is None:
        message = 'Continue (Y/N)?: '

    answer = input(message).lower()

    while answer != 'y' and answer != 'n':
        answer = input("Please enter (Y/N): ").lower()

    return answer == 'y'

@contextlib.contextmanager
def chdir_and_return(directory):
    """
    If you need to run a particular command in specific
    directory but would like for the program to return
    to its original working directory, you can use
    this function to automatically maintain that state

    :param directory: The directory to return to
    """
    old = os.getcwd()

    if not os.path.exists(directory) or not \
        os.path.isdir(directory):
        raise ValueError("Invalid directory: {}".format(directory))

    try:
        os.chdir(directory)
        yield
    finally:
        os.chdir(old)

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

NCBI_ENDPOINT = "ftp.ncbi.nlm.nih.gov"

_tools_dirs = [
    os.path.normpath(os.path.expanduser("~/.tools")),
    os.path.normpath(os.path.expanduser("~/.bin")),
    os.path.normpath(os.path.expanduser("~/.tmp"))
]

_valid_directory_names = { 
        "NCBI" :'ncbi_blast',
        "DB" : 'pointfinder_db',
        "SEQ" : 'sequence_data'
        }

def get_discriminator_string_by_platform(platform):
    """
    This is specifically for the NCBI tools since
    they have the OS a particular package was intended
    for in the filename.
    """

    if "linux" in platform:
        return "linux"

    elif "darwin" in platform:
        return "macosx"

    elif "win" in platform:
        return "win"

    raise ValueError("I don't know what platform this is: {}"
        .format(platform))

def populate_syspath(directories):
    """
    This function will populate the sys path with the directories
    and sub directories from a list of directories.
    """

    path_parts = []
    # We can use shutil.which to get the correct binary we want
    # if we ensure that the path is correct. Unfortunately, I
    # didn't feel safe manipulating your path (esp if you're on
    # windows) so am resorting to this.
    for directory in directories:
        for root, dirs, files in os.walk(directory):
            path_parts.append(root)

    os.environ["PATH"] += os.pathsep + os.pathsep.join(path_parts)

def parse_md5_data(hash_data):
    """
    Usually md5 data usually consists of the hash followed by the
    filename, like so:

    a982lk3j09aosilejk3 some_file_with_integrity

    This function returns back the hash from the data
    """

    if not hash_data:
        raise ValueError("MD5 data cannot be empty")

    hash_info = hash_data.split()

    if len(hash_info) != 2:
        raise ValueError("Invalid hash data")

    return hash_info[0]

def parse_md5_file(path_to_file):
    """
    Reads in a md5 file and returns back the hash data

    :param path_to_file: The path to the file you want to
        read in.
    """

    if not os.path.exists(path_to_file):
        raise OSError("Path to MD5 file does not exist!")

    with open(path_to_file, 'r') as f:
        return parse_md5_data(f.read())

def verify_download_md5(file, hash_data):
    """
    Verifies a that a file matches the hash_data provided.

    :param file: The path to the file to be loaded
    :param hash_data: The actual hash data the file should match
    """
    log.info("Verifying: " + file)
    
    if not os.path.exists(file):
        raise OSError("Path {} does not exist".format(file))

    md5_digest = hashlib.md5()
    with open(file, 'rb') as f:

        for block in chunked_file_reader(f):
            md5_digest.update(block)

    return md5_digest.hexdigest() == hash_data

def untar_contents_to_dir(source_tar, tarfile_final_path):
    """
    Untars the _contents_ of a tar file to a directory. This means
    that if your tar file looks like:
    some_folder
    |__ dir_1
    |__ dir_2
    |__ dir_3

    The tarfile_final_path will not contain some_folder, but rather
    the contents of some_folder.

    :param source_tar: The source tar file
    :param tarfile_final_path: The directory to unpack to
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

def gunzip_file(source, destination):
    """
    Unzip from a source path to a destination path

    :param source: The source gzip file
    :param destination: The destination zip file
    """

    if not hasattr(source, 'read'):
        source = gzip.open(source, 'rb')

    if not hasattr(destination, 'write'):
        destination = open(destination, 'wb')

    shutil.copyfileobj(source, destination)

    source.close()
    destination.close()

def git_checkout(branch, stdout=None, stderr=None):
    """
    Checks out a particular branch in a git repo. It assumes
    that the current working directory of this script
    is in a git repo.

    :param branch: Assumes that this branch exists in your git repo
        either locally or remotely.
    :raises `subprocess.CalledProcessError`: When the subprocess fails
    :raises `ValueError` If branch is empty or not a string
    """

    if not branch or not isinstance(branch, str):
        raise ValueError("Invalid branch argument")

    args = ["git", "checkout", branch]
    subprocess.Popen(args, stdout=stdout, stderr=stderr)

def git_reset_hard_head(stdout=None, stderr=None):
    """
    Assumes that the current working directory of the script is
    in the git repo you want. Resets the local working directory
    to known state

    :param directory: The directory to do the reset
    """
    args = ['git', 'reset', '--hard', 'HEAD']
    subprocess.Popen(args, stdout=stdout, stderr=stderr)

def git_pull(branch, stdout=None, stderr=None):
    """
    Assumes that the current working directory of the script is
    in the git repo you want. Also assumes that we don't need
    to authenticate to run this.

    :param branch: The branch that you want checked out to latest
    :raises `subprocess.CalledProcessError`: When the subprocess fails
    :raises `ValueError` If branch is empty or not a string
    """
    if not branch or not isinstance(branch, str):
        raise ValueError("Invalid branch argument")

    git_checkout(branch, stdout=stdout, stderr=stderr)
    args = ["git", "pull"]
    subprocess.Popen(args, stdout=stdout, stderr=stderr)

BLAST_DIR = "blast/executables/blast+/LATEST/"
def fetch_ncbi_tools(tools_dir, retry_count=0):
    """
    Fetches the ncbi tools by logging into their FTP site
    anonymously, changing directories to the latest tools path,
    determining the correct package based on OS and downloading
    both the tar and the md5. This is all assuming you're using
    a x64 OS (sorry if that's not the case)

    :param tarfile_final_path: The final resting place for
        these tools

    :param retry_count: The current retry count in case the download
        from NCBI fails
    """
    actual_tools_needed = ['makeblastdb', 'blastn']
    
    if os.name == 'nt':
        for i in range(len(actual_tools_needed)):
            actual_tools_needed[i] += '.exe'
    
    if all(shutil.which(tool) for tool in actual_tools_needed):
        log.info("Looks like we have the blast tools we need " \
            "already!")

        return

    log.info("Needed blast tools not found, proceeding with fetch")

    if retry_count >= 3:
        raise RuntimeError("Tried to download the NCBI tools three "
            "times. Now I'm giving up sorry!")

    tarfile_final_path = os.path.join(tools_dir,
        _valid_directory_names.get("NCBI", "ncbi_blast"))
    os.makedirs(tarfile_final_path, exist_ok=True)
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

        search_string = get_discriminator_string_by_platform(sys.platform)
        package = next(file for file in available_versions \
                        if search_string in file)

        tar_file_path = os.path.join(tempdir, package)
        files = [package, package+".md5"]

        for i, file in enumerate(files):
            path = os.path.join(tempdir, file)

            with open(path, "wb") as handle:
                downloader = functools.partial(download_and_write,
                    file_handle=handle)
                
                log.info("Beginning to download: {}"
                    .format(file))

                ftp_conn.retrbinary("RETR "+file, downloader)

            files[i] = path

    local_file = files[0]
    local_file_md5 = files[1]
    md5_data = parse_md5_file(local_file_md5)

    if not verify_download_md5(local_file, md5_data):
        log.error("Download failed md5 verification, will try a total of "
            " 3 times before quitting.")

        # Back off factor of 2
        time.sleep(2*(retry_count+1))
        return fetch_ncbi_tools(tarfile_final_path, retry_count=retry_count+1)

    untar_contents_to_dir(tar_file_path, tarfile_final_path)
    tempdir_obj.cleanup()

_default_pointfinder_branch = 'master'
_default_commit_hash = "3f07218144693c977f976e6e26a3d68d99011a34"
def fetch_pointfinder_db(tools_dir):
    """
    Will fetch the point finder database that we use to BLAST
    a query genome against

    :param tools_dir: The tools directory where you want this installed
    """

    if not shutil.which("git"):
    # You don't have git installed, what?!?!
        raise RuntimeError("Please install git somewhere on your path " \
            "before continuing!")

    log.info("Fetching sequence database for mutation search..")
    database = "https://bitbucket.org/genomicepidemiology/" \
        "pointfinder_db.git"

    directory_name = _valid_directory_names.get("DB", "pointfinder_db")
    directory_path = os.path.join(tools_dir, directory_name)
    git_directory = os.path.join(directory_path, ".git")

    args = ["git", "clone", database, directory_path]
    to_clone = not (os.path.exists(directory_path) and os.path.exists(git_directory))

    out = subprocess.DEVNULL
    err = out

    if to_clone:
        subprocess.Popen(args, stdout=out, stderr=err)

    with chdir_and_return(directory_path):
        git_reset_hard_head(stdout=out, stderr=err)
        git_pull(_default_pointfinder_branch, stdout=out, stderr=err)
        git_checkout(_default_commit_hash, stdout=out, stderr=err)

    log.info("Successfully retrieved database!")

def fetch_sequence_data(tools_dir):
    """
    Fetches assemblies from the NCBI ftp site
    GCA_000299455.1_ASM29945v1_genomic.fna
    """
    
    log.info("Fetching demonstration sequences..")
    directory_name = _valid_directory_names.get("SEQ", "sequence_data")

    # (Root_dir, Filename)
    root = "/"
    to_retrieve = [
        ("genomes/all/GCF/000/299/455/GCF_000299455.1_ASM29945v1",
            "GCF_000299455.1_ASM29945v1_genomic.fna.gz"),
        ("/genomes/all/GCA/003/691/425/GCA_003691425.1_ASM369142v1",
            "GCA_003691425.1_ASM369142v1_genomic.fna.gz"),
        ("/genomes/all/GCF/004/368/015/GCF_004368015.1_ASM436801v1",
            "GCF_004368015.1_ASM436801v1_genomic.fna.gz")
    ]

    sequence_final_path = os.path.join(tools_dir,
        _valid_directory_names.get("SEQ", "sequence_data"))
    os.makedirs(sequence_final_path, exist_ok=True)

    tempdir_obj = tempfile.TemporaryDirectory()
    tempdir = tempdir_obj.name

    with ftplib.FTP(NCBI_ENDPOINT) as ftp_conn:
        ftp_conn.login()

        for (remote_path, filename) in to_retrieve:

            dest_path = os.path.join(sequence_final_path, filename).replace(
                ".gz", "")

            if os.path.exists(dest_path):
                log.info("Already have: {}".format(dest_path))
                continue

            ftp_conn.cwd(remote_path)
            local_file = os.path.join(tempdir, filename)

            with open(local_file, "wb") as handle:
                downloader = functools.partial(download_and_write,
                    file_handle=handle)

                log.info("Beginning to download: {}"
                    .format(filename))

                ftp_conn.retrbinary("RETR "+filename, downloader)

    for downloaded_file in os.listdir(tempdir):
        source_path = os.path.join(tempdir, downloaded_file)
        dest_path = os.path.join(sequence_final_path, downloaded_file)

        if ".gz" in dest_path:
            dest_path = dest_path.replace(".gz", "")
        try:
            gunzip_file(source_path, dest_path)
        except:
            pass

    tempdir_obj.cleanup()

def retrieve_necessary_deps():
    """
    Retrieves the necessary tools from the remote endpoints and
    selects an ideal place on the filesystem to place them.

    :raises ValueError: If I can't figure out a good place to put
        these tools
    """
    populate_syspath(_tools_dirs)

    chosen_tools_dir = ''
    for path in _tools_dirs:
        log.info("Checking for pre-installed tools: {}".format(path))

        if not os.path.exists(path):
            continue

        # I'm going to check to see if you already have any of
        # the tools that I need or a portion of them
        paths = list(
            map(
                lambda tail: os.path.join(path, tail), 
                _valid_directory_names.values()
                )
            )

        if any(os.path.exists(p) for p in paths):
            chosen_tools_dir = path
            break

    else:
        # I couldn't find anything, so I'm going to pick
        # whichever directory either doesn't exist and if they
        # all exist whichever has the least amount of things
        # in it
        counts = []
        for path in _tools_dirs:
            if not os.path.exists(path):
                chosen_tools_dir = path
                break

            counts.append([len(os.listdir(path)), path])

        else:
            counts.sort(reverse=True)
            chosen_tools_dir = counts.pop()[-1]

    if not chosen_tools_dir:
        raise ValueError("Unable to select an ideal spot for my extra tools!")

    fetch_ncbi_tools(chosen_tools_dir)
    fetch_sequence_data(chosen_tools_dir)
    fetch_pointfinder_db(chosen_tools_dir)

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
        log.info("Fetching custom dependencies")
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
        log.info("Fetching custom dependencies")
        retrieve_necessary_deps()

class PyTest(TestCommand):
    user_options = [('pytest-args=', 'a', "Arguments to pass into py.test")]

    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.pytest_args = ""

    def run_tests(self):
        import shlex

        try:
            import pytest
        except ImportError:
            install_package("pytest")
            import pytest

        errno = pytest.main(shlex.split(self.pytest_args))
        sys.exit(errno)

class RemoveExtraTools(develop):
    """
    Class that will help to cleanly uninstall anything that I have
    downloaded to your computer
    """

    def run(self):
        """
        Will go through all the places where I could have installed things
        and removes the directories
        """

        for root, tail in itertools.product(_tools_dirs,
            _valid_directory_names.values()):
            to_rm = os.path.join(root, tail)
            if os.path.exists(to_rm):
                log.info("Removing: {}".format(to_rm))

                try:
                    if os.path.isdir(to_rm):
                        shutil.rmtree(to_rm)
                    else:
                        os.unlink(to_rm)
                except OSError:
                    log.info("Failed to remove: {}".format(to_rm))

__version__ = '0.0.1-dev'
__author__ = 'Milan Patel'
__title__ = 'genomics_tools'

packages = [
            'genomics_tools',
            'genomics_tools.tools',
            'genomics_tools.genotyping'
]

# Absolutely contained, no external dependencies.
requires = [
]

test_requires = [
    "pytest>=4.4.1"
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
            "install" : PostInstallCommand,
            "test" : PyTest,
            "uninstall" : RemoveExtraTools
        },
        tests_requires=test_requires,
        entry_points={
            'console_scripts': [
                'genomics_tools=genomics_tools.__main__:_main'
                ]
            }
        )