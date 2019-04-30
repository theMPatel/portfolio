# Genomics Tools

This set of code represents some of the work that I did as a Fellow at the CDC in the Enteric Diseases Laboratory Branch. The purpose of the tools were to provide the ability to predict certain genetic traits of interest deemed important from the perspective of microbial public health surveillance.

One of the more important aspects of these tools is its ability to predict antibiotic resistance, specifically in foodborne illness causing pathogens (such as E. coli). I have packaged up this particular functionality into a sample for you to run and see the results for yourself.

The reason as to why I have chosen this particular sample to share with you lie partially in nostalgia: It was the first major software product for which I was responsible for the design, implementation and delivery. I have hopefully since come further as an engineer in capability but I felt that this is an accurate representation of the level of complexity I was able to manage coming right out of college.

One more thing, these tools were originally (and still are) run on a high performance compute cluster (HPC) that is on-site at CDC so what is presented here is not completely representative of what the larger project looks like. I'm excited to hear what thoughts you have on this piece of code and I hope using this will be an enjoyable experience!

## Setting things up

I have repackaged the original tools and turned it into a pip installable package (python's package distribution system). Prior to running the install script in this repository, you should ensure a few things:

1. Make sure you have python3 installed, preferably one that is greater than 3.5
    * [Windows](https://www.python.org/downloads/windows/)
        * __NOTE__: If you are on windows, make sure to add the python interpreter to your path.
        * You can usually accomplish by searching control panel, or using the Start search bar to search for "edit user environment variables" or something like that.
    * OSX: Use brew or comparable
        * Python should automatically be added to your path
    * Linux: Use apt or comparable
        * Python should automatically be added to your path

2. Make sure you have Git installed. Some of my dependencies are available only as a git repository
    * [Windows:](https://git-scm.com/downloads)
        * Git bash is nice to have too
    * OSX
        * Use brew or comparable
    * Linux
        * Use apt or comparable

3. It is recommended that you create a python virtual environment (venv) to install these tools into. I would prefer not to muck up your system python with dependencies and installation files from these tools. This is a choice that I leave up to you. Below is the command that has worked for me on all platforms (linux, macos, windows) to create a venv. Make sure that the parent directories for the final venv directory have already been created, otherwise you will get an error. It should leave you with a directory that has a self contained python environment.

    ```bash
    python3 -m venv "${path_to_home}"/.virtualenvs/genomics_tools
    ```

### Using the virtualenv
* Windows
    - If you are using git bash:
        ```bash
        source "${path_to_home}"/.virtualenvs/genomics_tools/Scripts/activate
        ```
    - If you are using the command prompt:
        ```bash
        %USERPROFILE%/.virtualenvs/genomics_tools/Scripts/activate.bat
        ```
* UNIX
    ```bash
    source "${path_to_home}"/.virtualenvs/genomics_tools/bin/activate
    ```
You should note that your console prompt has now been decorated with a (genomics_tools) token [or whichever name you had picked]. This is your proof that the virtualenv was successfully started.

---

At this point, my setup script can pretty much take over fetching the remaining tools that will be needed to make things run. For some transparency, I list below the things that I do which are outside the normal pip install process.

1. I will fetch some binaries from the FTP site of NCBI (National Center for Biotechnology Information) and place them somewhere I can find later when running the antibiotic resistance prediction pipeline.
2. I will fetch sequence data from the FTP site of NCBI so that I have something to run the tools against.
    * __NOTE__: Everything that I have picked are complete genomes meaning they have more or less been circularized. The significance of this is that there will only ever be one contig for which all results fall into. 
3. I will git clone the 'database' which forms the backend of the tools. These are actually just sequence and tabulated text files which contain information regarding which genes confer the corresponding resistances. This database is curated by the Center of Genomic Epidemiology.

### Unit tests

If you would like to see the results of the unit tests that have been setup for this repository you can run the below command.
```bash
./setup.py test
```
If you do not have the appropriate test dependencies, I will install them on the fly.

### Installing the tooling
Assuming you are in the directory where setup.py exists, you can run:
```bash
pip install -e .
```
This will install the tools in editable mode. If everything passes you can run and see the fruits of this endeavor!

### Running the tooling
```bash
genomics_tools --run
```

### Postscript

I have tested these tools on Ubuntu 18.04, macOS Mojave (10.14.4) and Windows 10. If these tools do not work for you, I would love to get some feedback on any errors you encounter so that I can make this more robust. Hope you have enjoyed this exercise!