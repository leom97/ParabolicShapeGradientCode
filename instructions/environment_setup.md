# Instructions for environment setup

The following instructions are tested on Ubuntu 22.04 and Python 3.9.
We used the package manager Anaconda throughout, where a pre-packaged version of FEniCS is available.

Perform the following steps to set up the environment:

- change the current directory to the directory of this file `./`
- create a virtual environment: `conda create --name article_env python==3.9`, and take note of the installation location
- activate the newly created environment: `conda activate article_env`
- run `source install_packages.sh` to install the required packages. This will require some interaction with the command line, and some time
- run `python -m site` from within the virtual environment to obtain the path to the python executable
- run `source apply_patches.sh` to suitably patch the already installed packages (*): this requires the knowledge of the location of the conda environment executable obtained at the previous point

(*) This will apply the commit https://bitbucket.org/fenics-project/ffc/pull-requests/98 to the file tools.py.

You can move on to `../code/getting_started.md` for an overview of the code and some indications on how to run it.

<!-- Maybe be now more explicit in th end of the instructions -->
<!-- also, instruct about the usage of the PDE run instructions! -->