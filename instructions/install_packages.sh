#!/bin/bash

conda install -c conda-forge mamba
mamba install --yes -c conda-forge fenics
mamba install -c conda-forge gmsh meshio
pip install pygmsh
pip install scipy
pip install ipython
pip install tqdm
