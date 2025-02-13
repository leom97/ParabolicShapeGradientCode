"""
This script is to numerically verify the estimate of Theorem 5.1.
It serves to generate the experiments of Section 6.

The order of convergence (OOC) is expected to asymptotically approach 2,
since we couple space and time discretization like dt = h^2.

Head to ./configuration.py to set up the data needed for this script
(in configuration.py, instructions are provided to reproduce Table 1).

Then, run the script normally through e.g. the command line.
"""

# %% Imports

from dolfin import *
import logging
import numpy as np
from tqdm import tqdm
import os
import sys

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '../')))
from utilities.radial_transfer import backend_radial_displacement
from utilities.shape_optimization import ShapeOptimizationProblem
from utilities.ooc_verification import get_spiky_radial_function, get_assembled_shape_gradient, W1i_norm, \
    string_to_vector_field
from shape_gradients_ooc.configuration import _f, _g

os.chdir(os.path.dirname(os.path.abspath(__file__)))

# %% Setting log and global parameters

parameters['allow_extrapolation'] = True  # I want a function to be taken from a mesh to a slightly different one
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%H:%M')
set_log_level(LogLevel.ERROR)
runs_path = "./"  # path to the directory containing this .py file, try to change to absolute path in case of problems
mesh_path = "../../../mesh_data"
from pathlib import Path
Path(mesh_path).mkdir(parents=True, exist_ok=True)

# %% Set-ups

# Please head to configuration.py to tweak settings and data

try:
    import importlib.util

    spec = importlib.util.spec_from_file_location("configuration", runs_path + "configuration.py")
    pd = importlib.util.module_from_spec(spec)
    sys.modules["configuration"] = pd
    spec.loader.exec_module(pd)

    geometry_dict = pd.geometry_dict
    pde_dict = pd.pde_dict
    cost_functional_dict = pd.cost_functional_dict
    experiment_dict = pd.experiment_dict
    smooth_displacements_dict = pd.smooth_displacements_dict
    h_tentative = pd.h_tentative
except:
    raise Exception("Couldn't load configuration file from path")

# the data for problem 2.1.1
f = _f(t=0)
g = _g(t=0)

if pde_dict["ode_scheme"] == "crank_nicolson":
    dt_power = 1
else:
    dt_power = 2

dt_multiplier = experiment_dict["dt_multiplier"]

h_actual = np.array([])
dt_actual = np.array([])

# Displacement parameters
amp = [0.1, 0.2]
spikes = range(10)

results = []

# %% Run

for h, k in zip(h_tentative, range(len(h_tentative))):

    logging.info("#######################################")
    logging.info(f"{(k + 1)}/{len(h_tentative)}")
    logging.info("#######################################")

    geometry_dict["domain"]["resolution"] = h
    problem = ShapeOptimizationProblem()  # let's create a problem class for easily holding various useful variables
    problem.problem_folder = mesh_path
    M2 = problem.create_optimal_geometry(geometry_dict)

    h_actual = np.append(h_actual, problem.exact_domain.mesh.hmax())

    V_vol = FunctionSpace(problem.exact_domain.mesh, "CG", 1)  # scalar, linear FEM on the volume
    V_sph = FunctionSpace(problem.exact_sphere.mesh, "CG", 1)  # scalar, linear FEM on the unit sphere
    V_def = VectorFunctionSpace(problem.exact_domain.mesh, "CG", 1)  # vector (2D), linear FEM on the volume
    Q = VectorFunctionSpace(problem.exact_domain.mesh, 'DG', 0)  # space for gradients of V_def functions

    pde_dict["N_steps"] = int(np.ceil(dt_multiplier * pde_dict["T"] / (h_actual[-1] ** (dt_power))))
    dt_actual = np.append(dt_actual, pde_dict["T"] / pde_dict["N_steps"])

    dj = get_assembled_shape_gradient(V_sph, V_def, V_vol, M2, problem.exact_domain, pde_dict,
                                      cost_functional_dict, f, g)
    evaluations = {"dj": [], "norms": []}

    logging.info("Evaluating the gradient")

    # First, test dj with non-smooth fields.
    with tqdm(total=len(amp) * len(spikes)) as pbar:
        for A in amp:
            for s in spikes:
                dq = get_spiky_radial_function(problem, V_sph, A, s)

                W = backend_radial_displacement(dq, M2, V_def)
                dj_dq = np.dot(dj[:], W.vector()[:])

                evaluations["dj"].append(dj_dq)
                evaluations["norms"].append(W1i_norm(W, Q, problem, V_vol))

                pbar.update(1)

    # Then, test dj by smooth displacement fields.
    for i in tqdm(range(len(smooth_displacements_dict["x"]))):
        W = string_to_vector_field(smooth_displacements_dict["x"][i], smooth_displacements_dict["y"][i], problem)

        dj_dq = np.dot(dj[:], W.vector()[:])

        evaluations["dj"].append(dj_dq)
        evaluations["norms"].append(W1i_norm(W, Q, problem, V_vol))

        pbar.update(1)

    evaluations["dj"] = np.array(evaluations["dj"])
    evaluations["norms"] = np.array(evaluations["norms"])

    results.append(evaluations)

# %% Post-processing

# Let's create a matrix of gradient evaluations, and of gradient norms
dj = []
norms = []

for e in results:
    dj.append(e["dj"])
    norms.append(e["norms"])

dj = np.array(dj)
norms = np.array(norms)
dual_errors = np.max(np.abs(dj - dj[-1, :]) / norms[-1], axis=1)
ooc = np.log(dual_errors[1:] / dual_errors[:-1]) / np.log(h_actual[1:] / h_actual[:-1])

logging.info(f"OOCs are {ooc}")

# %% Gather results (to save, e.g. by means of pickles)

results_dict = {
    "pde_dict": pde_dict,
    "geometry_dict": geometry_dict,
    "experiment_dict": experiment_dict,
    "dj": dj,
    "norms": norms,
    "dual_errors": dual_errors,
    "ooc": ooc
}

import pickle
from datetime import datetime
filename = f"results/data_{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.pkl"
with open(filename, "wb") as file:
    pickle.dump(results_dict, file)
