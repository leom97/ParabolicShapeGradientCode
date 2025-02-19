import numpy as np
from dolfin import *

"""

Configuration file for shape_gradient_ooc_verification.py.
Default values are provided, they yield first row of Table 1 in the article.

To obtain row 2, change cost_functional_dict["final_smoothing_lambda"] = "lambda t: 0*t + 1.
And for row 3, start from the default values and replace the geometry_dict["domain"] value
with the one commented below it.

"""

experiment_dict = {
    "N_it": 6,  # number of refinements for the spatial mesh
    "N_steps_multiplier": 15,  # used to couple dt and h: dt = 1/N_steps_multiplier * h^(1 or 2), 2 is for implicit Euler
    "h_multiplier": 0.9    # multiplies h_tentative = 2^{-i}, i = 0, ..., N_it - 1
}

# This dictionary describes the domain on which the shape gradient is computed (domain),
# and on which the spherical functions live (sphere).
# Comment the "annulus" domain and replace it with the "square_annulus"
# to reproduce the 3rd row of Table 1.
geometry_dict = {
    "domain": {"type": "annulus", "resolution": None, "ext_refinement": 1.0, "int_refinement": 1.0, "inner_radius": 1,
               "outer_radius": 2,
               "center": np.array([0, 0]), "reload_xdmf": False},
    # "domain": {"type": "square_annulus", "resolution": None, "ext_refinement": 1.0, "int_refinement": 1.0, "inner_radius": 2,
    #            "outer_radius": 2, "side_length": 1,
    #            "center": np.array([0, 0]), "reload_xdmf": False},
    "sphere": {"dimension": 2, "resolution": 0.5},
    "q_ex_lambda": 'lambda x: 0'
}

# The resolution of the domain mesh is dictated by the following vector: every entry correspond to a mesh width,
# which should get progressively finer
h_tentative = experiment_dict["h_multiplier"] * 1 / (2 ** np.arange(0, experiment_dict["N_it"]))

# Change "ode_scheme" to "implicit_euler" or "crank_nicolson".
# In the paper, only "implicit_euler" is used.
pde_dict = {
    "ode_scheme": "implicit_euler",
    "marker_dirichlet": [3],  # interior of annulus
    "marker_neumann": [2],  # exterior of annulus
    "T": 1,  # final time of simulation
    "N_steps": None
}

# Describes the function \eta inside of the cost function.
# To obtain the 2nd row of Table 1, set "final_smoothing_lambda" to "lambda t: 0*t + 1".
cost_functional_dict = {
    "final_smoothing_lambda": "lambda t: exp(-0.005/pow(t,2)) if t > DOLFIN_EPS else 0.0",
    "discretization": None,
    "H1_smoothing": 0
}

# Smooth_displacements_dict["x"][i] and smooth_displacements_dict["y"][i] represents a smooth vector field.
smooth_displacements_dict = {
    "x": ["x^3", ".5*x", ".5*y^3", ".2*y", ".2*x+.2*y^2", ".2*x+.2*y^2", ".2*x*y+.2*y*2", "-.5*x-.2*x*y", ".5*x+.4*x*y",
          "x"],
    "y": [".5*y", "y^3", ".2*x", ".5*x^3", ".5*x^3-.2*x", ".5*x^3-.2*x*y", ".5*x^3-.2*x*y", "-.4*x^2*y+.4*y^2",
          "-.3*y^3-.1*x*y", "y"]
}


# Pde data: used to set up the states v, w
class _f(UserExpression):  # Dirichlet datum
    def __init__(self, t, **kwargs):
        super().__init__(self, **kwargs)
        self.t = t

    def eval(self, values, x):
        values[0] = self.t ** 3 * np.cos(5 * self.t) * (x[0] ** 2 + x[1])

    def value_shape(self):
        return ()


class _g(UserExpression):  # Neumann datum
    def __init__(self, t, **kwargs):
        super().__init__(self, **kwargs)
        self.t = t

    def eval(self, values, x):
        values[0] = self.t ** 2 * np.sin(5 * self.t) * np.sin(x[0])

    def value_shape(self):
        return ()