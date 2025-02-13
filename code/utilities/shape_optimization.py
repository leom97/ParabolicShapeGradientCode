"""
Implements a class representing the shape optimization problem, with several methods to solve it, and visualize the
results. It is used in shape_optimization_main.py.
"""

# %% Imports

from dolfin import *
import logging

from utilities.meshing import AnnulusMesh, CircleMesh, SquareAnnulusMesh
from utilities.radial_transfer import compute_spherical_transfer_matrix, compute_radial_displacement_matrix, backend_radial_displacement

# %% Setting log and global parameters

parameters["refinement_algorithm"] = "plaza_with_parent_facets"  # must do this to be able to refine stuff associated with a mesh that is being refined
parameters['allow_extrapolation'] = True  # needed if I want a function to be taken from a mesh to a slightly different one

logging.basicConfig(level=logging.INFO,
                    format='%(asctime)s %(levelname)-8s %(message)s',
                    datefmt='%H:%M')

set_log_level(LogLevel.ERROR)
logging.getLogger('FFC').setLevel(logging.ERROR)
logging.getLogger('UFL').setLevel(logging.ERROR)


# %% Class definitions

class ShapeOptimizationProblem:

    def __init__(self):
        # Geometry data
        self.exact_domain = None  # it will be deformed by q_ex
        self.optimization_domain = None
        self.exact_sphere = None
        self.optimization_sphere = None
        self.q_ex = None  # radial function describing the analytical solution to the shape optimization problem
        self.W_ex = None  # associated radial displacement
        self.exact_geometry_dict = None
        self.optimization_geometry_dict = None
        self.simulated_geometry_dict = None
        self.problem_folder = None

        # PDE data: dictionaries and variables that contain info to simulate pdes
        self.exact_pde_dict = None
        self.optimization_pde_dict = None
        self.marker_neumann = None
        self.marker_dirichlet = None
        self.u_N = None  # expression for exterior neumann data
        self.u_D = None  # expression for exterior dirichlet data
        self.T = None
        self.u_ex = None
        self.exact_pde = None  # it will be of HeatEquation class
        self.exact_exterior_BC = None  # can be "N" or "D" (it corresponds to which one is the exact data we prescribe)
        self.pre_assembled_BCs = None  # it contains time-lists of functions, representing the time-dependent boundary conditions at the right instants

        # Function spaces
        self.V_vol = None
        self.V_def = None
        self.V_sph = None

        # Optimization
        self.optimization_dict = None
        self.cost_functional_dict = None
        self.j = None
        self.q_opt = None
        self.opt_results = None
        self.duration = None
        self.v_equation = None
        self.w_equation = None

    def get_domain(self, domain_dict, ground_truth):
        if ground_truth:
            subfolder = "/domain/ground_truth/"
        else:
            subfolder = "/domain/simulation/"
        xdmf_path = None
        if "reload_xdmf" in domain_dict.keys() and ground_truth == True:
            if domain_dict["reload_xdmf"]:
                xdmf_path = self.problem_folder + subfolder
        if domain_dict["type"] == "annulus":
            domain = AnnulusMesh(resolution=domain_dict["resolution"],
                                 path=self.problem_folder + subfolder,
                                 int_refinement=domain_dict["int_refinement"],
                                 ext_refinement=domain_dict["ext_refinement"],
                                 inner_radius=domain_dict["inner_radius"],
                                 outer_radius=domain_dict["outer_radius"],
                                 center=domain_dict["center"],
                                 xdmf_path=xdmf_path)
        elif domain_dict["type"] == "square_annulus":
            domain = SquareAnnulusMesh(resolution=domain_dict["resolution"],
                                       path=self.problem_folder + subfolder,
                                       inner_radius=domain_dict["inner_radius"],
                                       side_length=domain_dict["side_length"],
                                       int_refinement=domain_dict["int_refinement"],
                                       ext_refinement=domain_dict["ext_refinement"],
                                       xdmf_path=xdmf_path
                                       )
        elif domain_dict["type"] == "smoothed_square_annulus":
            domain = SmoothedSquareAnnulusMesh(resolution=domain_dict["resolution"],
                                               path=self.problem_folder + subfolder,
                                               inner_radius=domain_dict["inner_radius"],
                                               side_length=domain_dict["side_length"],
                                               int_refinement=domain_dict["int_refinement"],
                                               ext_refinement=domain_dict["ext_refinement"],
                                               xdmf_path=xdmf_path,
                                               smoothing_radius=domain_dict["smoothing_radius"]
                                               )
        else:
            raise ValueError("Domain unsupported")
        return domain

    def get_sphere(self, sphere_dict, ground_truth):
        if ground_truth:
            subfolder = "ground_truth/"
        else:
            subfolder = "simulation/"
        if sphere_dict["dimension"] == 2:
            sphere = CircleMesh(resolution=sphere_dict["resolution"],
                                path=self.problem_folder + "/sphere/" + subfolder)
        else:
            raise ValueError("Only 2D is supported for now")
        return sphere

    def create_optimal_geometry(self, exact_geometry_dict):

        """
        Create the domain which we will later want to reconstruct.
        """

        logging.info("Creating optimal geometry")

        self.exact_domain = self.get_domain(exact_geometry_dict["domain"], ground_truth=True)
        self.exact_sphere = self.get_sphere(exact_geometry_dict["sphere"], ground_truth=True)
        self.exact_geometry_dict = exact_geometry_dict

        p = self.exact_domain.center  # star shape point
        f_D = self.exact_domain.boundary_radial_function  # radial function to external boundary

        L1_vol = FiniteElement("Lagrange", self.exact_domain.mesh.ufl_cell(), 1)
        V_vol = FunctionSpace(self.exact_domain.mesh, L1_vol)
        L1_sph = FiniteElement("Lagrange", self.exact_sphere.mesh.ufl_cell(), 1)
        V_sph = FunctionSpace(self.exact_sphere.mesh, L1_sph)
        self.V_sph_ex = V_sph
        VD = VectorFunctionSpace(self.exact_domain.mesh, "Lagrange", 1)

        self.q_ex = Function(V_sph)

        circle_coords = self.exact_sphere.mesh.coordinates()[dof_to_vertex_map(V_sph), :]
        q_ex_lambda = eval(exact_geometry_dict["q_ex_lambda"])
        self.q_ex.vector()[:] = q_ex_lambda(circle_coords)

        M = compute_spherical_transfer_matrix(V_vol, V_sph, p=p)
        M2 = compute_radial_displacement_matrix(M, VD, p=p, f_D=f_D)
        self.W_ex = backend_radial_displacement(self.q_ex, M2, VD)

        ALE.move(self.exact_domain.mesh, self.W_ex)  # note, id(domain) = id(self.domain), so, both are moved

        return M2

