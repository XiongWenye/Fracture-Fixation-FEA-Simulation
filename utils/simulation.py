import numpy as np
import os
import pandas as pd

from utils.fea_core import create_mesh, assemble_global_stiffness
from utils.materials import get_material_type
from utils.analysis import (
    calculate_stress,
    calculate_average_stresses,
    calculate_fracture_gap_strain,
)
from utils.plot_utils import plot_analysis_results, plot_stress_model


def apply_bc_and_loads(nodes, K, F, applied_load, fixator_length):
    """
    Applies boundary conditions and loads to the global stiffness matrix and force vector.
    """

    left_nodes = np.where(np.isclose(nodes[:, 0], 0.0))[0]
    for node in left_nodes:

        for dof in range(2):
            idx = 2 * node + dof

            K[idx, :] = 0
            K[:, idx] = 0
            K[idx, idx] = 1.0
            F[idx] = 0.0

    right_nodes = np.where(np.isclose(nodes[:, 0], fixator_length))[0]
    if len(right_nodes) > 0:
        load_per_node = applied_load / len(right_nodes)
        for node in right_nodes:
            F[2 * node] = load_per_node

    return K, F


def run_simulation(simulation_params, output_dir):
    """
    Runs a finite element simulation for a bone fracture model with a fixator.
    """

    (E_fixator, label) = simulation_params
    E_bone, E_callus_initial, E_callus_final = 18e9, 0.1e9, 18e9
    nu = 0.3
    bone_length, bone_width, fracture_width = 0.2, 0.02, 0.01
    fixator_thickness, fixator_length = 0.005, bone_length + 0.04
    applied_load = 1000
    num_steps = 20
    nx, ny = 40, 12

    fracture_params = (fixator_length, bone_length, fracture_width)

    nodes, elements = create_mesh(fixator_length, bone_width, fixator_thickness, nx, ny)

    material_props = {
        "bone": {"E": E_bone, "nu": nu},
        "callus": {"E": E_callus_initial, "nu": nu},
        "fixator": {"E": E_fixator, "nu": nu},
    }

    results = {
        "avg_stress_fixator": [],
        "avg_stress_bone": [],
        "avg_stress_callus": [],
        "gap_strain": [],
    }

    for step in range(num_steps):
        print(f"  - Running step {step + 1}/{num_steps}")

        t = step / (num_steps - 1)
        material_props["callus"]["E"] = E_callus_initial + t * (
            E_callus_final - E_callus_initial
        )

        K = assemble_global_stiffness(
            nodes,
            elements,
            material_props,
            bone_width,
            fixator_thickness,
            fracture_params,
        )
        F = np.zeros(2 * nodes.shape[0])
        K, F = apply_bc_and_loads(nodes, K, F, applied_load, fixator_length)
        U = np.linalg.solve(K, F)
        stresses = calculate_stress(
            nodes,
            elements,
            U,
            material_props,
            bone_width,
            fixator_thickness,
            fracture_params,
        )

        plot_stress_model(
            nodes, elements, material_props, step, U, stresses, output_dir
        )

        avg_stresses = calculate_average_stresses(
            nodes, elements, stresses, bone_width, fixator_thickness, fracture_params
        )
        gap_strain = calculate_fracture_gap_strain(
            nodes, U, fixator_length, bone_length, fracture_width
        )

        results["avg_stress_fixator"].append(avg_stresses["fixator"])
        results["avg_stress_bone"].append(avg_stresses["bone"])
        results["avg_stress_callus"].append(avg_stresses["callus"])
        results["gap_strain"].append(gap_strain)

    plot_analysis_results(results, output_dir)

    df = pd.DataFrame(
        {
            "step": list(range(len(results["gap_strain"]))),
            "avg_stress_fixator": results["avg_stress_fixator"],
            "avg_stress_bone": results["avg_stress_bone"],
            "avg_stress_callus": results["avg_stress_callus"],
            "gap_strain": results["gap_strain"],
        }
    )
    df.to_csv(os.path.join(output_dir, f"simulation_results.csv"), index=False)
    return results
