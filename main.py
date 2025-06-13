import os
from utils.simulation import run_simulation
from utils.plot_utils import plot_parametric_comparison
import pandas as pd
import shutil


# 确保输出目录存在
def setup_directories(base_dir="output_advanced"):
    if os.path.exists(base_dir):
        shutil.rmtree(base_dir)
    os.makedirs(base_dir)
    return base_dir


if __name__ == "__main__":
    base_dir = setup_directories()

    parametric_studies = {
        "Flexible Fixator": (35.0e9, "Flexible"),
        "Standard Fixator": (70.0e9, "Standard"),
        "Rigid Fixator": (140.0e9, "Rigid"),
    }

    all_results = {}
    for name, params in parametric_studies.items():
        print(f"--- Running Simulation for: {name} ---")
        sim_output_dir = os.path.join(base_dir, params[1])
        os.makedirs(sim_output_dir)
        os.makedirs(os.path.join(sim_output_dir, "stress_images"))

        results = run_simulation(params, sim_output_dir)
        all_results[name] = results

    print("\n--- Generating Parametric Comparison Plot ---")
    plot_parametric_comparison(all_results, base_dir)

    # 保存最终结果
    summary_rows = []
    for label, res in all_results.items():
        summary_rows.append(
            {
                "group": label,
                "final_stress_fixator": res["avg_stress_fixator"][-1],
                "final_stress_bone": res["avg_stress_bone"][-1],
                "final_stress_callus": res["avg_stress_callus"][-1],
                "final_gap_strain": res["gap_strain"][-1],
            }
        )
    pd.DataFrame(summary_rows).to_csv(
        os.path.join(base_dir, "summary_final_results.csv"), index=False
    )

    print(f"\nSimulation completed. All results saved to '{base_dir}' directory.")
