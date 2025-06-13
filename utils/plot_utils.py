import os
import matplotlib.pyplot as plt
import numpy as np


def plot_analysis_results(results, output_dir):
    """绘制应力遮挡和间隙应变的演化图"""
    steps = range(len(results["gap_strain"]))

    # 绘制应力遮挡图
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(steps, results["avg_stress_fixator"], "b-o", label="Fixator ")
    ax.plot(steps, results["avg_stress_bone"], "k-^", label="Bone ")
    ax.plot(steps, results["avg_stress_callus"], "g-s", label="Callus ")
    ax.set_xlabel("Simulation Step ")
    ax.set_ylabel("Average Von Mises Stress (Pa) ")
    ax.set_title("Stress Shielding Effect ")
    ax.legend()
    ax.grid(True)
    plt.savefig(os.path.join(output_dir, "stress_shielding.png"))
    plt.close()

    # 绘制骨折间隙应变图
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(steps, np.array(results["gap_strain"]) * 100, "r-o")  # 转换为百分比
    ax.set_xlabel("Simulation Step ")
    ax.set_ylabel("Fracture Gap Strain (%) ")
    ax.set_title("Fracture Gap Strain Evolution ")
    ax.grid(True)
    plt.savefig(os.path.join(output_dir, "gap_strain.png"))
    plt.close()


def plot_parametric_comparison(all_results, output_dir="output_advanced"):
    """绘制参数化研究的对比图"""
    fig, axes = plt.subplots(1, 2, figsize=(18, 7))

    colors = plt.cm.viridis(np.linspace(0, 1, len(all_results)))

    # 对比最终应力分布
    ax = axes[0]
    for i, (label, results) in enumerate(all_results.items()):
        final_stresses = [
            results["avg_stress_fixator"][-1],
            results["avg_stress_bone"][-1],
            results["avg_stress_callus"][-1],
        ]
        ax.bar(
            np.arange(3) + i * 0.2,
            final_stresses,
            width=0.2,
            label=label,
            color=colors[i],
        )
    ax.set_xticks(np.arange(3) + 0.2)
    ax.set_xticklabels(["Fixator", "Bone", "Callus"])
    ax.set_ylabel("Final Average Stress (Pa)")
    ax.set_title("Final Stress Distribution vs. Fixator Stiffness")
    ax.legend()

    # 对比间隙应变演化
    ax = axes[1]
    for i, (label, results) in enumerate(all_results.items()):
        ax.plot(results["gap_strain"], "-o", label=label, color=colors[i])
    ax.set_xlabel("Simulation Step")
    ax.set_ylabel("Fracture Gap Strain")
    ax.set_title("Gap Strain Evolution vs. Fixator Stiffness")
    ax.legend()
    ax.grid(True)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "parametric_comparison.png"))
    plt.close()
