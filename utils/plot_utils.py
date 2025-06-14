import os
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata

from utils.analysis import calculate_von_mises


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


def plot_stress_model(
    nodes,
    elements,
    material_props,
    step,
    U=None,
    stresses=None,
    output_dir="output_advanced",
):
    # 静态变量用于收集每2步的图像
    if not hasattr(plot_stress_model, "collected_steps"):
        plot_stress_model.collected_steps = []
        plot_stress_model.collected_imgs = []

    if stresses is not None:
        # 计算von Mises应力
        vm_stress = calculate_von_mises(stresses)

        # 创建每个节点的应力值 (取相邻单元的平均)
        node_stresses = np.zeros(nodes.shape[0])
        node_counts = np.zeros(nodes.shape[0])
        for elem_idx, elem in enumerate(elements):
            for node in elem:
                node_stresses[node] += vm_stress[elem_idx]
                node_counts[node] += 1
        node_stresses /= node_counts

        x = nodes[:, 0]
        y = nodes[:, 1]
        z = node_stresses

        xi = np.linspace(x.min(), x.max(), 100)
        yi = np.linspace(y.min(), y.max(), 100)
        zi = griddata((x, y), z, (xi[None, :], yi[:, None]), method="cubic")

        # 只收集步数为2的倍数的图像
        if step % 2 == 0:
            plot_stress_model.collected_steps.append(step)
            plot_stress_model.collected_imgs.append((xi, yi, zi, z.min(), z.max()))

        # 仍然单独保存每步的应力云图（可选）
        fig, ax = plt.subplots(figsize=(12, 6))
        levels = np.linspace(z.min(), z.max(), 20)
        cs = ax.contourf(xi, yi, zi, levels=levels, cmap="jet", extend="both")
        cbar = fig.colorbar(cs)
        cbar.set_label("Von Mises Stress (Pa)")
        ax.set_title(f"Stress Distribution - Step {step}")
        ax.set_xlabel("Length (m)")
        ax.set_ylabel("Height (m)")
        plt.savefig(
            os.path.join(output_dir, "stress_images", f"stress_step_{step:02d}.png")
        )
        plt.close()

    # 在最后一次调用时，汇总并绘制所有收集到的图像
    if hasattr(plot_stress_model, "collected_steps") and step == 19:
        n = len(plot_stress_model.collected_imgs)
        ncols = 5
        nrows = 2
        fig, axes = plt.subplots(
            nrows, ncols, figsize=(6 * ncols, 5 * nrows), sharey=True
        )
        axes = axes.flatten()
        # 统一颜色范围
        vmin = min([img[3] for img in plot_stress_model.collected_imgs])
        vmax = max([img[4] for img in plot_stress_model.collected_imgs])
        for i, (xi, yi, zi, _, _) in enumerate(plot_stress_model.collected_imgs):
            ax = axes[i]
            levels = np.linspace(vmin, vmax, 20)
            cs = ax.contourf(
                xi,
                yi,
                zi,
                levels=levels,
                cmap="jet",
                extend="both",
                vmin=vmin,
                vmax=vmax,
            )
            ax.set_title(f"Step {plot_stress_model.collected_steps[i]}")
            ax.set_xlabel("Length (m)")
            if i % ncols == 0:
                ax.set_ylabel("Height (m)")
        # 隐藏多余的子图
        for j in range(i + 1, nrows * ncols):
            axes[j].axis("off")
        # 添加一个统一的颜色条
        fig.subplots_adjust(right=0.88)
        cbar_ax = fig.add_axes([0.90, 0.15, 0.02, 0.7])
        fig.colorbar(cs, cax=cbar_ax, label="Von Mises Stress (Pa)")
        plt.suptitle("Stress Distribution at Steps Multiple of 2")
        plt.savefig(os.path.join(output_dir, "stress_images", "stress_summary.png"))
        plt.close()
        # 清理静态变量
        del plot_stress_model.collected_steps
        del plot_stress_model.collected_imgs
