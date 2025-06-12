import numpy as np
import matplotlib.pyplot as plt
import os
import shutil
import pandas as pd


# 平面应力弹性矩阵
def plane_stress_matrix(E, nu):
    return E / (1 - nu**2) * np.array([[1, nu, 0], [nu, 1, 0], [0, 0, (1 - nu) / 2]])


# 计算单元刚度矩阵
def element_stiffness(nodes, element, D):
    # 4节点四边形等参元
    xi = np.array([-1 / np.sqrt(3), 1 / np.sqrt(3), 1 / np.sqrt(3), -1 / np.sqrt(3)])
    eta = np.array([-1 / np.sqrt(3), -1 / np.sqrt(3), 1 / np.sqrt(3), 1 / np.sqrt(3)])
    weights = np.array([1, 1, 1, 1])

    ke = np.zeros((8, 8))

    for gp in range(4):
        # 形函数导数
        dN_dxi = 0.25 * np.array(
            [-(1 - eta[gp]), (1 - eta[gp]), (1 + eta[gp]), -(1 + eta[gp])]
        )
        dN_deta = 0.25 * np.array(
            [-(1 - xi[gp]), -(1 + xi[gp]), (1 + xi[gp]), (1 - xi[gp])]
        )

        # 雅可比矩阵
        J = np.zeros((2, 2))
        for i in range(4):
            node_idx = element[i]
            x, y = nodes[node_idx, :]
            J[0, 0] += dN_dxi[i] * x
            J[0, 1] += dN_dxi[i] * y
            J[1, 0] += dN_deta[i] * x
            J[1, 1] += dN_deta[i] * y

        detJ = np.linalg.det(J)
        invJ = np.linalg.inv(J)

        # 形函数对x,y的导数
        dN_dx = invJ[0, 0] * dN_dxi + invJ[0, 1] * dN_deta
        dN_dy = invJ[1, 0] * dN_dxi + invJ[1, 1] * dN_deta

        # B矩阵
        B = np.zeros((3, 8))
        for i in range(4):
            B[0, 2 * i] = dN_dx[i]
            B[1, 2 * i + 1] = dN_dy[i]
            B[2, 2 * i] = dN_dy[i]
            B[2, 2 * i + 1] = dN_dx[i]

        ke += B.T @ D @ B * detJ * weights[gp]

    return ke


# 计算von Mises应力
def calculate_von_mises(stress):
    sxx = stress[:, 0]
    syy = stress[:, 1]
    txy = stress[:, 2]
    return np.sqrt(sxx**2 + syy**2 - sxx * syy + 3 * txy**2)


# 确保输出目录存在
def setup_directories(base_dir="output_advanced"):
    if os.path.exists(base_dir):
        shutil.rmtree(base_dir)
    os.makedirs(base_dir)
    return base_dir


# 创建网格
def create_mesh(fixator_length, bone_width, fixator_thickness, nx, ny):
    x = np.linspace(0, fixator_length, nx + 1)
    y = np.linspace(0, bone_width + 2 * fixator_thickness, ny + 1)
    nodes = np.zeros(((nx + 1) * (ny + 1), 2))
    for i in range(nx + 1):
        for j in range(ny + 1):
            nodes[i * (ny + 1) + j, 0] = x[i]
            nodes[i * (ny + 1) + j, 1] = y[j]
    elements = np.zeros((nx * ny, 4), dtype=int)
    for i in range(nx):
        for j in range(ny):
            n1 = i * (ny + 1) + j
            n2 = (i + 1) * (ny + 1) + j
            n3 = (i + 1) * (ny + 1) + j + 1
            n4 = i * (ny + 1) + j + 1
            elements[i * ny + j, :] = [n1, n2, n3, n4]
    return nodes, elements


# 确定材料类型
def get_material_type(x, y, bone_width, fixator_thickness, fracture_params):
    fixator_length, bone_length, fracture_width = fracture_params
    if y < fixator_thickness or y > bone_width + fixator_thickness:
        return "fixator"
    fracture_start = (
        (fixator_length - bone_length) / 2 + bone_length / 2 - fracture_width / 2
    )
    fracture_end = fracture_start + fracture_width
    if fracture_start <= x <= fracture_end:
        return "callus"
    else:
        return "bone"


# 平面应力弹性矩阵
def plane_stress_matrix(E, nu):
    return E / (1 - nu**2) * np.array([[1, nu, 0], [nu, 1, 0], [0, 0, (1 - nu) / 2]])


# 计算单元刚度矩阵
def element_stiffness(nodes, element, D):
    xi_vals = np.array([-1, 1, 1, -1]) / np.sqrt(3)
    eta_vals = np.array([-1, -1, 1, 1]) / np.sqrt(3)
    weights = np.array([1, 1, 1, 1])
    ke = np.zeros((8, 8))
    for gp in range(4):
        xi, eta = xi_vals[gp], eta_vals[gp]
        dN_dxi = 0.25 * np.array([-(1 - eta), (1 - eta), (1 + eta), -(1 + eta)])
        dN_deta = 0.25 * np.array([-(1 - xi), -(1 + xi), (1 + xi), (1 - xi)])
        elem_nodes = nodes[element, :]
        J = np.dot(np.array([dN_dxi, dN_deta]), elem_nodes)
        detJ = np.linalg.det(J)
        invJ = np.linalg.inv(J)
        dN_dxy = np.dot(invJ, np.array([dN_dxi, dN_deta]))
        B = np.zeros((3, 8))
        for i in range(4):
            B[0, 2 * i] = dN_dxy[0, i]
            B[1, 2 * i + 1] = dN_dxy[1, i]
            B[2, 2 * i] = dN_dxy[1, i]
            B[2, 2 * i + 1] = dN_dxy[0, i]
        ke += B.T @ D @ B * detJ * weights[gp]
    return ke


def assemble_global_stiffness(
    nodes, elements, material_props, bone_width, fixator_thickness, fracture_params
):
    """
    组装全局刚度矩阵 (已修正版本)
    """
    num_nodes = nodes.shape[0]
    K = np.zeros((2 * num_nodes, 2 * num_nodes))

    for elem_idx, elem in enumerate(elements):
        # 获取单元中心坐标以确定材料类型
        x_center = np.mean(nodes[elem, 0])
        y_center = np.mean(nodes[elem, 1])
        material_type = get_material_type(
            x_center, y_center, bone_width, fixator_thickness, fracture_params
        )

        # 获取材料属性
        E = material_props[material_type]["E"]
        nu = material_props[material_type]["nu"]
        D = plane_stress_matrix(E, nu)

        # 计算单元刚度矩阵
        ke = element_stiffness(nodes, elem, D)

        # ### 正确的组装逻辑 ###
        # 将8x8的单元刚度矩阵ke组装到全局刚度矩阵K中
        for i in range(4):  # 局部节点i
            for j in range(4):  # 局部节点j
                for di in range(2):  # 自由度 (0 for x, 1 for y)
                    for dj in range(2):
                        # 全局行索引
                        row = 2 * elem[i] + di
                        # 全局列索引
                        col = 2 * elem[j] + dj
                        # 局部索引
                        local_row = 2 * i + di
                        local_col = 2 * j + dj
                        # 累加到全局矩阵
                        K[row, col] += ke[local_row, local_col]
    return K


def apply_bc_and_loads(nodes, K, F, applied_load, fixator_length):
    """
    施加边界条件和载荷 (使用 np.isclose 提高稳健性)
    """
    # 使用isclose来处理浮点数比较
    left_nodes = np.where(np.isclose(nodes[:, 0], 0.0))[0]
    for node in left_nodes:
        # 固定x和y方向的自由度
        for dof in range(2):
            idx = 2 * node + dof
            # 罚函数法施加边界条件
            K[idx, :] = 0
            K[:, idx] = 0
            K[idx, idx] = 1.0
            F[idx] = 0.0

    # 在右侧施加拉力
    right_nodes = np.where(np.isclose(nodes[:, 0], fixator_length))[0]
    if len(right_nodes) > 0:
        load_per_node = applied_load / len(right_nodes)
        for node in right_nodes:
            F[2 * node] = load_per_node  # x方向拉力

    return K, F


def calculate_stress(
    nodes, elements, U, material_props, bone_width, fixator_thickness, fracture_params
):
    element_stresses = np.zeros((len(elements), 3))
    for elem_idx, elem in enumerate(elements):
        x_center = np.mean(nodes[elem, 0])
        y_center = np.mean(nodes[elem, 1])
        material_type = get_material_type(
            x_center, y_center, bone_width, fixator_thickness, fracture_params
        )
        E, nu = material_props[material_type]["E"], material_props[material_type]["nu"]
        D = plane_stress_matrix(E, nu)
        Ue = U[
            np.array([2 * n for n in elem] + [2 * n + 1 for n in elem])
            .reshape(2, 4)
            .T.flatten()
        ]

        xi, eta = 0, 0
        dN_dxi = 0.25 * np.array([-(1 - eta), (1 - eta), (1 + eta), -(1 + eta)])
        dN_deta = 0.25 * np.array([-(1 - xi), -(1 + xi), (1 + xi), (1 - xi)])
        elem_nodes = nodes[elem, :]
        J = np.dot(np.array([dN_dxi, dN_deta]), elem_nodes)
        invJ = np.linalg.inv(J)
        dN_dxy = np.dot(invJ, np.array([dN_dxi, dN_deta]))
        B = np.zeros((3, 8))
        for i in range(4):
            B[0, 2 * i] = dN_dxy[0, i]
            B[1, 2 * i + 1] = dN_dxy[1, i]
            B[2, 2 * i] = dN_dxy[1, i]
            B[2, 2 * i + 1] = dN_dxy[0, i]

        element_stresses[elem_idx, :] = D @ B @ Ue
    return element_stresses


# --- 新增分析函数 ---


def calculate_average_stresses(
    nodes, elements, stresses, bone_width, fixator_thickness, fracture_params
):
    """按材料类型计算平均Von Mises应力"""
    vm_stresses = calculate_von_mises(stresses)

    avg_stresses = {"bone": 0, "callus": 0, "fixator": 0}
    counts = {"bone": 0, "callus": 0, "fixator": 0}

    for i, elem in enumerate(elements):
        x_center = np.mean(nodes[elem, 0])
        y_center = np.mean(nodes[elem, 1])
        material_type = get_material_type(
            x_center, y_center, bone_width, fixator_thickness, fracture_params
        )

        avg_stresses[material_type] += vm_stresses[i]
        counts[material_type] += 1

    for key in avg_stresses:
        if counts[key] > 0:
            avg_stresses[key] /= counts[key]

    return avg_stresses


# ===================================================================
#  请用这个修正版的函数替换你的旧函数
# ===================================================================


def calculate_fracture_gap_strain(
    nodes, U, fixator_length, bone_length, fracture_width
):
    """计算骨折间隙的平均轴向应变 (已修正版本)"""
    fracture_start = (
        (fixator_length - bone_length) / 2 + bone_length / 2 - fracture_width / 2
    )
    fracture_end = fracture_start + fracture_width

    # 获取网格中所有不重复的x坐标
    unique_x_coords = np.unique(nodes[:, 0])

    # 找到离 fracture_start 最近且小于它的那一列节点的x坐标
    left_boundary_x = unique_x_coords[unique_x_coords < fracture_start].max()

    # 找到离 fracture_end 最近且大于它的那一列节点的x坐标
    right_boundary_x = unique_x_coords[unique_x_coords > fracture_end].min()

    # 根据找到的x坐标来选择节点
    left_gap_nodes = np.where(np.isclose(nodes[:, 0], left_boundary_x))[0]
    right_gap_nodes = np.where(np.isclose(nodes[:, 0], right_boundary_x))[0]

    # 过滤掉固定器部分的节点，只计算骨骼和骨痂区域的应变
    bone_width = 0.02
    fixator_thickness = 0.005

    left_gap_nodes_bone = [
        n
        for n in left_gap_nodes
        if fixator_thickness <= nodes[n, 1] <= bone_width + fixator_thickness
    ]
    right_gap_nodes_bone = [
        n
        for n in right_gap_nodes
        if fixator_thickness <= nodes[n, 1] <= bone_width + fixator_thickness
    ]

    if not left_gap_nodes_bone or not right_gap_nodes_bone:
        return 0.0

    # 提取这些节点的x方向位移
    ux_left = U[2 * np.array(left_gap_nodes_bone)]
    ux_right = U[2 * np.array(right_gap_nodes_bone)]

    # 计算实际的间隙宽度
    actual_gap_width = right_boundary_x - left_boundary_x

    # 计算平均位移差并除以实际间隙宽度得到应变
    avg_displacement_diff = np.mean(ux_right) - np.mean(ux_left)
    strain = avg_displacement_diff / actual_gap_width

    return strain


# --- 新增绘图函数 ---


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


# --- 主模拟函数（重构后）---


def run_simulation(simulation_params, output_dir):
    """运行单次模拟并返回分析结果"""
    # 解包参数
    (E_fixator, label) = simulation_params
    E_bone, E_callus_initial, E_callus_final = 18e9, 0.1e9, 18e9
    nu = 0.3
    bone_length, bone_width, fracture_width = 0.2, 0.02, 0.01
    fixator_thickness, fixator_length = 0.005, bone_length + 0.04
    applied_load = 1000
    num_steps = 20
    nx, ny = 40, 12

    fracture_params = (fixator_length, bone_length, fracture_width)

    # 创建网格
    nodes, elements = create_mesh(fixator_length, bone_width, fixator_thickness, nx, ny)

    # 初始化材料属性
    material_props = {
        "bone": {"E": E_bone, "nu": nu},
        "callus": {"E": E_callus_initial, "nu": nu},
        "fixator": {"E": E_fixator, "nu": nu},
    }

    # 存储每个时间步的分析结果
    results = {
        "avg_stress_fixator": [],
        "avg_stress_bone": [],
        "avg_stress_callus": [],
        "gap_strain": [],
    }

    # 时间步循环
    for step in range(num_steps):
        print(f"  - Running step {step + 1}/{num_steps}")
        # 更新骨痂弹性模量
        t = step / (num_steps - 1)
        material_props["callus"]["E"] = E_callus_initial + t * (
            E_callus_final - E_callus_initial
        )

        # 有限元求解
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

        # 计算并存储分析指标
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

    # 绘制本次模拟的分析图
    plot_analysis_results(results, output_dir)

    # --- 新增: 保存每组数据为CSV ---
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
    # --- 新增结束 ---

    return results


# --- 主程序入口 ---

if __name__ == "__main__":
    base_dir = setup_directories()

    # 定义参数化研究的参数
    # E_fixator for Aluminum is ~70 GPa
    parametric_studies = {
        "Flexible Fixator ": (35.0e9, "Flexible"),
        "Standard Fixator ": (70.0e9, "Standard"),
        "Rigid Fixator ": (140.0e9, "Rigid"),
    }

    all_results = {}

    for name, params in parametric_studies.items():
        print(f"--- Running Simulation for: {name} ---")

        # 为每次模拟创建单独的输出目录
        sim_output_dir = os.path.join(base_dir, params[1])
        os.makedirs(sim_output_dir)

        # 运行模拟并保存结果
        results = run_simulation(params, sim_output_dir)
        all_results[name] = results

    # 绘制参数化研究的最终对比图
    print("\n--- Generating Parametric Comparison Plot ---")
    plot_parametric_comparison(all_results, base_dir)

    print(f"\nSimulation completed. All results saved to '{base_dir}' directory.")

    # --- 新增: 汇总最终结果 ---
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
    # --- 新增结束 ---
