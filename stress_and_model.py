import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import LinearSegmentedColormap
from scipy.interpolate import griddata

# ...existing code...
import os

# 确保输出目录存在
os.makedirs("output", exist_ok=True)

# 材料参数
E_bone = 18.0e9  # 骨弹性模量 (Pa)
E_callus_initial = 0.1e9  # 初始骨痂弹性模量 (Pa)
E_callus_final = 18.0e9  # 最终骨痂弹性模量 (Pa)
E_fixator = 70.0e9  # 外固定器弹性模量 (Pa)
nu = 0.3  # 泊松比

# 几何参数
bone_length = 0.2  # 骨长度 (m)
bone_width = 0.02  # 骨宽度 (m)
fracture_width = 0.01  # 骨折区域宽度 (m)
fixator_thickness = 0.005  # 外固定器厚度 (m)
fixator_length = bone_length + 0.04  # 外固定器长度 (m)

# 载荷参数
applied_load = 1000  # 施加的载荷 (N)

# 模拟参数
num_steps = 20  # 时间步数
dt = 1.0  # 时间步长 (天)

# 有限元网格参数
nx = 40  # x方向单元数
ny = 12  # y方向单元数


# 创建网格
def create_mesh():
    # 计算节点坐标
    x = np.linspace(0, fixator_length, nx + 1)
    y = np.linspace(0, bone_width + 2 * fixator_thickness, ny + 1)

    # 创建节点坐标矩阵
    nodes = np.zeros(((nx + 1) * (ny + 1), 2))
    for i in range(nx + 1):
        for j in range(ny + 1):
            nodes[i * (ny + 1) + j, 0] = x[i]
            nodes[i * (ny + 1) + j, 1] = y[j]

    # 创建单元连接矩阵
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
def get_material_type(x, y):
    # 外固定器区域
    if y < fixator_thickness or y > bone_width + fixator_thickness:
        return "fixator"

    # 骨折区域
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


# 组装全局刚度矩阵
def assemble_global_stiffness(nodes, elements, material_props):
    num_nodes = nodes.shape[0]
    K = np.zeros((2 * num_nodes, 2 * num_nodes))

    for elem_idx, elem in enumerate(elements):
        # 获取单元中心坐标以确定材料类型
        x_center = np.mean(nodes[elem, 0])
        y_center = np.mean(nodes[elem, 1])
        material_type = get_material_type(x_center, y_center)

        # 获取材料属性
        E = material_props[material_type]["E"]
        nu = material_props[material_type]["nu"]
        D = plane_stress_matrix(E, nu)

        # 计算单元刚度矩阵
        ke = element_stiffness(nodes, elem, D)

        # 组装到全局矩阵
        for i in range(4):
            for j in range(4):
                for di in range(2):
                    for dj in range(2):
                        row = 2 * elem[i] + di
                        col = 2 * elem[j] + dj
                        K[row, col] += ke[2 * i + di, 2 * j + dj]

    return K


# 施加边界条件和载荷
def apply_bc_and_loads(nodes, K, F):
    # 固定左侧边界 (x=0)
    left_nodes = np.where(nodes[:, 0] == 0)[0]
    for node in left_nodes:
        F[2 * node] = 0  # ux = 0
        F[2 * node + 1] = 0  # uy = 0
        K[2 * node, :] = 0
        K[:, 2 * node] = 0
        K[2 * node, 2 * node] = 1
        K[2 * node + 1, :] = 0
        K[:, 2 * node + 1] = 0
        K[2 * node + 1, 2 * node + 1] = 1

    # 在右侧施加拉力 (x=fixator_length)
    right_nodes = np.where(nodes[:, 0] == fixator_length)[0]
    total_nodes = len(right_nodes)
    load_per_node = applied_load / total_nodes

    for node in right_nodes:
        F[2 * node] = load_per_node  # x方向拉力

    return K, F


# 计算应力
def calculate_stress(nodes, elements, U, material_props):
    num_elements = elements.shape[0]
    element_stresses = np.zeros((num_elements, 3))  # sigma_xx, sigma_yy, tau_xy

    for elem_idx, elem in enumerate(elements):
        # 获取单元中心坐标以确定材料类型
        x_center = np.mean(nodes[elem, 0])
        y_center = np.mean(nodes[elem, 1])
        material_type = get_material_type(x_center, y_center)

        # 获取材料属性
        E = material_props[material_type]["E"]
        nu = material_props[material_type]["nu"]
        D = plane_stress_matrix(E, nu)

        # 单元位移
        Ue = np.zeros(8)
        for i in range(4):
            Ue[2 * i] = U[2 * elem[i]]
            Ue[2 * i + 1] = U[2 * elem[i] + 1]

        # 在单元中心计算应力
        xi, eta = 0, 0  # 单元中心

        # 形函数导数
        dN_dxi = 0.25 * np.array([-(1 - eta), (1 - eta), (1 + eta), -(1 + eta)])
        dN_deta = 0.25 * np.array([-(1 - xi), -(1 + xi), (1 + xi), (1 - xi)])

        # 雅可比矩阵
        J = np.zeros((2, 2))
        for i in range(4):
            node_idx = elem[i]
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

        # 计算应力
        stress = D @ B @ Ue
        element_stresses[elem_idx, :] = stress

    return element_stresses


# 计算von Mises应力
def calculate_von_mises(stress):
    sxx = stress[:, 0]
    syy = stress[:, 1]
    txy = stress[:, 2]
    return np.sqrt(sxx**2 + syy**2 - sxx * syy + 3 * txy**2)


# 可视化函数
def plot_model(nodes, elements, material_props, step, U=None, stresses=None):
    fig, ax = plt.subplots(figsize=(12, 6))

    # 计算当前骨痂与最终骨痂的模量比，用于颜色渐变
    E_callus = material_props["callus"]["E"]
    healing_ratio = (E_callus - E_callus_initial) / (E_callus_final - E_callus_initial)

    for elem in elements:
        x = nodes[elem, 0]
        y = nodes[elem, 1]

        # 获取材料类型
        x_center = np.mean(x)
        y_center = np.mean(y)
        material_type = get_material_type(x_center, y_center)

        # 设置颜色
        if material_type == "bone":
            color = "lightgray"
        elif material_type == "callus":
            # 根据愈合程度改变颜色 (从浅绿到深绿)
            green_value = 0.3 + 0.7 * healing_ratio
            color = (0.0, green_value, 0.0)  # RGB格式
        else:  # fixator
            # 外固定器颜色可以随时间变浅，模拟刚度降低
            blue_value = 0.7 - 0.5 * healing_ratio
            color = (0.0, 0.0, max(0.2, blue_value))  # 保持最小蓝色值

        # 绘制单元
        ax.fill(
            x[[0, 1, 2, 3, 0]],
            y[[0, 1, 2, 3, 0]],
            color=color,
            edgecolor="k",
            linewidth=0.5,
        )

    # 如果有位移，绘制变形后的形状
    if U is not None:
        scale = 100  # 位移放大系数
        deformed_nodes = nodes.copy()
        deformed_nodes[:, 0] += scale * U[::2]
        deformed_nodes[:, 1] += scale * U[1::2]

        for elem in elements:
            x = deformed_nodes[elem, 0]
            y = deformed_nodes[elem, 1]
            ax.plot(x[[0, 1, 2, 3, 0]], y[[0, 1, 2, 3, 0]], "r-", linewidth=0.5)

    # 设置标题和标签
    ax.set_title(f"Fracture Healing Simulation - Step {step}")
    ax.set_xlabel("Length (m)")
    ax.set_ylabel("Height (m)")
    ax.grid(True)

    # 添加图例
    bone_patch = plt.Rectangle(
        (0, 0), 1, 1, fc="lightgray", edgecolor="k", label="Bone"
    )
    callus_patch = plt.Rectangle(
        (0, 0), 1, 1, fc="lightgreen", edgecolor="k", label="Callus"
    )
    fixator_patch = plt.Rectangle(
        (0, 0), 1, 1, fc="lightblue", edgecolor="k", label="Fixator"
    )
    ax.legend(handles=[bone_patch, callus_patch, fixator_patch])

    # 保存图像
    plt.savefig(f"output/model_images/model_step_{step:02d}.png")
    plt.close()

    # 如果有应力数据，绘制应力云图
    if stresses is not None:
        fig, ax = plt.subplots(figsize=(12, 6))

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

        # 创建应力云图
        x = nodes[:, 0]
        y = nodes[:, 1]
        z = node_stresses

        # 创建网格用于绘图
        xi = np.linspace(x.min(), x.max(), 100)
        yi = np.linspace(y.min(), y.max(), 100)
        zi = griddata((x, y), z, (xi[None, :], yi[:, None]), method="cubic")

        # 绘制云图
        levels = np.linspace(z.min(), z.max(), 20)
        cs = ax.contourf(xi, yi, zi, levels=levels, cmap="jet", extend="both")

        # 添加颜色条
        cbar = fig.colorbar(cs)
        cbar.set_label("Von Mises Stress (Pa)")

        # 绘制网格
        for elem in elements:
            x = nodes[elem, 0]
            y = nodes[elem, 1]
            ax.plot(x[[0, 1, 2, 3, 0]], y[[0, 1, 2, 3, 0]], "k-", linewidth=0.2)

        # 设置标题和标签
        ax.set_title(f"Stress Distribution - Step {step}")
        ax.set_xlabel("Length (m)")
        ax.set_ylabel("Height (m)")

        # 保存应力云图
        plt.savefig(f"output/stress_images/stress_step_{step:02d}.png")
        plt.close()


# 主模拟函数
def simulate_healing():
    # 创建网格
    nodes, elements = create_mesh()

    # 初始化材料属性
    material_props = {
        "bone": {"E": E_bone, "nu": nu},
        "callus": {"E": E_callus_initial, "nu": nu},
        "fixator": {"E": E_fixator, "nu": nu},
    }

    # 初始化结果存储
    model_images = []
    stress_images = []

    # 时间步循环
    for step in range(num_steps):
        print(f"Processing step {step + 1}/{num_steps}")

        # 更新骨痂弹性模量 (模拟愈合过程)
        t = step / (num_steps - 1)  # 归一化时间 [0,1]
        E_callus = E_callus_initial + t * (E_callus_final - E_callus_initial)
        material_props["callus"]["E"] = E_callus

        # 组装刚度矩阵
        K = assemble_global_stiffness(nodes, elements, material_props)

        # 初始化载荷向量
        F = np.zeros(2 * nodes.shape[0])

        # 施加边界条件和载荷
        K, F = apply_bc_and_loads(nodes, K, F)

        # 求解位移
        U = np.linalg.solve(K, F)

        # 计算应力
        stresses = calculate_stress(nodes, elements, U, material_props)

        # 可视化模型和应力
        plot_model(nodes, elements, material_props, step, U, stresses)

        # 保存图像路径
        model_images.append(f"output/model_images/model_step_{step:02d}.png")
        stress_images.append(f"output/stress_images/stress_step_{step:02d}.png")

    # 创建动画
    def create_animation(image_files, output_name, fps=2):
        fig = plt.figure(figsize=(12, 6))
        ims = []

        for file in image_files:
            img = plt.imread(file)
            im = plt.imshow(img, animated=True)
            ims.append([im])

        ani = animation.ArtistAnimation(fig, ims, interval=500, blit=True)
        ani.save(f"./assets/{output_name}.gif", writer="pillow", fps=fps)
        plt.close()

    # 生成模型动画
    create_animation(model_images, "healing_process")

    # 生成应力动画
    create_animation(stress_images, "stress_distribution")


# 运行模拟
if __name__ == "__main__":
    simulate_healing()
    print("Simulation completed. Results saved to 'output' directory.")
