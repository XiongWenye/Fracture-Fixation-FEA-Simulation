import numpy as np

from utils.materials import get_material_type
from utils.fem_core import plane_stress_matrix


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


# 计算von Mises应力
def calculate_von_mises(stress):
    sxx, syy, txy = stress[:, 0], stress[:, 1], stress[:, 2]
    return np.sqrt(sxx**2 + syy**2 - sxx * syy + 3 * txy**2)


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
