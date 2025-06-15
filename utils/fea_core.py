import numpy as np

from utils.materials import get_material_type


def plane_stress_matrix(E, nu):
    return E / (1 - nu**2) * np.array([[1, nu, 0], [nu, 1, 0], [0, 0, (1 - nu) / 2]])


def element_stiffness(nodes, element, D):

    xi = np.array([-1 / np.sqrt(3), 1 / np.sqrt(3), 1 / np.sqrt(3), -1 / np.sqrt(3)])
    eta = np.array([-1 / np.sqrt(3), -1 / np.sqrt(3), 1 / np.sqrt(3), 1 / np.sqrt(3)])
    weights = np.array([1, 1, 1, 1])

    ke = np.zeros((8, 8))

    for gp in range(4):

        dN_dxi = 0.25 * np.array(
            [-(1 - eta[gp]), (1 - eta[gp]), (1 + eta[gp]), -(1 + eta[gp])]
        )
        dN_deta = 0.25 * np.array(
            [-(1 - xi[gp]), -(1 + xi[gp]), (1 + xi[gp]), (1 - xi[gp])]
        )

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

        dN_dx = invJ[0, 0] * dN_dxi + invJ[0, 1] * dN_deta
        dN_dy = invJ[1, 0] * dN_dxi + invJ[1, 1] * dN_deta

        B = np.zeros((3, 8))
        for i in range(4):
            B[0, 2 * i] = dN_dx[i]
            B[1, 2 * i + 1] = dN_dy[i]
            B[2, 2 * i] = dN_dy[i]
            B[2, 2 * i + 1] = dN_dx[i]

        ke += B.T @ D @ B * detJ * weights[gp]

    return ke


def assemble_global_stiffness(
    nodes, elements, material_props, bone_width, fixator_thickness, fracture_params
):
    num_nodes = nodes.shape[0]
    K = np.zeros((2 * num_nodes, 2 * num_nodes))

    for elem_idx, elem in enumerate(elements):

        x_center = np.mean(nodes[elem, 0])
        y_center = np.mean(nodes[elem, 1])
        material_type = get_material_type(
            x_center, y_center, bone_width, fixator_thickness, fracture_params
        )

        E = material_props[material_type]["E"]
        nu = material_props[material_type]["nu"]
        D = plane_stress_matrix(E, nu)

        ke = element_stiffness(nodes, elem, D)

        for i in range(4):
            for j in range(4):
                for di in range(2):
                    for dj in range(2):

                        row = 2 * elem[i] + di

                        col = 2 * elem[j] + dj

                        local_row = 2 * i + di
                        local_col = 2 * j + dj

                        K[row, col] += ke[local_row, local_col]
    return K


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
