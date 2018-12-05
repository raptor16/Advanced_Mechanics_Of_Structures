"""
elementmatrix.m (element matrix calculation)
assembly.m (assemble global stiffness matrix)
solve.m (incorporate BCs and solve for displacements)
postprocess.m (calculate stresses and strains)
maindriver.m (read inputs, call assembly(), call solve() and call postprocess())
"""

import math
import numpy as np

def element_matrix(connec, nodal_coordinates, list_of_E, list_of_A, list_of_L): # pascal for E and meters for any lengths

    list_of_k = []
    list_of_l = []
    list_of_m = []
    i = 0
    for node_pair in connec:
        node_1_index = node_pair[0]-1
        node_2_index = node_pair[1]-1

        x2 = nodal_coordinates[node_2_index][0]  # 0 is the x-coord
        x1 = nodal_coordinates[node_1_index][0]
        y2 = nodal_coordinates[node_2_index][1]
        y1 = nodal_coordinates[node_1_index][1]

        length = math.sqrt((y2-y1)**2 + (x2-x1)**2)  # error if length is 0
        l = (x2 - x1)/length
        list_of_l.append(l)
        m = (y2 - y1)/length
        list_of_m.append(m)

        k = np.matrix([[l**2, l*m, -l**2, -l*m],
             [l*m, m**2, -l*m, -m**2],
             [-l**2, -l*m, l**2, l*m],
             [-l*m, -m**2, l*m, m**2]])
        k = list_of_E[i] * list_of_A[i]  / list_of_L[i] * k

        list_of_k.append(k)
        i += 1
    return list_of_k, list_of_l, list_of_m

def assembly(connec, nodal_coordinates, list_of_k):

    connec_dof_index = construct_connec_dof(connec)
    length_nodal_coord = len(nodal_coordinates)
    K = np.matrix([[0.] * length_nodal_coord*2] * length_nodal_coord*2) # change 8 to  # of nodes * 2

    k_counter = 0
    for k in list_of_k:
        for row_number in range (0, len(k)):
            for col_number in range (0, np.shape(k)[1]):
                K[connec_dof_index[k_counter, row_number], connec_dof_index[k_counter, col_number]] += \
                    k[row_number, col_number]
        k_counter += 1
    return K


def construct_connec_dof(connec):
    connec_dof_index = []
    for node_pair in connec:
        node_1 = node_pair[0]
        node_2 = node_pair[1]
        connec_dof_row_index = [node_1*2-1-1, node_1*2-1, node_2*2-1-1, node_2*2-1]  # subtract 1 again to adjust index
        connec_dof_index.append(connec_dof_row_index)
    connec_dof_index = np.matrix(connec_dof_index)
    return connec_dof_index

def solve(K, displacements, forces_applied):

    unknown_BC_indices = []
    known_BC_indices = []
    for i in range(0, len(displacements)):
        if (np.isnan(displacements[i])):
            unknown_BC_indices.append(i)
        else:
            known_BC_indices.append(i)

    K_free = np.matrix([[0.] * len(unknown_BC_indices)] * len(unknown_BC_indices))
    K_BC = np.matrix([[0.] * len(known_BC_indices)] * len(unknown_BC_indices))

    for row in range(len(unknown_BC_indices)):
        for col in range(len(unknown_BC_indices)):
            K_free[row, col] = K[unknown_BC_indices[row], unknown_BC_indices[col]]
    for row in range(len(unknown_BC_indices)):
        for col in range(len(known_BC_indices)):
            K_BC[row, col] = K[unknown_BC_indices[row], known_BC_indices[col]]

    BC = []
    for index in known_BC_indices:
        BC.append(displacements[index])

    F_free = [0.] * len(unknown_BC_indices)
    for i in range(len(unknown_BC_indices)):
        F_free[i] = forces_applied[unknown_BC_indices[i]]

    # K_free * disp_unk + K_BC * BC = F_free
    disp_unk = np.linalg.inv(K_free) * (np.transpose(np.matrix(F_free)) - np.matrix(K_BC) * np.transpose(np.matrix(BC)))

    j = 0
    for i in range(len(displacements)):
        if(np.isnan(displacements[i])):
            displacements[i] = np.asscalar(disp_unk[j])
            j += 1

    forces_reaction = K * np.transpose(np.matrix(displacements))
    return displacements, forces_reaction

def post_process(connec, displacements, list_of_E, list_of_L, list_of_l, list_of_m):

    stress = [0.] * len(connec)
    strain = [0.] * len(connec)
    d_element = np.matrix([[0.] * 4] * len(connec))
    connec_dof_index = construct_connec_dof(connec)

    for element in range(0, len(connec)):
        d_element[element, 0] = displacements[connec_dof_index[element, 0]]  # index adjustment necessitates the addition - 1
        d_element[element, 1] = displacements[connec_dof_index[element, 1]]  # index adjustment necessitates the addition - 1
        d_element[element, 2] = displacements[connec_dof_index[element, 2]]  # index adjustment necessitates the addition - 1
        d_element[element, 3] = displacements[connec_dof_index[element, 3]]  # index adjustment necessitates the addition - 1

    for element in range(len(connec)):
        tmp = np.matrix([-list_of_l[element], -list_of_m[element], list_of_l[element], list_of_m[element]])
        stress[element] = list_of_E[element] / list_of_L[element] * np.asscalar((tmp * np.transpose(np.matrix(d_element[element]))))
        strain[element] = stress[element] / list_of_E[element]
    return stress, strain


"""
connec: connectivity matrix, where the index represents the element numbers
nodal_coordinates: list of (x,y) of each nodes where the index represents the node number
"""
def main_driver(connec, nodal_coordinates, list_of_E, list_of_A, displacements, forces_applied):

    list_of_L = get_lengths(connec, nodal_coordinates)
    list_of_k, list_of_l, list_of_m = element_matrix(connec, nodal_coordinates, list_of_E, list_of_A, list_of_L)  # array of E's, Array of
    K = assembly(connec, nodal_coordinates, list_of_k)
    d, F_reaction = solve(K, displacements, forces_applied)
    stress, strain = post_process(connec, d, list_of_E, list_of_L, list_of_l, list_of_m)

    return stress


def get_lengths(connec, nodal_coordinates):
    list_of_L = []

    for i in range(len(connec)):
        node_1 = connec[i][0] - 1
        node_2 = connec[i][1] - 1
        dist_squared = (nodal_coordinates[node_2][0] - nodal_coordinates[node_1][0])**2 +\
               (nodal_coordinates[node_2][1] - nodal_coordinates[node_1][1])**2
        list_of_L.append(dist_squared**0.5)
    return list_of_L


"""
# INPUT FOR 3-bar truss
sample_connec = [[1, 4],
                 [2, 4],
                 [3, 4]]
nodal_coordinates = [[0, 0],
                     [1, 0],
                     [2, 0],
                     [1, -1]]

list_of_E = [70 * 10 ** 9] * 3
list_of_A = [1 * 10 ** -4] * 3
displacements = [0, 0, 0, 0, 0, 0, np.NAN, np.NAN]  # Use NAN for unknown boundary conditions
forces = [0, 0, 0, 0, 0, 0, 100, 0]

"""


# INPUT FOR 10-bar truss
sample_connec = [[1, 2],
                 [2, 3],
                 [3, 6],
                 [6, 5],
                 [5, 4],
                 [4, 2],
                 [1, 5],
                 [2, 5],
                 [2, 6],
                 [3, 5]]

nodal_coordinates = [[0, 1],
                     [1, 1],
                     [2, 1],
                     [0, 0],
                     [1, 0],
                     [2, 0]]

list_of_E = [70 * 10 ** 9] * 10
list_of_A = [1 * 10 ** -4] * 10
displacements = [0, 0, np.NAN, np.NAN, np.NAN, np.NAN, 0, 0, np.NAN, np.NAN, np.NAN, np.NAN]  # Use NAN for unknown boundary conditions
forces = [0, 0, 0, 0, 0, 0, 0, 0, 0, -100, 0, -100]


main_driver(sample_connec, nodal_coordinates, list_of_E, list_of_A, displacements, forces)
