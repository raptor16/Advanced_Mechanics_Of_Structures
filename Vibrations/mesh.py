"""
Pre-made Meshes for AER501 Assignment 2 Part B, Ported to Python 2.7 from MATLAB

Mesh1.m -> mesh.Mesh1
Mesh2.m -> mesh.Mesh2
Mesh3.m -> mesh.Mesh3
Mesh4.m -> mesh.Mesh4
"""
# Let K and M denote the global stiffness and mass matrices after assembly.
# The rows and columns corresponding to the zero displacement BCs can be
# deleted as follows:
#
# K = K(dof_active, dof_active)
# M = M(dof_active, dof_active)
#
# Using the mesh information and Fig 1 in Assignment2.pdf you can identify
# the DOF where the external force is applied and the DOF whose displacement
# response is of interest.

import numpy as np

class _Mesh(object):
    """Container class for mesh data"""
    NOD = None # connectivity matrix
    X = None # vector containing x-coordinates of each node
    Y = None # vector containing y-coordinates of each node
    constr_node = None # Indices of constrained nodes

    # indices of active and inactive DOFs
    dof_restr = None
    dof_active = None

def _active_inactive_dofs(constr_node_2, length_x):
    dof_restr = np.array([0, 1, 2, 3 * constr_node_2, 3 * constr_node_2 + 1, 3 * constr_node_2 + 2]) # indices of restrained DOF
    dof_active = np.concatenate((np.arange(3, 3 * constr_node_2), np.arange(3 * constr_node_2 + 3, 3 * length_x))) # indices of unrestrained DOF
    return dof_restr, dof_active


# Mesh1
#
#  1     3     5     7
# |>o-----o-----o-----o
#    \    |\    |\    |
#      \  |  \  |  \  |
#        \|    \|    \|
# |>o-----o-----o-----o
#   0     2     4     6
#
#
# Mesh with one element per member
#
Mesh1 = _Mesh()
Mesh1.NOD = np.array([
    [0, 2],
    [2, 3],
    [1, 3],
    [1, 2],
    [2, 4],
    [4, 5],
    [3, 5],
    [3, 4],
    [4, 6],
    [6, 7],
    [5, 7],
    [5, 6],
])
Mesh1.X = np.array([0,0,1,1,2,2,3,3])
Mesh1.Y = np.array([0,1,0,1,0,1,0,1])
Mesh1.constr_node = np.array([0, 1])
Mesh1.dof_restr = np.arange(0, 6)
Mesh1.dof_active = np.arange(6, 24)


# Mesh2
#
# |>o-----o-----o-----o
#    \    |\    |\    |
#      \  |  \  |  \  |
#        \|    \|    \|
# |>o-----o-----o-----o
#
# Mesh with two elements per member
#
Mesh2 = _Mesh()
Mesh2.NOD = np.array([
    [0, 1],
    [1, 2],
    [2, 3],
    [3, 4],
    [4, 5],
    [5, 6],
    [6, 7],
    [7, 2],
    [2, 8],
    [8, 9],
    [9, 10],
    [10, 11],
    [11, 12],
    [12, 4],
    [4, 13],
    [13, 9],
    [9, 14],
    [14, 15],
    [15, 16],
    [16, 17],
    [17, 18],
    [18, 11],
    [11, 19],
    [19, 15],
])
Mesh2.X = np.array([ 0.00000, 0.50000, 1.00000, 1.00000, 1.00000, 0.50000, 0.00000, 0.50000, 1.50000, 2.00000, 2.00000, 2.00000, 1.50000, 1.50000, 2.50000, 3.00000, 3.00000, 3.00000, 2.50000, 2.50000,])
Mesh2.Y = np.array([ 0.00000, 0.00000, 0.00000, 0.50000, 1.00000, 1.00000, 1.00000, 0.50000, 0.00000, 0.00000, 0.50000, 1.00000, 1.00000, 0.50000, 0.00000, 0.00000, 0.50000, 1.00000, 1.00000, 0.50000,])
Mesh2.constr_node = np.array([0, 6])
Mesh2.dof_restr, Mesh2.dof_active = _active_inactive_dofs(Mesh2.constr_node[1], len(Mesh2.X))


# Mesh3
#
# |>o-----o-----o-----o
#    \    |\    |\    |
#      \  |  \  |  \  |
#        \|    \|    \|
# |>o-----o-----o-----o
#
#
# Mesh with three elements per member
#
Mesh3 = _Mesh()
Mesh3.NOD = np.array([
    [0, 1],
    [1, 2],
    [2, 3],
    [3, 4],
    [4, 5],
    [5, 6],
    [6, 7],
    [7, 8],
    [8, 9],
    [9, 10],
    [10, 11],
    [11, 3],
    [3, 12],
    [12, 13],
    [13, 14],
    [14, 15],
    [15, 16],
    [16, 17],
    [17, 18],
    [18, 19],
    [19, 6],
    [6, 20],
    [20, 21],
    [21, 14],
    [14, 22],
    [22, 23],
    [23, 24],
    [24, 25],
    [25, 26],
    [26, 27],
    [27, 28],
    [28, 29],
    [29, 17],
    [17, 30],
    [30, 31],
    [31, 24],
])
Mesh3.X = np.array([ 0.00000, 0.33333, 0.66667, 1.00000, 1.00000, 1.00000, 1.00000, 0.66667, 0.33333, 0.00000, 0.33333, 0.66667, 1.33333, 1.66667, 2.00000, 2.00000, 2.00000, 2.00000, 1.66667, 1.33333, 1.33333, 1.66667, 2.33333, 2.66667, 3.00000, 3.00000, 3.00000, 3.00000, 2.66667, 2.33333, 2.33333, 2.66667,])
Mesh3.Y = np.array([ 0.00000, 0.00000, 0.00000, 0.00000, 0.33333, 0.66667, 1.00000, 1.00000, 1.00000, 1.00000, 0.66667, 0.33333, 0.00000, 0.00000, 0.00000, 0.33333, 0.66667, 1.00000, 1.00000, 1.00000, 0.66667, 0.33333, 0.00000, 0.00000, 0.00000, 0.33333, 0.66667, 1.00000, 1.00000, 1.00000, 0.66667, 0.33333,])
Mesh3.constr_node = np.array([0, 9])
Mesh3.dof_restr, Mesh3.dof_active = _active_inactive_dofs(Mesh3.constr_node[1], len(Mesh3.X))


# Mesh4
#
# |>o-----o-----o-----o
#    \    |\    |\    |
#      \  |  \  |  \  |
#        \|    \|    \|
# |>o-----o-----o-----o
#
# Mesh with four elements per member
#
Mesh4 = _Mesh()
Mesh4.NOD = np.array([
    [0, 1],
    [1, 2],
    [2, 3],
    [3, 4],
    [4, 5],
    [5, 6],
    [6, 7],
    [7, 8],
    [8, 9],
    [9, 10],
    [10, 11],
    [11, 12],
    [12, 13],
    [13, 14],
    [14, 15],
    [15, 4],
    [4, 16],
    [16, 17],
    [17, 18],
    [18, 19],
    [19, 20],
    [20, 21],
    [21, 22],
    [22, 23],
    [23, 24],
    [24, 25],
    [25, 26],
    [26, 8],
    [8, 27],
    [27, 28],
    [28, 29],
    [29, 19],
    [19, 30],
    [30, 31],
    [31, 32],
    [32, 33],
    [33, 34],
    [34, 35],
    [35, 36],
    [36, 37],
    [37, 38],
    [38, 39],
    [39, 40],
    [40, 23],
    [23, 41],
    [41, 42],
    [42, 43],
    [43, 33],
])
Mesh4.X = np.array([ 0.00000, 0.25000, 0.50000, 0.75000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 0.75000, 0.50000, 0.25000, 0.00000, 0.25000, 0.50000, 0.75000, 1.25000, 1.50000, 1.75000, 2.00000, 2.00000, 2.00000, 2.00000, 2.00000, 1.75000, 1.50000, 1.25000, 1.25000, 1.50000, 1.75000, 2.25000, 2.50000, 2.75000, 3.00000, 3.00000, 3.00000, 3.00000, 3.00000, 2.75000, 2.50000, 2.25000, 2.25000, 2.50000, 2.75000,])
Mesh4.Y = np.array([ 0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.25000, 0.50000, 0.75000, 1.00000, 1.00000, 1.00000, 1.00000, 1.00000, 0.75000, 0.50000, 0.25000, 0.00000, 0.00000, 0.00000, 0.00000, 0.25000, 0.50000, 0.75000, 1.00000, 1.00000, 1.00000, 1.00000, 0.75000, 0.50000, 0.25000, 0.00000, 0.00000, 0.00000, 0.00000, 0.25000, 0.50000, 0.75000, 1.00000, 1.00000, 1.00000, 1.00000, 0.75000, 0.50000, 0.25000,])
Mesh4.constr_node = np.array([0, 12])
Mesh4.dof_restr, Mesh4.dof_active = _active_inactive_dofs(Mesh4.constr_node[1], len(Mesh4.X))
