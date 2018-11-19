import element_mats as em
import mesh
import numpy as np
import math
from scipy.linalg import eigh
import matplotlib.pyplot as plt


"""
 ElementStiMat.m (element stiness matrix calculation)
 ElementMassMat.m (element mass matrix calculation)
 assembly.m (assemble global stiness and mass matrices)
 freevibration.m (calculate natural frequencies and mode shapes)
 directfrfanalysis.m (calculate frequency response using the full-order dy-
namic stiness matrix)
 modalfrfanalysis.m (calculate frequency response using a modal analysis
approach)
 maindriver.m (read inputs, call assembly, incorporate BCs, call free vibra-
tion and frequency response analysis functions)
"""

def element_stiffness(EA, EI, x1, y1, x2, y2):
    """
    Calculates the element stiffness matrix

    :param EA:
    :param EI:
    :param x1:
    :param y1:
    :param x2:
    :param y2:

    :return:
    """
    return em.stiffness(EA, EI, x1, y1, x2, y2)


def element_mass_matrix(rho, x1, y1, x2, y2):
    """
    Calculates the element mass matrix
    :return:
    """
    return em.mass(rho, x1, y1, x2, y2)

def assembly(list_of_matrix, connec_mat, num_nodes):
    """
    Assembles the global stiffness matrix and mass matrices
    :return:
    """
    elements = list_of_matrix.shape[0]
    dof = num_nodes * 3
    shape = (elements, dof, dof)
    Mat = np.zeros(shape)  # This is a sparse matrix

    """
    Mass or Stiffness matrix will be a 6x6
      node1x  node1y node1z node2x node2y node2z
    [[ a      b      c      d      e      f]    node1x   
     [ g      h      i      j      k      l]    node1y 
     [ m      n      o      p      q      r]    node1z 
     [ s      t      u      v      w      x]    node2x 
     [ y      z      A      B      C      D]    node2y 
     [ E      F      G      H      I      J]]   node2z 
    
    I ran out of letters those capitals don't mean matrix.
    
    a: x, x
    b: y, x
    g: x, y
    h: y, y
    
    0->0 , +1, +2
    [0 1 2]
    
    [2 3 4] Node 1 at x y z hypothesis: *3-1
    
    [5 6 7] Node 2 at x y z 2*3-1 = 5
    """

    y_offset = 1
    z_offset = 2
    for i in range(elements):
        node1 = connec_mat[i][0]
        node2 = connec_mat[i][1]

        dof_starting_point1 = node1 * 3  # add by 1 for y and add by 2 for z of that node
        dof_starting_point2 = node2 * 3

        # TODO: use symmetric matrix
        Mat[i][dof_starting_point1][dof_starting_point1] = list_of_matrix[i][0][0]
        Mat[i][dof_starting_point1][dof_starting_point1 + y_offset] = list_of_matrix[i][0][1]
        Mat[i][dof_starting_point1][dof_starting_point1 + z_offset] = list_of_matrix[i][0][2]
        Mat[i][dof_starting_point1][dof_starting_point2] = list_of_matrix[i][0][3]
        Mat[i][dof_starting_point1][dof_starting_point2 + y_offset] = list_of_matrix[i][0][4]
        Mat[i][dof_starting_point1][dof_starting_point2 + z_offset] = list_of_matrix[i][0][5]

        Mat[i][dof_starting_point1 + y_offset][dof_starting_point1] = list_of_matrix[i][1][0]
        Mat[i][dof_starting_point1 + y_offset][dof_starting_point1 + y_offset] = list_of_matrix[i][1][1]
        Mat[i][dof_starting_point1 + y_offset][dof_starting_point1 + z_offset] = list_of_matrix[i][1][2]
        Mat[i][dof_starting_point1 + y_offset][dof_starting_point2] = list_of_matrix[i][1][3]
        Mat[i][dof_starting_point1 + y_offset][dof_starting_point2 + y_offset] = list_of_matrix[i][1][4]
        Mat[i][dof_starting_point1 + y_offset][dof_starting_point2 + z_offset] = list_of_matrix[i][1][5]

        Mat[i][dof_starting_point1 + z_offset][dof_starting_point1] = list_of_matrix[i][2][0]
        Mat[i][dof_starting_point1 + z_offset][dof_starting_point1 + y_offset] = list_of_matrix[i][2][1]
        Mat[i][dof_starting_point1 + z_offset][dof_starting_point1 + z_offset] = list_of_matrix[i][2][2]
        Mat[i][dof_starting_point1 + z_offset][dof_starting_point2] = list_of_matrix[i][2][3]
        Mat[i][dof_starting_point1 + z_offset][dof_starting_point2 + y_offset] = list_of_matrix[i][2][4]
        Mat[i][dof_starting_point1 + z_offset][dof_starting_point2 + z_offset] = list_of_matrix[i][2][5]

        Mat[i][dof_starting_point2][dof_starting_point1] = list_of_matrix[i][3][0]
        Mat[i][dof_starting_point2][dof_starting_point1 + y_offset] = list_of_matrix[i][3][1]
        Mat[i][dof_starting_point2][dof_starting_point1 + z_offset] = list_of_matrix[i][3][2]
        Mat[i][dof_starting_point2][dof_starting_point2] = list_of_matrix[i][3][3]
        Mat[i][dof_starting_point2][dof_starting_point2 + y_offset] = list_of_matrix[i][3][4]
        Mat[i][dof_starting_point2][dof_starting_point2 + z_offset] = list_of_matrix[i][3][5]

        Mat[i][dof_starting_point2 + y_offset][dof_starting_point1] = list_of_matrix[i][4][0]
        Mat[i][dof_starting_point2 + y_offset][dof_starting_point1 + y_offset] = list_of_matrix[i][4][1]
        Mat[i][dof_starting_point2 + y_offset][dof_starting_point1 + z_offset] = list_of_matrix[i][4][2]
        Mat[i][dof_starting_point2 + y_offset][dof_starting_point2] = list_of_matrix[i][4][3]
        Mat[i][dof_starting_point2 + y_offset][dof_starting_point2 + y_offset] = list_of_matrix[i][4][4]
        Mat[i][dof_starting_point2 + y_offset][dof_starting_point2 + z_offset] = list_of_matrix[i][4][5]

        Mat[i][dof_starting_point2 + z_offset][dof_starting_point1] = list_of_matrix[i][5][0]
        Mat[i][dof_starting_point2 + z_offset][dof_starting_point1 + y_offset] = list_of_matrix[i][5][1]
        Mat[i][dof_starting_point2 + z_offset][dof_starting_point1 + z_offset] = list_of_matrix[i][5][2]
        Mat[i][dof_starting_point2 + z_offset][dof_starting_point2] = list_of_matrix[i][5][3]
        Mat[i][dof_starting_point2 + z_offset][dof_starting_point2 + y_offset] = list_of_matrix[i][5][4]
        Mat[i][dof_starting_point2 + z_offset][dof_starting_point2 + z_offset] = list_of_matrix[i][5][5]
        #print "is Mat[i] symm", check_symmetric(Mat[i])

        """
        # BROKEN - Not?
        # Must do for node1x, node1y, node1z, etc.
        for j in range(6):
            for k in range(6):
                j_offset = j % 3  # 3 for x, y, and z
                k_offset = k % 3
                if (j >= 3):
                    dof_starting_pointx = dof_starting_point2
                else:
                    dof_starting_pointx = dof_starting_point1
                if (k >= 3):
                    dof_starting_pointy = dof_starting_point2
                else:
                    dof_starting_pointy = dof_starting_point1

                Mat[i][dof_starting_pointx + j_offset][dof_starting_pointy + k_offset] = list_of_matrix[i][j][k]
        """

    shape = (dof, dof)
    assembled_Mat = np.zeros(shape)
    for i in range(elements):
        assembled_Mat += Mat[i]
    return assembled_Mat


def check_symmetric(a, tol=1e-8):
    return np.allclose(a, a.T, atol=tol)


def eigensolver(M, K):
    # K \phi = \lambda M \phi ; \lamda is the eigenvalues and \phi is the eigenvectors
    eigval, eigvec = eigh(K, M, eigvals_only=False)
    omega = eigval ** 0.5
    # eigenvalues are the sqaure of natural frequencies and eigenvectors are the mode shape

    return eigval, eigvec, omega


def plot_mode_shape(eigvec):

    plt.ion()
    plt.plot(eigvec)

    plt.pause(1000)
    plt.show()


def free_vibration(M, K, local_mesh):
    """
    Question 2: Calculates natural frequencies and mode shapes
    :return:
    """
    # When solving the eigenproblem, take all the non restricted nodes aka active dof -- we want a good slicer for this
    # Therefore we eigensolve only 6:
    # We'll do this by brute force slicing
    reduced_M = M[local_mesh.dof_active, :][:, local_mesh.dof_active]
    reduced_K = K[local_mesh.dof_active, :][:, local_mesh.dof_active]
    eigval, eigvec, omega = eigensolver(M[local_mesh.dof_active, :][:, local_mesh.dof_active],
                                        K[local_mesh.dof_active, :][:, local_mesh.dof_active])

    hertz = omega / (2 * math.pi)

    return eigval, eigvec, omega, hertz, reduced_M, reduced_K

def direct_freq_analysis():
    """
    Calculates frequency response using the full-order dynamic stiffness matrix
    :return:
    """

    return NotImplemented

def modal_freq_analyis():
    """
    Calculates frequency response using modal analysis approach.
    :return:
    """

    return NotImplemented

def get_mesh_info (local_mesh):
    connec = local_mesh.NOD
    num_nodes = local_mesh.X.shape[0]
    X = local_mesh.X
    Y = local_mesh.Y
    num_elements = local_mesh.NOD.shape[0]  # number of elements

    return connec, num_nodes, X, Y, num_elements

def plot_excitation_freq(X, Y):
    plt.ion()
    plt.semilogy(X, Y)
    plt.pause(1000)
    plt.show()


def main():
    L = 1               # unit: m
    EI = 1.286 * 10**4   # Nm^2
    EA = 6.987 * 10**6   # N
    rhoA = 2.74           # kg/m
    x1, y1, x2, y2 = 0, 0, 0, 0


    # Get Mesh1
    local_mesh = mesh.Mesh1
    connec, num_nodes, X, Y, e = get_mesh_info(local_mesh)

    shape = (e, 6, 6)  # 6x6 mass matrix for e elements
    m_e = np.empty(shape)
    k_e = np.empty(shape)

    for i in range(e):
        node1 = connec[i][0]
        node2 = connec[i][1]
        m_e[i] = element_mass_matrix(rhoA, X[node1], Y[node1], X[node2], Y[node2])
        k_e[i] = element_stiffness(EA, EI, X[node1], Y[node1], X[node2], Y[node2])

    M = assembly(m_e, connec, num_nodes)
    K = assembly(k_e, connec, num_nodes)

    eigval, eigvec, omega, hertz, reduced_M, reduced_K = free_vibration(M, K, local_mesh)
    print "hertz", hertz

    # QUESTION 3
    """
       [K - \omega ^ 2 * M  + j * \omega * C] u_0 = f
               dynamic stiffness matrix
    """
    dof = 3
    fshape = (12*dof)
    f = np.zeros(fshape)
    y_offset = 1
    f[2*dof + y_offset] = 1  # Force at node 2 at the y dof is 1
    reduced_f = f[local_mesh.dof_active]
    print "reduced_f.shape", reduced_f.shape

    C = 10.0 * reduced_M
    j = 1j

    Omega_list = np.linspace(0, 250, 250)
    u0_at_R = np.zeros(Omega_list.shape[0])
    R_node = 6
    R_dof = R_node * dof - 6 + y_offset # minus 6 because there's two nodes that are dead.
    print "R_dof", R_dof

    # Omega = 250 * 2 * math.pi
    for i in range (Omega_list.shape[0]):
        Omega = Omega_list[i] * 2 * math.pi # Convert to Radians
        dynamic_stiffness_matrix = reduced_K - Omega ** 2 * reduced_M + j * Omega * C   # 18 x 18
        u0 = np.linalg.solve(dynamic_stiffness_matrix, reduced_f)   # 18 x 1
        u0_at_R[i] = abs(u0[R_dof])
    print "u0_at_R", u0_at_R
    plot_excitation_freq(Omega_list, u0_at_R)
    # Find the u0 at R_dof and plot it against the frequencies from 0 to 250.
    # Then plot u0_at_R against the frequencies.



    # Get Mesh2
    connec, num_nodes, X, Y, e = get_mesh_info(mesh.Mesh2)
    shape = (e, 6, 6)  # 6x6 mass matrix for e elements
    m_e = np.empty(shape)
    k_e = np.empty(shape)
    for i in range(e):
        node1 = connec[i][0]
        node2 = connec[i][1]
        m_e[i] = element_mass_matrix(rhoA, X[node1], Y[node1], X[node2], Y[node2])
        k_e[i] = element_stiffness(EA, EI, X[node1], Y[node1], X[node2], Y[node2])
    M = assembly(m_e, connec, num_nodes)
    K = assembly(k_e, connec, num_nodes)
    eigval, eigvec, omega, hertz, reduced_M, reduced_K = free_vibration(M, K, mesh.Mesh2)



    # Get Mesh3
    connec, num_nodes, X, Y, e = get_mesh_info(mesh.Mesh3)
    shape = (e, 6, 6)  # 6x6 mass matrix for e elements
    m_e = np.empty(shape)
    k_e = np.empty(shape)
    for i in range(e):
        node1 = connec[i][0]
        node2 = connec[i][1]
        m_e[i] = element_mass_matrix(rhoA, X[node1], Y[node1], X[node2], Y[node2])
        k_e[i] = element_stiffness(EA, EI, X[node1], Y[node1], X[node2], Y[node2])
    M = assembly(m_e, connec, num_nodes)
    K = assembly(k_e, connec, num_nodes)
    eigval, eigvec, omega, hertz, reduced_M, reduced_K = free_vibration(M, K, mesh.Mesh3)



    # Get Mesh4
    connec, num_nodes, X, Y, e = get_mesh_info(mesh.Mesh4)
    shape = (e, 6, 6)  # 6x6 mass matrix for e elements
    m_e = np.empty(shape)
    k_e = np.empty(shape)
    for i in range(e):
        node1 = connec[i][0]
        node2 = connec[i][1]
        m_e[i] = element_mass_matrix(rhoA, X[node1], Y[node1], X[node2], Y[node2])
        k_e[i] = element_stiffness(EA, EI, X[node1], Y[node1], X[node2], Y[node2])
    M = assembly(m_e, connec, num_nodes)
    K = assembly(k_e, connec, num_nodes)
    eigval, eigvec, omega, hertz, reduced_M, reduced_K = free_vibration(M, K, mesh.Mesh4)


    #plot_mode_shape(eigvec)

    # ^^ RINSE WASH REPEAT FOR THE REST OF THE MESHES ^^



    # incorporates BCs
    direct_freq_analysis()
    modal_freq_analyis()


"""
Part A
Study the convergence of the first 5 natural frequencies and mode shapes of
the structure when the finite element spatial mesh is refined.
"""

main()
