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


def plot_mode_shape(eigvec, eigvec1, eigvec2, eigvec3):

    plt.ion()
    plt.subplot(2, 2, 1)
    plt.plot(eigvec, label="Eigenvector(Mode Shape 1)")
    plt.title("Eigenvector(Mode Shape) 1")
    plt.subplot(2, 2, 2)
    plt.plot(eigvec1, label="Eigenvector(Mode Shape 1)")
    plt.title("Eigenvector(Mode Shape) 2")
    plt.subplot(2, 2, 3)
    plt.plot(eigvec2, label="Eigenvector(Mode Shape 1)")
    plt.title("Eigenvector(Mode Shape) 3")
    plt.subplot(2, 2, 4)
    plt.plot(eigvec3, label="Eigenvector(Mode Shape 1)")
    plt.title("Eigenvector(Mode Shape) 4")
    # plt.legend()
    plt.pause(1000)
    plt.show()


def free_vibration(M, K, local_mesh):
    """
    Question 2: Calculates natural frequencies and mode shapes
    :return:
    """
    reduced_M = M[local_mesh.dof_active, :][:, local_mesh.dof_active]
    reduced_K = K[local_mesh.dof_active, :][:, local_mesh.dof_active]
    eigval, eigvec, omega = eigensolver(M[local_mesh.dof_active, :][:, local_mesh.dof_active],
                                        K[local_mesh.dof_active, :][:, local_mesh.dof_active])

    hertz = omega / (2 * math.pi)

    return eigval, eigvec, omega, hertz, reduced_M, reduced_K

def direct_freq_analysis(rhoA, EA, EI, local_mesh, dof, y_offset, reduced_M, reduced_K, Omega_list):
    """
    Calculates frequency response using the full-order dynamic stiffness matrix
    :return:
    """

    # QUESTION 3 A
    #
    #    [K - \omega ^ 2 * M  + j * \omega * C] u_0 = f
    #            dynamic stiffness matrix
    #
    # TODO: cleanup the static constants (eqivalent for Python) for variables dof and y_offset
    j = 1j
    connec, num_nodes, X, Y, e = get_mesh_info(local_mesh)

    members = connec.shape[0]
    fshape = (members * dof)
    f = np.zeros(fshape)
    P_dof = np.where((X == 1) & (Y == 0))[0][0] * dof + y_offset
    f[P_dof] = 1  # Force at node 2 at the y dof is 1
    reduced_f = f[local_mesh.dof_active]
    C = 10.0 * reduced_M

    u0_at_R = np.zeros(Omega_list.shape[0])
    # R_node = 6
    R_dof = np.where((X == 3) & (Y == 0))[0][0] * dof - len(local_mesh.dof_restr) + y_offset
    # R_dof = R_node * dof - len(local_mesh.dof_restr) + y_offset  # minus 6 because there's two nodes that are dead.

    # Omega = 250 * 2 * math.pi
    for i in range(Omega_list.shape[0]):
        Omega = Omega_list[i] * 2 * math.pi  # Convert to Radians
        dynamic_stiffness_matrix = reduced_K - Omega ** 2 * reduced_M + j * Omega * C  # 18 x 18
        u0 = np.linalg.solve(dynamic_stiffness_matrix, reduced_f)  # 18 x 1
        u0_at_R[i] = abs(u0[R_dof])

    return u0_at_R, reduced_f, R_dof

def modal_freq_analyis(eigvec, eigval, reduced_f, R_dof, Omega_list):
    """
    Calculates frequency response using modal analysis approach.
    :return:
    """
    # Question 3 B
    gamma2 = 10
    j = 1j
    num_of_eig = len(eigval)
    num_of_omega = len(Omega_list)
    q_shape = (num_of_omega, num_of_eig)
    u0_modal = np.zeros(q_shape, dtype=complex)
    print "eigvec", eigvec.shape

    for k in range(len(Omega_list)):
        for i in range(num_of_eig):
            Omega = Omega_list[k] * 2 * math.pi  # Convert to Radians
            q_at_Omega = np.dot(eigvec[:, i], reduced_f) / (eigval[i] - Omega**2 + j* Omega * gamma2)
            u0_modal[k] = np.add(u0_modal[k], q_at_Omega * eigvec[:, i])
    # plot_modal(Omega_list, u0)
    u0_modal = np.absolute(u0_modal)
    plot_modal(Omega_list, u0_modal[:, R_dof])


def get_mesh_info (local_mesh):
    connec = local_mesh.NOD
    num_nodes = local_mesh.X.shape[0]
    X = local_mesh.X
    Y = local_mesh.Y
    num_elements = local_mesh.NOD.shape[0]  # number of elements

    return connec, num_nodes, X, Y, num_elements

def plot_natural_freq(nat_freq, f1, f2, f3):
    # Part 2
    #n = 4
    #x = np.linspace(0, n, n)

    plt.ion()
    plt.plot(nat_freq, label="Mesh1 - 1 Element per Beam ")
    plt.plot(f1, label="Mesh2 - 2 Elements per Beam")
    plt.plot(f2, label="Mesh3 - 3 Elements per Beam")
    plt.plot(f3, label="Mesh4 - 4 Elements per Beam")
    plt.legend()
    plt.title("Natural Frequency vs. Frequency Mode")
    plt.xlabel("Eigenmodes")
    plt.ylabel("Natural frequency in Hz")
    plt.pause(2000)
    plt.show()


def plot_excitation_freq(X, Y, Y1, Y2, Y3):
    """

    :param X: The linspace of frequencies that is to be calculated for each Mesh's u0 at R
    :param Y: u0 at R for Mesh 1
    :param Y1: u0 at R for Mesh 2
    :param Y2: u0 at R for Mesh 3
    :param Y3: u0 at R for Mesh 4
    :return: None -- a figure will popup with the plots.
    """
    plt.ion()
    plt.semilogy(X, Y, label="Mesh 1")
    plt.semilogy(X, Y1, label="Mesh 2")
    plt.semilogy(X, Y2, label="Mesh 3")
    plt.semilogy(X, Y3, label="Mesh 4")
    plt.legend()
    plt.title("Excitation Frequency ")
    plt.xlabel("Frequency in Hz")
    plt.ylabel("u0 in m)")
    plt.pause(1000)
    plt.show()

def plot_modal (X, Y):
    plt.ion()
    plt.semilogy(X, Y, label="Mesh 4")
    plt.legend()
    plt.title("Modal Analysis")
    plt.xlabel("Frequency in Hz")
    plt.ylabel("u0 in m)")
    plt.pause(1000)
    plt.show()

def main():
    EI = 1.286 * 10**4   # Nm^2
    EA = 6.987 * 10**6   # N
    rhoA = 2.74           # kg/m
    dof = 3
    y_offset = 1
    Omega_list = np.linspace(0, 250, 250)  # Omega_list (big Omega) is the X axis Freq from 0 to 250 in Hz

    ### MESH 1 ###
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
    print "Mesh 1 Frequencies in Hz, first 12", hertz[0:12]

    u0_at_R, _, _ = direct_freq_analysis(rhoA, EA, EI, local_mesh, dof, y_offset, reduced_M, reduced_K, Omega_list)
    ### END MESH 1 ###

    ### MESH 2 ###
    local_mesh = mesh.Mesh2
    connec1, num_nodes1, X1, Y1, e1 = get_mesh_info(local_mesh)
    shape1 = (e1, 6, 6)  # 6x6 mass matrix for e elements
    m_e1 = np.empty(shape1)
    k_e1 = np.empty(shape1)
    for i in range(e1):
        node1 = connec1[i][0]
        node2 = connec1[i][1]
        m_e1[i] = element_mass_matrix(rhoA, X1[node1], Y1[node1], X1[node2], Y1[node2])
        k_e1[i] = element_stiffness(EA, EI, X1[node1], Y1[node1], X1[node2], Y1[node2])
    M1 = assembly(m_e1, connec1, num_nodes1)
    K1 = assembly(k_e1, connec1, num_nodes1)
    eigval1, eigvec1, omega1, hertz1, reduced_M, reduced_K = free_vibration(M1, K1, local_mesh)
    print "Mesh 2 Frequencies in Hz, first 12", hertz1[0:12]

    u0_at_R1, _, _ = direct_freq_analysis(rhoA, EA, EI, local_mesh, dof, y_offset, reduced_M, reduced_K, Omega_list)


    ### END MESH 2


    ### MESH 3 ###
    # Get Mesh3
    local_mesh = mesh.Mesh3
    connec2, num_nodes2, X2, Y2, e2 = get_mesh_info(local_mesh)
    shape2 = (e2, 6, 6)  # 6x6 mass matrix for e elements
    m_e2 = np.empty(shape2)
    k_e2 = np.empty(shape2)
    for i in range(e2):
        node1 = connec2[i][0]
        node2 = connec2[i][1]
        m_e2[i] = element_mass_matrix(rhoA, X2[node1], Y2[node1], X2[node2], Y2[node2])
        k_e2[i] = element_stiffness(EA, EI, X2[node1], Y2[node1], X2[node2], Y2[node2])
    M2 = assembly(m_e2, connec2, num_nodes2)
    K2 = assembly(k_e2, connec2, num_nodes2)
    eigval2, eigvec2, omega2, hertz2, reduced_M, reduced_K = free_vibration(M2, K2, local_mesh)
    print "Mesh 3 Frequencies in Hz, first 12", hertz2[0:12]

    u0_at_R2, _, _ = direct_freq_analysis(rhoA, EA, EI, local_mesh, dof, y_offset, reduced_M, reduced_K, Omega_list)
    ### END MESH 3


    ###### MESH 4
    local_mesh = mesh.Mesh4
    connec4, num_nodes4, X4, Y4, e4 = get_mesh_info(local_mesh)
    shape4 = (e4, 6, 6)  # 6x6 mass matrix for e elements
    m_e4 = np.empty(shape4)
    k_e4 = np.empty(shape4)
    for i in range(e4):
        node1 = connec4[i][0]
        node2 = connec4[i][1]
        m_e4[i] = element_mass_matrix(rhoA, X4[node1], Y4[node1], X4[node2], Y4[node2])
        k_e4[i] = element_stiffness(EA, EI, X4[node1], Y4[node1], X4[node2], Y4[node2])
    M4 = assembly(m_e4, connec4, num_nodes4)
    K4 = assembly(k_e4, connec4, num_nodes4)
    eigval4, eigvec4, omega4, hertz4, reduced_M, reduced_K = free_vibration(M4, K4, local_mesh)
    print "Mesh 4 Frequencies in Hz, first 12", hertz4[0:12]
    u0_at_R4, reduced_f, R_dof = direct_freq_analysis(rhoA, EA, EI, local_mesh, dof, y_offset, reduced_M, reduced_K, Omega_list)
    modal_freq_analyis(eigvec4, eigval4, reduced_f, R_dof, Omega_list)

    ### END MESH 4 ###


    # Q2
    plot_natural_freq(hertz, hertz1, hertz2, hertz4)
    # Q3
    plot_excitation_freq(Omega_list, u0_at_R, u0_at_R1, u0_at_R2, u0_at_R4)
    # plot_mode_shape(eigvec, eigvec1, eigvec2, eigvec4) # Looks pretty but is absolutely useless.



"""
Part A
Study the convergence of the first 5 natural frequencies and mode shapes of
the structure when the finite element spatial mesh is refined.
"""

main()
