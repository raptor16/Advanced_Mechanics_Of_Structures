import math
import numpy as np
from scipy.linalg import eigh
import matplotlib.pyplot as plt

def calculate_mass_matrix(rhoA, L):
    # convert to floats just in case the inputs aren't floats.
    rhoA, L = float(rhoA), float(L)

    mat = np.array([[2, 1],
                   [1, 2]])
    mat = mat.astype(float)
    M = (rhoA * L / 6) * mat

    return M


def calculate_stiffness_matrix(EA, L):
    # convert to floats just in case the inputs aren't floats.

    EA, L = float(EA), float(L)

    mat = np.array([[1, -1],
                    [-1, 1]])
    mat = mat.astype(float)
    M = (EA / L) * mat

    return M


def assembler(list_of_matrix, connec_mat, num_nodes):
    # read the connec_mat and place the m_e corresponding the e^th index of the connec_mat
    # just remember to sum it up when placing it in the new matrix
    # or better yet, construct a nxn 0 matrix where n is the number of nodes
    # the indices of the mat corresponding in the e^th connec_mat
    # then just sum all those matrices up.

    # make a 3d matrix for a list of 'resized' matrix [len(list_of_matrix), n, n]
    elements = list_of_matrix.shape[0]
    shape = (elements, num_nodes, num_nodes)
    Mat = np.zeros(shape)       # This is a sparse matrix
    for i in range(elements):
        node1 = connec_mat[i][0]
        node2 = connec_mat[i][1]
        Mat[i][node1][node1] = list_of_matrix[i][0][0]
        Mat[i][node1][node2] = list_of_matrix[i][0][1]
        Mat[i][node2][node1] = list_of_matrix[i][1][0]
        Mat[i][node2][node2] = list_of_matrix[i][1][1]
    shape = (num_nodes, num_nodes)
    assembled_Mat = np.zeros(shape)
    for i in range(elements):
        assembled_Mat += Mat[i]

    return assembled_Mat


def eigensolver(M, K):
    # K \phi = \lambda M \phi ; \lamda is the eigenvalues and \phi is the eigenvectors
    eigval, eigvec = eigh(K, M, eigvals_only=False)
    omega = eigval ** 0.5
    # eigenvalues are the sqaure of natural frequencies and eigenvectors are the mode shape

    return eigval, eigvec, omega


def connec_generator(elements, L=1):
    num_nodes = elements + 1
    l_e = float(L) / float(elements)
    shape = (elements, 2)
    connec_mat = np.empty(shape, dtype=int)

    for i in range (connec_mat.shape[0]):
        for j in range (connec_mat.shape[1]):
            connec_mat[i][j] = int(i + j)
    return connec_mat, num_nodes, l_e


def plot_frequency(nat_freq, f1, f2, f3, f4):

    n = nat_freq.shape[0]
    x = np.linspace(0, n, n)

    n1 = f1.shape[0]
    x1 = np.linspace(0, n1, n1)

    n2 = f2.shape[0]
    x2 = np.linspace(0, n2, n2)

    n3 = f3.shape[0]
    x3 = np.linspace(0, n3, n3)

    n4 = f4.shape[0]
    x4 = np.linspace(0, n4, n4)

    plt.ion()
    plt.figure()
    plt.plot(x, nat_freq, label="5 elements")
    plt.plot(x1, f1, label="10 elements")
    plt.plot(x2, f2, label="30 elements")
    plt.plot(x3, f3, label="60 elements")
    plt.plot(x4, f4, label="100 elements")
    plt.legend()
    plt.title("Natural Frequency vs. Frequency Mode")
    plt.xlabel("Frequency Modes")
    plt.ylabel("Natural frequency in Hz")
    plt.pause(2000)
    plt.show()

    print "done"


def plot_mode_shape(elements, eigvec, eigvec1, eigvec2, eigvec3, eigvec4, L=1):

    x = np.linspace(0, 1, elements[0])
    x1 = np.linspace(0, 1, elements[1])
    x2 = np.linspace(0, 1, elements[2])
    x3 = np.linspace(0, 1, elements[3])
    x4 = np.linspace(0, 1, elements[4])
    """
    for j = 1:5
        plot(linspace(0,1,num_elements(i)), eig_vectors(:,j));
    end
    """



    plt.ion()
    plt.figure()
    for i in range (len(elements)):
        plt.plot(x, eigvec[:, i], label="Eigenvector "+str(i))
        plt.title("Mode Shape for N=6")
        plt.xlabel("Length of rod in m")
        plt.ylabel("Shape Function at Node m")
        plt.legend()
    plt.figure()
    for i in range(len(elements)):
        plt.plot(x1, eigvec1[:, i], label="Eigenvector " + str(i))
        plt.title("Mode Shape for N=11")
        plt.xlabel("Length of rod in m")
        plt.ylabel("Shape Function at Node m")
        plt.legend()
    plt.figure()
    for i in range(len(elements)):
        plt.plot(x2, eigvec2[:, i], label="Eigenvector " + str(i))
        plt.title("Mode Shape for N=31")
        plt.xlabel("Length of rod in m")
        plt.ylabel("Shape Function at Node m")
        plt.legend()
    plt.figure()
    for i in range(len(elements)):
        plt.plot(x3, eigvec3[:, i], label="Eigenvector " + str(i))
        plt.title("Mode Shape for N=61")
        plt.xlabel("Length of rod in m")
        plt.ylabel("Shape Function at Node m")
        plt.legend()
    plt.figure()
    for i in range(len(elements)):
        plt.plot(x4, eigvec4[:, i], label="Eigenvector " + str(i))
        plt.title("Mode Shape for N=101")
        plt.xlabel("Length of rod in m")
        plt.ylabel("Shape Function at Node m")
        plt.legend()

    plt.pause(1000)
    plt.show()


def main(rhoA=1, EA=1, L=1.):

    """
    split the rod into 6 nodes i.e. 5 elements of equal length
    0 -- 1 -- 2 -- 3 -- 4 -- 5
      e1   e2   e3   e4   e5

    connec_mat will represent the following table:
    
    element | node 1 | node 2
    e1      | 0      | 1
    e2      | 1      | 2
    e3      | 2      | 3
    e4      | 3      | 4
    e5      | 4      | 5
    """
    # TODO: Use classes for each mesh.

    elements = [5, 10, 30, 60, 100]

    # Mesh 0
    connec_mat, num_nodes, l_e = connec_generator(elements[0], L)

    # Mesh 1
    connec_mat1, num_nodes1, l_e1 = connec_generator(elements[1], L)

    # Mesh 2
    connec_mat2, num_nodes2, l_e2 = connec_generator(elements[2], L)

    # Mesh 3
    connec_mat3, num_nodes3, l_e3 = connec_generator(elements[3], L)

    # Mesh 4
    connec_mat4, num_nodes4, l_e4 = connec_generator(elements[4], L)

    # since rhoA and l_e is the same for each element then we only need to run calculate_mass_matrix once
    m = calculate_mass_matrix(rhoA, l_e)        # 2x2 mass matrix
    # since EA and l_e is the same for each element then we only need to run calculate_stiffness_matrix once
    k = calculate_stiffness_matrix(EA, l_e)     # 2x2 stiffness matrix

    m1 = calculate_mass_matrix(rhoA, l_e1)
    k1 = calculate_stiffness_matrix(EA, l_e1)
    m2 = calculate_mass_matrix(rhoA, l_e2)
    k2 = calculate_stiffness_matrix(EA, l_e2)
    m3 = calculate_mass_matrix(rhoA, l_e3)
    k3 = calculate_stiffness_matrix(EA, l_e3)
    m4 = calculate_mass_matrix(rhoA, l_e4)
    k4 = calculate_stiffness_matrix(EA, l_e4)

    shape = (elements[0], 2, 2)
    list_of_M = np.empty(shape)
    list_of_K = np.empty(shape)

    shape1 = (elements[1], 2, 2)
    shape2 = (elements[2], 2, 2)
    shape3 = (elements[3], 2, 2)
    shape4 = (elements[4], 2, 2)

    list_of_M1 = np.empty(shape1)
    list_of_K1 = np.empty(shape1)
    list_of_M2 = np.empty(shape2)
    list_of_K2 = np.empty(shape2)
    list_of_M3 = np.empty(shape3)
    list_of_K3 = np.empty(shape3)
    list_of_M4 = np.empty(shape4)
    list_of_K4 = np.empty(shape4)

    for i in range(elements[0]):
        list_of_M[i] = m
        list_of_K[i] = k

    for i in range(elements[1]):
        list_of_M1[i] = m1
        list_of_K1[i] = k1

    for i in range(elements[2]):
        list_of_M2[i] = m2
        list_of_K2[i] = k2

    for i in range(elements[3]):
        list_of_M3[i] = m3
        list_of_K3[i] = k3

    for i in range(elements[4]):
        list_of_M4[i] = m4
        list_of_K4[i] = k4

    """
      node1 node2
    [[  a    b  ]   node 1
     [  c    d  ]]  node 2 

     node2 node3
    [[  e    f  ]   node 2
     [  g    h  ]]  node 3 

     node3 node4
    [[  i    j  ]   node 3
     [  k    l  ]]  node 4

     when assembled will looks more like this

     [[  a    b       0    0 ]
      [  c  d + e     f    0 ]
      [  0    g     h + i  k ]
      [  0    0       k    l ]]

      etc. -- in the end you'll get a 6x6

    Because it is just a 1d rod that's been 'chopped' up, we will see a trend where the middle part of the diagonal
    is just going to be a summation of the common nodes which makes sense physically. Now, we can use this pattern
    to make our code for assembly more efficient by just brute forcing this known matrix instead of looping through it.
    But for the sake of scalability and because we need to reuse it anyway for part b, I'm just gonna write scalable 
    code.  

    """

    M = assembler(list_of_M, connec_mat, num_nodes)
    K = assembler(list_of_K, connec_mat, num_nodes)

    M1 = assembler(list_of_M1, connec_mat1, num_nodes1)
    K1 = assembler(list_of_K1, connec_mat1, num_nodes1)
    M2 = assembler(list_of_M2, connec_mat2, num_nodes2)
    K2 = assembler(list_of_K2, connec_mat2, num_nodes2)
    M3 = assembler(list_of_M3, connec_mat3, num_nodes3)
    K3 = assembler(list_of_K3, connec_mat3, num_nodes3)
    M4 = assembler(list_of_M4, connec_mat4, num_nodes4)
    K4 = assembler(list_of_K4, connec_mat4, num_nodes4)

    eigval, eigvec, omega = eigensolver(M[1:, 1:], K[1:, 1:])

    eigval1, eigvec1, omega1 = eigensolver(M1[1:, 1:], K1[1:, 1:])
    eigval2, eigvec2, omega2 = eigensolver(M2[1:, 1:], K2[1:, 1:])
    eigval3, eigvec3, omega3 = eigensolver(M3[1:, 1:], K3[1:, 1:])
    eigval4, eigvec4, omega4 = eigensolver(M4[1:, 1:], K4[1:, 1:])


    nat_freq_hertz = omega / (2 * math.pi)

    nat_freq_hertz1 = omega1 / (2 * math.pi)
    nat_freq_hertz2 = omega2 / (2 * math.pi)
    nat_freq_hertz3 = omega3 / (2 * math.pi)
    nat_freq_hertz4 = omega4 / (2 * math.pi)

    print "nat_freq_hertz first 5", nat_freq_hertz[0:5]
    print "nat_freq_hertz1", nat_freq_hertz1[0:5]
    print "nat_freq_hertz2", nat_freq_hertz2[0:5]
    print "nat_freq_hertz3", nat_freq_hertz3[0:5]
    print "nat_freq_hertz4", nat_freq_hertz4[0:5]

    plot_frequency(nat_freq_hertz, nat_freq_hertz1[0:5], nat_freq_hertz2[0:5], nat_freq_hertz3[0:5], nat_freq_hertz4[0:5])
    plot_mode_shape(elements, eigvec, eigvec1, eigvec2, eigvec3, eigvec4)


if __name__ == '__main__':
    connec_generator(5, 1)
    EA = 6.987 * 10**6   # N
    rhoA = 2.74         # kg/m
    L = 1               # m
    main(rhoA, EA, L)
