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

    print "EA", EA
    print "L", L

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
    print "eigval", eigval
    print "eigvec", eigvec
    print "nat freq.", eigval**0.5
    omega = eigval ** 0.5
    # eigenvalues are the sqaure of natural frequencies and eigenvectors are the mode shape

    return eigval, eigvec, omega


def plot_frequency(nat_freq):
    n = nat_freq.shape[0]
    x = np.linspace(0, 1, n)

    plt.ion()
    plt.show()
    plt.plot(x, nat_freq)
    plt.draw()
    plt.pause(5000)

    print "done"


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
    elements = 5
    num_nodes = 6
    l_e = float(L)/float(elements)

    connec_mat = np.array([[0, 1],
                           [1, 2],
                           [2, 3],
                           [3, 4],
                           [4, 5]])
    # since rhoA and l_e is the same for each element then we only need to run calculate_mass_matrix once
    m = calculate_mass_matrix(rhoA, l_e)        # 2x2 mass matrix
    # since EA and l_e is the same for each element then we only need to run calculate_stiffness_matrix once
    k = calculate_stiffness_matrix(EA, l_e)     # 2x2 stiffness matrix

    shape = (elements, 2, 2)
    list_of_M = np.empty(shape)
    list_of_K = np.empty(shape)

    for i in range(elements):
        list_of_M[i] = m
        list_of_K[i] = k

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
    print "M", M
    print "K", K
    eigval, eigvec, omega = eigensolver(M[1:, 1:], K[1:, 1:])
    nat_freq_hertz = omega / (2 * math.pi)
    print "nat_freq_hertz", nat_freq_hertz
    plot_frequency(nat_freq_hertz)


if __name__ == '__main__':
    EA = 6.987 * 10**6   # N
    rhoA = 2.74         # kg/m
    L = 1               # m
    main(rhoA, EA, L)
