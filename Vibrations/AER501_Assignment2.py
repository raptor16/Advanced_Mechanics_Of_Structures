import element_mats as em
import mesh

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

def assembly():
    """
    Assembles the global stiffness matrix and mass matrices
    :return:
    """

def free_vibration():
    """
    Question 2: Calculates natural frequencies and mode shapes
    :return:
    """

def direct_freq_analysis():
    """
    Calculates frequency response using the full-order dynamic stiffness matrix
    :return:
    """

def modal_freq_analyis():
    """
    Calculates frequency response using modal analysis approach.
    :return:
    """


def main():
    L = 1               # unit: m
    EI = 1.286 * 10e4   # Nm^2
    EA = 6.987 * 10e6   # N
    rhoA = 2.74           # kg/m
    x1, y1, x2, y2 = 0, 0, 0, 0

    element_stiffness(EA, EI, x1, y1, x2, y2)
    element_mass_matrix(rhoA, x1, y1, x2, y2)

    assembly()
    # incorporates BCs
    free_vibration()
    direct_freq_analysis()
    modal_freq_analyis()


"""
Part A
Study the convergence of the rst 5 natural frequencies and mode shapes of
the structure when the finite element spatial mesh is rened.
"""

main()
