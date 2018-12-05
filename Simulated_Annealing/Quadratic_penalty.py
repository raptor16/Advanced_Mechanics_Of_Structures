import numpy as np
import truss
from Assignment3 import get_perturbed_values
from Assignment3 import decision
from Assignment3 import schedule
from matplotlib import pyplot as plt

def SA_3(x0, lb, ub, epsilon=2, max_iter=5000, t_start=1000, c=0.98, n=2, rho=1):
    """
    :param x0: initial guess
    :param lb: vector lower bound
    :param ub: vector upper bound
    :param epsilon:  step-size controlling magnitude of perturbation to design variables
    :param max_iter: maximum number of iterations
    :param t_start: starting temperature
    :param c: cooling schedule parameter
    :return:
    """
    # initialize
    iteration = 0
    T = t_start
    x, xopt = x0, x0
    shape = (5000,)
    fopt_graph = np.empty(shape)

    while(iteration < max_iter):
        x_prime = get_perturbed_values(x, lb, ub, epsilon)
        x = delta_E_acceptance(T, x, x_prime, rho)
        fopt_graph[iteration] = calculate_pi_func(x, rho)
        if calculate_pi_func(x, rho) < calculate_pi_func(xopt, rho):
            xopt = x
            fopt = calculate_pi_func(xopt, rho)
        T = schedule(c, iteration, t_start)
        print iteration
        iteration += 1

    # fopt = calculate_pi_func(xopt)
    return xopt, fopt, fopt_graph


def delta_E_acceptance(T, x, x_prime, rho):
    """
    Generates a new x value
    :param T:
    :param x:
    :param x_prime:
    :param n:
    :return:
    """
    delta_E = compute_delta_E(x, x_prime, rho)

    if delta_E > 0:
        # accept move
        ret = x_prime
    else:
        if decision(delta_E, T):
            ret = x_prime
        else:
            ret = x
    return ret


def compute_delta_E(x_areas, x_areas_prime, rho):
    delta_E = calculate_pi_func(x_areas, rho) - calculate_pi_func(x_areas_prime, rho)
    return delta_E


def is_violate_constraints_stress(sigmas):
    sigma_max = 270 * 10**6 # Newton per Meter Square
    a_large_number = 100
    return np.greater(sigmas, sigma_max) * a_large_number


#################### TRUSS   ##########################

def get_weight(areas):
    """
    Gets the objective function with the added ONE penalty
    :param areas: np array of areas
    :return: total weight scalar of the truss
    """
    lengths = get_lengths()
    density = 2700 # 6061 Aluminium alloy in kg/m3
    weights = areas * lengths * density
    # stress = get_stress(areas)
    # weights = weights + is_violate_constraints_stress(stress)
    total_weight = np.sum(weights)

    return total_weight


def get_stress(list_of_A):
    """
    Gets the stress from the truss solver in A1
    :param list_of_A: list of areas
    :return: stress in N/m^2
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

    displacements = [0, 0, np.NAN, np.NAN, np.NAN, np.NAN, 0, 0, np.NAN, np.NAN, np.NAN,
                     np.NAN]  # Use NAN for unknown boundary conditions
    forces = [0, 0, 0, 0, 0, 0, 0, 0, 0, -100, 0, -100]

    stress = truss.main_driver(sample_connec, nodal_coordinates, list_of_E, list_of_A, displacements, forces)
    stress = np.array(stress)

    return stress


def get_lengths():
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

    list_of_L = truss.get_lengths(sample_connec, nodal_coordinates)
    list_of_L = np.array(list_of_L)
    return list_of_L

###########################################################

def quadratic_penalty(sigmas):
    """

    :param sigmas: stress vector
    :return: penalty function
    """
    # h(x) is 0 because there are no equality constraints.
    # g(x) is the inequality constraint g(x) <= 0. In this case, sigma(x) - sigma_max

    sigma_max = 270 * 10**6 # Newton per Meter Square

    zeros = np.zeros(sigmas.shape)
    g_of_x = sigmas - sigma_max

    phi = np.sum(np.power(np.maximum(zeros, g_of_x), 2))

    return phi

"""
def calculate_rho_scalar(iteration, C=0.5, alpha=1.5):
    rho = np.power(C * iteration, alpha)
    return rho
"""


def calculate_pi_func(x_areas, rho):
    alpha = 1.5
    C = 0.5
    # rho = calculate_rho_scalar(iteration, C, alpha)

    sigmas = get_stress(x_areas)
    phi = quadratic_penalty(sigmas)

    obj_func = get_weight(x_areas)
    pi = obj_func + rho * phi
    return pi


def plot_convergence(fopt, c):
    """
    helper function for plotting the convergence of fopt
    :param fopt: y axis being plotted
    :param c: used for labeling the graph, this is the cooling schedule parameter
    :return: Nothing
    """
    num_of_iterations = 5000
    x = np.linspace(1, num_of_iterations, 5000)
    label_text = "c=" + str(c)
    plt.figure()
    plt.plot(x, fopt, label=label_text)
    plt.xlabel("iterations")
    plt.ylabel("optimal f")
    plt.title("10-bar truss total weight vs. Number of Iterations")
    plt.legend()
    plt.show()

    plt.figure()
    plt.semilogy(x, fopt, label=label_text)
    plt.xlabel("iterations")
    plt.ylabel("optimal f")
    plt.title("10-bar truss total weight  vs. Number of Iterations")
    plt.legend()
    plt.show()

if __name__ == '__main__':

    ######################## 10-bar Truss ######################
    n = 10
    # Design variable is cross-sectional area. There are 10 of them -- list_of_A

    #x_areas = np.array([1 * 10 ** -4, 1 * 10 ** -4, 1 * 10 ** -4, 1 * 10 ** -4, 1 * 10 ** -4,
    #          1 * 10 ** -4, 1 * 10 ** -4, 1 * 10 ** -4, 1 * 10 ** -4, 1 * 10 ** -4])
    x_areas_list = [1. * 10 ** -4] * 10
    x_areas = np.array(x_areas_list)
    # Minimize the weight of the structure subject to stress constraints
    lb_list = [0.] * 10
    ub_list = [0.0001] * 10
    lb = np.array(lb_list)
    ub = np.array(ub_list)
    epsilon = 0.3 * ub[0]
    max_iter = 5000
    t_start = 1000
    c = 0.996
    average_n = 5
    shape = (5000,)
    total_fopt_graph = np.zeros(shape)

    rho = 100
    for i in range(average_n):
        xopt, fopt, fopt_graph = SA_3(x_areas, lb, ub, epsilon, max_iter, t_start, c, n, rho)
        total_fopt_graph += fopt_graph
    average = total_fopt_graph / average_n
    print "xopt", xopt
    print "fopt", fopt
    plot_convergence(average, c)
