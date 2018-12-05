import move
import numpy as np
import random
import truss

def SA(x0, lb, ub, epsilon=2, max_iter=5000, t_start=1000, c=0.98, n=2):
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

    while(iteration < max_iter):
        x_prime = get_perturbed_values(x, lb, ub, epsilon)
        x = delta_E_acceptance(T, x, x_prime, n)
        if bump(x) < bump(xopt):
            xopt = x
        T = schedule(c, iteration, t_start)
        print iteration
        iteration += 1

    fopt = bump(xopt, n)
    return xopt, fopt


def get_perturbed_values(x, lb, ub, epsilon):
    """
    Wrapper function for move function in move.py
    :param x: design variable
    :param lb: lower bounds
    :param ub: upper bounds
    :param epsilon: step-size controlling magnitude of perturbation to design variabless
    :return: the perturbed x'
    """
    x_prime = move.move(x, lb, ub, epsilon)
    return x_prime


def schedule(c, iteration, t_start = 1000):
    """
    returns the temperature based on the exponential cooling schedule T(t) = t_start where 0<c<1.
    :param c: cooling schedule parameter
    :param t: time (scalar)
    :param t_start: starting temperature nx1 column vector
    :return:
    """
    T = t_start * np.power(c, iteration)

    return T


def delta_E_acceptance(T, x, x_prime, n=2):
    """
    Generates a new x value
    :param T:
    :param x:
    :param x_prime:
    :param n:
    :return:
    """
    delta_E = compute_delta_E(x, x_prime, n)

    if delta_E > 0:
        # accept move
        ret = x_prime
    else:
        if decision(delta_E, T):
            ret = x_prime
        else:
            ret = x
    return ret

#### HELPERS

def compute_delta_E(x, x_prime, n=2):
    delta_E = bump(x, n) - bump(x_prime, n)

    return delta_E


def decision(delta_E, T):
    probability = np.exp(delta_E / T)
    return random.random() < probability


def is_violate_constraints_bump(x, n=2):
    return np.prod(x) < 0.75 or np.sum(x) > (15 * n / 2)


def bump(x, n=2):
    """
    a funtion that returns the objective function given the design variable x
    :param x: design variable (nx1) column vector
    :return: objective function
    """
    a_large_number = 100

    numerator = - np.abs(np.sum(np.power(np.cos(x), 4)) - 2 * np.prod(np.power(np.cos(x), 2)))
    denominator = np.sqrt(np.sum(np.arange(1, len(x)+1) * np.power(x, 2)))
    if is_violate_constraints_bump(x, n):
        return numerator/denominator + a_large_number

    return numerator/denominator


################### STRESS STUFF ####################


def SA_3(x0, lb, ub, epsilon=2, max_iter=5000, t_start=1000, c=0.98, n=2):
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

    while(iteration < max_iter):
        x_prime = get_perturbed_values(x, lb, ub, epsilon)
        # get_stress
        x = delta_E_acceptance(T, x, x_prime, n)
        if bump(x) < bump(xopt):
            xopt = x
        T = schedule(c, iteration, t_start)
        print iteration
        iteration += 1

    fopt = bump(xopt, n)
    return xopt, fopt


def delta_E_acceptance_3(T, x, x_prime, n=2):
    """
    Generates a new x value
    :param T:
    :param x:
    :param x_prime:
    :param n:
    :return:
    """
    delta_E = compute_delta_E(x, x_prime, n)

    if delta_E > 0:
        # accept move
        ret = x_prime
    else:
        if decision(delta_E, T):
            ret = x_prime
        else:
            ret = x
    return ret


def is_violate_constraints_stress(sigmas, n=2):
    sigma_max = 270 # MegaPascal
    print np.greater(sigmas, sigma_max)
    return np.greater(sigmas, sigma_max)


###################################################################

def get_weight(areas, lengths, density):
    return areas * lengths * density


def get_stress():
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
    displacements = [0, 0, np.NAN, np.NAN, np.NAN, np.NAN, 0, 0, np.NAN, np.NAN, np.NAN,
                     np.NAN]  # Use NAN for unknown boundary conditions
    forces = [0, 0, 0, 0, 0, 0, 0, 0, 0, -100, 0, -100]

    stress = truss.main_driver(sample_connec, nodal_coordinates, list_of_E, list_of_A, displacements, forces)
    stress = np.array(stress)
    return stress

###########################################################


if __name__ == '__main__':

    ######################## Bump ######################
    """
    n = 2
    sample_x = np.array([1, 1])
    epsilon = 2
    lb = np.array([0, 0])
    ub = np.array([10, 10])
    t_start = 1000
    c = 0.99
    a_large_number = 100
    max_iter = 5000
    xopt, fopt = SA(sample_x, lb, ub, epsilon, max_iter, t_start, c, n)
    print "xopt", xopt
    print "fopt", fopt
    """

    ######################## 10-bar Truss ######################
    n = 2 # dimensions
    # Design variable is cross-sectional area. There are 10 of them -- list_of_A
    x_areas = np.array([1 * 10 ** -4, 1 * 10 ** -4, 1 * 10 ** -4, 1 * 10 ** -4, 1 * 10 ** -4,
              1 * 10 ** -4, 1 * 10 ** -4, 1 * 10 ** -4, 1 * 10 ** -4, 1 * 10 ** -4])
    # Minimize the weight of the structure subject to stress constraints


    ############## TRUSS STUFF  ###############



    ######### GET WEIGHT ###########
    list_of_L = truss.get_lengths(sample_connec, nodal_coordinates)
    list_of_L = np.array(list_of_L)
    density = 2700 # 6063 Aluminium alloy in kg/m
    weights = get_weight(x_areas, list_of_L, density)


    # Design Variable is Cross Section
    # Objective is Weight ** Weight is just density * cross section * L


    ###### BACK TO ASSIGN 3 ######
    is_violate_constraints_stress(stress, n)


