import move
import numpy as np
import random

def SA(x0, lb, ub, epsilon=3, max_iter=5000, t_start=1000, c=0.99, n=2):
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
    random.seed(42)
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


if __name__ == '__main__':
    n = 2
    sample_x = np.array([1., 1.])
    epsilon = 0.2 * 10
    lb = np.array([0., 0.])
    ub = np.array([10., 10.])
    t_start = 1000
    c = 0.997
    a_large_number = 100
    max_iter = 5000
    xopt, fopt = SA(sample_x, lb, ub, epsilon, max_iter, t_start, c, n)
    print "xopt", xopt
    print "fopt", fopt
