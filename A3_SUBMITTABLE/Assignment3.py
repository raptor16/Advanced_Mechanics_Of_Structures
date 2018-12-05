import move
import numpy as np
import random
from matplotlib import pyplot as plt

def SA(x0, lb, ub, epsilon=3, max_iter=5000, t_start=1000, c=0.99, n=2):
    """
    :param x0: initial guess
    :param lb: vector lower bound
    :param ub: vector upper bound
    :param epsilon:  step-size controlling magnitude of perturbation to design variables
    :param max_iter: maximum number of iterations
    :param t_start: starting temperature
    :param c: cooling schedule parameter
    :return: xopt is the optimized design variable values, fopt is the optimized objective function values,
        evaluated objective function at x, primarily used for graphing
    """
    # initialize
    iteration = 0
    T = t_start
    x, xopt = x0, x0
    shape = (5000, )
    fopt_graph = np.empty(shape)

    while(iteration < max_iter):
        x_prime = get_perturbed_values(x, lb, ub, epsilon)
        x = delta_E_acceptance(T, x, x_prime, n)
        if not is_violate_constraints_bump(x, n):
            fopt_graph[iteration] = bump(x)
        if bump(x) < bump(xopt):
            xopt = x
        T = schedule(c, iteration, t_start)
        print iteration
        iteration += 1

    fopt = bump(xopt, n)
    return xopt, fopt, fopt_graph


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
    :return: Temperature based on the cooling schedule function used
    """
    T = t_start * np.power(c, iteration)

    return T


def delta_E_acceptance(T, x, x_prime, n=2):
    """
    Generates a x value depending on if a move is accepted or rejected
    :param T: caculated Temperature from cooling schedule function
    :param x: the array of design variables
    :param x_prime: the array of design variables perturbed
    :param n: this is not necessary; it is the dimensions
    :return: the design variable update based on the move
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
    """
    helper for calculating the deltaE
    :param x: design variable
    :param x_prime: design variable perturbed
    :param n: not needed, dimension
    :return:
    """
    delta_E = bump(x, n) - bump(x_prime, n)

    return delta_E


def decision(delta_E, T):
    """
    helper for 'rolling the dice' based on probability
    :param delta_E: deltaE = f(x) - f(x')
    :param T: Temperature from cooling shcedule
    :return: true or false depending on the decision absed on the probability
    """
    probability = np.exp(delta_E / T)
    random.seed(42)
    return random.random() < probability


def is_violate_constraints_bump(x, n=2):
    """
    checks if any of the bump function constraints are being violated
    :param x: design variable
    :param n: dimension, optional
    :return:
    """
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

########################################################################

def plot_convergence(fopt, c):
    num_of_iterations = 5000
    x = np.linspace(1, num_of_iterations, 5000)
    print x
    label_text = "c=" + str(c)
    plt.plot(x, fopt, label=label_text)
    plt.xlabel("iterations")
    plt.ylabel("optimal f")
    plt.title("Bump Function vs. Number of Iterations")
    plt.legend()
    plt.show()

########################################################################

if __name__ == '__main__':
    n = 2
    sample_x = np.array([1., 1.])
    epsilon = 0.2 * 10
    lb = np.array([0., 0.])
    ub = np.array([10., 10.])
    t_start = 1000
    c = 0.996
    a_large_number = 100
    max_iter = 5000

    average_n = 5

    shape = (5000, )
    total_fopt_graph = np.zeros(shape)
    """
    for i in range(average_n):
        xopt, fopt, fopt_graph = SA(sample_x, lb, ub, epsilon, max_iter, t_start, c, n)
        total_fopt_graph += fopt_graph
    average = total_fopt_graph / average_n
    print "xopt", xopt
    print "fopt", fopt
    plot_convergence(average, c)
    """
    c_arr = np.array([0.97, 0.975, 0.98, 0.985, 0.99, 0.995, 0.996, 0.997, 0.998, 0.999])
    fopt_arr = np.empty(c_arr.shape)
    print c_arr.shape
    print fopt_arr.shape
    for i in range(len(c_arr)):
        _, fopt1, _ = SA(sample_x, lb, ub, epsilon, max_iter, t_start, c_arr[i], n)
        fopt_arr[i] = fopt1
    plt.figure()
    plt.plot(c_arr, fopt_arr)
    plt.xlabel("c values")
    plt.ylabel("fopt")
    plt.title("fopt vs c parameter")
    plt.show()