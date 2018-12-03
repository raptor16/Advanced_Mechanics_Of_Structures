import move
import numpy as np


def SA(x0, lb, ub, epsilon, max_iter, t_start, c):
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

    xopt, fopt = 0, 0
    return xopt, fopt

def obj_func(x):
    """
    a funtion that returns the objective function given the design variable x
    :param x: design variable
    :return: objective function
    """

    return NotImplemented


def schedule(c, t, t_start = 1000):
    """
    returns the temperature based on the exponential cooling schedule T(t) = t_start where 0<c<1.
    :param c: cooling schedule parameter
    :param t: temperature
    :param t_start: starting temperature
    :return:
    """

    return NotImplemented


if __name__ == '__main__':
    print move.move(np.array([0, 0]), np.array([-1, -1]), np.array([1, 1]), 5)
    #x = np.array([1, 2, 3, 4, 5, 6 ,6, 234 ])
    #lb = np.array([1, 2, 3, 4, 2, 0, 1, 2])
    #ub = np.array([10, 10, 10, 10, 90, 234, 20, 700])
    #move.move(x, lb, ub, 0.1)
    #print move.move(x, lb, ub, 0.1)