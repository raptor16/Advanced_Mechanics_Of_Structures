import numpy as np

def penalty_function(x):
    phi = 1
    if(is_x_valid(x)):
        phi = 0
    return phi


def penalty_parameter(C, iteration_number):
    alpha = 1.4-10
    rho = np.power(C * iteration_number, alpha)
    return rho


def objective_function_update(C, iteration_number, x):
    """
    Updates the original objective bump function
    :param C:
    :param iteration_number:
    :param x:
    :return:
    """
    rho = penalty_parameter(C, iteration_number)
    f = bump(x) + rho * penalty_function(x)
    return f


def is_x_valid(x):
    """
    return if x is valid or not.
    :param x:
    :return:
    """
    return NotImplemented

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

def calculate_rho_scakar(iteration, C=0.5, alpha=1.5):
    rho = np.power(C * iteration, alpha)
    return rho


def calculate_pi_func(obj_func, rho, phi):
    pi = obj_func + rho * phi
    return pi
