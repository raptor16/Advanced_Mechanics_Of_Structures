import numpy as np

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