import numpy in np


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


def one_pass_penalty(f):
    if(is_constriant_violated):
        return f + a_large_number
    else:
        return f


def is_constriant_violated():
    return NotImplemented