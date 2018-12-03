import numpy as np
import random

def move (x, lb, ub, epsilon):
    """
    Randomly perturns the design variable vector x such that the perturbed vector satisfies the bound constraints.

    All 'vectors' are np.array of shape (n, )

    :param x: design variable vector
    :param lb: vector containing lower bounds on the design variables
    :param ub: vector contatining upper bounds on the design variables
    :param epsilon: parameter controlling the magnitude of the pertubation. It is recommended that this parameter is
    set to a value between 0.1 and 0.3 if the design variables are normalized to [0, 1]. Try some typical values and see
    what the impact this parameter makes on the convergence trends.
    :return: z: the perturbed design variable vector satisfying the bound constraints
    """
    # typecast
    x, lb, ub = x.astype(float), lb.astype(float), ub.astype(float)

    n = x.shape[0]
    z = np.empty(n, dtype=float)
    flag = 0
    while(flag == 0):
        rand = random.uniform(0, 1)
        print rand
        ind = int(np.ceil(rand * n) -1)
        z = x.copy()
        z[ind] = x[ind] + epsilon * (-1 + rand * 2)
        print z
        if(z[ind] < lb[ind] or z[ind] > ub[ind]):
            flag = 0
        else:
            flag = 1
    return z
