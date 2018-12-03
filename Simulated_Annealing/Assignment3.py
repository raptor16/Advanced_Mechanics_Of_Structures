import move
import numpy as np

if __name__ == '__main__':
    print move.move(np.array([0, 0]), np.array([-1, -1]), np.array([1, 1]), 5)
    #x = np.array([1, 2, 3, 4, 5, 6 ,6, 234 ])
    #lb = np.array([1, 2, 3, 4, 2, 0, 1, 2])
    #ub = np.array([10, 10, 10, 10, 90, 234, 20, 700])
    #move.move(x, lb, ub, 0.1)
    #print move.move(x, lb, ub, 0.1)