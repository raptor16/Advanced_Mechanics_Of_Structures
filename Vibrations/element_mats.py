"""
Pre-made Element Matrix for AER501 Assignment 2 Part B, Ported to Python 2.7 from MATLAB

ElementMassMat.m -> element_mats.mass()
ElementStiffMat.m -> element_mats.stiffness()
"""

from __future__ import division

import numpy as np

def mass(rho, x1, y1, x2, y2):
    """
    This routine computes the element mass matrix for
    a two-dimensional frame element

    function [MM] = ElementMassMat(rho, x1, y1, x2, y2)

    INPUTS
    ------
    rho    - Mass per unit length
    x1, y1 - X and Y coordinates of node 1
    x2, y2 - X and Y coordinates of node 2

    OUTPUT
    ------
    MM     - 6 x 6 element mass matrix
    """
    ELL = np.sqrt((x2-x1)**2 + (y2-y1)**2)

    C=(x2-x1)/ELL
    S=(y2-y1)/ELL

    MM = np.zeros((6,6))

    # Consistent Mass matrix

    MM[0,0]=140*C*C+156*S*S
    MM[0,1]=-16*C*S
    MM[0,2]=-22*ELL*S
    MM[0,3]=70*C*C+54*S*S
    MM[0,4]=-MM[0,1]
    MM[0,5]=13*ELL*S

    MM[1,0]=MM[0,1]
    MM[1,1]=156*C*C+140*S*S
    MM[1,2]=22*ELL*C
    MM[1,3]=-MM[0,1]
    MM[1,4]=54*C*C+70*S*S
    MM[1,5]=-13*ELL*C

    MM[2,0]=MM[0,2]
    MM[2,1]=MM[1,2]
    MM[2,2]=4*ELL*ELL
    MM[2,3]=-MM[0,5]
    MM[2,4]=-MM[1,5]
    MM[2,5]=-3*ELL*ELL

    MM[3,0]=MM[0,3]
    MM[3,1]=MM[1,3]
    MM[3,2]=MM[2,3]
    MM[3,3]=MM[0,0]
    MM[3,4]=MM[0,1]
    MM[3,5]=-MM[0,2]

    MM[4,0]=MM[0,4]
    MM[4,1]=MM[1,4]
    MM[4,2]=MM[2,4]
    MM[4,3]=MM[3,4]
    MM[4,4]=MM[1,1]
    MM[4,5]=-MM[1,2]

    MM[5,0]=MM[0,5]
    MM[5,1]=MM[1,5]
    MM[5,2]=MM[2,5]
    MM[5,3]=MM[3,5]
    MM[5,4]=MM[4,5]
    MM[5,5]=MM[2,2]

    FAC=rho*ELL/420.

    MM = FAC*MM

    return MM


def stiffness(EA, EI, x1, y1, x2, y2):
    """
    This routine computes the element stiffness matrix for
    a two-dimensional frame element

    function [KM] = ElementStiffMat(EA, EI, X1, Y1, X2, Y2);

    INPUTS
    ------
    EA      - axial rigidity
    EI      - flexural rigidity
    X1, Y1  - X and Y coordinates of node 1
    X2, Y2  - X and Y coordinates of node 2

    OUTPUT
    ------
    KM      - 6 x 6 element stiffness matrix
    """

    L = np.sqrt((x2-x1)**2 + (y2-y1)**2)

    C=(x2-x1)/L
    S=(y2-y1)/L

    E1=EA/L
    E2=12*EI/(L*L*L)
    E3=EI/L
    E4=6*EI/(L*L)

    KM = np.zeros((6,6))

    # Define elements of K matrix

    KM[0,0]=C*C*E1+S*S*E2
    KM[3,3]=KM[0,0]
    KM[0,1]=S*C*(E1-E2)
    KM[1,0]=KM[0,1]
    KM[3,4]=KM[0,1]
    KM[4,3]=KM[3,4]
    KM[0,2]=-S*E4
    KM[2,0]=KM[0,2]
    KM[0,5]=KM[0,2]
    KM[5,0]=KM[0,5]
    KM[2,3]=S*E4
    KM[3,2]=KM[2,3]
    KM[3,5]=KM[2,3]
    KM[5,3]=KM[3,5]
    KM[0,3]=-KM[0,0]
    KM[3,0]=KM[0,3]
    KM[0,4]=S*C*(-E1+E2)
    KM[4,0]=KM[0,4]
    KM[1,3]=KM[0,4]
    KM[3,1]=KM[1,3]
    KM[1,1]=S*S*E1+C*C*E2
    KM[4,4]=KM[1,1]
    KM[1,4]=-KM[1,1]
    KM[4,1]=KM[1,4]
    KM[1,2]=C*E4
    KM[2,1]=KM[1,2]
    KM[1,5]=KM[1,2]
    KM[5,1]=KM[1,5]
    KM[2,2]=4*E3
    KM[5,5]=KM[2,2]
    KM[2,4]=-C*E4
    KM[4,2]=KM[2,4]
    KM[4,5]=KM[2,4]
    KM[5,4]=KM[4,5]
    KM[2,5]=2*E3
    KM[5,2]=KM[2,5]

    return KM
