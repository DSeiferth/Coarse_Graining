import numpy as np
import math

def entropyProduction(matrix, p):
    P1=0.
    P2=0.
    N=len(matrix)
    for i in range(N):
        for j in range(N):
            P1=P1+(matrix[j][i]*p[j]-matrix[i][j]*p[i])*math.log(p[j]/p[i])
            if matrix[i][j]!=0:
                P2=P2+((matrix[j][i]*p[j]-matrix[i][j]*p[i]))*math.log(matrix[j][i]/matrix[i][j])
    P1=0.5*P1
    P2=0.5*P2
    print('P1=dS/dt= '+str(P1))
    print('P2(copling to set of thermodynamic forces)= '+str(P2))
    return P1+P2
    
    
def affinity(matrix, p, state1, state2):
    return np.log(matrix[state1][state2]*p[state1])-np.log(matrix[state2][state1]*p[state2])
