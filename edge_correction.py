import numpy as np
import sympy.physics.wigner as wg
import fourPCF_estimator as fe
from sympy.physics.wigner import wigner_3j, wigner_9j

def l_modes(odd_mode=True, l_max=5):
    modes = []
    for l1 in range(l_max+1):
        for l2 in range(l_max+1):  
            for l3 in range(abs(l1-l2), min(l1+l2, l_max)+1):
                if (odd_mode == False) or ((l1+l2+l3)%2 == 1):
                    modes.append([l1, l2, l3])
    return np.array(modes)

def D(l, L, lp):
      return np.sqrt((2*l + 1)*(2*L + 1)*(2*lp + 1))

def wigner_9j_prod(l_modes):
    num_modes = l_modes.shape[0]
    prod = np.zeros((num_modes, num_modes, num_modes))

    for i, j, k in zip(range(num_modes), range(num_modes), range(num_modes)):
        l1, l2, l3 = l_modes[i][0], l_modes[i][1], l_modes[i][2]
        L1, L2, L3 = l_modes[j][0], l_modes[j][1], l_modes[j][2]
        lp1, lp2, lp3 = l_modes[k][0], l_modes[k][1], l_modes[k][2]
        prod[i, j, k] = (-1)**(lp1 + lp2 + lp3)/(4*np.pi)**(3/2)*D(l1, L1, lp1)*D(l2, L2, lp2)*D(l3, L3, lp3)  \
            *np.float64(wigner_3j(l1, L1, lp1, 0, 0, 0)) \
            *np.float64(wigner_3j(l2, L2, lp2, 0, 0, 0)) \
            *np.float64(wigner_3j(l3, L3, lp3, 0, 0, 0)) \
            *np.float64(wigner_9j(l1, L1, lp1, l2, L2, lp2, l3, L3, lp3))
    return prod

def R_frac(Random_catalog, Random_catalog_weights, bins, space_vol, l_modes):
    num_modes = l_modes.shape[0]
    frac = np.zeros(num_modes)
    R000 = fe.estimator(0, 0, 0, Random_catalog, Random_catalog_weights, bins, None, None, 0, space_vol)

    for i in range(num_modes):
        Ri = fe.estimator(l_modes[i][0], l_modes[i][1], l_modes[i][2], Random_catalog, Random_catalog_weights, bins, None, None, 0, space_vol)
        frac[i] = Ri/R000
    return frac

def coupling_matrix(l_modes):
    num_modes = l_modes.shape[0]
    M = np.zeros((num_modes, num_modes))
    prod = wigner_9j_prod(l_modes)
    frac = R_frac(l_modes)
    
    for i, j, k in zip(range(num_modes), range(num_modes), range(num_modes)):
        M[i][j] += frac[k]*prod[i][j][k]
    return M