from integrals_m import intis
import numpy as np
from molecule import System, Nucleus
import matplotlib.pyplot as plt
HeHH = System(nuclei=[Nucleus("He", position=[0, 0, 0]),
                      Nucleus("H",  position=[1.4632, 0, 0])
                      ], electrons=2, xyz=False)
H2_far = System(nuclei=[Nucleus("H", position=[0, 0, 0]),
                        Nucleus("H", position=[5.0, 0, 0])], electrons=2, xyz=False)



Vnne, Vdu, S, T, Vne, X, P = intis(H2_far)

np.set_printoptions(precision=4, suppress=True, linewidth=120)
req = True


max_iter = 100

#print(Vdu) #(K,K,K,K) -> (u,v,s,l)
N = 1

E = []
P = np.zeros_like(S)
E_old = 0.0
for it in range(max_iter):


    n_occ = HeHH.electrons // 2



    trans_Vdu = np.transpose(Vdu,(0,3,2,1)) # -> (u,l,s,v)
    diff = Vdu-1/2*trans_Vdu
    G = np.sum(P[None,None,:,:]*diff,axis=(2,3)) #(K,K)

    H = T + Vne
    F = G + H
    F_=X.T@F@X
    eps, C_prime = np.linalg.eigh(F_)
    C = X@C_prime


    E_el = 0.5 * np.sum(P * (H + F))
    E_tot = E_el + Vnne
    print('E_tot:', E_tot)

    C_occ = C[:, :n_occ]
    P_new = 2 * np.sum(C_occ[:, None, :] * C_occ[None, :, :], axis=2)

    if np.max(np.abs(P_new - P)) < 1e-6 and abs(E_el - E_old) < 1e-8:
        break
    P, E_old = P_new, E_el
    E.append(E_el)

print(E)
plt.scatter(range(len(E)), E)
plt.show()


