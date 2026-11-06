from molecule import System, Nucleus
import numpy as np
import math

HeH = System(nuclei = [Nucleus("He" ,position = [0,0,0]), Nucleus("H", position = [1.4632,0,0])], electrons = 2, xyz = False)
HeHH = System(nuclei = [Nucleus("He" ,position = [0,0,0]),
                        Nucleus("H", position = [1.4632,0,0])], electrons = 2, xyz = False)
system = HeHH
#overlap

alphas = np.array([base.alpha for nuc in system.nuclei for base in nuc.basis])
Ns = np.array([base.N for nuc in system.nuclei for base in nuc.basis])
centers = np.array([base.center for nuc in system.nuclei for base in nuc.basis])
coefs = np.array([base.coef for nuc in system.nuclei for base in nuc.basis])
Nalpha = alphas*Ns #(3,2)

bases = np.array([base for nuc in system.nuclei for base in nuc.basis])



def fo(t):
    return 0.5 *np.sqrt(np.pi/t)*math.erf(np.sqrt(t))

def one_ele(center, Nsa, alpha, coef, nuclei):
    w = coef * Nsa
    m = []
    T = []
    Vne = []

    for i in range(len(alpha[0])):
        for j in range(len(alpha[1])):
            p = alpha[0][i] + alpha[1][j]
            pl = alpha[0][i] * alpha[1][j]
            d = np.linalg.norm(center[0] - center[1]) ** 2  # (1,)
            Sab_p = (np.pi/p)** 1.5* np.exp(-alpha[0][i] * alpha[1][j]/p*d)
            m.append(Sab_p )
            T.append((pl/p*(3-2*pl/p*d))*Sab_p)

            for k in range(len(nuclei)):
                dee = np.sum((center[0] - nuclei[k]) ** 2)
                Vne.append(-2 * np.pi/p*nuclei[k]*np.exp(-pl / p * d)*fo(p * dee))





    l = np.outer(w[0], w[1]).ravel()
    return np.sum(l*np.array(m)), np.sum(l*np.array(T)), np.sum(l*np.array(Vne))


k,T =[],[]



for i in range(len(centers)):
    for j in range(len(centers)):
        Sab, Tab = one_ele(np.vstack([centers[i], centers[j]]), np.vstack([Ns[i], Ns[j]]), np.vstack([alphas[i], alphas[j]]), np.vstack([coefs[i], coefs[j]]), centers)
        k.append(Sab)
        T.append(Tab)

print(np.array(k).reshape(len(centers),len(centers)))
print(np.array(T).reshape(len(centers),len(centers)))