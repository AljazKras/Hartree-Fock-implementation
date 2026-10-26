from molecule import System, Nucleus
import numpy as np


HeH = System(nuclei = [Nucleus("He" ,position = [0,0,0]), Nucleus("H", position = [1.4632,0,0])], electrons = 2, xyz = False)
HeHH = System(nuclei = [Nucleus("He" ,position = [0,0,0]), Nucleus("H", position = [1.4632,0,0]), Nucleus("H", position = [-1.4632,0,0])], electrons = 2, xyz = False)
system = HeHH
#overlap

alphas = np.array([base.alpha for nuc in system.nuclei for base in nuc.basis])
Ns = np.array([base.N for nuc in system.nuclei for base in nuc.basis])
centers = np.array([base.center for nuc in system.nuclei for base in nuc.basis])
coefs = np.array([base.coef for nuc in system.nuclei for base in nuc.basis])
Nalpha = alphas*Ns #(3,2)

bases = np.array([base for nuc in system.nuclei for base in nuc.basis])




def overlap(center, Nsa, alpha, coef):
    w = coef * Nsa
    m = []
    for i in range(len(Ns[0])):
        for j in range(len(Ns[1])):
            p = alpha[0][i] + alpha[1][j]
            d = np.linalg.norm(center[0] - center[1]) ** 2  # (1,)
            m.append((np.pi/p)** 1.5* np.exp(-alpha[0][i] * alpha[1][j]/p*d))

    l = np.outer(w[0], w[1]).ravel()
    return np.sum(l*np.array(m))

k =[]
for i in range(3):
    for j in range(3):
        k.append(overlap(np.vstack([centers[i], centers[j]]), np.vstack([Ns[i], Ns[j]]), np.vstack([alphas[i], alphas[j]]), np.vstack([coefs[i], coefs[j]])))

print(np.array(k).reshape(3,3))