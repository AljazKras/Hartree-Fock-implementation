from molecule import System, Nucleus
import numpy as np


HeH = System(nuclei = [Nucleus("He" ,position = [0,0,0]), Nucleus("H", position = [1.4632,0,0])], electrons = 2, xyz = False)
system = HeH
#overlap

alphas = np.array([base.alpha for nuc in system.nuclei for base in nuc.basis])
Ns = np.array([base.N for nuc in system.nuclei for base in nuc.basis])
centers = np.array([base.center for nuc in system.nuclei for base in nuc.basis])
coefs = np.array([base.coef for nuc in system.nuclei for base in nuc.basis])
Nalpha = alphas*Ns #(3,2)
w = coefs * Ns
m = []

for i in range(len(Nalpha[0])):
    for j in range(len(Nalpha[1])):
        p = alphas[0][i] + alphas[1][j] # (3,)
        d =  np.linalg.norm(centers[0] - centers[1])**2 #(1,)
        m.append((np.pi/p)**1.5 *np.exp(-alphas[0][i]*alphas[1][j] /p * d))#(3,)

l = np.outer(w[0], w[1]).ravel()
print(np.sum(l*np.array(m)))