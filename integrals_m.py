from molecule import System, Nucleus
import numpy as np
import math
import scipy.special
HeHH = System(nuclei = [Nucleus("He" ,position = [0,0,0]),
                        Nucleus("H", position = [1.4632,0,0])], electrons = 2, xyz = False)

def intis(system):


    # first only for 1 orbital per nucleus
    alphas = np.array([base.alpha for nuc in system.nuclei for base in nuc.basis]) #(n,3)
    N = np.array([base.N for nuc in system.nuclei for base in nuc.basis]) #(n,3)
    centers = np.array([base.center for nuc in system.nuclei for base in nuc.basis]) #(n,3) - (n,1)
    coefs = np.array([base.coef for nuc in system.nuclei for base in nuc.basis]) #(n,3)
    charges = np.array([nuc.charge for nuc in system.nuclei])


    #Sab per nucleus so it's nuc x nuc combinations
    w = N*coefs #(n,3)

    dist = np.linalg.norm((centers[:,None,:] - centers[None,:,:]),axis = 2) #(n,n) KxK combinations
    sum_alpha = alphas[:,None,:,None]+ alphas[None,:,None,:] # (K3)^n combinations HEH 36 -> S11, S12, S21, S21 4x3x3
    mul_alpha = alphas[:,None,:,None]  * alphas[None,:,None,:]
    mul_w = w[:,None,:,None] * w[None,:,None,:]


    dist_expanded = dist[:,:,None, None]
    Sb = (np.pi/sum_alpha)**(3/2)*np.exp(-mul_alpha/sum_alpha*dist_expanded**2)
    S = np.sum(mul_w*(Sb), axis = (2,3))
    T = np.sum(mul_w*mul_alpha/sum_alpha* (3-2*mul_alpha/sum_alpha*dist_expanded**2)*Sb, axis = (2,3))

    center_gauss = ((alphas[:,None,:,None,None]*centers[:,None,None,None,:] + alphas[None,:,None,:,None]*centers[None,:,None,None,:])/(alphas[:,None,:,None,None] +alphas[None,:,None,:,None]))
    t = sum_alpha[:, :, :, :, None] * np.sum((center_gauss[:, :, :, :, None, :] - centers[None, None, None, None, :, :])**2,axis=-1)

    #limit for the error function is 1
    def fo(t):
        Fo = np.ones_like(t)

        mask_t = t > 1e-10

        Fo[mask_t] = 1/2*np.sqrt(np.pi/(t[mask_t]))*scipy.special.erf(np.sqrt(t[mask_t]))
        return Fo

    Vne = -np.sum(mul_w[..., None]*(2*np.pi/sum_alpha[...,None])*charges[None,None,None,None,:]*np.exp(-mul_alpha[...,None]/sum_alpha[...,None]*dist_expanded[...,None]**2)*fo(t), axis = (2,3,-1))

    #2 electron integrals 11/22, 12/22, 11/12, 12/12 22/11, 22/12
    # 8 axes
    # axes: (mu, nu, p, q, lambda, sigma, r, s, xyz)
    #
    mul_w_expanded1 = mul_w[None,None,None,None,...] # (1,1,1,1,K,K,3,3)
    mul_w_expanded2 = mul_w[...,None,None,None,None] # (K,K,3,3,1,1,1,1)
    mul_w_expanded_mul = mul_w_expanded1*mul_w_expanded2 #(K,K,3,3,K,K,3,3)


    sum_alpha_1 = sum_alpha[...,None,None,None,None]
    sum_alpha_2 = sum_alpha[None,None,None,None,...]

    mul_alpha_1 = mul_alpha[...,None,None,None,None]
    mul_alpha_2 = mul_alpha[None,None,None,None,...]

    kab = np.exp(-mul_alpha_1 / sum_alpha_1 * dist_expanded[..., None, None, None, None]**2)
    kcd = np.exp(-mul_alpha_2 / sum_alpha_2 * dist_expanded[None, None, None, None, ...]**2)
    #double electron integral
    d = np.sum(mul_w_expanded_mul * 2 * np.pi**2.5 / (sum_alpha_1 * sum_alpha_2 * np.sqrt(sum_alpha_1 + sum_alpha_2)) * kab * kcd * fo(sum_alpha_1 * sum_alpha_2 / (sum_alpha_1 + sum_alpha_2) * np.sum((center_gauss[:, :, :, :, None, None, None, None, :] - center_gauss[None, None, None, None, :, :, :, :, :])**2, axis=-1)), axis=(2, 3, 6, 7))

    s, U = np.linalg.eigh(S)
    X = U@np.diag(s**-0.5)@U.T
    P = np.zeros_like(S)
    charg_ab = charges[:,None]*charges[None,:]
    Vnn = np.zeros_like(S)
    mask_v = dist != 0
    Vnn[mask_v] = 1

    Vnn[mask_v] = (charg_ab[mask_v]/dist[mask_v])
    Vnne = np.sum(np.tril(Vnn)) # pretty elegant but useless haha

    return Vnne, d, S, T, Vne, X, P


