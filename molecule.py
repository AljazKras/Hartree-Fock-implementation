import numpy as np

A_B = 1.8897261
A_N = {"H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5,"C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10}
rng = np.random.default_rng()

ZETA_SZABO = {"H": {"1s": 1.24}, "He": {"1s": 2.0925}}
STO3G = {
    "1s": {"alpha": [2.22766, 0.405771, 0.109818],
           "coef":  [0.154329, 0.535328, 0.444635]},
    "2s": {"alpha": [0.994203, 0.231031, 0.0751386],
           "coef":  [-0.0999672, 0.399515, 0.700115]},
    "2p": {"alpha": [0.994203, 0.231031, 0.0751386],
           "coef":  [0.155916, 0.607684, 0.391957]},
}
class Basis:
    def __init__(self, center, element, orbital):

        self.center = np.array(center)
        self.alpha= np.array(STO3G[orbital]["alpha"]) * ZETA_SZABO[element][orbital] ** 2
        self.coef= np.array(STO3G[orbital]["coef"])
        self.N = (2 * self.alpha / np.pi) ** 0.75

# not the cleanest code
class Nucleus:
    def __init__(self, element, charge= None, position: np.ndarray = None):
        self.element = element
        if charge == None:
            self.charge = A_N[element]
        else:
            self.charge = charge




        if position is not None:
            self.position = np.array(position, dtype=float)
        else:
            self.position = np.array(rng.normal(0,[1,1,1]), dtype=float)

        orbitals = ["1s"]
        self.basis = [Basis(self.position, self.element, orbital) for orbital in orbitals]

class System:

    def __init__(self,nuclei, electrons = None, charge = 0, xyz = False):
        if not xyz:
            self.nuclei = list(nuclei)
            c_nuc = sum(nucleus.charge for nucleus in self.nuclei)
            if electrons is None:
                self.electrons = c_nuc - charge
                self.charge = charge
            else:
                self.electrons = electrons
                self.charge = c_nuc - electrons

        elif xyz:
            # it's actually ORCA input
            print("Enter xyz:")
            lines = []
            while True:
                line = input()
                if line.strip() == "*" or line == "":
                    break
                lines.append(line.strip().split())

            self.charge = int(lines[0][2])
            self.electrons = np.sum([A_N[l[0]] for l in lines[1:]], axis = 0) - self.charge
            self.nuclei = [Nucleus(element = line[0], charge = A_N[line[0]], position= np.array(line[1:4], dtype=float) * A_B) for line in lines[1:]]

        if self.electrons % 2 != 0:
            raise ValueError("The number of electrons is not even.")


#xx = System(nuclei = [Nucleus("He", 2, )], electrons = 4, xyz = True)
#print(xx.nuclei[0].position)