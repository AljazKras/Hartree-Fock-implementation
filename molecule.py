import numpy as np

A_B = 1.8897261
A_N = {"H": 1, "He": 2, "Li": 3, "Be": 4, "B": 5,"C": 6, "N": 7, "O": 8, "F": 9, "Ne": 10}
rng = np.random.default_rng()

# not the cleanest code
class Nucleus:
    def __init__(self, element, charge, position: np.array = None):
        self.name = element
        self.charge = charge
        if position is not None:
            self.position = np.array(position, dtype=float)
        else:
            self.position = np.array(rng.normal(0,[1,1,1]), dtype=float)


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







xx = System(nuclei = [Nucleus("He", 2, )], electrons = 4, xyz = True)
print(xx.nuclei[0].position)