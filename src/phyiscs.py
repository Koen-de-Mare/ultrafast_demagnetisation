import math
import random

import matplotlib.pyplot as plt

# UNITS:
# length: nm
# temperature: K
# energy: eV
# time: fs

Ds: float = 10000.0     # density of states, (eV^-1 nm^-3)
diffusivity: float = 0.01  # (nm^2 fs^-1)
penetration_depth: float = 1.0  # (nm)
E_nt: float = 1.0  # (eV)
v_fermi: float = 0.1  # (nm/fs)
tau_ee: float = 50.0  # (fs)
alpha: float = 1.645  # related to heat capacity of thermal electrons (1)

class ExcitedElectron:
    def __init__(self):
        self.z: float = 0.0
        self.vz: float = 0.0


class ElectronState:
    def __init__(self):
        self.sliceLength: float = 0
        self.num_slices: int = 0

        self.muList: list = []
        self.TList: list = []

        self.excitedList: list = []

    def excited_density(self):
        excited_per_slice = [0] * self.num_slices
        for i in range(len(self.excitedList)):
            n = math.floor(self.excitedList[i].z / self.sliceLength)
            excited_per_slice[n] += 1

        excited_density = list(map((lambda x: x / (Ds * self.sliceLength)), excited_per_slice))

        return excited_density

    def plot(self):
        steps = [x * self.sliceLength for x in range(0, self.num_slices)]
        excited_density = self.excited_density()

        plt.plot(steps, self.muList, steps, excited_density)
        plt.ylabel("mu")
        plt.xlabel("z")
        plt.axis([0, self.num_slices * self.sliceLength, -1, 1])

        plt.show()

    def clone_empty(self):  # -> ElectronState
        result = ElectronState()

        num_slices: int = self.num_slices

        result.sliceLength = self.sliceLength
        result.num_slices = num_slices

        result.muList = [0.0] * num_slices
        result.TList = [0.0] * num_slices

        result.excitedList = []

        return result

    def advance(self, time: object, power: object):
        dt: float = 0.5
        accumulator = self

        for _i in range(math.floor(time / dt)):
            accumulator = accumulator.step(dt, power)

        return accumulator

    def step(self, dt: float, power):
        result: ElectronState = self.clone_empty()

        num_slices = self.num_slices

        # diffusive transport of thermalised electrons -----------------------------------------------------------------
        # implements the following differential equation:
        # d_mu / d_t = diffusivity * d^2_mu / d_x^2
        for i in range(1, num_slices - 1):
            d2mu_dx2: float = \
                (self.muList[i-1] + self.muList[i+1] - 2 * self.muList[i]) / (self.sliceLength * self.sliceLength)
            result.muList[i] = self.muList[i] + dt * diffusivity * d2mu_dx2

        # first and last slices as special cases
        # equation derived from local conservation of electrons
        result.muList[0] = \
            self.muList[0] + \
            dt * diffusivity * (self.muList[1] - self.muList[0]) / (self.sliceLength * self.sliceLength)
        result.muList[num_slices - 1] = \
            self.muList[num_slices - 1] + \
            dt * diffusivity * (self.muList[num_slices - 2] - self.muList[num_slices - 1]) / \
            (self.sliceLength * self.sliceLength)

        # ballistic transport of nonthermal electrons ------------------------------------------------------------------
        excited_list = []
        for i in range(len(self.excitedList)):
            current_electron = self.excitedList[i]
            new_electron = ExcitedElectron()

            new_electron.z = current_electron.z + current_electron.vz
            new_electron.vz = current_electron.vz

            if new_electron.z < 0.0:
                new_electron.z = -new_electron.z
                new_electron.vz = -new_electron.vz
            if new_electron.z > num_slices * self.sliceLength:
                new_electron.z = 2 * num_slices * self.sliceLength - new_electron.z
                new_electron.vz = -new_electron.vz

            excited_list.append(new_electron)

        # decay of non-equilibrium electrons ---------------------------------------------------------------------------
        p_decay = 1.0 - math.exp(-dt / tau_ee)  # probability of a given nonequilibrium electron thermalising
        for i in range(len(excited_list)):
            if random.random() < p_decay:
                n: int = math.floor(excited_list[i].z / self.sliceLength)
                result.muList[n] += 1.0 / (Ds * self.sliceLength)
            else:
                result.excitedList.append(excited_list[i])

        # excitation of non-equilibrium electrons ----------------------------------------------------------------------
        for i in range(num_slices):
            p_i = power * math.exp(- (i + 0.5) * self.sliceLength / penetration_depth)
            num_excitations_i: int = math.trunc(self.sliceLength * p_i / (E_nt - self.muList[i]))

            result.muList[i] -= num_excitations_i / (Ds * self.sliceLength)

            for j in range(num_excitations_i):
                new_electron = ExcitedElectron()

                new_electron.z = self.sliceLength * (i + random.random())
                new_electron.vz = v_fermi * (2.0 * random.random() - 1)

                result.excitedList.append(new_electron)

        return result


def make_equilibrium(temp: float, h: float, num_slices: int) -> ElectronState:
    result = ElectronState()

    result.sliceLength = h
    result.num_slices = num_slices

    result.muList = [0.0] * num_slices
    result.TList = [temp] * num_slices

    result.excitedList = []

    return result
