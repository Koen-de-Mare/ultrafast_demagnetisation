import math
import random

import matplotlib.pyplot as plt

# UNITS:
# length: nm
# temperature: K
# energy: eV
# time: fs

# NOTES:
# as the model is 1D some peculiarities occur; the unit of number of excited electrons is actually nm^-2

# physical constants:
kB: float = 0.00008617333262145  # Boltzmann constant (eV K^-1)
alpha: float = 1.645  # related to heat capacity of thermal electrons (1)

# system parameters:
Ds: float = 100000.0  # (eV^-1 nm^-3) density of states
diffusivity: float = 0.025  # (nm^2 fs^-1) constant of diffusion of thermal electrons
penetration_depth: float = 0.5  # (nm)
E_nt: float = 1.0  # (eV)
v_fermi: float = 0.1  # (nm fs^1)
tau_ee: float = 200.0  # (fs)
conductivity: float = 0.005  # (eV fs^-1 nm^-1 K^-1) heat conductivity for the thermalised electron system
electric_conductivity: float = 0.5  # (nm^-1 fs^-1) weird units as vacuum permitivity and electron charge are included

# simulation parameters:
diffusive_transport: bool = False
electrons_per_entity: float = 0.1  # (nm^-2)


class ExcitedElectron:
    def __init__(self):
        self.z: float = 0.0  # (nm)
        self.vz: float = 0.0  # (nm)


class ElectronState:
    def __init__(self):
        self.sliceLength: float = 0  # (nm)
        self.num_slices: int = 0  # (1)

        self.muList: list = []  # (eV)
        self.thermalEnergyList: list = []  # (eV nm^-3), note that this is energy density per volume, not per slice

        self.excitedList: list = []

        self.accumulated_energy: float = 0.0  # (eV nm^-2)

    def excited_density(self):
        excited_per_slice = [0] * self.num_slices
        for i in range(len(self.excitedList)):
            n = math.floor(self.excitedList[i].z / self.sliceLength)
            excited_per_slice[n] += 1

        excited_density = list(map((lambda x: x * electrons_per_entity / self.sliceLength), excited_per_slice))

        return excited_density  # (nm^-3)

    def electron_density(self):
        excited_density = self.excited_density()

        electron_density = [0.0] * self.num_slices
        for i in range(self.num_slices):
            electron_density[i] = Ds * self.muList[i] + excited_density[i]

        return electron_density  # (nm^-3)

    def electric_field_lattice_sites(self):
        # note that the units are off
        # elementary charge and vacuum permitivity have not been included
        # as they can just be absorbed into the conductivity.
        # electric field is calculated in the center of a slice

        electron_density = self.electron_density()

        e_x = [0.0] * self.num_slices

        e_x[0] = -0.5 * electron_density[0]

        for i in range(1, self.num_slices):
            e_x[i] = e_x[i - 1] - 0.5 * self.sliceLength * (electron_density[i - 1] + electron_density[i])

        return e_x  # (nm^-2)

    def electric_field_intersite(self):
        # like in electric_field_lattice_sites() the units are strange
        # electric field at edge between slices
        # for this reason the result has a length of (self.num_slices - 1)

        electron_density = self.electron_density()

        e_x = [0.0] * (self.num_slices - 1)

        e_x[0] = -electron_density[0]
        for i in range(1, self.num_slices - 1):
            e_x[i] = e_x[i - 1] - electron_density[i] * self.sliceLength

        return e_x  # (nm^-2)

    def electron_temperature_distribution(self):
        temperatures = [0.0] * self.num_slices

        for i in range(self.num_slices):
            temperatures[i] = math.sqrt(
                max(0.0, (self.thermalEnergyList[i] / Ds - 0.5 * self.muList[i] * self.muList[i]) / (alpha * kB * kB))
            )

        return temperatures  # (K)

    def plot(self):
        steps = [x * self.sliceLength for x in range(0, self.num_slices)]
        excited_density = self.excited_density()

        plt.figure(figsize=(9, 9))

        plt.subplot(211)
        plt.plot(steps, self.muList, steps, list(map(lambda x: x / Ds, excited_density)))
        plt.ylabel("mu (eV)")
        plt.xlabel("z (nm)")
        plt.axis([0, self.num_slices * self.sliceLength, -0.1, 0.1])

        plt.subplot(212)
        plt.plot(steps, self.electron_temperature_distribution())
        plt.ylabel("Te (K)")
        plt.xlabel("z (nm)")
        plt.axis([0, self.num_slices * self.sliceLength, 0, 1000])

        plt.show()

    def energy(self) -> float:
        temperatures = self.electron_temperature_distribution()
        energy_accumulator: float = 0.0

        for i in range(self.num_slices):
            e_i: float = self.sliceLength * Ds * (
                    0.5 * self.muList[i] * self.muList[i] + alpha * kB * kB * temperatures[i] * temperatures[i]
            )
            energy_accumulator += e_i

        nonthermal_energy: float = len(self.excitedList) * electrons_per_entity * E_nt

        return energy_accumulator + nonthermal_energy  # (eV)

    def extra_electrons(self) -> float:
        # calculates the number of excess electrons to check if conservation laws are met

        num_excited = len(self.excitedList) * electrons_per_entity

        mu_accumulated = 0.0
        for i in range(self.num_slices):
            mu_accumulated += self.muList[i]

        return num_excited + self.sliceLength * mu_accumulated * Ds  # (nm^-2)

    def advance(self, time: float, power: float):
        dt: float = 0.5
        accumulator = self

        for _i in range(math.floor(time / dt)):
            accumulator = accumulator.step(dt, power)

        return accumulator

    def step(self, dt: float, power):
        # for brevity
        num_slices = self.num_slices
        sliceLength = self.sliceLength

        # make a new ElectronState for writing the result into
        result = ElectronState()
        result.sliceLength = sliceLength
        result.num_slices = num_slices
        result.muList = self.muList.copy()
        result.thermalEnergyList = self.thermalEnergyList.copy()
        result.excitedList = []
        result.accumulated_energy = self.accumulated_energy

        # check for stability of diffusive transport
        assert (dt < sliceLength * sliceLength / (2.25 * diffusivity))

        # HEAT DIFFUSION -----------------------------------------------------------------------------------------------
        temperatures = self.electron_temperature_distribution()

        for i in range(1, num_slices - 1):
            d2T_dx2: float = \
                (temperatures[i - 1] + temperatures[i + 1] - 2 * temperatures[i]) / (
                        self.sliceLength * self.sliceLength)
            result.thermalEnergyList[i] = self.thermalEnergyList[i] + dt * conductivity * d2T_dx2

        # first and last slices as special cases
        # equations derived from conservation of energy
        result.thermalEnergyList[0] = self.thermalEnergyList[0] + \
            dt * conductivity * (temperatures[1] - temperatures[0]) / \
                (self.sliceLength * self.sliceLength)
        result.thermalEnergyList[num_slices - 1] = self.thermalEnergyList[num_slices - 1] + \
            dt * conductivity * (temperatures[num_slices - 2] - temperatures[num_slices - 1]) / \
                (self.sliceLength * self.sliceLength)

        # TRANSPORT OF THERMALISED ELECTRONS ---------------------------------------------------------------------------
        # implements the following differential equation:
        # d_mu / d_t = diffusivity * d^2_mu / d_x^2
        # The value in self.muList are used to calculate how many electrons to transport,
        # conservation of energy uses the values in result.muList .
        # Because for every slice transport to the left is evaluated before transport to the right,
        # this may break symmetry slightly but is a compromise to achieve conservation of energy.

        if not diffusive_transport:
            electric_field_is = self.electric_field_intersite()

        for i in range(0, num_slices-1):
            #transport in units of (nm^-3)
            if diffusive_transport:
                transport = Ds * dt * diffusivity * (self.muList[i + 1] - self.muList[i]) / (sliceLength * sliceLength)
            else:
                transport = dt * electric_field_is[i] * electric_conductivity

            mu_a_0 = result.muList[i]
            mu_b_0 = result.muList[i + 1]

            E_a_0 = result.thermalEnergyList[i]
            E_b_0 = result.thermalEnergyList[i + 1]

            E_mu_a_0 = Ds * 0.5 * mu_a_0 * mu_a_0
            E_mu_b_0 = Ds * 0.5 * mu_b_0 * mu_b_0

            E_th_a_0 = E_a_0 - E_mu_a_0
            E_th_b_0 = E_b_0 - E_mu_b_0
            assert (E_th_a_0 >= 0.0)
            assert (E_th_b_0 >= 0.0)

            mu_a_1 = mu_a_0 + transport / Ds
            mu_b_1 = mu_b_0 - transport / Ds

            E_mu_a_1 = Ds * 0.5 * mu_a_1 * mu_a_1
            E_mu_b_1 = Ds * 0.5 * mu_b_1 * mu_b_1

            delta_E = E_mu_a_0 + + E_mu_b_0 - E_mu_a_1 - E_mu_b_1
            assert(delta_E >= -0.01)

            E_th_a_1 = E_th_a_0 + 0.5 * delta_E
            E_th_b_1 = E_th_b_0 + 0.5 * delta_E

            E_a_1 = E_mu_a_1 + E_th_a_1
            E_b_1 = E_mu_b_1 + E_th_b_1

            result.muList[i] = mu_a_1
            result.muList[i + 1] = mu_b_1

            result.thermalEnergyList[i] = E_a_1
            result.thermalEnergyList[i + 1] = E_b_1

        # BALLISTIC TRANSPORT OF NON-THERMAL ELECTRONS -----------------------------------------------------------------
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

        # DECAY OF NON-THERMAL ELECTRONS -------------------------------------------------------------------------------
        p_decay = 1.0 - math.exp(-dt / tau_ee)  # probability of a given nonequilibrium electron thermalising
        for i in range(len(excited_list)):
            if random.random() < p_decay:
                n: int = math.floor(excited_list[i].z / self.sliceLength)
                result.thermalEnergyList[n] += electrons_per_entity * (E_nt - result.muList[n]) / self.sliceLength
                result.muList[n] += electrons_per_entity / (Ds * self.sliceLength)
            else:
                result.excitedList.append(excited_list[i])

        # EXCITATION OF NON-THERMAL ELECTRONS --------------------------------------------------------------------------
        for i in range(num_slices):
            p_i = power * math.exp(- (i + 0.5) * self.sliceLength / penetration_depth)  # (eV fs^-1 nm^-3)
            num_excitations_i: int = round(
                self.sliceLength * dt * p_i / (E_nt - self.muList[i]) / electrons_per_entity
            )  # (1)

            result.accumulated_energy += num_excitations_i * electrons_per_entity * (E_nt - self.muList[i])

            for j in range(num_excitations_i):
                new_electron = ExcitedElectron()

                new_electron.z = self.sliceLength * (i + random.random())
                new_electron.vz = v_fermi * (2.0 * random.random() - 1)

                result.thermalEnergyList[i] -= electrons_per_entity * result.muList[i] / self.sliceLength
                result.muList[i] -= electrons_per_entity / (Ds * self.sliceLength)

                result.excitedList.append(new_electron)

        return result


def make_equilibrium(temp: float, h: float, num_slices: int) -> ElectronState:
    result = ElectronState()

    result.sliceLength = h
    result.num_slices = num_slices

    result.muList = [0.0] * num_slices
    result.thermalEnergyList = [Ds * alpha * kB * kB * temp * temp] * num_slices

    result.excitedList = []

    result.accumulated_energy = Ds * alpha * kB * kB * temp * temp * num_slices * h

    return result
