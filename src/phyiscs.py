import math
import random

import matplotlib.pyplot as plt

# FUNDAMENTAL UNITS:
# length: nm
# temperature: K
# energy: eV
# time: fs
# charge: e, the elementary charge.
#     The charge of an electron is -1.0 e, but please use the electron_charge constant for clarity

# DERIVED UNITS
# electric field: (eV nm^-1 e^-1)
# voltage: (eV e^-1)

# physical constants:
kB: float = 0.00008617333262145  # Boltzmann constant (eV K^-1)
alpha: float = 1.645  # related to heat capacity of thermal electrons (1)
epsilon_0: float = 0.055263  # vacuum permeability (e^2 eV^-1 nm^-1)
electron_charge = -1.0  # (e)

# system parameters:
Ds_up: float = 80.0  # (eV^-1 nm^-3) density of states
Ds_dn: float = 20.0  # (eV^-1 nm^-3) density of states
Ds = Ds_up + Ds_dn
diffusivity: float = 0.025  # (nm^2 fs^-1) constant of diffusion of thermal electrons
penetration_depth: float = 15.0  # (nm)
E_nt: float = 1.0  # (eV)
v_fermi: float = 1.0  # (nm fs^1)
tau_ee_up: float = 200.0  # (fs) previously 200.0
tau_ee_dn: float = 10.0  # (fs), previously 10.0

# the values of the parameters below should both be around 90, taken from the bulk conductivity of Ni,
# for the driven transport this would lead to numerical instability.
# The largest value that approximately achieves instant screening is chosen.
electric_conductivity_driven: float = 0.1  # (nm^-1 fs^-1 eV^-1) written as C in the documentation
electric_conductivity_diffusive: float = 90.0  # (nm^-1 fs^-1 eV^-1) written as C in the documentation


# simulation parameters:
electrons_per_packet: float = 0.001  # (nm^-2), 0.0001 for animations


class ExcitedElectron:
    def __init__(self):
        self.z: float = 0.0  # (nm)
        self.vz: float = 0.0  # (nm)
        self.is_up = True


class ElectronState:
    def __init__(self):
        self.sliceLength: float = 0  # (nm)
        self.num_slices: int = 0  # (1)

        self.mu0List_up: list = []  # (eV)
        self.mu0List_dn: list = []  # (eV)

        self.excitedList: list = []

        self.spin_current: list = []  # (nm^-2 fs^-1)

        self.accumulated_energy: float = 0.0  # (eV nm^-2)

    def excited_density(self):
        excited_per_slice = [0] * self.num_slices
        for i in range(len(self.excitedList)):
            n = math.floor(self.excitedList[i].z / self.sliceLength)
            excited_per_slice[n] += 1

        excited_density = list(map((lambda x: x * electrons_per_packet / self.sliceLength), excited_per_slice))

        return excited_density  # (nm^-3)

    def excited_density_up(self):
        excited_per_slice = [0] * self.num_slices
        for i in range(len(self.excitedList)):
            if self.excitedList[i].is_up:
                n = math.floor(self.excitedList[i].z / self.sliceLength)
                excited_per_slice[n] += 1

        excited_density_up = list(map((lambda x: x * electrons_per_packet / self.sliceLength), excited_per_slice))

        return excited_density_up  # (nm^-3)

    def excited_density_dn(self):
        excited_per_slice = [0] * self.num_slices
        for i in range(len(self.excitedList)):
            if not self.excitedList[i].is_up:
                n = math.floor(self.excitedList[i].z / self.sliceLength)
                excited_per_slice[n] += 1

        excited_density_dn = list(map((lambda x: x * electrons_per_packet / self.sliceLength), excited_per_slice))

        return excited_density_dn  # (nm^-3)

    def electron_density(self):
        excited_density = self.excited_density()

        electron_density = [0.0] * self.num_slices
        for i in range(self.num_slices):
            electron_density[i] = (Ds_up * self.mu0List_up[i] + Ds_dn * self.mu0List_dn[i]) + excited_density[i]

        return electron_density  # (nm^-3)

    def charge_density(self):
        electron_density = self.electron_density()
        charge_density = list(map((lambda x: x * electron_charge), electron_density))
        return charge_density  # (e nm^-3)

    def electric_field_intersite(self):
        # electric field on the interface between slice [n] and [n+1]
        charge_density = self.charge_density()

        e_x = [0.0] * (self.num_slices - 1)

        e_x[0] = charge_density[0] * self.sliceLength / epsilon_0
        for i in range(1, self.num_slices - 1):
            e_x[i] = e_x[i - 1] + charge_density[i] * self.sliceLength / epsilon_0

        return e_x  # (eV nm^-1 e^-1)

    def potential(self):
        electric_field_is = self.electric_field_intersite()
        potential = [0.0] * self.num_slices

        for i in range(self.num_slices-1):
            potential[i + 1] = potential[i] - electric_field_is[i] * self.sliceLength

        return potential  # (eV e^-1)

    def electrochemical_potential_up(self):
        potential = self.potential()
        ec_potential = [0.0] * self.num_slices

        for i in range(self.num_slices):
            ec_potential[i] = self.mu0List_up[i] + electron_charge * potential[i]

        return ec_potential  # (eV)

    def electrochemical_potential_dn(self):
        potential = self.potential()
        ec_potential = [0.0] * self.num_slices

        for i in range(self.num_slices):
            ec_potential[i] = self.mu0List_dn[i] + electron_charge * potential[i]

        return ec_potential  # (eV)

    def plot(self):
        steps = [x * self.sliceLength for x in range(0, self.num_slices)]
        excited_density_up = self.excited_density_up()
        excited_density_dn = self.excited_density_dn()

        plt.figure(figsize=(9, 9))

        plt.subplot(211)
        plt.plot(steps, self.mu0List_up, steps, list(map(lambda x: x / Ds_up, excited_density_up)))
        plt.ylabel("mu_0_up (eV)")
        plt.xlabel("z (nm)")
        plt.axis([0, self.num_slices * self.sliceLength, -0.1, 0.1])

        plt.subplot(212)
        plt.plot(steps, self.mu0List_dn, steps, list(map(lambda x: x / Ds_dn, excited_density_dn)))
        plt.ylabel("mu_0_dn (eV)")
        plt.xlabel("z (nm)")
        plt.axis([0, self.num_slices * self.sliceLength, -0.1, 0.1])

        plt.show()

    def extra_electrons_up(self) -> float:
        # calculates the number of excess electrons to check if conservation laws are met

        num_excited = 0.0
        for i in range(len(self.excitedList)):
            if self.excitedList[i].is_up:
                num_excited += 1

        mu_accumulated = 0.0
        for i in range(self.num_slices):
            mu_accumulated += self.mu0List_up[i]

        return num_excited * electrons_per_packet + self.sliceLength * mu_accumulated * Ds_up  # (nm^-2)

    def extra_electrons_dn(self) -> float:
        # calculates the number of excess electrons to check if conservation laws are met

        num_excited = 0.0
        for i in range(len(self.excitedList)):
            if not self.excitedList[i].is_up:
                num_excited += 1

        mu_accumulated = 0.0
        for i in range(self.num_slices):
            mu_accumulated += self.mu0List_dn[i]

        return num_excited * electrons_per_packet + self.sliceLength * mu_accumulated * Ds_dn  # (nm^-2)

    def extra_electrons(self) -> float:
        return self.extra_electrons_up() + self.extra_electrons_dn()

    def advance(self, time: float, power: float):
        dt: float = 0.5
        accumulator = self

        for _i in range(math.floor(time / dt)):
            accumulator = accumulator.step(dt, power)

        return accumulator

    def step(self, dt: float, fluence: float):
        # for brevity
        num_slices = self.num_slices
        sliceLength = self.sliceLength

        # make a new ElectronState for writing the result into
        result = ElectronState()
        result.sliceLength = sliceLength
        result.num_slices = num_slices
        result.mu0List_up = self.mu0List_up.copy()
        result.mu0List_dn = self.mu0List_dn.copy()
        result.excitedList = []
        result.spin_current = [0.0] * num_slices
        result.accumulated_energy = self.accumulated_energy

        # TRANSPORT OF THERMALISED ELECTRONS ---------------------------------------------------------------------------
        # The value in self.muList are used to calculate how many electrons to transport,
        # conservation of energy uses the values in result.muList .
        # Because for every slice transport to the left is evaluated before transport to the right,
        # this may break symmetry slightly but is a compromise to achieve conservation of energy.

        electric_field_is = self.electric_field_intersite()

        for i in range(0, num_slices-1):
            # both types of transport in units of (nm^-3)

            # up
            transport_diffusive_up = (Ds_up / Ds) * \
                dt * electric_conductivity_diffusive * (self.mu0List_up[i] - self.mu0List_up[i + 1]) / (sliceLength * sliceLength)
            transport_driven_up = (Ds_up / Ds) * \
                electric_conductivity_driven * electron_charge * electric_field_is[i] * dt / sliceLength

            transport_up = transport_diffusive_up + transport_driven_up  # (nm^-3)

            result.mu0List_up[i] = result.mu0List_up[i] - transport_up / Ds_up
            result.mu0List_up[i + 1] = result.mu0List_up[i + 1] + transport_up / Ds_up

            spin_current_up = transport_up * sliceLength / dt

            result.spin_current[i] += 0.5 * spin_current_up
            result.spin_current[i + 1] += 0.5 * spin_current_up

            # dn
            transport_diffusive_dn = (Ds_dn / Ds) * \
                dt * electric_conductivity_diffusive * (self.mu0List_dn[i] - self.mu0List_dn[i + 1]) / (sliceLength * sliceLength)
            transport_driven_dn = (Ds_dn / Ds) * \
                electric_conductivity_driven * electron_charge * electric_field_is[i] * dt / sliceLength

            transport_dn = transport_diffusive_dn + transport_driven_dn  # (nm^-3)

            result.mu0List_dn[i] = result.mu0List_dn[i] - transport_dn / Ds_dn
            result.mu0List_dn[i + 1] = result.mu0List_dn[i + 1] + transport_dn / Ds_dn

            spin_current_dn = -transport_dn * sliceLength / dt

            result.spin_current[i] += 0.5 * spin_current_dn
            result.spin_current[i + 1] += 0.5 * spin_current_dn

        # BALLISTIC TRANSPORT OF NON-THERMAL ELECTRONS -----------------------------------------------------------------
        excited_list = []
        for i in range(len(self.excitedList)):
            current_electron = self.excitedList[i]
            new_electron = ExcitedElectron()

            new_electron.z = current_electron.z + current_electron.vz
            new_electron.vz = current_electron.vz
            new_electron.is_up = current_electron.is_up

            initial_slice: int = math.floor(current_electron.z / self.sliceLength)

            if new_electron.z < 0.0:
                new_electron.z = -new_electron.z
                new_electron.vz = -new_electron.vz
            if new_electron.z > num_slices * self.sliceLength:
                new_electron.z = 2 * num_slices * self.sliceLength - new_electron.z
                new_electron.vz = -new_electron.vz

            final_slice: int = math.floor(new_electron.z / self.sliceLength)

            # contribute to spin current
            # the sign is taken care of by the factor spin / (final_slice - initial_slice)
            for n in range(min(initial_slice, final_slice), max(initial_slice, final_slice)):
                if current_electron.is_up:
                    spin = 1
                else:
                    spin = -1
                self.spin_current[n] += spin * electrons_per_packet / (final_slice - initial_slice) / dt

            excited_list.append(new_electron)

        # DECAY OF NON-THERMAL ELECTRONS -------------------------------------------------------------------------------
        p_decay_up = 1.0 - math.exp(-dt / tau_ee_up)  # probability of a given nonequilibrium electron thermalising
        p_decay_dn = 1.0 - math.exp(-dt / tau_ee_dn)  # probability of a given nonequilibrium electron thermalising
        for i in range(len(excited_list)):
            if (excited_list[i].is_up and random.random() < p_decay_up) or \
                    ((not excited_list[i].is_up) and random.random() < p_decay_dn):
                n: int = math.floor(excited_list[i].z / self.sliceLength)
                if excited_list[i].is_up:
                    result.mu0List_up[n] += electrons_per_packet / (Ds_up * self.sliceLength)
                else:
                    result.mu0List_dn[n] += electrons_per_packet / (Ds_dn * self.sliceLength)
            else:
                result.excitedList.append(excited_list[i])

        # EXCITATION OF NON-THERMAL ELECTRONS --------------------------------------------------------------------------
        for i in range(num_slices):

            p_i = fluence * (
                    math.exp(-i * sliceLength / penetration_depth) -
                    math.exp(-(i+1) * sliceLength / penetration_depth)
            )  # (eV fs^-1 nm^-2)

            # for simplicity assuming mu approximately 0
            num_excitations_i = p_i * dt / E_nt / electrons_per_packet  # (1)

            num_excitations_i_up = round(num_excitations_i * Ds_up / Ds)
            num_excitations_i_dn = round(num_excitations_i * Ds_dn / Ds)

            result.mu0List_up[i] -= num_excitations_i_up * electrons_per_packet / (Ds_up * self.sliceLength)
            result.mu0List_dn[i] -= num_excitations_i_dn * electrons_per_packet / (Ds_dn * self.sliceLength)

            result.accumulated_energy += (num_excitations_i_up + num_excitations_i_dn) * electrons_per_packet * E_nt

            for j in range(num_excitations_i_up):
                new_electron = ExcitedElectron()

                new_electron.z = self.sliceLength * (i + random.random())
                new_electron.vz = v_fermi * (2.0 * random.random() - 1)
                new_electron.is_up = True

                result.excitedList.append(new_electron)

            for j in range(num_excitations_i_dn):
                new_electron = ExcitedElectron()

                new_electron.z = self.sliceLength * (i + random.random())
                new_electron.vz = v_fermi * (2.0 * random.random() - 1)
                new_electron.is_up = False

                result.excitedList.append(new_electron)

        return result


def make_equilibrium(temp: float, h: float, num_slices: int) -> ElectronState:
    result = ElectronState()

    result.sliceLength = h
    result.num_slices = num_slices

    result.mu0List_up = [0.0] * num_slices
    result.mu0List_dn = [0.0] * num_slices

    result.excitedList = []

    result.spin_current = [0.0] * num_slices

    result.accumulated_energy = Ds * alpha * kB * kB * temp * temp * num_slices * h

    return result
