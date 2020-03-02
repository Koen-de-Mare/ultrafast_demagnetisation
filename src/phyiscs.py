import math
import random

import matplotlib.pyplot as plt

from src.simulation_parameters_container import SimulationParameters
from src.system_parameters_container import SystemParameters

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


class ExcitedElectron:
    def __init__(self):
        self.z: float = 0.0  # (nm)
        self.vz: float = 0.0  # (nm)
        self.is_up = True


class SimulationState:
    def __init__(self):
        self.syspar = SystemParameters()
        self.simpar = SimulationParameters()

        self.num_slices: int = 0  # (1)

        self.mu0List_up: list = []  # (eV)
        self.mu0List_dn: list = []  # (eV)

        self.excitedList: list = []

        self.spin_current: list = []  # (nm^-2 fs^-1)

    def excited_density(self):
        excited_per_slice = [0] * self.num_slices
        for i in range(len(self.excitedList)):
            n = math.floor(self.excitedList[i].z / self.simpar.sliceLength)
            excited_per_slice[n] += 1

        excited_density = list(map((lambda x: x * self.simpar.electrons_per_packet / self.simpar.sliceLength), excited_per_slice))

        return excited_density  # (nm^-3)

    def excited_density_up(self):
        excited_per_slice = [0] * self.num_slices
        for i in range(len(self.excitedList)):
            if self.excitedList[i].is_up:
                n = math.floor(self.excitedList[i].z / self.simpar.sliceLength)
                excited_per_slice[n] += 1

        excited_density_up = list(map((lambda x: x * self.simpar.electrons_per_packet / self.simpar.sliceLength), excited_per_slice))

        return excited_density_up  # (nm^-3)

    def excited_density_dn(self):
        excited_per_slice = [0] * self.num_slices
        for i in range(len(self.excitedList)):
            if not self.excitedList[i].is_up:
                n = math.floor(self.excitedList[i].z / self.simpar.sliceLength)
                excited_per_slice[n] += 1

        excited_density_dn = list(map((lambda x: x * self.simpar.electrons_per_packet / self.simpar.sliceLength), excited_per_slice))

        return excited_density_dn  # (nm^-3)

    def electron_density(self):
        excited_density = self.excited_density()

        electron_density = [0.0] * self.num_slices
        for i in range(self.num_slices):
            electron_density[i] = (self.syspar.Ds_up * self.mu0List_up[i] + self.syspar.Ds_dn * self.mu0List_dn[i]) + excited_density[i]

        return electron_density  # (nm^-3)

    def charge_density(self):
        electron_density = self.electron_density()
        charge_density = list(map((lambda x: x * electron_charge), electron_density))
        return charge_density  # (e nm^-3)

    def electric_field_intersite(self):
        # electric field on the interface between slice [n] and [n+1]
        charge_density = self.charge_density()

        e_x = [0.0] * (self.num_slices - 1)

        e_x[0] = charge_density[0] * self.simpar.sliceLength / epsilon_0
        for i in range(1, self.num_slices - 1):
            e_x[i] = e_x[i - 1] + charge_density[i] * self.simpar.sliceLength / epsilon_0

        return e_x  # (eV nm^-1 e^-1)

    def potential(self):
        electric_field_is = self.electric_field_intersite()
        potential = [0.0] * self.num_slices

        for i in range(self.num_slices-1):
            potential[i + 1] = potential[i] - electric_field_is[i] * self.simpar.sliceLength

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
        steps = [x * self.simpar.sliceLength for x in range(0, self.num_slices)]
        excited_density_up = self.excited_density_up()
        excited_density_dn = self.excited_density_dn()

        plt.figure(figsize=(9, 9))

        plt.subplot(211)
        plt.plot(steps, self.mu0List_up, steps, list(map(lambda x: x / self.syspar.Ds_up, excited_density_up)))
        plt.ylabel("mu_0_up (eV)")
        plt.xlabel("z (nm)")
        plt.axis([0, self.num_slices * self.simpar.sliceLength, -0.1, 0.1])

        plt.subplot(212)
        plt.plot(steps, self.mu0List_dn, steps, list(map(lambda x: x / self.syspar.Ds_dn, excited_density_dn)))
        plt.ylabel("mu_0_dn (eV)")
        plt.xlabel("z (nm)")
        plt.axis([0, self.num_slices * self.simpar.sliceLength, -0.1, 0.1])

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

        return num_excited * self.simpar.electrons_per_packet + self.simpar.sliceLength * mu_accumulated * self.syspar.Ds_up  # (nm^-2)

    def extra_electrons_dn(self) -> float:
        # calculates the number of excess electrons to check if conservation laws are met

        num_excited = 0.0
        for i in range(len(self.excitedList)):
            if not self.excitedList[i].is_up:
                num_excited += 1

        mu_accumulated = 0.0
        for i in range(self.num_slices):
            mu_accumulated += self.mu0List_dn[i]

        return num_excited * self.simpar.electrons_per_packet + self.simpar.sliceLength * mu_accumulated * self.syspar.Ds_dn  # (nm^-2)

    def extra_electrons(self) -> float:
        return self.extra_electrons_up() + self.extra_electrons_dn()

    def advance(self, time: float, power: float):
        accumulator = self

        for _i in range(math.floor(time / self.simpar.dt)):
            accumulator = accumulator.step(power)

        return accumulator

    def step(self, fluence: float):
        # for brevity
        num_slices = self.num_slices
        sliceLength = self.simpar.sliceLength

        # make a new SimulationState for writing the result into
        result = SimulationState()
        result.simpar = self.simpar
        result.syspar = self.syspar
        result.num_slices = num_slices
        result.mu0List_up = self.mu0List_up.copy()
        result.mu0List_dn = self.mu0List_dn.copy()
        result.excitedList = []
        result.spin_current = [0.0] * num_slices

        # TRANSPORT OF THERMALISED ELECTRONS ---------------------------------------------------------------------------
        # The value in self.muList are used to calculate how many electrons to transport,
        # conservation of energy uses the values in result.muList .
        # Because for every slice transport to the left is evaluated before transport to the right,
        # this may break symmetry slightly but is a compromise to achieve conservation of energy.

        electric_field_is = self.electric_field_intersite()

        for i in range(0, num_slices-1):
            # both types of transport in units of (nm^-3)

            # up
            transport_diffusive_up = (self.syspar.Ds_up / self.syspar.Ds) * \
                self.simpar.dt * self.syspar.electric_conductivity_diffusive * (self.mu0List_up[i] - self.mu0List_up[i + 1]) / (sliceLength * sliceLength)
            transport_driven_up = (self.syspar.Ds_up / self.syspar.Ds) * \
                self.syspar.electric_conductivity_driven * electron_charge * electric_field_is[i] * self.simpar.dt / sliceLength

            transport_up = transport_diffusive_up + transport_driven_up  # (nm^-3)

            result.mu0List_up[i] = result.mu0List_up[i] - transport_up / self.syspar.Ds_up
            result.mu0List_up[i + 1] = result.mu0List_up[i + 1] + transport_up / self.syspar.Ds_up

            spin_current_up = transport_up * sliceLength / self.simpar.dt

            result.spin_current[i] += 0.5 * spin_current_up
            result.spin_current[i + 1] += 0.5 * spin_current_up

            # dn
            transport_diffusive_dn = (self.syspar.Ds_dn / self.syspar.Ds) * \
                self.simpar.dt * self.syspar.electric_conductivity_diffusive * (self.mu0List_dn[i] - self.mu0List_dn[i + 1]) / (sliceLength * sliceLength)
            transport_driven_dn = (self.syspar.Ds_dn / self.syspar.Ds) * \
                self.syspar.electric_conductivity_driven * electron_charge * electric_field_is[i] * self.simpar.dt / sliceLength

            transport_dn = transport_diffusive_dn + transport_driven_dn  # (nm^-3)

            result.mu0List_dn[i] = result.mu0List_dn[i] - transport_dn / self.syspar.Ds_dn
            result.mu0List_dn[i + 1] = result.mu0List_dn[i + 1] + transport_dn / self.syspar.Ds_dn

            spin_current_dn = -transport_dn * sliceLength / self.simpar.dt

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

            initial_slice: int = math.floor(current_electron.z / sliceLength)

            if new_electron.z < 0.0:
                new_electron.z = -new_electron.z
                new_electron.vz = -new_electron.vz
            if new_electron.z > num_slices * sliceLength:
                new_electron.z = 2 * num_slices * sliceLength - new_electron.z
                new_electron.vz = -new_electron.vz

            final_slice: int = math.floor(new_electron.z / sliceLength)

            # contribute to spin current
            # the sign is taken care of by the factor spin / (final_slice - initial_slice)
            for n in range(min(initial_slice, final_slice), max(initial_slice, final_slice)):
                if current_electron.is_up:
                    spin = 1
                else:
                    spin = -1
                self.spin_current[n] += spin * self.simpar.electrons_per_packet / (final_slice - initial_slice) / self.simpar.dt

            excited_list.append(new_electron)

        # DECAY OF NON-THERMAL ELECTRONS -------------------------------------------------------------------------------
        p_decay_up = 1.0 - math.exp(-self.simpar.dt / self.syspar.tau_ee_up)  # probability of a given nonequilibrium electron thermalising
        p_decay_dn = 1.0 - math.exp(-self.simpar.dt / self.syspar.tau_ee_dn)  # probability of a given nonequilibrium electron thermalising
        for i in range(len(excited_list)):
            if (excited_list[i].is_up and random.random() < p_decay_up) or \
                    ((not excited_list[i].is_up) and random.random() < p_decay_dn):
                n: int = math.floor(excited_list[i].z / sliceLength)
                if excited_list[i].is_up:
                    result.mu0List_up[n] += self.simpar.electrons_per_packet / (self.syspar.Ds_up * sliceLength)
                else:
                    result.mu0List_dn[n] += self.simpar.electrons_per_packet / (self.syspar.Ds_dn * sliceLength)
            else:
                result.excitedList.append(excited_list[i])

        # EXCITATION OF NON-THERMAL ELECTRONS --------------------------------------------------------------------------
        for i in range(num_slices):

            p_i = fluence * (
                    math.exp(-i * sliceLength / self.syspar.penetration_depth) -
                    math.exp(-(i+1) * sliceLength / self.syspar.penetration_depth)
            )  # (eV fs^-1 nm^-2)

            # for simplicity assuming mu approximately 0
            num_excitations_i = p_i * self.simpar.dt / self.syspar.E_nt / self.simpar.electrons_per_packet  # (1)

            num_excitations_i_up = round(num_excitations_i * self.syspar.Ds_up / self.syspar.Ds)
            num_excitations_i_dn = round(num_excitations_i * self.syspar.Ds_dn / self.syspar.Ds)

            result.mu0List_up[i] -= num_excitations_i_up * self.simpar.electrons_per_packet / (self.syspar.Ds_up * sliceLength)
            result.mu0List_dn[i] -= num_excitations_i_dn * self.simpar.electrons_per_packet / (self.syspar.Ds_dn * sliceLength)

            for j in range(num_excitations_i_up):
                new_electron = ExcitedElectron()

                new_electron.z = sliceLength * (i + random.random())
                new_electron.vz = self.syspar.v_fermi * (2.0 * random.random() - 1)
                new_electron.is_up = True

                result.excitedList.append(new_electron)

            for j in range(num_excitations_i_dn):
                new_electron = ExcitedElectron()

                new_electron.z = sliceLength * (i + random.random())
                new_electron.vz = self.syspar.v_fermi * (2.0 * random.random() - 1)
                new_electron.is_up = False

                result.excitedList.append(new_electron)

        return result


def make_equilibrium(system_parameters, simulation_parameters) -> SimulationState:
    num_slices = math.floor(system_parameters.length / simulation_parameters.sliceLength)

    result = SimulationState()

    result.syspar = system_parameters
    result.simpar = simulation_parameters

    result.num_slices = num_slices

    result.mu0List_up = [0.0] * num_slices
    result.mu0List_dn = [0.0] * num_slices

    result.excitedList = []

    result.spin_current = [0.0] * num_slices

    return result
