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
e = 1.0  # (e)

target_sliceWidth = 0.05  # (nm)

class SimulationState:
    def __init__(self):
        self.dt = 0.0

        # slice properties
        self.h: list = []
        self.Ds_up: list = []
        self.Ds_dn: list = []
        self.tau_sf: list = []

        # interface properties
        self.alpha_up: list = []
        self.alpha_dn: list = []
        self.beta_up: list = []
        self.beta_dn: list = []

        self.num_slices: int = 0  # (1)

        self.mu0List_up: list = []  # (eV)
        self.mu0List_dn: list = []  # (eV)

    def add_layer(self, material, layer_thickness, C_interface):

        new_layers = math.ceil(layer_thickness / target_sliceWidth)
        sliceWidth = layer_thickness / new_layers

        # slice properties
        self.h.extend([sliceWidth] * new_layers)
        self.Ds_up.extend([material.Ds_up] * new_layers)
        self.Ds_dn.extend([material.Ds_dn] * new_layers)
        self.tau_sf.extend([material.tau_sf] * new_layers)

        # interlayer interface
        self.alpha_up.append(-C_interface)
        self.alpha_dn.append(-C_interface)
        self.beta_up.append(0.0)
        self.beta_dn.append(0.0)

        # bulk interface properties
        self.alpha_up.extend([-material.conductivity_up / sliceWidth] * (new_layers - 1))
        self.alpha_dn.extend([-material.conductivity_dn / sliceWidth] * (new_layers - 1))
        self.beta_up.extend([-material.conductivity_up * e] * (new_layers - 1))
        self.beta_dn.extend([-material.conductivity_dn * e] * (new_layers - 1))

        self.num_slices += new_layers

        self.mu0List_up.extend([0.0] * new_layers)
        self.mu0List_dn.extend([0.0] * new_layers)

        self.e_external: float = 0.0

    def make_steps(self):
        steps = []
        z = 0.0
        for n in range(self.num_slices):
            steps.append(z)
            z += self.h[n]

        return steps

    def plot(self):
        steps = self.make_steps()

        plt.figure(figsize=(9, 9))

        plt.plot(steps, self.mu0List_up, steps, self.mu0List_dn)
        plt.ylabel("mu_0 (eV)")
        plt.xlabel("z (nm)")
        plt.axis([0, steps[self.num_slices-1], -1, 1])

        plt.show()

    def plot2(self):
        steps = self.make_steps()
        zmax = steps[len(steps) - 1]

        fig, [ax1, ax2, ax3] = plt.subplots(3, 1, figsize=(10, 8))

        ax1.set_xlim(0, zmax)
        ax1.set_ylim(-0.2, 0.2)
        ax1.set_ylabel("mu_0 \n (eV)")
        ax1.axhline(y=0, color='b')
        ax1.plot(steps, self.mu0List_up, steps, self.mu0List_dn)

        ax2.set_xlim(0, zmax)
        ax2.set_ylim(-20.0, 20.0)
        ax2.set_ylabel("V \n (eV)")
        ax2.axhline(y=0, color='b')
        ax2.plot(steps, self.potential())

        ax3.set_xlim(0, zmax)
        ax3.set_ylim(-1.0, 1.0)
        ax3.set_ylabel("rho \n (e nm^-3)")
        ax3.axhline(y=0, color='b')
        ax3.plot(steps, self.charge_density())

        plt.xlabel("depth (nm)")

        plt.show()

    def electron_density(self):
        electron_density = [0.0] * self.num_slices
        for n in range(self.num_slices):
            electron_density[n] = (self.Ds_up[n] * self.mu0List_up[n] + self.Ds_dn[n] * self.mu0List_dn[n])

        return electron_density  # (nm^-3)

    def charge_density(self):
        electron_density = self.electron_density()
        charge_density = list(map((lambda x: x * electron_charge), electron_density))
        return charge_density  # (e nm^-3)

    def electric_field(self):
        # electric field on the interface between slice [n] and [n+1]
        rho = self.charge_density()

        e_x = [0.0] * (self.num_slices - 1)

        e_x[0] = rho[0] * self.h[0] / epsilon_0
        for n in range(self.num_slices - 2):
            e_x[n+1] = e_x[n] + self.h[n+1] * rho[n+1] / epsilon_0

        # shifts the electric field such that no external electric field is applied
        correction = - 0.5 * e_x[self.num_slices - 2]
        for i in range(self.num_slices - 1):
            e_x[i] += correction + self.e_external

        return e_x  # (eV nm^-1 e^-1)

    def potential(self):
        electric_field_is = self.electric_field()
        potential = [0.0] * self.num_slices

        for n in range(self.num_slices-1):
            potential[n+1] = potential[n] - 0.5 * (self.h[n] + self.h[n+1]) * electric_field_is[n]

        return potential  # (eV e^-1)

    def electrochemical_potential_up(self):
        potential = self.potential()
        ec_potential = [0.0] * self.num_slices

        for n in range(self.num_slices):
            ec_potential[n] = self.mu0List_up[n] + electron_charge * potential[n]

        return ec_potential  # (eV)

    def electrochemical_potential_dn(self):
        potential = self.potential()
        ec_potential = [0.0] * self.num_slices

        for n in range(self.num_slices):
            ec_potential[n] = self.mu0List_dn[n] + electron_charge * potential[n]

        return ec_potential  # (eV)

    def step(self):
        num_slices = self.num_slices

        electric_field = self.electric_field()

        J_up: list = [0.0] * (num_slices - 1)
        J_dn: list = [0.0] * (num_slices - 1)

        for n in range(num_slices - 1):
            J_up[n] = \
                self.alpha_up[n] * (self.mu0List_up[n + 1] - self.mu0List_up[n]) + \
                self.beta_up[n] * electric_field[n]
            J_dn[n] = \
                self.alpha_dn[n] * (self.mu0List_dn[n + 1] - self.mu0List_dn[n]) + \
                self.beta_dn[n] * electric_field[n]

        dNup_dt = [0.0] * num_slices
        dNdn_dt = [0.0] * num_slices

        for n in range(num_slices):
            dNup_dt[n] += (self.mu0List_dn[n] - self.mu0List_up[n]) / self.tau_sf[n]
            dNdn_dt[n] += (self.mu0List_up[n] - self.mu0List_dn[n]) / self.tau_sf[n]

        for n in range(1, num_slices):
            dNup_dt[n] += J_up[n - 1] / self.h[n]
            dNdn_dt[n] += J_dn[n - 1] / self.h[n]

        for n in range(0, num_slices - 1):
            dNup_dt[n] -= J_up[n] / self.h[n]
            dNdn_dt[n] -= J_dn[n] / self.h[n]

        for n in range(num_slices):
            self.mu0List_up[n] += self.dt * dNup_dt[n] / self.Ds_up[n]
            self.mu0List_dn[n] += self.dt * dNdn_dt[n] / self.Ds_dn[n]

def single_layer(material, layer_thickness):
    result = SimulationState()

    new_layers = math.ceil(layer_thickness / target_sliceWidth)
    sliceWidth = layer_thickness / new_layers

    # slice properties
    result.h.extend([sliceWidth] * new_layers)
    result.Ds_up.extend([material.Ds_up] * new_layers)
    result.Ds_dn.extend([material.Ds_dn] * new_layers)
    result.tau_sf.extend([material.tau_sf] * new_layers)

    # bulk interface properties
    result.alpha_up.extend([-material.conductivity_up / sliceWidth] * (new_layers-1))
    result.alpha_dn.extend([-material.conductivity_dn / sliceWidth] * (new_layers-1))
    result.beta_up.extend([-material.conductivity_up * e] * (new_layers-1))
    result.beta_dn.extend([-material.conductivity_dn * e] * (new_layers-1))

    result.num_slices += new_layers

    result.mu0List_up.extend([0.0] * new_layers)
    result.mu0List_dn.extend([0.0] * new_layers)

    return result