import math

import matplotlib.animation as animation

from src.material import Material
from src.phyiscs import *

material_up = Material()
material_dn = Material()

sigma_a = 1
sigma_b = 10

Ds_tot = 5.0
Ds_majority = Ds_tot * 0.5
Ds_minority = Ds_tot * 0.5

material_up.Ds_up = Ds_majority
material_up.Ds_dn = Ds_minority
material_up.conductivity_up = sigma_a * material_up.Ds_up / Ds_tot
material_up.conductivity_dn = sigma_a * material_up.Ds_dn / Ds_tot
material_up.tau_sf = 10.0

material_dn.Ds_up = Ds_minority
material_dn.Ds_dn = Ds_majority
material_dn.conductivity_up = sigma_b * material_dn.Ds_up / Ds_tot
material_dn.conductivity_dn = sigma_b * material_dn.Ds_dn / Ds_tot
material_dn.tau_sf = 10.0

system = single_layer(material_up, 5.0)
system.add_layer(material_dn, 5.0, 100.0)

system.e_external = -1.0  # negative for electrons moving toward +z

dt = 0.0001  # (fs)
system.dt = dt

num_frames: int = round(50.0 / dt)

mu0_up_lists = []  # (eV)
mu0_dn_lists = []  # (eV)

t: float = 0.0

for i in range(num_frames):
    print("progress: {} %".format(100 * i / num_frames))

    mu0_up_lists.append(system.mu0List_up)
    mu0_dn_lists.append(system.mu0List_dn)

    system.step()
    t += dt

    system.mu0List_up[0] = 0
    system.mu0List_dn[0] = 0
    system.mu0List_up[system.num_slices-1] = 0
    system.mu0List_dn[system.num_slices-1] = 0

    if i % 25 == 1:
        system.plot2()