import math

import matplotlib.animation as animation

from src.material import Material
from src.phyiscs import *

material_up = Material()
material_dn = Material()

sigma = 1

Ds_tot = 100.0
Ds_majority = Ds_tot * 0.9
Ds_minority = Ds_tot * 0.1

material_up.Ds_up = Ds_majority
material_up.Ds_dn = Ds_minority
material_up.conductivity_up = sigma * material_up.Ds_up / Ds_tot
material_up.conductivity_dn = sigma * material_up.Ds_dn / Ds_tot
material_up.tau_sf = 5.0

material_dn.Ds_up = Ds_minority
material_dn.Ds_dn = Ds_majority
material_dn.conductivity_up = sigma * material_dn.Ds_up / Ds_tot
material_dn.conductivity_dn = sigma * material_dn.Ds_dn / Ds_tot
material_dn.tau_sf = 5.0

system = single_layer(material_up, 10.0)
system.add_layer(material_dn, 10.0, 0.00001)

# pull out of equilibrium
for n in range(20):
    system.mu0List_up[10 + n] = 1.0
    system.mu0List_dn[350 + n] = -1.0


# wipe beta
#for n in range(system.num_slices-1):
#    system.beta_up[n] = 0.0
#    system.beta_dn[n] = 0.0

dt = 0.01  # (fs)
system.dt = dt

num_frames: int = round(10.0 / dt)

mu0_up_lists = []  # (eV)
mu0_dn_lists = []  # (eV)

t: float = 0.0

for i in range(num_frames):
    print("progress: {} %".format(100 * i / num_frames))

    mu0_up_lists.append(system.mu0List_up)
    mu0_dn_lists.append(system.mu0List_dn)

    system.step()
    t += dt

    #if i % 50 == 1:
    #    system.plot()

steps = system.make_steps()
zmax = steps[len(steps) - 1]

# code based on:
# https://stackoverflow.com/questions/49165233/two-lines-matplotib-animation

fig, ax1 = plt.subplots(1, 1, figsize=(10, 8))

ax1.set_xlim(0, zmax)
ax1.set_ylim(-1.0, 1.0)
ax1.set_ylabel("mu \n (eV)")
line1, = ax1.plot([], [], lw=3)

plt.xlabel("depth (nm)")


def init():
    line1.set_data([], [])
    #line2.set_data([], [])
    return [line1]


def animate(n):
    line1.set_data(steps, mu0_up_lists[n])
    #line2.set_data(steps, mu0_dn_lists[n])
    return [line1]


anim = animation.FuncAnimation(fig, animate, init_func=init, frames=num_frames, interval=dt*10, blit=True)

show = True
if show:
    plt.show()
else:
    writer = animation.FFMpegWriter(fps=30.0)
    anim.save("output/animation.mp4", writer=writer)
