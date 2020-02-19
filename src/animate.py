import math

import matplotlib.animation as animation

from src.phyiscs import *

#system = make_equilibrium(300, 1.25, 160)
#system = make_equilibrium(300, 2.5, 80)
system = make_equilibrium(300, 5.0, 40)
#system = make_equilibrium(300, 10.0, 20)

dt: float = 0.5
num_frames: int = round(150 / dt)

include_ballistic = False  # include magnetisation of ballistic electrons in moke

# pulse properties
t_pulse: float = 20.0
pulse_duration: float = 10.0
pulse_energy: float = 30.0  # really energy per area (eV nm^-2)

mu0_up_lists = []  # (eV)
mu0_dn_lists = []  # (eV)
excited_up_lists = []  # (nm^-3)
excited_dn_lists = []  # (nm^-3)
spin_current_lists = []  # (nm^-2 fs^-1)

t: float = 0.0

for i in range(num_frames):
    print("progress: {} %".format(100 * i / num_frames))

    mu0_up_lists.append(system.mu0List_up)
    mu0_dn_lists.append(system.mu0List_dn)
    excited_up_lists.append(system.excited_density_up())
    excited_dn_lists.append(system.excited_density_dn())
    spin_current_lists.append(system.spin_current)

    print("extra up: {}".format(system.extra_electrons_up()))
    print("extra dn: {}".format(system.extra_electrons_up()))
    print("")

    fluence: float = pulse_energy * math.exp(
        0.5 * (t - t_pulse)*(t_pulse - t) / (pulse_duration * pulse_duration)) / (pulse_duration * math.sqrt(2.0 * math.pi))
    # (eV fs^-1 nm^-2)

    system = system.step(dt, fluence)
    t += dt

print(system.extra_electrons())

steps = [x * system.sliceLength for x in range(0, system.num_slices)]

# code based on:
# https://stackoverflow.com/questions/49165233/two-lines-matplotib-animation

fig, (ax1, ax2, ax3, ax4) = plt.subplots(4, 1, figsize=(10, 8))

ax1.set_xlim(0, system.sliceLength * system.num_slices)
ax1.set_ylim(-0.01, 0.01)
ax1.set_ylabel("mu up \n (eV)")
ax1.axhline(y=0, color='b')
line1, = ax1.plot([], [], lw=3)
line2, = ax1.plot([], [], lw=3)

ax2.set_xlim(0, system.sliceLength * system.num_slices)
ax2.set_ylim(-0.01, 0.01)
ax2.set_ylabel("mu dn \n (eV)")
ax2.axhline(y=0, color='b')
line3, = ax2.plot([], [], lw=3)
line4, = ax2.plot([], [], lw=3)

ax3.set_xlim(0, system.sliceLength * system.num_slices)
ax3.set_ylim(-0.5, 0.5)
ax3.set_ylabel("magnetisation \n (nm^-3)")
ax3.axhline(y=0, color='b')
line5, = ax3.plot([], [], lw=3)


ax4.set_xlim(0, system.sliceLength * system.num_slices)
ax4.set_ylim(-0.1, 0.1)
ax4.set_ylabel("spin current \n (nm^-2 fs^-1)")
ax4.axhline(y=0, color='b')
line6, = ax4.plot([], [], lw=3)

plt.xlabel("depth (nm)")


def init():
    line1.set_data([], [])
    line2.set_data([], [])
    line3.set_data([], [])
    line4.set_data([], [])
    line5.set_data([], [])
    line6.set_data([], [])
    return [line1, line2, line3, line4, line5, line6]


def animate(n):
    line1.set_data(steps, mu0_up_lists[n])
    line2.set_data(steps, list(map((lambda x: x / Ds_up), excited_up_lists[n])))
    line3.set_data(steps, mu0_dn_lists[n])
    line4.set_data(steps, list(map((lambda x: x / Ds_dn), excited_dn_lists[n])))

    magnetisation = [0.0] * system.num_slices
    for i in range(system.num_slices):
        magnetisation_thermal = \
            mu0_up_lists[n][i] * Ds_up - mu0_dn_lists[n][i] * Ds_dn
        magnetisation_nonthermal = \
            excited_up_lists[n][i] - excited_dn_lists[n][i]

        magnetisation[i] = magnetisation_thermal
        if include_ballistic:
            magnetisation[i] += magnetisation_nonthermal

    line5.set_data(steps, magnetisation)
    line6.set_data(steps, spin_current_lists[n])
    return [line1, line2, line3, line4, line5, line6]


anim = animation.FuncAnimation(fig, animate, init_func=init, frames=num_frames, interval=dt*100, blit=True)

show = True
if show:
    plt.show()
else:
    writer = animation.FFMpegWriter(fps=30.0)
    anim.save("output/animation.mp4", writer=writer)
