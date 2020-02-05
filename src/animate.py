import math

import matplotlib.animation as animation

from src.phyiscs import *

#system = make_equilibrium(300, 0.125, 160)
#system = make_equilibrium(300, 0.25, 80)
#system = make_equilibrium(300, 0.5, 40)
system = make_equilibrium(300, 1.0, 20)

dt: float = 0.5
num_frames: int = round(150 / dt)

# pulse properties
t_pulse: float = 20.0
pulse_duration: float = 10.0
pulse_energy: float = 2000.0  # really energy per area (eV nm^-2)

mu_lists = []
excited_lists = []
temperature_lists = []

t: float = 0.0

for i in range(num_frames):
    mu_lists.append(system.muList)
    excited_lists.append(list(map(lambda x: x / Ds, system.excited_density())))
    temperature_lists.append(system.electron_temperature_distribution())

    fluence: float = pulse_energy * math.exp(
        0.5 * (t - t_pulse)*(t_pulse - t) / (pulse_duration * pulse_duration)) / (pulse_duration * math.sqrt(2.0 * math.pi))
    # (eV fs^-1 nm^-2)

    system = system.step(dt, fluence)
    t += dt

print(system.extra_electrons())

steps = [x * system.sliceLength for x in range(0, system.num_slices)]

# code based on:
# https://stackoverflow.com/questions/49165233/two-lines-matplotib-animation

fig, (ax1, ax2) = plt.subplots(2, 1)

ax1.set_xlim(0, system.sliceLength * system.num_slices)
ax1.set_ylim(-0.01, 0.01)
line1, = ax1.plot([], [], lw=3)
line2, = ax1.plot([], [], lw=3)

ax2.set_xlim(0, system.sliceLength * system.num_slices)
ax2.set_ylim(0, 1000)
line3, = ax2.plot([], [], lw=3)

def init():
    line1.set_data([], [])
    line2.set_data([], [])
    line3.set_data([], [])
    return [line1, line2, line3]


def animate(n):
    line1.set_data(steps, mu_lists[n])
    line2.set_data(steps, excited_lists[n])
    line3.set_data(steps, temperature_lists[n])
    return [line1, line2, line3]


anim = animation.FuncAnimation(fig, animate, init_func=init, frames=num_frames, interval=dt*100, blit=True)

show = True
if show:
    plt.show()
else:
    writer = animation.FFMpegWriter(fps=30.0)
    anim.save("output/animation.mp4", writer=writer)
