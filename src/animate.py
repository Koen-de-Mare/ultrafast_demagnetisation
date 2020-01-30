import math

import matplotlib.animation as animation

from src.phyiscs import *

system = make_equilibrium(300, 0.25, 20)

num_frames: int = 500
dt: float = 0.5

# pulse properties
t_pulse: float = 20.0
pulse_duration: float = 10.0
pulse_energy: float = 10000.0

mu_lists = []
excited_lists = []

t: float = 0.0

for i in range(num_frames):
    mu_lists.append(system.muList)
    excited_lists.append(system.excited_density())

    power: float = pulse_energy * math.exp(
        0.5 * (t - t_pulse)*(t_pulse - t) / (pulse_duration * pulse_duration)) / (pulse_duration * math.sqrt(2.0 * math.pi))

    system = system.step(dt, power)
    t += dt

print(system.extra_electrons())

steps = [x * system.sliceLength for x in range(0, system.num_slices)]

# code based on:
# https://stackoverflow.com/questions/49165233/two-lines-matplotib-animation

fig = plt.figure()
ax = plt.axes(xlim=(0, system.sliceLength * system.num_slices), ylim=(-0.03, 0.03))
line1, = ax.plot([], [], lw=3)
line2, = ax.plot([], [], lw=3)


def init():
    line1.set_data([], [])
    line2.set_data([], [])
    return [line1, line2]


def animate(n):
    line1.set_data(steps, mu_lists[n])
    line2.set_data(steps, excited_lists[n])
    return [line1, line2]


anim = animation.FuncAnimation(fig, animate, init_func=init, frames=num_frames, interval=50, blit=True)

plt.show()
