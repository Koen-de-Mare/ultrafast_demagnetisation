import math

import matplotlib.pyplot as plt

from src.phyiscs import *

#system = make_equilibrium(300, 1.25, 160)
#system = make_equilibrium(300, 2.5, 80)
#system = make_equilibrium(300, 5.0, 40)
system = make_equilibrium(300, 10.0, 20)

dt: float = 0.5
num_frames: int = round(300 / dt)

# pulse properties
t_pulse: float = 20.0
pulse_duration: float = 10.0
pulse_energy: float = 30.0  # really energy per area (eV nm^-2)

moke_depth = 25.0  # (nm)
moke_list = []
t_list = []

t: float = 0.0

for i in range(num_frames):
    print("progress: {} %".format(100 * i / num_frames))

    uplist = system.mu0List_up
    dnlist = system.mu0List_dn

    moke_signal = 0.0  # (nm^-2)
    for j in range(system.num_slices):
        sensitivity = math.exp(-1 * j * system.sliceLength / moke_depth)
        moke_signal += (Ds_up * uplist[j] - Ds_dn * dnlist[j]) * system.sliceLength

    moke_list.append(moke_signal)
    t_list.append(t)

    print("extra up: {}".format(system.extra_electrons_up()))
    print("extra dn: {}".format(system.extra_electrons_up()))
    print("")

    fluence: float = pulse_energy * math.exp(
        0.5 * (t - t_pulse)*(t_pulse - t) / (pulse_duration * pulse_duration)) / (pulse_duration * math.sqrt(2.0 * math.pi))
    # (eV fs^-1 nm^-2)

    system = system.step(dt, fluence)
    t += dt

print(system.extra_electrons())

plt.plot(t_list, moke_list)
plt.xlabel("time (fs)")
plt.ylabel("magetisation (nm^-2)")
plt.show()