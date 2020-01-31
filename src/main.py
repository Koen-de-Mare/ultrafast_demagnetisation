from src.phyiscs import *

system = make_equilibrium(300, 0.25, 25)
print(system.extra_electrons())
print("extra_electrons: {}".format(system.extra_electrons()))
print("energy: {}".format(system.energy()))
print("input energy: {}".format(system.accumulated_energy))
system.plot()

system = system.step(1.0, 2500.0)
print(system.extra_electrons())
print("extra_electrons: {}".format(system.extra_electrons()))
print("energy: {}".format(system.energy()))
print("input energy: {}".format(system.accumulated_energy))
system.plot()

system = system.advance(25.0, 0.0)
print("extra_electrons: {}".format(system.extra_electrons()))
print("energy: {}".format(system.energy()))
print("input energy: {}".format(system.accumulated_energy))
system.plot()

system = system.advance(25.0, 0.0)
print("extra_electrons: {}".format(system.extra_electrons()))
print("energy: {}".format(system.energy()))
print("input energy: {}".format(system.accumulated_energy))
system.plot()

system = system.advance(25.0, 0.0)
print("extra_electrons: {}".format(system.extra_electrons()))
print("energy: {}".format(system.energy()))
print("input energy: {}".format(system.accumulated_energy))
system.plot()

system = system.advance(25.0, 0.0)
print("extra_electrons: {}".format(system.extra_electrons()))
print("energy: {}".format(system.energy()))
print("input energy: {}".format(system.accumulated_energy))
system.plot()