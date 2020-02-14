from src.phyiscs import *

system = make_equilibrium(300, 5.0, 40)

print("extra_electrons: {}".format(system.extra_electrons()))
system.plot()

system = system.step(1.0, 100.0)
print("extra_electrons: {}".format(system.extra_electrons()))
system.plot()


system = system.advance(50.0, 0.0)
print("extra_electrons: {}".format(system.extra_electrons()))
system.plot()

system = system.advance(50.0, 0.0)
print("extra_electrons: {}".format(system.extra_electrons()))
system.plot()

system = system.advance(50.0, 0.0)
print("extra_electrons: {}".format(system.extra_electrons()))
system.plot()

system = system.advance(50.0, 0.0)
print("extra_electrons: {}".format(system.extra_electrons()))
system.plot()