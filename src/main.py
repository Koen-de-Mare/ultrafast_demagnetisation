from src.phyiscs import *

system = make_equilibrium(300, 0.25, 25)
system.plot()

system = system.step(1.0, 10000.0)
print(system.extra_electrons())
system.plot()

system = system.advance(25.0, 0.0)
print(system.extra_electrons())
system.plot()

system = system.advance(25.0, 0.0)
print(system.extra_electrons())
system.plot()

system = system.advance(25.0, 0.0)
print(system.extra_electrons())
system.plot()

print("long time")

system = system.advance(2000.0, 0.0)
print(system.extra_electrons())
system.plot()
