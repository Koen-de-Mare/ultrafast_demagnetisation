from src.phyiscs import *

system = make_equilibrium(300, 0.25, 25)

#system.plot()
system = system.step(1.0, 10000.0)
system.plot()
system = system.advance(25.0)
system.plot()
system = system.advance(25.0)
system.plot()
system = system.advance(25.0)
system.plot()
system = system.advance(25.0)
system.plot()
system = system.advance(25.0)
system.plot()
system = system.advance(25.0)
system.plot()
