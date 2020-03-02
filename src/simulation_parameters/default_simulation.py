from src.simulation_parameters_container import SimulationParameters

def make_default_simulation() -> SimulationParameters:
    simpar = SimulationParameters()

    simpar.electrons_per_packet: float = 0.001  # (nm^-2), 0.0001 for animations

    simpar.dt = 0.5  # (fs)
    simpar.sliceLength = 10.0  # (nm)

    return simpar
