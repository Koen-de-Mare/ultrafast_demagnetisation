from src.system_parameters_container import SystemParameters

def make_default_system() -> SystemParameters:

    syspar = SystemParameters()

    # system parameters:
    syspar.Ds_up: float = 80.0  # (eV^-1 nm^-3) density of states
    syspar.Ds_dn: float = 20.0  # (eV^-1 nm^-3) density of states
    syspar.Ds = syspar.Ds_up + syspar.Ds_dn

    #diffusivity: float = 0.025  # (nm^2 fs^-1) constant of diffusion of thermal electrons

    syspar.penetration_depth: float = 15.0  # (nm)
    syspar.E_nt: float = 1.0  # (eV)
    syspar.v_fermi: float = 1.0  # (nm fs^1)
    syspar.tau_ee_up: float = 200.0  # (fs)
    syspar.tau_ee_dn: float = 10.0  # (fs)

    # the values of the parameters below should both be around 90, taken from the bulk conductivity of Ni,
    # for the driven transport this would lead to numerical instability.
    # The largest value that approximately achieves instant screening is chosen.
    syspar.electric_conductivity_driven: float = 0.1  # (nm^-1 fs^-1 eV^-1) written as C in the documentation
    syspar.electric_conductivity_diffusive: float = 90.0  # (nm^-1 fs^-1 eV^-1) written as C in the documentation

    syspar.length = 200.0  # (nm)

    return syspar