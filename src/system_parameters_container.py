class SystemParameters:
    def __init__(self):
        Ds_up: float = 0.0  # (eV^-1 nm^-3) density of states
        Ds_dn: float = 0.0  # (eV^-1 nm^-3) density of states
        Ds: float = 0.0

        penetration_depth: float = 0.0  # (nm)
        E_nt: float = 0.0  # (eV)
        v_fermi: float = 0.0  # (nm fs^1)

        tau_ee_up: float = 0.0  # (fs)
        tau_ee_dn: float = 0.0  # (fs)

        electric_conductivity_driven: float = 0.0  # (nm^-1 fs^-1 eV^-1) written as C in the documentation
        electric_conductivity_diffusive: float = 0.0  # (nm^-1 fs^-1 eV^-1) written as C in the documentation

        length = 0.0  # (nm)