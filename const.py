import scipy.constants as sp

angstrom = 1e-10 / 4.0 / sp.pi / sp.epsilon_0 / sp.hbar**2 * sp.m_e * sp.e**2

atom_number = {
    1: "H",
    2: "He",
    3: "Li",
    4: "Be",
    5: "B",
    6: "C",
    7: "N",
}