import numpy as np


def Hurkx(electric_field_cm, tunneling_mass_ratio=0.03, temperature=300):
    """
    這是用來計算 Hurkx model 的 Gamma 值，用以比較與 pure SRH 的大小關係。
    :param electric_field: 單位為 (V/cm)，務必注意單位。
    :param tunneling_mass_ratio: 預設為 0.03
    :param temperature: 預設為 300(K)
    :return: 2 * (3pi)^2 * |F|/Fg * exp[(F/Fg)^2]
    """

    k = 1.38e-23  # [J/K]
    q = 1.6e-19  # [C]
    h_bar = 6.626e-34 / (2 * np.pi)
    m0 = 9.11e-31  # [kg]
    f_gamma = (24 * tunneling_mass_ratio * m0 * (k * temperature) ** 3) ** 0.5 / (q * h_bar)
    electric_field = electric_field_cm * 100  # [V/cm] --> [V/m]
    print(electric_field / f_gamma)
    gamma = 2 * np.sqrt(3 * np.pi) * abs(electric_field) / f_gamma * np.exp((electric_field / f_gamma) ** 2)
    return gamma