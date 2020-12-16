import numpy as np


def alpha(electric_field, position, interface, form_InP, temperature=300):
    """
    :param electric_field: [float], rather than [np.darray]
    :param position: the unit should be the same as interface
    :param interface: [x_InP_InGaAsP, x_InGaAsP_InGaAs, x_InGaAs_InP]
    :param form_InP: 'Okuto' or 'vanOverst'
    :param temperature: default value is 300(K)
    :return: ionization coefficient of electron
    """
    # InP / InGaAsP / InGaAs / InP
    # interface = [x_InP_InGaAsP, x_InGaAsP_InGaAs, x_InGaAs_InP]
    if interface[0] < position < interface[1]:  # InGaAsP
        if electric_field != 0:
            return 9.051e5 * np.exp(-1.73287e6 / abs(electric_field))
        else:
            return 0
    elif interface[1] < position < interface[2]:  # InGaAs
        if electric_field != 0:
            return 1.0e9 * np.exp(-3.6e6 / abs(electric_field))
        else:
            return 0
    else:  # InP
        if electric_field != 0:
            if form_InP == 'Okuto':
                a = 0.2885
                b = 9.7885e5
                c = 3.2237e-4
                d = 9.1731e-4
                gamma = 1.052
                delta = 1.8
                T0 = 300
                return a * (1 + c * (temperature - T0)) * \
                       abs(electric_field) ** gamma * \
                       np.exp(- ((b * (1 + d * (temperature - T0))) / abs(electric_field)) ** delta)
            elif form_InP == 'vanOverst':
                if electric_field < 3.8e5:
                    return 1.12e7 * np.exp(-3.11e6 / abs(electric_field))
                else:
                    return 2.93e6 * np.exp(-2.64e6 / abs(electric_field))
            else:
                raise BaseException("Wrong form")
        else:
            return 0


def beta(electric_field, position, interface, form_InP, temperature=300):
    """
    :param electric_field: [float], rather than [np.darray]
    :param position: the unit should be the same as interface
    :param interface: [x_InP_InGaAsP, x_InGaAsP_InGaAs, x_InGaAs_InP]
    :param form_InP: 'Okuto' or 'vanOverst'
    :param temperature: default value is 300(K)
    :return: ionization coefficient of hole
    """
    if interface[0] < position < interface[1]:  # InGaAsP
        if electric_field != 0:
            return 4.755e6 * np.exp(-2.449e6 / abs(electric_field))
        else:
            return 0
    elif interface[1] < position < interface[2]:  # InGaAs
        if electric_field != 0:
            return 1.3e8 * np.exp(-2.7e6 / abs(electric_field))
        else:
            return 0
    else:  # InP
        if electric_field != 0:
            if form_InP == 'Okuto':
                a = 4.6632
                b = 9.1762e+5
                c = 3.7237e-4
                d = 1.0244e-3
                gamma = 0.8623
                delta = 1.8149
                T0 = 300
                return a * (1 + c * (temperature - T0)) * abs(electric_field) ** gamma * np.exp(- ((b * (1 + d * (temperature - T0))) / abs(electric_field)) ** delta)
            elif form_InP == 'vanOverst':
                if electric_field < 3.8e5:
                    return 4.79e6 * np.exp(-2.55e6 / abs(electric_field))
                else:
                    return 1.62e6 * np.exp(-2.11e6 / abs(electric_field))
            else:
                raise BaseException("Wrong form")
        else:
            return 0