import numpy as np
import physics as phys
import utils

q = 1.6e-19  # [C]
me = 9.11e-31  # [kg]
hbar = 1.054e-34  # [J-s]
eps_InP = 12.5 * 8.85e-14  # [F/cm]
eps_InGaAs = 13.9 * 8.85e-14  # [F/cm] In 0.53 Ga 0.47 As
Eg_InP = 1.35 * q  # [J]
Eg_InGaAs = 0.742 * q  # [J]  # TCAD 應該是 0.7428


def Jt_E_InGaAs(F, ND):
    # 前置參數（來自 Jt_InGaAs）
    mr = 1
    ratio = 0.06
    alpha = 1
    gamma = 1  # [Edx ~ E * (W * gamma)] 之前不知為何設定為 0.5，目前覺得這個修正 gamma 應該不需存在。
    F = F * 100  # 把單位從 V/cm 轉成 V/m

    # 其他參數
    me = 9.11e-31  # [kg]
    mc = 0.04 * me  # [kg]
    mv = 0.04 * me  # [kg]
    # meff = 2 * mc * mv / (mc + mv) * mr  # [kg]
    meff = 0.04 * mr * me
    # Eg = 0.718 * q  # [J]  [From TCAD]
    Eg = Eg_InGaAs  # [J] [https://www.batop.de/information/Eg_InGaAs.html#]
    w = ratio * eps_InGaAs / (q * ND)
    TunnelingCurrent = (2 * meff / Eg) ** 0.5 * (
                q ** 3 * F ** (alpha + 1) * w * gamma / ((2 * np.pi) ** 3 * hbar ** 2)) * \
                       np.exp(- np.pi / (4 * q * hbar * F) * (2 * meff * Eg ** 3) ** 0.5)
    return TunnelingCurrent * 1e-4  # [A / cm^2]



def Jt_InGaAs(V, ND, alpha, mr, ratio):
    # (J.J. Liou 1980)
    gamma = 1  # [Edx ~ E * (W * gamma)] 之前不知為何設定為 0.5，目前覺得這個修正 gamma 應該不需存在。
    me = 9.11e-31  # [kg]
    mc = 0.04 * me  # [kg]
    mv = 0.04 * me  # [kg]
    # meff = 2 * mc * mv / (mc + mv) * mr  # [kg]
    meff = 0.04 * mr * me
    # Eg = 0.718 * q  # [J]  [From TCAD]
    Eg = Eg_InGaAs  # [J] [https://www.batop.de/information/Eg_InGaAs.html#]
    F = (2 * q * V * ND / eps_InGaAs) ** 0.5 * 100  # [V/m]
    w = ratio * (2 * eps_InGaAs * V / (q * ND)) ** 0.5 / 100  # [m]
    hbar = 1.054e-34  # [J-s]
    #print('A: %.3e,  TCAD: %.3e ' % ((2 * meff / Eg) ** 0.5 * q ** 2 / ((2 * np.pi) ** 3 * hbar ** 2), 7.271e19))
    #print('B: %.3e,  TCAD: %.3e' % (np.pi / (4 * q * hbar) * (2 * meff * Eg ** 3) ** 0.5, 5.14e6))

    A_TCAD = 7.271e19 * 1e4  # [m-2s-1V-2]
    B_TCAD = 5.14e6 * 100  # [V/m]

    TunnelingCurrent = (2 * meff / Eg) ** 0.5 * (q ** 3 * F ** (alpha + 1) * w * gamma / ((2 * np.pi) ** 3 * hbar ** 2)) * \
                       np.exp(- np.pi / (4 * q * hbar * F) * (2 * meff * Eg ** 3) ** 0.5)

    TunnelingCurrent_TCAD = A_TCAD * q * w * F ** (alpha + 1) * np.exp(- B_TCAD / F) * 1e-4  # [A/cm2]
    return TunnelingCurrent * 1e-4  # [A / cm^2]


def Jt_InP(V, ND, alpha, mr, ratio):
    # (J.J. Liou 1980)
    gamma = 1  # [Edx ~ E * (W * gamma)] 之前不知為何設定為 0.5，目前覺得這個修正 gamma 應該不需存在。
    me = 9.11e-31  # [kg]
    mc = 0.1149 * me  # [kg]  正確是 0.1149，但 Ando 似乎是 0.065
    mv = 0.1149 * me  # [kg]
    #meff = 0.1149 * mr * me
    meff = 2 * mc * mv / (mc + mv)
    # Eg = 0.718 * q  # [J]  [From TCAD]
    Eg = Eg_InP  # [J] [https://www.batop.de/information/Eg_InGaAs.html#]
    F = (2 * q * V * ND / eps_InP) ** 0.5 * 100  # [V/m]
    w = ratio * (2 * eps_InP * V / (q * ND)) ** 0.5 / 100  # [m]
    hbar = 1.054e-34  # [J-s]
    #print((2 * meff / Eg) ** 0.5 * q ** 2 / ((2 * np.pi) ** 3 * hbar ** 2) * 0.4)
    TunnelingCurrent = (2 * meff / Eg) ** 0.5 * (q ** 3 * F ** (alpha + 1) * w * gamma / ((2 * np.pi) ** 3 * hbar ** 2)) * \
                       np.exp(- np.pi / (4 * q * hbar * F) * (2 * meff * Eg ** 3) ** 0.5)
    return TunnelingCurrent * 1e-4  # [A / cm^2]


def G_BTB_InGaAs(E_Vcm, T, mr):
    if type(E_Vcm) is np.ndarray:
        G = []
        for F in E_Vcm:
            if F == 0:
                G.append(0)
            else:
                E = F * 100
                meff = 0.04 * mr * me
                # TCAD: A = 7.271e19 [cm-2s-1V-2]
                A = (2 * meff / (q * phys.Eg_InGaAs(T))) ** 0.5 * q ** 2 / ((2 * np.pi) ** 3 * hbar ** 2)  # [m-2s-1V-2]
                # A = 7.271e19 * 1e4
                # TCAD: B = 5.14e6 [V/cm]
                B = np.pi / (4 * q * hbar) * (2 * meff * (q * phys.Eg_InGaAs(T)) ** 3) ** 0.5  # [V/m]
                # B = 5.14e6 * 100  # [V/m]
                G.append(A * E ** 2 * np.exp(- B / E))
        return np.asarray(G)
    else:
        if E_Vcm == 0:
            return 0
        else:
            E = E_Vcm * 100  # [V/m]
            meff = 0.04 * mr * me
            # TCAD: A = 7.271e19 [cm-2s-1V-2]
            A = (2 * meff / phys.Eg_InGaAs(T)) ** 0.5 * q ** 2 / ((2 * np.pi) ** 3 * hbar ** 2)  # [m-2s-1V-2]
            # A = 7.271e19 * 1e4
            # TCAD: B = 5.14e6 [V/cm]
            B = np.pi / (4 * q * hbar) * (2 * meff * phys.Eg_InGaAs(T) ** 3) ** 0.5  # [V/m]
            # B = 5.14e6 * 100  # [V/M]
            return A * E ** 2 * np.exp(- B / E)  # E[V/m] ---> G = [m-3s-1]


def G_BTB_InP(E_Vcm, T, m_ratio):
    if type(E_Vcm) is np.ndarray:
        G = []
        for F in E_Vcm:
            if F == 0:
                G.append(0)
            else:
                E = F * 100
                meff = m_ratio * me  # 0.1149 * me
                # TCAD: A = 7.271e19 [cm-2s-1V-2]
                A = (2 * meff / (q * phys.Eg_InP(T))) ** 0.5 * q ** 2 / ((2 * np.pi) ** 3 * hbar ** 2)  # [m-2s-1V-2]
                # TCAD: B = 5.14e6 [V/cm]
                B = np.pi / (4 * q * hbar) * (2 * meff * (q * phys.Eg_InP(T)) ** 3) ** 0.5  # [V/m]
                G.append(A * E ** 2 * np.exp(- B / E))  # E[V/m] ---> G = [m-3s-1]
        return np.asarray(G)
    else:
        if E_Vcm == 0:
            return 0
        else:
            E = E_Vcm * 100
            meff = m_ratio * me  # 0.1149 * me
            # TCAD: A = 7.271e19 [cm-2s-1V-2]
            A = (2 * meff / phys.Eg_InP(T)) ** 0.5 * q ** 2 / ((2 * np.pi) ** 3 * hbar ** 2)  # [m-2s-1V-2]
            # TCAD: B = 5.14e6 [V/cm]
            B = np.pi / (4 * q * hbar) * (2 * meff * phys.Eg_InP(T) ** 3) ** 0.5  # [V/m]
            return A * E ** 2 * np.exp(- B / E)  # E[V/m] ---> G = [m-3s-1]


def J_BTB_InP(x_cm_array, E_Vcm_array, T, m_ratio):
    return utils.ydx(x_cm_array, q * G_BTB_InP(E_Vcm_array, T, m_ratio), 0, len(x_cm_array) - 1)


def J_BTB_InGaAs(x_cm_array, E_Vcm_array, T, m_ratio):
    return utils.ydx(x_cm_array, q * G_BTB_InGaAs(E_Vcm_array, T, m_ratio / 0.04), 0, len(x_cm_array) - 1)


def Em_InGaAs(V, ND):
    F = (2 * q * V * ND / eps_InGaAs) ** 0.5  # [V/cm]
    return F


def Em_InP(V, ND):
    F = (2 * q * V * ND / eps_InP) ** 0.5  # [V/cm]
    return F


def W_InGaAs(V, ND):
    w = (2 * eps_InGaAs * V / (q * ND)) ** 0.5
    return w  # [cm]


def W_InP(V, ND):
    w = (2 * eps_InP * V / (q * ND)) ** 0.5
    return w  # [cm]
