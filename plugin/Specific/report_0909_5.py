import numpy as np
import os
import csv
import physics as phys
import Experiment as Exp
import ExpInterface as EI
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import matplotlib.pylab as pylab
import DataAnalysis as Data
import utils
from scipy.optimize import curve_fit
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (20, 9.3),
          'axes.labelsize': 'x-large',
          'axes.titlesize':'x-large',
          'xtick.labelsize':'x-large',
          'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)
plt.rcParams.update({'font.size': 13})

# 這是分析變溫IV的程式。目的是要寫一個適用 TCAD 以及實驗數據的程式。
# 作法是獨立出讀取資料的程式。

kB = 1.38e-23  # [J/k]
e = 1.6e-19  # [C]
me = 9.11e-31  # [kg]
eps_InP = 12.5 * 8.85e-14  # [F/cm]
eps_InGaAs = 13.9 * 8.85e-14  # [F/cm] In 0.53 Ga 0.47 As
eps_InGaAsP = 13.436 * 8.85e-14  # [F/cm] Approximated by In 0.53 Ga 0.47 As 0.65 P 0.35
hbar = 1.054e-34  # [J-s]
Vpt = 22.5  # [V]
Vpt_dec = 25  # [V]
Location = 'data/TAPD.csv'
A = np.pi * (120e-4 ** 2)  # [cm-2]
EtValue = {'InP': -0.1, 'InGaAs': 0.3}
EtY = [-1.15, 1.15]
Temperature = np.arange(240, 340, 10)  # [K]
V0 = 2
Vf = 33
step = 0.1
iterations = (Vf - V0) / step
Voltage = np.asarray([-V0 - step * i for i in range(int(iterations))])
Tmin = 240
Tmax = 330
T_analysis = np.arange(Tmin, Tmax+10, 10)
invT_list = np.arange(1/Tmax, 1/Tmin + 1e-5, 1e-5)
add_V = [-33.5, -34.5, -35.5, -36.5, -37.5, -38.5]
Merge_V = list(Voltage) + list(add_V)
Voltage_InP = [V for V in Merge_V if abs(V) < abs(Vpt)]
Voltage_InGaAs = [V for V in Merge_V if abs(V) >= abs(Vpt_dec)]

# 讀取IV，這裡必須給出 RawIV，不論TCAD還是實驗。
RawIV = dict()
LocationNew = {T: 'data/raw data/T' + str(T) + '/T.csv' for T in Temperature}
fieldnames = ['DataName', 'V1', 'I1', 'I2', 'I11', 'I4']
for T in Temperature:
    tempVoltage = []
    tempCurrent = []
    with open(LocationNew[T], newline='', encoding='utf-8-sig') as csv_file:
        rows = csv.DictReader(csv_file, fieldnames=fieldnames)
        for _index, row in enumerate(rows):
            if row['DataName'] == 'DataValue':
                tempVoltage.append(float(row['V1']))
                tempCurrent.append(float(row['I11']))
        RawIV[T] = Data.CSV(data=None, xlabel=None, ylabel=None, x=tempVoltage, y=tempCurrent)

# 讀取 InP 最大電場與偏壓的分佈
Ef_InP = Data.CSV('data/SRH_off_InP_ElectricField.csv', xlabel='MaxInP_ElectricField', ylabel='MaxInP_ElectricField')
Ef_InGaAs = Data.CSV('data/SRH_off_ElectricField.csv',
                     xlabel='MaxWindow((-3.5:0):(-3.25:3)) ElectricField(IV_n6504_des_6)',
                     ylabel='MaxWindow((-3.5:0):(-3.25:3)) ElectricField(IV_n6504_des_6)')


def Em_InP(V):
    return utils.find(Ef_InP.X, Ef_InP.Y, -abs(V), 'linear')


def Em_InGaAs(V):
    return utils.find(Ef_InGaAs.X, Ef_InGaAs.Y, -abs(V), 'linear')


def dm_InP(E, ND, ND_c, d_mul, d_charge):
    d = E * eps_InP / (e * ND)
    if d <= d_mul:
        return d
    else:
        E2 = E - (e * ND * d_mul) / eps_InP
        d2 = E2 * eps_InP / (e * ND_c)
        if d2 <= d_charge:
            return d_mul + d2
        else:
            return d_mul + d_charge


def dm_InGaAs(E, ND_abs, d_abs):
    ND_g = 2e15
    d_g = 0.12e-4
    d = E * eps_InGaAs / (e * ND_abs)
    if d <= d_abs:
        return d
    else:
        return d_abs


def GammaDx_InP(X, Eti, mt300, C):
    T, Emax_Vcm = X
    alpha = 1
    tp = 0.1
    tn = 0.226
    prefactor = 1

    me = 9.11e-31
    Nc300 = 5.716e17  # [cm-3]
    Nv300 = 1.143e19  # [cm-3]
    tau_p0 = 300e-9  # [s]
    tau_n0 = 2.89e-9  # [s]
    tau_p = tp * tau_p0 * (T / 300) ** alpha
    tau_n = tn * tau_n0 * (T / 300) ** alpha
    ND = 5e16  # [cm-3]
    Ncharge = 7.8e16  # [cm-3]
    d_mul = 0.42e-4  # [cm]
    d_ch = 0.2e-4  # [cm]

    mt = C * (T - 300) + mt300
    ni = np.sqrt(Nc300 * Nv300) * (T / 300) ** 1.5 * np.exp(- e * phys.Eg_InP(T) / (2 * kB * T))
    G_SRH = ni / (2 * np.sqrt(tau_p * tau_n) * np.cosh(e * Eti / (kB * T) + 0.5 * np.log(tau_p / tau_n)))
    if type(mt) is np.ndarray:
        for i, element in enumerate(mt):
            if element < 0:
                mt[i] = 2
    dM = dm_InP(Emax_Vcm, ND, Ncharge, d_mul, d_ch)  # 0.42e-4  # [cm]
    F_Gamma = np.sqrt(24 * (mt * me) * (kB * T) ** 3) / (e * hbar) / 100  # [V/cm]
    E1 = Emax_Vcm
    if dM <= d_mul:
        E2 = E1 - (e * ND * dM) / eps_InP
        d_Gamma_1 = (np.sqrt(3 * np.pi) * eps_InP * F_Gamma) / (e * ND) * \
                    (np.exp((E1 / F_Gamma) ** 2) - np.exp(E2 / F_Gamma ** 2))  # [cm]
        return d_Gamma_1
    else:
        E2 = E1 - (e * ND * d_mul) / eps_InP
        E3 = E2 - (e * Ncharge * (dM - d_mul)) / eps_InP
        d_Gamma_1 = (np.sqrt(3 * np.pi) * eps_InP * F_Gamma) / (e * ND) * \
                    (np.exp((E1 / F_Gamma) ** 2) - np.exp(E2 / F_Gamma ** 2))  # [cm]
        d_Gamma_2 = (np.sqrt(3 * np.pi) * eps_InP * F_Gamma) / (e * Ncharge) * \
                    (np.exp((E2 / F_Gamma) ** 2) - np.exp(E3 / F_Gamma ** 2))  # [cm]
        return d_Gamma_1 + d_Gamma_2


def GammaDx_InGaAs(X, Eti, mt300, C):
    T, Emax_Vcm = X
    prefactor = 1
    tp = 1
    tn = 1
    alpha = 1.5
    me = 9.11e-31
    Nc300 = 2.53956e17  # [cm-3]
    Nv300 = 7.51e18  # [cm-3]
    tau_p0 = 8e-9  # [s]
    tau_n0 = 0.25e-9  # [s]
    tau_p = tp * tau_p0 * (T / 300) ** alpha
    tau_n = tn * tau_n0 * (T / 300) ** alpha
    ND_abs = 7.5e14  # [cm-3]
    d_InGaAs = 3e-4  # [cm]

    ni = np.sqrt(Nc300 * Nv300) * (T / 300) ** 1.5 * np.exp(- e * phys.Eg_InGaAs(T) / (2 * kB * T))
    G_SRH = ni / (2 * np.sqrt(tau_p * tau_n) * np.cosh(e * Eti / (kB * T) + 0.5 * np.log(tau_p / tau_n)))

    dM = dm_InGaAs(Emax_Vcm, ND_abs, d_InGaAs)  # [cm]
    mt = C * (T - 300) + mt300
    F_Gamma = np.sqrt(24 * (mt * me) * (kB * T) ** 3) / (e * hbar) / 100  # [V/cm]
    E1 = Emax_Vcm
    E2 = 0
    d_Gamma = (np.sqrt(3 * np.pi) * eps_InGaAs * F_Gamma) / (e * ND_abs) * (
                np.exp((E1 / F_Gamma) ** 2) - np.exp((E2 / F_Gamma) ** 2))  # [cm]
    return d_Gamma


def GammaDx(X, Eti, mt300, C, V, Vpt, Vpt_dec):
    if abs(V) < abs(Vpt):
        return GammaDx_InP(X, Eti, mt300, C)
    elif abs(V) >= Vpt_dec:
        return GammaDx_InGaAs(X, Eti, mt300, C)
    else:
        return np.nan


def Em(V, Vpt):
    if abs(V) < abs(Vpt):
        return Em_InP(V)
    else:
        return Em_InGaAs(V)


def dm(V, Vpt):
    ND = 5e16  # [cm-3]
    Ncharge = 7.8e16  # [cm-3]
    ND_abs = 7.5e14  # [cm-3]
    d_mul = 0.42e-4  # [cm]
    d_ch = 0.2e-4  # [cm]
    d_InGaAs = 3e-4  # [cm]
    if abs(V) < abs(Vpt):
        return dm_InP(Em(V, Vpt), ND, Ncharge, d_mul, d_ch)
    else:
        return dm_InGaAs(Em(V, Vpt), ND_abs, d_InGaAs)


def mt_InP(x_cm, T, N_list, d_list):
    N_mul, N_ch = N_list
    d_mul_cm, d_charge_cm = d_list
    x_um = x_cm * 1e4
    if x_um < 0.3:
        return 1  # 2.5 * (4.9e-3 * np.exp(x_um / 0.62 * 3.7) + 1 * x_um)
    else:
        return 1  # 2.5 * (4.9e-3 * np.exp(0.3 / 0.62 * 3.7) + 1 * 0.3)


def mt_InGaAs(x_cm, T, N_list, d_list):
    N_abs, = N_list
    d_abs, = d_list
    if x_cm <= d_abs:
        pass
    else:
        pass


def ElectricField_InP(w_cm, x, N_list, d_list):
    N_mul, N_ch = N_list
    d_mul_cm, d_charge_cm = d_list
    d_InP = [dm_InP(Em_InP(element), 5e16, 7.8e16, 0.42e-4, 0.2e-4) for element in Voltage_InP]
    v = utils.find(d_InP, Voltage_InP, w_cm, 'linear', extrapolation='Yes')
    Ej = Em_InP(v)
    E = []
    if type(x) is np.ndarray:
        for position in x:
            if position <= d_mul_cm:
                Ex = Ej - (e * N_mul * position) / eps_InP
            elif d_mul_cm < position <= d_mul_cm + d_charge_cm:
                Ex = Ej - (e * N_mul * d_mul_cm) / eps_InP - (e * N_ch * (position - d_mul_cm)) / eps_InP
            else:
                raise BaseException("Wrong x: %.2e" % position)

            if Ex <= 0:
                E.append(0)
            else:
                E.append(Ex)
        return np.asarray(E)
    else:
        if x <= d_mul_cm:
            Ex = Ej - (e * N_mul * x) / eps_InP
        elif d_mul_cm < x <= d_mul_cm + d_charge_cm:
            Ex = Ej - (e * N_mul * d_mul_cm) / eps_InP - (e * N_ch * (x - d_mul_cm)) / eps_InP
        else:
            raise BaseException("Wrong x: %.2e" % x)

        if Ex <= 0:
            return 0
        else:
            return Ex


def ElectricField_InGaAs(w_cm, x, N_list, d_list):
    N_abs, = N_list
    d_abs, = d_list
    d_InGaAs = [dm_InGaAs(Em_InGaAs(element), 7.5e14, 3e-4) for element in Voltage_InGaAs]
    v = utils.find(d_InGaAs, Voltage_InGaAs, w_cm, 'linear')
    Ej = Em_InGaAs(v)
    E = []
    if type(x) is np.ndarray:
        for position in x:
            if position <= d_abs:
                Ex = Ej - (e * N_abs * position) / eps_InGaAs
            else:
                raise BaseException("Wrong x: %.2e" % position)

            if Ex <= 0:
                E.append(0)
            else:
                E.append(Ex)
        return np.asarray(E)
    else:
        if x <= d_abs:
            Ex = Ej - (e * N_abs * x) / eps_InGaAs
        else:
            raise BaseException("Wrong x: %.2e" % x)

        if Ex <= 0:
            return 0
        else:
            return Ex


def intergral_GammaInP(x_cm, mt_kg, T, N_list, d_list):
    w_cm = x_cm[-1]
    F_Gamma = np.sqrt(24 * mt_kg * (kB * T) ** 3) / (e * hbar) / 100  # [V/cm]
    ratio = ElectricField_InP(w_cm, x_cm, N_list, d_list) / F_Gamma
    y = 2 * np.sqrt(3 * np.pi) * ratio * np.exp(ratio ** 2)
    return utils.ydx(x_cm, y, 0, len(x_cm) - 1)


def intergral_GammaInGaAs(x_cm, mt_kg, T, N_list, d_list):
    w_cm = x_cm[-1]
    F_Gamma = np.sqrt(24 * mt_kg * (kB * T) ** 3) / (e * hbar) / 100  # [V/cm]
    ratio = ElectricField_InGaAs(w_cm, x_cm, N_list, d_list) / F_Gamma
    y = 2 * np.sqrt(3 * np.pi) * ratio * np.exp(ratio ** 2)
    return utils.ydx(x_cm, y, 0, len(x_cm) - 1)