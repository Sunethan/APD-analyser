import numpy as np
import os
import csv
import physics as phys
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import matplotlib.pylab as pylab
import DataAnalysis as Data
import utils
import GenerationRate.BandToBandTunneling as BTB
from scipy.optimize import curve_fit
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (20, 9.3),
          'axes.labelsize': 'x-large',
          'axes.titlesize':'x-large',
          'xtick.labelsize':'x-large',
          'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)
plt.rcParams.update({'font.size': 9})

# 物理常數
kB = 1.38e-23  # [J/k]
me = 9.11e-31  # [kg]
e = 1.6e-19  # [C]
eps_InP = 12.5 * 8.85e-14  # [F/cm]
eps_InGaAs = 13.9 * 8.85e-14  # [F/cm] In 0.53 Ga 0.47 As
eps_InGaAsP = 13.436 * 8.85e-14  # [F/cm] Approximated by In 0.53 Ga 0.47 As 0.65 P 0.35
h_bar = 1.054e-34  # [J-s]
Lifetime_p = {'InP': 0.5, 'InGaAs': 10, 'InGaAsP': 1000}  # [ns]  300, 8
Lifetime_n = {'InP': 0.1, 'InGaAs': 500, 'InGaAsP': 10}  # [ns]  2.89, 0.25
Lifetime_ref = {'InP': ['Parks(1996)', {'n': 1, 'p': 1}], 'InGaAs': ['Datta(1999)', {'n': 0.25, 'p': 8}]}
Eti = {'InP': -0.025, 'InGaAs': 0.16}

# 元件參數
A = np.pi * (120e-4 ** 2)  # [cm2]

# 溫度範圍
Temperature = np.arange(240, 340, 10)  # [K]
T_min = 240
T_max = 330
T_analysis = np.arange(T_min, T_max + 10, 10)
T_analysis_IT = np.arange(T_min + 20, T_max + 10, 10)
invT_list = np.arange(1 / T_max, 1 / T_min + 1e-5, 1e-5)

# 繪圖參數
count = 6
ColorSet10 = ['orangered', 'yellowgreen', 'goldenrod', 'darkviolet', 'darkorange',
              'brown', 'b', 'r', 'fuchsia', 'g']
LineSet2 = ['-', '-.']
ColorModel = {'SRH': 'r', 'TAT': 'b'}


class CurrentFitting(object):
    def __init__(self, voltage_settings, mode, electric_field, doping, grading, effective_mass, structure, others):
        self.v_min, self.v_max, v_max_range, self.Vpt, self.V1, self.V2 = voltage_settings
        self.method, self.mechanism, self.material = mode
        location_electric_field, label_electric_field = electric_field
        location_doping, label_doping = doping
        self.d_mul, self.interface_um = structure  # interface_um = [-3.62, -3.5, -0.5]
        self.effective_mass_InP = effective_mass['InP']
        self.effective_mass_InGaAs = effective_mass['InGaAs']
        self.TCAD_IV, self.TCAD_lifetime, self.TCAD_check = others

        # 設定電壓範圍
        v_step = 0.1
        iterations = (self.v_max['InGaAs'] - self.v_min['InP']) / v_step
        self.voltage = np.asarray([round(-self.v_min['InP'] - v_step * i, 1) for i in range(int(iterations))])
        self.V_InP = np.asarray([element for element in self.voltage
                                 if abs(self.v_min['InP']) <= abs(element) <= self.v_max['InP']])
        self.V_InGaAs = np.asarray([element for element in self.voltage
                                    if abs(self.v_min['InGaAs']) <= abs(element) <= self.v_max['InGaAs']])
        if v_max_range == 'All':
            self.Temperature_v_max = {240: 32.8, 250: 34, 260: 34.6, 270: 35.8, 280: 37.0,
                                      290: 37.6, 300: 38.8, 310: 39.4, 320: 40, 330: 41.2}
            for T in T_analysis:
                self.Temperature_v_max[T] = self.Temperature_v_max[T] - 0.3
        elif v_max_range == 'Partial':
            self.Temperature_v_max = {T: self.v_max['InGaAs'] for T in T_analysis}  #
        else:
            raise BaseException("Wrong InGaAs analysis range: %s" % v_max_range)

        # 製作 guess & bound
        SRH_InP_guess_IV = {T: [Eti['InP'], 1, 1] for T in T_analysis}
        SRH_InP_bound_IV = {T: ([-phys.Eg_InP(300) / 2, 1, 1], [phys.Eg_InP(300) / 2, 10, 10]) for T in T_analysis}
        SRH_InGaAs_guess_IV = {T: [Eti['InGaAs'], 1, 1] for T in T_analysis}
        SRH_InGaAs_bound_IV = {T: ([0, 0.1, 0.1], [phys.Eg_InGaAs(T) / 2, 10, 10]) for T in T_analysis}

        TAT_InP_guess_IV = {T: [Eti['InP'], 1, 1] for T in T_analysis}
        TAT_InP_bound_IV = {T: ([-0.1, 1, 1], [0.1, 1.5, 1.5]) for T in T_analysis}
        TAT_InGaAs_guess_IV = {T: [Eti['InGaAs'], 1, 1] for T in T_analysis}
        TAT_InGaAs_bound_IV = {T: ([0, 0.5, 0.85], [phys.Eg_InGaAs(T) / 2, 1.5, 1.5]) for T in T_analysis}

        # 製作 guess & bounds for IT fitting  (Eti, tp, tn, alpha_p, alpha_n)
        SRH_InP_guess_IT = {V: [Eti['InP'], 1, 1, 10, 1] for V in self.V_InP}
        SRH_InP_bound_IT = {V: ([-phys.Eg_InP(300) / 2, 1, 1, 0.1, 0.1], [phys.Eg_InP(300) / 2, 3, 3, 10, 10]) for V in self.V_InP}
        SRH_InGaAs_guess_IT = {V: [Eti['InGaAs'], 1, 1, 5, 5] for V in self.V_InGaAs}
        SRH_InGaAs_bound_IT = {V: ([-phys.Eg_InGaAs(300) / 2, 1e-1, 1, 0, 0],
                                   [phys.Eg_InGaAs(300) / 2, 1, 10, 8, 8]) for V in self.V_InGaAs}

        TAT_InP_guess_IT = {V: [Eti['InP'], 1, 1, 4, 4] for V in self.V_InP}
        TAT_InP_bound_IT = {V: ([- phys.Eg_InP(300) / 2, 0.8, 0.8, 1, 1], [0, 1.5, 1.5, 8, 8]) for V in self.V_InP}
        TAT_InGaAs_guess_IT = {V: [Eti['InGaAs'], 1, 1, 5, 5] for V in self.V_InGaAs}
        TAT_InGaAs_bound_IT = {V: ([-phys.Eg_InGaAs(300) / 2, 1e-1, 1, 0, 0],
                                   [phys.Eg_InGaAs(300) / 2, 1, 10, 8, 8]) for V in self.V_InGaAs}

        self.guess = {'InP': {'SRH': {'IV': SRH_InP_guess_IV, 'IT': SRH_InP_guess_IT},
                              'TAT': {'IV': TAT_InP_guess_IV, 'IT': TAT_InP_guess_IT}},
                      'InGaAs': {'SRH': {'IV': SRH_InGaAs_guess_IV, 'IT': SRH_InGaAs_guess_IT},
                                 'TAT': {'IV': TAT_InGaAs_guess_IV, 'IT': TAT_InGaAs_guess_IT}}}
        self.bound = {'InP': {'SRH': {'IV': SRH_InP_bound_IV, 'IT': SRH_InP_bound_IT},
                              'TAT': {'IV': TAT_InP_bound_IV, 'IT': TAT_InP_bound_IT}},
                      'InGaAs': {'SRH': {'IV': SRH_InGaAs_bound_IV, 'IT': SRH_InGaAs_bound_IT},
                                 'TAT': {'IV': TAT_InGaAs_bound_IV, 'IT': TAT_InGaAs_bound_IT}}}

        # 讀取IV，這裡必須給出 RawIV，不論TCAD還是實驗。
        self.RawIV = dict()
        LocationRawIV = {T: '/Users/Ethan/PycharmProjects/Python/Avalanche Photodiode/Report/1105/data' \
                            '/raw data/T' + str(T) + '/T.csv' for T in T_analysis}
        fieldnames = ['DataName', 'V1', 'I1', 'I2', 'I11', 'I4']
        for T in T_analysis:
            tempVoltage = []
            tempCurrent = []
            with open(LocationRawIV[T], newline='', encoding='utf-8-sig') as csv_file:
                rows = csv.DictReader(csv_file, fieldnames=fieldnames)
                for _index, row in enumerate(rows):
                    if row['DataName'] == 'DataValue':
                        tempVoltage.append(float(row['V1']))
                        tempCurrent.append(float(row['I11']))
                self.RawIV[T] = Data.CSV(data=None, xlabel=None, ylabel=None, x=tempVoltage, y=tempCurrent)

        # 讀取 InP & InGaAs 最大電場與偏壓的分佈
        self.Ef_InP = Data.CSV(location_electric_field['InP'],
                               label_electric_field['InP'][grading], label_electric_field['InP'][grading])
        self.Ef_InGaAs = Data.CSV(location_electric_field['InGaAs'],
                                  label_electric_field['InGaAs'][grading], label_electric_field['InGaAs'][grading])
        self.DopingProfile = Data.DopingProfile(location_doping, label_doping[grading], label_doping[grading])

        #
        self.material_voltage = {'InP': self.V_InP, 'InGaAs': self.V_InGaAs}
        self.weight = {'InP': 1 / abs(self.V_InP), 'InGaAs': 1 / abs(self.V_InGaAs)}
        self.result = dict()
        for item in self.method:
            if item == 'IV':
                self.result['IV'] = {item: {model: {T: self.FitIV(T, item, model, self.guess[item][model]['IV'][T],
                                                                   self.bound[item][model]['IV'][T], fitsigma=1.5)
                                                    for T in T_analysis} for model in self.mechanism}
                                     for item in self.material}
                self.Lifetime = {item: {model: {T: self.result['IV'][item][model][T][2] for T in T_analysis}
                                        for model in self.mechanism} for item in self.material}
                self.Lifetime['InGaAsP'] = {model: {T: [Lifetime_p['InGaAsP'], Lifetime_n['InGaAsP']]
                                                    for T in T_analysis} for model in self.mechanism}
            if item == 'IT':
                self.result['IT'] = {item: {model: {V: self.FitIT(V, item, model, self.guess[item][model]['IT'][V],
                                                                   self.bound[item][model]['IT'][V], fitsigma=1)
                                                    for V in self.material_voltage[item]} for model in self.mechanism}
                                     for item in self.material}
        '''
        self.BTB = {item: {T: self.PlotIV(T, item, 'BTB', ['All', self.effective_mass_InP]) for T in T_analysis} for
                    item in self.material}
        '''

    def read_data(self, temperature):
        return self.RawIV[temperature]

    def read_result(self):
        return self.result

    def dm_InP(self, E_Vcm, ND, ND_c, d_mul, d_charge):
        d = E_Vcm * eps_InP / (e * ND)  # [cm]
        if type(d) is np.ndarray:
            dm_list = []
            for i, x in enumerate(d):
                if x <= d_mul:
                    dm_list.append(x)
                else:
                    E2 = E_Vcm[i] - (e * ND * d_mul) / eps_InP
                    d2 = E2 * eps_InP / (e * ND_c)
                    if d2 <= d_charge:
                        dm_list.append(d_mul + d2)
                    else:
                        dm_list.append(d_mul + d_charge)
            return np.asarray(dm_list)  # [cm]
        else:
            if d <= d_mul:
                return d  # [cm]
            else:
                E2 = E_Vcm - (e * ND * d_mul) / eps_InP
                d2 = E2 * eps_InP / (e * ND_c)
                if d2 <= d_charge:
                    return d_mul + d2  # [cm]
                else:
                    return d_mul + d_charge  # [cm]

    def dm_InGaAs(self, E, ND_abs, d_abs):
        d = E * eps_InGaAs / (e * ND_abs)
        if type(d) is np.ndarray:
            dm_list = []
            for x in d:
                if x <= d_abs:
                    dm_list.append(x)
                else:
                    dm_list.append(d_abs)
            return np.asarray(dm_list)
        else:
            if d <= d_abs:
                return d
            else:
                return d_abs

    def Em_InP(self, V):
        return utils.find(self.Ef_InP.X, self.Ef_InP.Y, -abs(V), 'linear')

    def Em_InGaAs(self, V):
        return utils.find(self.Ef_InGaAs.X, self.Ef_InGaAs.Y, -abs(V), 'linear')

    def FitIV(self, T, material, type, guess, bound, fitsigma):
        """
        :param T:
        :param material:
        :return: V, I, popt
        """
        if material == 'InP':
            V_InP = np.asarray([V for V in self.RawIV[T].X if -self.v_min['InP'] >= V > -self.v_max['InP']])
            F_InP = np.asarray([self.Em_InP(V) for V in V_InP])
            I_InP = np.asarray([I for i, I in enumerate(self.RawIV[T].Y) if self.RawIV[T].X[i] in V_InP])
            def lifetime(tp, tn):
                alpha = 1.5
                tau_p0 = Lifetime_p['InP'] * 1e-9  # [s]
                tau_n0 = Lifetime_n['InP'] * 1e-9  # [s]
                tau_p = tp * tau_p0 * (T / 300) ** alpha
                tau_n = tn * tau_n0 * (T / 300) ** alpha
                return tau_p, tau_n
            if type == 'TAT':
                def TAT_InP_IV(X, Eti, tp, tn):
                    Emax_Vcm, T = X
                    alpha = 1.5
                    # tp = 1
                    # tn = 0.1
                    mt = self.effective_mass_InP
                    prefactor = 1

                    me = 9.11e-31
                    Nc300 = 5.716e17  # [cm-3]
                    Nv300 = 1.143e19  # [cm-3]
                    tau_p0 = Lifetime_p['InP'] * 1e-9  # [s]
                    tau_n0 = Lifetime_n['InP'] * 1e-9  # [s]
                    tau_p = tp * tau_p0 * (T / 300) ** alpha
                    tau_n = tn * tau_n0 * (T / 300) ** alpha
                    ND = 5e16  # [cm-3]
                    Ncharge = 7.8e16  # [cm-3]
                    d_mul = self.d_mul  # 0.42e-4  # [cm]
                    d_ch = 0.2e-4  # [cm]

                    ni = np.sqrt(Nc300 * Nv300) * (T / 300) ** 1.5 * np.exp(-e * phys.Eg_InP(T) / (2 * kB * T))
                    G_SRH = ni / (2 * np.sqrt(tau_p * tau_n) * np.cosh(e * Eti / (kB * T) + 0.5 * np.log(tau_p / tau_n)))

                    dM = self.dm_InP(Emax_Vcm, ND, Ncharge, d_mul, d_ch)  # 0.42e-4  # [cm]
                    F_Gamma = np.sqrt(24 * (mt * me) * (kB * T) ** 3) / (e * h_bar) / 100  # [V/cm]
                    E1 = Emax_Vcm
                    log10_Current = []
                    for i, x in enumerate(dM):
                        if x <= d_mul:
                            E2 = E1[i] - (e * ND * x) / eps_InP
                            d_Gamma_1 = (np.sqrt(3 * np.pi) * eps_InP * F_Gamma) / (e * ND) * \
                                        (np.exp((E1[i] / F_Gamma) ** 2) - np.exp(E2 / F_Gamma ** 2))  # [cm]
                            log10_Current.append(
                                np.log10(A * e) + np.log10(prefactor * G_SRH) + np.log10(x + d_Gamma_1))
                        else:
                            E2 = E1[i] - (e * ND * d_mul) / eps_InP
                            E3 = E2 - (e * Ncharge * (x - d_mul)) / eps_InP
                            d_Gamma_1 = (np.sqrt(3 * np.pi) * eps_InP * F_Gamma) / (e * ND) * \
                                        (np.exp((E1[i] / F_Gamma) ** 2) - np.exp(E2 / F_Gamma ** 2))  # [cm]
                            d_Gamma_2 = (np.sqrt(3 * np.pi) * eps_InP * F_Gamma) / (e * Ncharge) * \
                                        (np.exp((E2 / F_Gamma) ** 2) - np.exp(E3 / F_Gamma ** 2))  # [cm]
                            log10_Current.append(
                                np.log10(A * e) + np.log10(prefactor * G_SRH) + np.log10(
                                    x + d_Gamma_1 + d_Gamma_2))
                    return np.asarray(log10_Current)

                TAT_InP_popt, TAT_InP_pcov = curve_fit(TAT_InP_IV, (F_InP, T), np.log10(I_InP), p0=guess, bounds=bound,
                                                       sigma=abs(np.log10(I_InP)) ** fitsigma)
                print('[TAT]   InP (%sK)  Eti:  %.3f,  tp: %.3e,  tn: %.3e' %
                      (T, TAT_InP_popt[0], TAT_InP_popt[1], TAT_InP_popt[2]))
                Eti = TAT_InP_popt[0]
                mt = self.effective_mass_InP
                tau_p, tau_n = lifetime(TAT_InP_popt[1], TAT_InP_popt[2])
                return V_InP, 10 ** TAT_InP_IV((F_InP, T), *TAT_InP_popt), [tau_p, tau_n], Eti, mt
            elif type == 'SRH':
                def SRH_InP(X, Eti, tp, tn):
                    """
                    使用 -U ~ ni * cosh(-(Eti+ln(tp/tn))/kT) 之近似公式，而不需要使用 |Eti| >> kT 之公式。
                    內建正確的 lifetime。
                    :param X: (T, Emax_Vcm)
                    :param Eti: eV
                    :return: np.log10(I)
                    """
                    Emax_Vcm, T = X
                    alpha = 1.5  # 1
                    # tp = 1  # 0.1
                    # tn = 1  # 0.226
                    prefactor = 1

                    me = 9.11e-31
                    Nc300 = 5.716e17  # [cm-3]
                    Nv300 = 1.143e19  # [cm-3]
                    tau_p0 = Lifetime_p['InP'] * 1e-9  # [s]
                    tau_n0 = Lifetime_n['InP'] * 1e-9  # [s]
                    tau_p = tp * tau_p0 * (T / 300) ** alpha
                    tau_n = tn * tau_n0 * (T / 300) ** alpha
                    ND = 5e16  # [cm-3]
                    Ncharge = 7.8e16  # [cm-3]
                    d_mul = self.d_mul  # 0.42e-4  # [cm]
                    d_ch = 0.2e-4  # [cm]

                    ni = np.sqrt(Nc300 * Nv300) * (T / 300) ** 1.5 * np.exp(- e * phys.Eg_InP(T) / (2 * kB * T))
                    G_SRH = ni / (
                                2 * np.sqrt(tau_p * tau_n) * np.cosh(e * Eti / (kB * T) + 0.5 * np.log(tau_p / tau_n)))
                    dM = self.dm_InP(Emax_Vcm, ND, Ncharge, d_mul, d_ch)  # 0.42e-4  # [cm]
                    return np.log10(A * e) + np.log10(prefactor * G_SRH) + np.log10(dM)

                popt_SRH_InP, pcov_SRH_InP = curve_fit(SRH_InP, (F_InP, T), np.log10(I_InP), p0=guess, bounds=bound,
                                                       sigma=abs(np.log10(I_InP)) ** fitsigma)
                print('[SRH]   InP (%sK)  Eti:  %.3f,  tp: %.3e,  tn: %.3e' %
                      (T, popt_SRH_InP[0], popt_SRH_InP[1], popt_SRH_InP[2]))
                Eti = popt_SRH_InP[0]
                mt = self.effective_mass_InP
                tau_p, tau_n = lifetime(popt_SRH_InP[1], popt_SRH_InP[2])
                return V_InP, 10 ** SRH_InP((F_InP, T), *popt_SRH_InP), [tau_p, tau_n], Eti, mt
            else:
                raise BaseException("Wrong type: %s" % type)
        elif material == 'InGaAs':
            V_InGaAs = np.asarray([V for V in self.RawIV[T].X
                                   if -self.Temperature_v_max[T] <= V <= -self.v_min['InGaAs']])
            F_InGaAs = np.asarray([self.Em_InGaAs(V) for V in V_InGaAs])
            I_InGaAs = np.asarray([I for i, I in enumerate(self.RawIV[T].Y) if self.RawIV[T].X[i] in V_InGaAs])
            def lifetime(tp, tn):
                alpha = 1.5
                tau_p0 = Lifetime_p['InGaAs'] * 1e-9  # [s]
                tau_n0 = Lifetime_n['InGaAs'] * 1e-9  # [s]
                tau_p = tp * tau_p0 * (T / 300) ** alpha
                tau_n = tn * tau_n0 * (T / 300) ** alpha
                return tau_p, tau_n
            if type == 'TAT':
                def TAT_InGaAs_IV(X, Eti, tp, tn):
                    Emax_Vcm, T = X
                    prefactor = 1
                    # tp = 1
                    # tn = 1
                    mt = self.effective_mass_InGaAs
                    alpha = 1.5
                    me = 9.11e-31
                    Nc300 = 2.53956e17  # [cm-3]
                    Nv300 = 7.51e18  # [cm-3]
                    tau_p0 = Lifetime_p['InGaAs'] * 1e-9  # [s]
                    tau_n0 = Lifetime_n['InGaAs'] * 1e-9  # [s]
                    tau_p = tp * tau_p0 * (T / 300) ** alpha
                    tau_n = tn * tau_n0 * (T / 300) ** alpha
                    ND_abs = 7.5e14  # [cm-3]
                    d_InGaAs = 3e-4  # [cm]

                    ni = np.sqrt(Nc300 * Nv300) * (T / 300) ** 1.5 * np.exp(-e * phys.Eg_InGaAs(T) / (2 * kB * T))
                    G_SRH = ni / (2 * np.sqrt(tau_p * tau_n) * np.cosh(e * Eti / (kB * T) + 0.5 * np.log(tau_p / tau_n)))

                    dM = self.dm_InGaAs(Emax_Vcm, ND_abs, d_InGaAs)  # [cm]
                    F_Gamma = np.sqrt(24 * (mt * me) * (kB * T) ** 3) / (e * h_bar) / 100  # [V/cm]
                    E1 = Emax_Vcm
                    E2 = 0
                    d_Gamma = (np.sqrt(3 * np.pi) * eps_InGaAs * F_Gamma) / (e * ND_abs) * \
                              (np.exp((E1 / F_Gamma) ** 2) - np.exp((E2 / F_Gamma) ** 2))  # [cm]
                    return np.log10(A * e) + np.log10(prefactor * G_SRH) + np.log10(dM + d_Gamma)

                if len(V_InGaAs) == 0:
                    return V_InGaAs, [], [0, 0], None, None
                else:
                    TAT_InGaAs_popt, TAT_InGaAs_pcov = curve_fit(TAT_InGaAs_IV, (F_InGaAs, T), np.log10(I_InGaAs),
                                                                 p0=guess,
                                                                 bounds=bound,
                                                                 sigma=abs(np.log10(I_InGaAs)) ** fitsigma)
                    print('[TAT]   InGaAs (%sK)  Eti:  %.3f,  tp: %.3e,  tn: %.3e' %
                          (T, TAT_InGaAs_popt[0], TAT_InGaAs_popt[1], TAT_InGaAs_popt[2]))

                    Eti = TAT_InGaAs_popt[0]
                    mt = self.effective_mass_InGaAs
                    tau_p, tau_n = lifetime(TAT_InGaAs_popt[1], TAT_InGaAs_popt[2])

                    return V_InGaAs, 10 ** TAT_InGaAs_IV((F_InGaAs, T), *TAT_InGaAs_popt), [tau_p, tau_n], Eti, mt
            elif type == 'SRH':
                def SRH_InGaAs_IV(X, Eti, tp, tn):
                    Emax_Vcm, T = X
                    prefactor = 1
                    # tp = 1
                    # tn = 1
                    alpha = 1.5
                    me = 9.11e-31
                    Nc300 = 2.53956e17  # [cm-3]
                    Nv300 = 7.51e18  # [cm-3]
                    tau_p0 = Lifetime_p['InGaAs'] * 1e-9  # [s]
                    tau_n0 = Lifetime_n['InGaAs'] * 1e-9  # [s]
                    tau_p = tp * tau_p0 * (T / 300) ** alpha
                    tau_n = tn * tau_n0 * (T / 300) ** alpha
                    ND_abs = 7.5e14  # [cm-3]
                    d_InGaAs = 3e-4  # [cm]

                    ni = np.sqrt(Nc300 * Nv300) * (T / 300) ** 1.5 * np.exp(-e * phys.Eg_InGaAs(T) / (2 * kB * T))
                    G_SRH = ni / (2 * np.sqrt(tau_p * tau_n) * np.cosh(e * Eti / (kB * T) + 0.5 * np.log(tau_p / tau_n)))

                    dM = self.dm_InGaAs(Emax_Vcm, ND_abs, d_InGaAs)  # [cm]
                    return np.log10(A * e) + np.log10(prefactor * G_SRH) + np.log10(dM)

                if len(V_InGaAs) == 0:
                    return V_InGaAs, [], [0, 0], None
                else:
                    SRH_InGaAs_popt, SRH_InGaAs_pcov = curve_fit(SRH_InGaAs_IV, (F_InGaAs, T), np.log10(I_InGaAs),
                                                                 p0=guess, bounds=bound,
                                                                 sigma=abs(np.log10(I_InGaAs)) ** fitsigma)
                    print('[SRH]   InGaAs (%sK)  Eti:  %.3f,   tp: %.3e,  tn: %.3e' %
                          (T, SRH_InGaAs_popt[0], SRH_InGaAs_popt[1], SRH_InGaAs_popt[2]))
                    Eti = SRH_InGaAs_popt[0]
                    tau_p, tau_n = lifetime(SRH_InGaAs_popt[1], SRH_InGaAs_popt[2])
                    return V_InGaAs, 10 ** SRH_InGaAs_IV((F_InGaAs, T), *SRH_InGaAs_popt), [tau_p, tau_n], Eti
            else:
                raise BaseException("Wrong type: %s" % type)
        else:
            raise BaseException("Wrong material: %s" % material)

    def FitIT(self, V, material, type, guess, bound, fitsigma):
        if material == 'InP':
            I_InP = np.asarray([utils.find(self.RawIV[T].X, self.RawIV[T].Y, V, 'log') for T in T_analysis_IT])
            if type == 'TAT':
                def TAT_InP_IT(X, Eti, tp, tn, alpha_p, alpha_n):
                    T, Emax_Vcm = X
                    mt = self.effective_mass_InP
                    prefactor = 1

                    me = 9.11e-31
                    Nc300 = 5.716e17  # [cm-3]
                    Nv300 = 1.143e19  # [cm-3]
                    tau_p0 = Lifetime_p['InP'] * 1e-9  # [s]
                    tau_n0 = Lifetime_n['InP'] * 1e-9  # [s]
                    tau_p = tp * tau_p0 * (T / 300) ** alpha_p
                    tau_n = tn * tau_n0 * (T / 300) ** alpha_n
                    ND = 5e16  # [cm-3]
                    Ncharge = 7.8e16  # [cm-3]
                    d_mul = self.d_mul  # 0.42e-4  # [cm]
                    d_ch = 0.2e-4  # [cm]

                    ni = np.sqrt(Nc300 * Nv300) * (T / 300) ** 1.5 * np.exp(-e * phys.Eg_InP(T) / (2 * kB * T))
                    G_SRH = ni / (2 * np.sqrt(tau_p * tau_n) * np.cosh(e * Eti / (kB * T) + 0.5 * np.log(tau_p / tau_n)))

                    dM = self.dm_InP(Emax_Vcm, ND, Ncharge, d_mul, d_ch)  # 0.42e-4  # [cm]
                    F_Gamma = np.sqrt(24 * (mt * me) * (kB * T) ** 3) / (e * h_bar) / 100  # [V/cm]
                    E1 = Emax_Vcm
                    if dM <= d_mul:
                        E2 = E1 - (e * ND * dM) / eps_InP
                        d_Gamma_1 = (np.sqrt(3 * np.pi) * eps_InP * F_Gamma) / (e * ND) * \
                                    (np.exp((E1 / F_Gamma) ** 2) - np.exp(E2 / F_Gamma ** 2))  # [cm]
                        return np.log10(A * e) + np.log10(prefactor * G_SRH) + np.log10(dM + d_Gamma_1)
                    else:
                        E2 = E1 - (e * ND * d_mul) / eps_InP
                        E3 = E2 - (e * Ncharge * (dM - d_mul)) / eps_InP
                        d_Gamma_1 = (np.sqrt(3 * np.pi) * eps_InP * F_Gamma) / (e * ND) * \
                                    (np.exp((E1 / F_Gamma) ** 2) - np.exp(E2 / F_Gamma ** 2))  # [cm]
                        d_Gamma_2 = (np.sqrt(3 * np.pi) * eps_InP * F_Gamma) / (e * Ncharge) * \
                                    (np.exp((E2 / F_Gamma) ** 2) - np.exp(E3 / F_Gamma ** 2))  # [cm]
                        return np.log10(A * e) + np.log10(prefactor * G_SRH) + np.log10(dM + d_Gamma_1 + d_Gamma_2)

                popt, pcov = curve_fit(TAT_InP_IT, (T_analysis_IT, self.Em_InP(V)), np.log10(I_InP), p0=guess,
                                       bounds=bound, sigma=abs(np.log10(I_InP)) ** fitsigma)
                Eti, tp, tn, alpha_p, alpha_n = popt
                print('[TAT]   InP (%.1f)  Eti:  %.3f,  tp: %.3e,  tn: %.3e,  alpha(p): %.3e,  alpha(n): %.3e' %
                      (V, Eti, tp, tn, alpha_p, alpha_n))
                return T_analysis_IT, 10 ** TAT_InP_IT((T_analysis_IT, self.Em_InP(V)), *popt), \
                       Eti, [tp, tn, alpha_p, alpha_n]
            elif type == 'SRH':
                def SRH_InP_IT(X, Eti, tp, tn, alpha_n, alpha_p):
                    T, Emax_Vcm = X
                    # tp = 1
                    # tn = 1
                    prefactor = 1

                    Nc300 = 5.716e17  # [cm-3]
                    Nv300 = 1.143e19  # [cm-3]
                    tau_p0 = Lifetime_p['InP'] * 1e-9  # [s]
                    tau_n0 = Lifetime_n['InP'] * 1e-9  # [s]
                    tau_p = tp * tau_p0 * (T / 300) ** alpha_p
                    tau_n = tn * tau_n0 * (T / 300) ** alpha_n
                    ND = 5e16  # [cm-3]
                    Ncharge = 7.8e16  # [cm-3]
                    d_mul = self.d_mul  # 0.42e-4  # [cm]
                    d_ch = 0.2e-4  # [cm]

                    ni = np.sqrt(Nc300 * Nv300) * (T / 300) ** 1.5 * np.exp(-e * phys.Eg_InP(T) / (2 * kB * T))
                    G_SRH = ni / (2 * np.sqrt(tau_p * tau_n) * np.cosh(e * Eti / (kB * T) + 0.5 * np.log(tau_p / tau_n)))

                    dM = self.dm_InP(Emax_Vcm, ND, Ncharge, d_mul, d_ch)  # 0.42e-4  # [cm]
                    return np.log10(A * e) + np.log10(prefactor * G_SRH) + np.log10(dM)

                popt, pcov = curve_fit(SRH_InP_IT, (T_analysis_IT, self.Em_InP(V)), np.log10(I_InP), p0=guess,
                                       bounds=bound, sigma=abs(np.log10(I_InP)) ** fitsigma)
                Eti, tp, tn, alpha_p, alpha_n = popt
                print('[SRH]   InP (%.1f)  Eti:  %.3f,  tp: %.3e,  tn: %.3e,  alpha(p): %.3e,  alpha(n): %.3e' %
                      (V, Eti, tp, tn, alpha_p, alpha_n))
                return T_analysis_IT, 10 ** SRH_InP_IT((T_analysis_IT, self.Em_InP(V)), *popt), \
                       Eti, [tp, tn, alpha_p, alpha_n]
            else:
                raise BaseException("Wrong type: %s" % type)
        elif material == 'InGaAs':
            I_InGaAs = np.asarray([utils.find(self.RawIV[T].X, self.RawIV[T].Y, V, 'log') for T in T_analysis_IT])

            # 檢查電流是否隨著溫度遞增
            if abs(V) > abs(self.Temperature_v_max[T_analysis_IT[0]]):
                raise BaseException("Voltage is too large: %s > Vmax(InGaAs,240K) = %s" %
                                    (abs(V), abs(self.Temperature_v_max[T_analysis_IT[0]])))

            if type == 'TAT':
                def TAT_InGaAs_IT(X, Eti, tp, tn, alpha_p, alpha_n):
                    T, Emax_Vcm = X
                    prefactor = 1
                    mt = self.effective_mass_InGaAs
                    me = 9.11e-31
                    Nc300 = 2.53956e17  # [cm-3]
                    Nv300 = 7.51e18  # [cm-3]
                    tau_p0 = Lifetime_p['InGaAs'] * 1e-9  # [s]
                    tau_n0 = Lifetime_n['InGaAs'] * 1e-9  # [s]
                    tau_p = tp * tau_p0 * (T / 300) ** alpha_p
                    tau_n = tn * tau_n0 * (T / 300) ** alpha_n
                    ND_abs = 7.5e14  # [cm-3]
                    d_InGaAs = 3e-4  # [cm]

                    ni = np.sqrt(Nc300 * Nv300) * (T / 300) ** 1.5 * np.exp(- e * phys.Eg_InGaAs(T) / (2 * kB * T))
                    G_SRH = ni / (2 * np.sqrt(tau_p * tau_n) * np.cosh(e * Eti / (kB * T) + 0.5 * np.log(tau_p / tau_n)))

                    dM = self.dm_InGaAs(Emax_Vcm, ND_abs, d_InGaAs)  # [cm]
                    F_Gamma = np.sqrt(24 * (mt * me) * (kB * T) ** 3) / (e * h_bar) / 100  # [V/cm]
                    E1 = Emax_Vcm
                    E2 = 0
                    d_Gamma = (np.sqrt(3 * np.pi) * eps_InGaAs * F_Gamma) / (e * ND_abs) * \
                              (np.exp((E1 / F_Gamma) ** 2) - np.exp((E2 / F_Gamma) ** 2))  # [cm]
                    return np.log10(A * e) + np.log10(prefactor * G_SRH) + np.log10(dM + d_Gamma)

                popt, pcov = curve_fit(TAT_InGaAs_IT, (T_analysis_IT, self.Em_InGaAs(V)), np.log10(I_InGaAs),
                                       p0=guess, bounds=bound, sigma=abs(np.log10(I_InGaAs)) ** fitsigma)
                Eti, tp, tn, alpha_p, alpha_n = popt
                print('[TAT]   InGaAs (%.1f)  Eti:  %.3f,  tp: %.3e,  tn: %.3e,  alpha(p): %.3e,  alpha(n): %.3e' %
                      (V, Eti, tp, tn, alpha_p, alpha_n))
                return T_analysis_IT, 10 ** TAT_InGaAs_IT((T_analysis_IT, self.Em_InGaAs(V)), *popt), \
                       Eti, [tp, tn, alpha_p, alpha_n]
            elif type == 'SRH':
                def SRH_InGaAs_IT(X, Eti, tp, tn, alpha_p, alpha_n):
                    T, Emax_Vcm = X
                    prefactor = 1
                    Nc300 = 2.53956e17  # [cm-3]
                    Nv300 = 7.51e18  # [cm-3]
                    tau_p0 = Lifetime_p['InGaAs'] * 1e-9  # [s]
                    tau_n0 = Lifetime_n['InGaAs'] * 1e-9  # [s]
                    tau_p = tp * tau_p0 * (T / 300) ** alpha_p
                    tau_n = tn * tau_n0 * (T / 300) ** alpha_n
                    ND_abs = 7.5e14  # [cm-3]
                    d_InGaAs = 3e-4  # [cm]

                    ni = np.sqrt(Nc300 * Nv300) * (T / 300) ** 1.5 * np.exp(- e * phys.Eg_InGaAs(T) / (2 * kB * T))
                    G_SRH = ni / (2 * np.sqrt(tau_p * tau_n) * np.cosh(e * Eti / (kB * T) + 0.5 * np.log(tau_p / tau_n)))

                    dM = self.dm_InGaAs(Emax_Vcm, ND_abs, d_InGaAs)  # [cm]
                    return np.log10(A * e) + np.log10(prefactor * G_SRH) + np.log10(dM)

                popt, pcov = curve_fit(SRH_InGaAs_IT, (T_analysis_IT, self.Em_InGaAs(V)), np.log10(I_InGaAs),
                                       p0=guess, bounds=bound, sigma=abs(np.log10(I_InGaAs)) ** fitsigma)
                Eti, tp, tn, alpha_p, alpha_n = popt
                print('[SRH]   InGaAs (%.1f)  Eti:  %.3f,  tp: %.3e,  tn: %.3e,  alpha(p): %.3e,  alpha(n): %.3e' %
                      (V, Eti, tp, tn, alpha_p, alpha_n))
                return T_analysis_IT, 10 ** SRH_InGaAs_IT((T_analysis_IT, self.Em_InGaAs(V)), *popt), \
                       Eti, [tp, tn, alpha_p, alpha_n]
            else:
                raise BaseException("Wrong type: %s" % type)
        else:
            raise BaseException("Wrong material: %s" % material)

    def d_Gamma(self, T, material):
        if material == 'InP':
            V_InP = np.asarray([V for V in self.RawIV[T].X if -self.v_min['InP'] >= V > -self.v_max['InP']])
            F_InP = np.asarray([self.Em_InP(V) for V in V_InP])
            def Gamma_InP(X):
                Emax_Vcm, T = X
                mt = self.effective_mass_InP

                me = 9.11e-31
                ND = 5e16  # [cm-3]
                Ncharge = 7.8e16  # [cm-3]
                d_mul = 0.42e-4  # [cm]
                d_ch = 0.2e-4  # [cm]

                dM = self.dm_InP(Emax_Vcm, ND, Ncharge, d_mul, d_ch)  # 0.42e-4  # [cm]
                F_Gamma = np.sqrt(24 * (mt * me) * (kB * T) ** 3) / (e * h_bar) / 100  # [V/cm]
                E1 = Emax_Vcm
                d_Gamma = []
                for i, x in enumerate(dM):
                    if x <= d_mul:
                        E2 = E1[i] - (e * ND * x) / eps_InP
                        d_Gamma_1 = (np.sqrt(3 * np.pi) * eps_InP * F_Gamma) / (e * ND) * \
                                    (np.exp((E1[i] / F_Gamma) ** 2) - np.exp(E2 / F_Gamma ** 2))  # [cm]
                        d_Gamma.append(d_Gamma_1)
                    else:
                        E2 = E1[i] - (e * ND * d_mul) / eps_InP
                        E3 = E2 - (e * Ncharge * (x - d_mul)) / eps_InP
                        d_Gamma_1 = (np.sqrt(3 * np.pi) * eps_InP * F_Gamma) / (e * ND) * \
                                    (np.exp((E1[i] / F_Gamma) ** 2) - np.exp(E2 / F_Gamma ** 2))  # [cm]
                        d_Gamma_2 = (np.sqrt(3 * np.pi) * eps_InP * F_Gamma) / (e * Ncharge) * \
                                    (np.exp((E2 / F_Gamma) ** 2) - np.exp(E3 / F_Gamma ** 2))  # [cm]
                        d_Gamma.append(d_Gamma_1 + d_Gamma_2)
                return np.asarray(d_Gamma)
            return V_InP, Gamma_InP((F_InP, T))
        elif material == 'InGaAs':
            V_InGaAs = np.asarray([V for V in self.RawIV[T].X
                                   if -self.Temperature_v_max[T] <= V <= -self.v_min['InGaAs']])
            F_InGaAs = np.asarray([self.Em_InGaAs(V) for V in V_InGaAs])
            def Gamma_InGaAs(X):
                Emax_Vcm, T = X
                mt = self.effective_mass_InGaAs
                me = 9.11e-31
                ND_abs = 7.5e14  # [cm-3]
                F_Gamma = np.sqrt(24 * (mt * me) * (kB * T) ** 3) / (e * h_bar) / 100  # [V/cm]
                E1 = Emax_Vcm
                E2 = 0
                d_Gamma = (np.sqrt(3 * np.pi) * eps_InGaAs * F_Gamma) / (e * ND_abs) * \
                          (np.exp((E1 / F_Gamma) ** 2) - np.exp((E2 / F_Gamma) ** 2))  # [cm]
                return d_Gamma
            return V_InGaAs, Gamma_InGaAs((F_InGaAs, T))

    def PlotIV(self, T, material, type, parameter):
        if material == 'InP':
            if type == 'BTB':
                V_range, effective_mass = parameter
                if V_range == 'All':
                    V_BTB = self.RawIV[T].X
                elif V_range == 'InP':
                    V_BTB = np.asarray([V for V in self.RawIV[T].X if -self.v_min['InP'] >= V > -self.v_max['InP']])
                elif V_range == 'InGaAs':
                    V_BTB = np.asarray([V for V in self.RawIV[T].X if -self.Temperature_v_max[T] <= V <= -self.v_min['InGaAs']])
                else:
                    raise BaseException("Wrong parameter: %s" % parameter)
                I_BTB = []
                for V in V_BTB:
                    Ej = self.Em_InP(V)
                    E_Vcm_array = self.DopingProfile.FieldTransform(Ej, self.interface_um)
                    I_BTB.append(BTB.J_BTB_InP(self.DopingProfile.X * 1e-4, E_Vcm_array, T, effective_mass) * A)
                return V_BTB, I_BTB
        elif material == 'InGaAs':
            V_InGaAs = np.asarray([V for V in self.RawIV[T].X if -self.Temperature_v_max[T] <= V <= -self.v_min['InGaAs']])
            F_InGaAs = np.asarray([self.Em_InGaAs(V) for V in V_InGaAs])
            if type == 'SRH':
                def SRH_InGaAs_IV(X, Eti, tp, tn):
                    Emax_Vcm, T = X
                    prefactor = 1
                    # tp = 1
                    # tn = 1
                    alpha = 1.5
                    me = 9.11e-31
                    Nc300 = 2.53956e17  # [cm-3]
                    Nv300 = 7.51e18  # [cm-3]
                    tau_p0 = Lifetime_p['InGaAs'] * 1e-9  # [s]
                    tau_n0 = Lifetime_n['InGaAs'] * 1e-9  # [s]
                    tau_p = tp * tau_p0 * (T / 300) ** alpha
                    tau_n = tn * tau_n0 * (T / 300) ** alpha
                    ND_abs = 7.5e14  # [cm-3]
                    d_InGaAs = 3e-4  # [cm]

                    ni = np.sqrt(Nc300 * Nv300) * (T / 300) ** 1.5 * np.exp(- e * phys.Eg_InGaAs(T) / (2 * kB * T))
                    G_SRH = ni / (
                            2 * np.sqrt(tau_p * tau_n) * np.cosh(e * Eti / (kB * T) + 0.5 * np.log(tau_p / tau_n)))

                    dM = self.dm_InGaAs(Emax_Vcm, ND_abs, d_InGaAs)  # [cm]
                    return np.log10(A * e) + np.log10(prefactor * G_SRH) + np.log10(dM)
                return V_InGaAs, 10 ** SRH_InGaAs_IV((F_InGaAs, T), *parameter)
            elif type == 'TAT':
                def TAT_InGaAs_IV(X, Eti, mt, tp, tn):
                    Emax_Vcm, T = X
                    prefactor = 1
                    # tp = 1
                    # tn = 1
                    alpha = 1.5
                    me = 9.11e-31
                    Nc300 = 2.53956e17  # [cm-3]
                    Nv300 = 7.51e18  # [cm-3]
                    tau_p0 = Lifetime_p['InGaAs'] * 1e-9  # [s]
                    tau_n0 = Lifetime_n['InGaAs'] * 1e-9  # [s]
                    tau_p = tp * tau_p0 * (T / 300) ** alpha
                    tau_n = tn * tau_n0 * (T / 300) ** alpha
                    ND_abs = 7.5e14  # [cm-3]
                    d_InGaAs = 3e-4  # [cm]

                    ni = np.sqrt(Nc300 * Nv300) * (T / 300) ** 1.5 * np.exp(- e * phys.Eg_InGaAs(T) / (2 * kB * T))
                    G_SRH = ni / (2 * np.sqrt(tau_p * tau_n) * np.cosh(e * Eti / (kB * T) + 0.5 * np.log(tau_p / tau_n)))

                    dM = self.dm_InGaAs(Emax_Vcm, ND_abs, d_InGaAs)  # [cm]
                    F_Gamma = np.sqrt(24 * (mt * me) * (kB * T) ** 3) / (e * h_bar) / 100  # [V/cm]
                    E1 = Emax_Vcm
                    E2 = 0
                    d_Gamma = (np.sqrt(3 * np.pi) * eps_InGaAs * F_Gamma) / (e * ND_abs) * (
                            np.exp((E1 / F_Gamma) ** 2) - np.exp((E2 / F_Gamma) ** 2))  # [cm]
                    return np.log10(A * e) + np.log10(prefactor * G_SRH) + np.log10(dM + d_Gamma)
                return V_InGaAs, 10 ** TAT_InGaAs_IV((F_InGaAs, T), *parameter)
            elif type == 'BTB':
                V_range, effective_mass = parameter
                if V_range == 'All':
                    V_BTB = self.RawIV[T].X
                elif V_range == 'InP':
                    V_BTB = np.asarray([V for V in self.RawIV[T].X if -self.v_min['InP'] >= V > -self.v_max['InP']])
                elif V_range == 'InGaAs':
                    V_BTB = np.asarray([V for V in self.RawIV[T].X
                                        if -self.Temperature_v_max[T] <= V <= -self.v_min['InGaAs']])
                else:
                    raise BaseException("Wrong parameter: %s" % parameter)
                I_BTB = []
                for V in V_BTB:
                    Ej = self.Em_InGaAs(V)
                    E_Vcm_array = self.DopingProfile.FieldTransform(Ej, self.interface_um)
                    I_BTB.append(BTB.J_BTB_InGaAs(self.DopingProfile.X * 1e-4, E_Vcm_array, T, effective_mass) * A)
                return V_BTB, I_BTB
            else:
                raise BaseException("Wrong type: %s" % type)
        else:
            raise BaseException("Wrong material: %s" % material)

    def PlotEm(self, T, material):
        V_net = np.asarray([V for V in self.RawIV[T].X if -self.v_min['InP'] >= V > -self.v_max['InGaAs']])
        if material == 'InP':
            F_InP = np.asarray([self.Em_InP(V) for V in V_net])
            return V_net, F_InP
        elif material == 'InGaAs':
            F_InGaAs = np.asarray([self.Em_InGaAs(V) for V in V_net])
            return V_net, F_InGaAs
        else:
            raise BaseException("Wrong material: %s" % material)

    def Plotdm(self, T, material):
        if material == 'InP':
            V_InP = np.asarray([V for V in self.RawIV[T].X if -self.v_min['InP'] >= V > -self.v_max['InP']])
            F_InP = np.asarray([self.Em_InP(V) for V in V_InP])
            return V_InP, self.dm_InP(F_InP, 5e16, 7.8e16, 0.42e-4, 0.2e-4)
        elif material == 'InGaAs':
            V_InGaAs = np.asarray([V for V in self.RawIV[T].X
                                   if -self.Temperature_v_max[T] <= V <= -self.v_min['InGaAs']])
            F_InGaAs = np.asarray([self.Em_InGaAs(V) for V in V_InGaAs])
            return V_InGaAs, self.dm_InGaAs(F_InGaAs, 7.5e14, 3e-4)
        else:
            raise BaseException("Wrong material: %s" % material)

    def plot_current_electric_field(self, number, figtitle):
        LabelMarker = 0
        temp = None
        fig_IV = plt.figure(number)
        plt.subplot(1, 2, 1)
        for i, T in enumerate(T_analysis):
            '''
            for material in self.material:
                plt.plot(self.BTB[material][T][0], self.BTB[material][T][1], linestyle=':', color=ColorSet10[i])
            '''
            for j, model in enumerate(self.mechanism):
                for material in self.material:
                    if LabelMarker < 2 and temp != model:
                        temp = model
                        plt.plot(self.result['IV'][material][model][T][0], self.result['IV'][material][model][T][1],
                                 color=ColorSet10[i], linestyle=LineSet2[j], label=model)
                        LabelMarker += 1
                    else:
                        plt.plot(self.result['IV'][material][model][T][0], self.result['IV'][material][model][T][1],
                                 color=ColorSet10[i], linestyle=LineSet2[j])
            raw_V = [element for i, element in enumerate(self.RawIV[T].X) if i % count == 0]
            raw_I = [element for i, element in enumerate(self.RawIV[T].Y) if i % count == 0]
            plt.plot(raw_V, raw_I, linestyle='none', marker='o', fillstyle='none', color=ColorSet10[i], label=T)
        plt.vlines(x=-self.Vpt, ymin=1e-13, ymax=1e-6, color='k', linestyle=':', label=r'V$_{pt}$')
        plt.grid()
        plt.yscale('log')
        plt.xlabel('Voltage (V)')
        plt.ylabel('Current (A)')
        plt.title(r'Fixed-temperature fitting')
        plt.xlim((-45, 1))
        plt.ylim((1e-13, 1e-6))
        plt.legend(loc='best')
        plt.subplot(1, 2, 2)
        plt.plot(*self.PlotEm(300, 'InP'), color='b', label='InP')
        plt.plot(*self.PlotEm(300, 'InGaAs'), color='r', label='InGaAs')
        plt.vlines(x=-self.Vpt, ymin=0, ymax=7e5, color='k',
                   linestyle=':', label=r'V$_{pt}=%s$' % -self.Vpt)
        plt.grid()
        plt.xlabel('Voltage (V)')
        plt.ylabel('Electric field (V/cm)')
        plt.title('Maximum electric field distribution')
        plt.legend(loc='best')
        plt.ylim((-0.3e5, 8e5))
        plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
        fig_IV.suptitle(figtitle, fontsize=18)

    def plot_trap_gammaDx(self, number, figtitle):

        Eti_avg = {'InP': {'SRH': 0, 'TAT': 0}, 'InGaAs': {'SRH': 0, 'TAT': 0}}
        if 'IT' in self.method:
            for material in self.material:
                for model in self.mechanism:
                    smooth_Eti = utils.center_smooth([self.result['IT'][material][model][V][2] for V in self.material_voltage[material]], 0)
                    for i, V in enumerate(self.material_voltage[material]):
                        Eti_avg[material][model] += smooth_Eti[i] / len(self.material_voltage[material])

        fig_Eti_TAT = plt.figure(number)
        plt.subplot(1, 2, 1)
        for i, item in enumerate(self.material):
            for j, model in enumerate(self.mechanism):
                plt.plot(T_analysis, [self.result['IV'][item][model][T][3] for T in T_analysis], marker='o',
                         fillstyle='none', linestyle=LineSet2[j], color=ColorSet10[-1 - i], label=r'[IV][%s][%s]' % (item, model))
                if 'IT' in self.method:
                    plt.plot(T_analysis, [Eti_avg[item][model] for T in T_analysis], color=ColorSet10[-3 - i],
                             linestyle=LineSet2[j], label=r'[IT][%s][%s]' % (item, model))
        plt.ylabel(r'E$_t$-E$_i$ (eV)')
        plt.legend(loc='best')
        plt.title('Trap level distribution')
        plt.ylim((-0.1, 0.18))
        plt.grid()

        plt.subplot(1, 2, 2)
        for i, T in enumerate(T_analysis):
            for j, item in enumerate(self.material):
                if item == 'InP':
                    plt.plot(*self.d_Gamma(T, item), color=ColorSet10[i], label='%s K' % T)
                else:
                    plt.plot(*self.d_Gamma(T, item), color=ColorSet10[i])
        plt.plot(*self.Plotdm(300, 'InP'), color='k', linestyle='-.', label='Depletion\nwidth')
        plt.plot(*self.Plotdm(330, 'InGaAs'), color='k', linestyle='-.')
        plt.grid()
        plt.xlabel('Voltage (V)')
        plt.ylabel(r'$d$ (cm)')
        plt.title(r'$d_{TAT}\equiv\int\Gamma(x)dx$')
        #plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
        #plt.xlim((-25, 2))
        plt.yscale('log')
        plt.legend(loc='best')
        fig_Eti_TAT.suptitle(figtitle, fontsize=18)

    def plot_lifetime_IV(self, number, figtitle):
        fig_Lifetime = plt.figure(number)
        for i, item in enumerate(self.material):
            plt.subplot(2, 2, i + 1)
            for model in self.mechanism:
                plt.plot(T_analysis, [1e9 * self.Lifetime[item][model][T][0] for T in T_analysis], marker='o',
                         fillstyle='none', linestyle='-', label=r'$\tau_p$[%s]' % model)
            plt.plot(300, Lifetime_ref[item][1]['p'], color='k', marker='o', fillstyle='none',
                     markersize=12, linestyle='none', label=Lifetime_ref[item][0])
            plt.grid(True)
            plt.legend(loc='best')
            plt.title(item)
            plt.ylabel(r'$\tau_p$ (ns)')

            plt.subplot(2, 2, i + 3)
            for model in self.mechanism:
                plt.plot(T_analysis, [1e9 * self.Lifetime[item][model][T][1] for T in T_analysis], marker='o',
                         fillstyle='none', linestyle='-', label=r'$\tau_n$[%s]' % model)
            plt.plot(300, Lifetime_ref[item][1]['n'], color='k', marker='o', fillstyle='none',
                     markersize=12, linestyle='none', label=Lifetime_ref[item][0])
            plt.grid(True)
            plt.legend(loc='best')
            plt.ylabel(r'$\tau_n$ (ns)')
            plt.xlabel('Temperature (K)')
        fig_Lifetime.suptitle(figtitle, fontsize=18)

    def plot_trap_IT_fitting(self, number, figtitle):
        fig_IT = plt.figure(number)
        plt.subplot(1, 2, 1)
        ColorModel = {'SRH': 'r', 'TAT': 'b'}
        for material in self.material:
            if material == 'InP':
                for model in self.mechanism:
                    plt.plot(self.material_voltage[material],
                             [self.result['IT'][material][model][V][2] for V in self.material_voltage[material]],
                             color=ColorModel[model], label=model)
            else:
                for model in self.mechanism:
                    plt.plot(self.material_voltage[material],
                             [self.result['IT'][material][model][V][2] for V in self.material_voltage[material]],
                             color=ColorModel[model])
        plt.grid()
        plt.xlabel('Voltage (V)')
        plt.ylabel(r'E$_t$-E$_i$ (eV)')
        plt.title('Trap level distribution')
        plt.legend(loc='best')

        fig_IT.suptitle(figtitle, fontsize=18)

    def plot_trap_IV_IT(self, number, figtitle):

        Eti_avg = {'InP': {'SRH': 0, 'TAT': 0}, 'InGaAs': {'SRH': 0, 'TAT': 0}}
        if 'IT' in self.method:
            for material in self.material:
                for model in self.mechanism:
                    smooth_Eti = utils.center_smooth([self.result['IT'][material][model][V][2] for V in self.material_voltage[material]], 0)
                    for i, V in enumerate(self.material_voltage[material]):
                        Eti_avg[material][model] += smooth_Eti[i] / len(self.material_voltage[material])

        fig_trap_IV_IT = plt.figure(number)
        plt.subplot(1, 2, 1)
        for i, item in enumerate(self.material):
            for j, model in enumerate(self.mechanism):
                plt.plot(T_analysis, [self.result['IV'][item][model][T][3] for T in T_analysis], marker='o',
                         fillstyle='none', linestyle=LineSet2[j], color=ColorSet10[-1 - i], label=r'[IV][%s][%s]' % (item, model))
                if 'IT' in self.method:
                    plt.plot(T_analysis, [Eti_avg[item][model] for T in T_analysis], color=ColorSet10[-3 - i],
                             linestyle=LineSet2[j], label=r'[IT][%s][%s]' % (item, model))
        plt.xlabel('Temperature (K)')
        plt.ylabel(r'E$_t$-E$_i$ (eV)')
        plt.legend(loc='best', ncol=2)
        plt.title('Trap level distribution')
        plt.ylim((-0.1, 0.18))
        plt.grid()

        plt.subplot(1, 2, 2)
        ColorModel = {'SRH': 'r', 'TAT': 'b'}
        for material in self.material:
            if material == 'InP':
                for model in self.mechanism:
                    plt.plot(self.material_voltage[material],
                             [self.result['IT'][material][model][V][2] for V in self.material_voltage[material]],
                             color=ColorModel[model], label=model)
            else:
                for model in self.mechanism:
                    plt.plot(self.material_voltage[material],
                             [self.result['IT'][material][model][V][2] for V in self.material_voltage[material]],
                             color=ColorModel[model])
        plt.grid()
        plt.xlabel('Voltage (V)')
        plt.ylabel(r'E$_t$-E$_t$ (eV)')
        plt.title('Trap level distribution')
        plt.legend(loc='best')
        fig_trap_IV_IT.suptitle(figtitle, fontsize=18)

    def plot_lifetime_IT_fitting(self, number, figtitle):
        fig_IT_info = plt.figure(number)
        plt.subplot(2, 2, 1)
        for material in self.material:
            for model in self.mechanism:
                if material == 'InP':
                    plt.plot(self.V_InP, [Lifetime_p['InP'] * self.result['IT'][material][model][V][3][0] for V in self.V_InP],
                             color=ColorModel[model], label=model)
                else:
                    plt.plot(self.V_InGaAs,
                             [Lifetime_p['InGaAs'] * self.result['IT'][material][model][V][3][0] for V in self.V_InGaAs],
                             color=ColorModel[model])
        plt.grid()
        plt.yscale('log')
        plt.xlabel(r'Voltage (V)')
        plt.ylabel(r'$\tau_{p, 300}$ (ns)')
        plt.legend(loc='best')

        plt.subplot(2, 2, 2)
        for material in self.material:
            for model in self.mechanism:
                if material == 'InP':
                    plt.plot(self.V_InP, [self.result['IT'][material][model][V][3][2] for V in self.V_InP],
                             color=ColorModel[model], label=model)
                else:
                    plt.plot(self.V_InGaAs, [self.result['IT'][material][model][V][3][2] for V in self.V_InGaAs],
                             color=ColorModel[model])
        plt.grid()
        plt.yscale('log')
        plt.xlabel(r'Voltage (V)')
        plt.ylabel(r'$\alpha_p$ (-)')
        plt.legend(loc='best')

        plt.subplot(2, 2, 3)
        for material in self.material:
            for model in self.mechanism:
                if material == 'InP':
                    plt.plot(self.V_InP, [Lifetime_n['InP'] * self.result['IT'][material][model][V][3][1] for V in self.V_InP],
                             color=ColorModel[model], label=model)
                else:
                    plt.plot(self.V_InGaAs,
                             [Lifetime_n['InGaAs'] * self.result['IT'][material][model][V][3][1] for V in self.V_InGaAs],
                             color=ColorModel[model])
        plt.grid()
        plt.yscale('log')
        plt.xlabel(r'Voltage (V)')
        plt.ylabel(r'$\tau_{n,300}$ (ns)')
        plt.legend(loc='best')

        plt.subplot(2, 2, 4)
        for material in self.material:
            for model in self.mechanism:
                if material == 'InP':
                    plt.plot(self.V_InP, [self.result['IT'][material][model][V][3][3] for V in self.V_InP],
                             color=ColorModel[model], label=model)
                else:
                    plt.plot(self.V_InGaAs, [self.result['IT'][material][model][V][3][3] for V in self.V_InGaAs],
                             color=ColorModel[model])
        plt.grid()
        plt.yscale('log')
        plt.xlabel(r'Voltage (V)')
        plt.ylabel(r'$\alpha_n$ (-)')
        plt.legend(loc='best')
        fig_IT_info.suptitle(figtitle, fontsize=18)

    def plot_lifetime_IT_IV_comparison(self, number, figtitle):
        def lifetime_p(T, tp, alpha_p, material):
            tau_p = tp * Lifetime_p[material] * 1e-9 * (T / 300) ** alpha_p
            return tau_p

        def lifetime_n(T, tn, alpha_n, material):
            tau_n = tn * Lifetime_n[material] * 1e-9 * (T / 300) ** alpha_n
            return tau_n

        def lifetime_approx(lifetime_array_at_T):
            lifetime_avg = 0
            num = len(lifetime_array_at_T)
            for i in range(num):
                lifetime_avg += (lifetime_array_at_T[i]) / num
            return lifetime_avg

        ColorModel = {'SRH': 'r', 'TAT': 'b'}
        self.Lifetime_p_DB = {'InP': {'SRH': {}, 'TAT': {}},
                         'InGaAs': {'SRH': {}, 'TAT': {}}}
        self.Lifetime_n_DB = {'InP': {'SRH': {}, 'TAT': {}},
                         'InGaAs': {'SRH': {}, 'TAT': {}}}
        fig_lifetime = plt.figure(number)
        for i, material in enumerate(self.material):
            plt.subplot(2, 2, i + 1)
            for j, model in enumerate(self.mechanism):
                LabelMarker1 = 0
                for V in self.material_voltage[material]:
                    tp = self.result['IT'][material][model][V][3][0]
                    alpha_p = self.result['IT'][material][model][V][3][2]
                    self.Lifetime_p_DB[material][model][V] = lifetime_p(T_analysis_IT, tp, alpha_p, material)  # [s]

                    if 'IV' not in self.method:
                        if LabelMarker1 == 0:
                            plt.plot(T_analysis_IT, 1e9 * self.Lifetime_p_DB[material][model][V],
                                     color=ColorModel[model], linestyle=':', label=model)
                            LabelMarker1 = 1
                        else:
                            plt.plot(T_analysis_IT, 1e9 * self.Lifetime_p_DB[material][model][V],
                                     color=ColorModel[model], linestyle=':')

                if 'IV' in self.method:
                    plt.plot(T_analysis, [1e9 * self.Lifetime[material][model][T][0] for T in T_analysis],
                             marker='o', color=ColorModel[model], fillstyle='none', linestyle='--',
                             label=r'[IV][%s]' % model)
                    Lifetime_avg_p = np.asarray([lifetime_approx([self.Lifetime_p_DB[material][model][V][i] for V
                                                                  in self.material_voltage[material]])
                                                 for i in range(len(T_analysis_IT))])
                    plt.plot(T_analysis_IT, 1e9 * Lifetime_avg_p, color=ColorModel[model], marker='o',
                             fillstyle='none', label=r'[IT][%s]' % model)
            #plt.plot(300, Lifetime_ref[material][1]['p'], color='k', marker='o', fillstyle='none',
            #         markersize=12, linestyle='none', label=Lifetime_ref[material][0])
            if material in self.TCAD_lifetime and self.TCAD_check == 'Yes':
                plt.plot([T for T in T_analysis if T in self.TCAD_lifetime[material]],
                         [1e9 * self.TCAD_lifetime[material][T]['p'] for T in T_analysis if T in self.TCAD_lifetime[material]],
                         marker='s', markersize=8, color='g', label='TCAD')
            plt.grid()
            plt.ylabel(r'$\tau_p$ (ns)')
            #plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
            plt.title('%s lifetime' % material)
            plt.legend(loc='best')

            plt.subplot(2, 2, i + 3)
            for j, model in enumerate(self.mechanism):
                LabelMarker2 = 0
                for V in self.material_voltage[material]:
                    tn = self.result['IT'][material][model][V][3][1]
                    alpha_n = self.result['IT'][material][model][V][3][3]
                    self.Lifetime_n_DB[material][model][V] = lifetime_n(T_analysis_IT, tn, alpha_n, material)  # [s]

                    if 'IV' not in self.method:
                        if LabelMarker2 == 0:
                            plt.plot(T_analysis_IT, 1e9 * self.Lifetime_n_DB[material][model][V],
                                     color=ColorModel[model], linestyle=':', label=model)
                            LabelMarker2 = 1
                        else:
                            plt.plot(T_analysis_IT, 1e9 * self.Lifetime_n_DB[material][model][V],
                                     color=ColorModel[model], linestyle=':')

                if 'IV' in self.method:
                    plt.plot(T_analysis, [1e9 * self.Lifetime[material][model][T][1] for T in T_analysis],
                             color=ColorModel[model], marker='o', fillstyle='none', linestyle='-.',
                             label=r'[IV][%s]' % model)
                    Lifetime_avg_n = np.asarray([lifetime_approx([self.Lifetime_n_DB[material][model][V][i] for V
                                                                  in self.material_voltage[material]])
                                                 for i in range(len(T_analysis_IT))])
                    plt.plot(T_analysis_IT, 1e9 * Lifetime_avg_n, color=ColorModel[model], marker='o',
                             fillstyle='none', label=r'[IT][%s]' % model)
            #plt.plot(300, Lifetime_ref[material][1]['n'], color='k', marker='o', fillstyle='none',
            #         markersize=12, linestyle='none', label=Lifetime_ref[material][0])
            if material in self.TCAD_lifetime and self.TCAD_check == 'Yes':
                plt.plot([T for T in T_analysis if T in self.TCAD_lifetime[material]],
                         [1e9 * self.TCAD_lifetime[material][T]['n'] for T in T_analysis if
                          T in self.TCAD_lifetime[material]],
                         marker='s', markersize=8, color='g', label='TCAD')
            plt.grid()
            plt.xlabel('Temperature (K)')
            plt.ylabel(r'$\tau_n$ (ns)')
            #plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
            plt.legend(loc='best')
        fig_lifetime.suptitle(figtitle, fontsize=18)

    def plot_current_fitting_comparison(self, number, figtitle):
        LabelMarker = 0
        temp = None
        fig_IV_IT = plt.figure(number)
        plt.subplot(1, 2, 1)
        for i, T in enumerate(T_analysis):
            '''
            for material in self.material:
                plt.plot(self.BTB[material][T][0], self.BTB[material][T][1], linestyle=':', color=ColorSet10[i])
            '''
            for j, model in enumerate(self.mechanism):
                for material in self.material:
                    if LabelMarker < 2 and temp != model:
                        temp = model
                        plt.plot(self.result['IV'][material][model][T][0], self.result['IV'][material][model][T][1],
                                 color=ColorSet10[i], linestyle=LineSet2[j], label=model)
                        LabelMarker += 1
                    else:
                        plt.plot(self.result['IV'][material][model][T][0], self.result['IV'][material][model][T][1],
                                 color=ColorSet10[i], linestyle=LineSet2[j])
            raw_V = [element for i, element in enumerate(self.RawIV[T].X) if i % count == 0]
            raw_I = [element for i, element in enumerate(self.RawIV[T].Y) if i % count == 0]
            plt.plot(raw_V, raw_I, linestyle='none', marker='o', fillstyle='none', color=ColorSet10[i], label=T)
        plt.vlines(x=-self.Vpt, ymin=1e-13, ymax=1e-6, color='k', linestyle=':', label=r'V$_{pt}$')
        plt.vlines(x=self.V1, ymin=1e-13, ymax=1e-9, color='k', linestyle='-.')
        plt.vlines(x=self.V2, ymin=1e-10, ymax=1e-6, color='k', linestyle='-.')
        plt.grid()
        plt.yscale('log')
        plt.xlabel('Voltage (V)')
        plt.ylabel('Current (A)')
        plt.title(r'Fixed-temperature fitting')
        plt.xlim((-45, 1))
        plt.ylim((1e-13, 1e-6))
        plt.legend(loc='best')
        plt.subplot(1, 2, 2)
        IT_InP = np.asarray([utils.find(self.RawIV[T].X, self.RawIV[T].Y, self.V1, 'log') for T in T_analysis_IT])
        IT_InGaAs = np.asarray([utils.find(self.RawIV[T].X, self.RawIV[T].Y, self.V2, 'log') for T in T_analysis_IT])
        plt.plot(1 / T_analysis_IT, np.log10(IT_InP), linestyle='none', marker='o', markersize='15',
                 fillstyle='none', label=r'$%s$ (V)' % self.V1)
        plt.plot(1 / T_analysis_IT, np.log10(IT_InGaAs), linestyle='none', marker='o', markersize='15',
                 fillstyle='none', label=r'$%s$ (V)' % self.V2)
        for material in self.material:
            for i, model in enumerate(self.mechanism):
                if material == 'InP':
                    plt.plot(1 / self.result['IT'][material][model][self.V1][0],
                             np.log10(self.result['IT'][material][model][self.V1][1]),
                             color=ColorModel[model], linestyle=LineSet2[i], label=model)
                else:
                    plt.plot(1 / self.result['IT'][material][model][self.V2][0],
                             np.log10(self.result['IT'][material][model][self.V2][1]),
                             linestyle=LineSet2[i], color=ColorModel[model])
        plt.grid()
        plt.xlabel(r'1/T (K$^{-1}$)')
        plt.ylabel('log10(I)')
        plt.legend(loc='best')
        plt.ylim((-13, -6))
        #plt.title(r'Fixed-voltage fitting')
        plt.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
        fig_IV_IT.suptitle(figtitle, fontsize=18)

    def plot_TCAD_python_fitting(self, number, figtitle):
        plt.figure(number)
        for i, T in enumerate(T_analysis):
            raw_V = [element for i, element in enumerate(self.RawIV[T].X) if i % count == 0]
            raw_I = [element for i, element in enumerate(self.RawIV[T].Y) if i % count == 0]
            plt.plot(raw_V, raw_I, linestyle='none', marker='o', fillstyle='none', color=ColorSet10[i], label=T)
            '''
            for material in self.material:
                plt.plot(self.BTB[material][T][0], self.BTB[material][T][1], linestyle=':', color=ColorSet10[i])
            '''
            for j, model in enumerate(self.mechanism):
                for material in self.material:
                    plt.plot(self.result['IV'][material][model][T][0], self.result['IV'][material][model][T][1],
                             color=ColorSet10[i], linestyle=LineSet2[j])
            if T in self.TCAD_IV:
                plt.plot(self.TCAD_IV[T].X, abs(self.TCAD_IV[T].Y), linestyle='--', color=ColorSet10[i], label='TCAD(%s)' % T)
        plt.vlines(x=-self.Vpt, ymin=1e-13, ymax=1e-6, color='k', linestyle='-.', label=r'V$_{pt}$')
        plt.grid()
        plt.yscale('log')
        plt.xlabel('Voltage (V)')
        plt.ylabel('Current (A)')
        plt.title(figtitle)
        plt.xlim((-45, 1))
        plt.ylim((1e-13, 1e-6))
        plt.legend(loc='best')

    def plot_TCAD_fitting(self, number, figtitle):
        fig_TCAD = plt.figure(number)
        for i, T in enumerate(T_analysis):
            raw_V = [element for i, element in enumerate(self.RawIV[T].X) if i % count == 0]
            raw_I = [element for i, element in enumerate(self.RawIV[T].Y) if i % count == 0]
            plt.plot(raw_V, raw_I, linestyle='none', marker='o', fillstyle='none', color=ColorSet10[i], label=T)
            if T in self.TCAD_IV:
                V_TCAD = [V for V in self.TCAD_IV[T].X if abs(V) <= 23]
                I_TCAD = [I for j, I in enumerate(abs(self.TCAD_IV[T].Y)) if abs(self.TCAD_IV[T].X[j]) <= 23]
                plt.plot(V_TCAD, I_TCAD, linestyle='-', color=ColorSet10[i])
        #plt.vlines(x=-self.Vpt, ymin=1e-13, ymax=1e-6, color='k', linestyle='-.', label=r'V$_{pt}$')
        plt.grid()
        plt.yscale('log')
        plt.xlabel('Voltage (V)')
        plt.ylabel('Current (A)')
        plt.title(figtitle)
        plt.xlim((-45, 1))
        plt.ylim((1e-14, 1e-6))
        plt.legend(loc='best')
        fig_TCAD.suptitle(figtitle, fontsize=18)

    def plot_double_IV(self, number, figtitle):
        LabelMarker = 0
        temp = None
        fig_double_IV = plt.figure(number)
        plt.subplot(1, 2, 1)
        for i, T in enumerate(T_analysis):
            '''
            for material in self.material:
                plt.plot(self.BTB[material][T][0], self.BTB[material][T][1], linestyle=':', color=ColorSet10[i])
            '''
            for j, model in enumerate(self.mechanism):
                for material in self.material:
                    if LabelMarker < 2 and temp != model:
                        temp = model
                        plt.plot(self.result['IV'][material][model][T][0], self.result['IV'][material][model][T][1],
                                 color=ColorSet10[i], linestyle=LineSet2[j], label=model)
                        LabelMarker += 1
                    else:
                        plt.plot(self.result['IV'][material][model][T][0], self.result['IV'][material][model][T][1],
                                 color=ColorSet10[i], linestyle=LineSet2[j])
            raw_V = [element for i, element in enumerate(self.RawIV[T].X) if i % count == 0]
            raw_I = [element for i, element in enumerate(self.RawIV[T].Y) if i % count == 0]
            plt.plot(raw_V, raw_I, linestyle='none', marker='o', fillstyle='none', color=ColorSet10[i], label=T)
        plt.vlines(x=-self.Vpt, ymin=1e-13, ymax=1e-6, color='k', linestyle=':', label=r'V$_{pt}$')
        plt.grid()
        plt.yscale('log')
        plt.xlabel('Voltage (V)')
        plt.ylabel('Current (A)')
        plt.title(r'Fixed-temperature fitting')
        plt.xlim((-45, 1))
        plt.ylim((1e-13, 1e-6))
        plt.legend(loc='best')
        plt.subplot(1, 2, 2)
        for i, T in enumerate(T_analysis):
            '''
            for material in self.material:
                plt.plot(self.BTB[material][T][0], self.BTB[material][T][1], linestyle=':', color=ColorSet10[i])
            '''
            for j, model in enumerate(self.mechanism):
                for material in self.material:
                    if LabelMarker < 2 and temp != model:
                        temp = model
                        plt.plot(self.result['IV'][material][model][T][0], self.result['IV'][material][model][T][1],
                                 color=ColorSet10[i], linestyle=LineSet2[j], label=model)
                        LabelMarker += 1
                    else:
                        plt.plot(self.result['IV'][material][model][T][0], self.result['IV'][material][model][T][1],
                                 color=ColorSet10[i], linestyle=LineSet2[j])
            raw_V = [element for i, element in enumerate(self.RawIV[T].X) if i % count == 0]
            raw_I = [element for i, element in enumerate(self.RawIV[T].Y) if i % count == 0]
            plt.plot(raw_V, raw_I, linestyle='none', marker='o', fillstyle='none', color=ColorSet10[i], label=T)
        plt.vlines(x=-self.Vpt, ymin=1e-13, ymax=1e-6, color='k', linestyle=':', label=r'V$_{pt}$')
        plt.grid()
        plt.yscale('log')
        plt.xlabel('Voltage (V)')
        plt.ylabel('Current (A)')
        plt.title(r'Fixed-temperature fitting')
        plt.xlim((-45, 1))
        plt.ylim((1e-13, 1e-6))
        plt.legend(loc='best')
        fig_double_IV.suptitle(figtitle, fontsize=18)

    def plot_any(self, number, figtitle):
        anyPlot = plt.figure(number)
        plt.subplot(1, 2, 1)
        for i, T in enumerate(T_analysis):
            raw_V = [element for i, element in enumerate(self.RawIV[T].X) if i % count == 0]
            raw_I = [element for i, element in enumerate(self.RawIV[T].Y) if i % count == 0]
            plt.plot(raw_V, raw_I, linestyle='none', marker='o', fillstyle='none', color=ColorSet10[i], label=T)
            if T in self.TCAD_IV:
                V_TCAD = [V for V in self.TCAD_IV[T].X if abs(V) <= 23]
                I_TCAD = [I for j, I in enumerate(abs(self.TCAD_IV[T].Y)) if abs(self.TCAD_IV[T].X[j]) <= 23]
                plt.plot(V_TCAD, I_TCAD, linestyle='-', color=ColorSet10[i])
        # plt.vlines(x=-self.Vpt, ymin=1e-13, ymax=1e-6, color='k', linestyle='-.', label=r'V$_{pt}$')
        plt.grid()
        plt.yscale('log')
        plt.xlabel('Voltage (V)')
        plt.ylabel('Current (A)')
        plt.title(figtitle)
        plt.xlim((-30, 0.5))
        plt.ylim((1e-14, 2e-9))
        plt.legend(loc='best')

        def lifetime_p(T, tp, alpha_p, material):
            tau_p = tp * Lifetime_p[material] * 1e-9 * (T / 300) ** alpha_p
            return tau_p

        def lifetime_n(T, tn, alpha_n, material):
            tau_n = tn * Lifetime_n[material] * 1e-9 * (T / 300) ** alpha_n
            return tau_n

        def lifetime_approx(lifetime_array_at_T):
            lifetime_avg = 0
            num = len(lifetime_array_at_T)
            for i in range(num):
                lifetime_avg += (lifetime_array_at_T[i]) / num
            return lifetime_avg

        ColorModel = {'SRH': 'r', 'TAT': 'b'}
        self.Lifetime_p_DB = {'InP': {'SRH': {}, 'TAT': {}},
                         'InGaAs': {'SRH': {}, 'TAT': {}}}
        self.Lifetime_n_DB = {'InP': {'SRH': {}, 'TAT': {}},
                         'InGaAs': {'SRH': {}, 'TAT': {}}}
        material = 'InP'
        plt.subplot(2, 2, 2)
        for j, model in enumerate(self.mechanism):
            LabelMarker1 = 0
            for V in self.material_voltage[material]:
                tp = self.result['IT'][material][model][V][3][0]
                alpha_p = self.result['IT'][material][model][V][3][2]
                self.Lifetime_p_DB[material][model][V] = lifetime_p(T_analysis_IT, tp, alpha_p, material)  # [s]

                if 'IV' not in self.method:
                    if LabelMarker1 == 0:
                        plt.plot(T_analysis_IT, 1e9 * self.Lifetime_p_DB[material][model][V],
                                 color=ColorModel[model], linestyle=':', label=model)
                        LabelMarker1 = 1
                    else:
                        plt.plot(T_analysis_IT, 1e9 * self.Lifetime_p_DB[material][model][V],
                                 color=ColorModel[model], linestyle=':')

            if 'IV' in self.method:
                plt.plot(T_analysis, [1e9 * self.Lifetime[material][model][T][0] for T in T_analysis],
                         marker='o', color=ColorModel[model], fillstyle='none', linestyle='--',
                         label=r'[IV][%s]' % model)
                Lifetime_avg_p = np.asarray([lifetime_approx([self.Lifetime_p_DB[material][model][V][i] for V
                                                              in self.material_voltage[material]])
                                             for i in range(len(T_analysis_IT))])
                plt.plot(T_analysis_IT, 1e9 * Lifetime_avg_p, color=ColorModel[model], marker='o',
                         fillstyle='none', label=r'[IT][%s]' % model)
        plt.plot(300, Lifetime_ref[material][1]['p'], color='k', marker='o', fillstyle='none',
                 markersize=12, linestyle='none', label=Lifetime_ref[material][0])
        if material in self.TCAD_lifetime and self.TCAD_check == 'Yes':
            plt.plot([T for T in T_analysis if T in self.TCAD_lifetime[material]],
                     [1e9 * self.TCAD_lifetime[material][T]['p'] for T in T_analysis if
                      T in self.TCAD_lifetime[material]],
                     marker='s', markersize=8, color='g', label='TCAD')
        plt.grid()
        plt.ylabel(r'$\tau_p$ (ns)')
        # plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
        plt.title('%s lifetime' % material)
        plt.legend(loc='upper left', ncol=3, fontsize=10)
        plt.ylim((0, 25))

        plt.subplot(2, 2, 4)
        for j, model in enumerate(self.mechanism):
            LabelMarker2 = 0
            for V in self.material_voltage[material]:
                tn = self.result['IT'][material][model][V][3][1]
                alpha_n = self.result['IT'][material][model][V][3][3]
                self.Lifetime_n_DB[material][model][V] = lifetime_n(T_analysis_IT, tn, alpha_n, material)  # [s]

                if 'IV' not in self.method:
                    if LabelMarker2 == 0:
                        plt.plot(T_analysis_IT, 1e9 * self.Lifetime_n_DB[material][model][V],
                                 color=ColorModel[model], linestyle=':', label=model)
                        LabelMarker2 = 1
                    else:
                        plt.plot(T_analysis_IT, 1e9 * self.Lifetime_n_DB[material][model][V],
                                 color=ColorModel[model], linestyle=':')

            if 'IV' in self.method:
                plt.plot(T_analysis, [1e9 * self.Lifetime[material][model][T][1] for T in T_analysis],
                         color=ColorModel[model], marker='o', fillstyle='none', linestyle='-.',
                         label=r'[IV][%s]' % model)
                Lifetime_avg_n = np.asarray([lifetime_approx([self.Lifetime_n_DB[material][model][V][i] for V
                                                              in self.material_voltage[material]])
                                             for i in range(len(T_analysis_IT))])
                plt.plot(T_analysis_IT, 1e9 * Lifetime_avg_n, color=ColorModel[model], marker='o',
                         fillstyle='none', label=r'[IT][%s]' % model)
        plt.plot(300, Lifetime_ref[material][1]['n'], color='k', marker='o', fillstyle='none',
                 markersize=12, linestyle='none', label=Lifetime_ref[material][0])
        if material in self.TCAD_lifetime and self.TCAD_check == 'Yes':
            plt.plot([T for T in T_analysis if T in self.TCAD_lifetime[material]],
                     [1e9 * self.TCAD_lifetime[material][T]['n'] for T in T_analysis if
                      T in self.TCAD_lifetime[material]],
                     marker='s', markersize=8, color='g', label='TCAD')
        plt.grid()
        plt.xlabel('Temperature (K)')
        plt.ylabel(r'$\tau_n$ (ns)')
        # plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
        plt.legend(loc='upper left', ncol=3, fontsize=10)
        plt.ylim((0, 1.75))
        anyPlot.suptitle(figtitle, fontsize=18)

    def plot_Python_only(self, number, figtitle):
        LabelMarker = 0
        temp = None
        fig_Python_only = plt.figure(number)
        plt.subplot(1, 2, 1)
        for i, T in enumerate(T_analysis):
            '''
            for material in self.material:
                plt.plot(self.BTB[material][T][0], self.BTB[material][T][1], linestyle=':', color=ColorSet10[i])
            '''
            for j, model in enumerate(self.mechanism):
                for material in self.material:
                    if LabelMarker < 2 and temp != model:
                        temp = model
                        plt.plot(self.result['IV'][material][model][T][0], self.result['IV'][material][model][T][1],
                                 color=ColorSet10[i], linestyle=LineSet2[j], label=model)
                        LabelMarker += 1
                    else:
                        plt.plot(self.result['IV'][material][model][T][0], self.result['IV'][material][model][T][1],
                                 color=ColorSet10[i], linestyle=LineSet2[j])
            raw_V = [element for i, element in enumerate(self.RawIV[T].X) if i % count == 0]
            raw_I = [element for i, element in enumerate(self.RawIV[T].Y) if i % count == 0]
            plt.plot(raw_V, raw_I, linestyle='none', marker='o', fillstyle='none', color=ColorSet10[i], label=T)
        #plt.vlines(x=-self.Vpt, ymin=1e-13, ymax=1e-6, color='k', linestyle=':', label=r'V$_{pt}$')
        #plt.vlines(x=self.V1, ymin=1e-13, ymax=1e-9, color='k', linestyle='-.')
        #plt.vlines(x=self.V2, ymin=1e-10, ymax=1e-6, color='k', linestyle='-.')
        plt.grid()
        plt.yscale('log')
        plt.xlabel('Voltage (V)')
        plt.ylabel('Current (A)')
        plt.title(r'$R_{TAT}\approx(1+\Gamma_p)R_{SRH}$ ; $m_h=0.86m_0$')
        plt.xlim((-30, 1))
        plt.ylim((1e-14, 1e-8))
        plt.legend(loc='best')
        plt.subplot(1, 2, 2)
        IT_InP = np.asarray([utils.find(self.RawIV[T].X, self.RawIV[T].Y, self.V1, 'log') for T in T_analysis_IT])
        IT_InGaAs = np.asarray([utils.find(self.RawIV[T].X, self.RawIV[T].Y, self.V2, 'log') for T in T_analysis_IT])
        plt.plot(1 / T_analysis_IT, np.log10(IT_InP), linestyle='none', marker='o', markersize='15',
                 fillstyle='none', label=r'$%s$ (V)' % self.V1)
        plt.plot(1 / T_analysis_IT, np.log10(IT_InGaAs), linestyle='none', marker='o', markersize='15',
                 fillstyle='none', label=r'$%s$ (V)' % self.V2)
        for material in self.material:
            for model in self.mechanism:
                if material == 'InP':
                    plt.plot(1 / self.result['IT'][material][model][self.V1][0],
                             np.log10(self.result['IT'][material][model][self.V1][1]),
                             color=ColorModel[model], label=model)
                else:
                    plt.plot(1 / self.result['IT'][material][model][self.V2][0],
                             np.log10(self.result['IT'][material][model][self.V2][1]),
                             color=ColorModel[model])
        plt.grid()
        plt.xlabel(r'1/T (K$^{-1}$)')
        plt.ylabel('log10(I)')
        plt.legend(loc='best')
        plt.ylim((-13, -6))
        plt.title(r'Fixed-voltage fitting')
        plt.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
        fig_Python_only.suptitle(figtitle, fontsize=18)

    def Lifetime_distribution(self, number, figtitle, xj, V, T):

        def Material_function(position):
            #  interface_um = [-3.62, -3.5, -0.5]
            if position < self.interface_um[0] or position > self.interface_um[2]:
                return 'InP'
            elif self.interface_um[0] <= position <= self.interface_um[1]:
                return 'InGaAsP'
            else:
                return 'InGaAs'

        def Gamma(material, carrier, F):
            # InGaAsP effective mass
            # http://www.ioffe.ru/SVA/NSM/Semicond/GaInAsP/bandstr.html#Masses
            # mn ~ 0.06，從圖近似
            # mp ~ 0.52494，從方程式近似(Goldberg Yu.A. & N.M. Schmidt (1999))
            mt = {'InP': {'n': 0.08, 'p': 0.86},
                  'InGaAs': {'n': 0.0463, 'p': 0.45},
                  'InGaAsP': {'n': 0.06, 'p': 0.52494}}
            F_Gamma = np.sqrt(24 * (mt[material][carrier] * me) * (kB * T) ** 3) / (e * h_bar) / 100  # [V/cm]
            return 2 * np.sqrt(3 * np.pi) * abs(F) / F_Gamma * np.exp((F / F_Gamma) ** 2)

        E_Vcm_array = self.DopingProfile.FieldTransform(self.Em_InP(V), [xj, *self.interface_um])
        Lifetime_n_array = np.asarray([self.Lifetime[Material_function(pos)]['TAT'][T][1] for pos in self.DopingProfile.X])
        Lifetime_p_array = np.asarray([self.Lifetime[Material_function(pos)]['TAT'][T][0] for pos in self.DopingProfile.X])
        Gamma_Distribution = {'n': np.asarray([Gamma(Material_function(self.DopingProfile.X[i]), 'n', Field)
                                               for i, Field in enumerate(E_Vcm_array)]),
                              'p': np.asarray([Gamma(Material_function(self.DopingProfile.X[i]), 'p', Field)
                                               for i, Field in enumerate(E_Vcm_array)])}
        fig_lifetime = plt.figure(number)
        plt.plot(self.DopingProfile.X, Lifetime_n_array, '-.b', label=r'$\tau_{n0}$')
        plt.plot(self.DopingProfile.X, Lifetime_p_array, '-.r', label=r'$\tau_{p0}$')
        plt.plot(self.DopingProfile.X, Lifetime_n_array / (1 + Gamma_Distribution['n']), '-b', label=r'$\tau_n$')
        plt.plot(self.DopingProfile.X, Lifetime_p_array / (1 + Gamma_Distribution['p']), '-r', label=r'$\tau_p$')
        plt.grid()
        plt.xlim((-4.25, -0.25))
        plt.legend(loc='best')
        plt.xlabel(r'Position ($\mu$ m)')
        plt.ylabel(r'Lifetime (s)')
        plt.yscale('log')
        fig_lifetime.suptitle(figtitle, fontsize=18)

    def plot_effective_mass_comparison(self, number, figtitle):
        mass_type = dict()
        Gamma_type = dict()
        if self.effective_mass_InP == 0.86:
            mass_type['InP'] = 'm_n'
            Gamma_type['InP'] = 'n'
        else:
            mass_type['InP'] = 'm_p'
            Gamma_type['InP'] = 'p'

        if self.effective_mass_InGaAs == 0.0463:
            mass_type['InGaAs'] = 'm_n'
            Gamma_type['InGaAs'] = 'n'
        else:
            mass_type['InGaAs'] = 'm_p'
            Gamma_type['InGaAs'] = 'p'

        LabelMarker = 0
        temp = None
        fig_effective_mass_comparison = plt.figure(number)
        plt.subplot(1, 2, 1)
        for i, T in enumerate(T_analysis):
            for j, model in enumerate(self.mechanism):
                if LabelMarker < 2 and temp != model:
                    temp = model
                    plt.plot(self.result['IV']['InP'][model][T][0], self.result['IV']['InP'][model][T][1],
                             color=ColorSet10[i], linestyle=LineSet2[j], label=model)
                    LabelMarker += 1
                else:
                    plt.plot(self.result['IV']['InP'][model][T][0], self.result['IV']['InP'][model][T][1],
                             color=ColorSet10[i], linestyle=LineSet2[j])
            raw_V = [element for i, element in enumerate(self.RawIV[T].X) if i % count == 0]
            raw_I = [element for i, element in enumerate(self.RawIV[T].Y) if i % count == 0]
            plt.plot(raw_V, raw_I, linestyle='none', marker='o', fillstyle='none', color=ColorSet10[i], label=T)
        plt.grid()
        plt.yscale('log')
        plt.xlabel('Voltage (V)')
        plt.ylabel('Current (A)')
        plt.title(r'$R_{TAT}\approx(1+\Gamma_%s)R_{SRH}$ ; $%s=%sm_0$' % (Gamma_type['InP'], mass_type['InP'], self.effective_mass_InP))
        plt.xlim((-30, 1))
        plt.ylim((1e-14, 1e-8))
        plt.legend(loc='best')

        LabelMarker = 0
        temp = None
        plt.subplot(1, 2, 2)
        for i, T in enumerate(T_analysis):
            for j, model in enumerate(self.mechanism):
                if LabelMarker < 2 and temp != model:
                    temp = model
                    plt.plot(self.result['IV']['InGaAs'][model][T][0], self.result['IV']['InGaAs'][model][T][1],
                             color=ColorSet10[i], linestyle=LineSet2[j], label=model)
                    LabelMarker += 1
                else:
                    plt.plot(self.result['IV']['InGaAs'][model][T][0], self.result['IV']['InGaAs'][model][T][1],
                             color=ColorSet10[i], linestyle=LineSet2[j])
            raw_V = [element for i, element in enumerate(self.RawIV[T].X) if i % count == 0]
            raw_I = [element for i, element in enumerate(self.RawIV[T].Y) if i % count == 0]
            plt.plot(raw_V, raw_I, linestyle='none', marker='o', fillstyle='none', color=ColorSet10[i], label=T)
        plt.grid()
        plt.yscale('log')
        plt.xlabel('Voltage (V)')
        plt.ylabel('Current (A)')
        plt.title(r'$R_{TAT}\approx(1+\Gamma_%s)R_{SRH}$ ; $%s=%sm_0$' % (Gamma_type['InGaAs'], mass_type['InGaAs'], self.effective_mass_InGaAs))
        plt.xlim((-45, -20))
        plt.ylim((1e-12, 1e-6))
        plt.legend(loc='best')
        fig_effective_mass_comparison.suptitle(figtitle, fontsize=18)