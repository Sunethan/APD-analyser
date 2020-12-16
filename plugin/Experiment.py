import numpy as np
import DataAnalysis as Data
import utils
from scipy.optimize import curve_fit
import csv
import physics as phys
import GenerationRate.BandToBandTunneling as BTB

Eg300_InGaAs = 0.744  # [eV]
Eg300_InP = 1.35  # [eV]
q = 1.6e-19  # [C]
kB = 1.38e-23


class IV_fitting(object):
    def __init__(self, settings, any_input=None, thin_film='Yes'):
        """
        :param settings:
        :param others:
        """
        self.kB, self.e, self.eps_InP, self.eps_InGaAs, self.eps_InGaAsP, self.hbar, self.Vpt, \
        self.Vpt_dec, self.Vabs_net, self.Location, self.A, self.EtY, self.Temperature, \
        self.V0, self.Vmax_InP, self.V0_TAT, self.step, self.iterations, self.Voltage, self.Tmin, \
        self.Tmax, self.T_analysis, self.T_analysis_IT, self.ABthick1, self.sdevice, self.snmesh, \
        self.Location_Ef, self.Label_InP, self.Label_InGaAs, self.Vf, self.Vf2, self.Location_Em_InP_NoAbthick, \
        self.Location_Em_InGaAs_NoAbthick, self.Label_InP_NoAbthick, self.Label_InGaAs_NoAbthick, \
        self.Location_Doping, self.Location_Doping_with_film, self.Label_Doping, self.Label_Doping_with_film, \
        self.interface_um, self.RawZn, self.Zn_xlabel, self.Zn_ylabel, self.d_mul = settings

        # 讀取IV，這裡必須給出 RawIV，不論TCAD還是實驗。
        self.RawIV = dict()
        LocationNew = {T: 'data/raw data/T' + str(T) + '/T.csv' for T in self.Temperature}
        fieldnames = ['DataName', 'V1', 'I1', 'I2', 'I11', 'I4']
        for T in self.Temperature:
            tempVoltage = []
            tempCurrent = []
            with open(LocationNew[T], newline='', encoding='utf-8-sig') as csv_file:
                rows = csv.DictReader(csv_file, fieldnames=fieldnames)
                for _index, row in enumerate(rows):
                    if row['DataName'] == 'DataValue':
                        tempVoltage.append(float(row['V1']))
                        tempCurrent.append(float(row['I11']))
                self.RawIV[T] = Data.CSV(data=None, xlabel=None, ylabel=None, x=tempVoltage, y=tempCurrent)

        # 讀取 InP & InGaAs 最大電場與偏壓的分佈
        self.film = thin_film
        if any_input is None:
            if self.film == 'Yes':
                if self.RawZn is None:
                    self.Ef_InP = Data.CSV(self.Location_Ef, xlabel=self.Label_InP, ylabel=self.Label_InP)
                    self.Ef_InGaAs = Data.CSV(self.Location_Ef, xlabel=self.Label_InGaAs, ylabel=self.Label_InGaAs)
                    self.DopingProfile = Data.DopingProfile(self.Location_Doping_with_film, self.Label_Doping_with_film,
                                                            self.Label_Doping_with_film)
                else:
                    self.Ef_InP = Data.CSV(self.Location_Ef, xlabel=self.Label_InP, ylabel=self.Label_InP)
                    self.Ef_InGaAs = Data.CSV(self.Location_Ef, xlabel=self.Label_InGaAs, ylabel=self.Label_InGaAs)
                    self.DopingProfile = Data.DopingProfile(self.Location_Doping_with_film, self.Label_Doping_with_film,
                                                            self.Label_Doping_with_film)
            elif self.film == 'No':
                self.Ef_InP = Data.CSV(self.Location_Em_InP_NoAbthick,
                                       xlabel=self.Label_InP_NoAbthick, ylabel=self.Label_InP_NoAbthick)
                self.Ef_InGaAs = Data.CSV(self.Location_Em_InGaAs_NoAbthick,
                                          xlabel=self.Label_InGaAs_NoAbthick, ylabel=self.Label_InGaAs_NoAbthick)
                self.DopingProfile = Data.DopingProfile(self.Location_Doping, self.Label_Doping, self.Label_Doping)
            else:
                raise BaseException("Wrong thin_film: %s" % self.film)
        else:
            Location_Em_InP, Location_Em_InGaAs, Label_InP, Label_InGaAs, DopingProfile, Label_doping, \
            self.Vpt, self.Vpt_dec, self.Vabs_net, self.Vmax_InP = any_input
            self.Ef_InP = Data.CSV(Location_Em_InP, Label_InP, Label_InP)
            self.Ef_InGaAs = Data.CSV(Location_Em_InGaAs, Label_InGaAs, Label_InGaAs)
            self.DopingProfile = Data.DopingProfile(DopingProfile, Label_doping, Label_doping)

    def RawData(self, T):
        return self.RawIV[T]

    def dm_InP(self, E_Vcm, ND, ND_c, d_mul, d_charge):
        d = E_Vcm * self.eps_InP / (self.e * ND)  # [cm]
        if type(d) is np.ndarray:
            dm_list = []
            for i, x in enumerate(d):
                if x <= d_mul:
                    dm_list.append(x)
                else:
                    E2 = E_Vcm[i] - (self.e * ND * d_mul) / self.eps_InP
                    d2 = E2 * self.eps_InP / (self.e * ND_c)
                    if d2 <= d_charge:
                        dm_list.append(d_mul + d2)
                    else:
                        dm_list.append(d_mul + d_charge)
            return np.asarray(dm_list)  # [cm]
        else:
            if d <= d_mul:
                return d  # [cm]
            else:
                E2 = E_Vcm - (self.e * ND * d_mul) / self.eps_InP
                d2 = E2 * self.eps_InP / (self.e * ND_c)
                if d2 <= d_charge:
                    return d_mul + d2  # [cm]
                else:
                    return d_mul + d_charge  # [cm]

    def dm_InGaAs(self, E, ND_abs, d_abs):
        """
        分段電場
        :param E:
        :param ND_abs:
        :param d_abs: [cm]
        :return:
        """
        ND_g = 2e15
        d_g = 0.12e-4
        if self.film == 'Yes':
            ND_abs_film = 2e17  # [cm-3]
            ABthick1_cm = self.ABthick1 * 1e-4  # [cm]
            d = E * self.eps_InGaAs / (self.e * ND_abs_film)  # [cm]
            if type(d) is np.ndarray:
                dm_list = []
                for i, x in enumerate(d):
                    if x <= ABthick1_cm:
                        dm_list.append(x)
                    else:
                        E2 = E[i] - (self.e * ND_abs_film * ABthick1_cm) / self.eps_InGaAs
                        d2 = E2 * self.eps_InGaAs / (self.e * ND_abs)
                        if d2 <= (d_abs - ABthick1_cm):
                            dm_list.append(ABthick1_cm + d2)
                        else:
                            dm_list.append(d_abs)
                return np.asarray(dm_list)
            else:
                if d <= ABthick1_cm:
                    return d
                else:
                    E2 = E - (self.e * ND_abs_film * ABthick1_cm) / self.eps_InGaAs
                    d2 = E2 * self.eps_InGaAs / (self.e * ND_abs)
                    if d2 <= (d_abs - ABthick1_cm):
                        return ABthick1_cm + d2
                    else:
                        return d_abs
        else:
            d = E * self.eps_InGaAs / (self.e * ND_abs)
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
            V_InP = np.asarray([V for V in self.RawIV[T].X if -self.V0_TAT >= V > -self.Vmax_InP])
            F_InP = np.asarray([self.Em_InP(V) for V in V_InP])
            I_InP = np.asarray([I for i, I in enumerate(self.RawIV[T].Y) if self.RawIV[T].X[i] in V_InP])
            def lifetime(tp, tn):
                alpha = 1.5
                tau_p0 = 300e-9  # [s]
                tau_n0 = 2.89e-9  # [s]
                tau_p = tp * tau_p0 * (T / 300) ** alpha
                tau_n = tn * tau_n0 * (T / 300) ** alpha
                return tau_p, tau_n
            if type == 'TAT':
                def TAT_InP_IV(X, Eti, tp, tn):
                    Emax_Vcm, T = X
                    alpha = 1.5
                    # tp = 1
                    # tn = 0.1
                    mt = 0.08
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
                    d_mul = self.d_mul  # 0.42e-4  # [cm]
                    d_ch = 0.2e-4  # [cm]

                    ni = np.sqrt(Nc300 * Nv300) * (T / 300) ** 1.5 * np.exp(- self.e * phys.Eg_InP(T) / (2 * kB * T))
                    G_SRH = ni / (2 * np.sqrt(tau_p * tau_n) * np.cosh(
                        self.e * Eti / (kB * T) + 0.5 * np.log(tau_p / tau_n)))

                    dM = self.dm_InP(Emax_Vcm, ND, Ncharge, d_mul, d_ch)  # 0.42e-4  # [cm]
                    F_Gamma = np.sqrt(24 * (mt * me) * (kB * T) ** 3) / (self.e * self.hbar) / 100  # [V/cm]
                    E1 = Emax_Vcm
                    log10_Current = []
                    for i, x in enumerate(dM):
                        if x <= d_mul:
                            E2 = E1[i] - (self.e * ND * x) / self.eps_InP
                            d_Gamma_1 = (np.sqrt(3 * np.pi) * self.eps_InP * F_Gamma) / (self.e * ND) * \
                                        (np.exp((E1[i] / F_Gamma) ** 2) - np.exp(E2 / F_Gamma ** 2))  # [cm]
                            log10_Current.append(
                                np.log10(self.A * self.e) + np.log10(prefactor * G_SRH) + np.log10(x + d_Gamma_1))
                        else:
                            E2 = E1[i] - (self.e * ND * d_mul) / self.eps_InP
                            E3 = E2 - (self.e * Ncharge * (x - d_mul)) / self.eps_InP
                            d_Gamma_1 = (np.sqrt(3 * np.pi) * self.eps_InP * F_Gamma) / (self.e * ND) * \
                                        (np.exp((E1[i] / F_Gamma) ** 2) - np.exp(E2 / F_Gamma ** 2))  # [cm]
                            d_Gamma_2 = (np.sqrt(3 * np.pi) * self.eps_InP * F_Gamma) / (self.e * Ncharge) * \
                                        (np.exp((E2 / F_Gamma) ** 2) - np.exp(E3 / F_Gamma ** 2))  # [cm]
                            log10_Current.append(
                                np.log10(self.A * self.e) + np.log10(prefactor * G_SRH) + np.log10(
                                    x + d_Gamma_1 + d_Gamma_2))
                    return np.asarray(log10_Current)

                TAT_InP_popt, TAT_InP_pcov = curve_fit(TAT_InP_IV, (F_InP, T), np.log10(I_InP), p0=guess, bounds=bound,
                                                       sigma=abs(np.log10(I_InP)) ** fitsigma)
                print('[TAT]   InP (%sK)  Eti:  %.3f,  tp: %.3e,  tn: %.3e' %
                      (T, TAT_InP_popt[0], TAT_InP_popt[1], TAT_InP_popt[2]))
                Eti = TAT_InP_popt[0]
                mt = 0.08
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
                    tau_p0 = 300e-9  # [s]
                    tau_n0 = 2.89e-9  # [s]
                    tau_p = tp * tau_p0 * (T / 300) ** alpha
                    tau_n = tn * tau_n0 * (T / 300) ** alpha
                    ND = 5e16  # [cm-3]
                    Ncharge = 7.8e16  # [cm-3]
                    d_mul = self.d_mul  # 0.42e-4  # [cm]
                    d_ch = 0.2e-4  # [cm]

                    ni = np.sqrt(Nc300 * Nv300) * (T / 300) ** 1.5 * np.exp(- self.e * phys.Eg_InP(T) / (2 * kB * T))
                    G_SRH = ni / (
                                2 * np.sqrt(tau_p * tau_n) * np.cosh(self.e * Eti / (kB * T) + 0.5 * np.log(tau_p / tau_n)))
                    dM = self.dm_InP(Emax_Vcm, ND, Ncharge, d_mul, d_ch)  # 0.42e-4  # [cm]
                    return np.log10(self.A * self.e) + np.log10(prefactor * G_SRH) + np.log10(dM)

                popt_SRH_InP, pcov_SRH_InP = curve_fit(SRH_InP, (F_InP, T), np.log10(I_InP), p0=guess, bounds=bound,
                                                       sigma=abs(np.log10(I_InP)) ** fitsigma)
                print('[SRH]   InP (%sK)  Eti:  %.3f,  tp: %.3e,  tn: %.3e' %
                      (T, popt_SRH_InP[0], popt_SRH_InP[1], popt_SRH_InP[2]))
                Eti = popt_SRH_InP[0]
                mt = 0.08
                tau_p, tau_n = lifetime(popt_SRH_InP[1], popt_SRH_InP[2])
                return V_InP, 10 ** SRH_InP((F_InP, T), *popt_SRH_InP), [tau_p, tau_n], Eti, mt
            else:
                raise BaseException("Wrong type: %s" % type)
        elif material == 'InGaAs':
            V_InGaAs = np.asarray([V for V in self.RawIV[T].X if -self.Vf2[T] <= V <= -self.Vabs_net])
            F_InGaAs = np.asarray([self.Em_InGaAs(V) for V in V_InGaAs])
            I_InGaAs = np.asarray([I for i, I in enumerate(self.RawIV[T].Y) if self.RawIV[T].X[i] in V_InGaAs])
            def lifetime(tp, tn):
                alpha = 1.5
                tau_p0 = 8e-9  # [s]
                tau_n0 = 0.25e-9  # [s]
                tau_p = tp * tau_p0 * (T / 300) ** alpha
                tau_n = tn * tau_n0 * (T / 300) ** alpha
                return tau_p, tau_n
            if type == 'TAT':
                def TAT_InGaAs_IV(X, Eti, tp, tn):
                    Emax_Vcm, T = X
                    prefactor = 1
                    # tp = 1
                    # tn = 1
                    mt = 0.042
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

                    ni = np.sqrt(Nc300 * Nv300) * (T / 300) ** 1.5 * np.exp(- self.e * phys.Eg_InGaAs(T) / (2 * kB * T))
                    G_SRH = ni / (2 * np.sqrt(tau_p * tau_n) * np.cosh(
                        self.e * Eti / (kB * T) + 0.5 * np.log(tau_p / tau_n)))

                    dM = self.dm_InGaAs(Emax_Vcm, ND_abs, d_InGaAs)  # [cm]
                    F_Gamma = np.sqrt(24 * (mt * me) * (kB * T) ** 3) / (self.e * self.hbar) / 100  # [V/cm]
                    E1 = Emax_Vcm
                    E2 = 0
                    d_Gamma = (np.sqrt(3 * np.pi) * self.eps_InGaAs * F_Gamma) / (self.e * ND_abs) * \
                              (np.exp((E1 / F_Gamma) ** 2) - np.exp((E2 / F_Gamma) ** 2))  # [cm]
                    return np.log10(self.A * self.e) + np.log10(prefactor * G_SRH) + np.log10(dM + d_Gamma)

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
                    mt = 0.042
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
                    tau_p0 = 8e-9  # [s]
                    tau_n0 = 0.25e-9  # [s]
                    tau_p = tp * tau_p0 * (T / 300) ** alpha
                    tau_n = tn * tau_n0 * (T / 300) ** alpha
                    ND_abs = 7.5e14  # [cm-3]
                    d_InGaAs = 3e-4  # [cm]

                    ni = np.sqrt(Nc300 * Nv300) * (T / 300) ** 1.5 * np.exp(- self.e * phys.Eg_InGaAs(T) / (2 * self.kB * T))
                    G_SRH = ni / (
                                2 * np.sqrt(tau_p * tau_n) * np.cosh(self.e * Eti / (kB * T) + 0.5 * np.log(tau_p / tau_n)))

                    dM = self.dm_InGaAs(Emax_Vcm, ND_abs, d_InGaAs)  # [cm]
                    return np.log10(self.A * self.e) + np.log10(prefactor * G_SRH) + np.log10(dM)

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
            I_InP = np.asarray([utils.find(self.RawIV[T].X, self.RawIV[T].Y, V, 'log') for T in self.T_analysis_IT])
            if type == 'TAT':
                def TAT_InP_IT(X, Eti, tp, tn, alpha_p, alpha_n):
                    T, Emax_Vcm = X
                    mt = 0.185
                    prefactor = 1

                    me = 9.11e-31
                    Nc300 = 5.716e17  # [cm-3]
                    Nv300 = 1.143e19  # [cm-3]
                    tau_p0 = 300e-9  # [s]
                    tau_n0 = 2.89e-9  # [s]
                    tau_p = tp * tau_p0 * (T / 300) ** alpha_p
                    tau_n = tn * tau_n0 * (T / 300) ** alpha_n
                    ND = 5e16  # [cm-3]
                    Ncharge = 7.8e16  # [cm-3]
                    d_mul = self.d_mul  # 0.42e-4  # [cm]
                    d_ch = 0.2e-4  # [cm]

                    ni = np.sqrt(Nc300 * Nv300) * (T / 300) ** 1.5 * np.exp(- self.e * phys.Eg_InP(T) / (2 * kB * T))
                    G_SRH = ni / (2 * np.sqrt(tau_p * tau_n) * np.cosh(
                        self.e * Eti / (kB * T) + 0.5 * np.log(tau_p / tau_n)))

                    dM = self.dm_InP(Emax_Vcm, ND, Ncharge, d_mul, d_ch)  # 0.42e-4  # [cm]
                    F_Gamma = np.sqrt(24 * (mt * me) * (kB * T) ** 3) / (self.e * self.hbar) / 100  # [V/cm]
                    E1 = Emax_Vcm
                    if dM <= d_mul:
                        E2 = E1 - (self.e * ND * dM) / self.eps_InP
                        d_Gamma_1 = (np.sqrt(3 * np.pi) * self.eps_InP * F_Gamma) / (self.e * ND) * \
                                    (np.exp((E1 / F_Gamma) ** 2) - np.exp(E2 / F_Gamma ** 2))  # [cm]
                        return np.log10(self.A * self.e) + np.log10(prefactor * G_SRH) + np.log10(dM + d_Gamma_1)
                    else:
                        E2 = E1 - (self.e * ND * d_mul) / self.eps_InP
                        E3 = E2 - (self.e * Ncharge * (dM - d_mul)) / self.eps_InP
                        d_Gamma_1 = (np.sqrt(3 * np.pi) * self.eps_InP * F_Gamma) / (self.e * ND) * \
                                    (np.exp((E1 / F_Gamma) ** 2) - np.exp(E2 / F_Gamma ** 2))  # [cm]
                        d_Gamma_2 = (np.sqrt(3 * np.pi) * self.eps_InP * F_Gamma) / (self.e * Ncharge) * \
                                    (np.exp((E2 / F_Gamma) ** 2) - np.exp(E3 / F_Gamma ** 2))  # [cm]
                        return np.log10(self.A * self.e) + np.log10(prefactor * G_SRH) + np.log10(dM + d_Gamma_1 + d_Gamma_2)

                popt, pcov = curve_fit(TAT_InP_IT, (self.T_analysis_IT, self.Em_InP(V)), np.log10(I_InP), p0=guess,
                                       bounds=bound, sigma=abs(np.log10(I_InP)) ** fitsigma)
                Eti, tp, tn, alpha_p, alpha_n = popt
                print('[TAT]   InP (%.1f)  Eti:  %.3f,  tp: %.3e,  tn: %.3e,  alpha(p): %.3e,  alpha(n): %.3e' %
                      (V, Eti, tp, tn, alpha_p, alpha_n))
                return self.T_analysis_IT, 10 ** TAT_InP_IT((self.T_analysis_IT, self.Em_InP(V)), *popt), \
                       Eti, [tp, tn, alpha_p, alpha_n]
            elif type == 'SRH':
                def SRH_InP_IT(X, Eti, tp, tn, alpha_n, alpha_p):
                    T, Emax_Vcm = X
                    alpha = 1.5
                    # tp = 1
                    # tn = 0.1
                    prefactor = 1

                    Nc300 = 5.716e17  # [cm-3]
                    Nv300 = 1.143e19  # [cm-3]
                    tau_p0 = 300e-9  # [s]
                    tau_n0 = 2.89e-9  # [s]
                    tau_p = tp * tau_p0 * (T / 300) ** alpha_p
                    tau_n = tn * tau_n0 * (T / 300) ** alpha_n
                    ND = 5e16  # [cm-3]
                    Ncharge = 7.8e16  # [cm-3]
                    d_mul = self.d_mul  # 0.42e-4  # [cm]
                    d_ch = 0.2e-4  # [cm]

                    ni = np.sqrt(Nc300 * Nv300) * (T / 300) ** 1.5 * np.exp(- self.e * phys.Eg_InP(T) / (2 * kB * T))
                    G_SRH = ni / (2 * np.sqrt(tau_p * tau_n) * np.cosh(
                        self.e * Eti / (kB * T) + 0.5 * np.log(tau_p / tau_n)))

                    dM = self.dm_InP(Emax_Vcm, ND, Ncharge, d_mul, d_ch)  # 0.42e-4  # [cm]
                    return np.log10(self.A * self.e) + np.log10(prefactor * G_SRH) + np.log10(dM)

                popt, pcov = curve_fit(SRH_InP_IT, (self.T_analysis_IT, self.Em_InP(V)), np.log10(I_InP), p0=guess,
                                       bounds=bound, sigma=abs(np.log10(I_InP)) ** fitsigma)
                Eti, tp, tn, alpha_p, alpha_n = popt
                print('[SRH]   InP (%.1f)  Eti:  %.3f,  tp: %.3e,  tn: %.3e,  alpha(p): %.3e,  alpha(n): %.3e' %
                      (V, Eti, tp, tn, alpha_p, alpha_n))
                return self.T_analysis_IT, 10 ** SRH_InP_IT((self.T_analysis_IT, self.Em_InP(V)), *popt), \
                       Eti, [tp, tn, alpha_p, alpha_n]
            else:
                raise BaseException("Wrong type: %s" % type)
        elif material == 'InGaAs':
            I_InGaAs = np.asarray([utils.find(self.RawIV[T].X, self.RawIV[T].Y, V, 'log') for T in self.T_analysis_IT])

            # 檢查電流是否隨著溫度遞增
            if abs(V) > abs(self.Vf2[self.T_analysis_IT[0]]):
                raise BaseException("Voltage is too large: %s > Vmax(InGaAs,240K) = %s" %
                                    (abs(V), abs(self.Vf2[self.T_analysis_IT[0]])))

            if type == 'TAT':
                def TAT_InGaAs_IT(X, Eti, tp, tn, alpha_p, alpha_n):
                    T, Emax_Vcm = X
                    prefactor = 1
                    mt = 0.042
                    me = 9.11e-31
                    Nc300 = 2.53956e17  # [cm-3]
                    Nv300 = 7.51e18  # [cm-3]
                    tau_p0 = 8e-9  # [s]
                    tau_n0 = 0.25e-9  # [s]
                    tau_p = tp * tau_p0 * (T / 300) ** alpha_p
                    tau_n = tn * tau_n0 * (T / 300) ** alpha_n
                    ND_abs = 7.5e14  # [cm-3]
                    d_InGaAs = 3e-4  # [cm]

                    ni = np.sqrt(Nc300 * Nv300) * (T / 300) ** 1.5 * np.exp(- self.e * phys.Eg_InGaAs(T) / (2 * kB * T))
                    G_SRH = ni / (2 * np.sqrt(tau_p * tau_n) * np.cosh(
                        self.e * Eti / (kB * T) + 0.5 * np.log(tau_p / tau_n)))

                    dM = self.dm_InGaAs(Emax_Vcm, ND_abs, d_InGaAs)  # [cm]
                    F_Gamma = np.sqrt(24 * (mt * me) * (kB * T) ** 3) / (self.e * self.hbar) / 100  # [V/cm]
                    E1 = Emax_Vcm
                    E2 = 0
                    d_Gamma = (np.sqrt(3 * np.pi) * self.eps_InGaAs * F_Gamma) / (self.e * ND_abs) * \
                              (np.exp((E1 / F_Gamma) ** 2) - np.exp((E2 / F_Gamma) ** 2))  # [cm]
                    return np.log10(self.A * self.e) + np.log10(prefactor * G_SRH) + np.log10(dM + d_Gamma)

                popt, pcov = curve_fit(TAT_InGaAs_IT, (self.T_analysis_IT, self.Em_InGaAs(V)), np.log10(I_InGaAs),
                                       p0=guess, bounds=bound, sigma=abs(np.log10(I_InGaAs)) ** fitsigma)
                Eti, tp, tn, alpha_p, alpha_n = popt
                print('[TAT]   InGaAs (%.1f)  Eti:  %.3f,  tp: %.3e,  tn: %.3e,  alpha(p): %.3e,  alpha(n): %.3e' %
                      (V, Eti, tp, tn, alpha_p, alpha_n))
                return self.T_analysis_IT, 10 ** TAT_InGaAs_IT((self.T_analysis_IT, self.Em_InGaAs(V)), *popt), \
                       Eti, [tp, tn, alpha_p, alpha_n]
            elif type == 'SRH':
                def SRH_InGaAs_IT(X, Eti, tp, tn, alpha_p, alpha_n):
                    T, Emax_Vcm = X
                    prefactor = 1
                    Nc300 = 2.53956e17  # [cm-3]
                    Nv300 = 7.51e18  # [cm-3]
                    tau_p0 = 8e-9  # [s]
                    tau_n0 = 0.25e-9  # [s]
                    tau_p = tp * tau_p0 * (T / 300) ** alpha_p
                    tau_n = tn * tau_n0 * (T / 300) ** alpha_n
                    ND_abs = 7.5e14  # [cm-3]
                    d_InGaAs = 3e-4  # [cm]

                    ni = np.sqrt(Nc300 * Nv300) * (T / 300) ** 1.5 * np.exp(- self.e * phys.Eg_InGaAs(T) / (2 * kB * T))
                    G_SRH = ni / (2 * np.sqrt(tau_p * tau_n) * np.cosh(
                        self.e * Eti / (kB * T) + 0.5 * np.log(tau_p / tau_n)))

                    dM = self.dm_InGaAs(Emax_Vcm, ND_abs, d_InGaAs)  # [cm]
                    return np.log10(self.A * self.e) + np.log10(prefactor * G_SRH) + np.log10(dM)

                popt, pcov = curve_fit(SRH_InGaAs_IT, (self.T_analysis_IT, self.Em_InGaAs(V)), np.log10(I_InGaAs),
                                       p0=guess, bounds=bound, sigma=abs(np.log10(I_InGaAs)) ** fitsigma)
                Eti, tp, tn, alpha_p, alpha_n = popt
                print('[SRH]   InGaAs (%.1f)  Eti:  %.3f,  tp: %.3e,  tn: %.3e,  alpha(p): %.3e,  alpha(n): %.3e' %
                      (V, Eti, tp, tn, alpha_p, alpha_n))
                return self.T_analysis_IT, 10 ** SRH_InGaAs_IT((self.T_analysis_IT, self.Em_InGaAs(V)), *popt), \
                       Eti, [tp, tn, alpha_p, alpha_n]
            else:
                raise BaseException("Wrong type: %s" % type)
        else:
            raise BaseException("Wrong material: %s" % material)

    def d_Gamma(self, T, material):
        if material == 'InP':
            V_InP = np.asarray([V for V in self.RawIV[T].X if -self.V0_TAT >= V > -self.Vmax_InP])
            F_InP = np.asarray([self.Em_InP(V) for V in V_InP])
            def Gamma_InP(X):
                Emax_Vcm, T = X
                mt = 0.185

                me = 9.11e-31
                ND = 5e16  # [cm-3]
                Ncharge = 7.8e16  # [cm-3]
                d_mul = 0.42e-4  # [cm]
                d_ch = 0.2e-4  # [cm]

                dM = self.dm_InP(Emax_Vcm, ND, Ncharge, d_mul, d_ch)  # 0.42e-4  # [cm]
                F_Gamma = np.sqrt(24 * (mt * me) * (kB * T) ** 3) / (self.e * self.hbar) / 100  # [V/cm]
                E1 = Emax_Vcm
                d_Gamma = []
                for i, x in enumerate(dM):
                    if x <= d_mul:
                        E2 = E1[i] - (self.e * ND * x) / self.eps_InP
                        d_Gamma_1 = (np.sqrt(3 * np.pi) * self.eps_InP * F_Gamma) / (self.e * ND) * \
                                    (np.exp((E1[i] / F_Gamma) ** 2) - np.exp(E2 / F_Gamma ** 2))  # [cm]
                        d_Gamma.append(d_Gamma_1)
                    else:
                        E2 = E1[i] - (self.e * ND * d_mul) / self.eps_InP
                        E3 = E2 - (self.e * Ncharge * (x - d_mul)) / self.eps_InP
                        d_Gamma_1 = (np.sqrt(3 * np.pi) * self.eps_InP * F_Gamma) / (self.e * ND) * \
                                    (np.exp((E1[i] / F_Gamma) ** 2) - np.exp(E2 / F_Gamma ** 2))  # [cm]
                        d_Gamma_2 = (np.sqrt(3 * np.pi) * self.eps_InP * F_Gamma) / (self.e * Ncharge) * \
                                    (np.exp((E2 / F_Gamma) ** 2) - np.exp(E3 / F_Gamma ** 2))  # [cm]
                        d_Gamma.append(d_Gamma_1 + d_Gamma_2)
                return np.asarray(d_Gamma)
            return V_InP, Gamma_InP((F_InP, T))
        elif material == 'InGaAs':
            V_InGaAs = np.asarray([V for V in self.RawIV[T].X if -self.Vf2[T] <= V <= -self.Vabs_net])
            F_InGaAs = np.asarray([self.Em_InGaAs(V) for V in V_InGaAs])
            def Gamma_InGaAs(X):
                Emax_Vcm, T = X
                mt = 0.042
                me = 9.11e-31
                ND_abs = 7.5e14  # [cm-3]
                F_Gamma = np.sqrt(24 * (mt * me) * (kB * T) ** 3) / (self.e * self.hbar) / 100  # [V/cm]
                E1 = Emax_Vcm
                E2 = 0
                d_Gamma = (np.sqrt(3 * np.pi) * self.eps_InGaAs * F_Gamma) / (self.e * ND_abs) * \
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
                    V_BTB = np.asarray([V for V in self.RawIV[T].X if -self.V0_TAT >= V > -self.Vmax_InP])
                elif V_range == 'InGaAs':
                    V_BTB = np.asarray([V for V in self.RawIV[T].X if -self.Vf2[T] <= V <= -self.Vabs_net])
                else:
                    raise BaseException("Wrong parameter: %s" % parameter)
                I_BTB = []
                for V in V_BTB:
                    Ej = self.Em_InP(V)
                    E_Vcm_array = self.DopingProfile.FieldTransform(Ej, self.interface_um)
                    I_BTB.append(BTB.J_BTB_InP(self.DopingProfile.X * 1e-4, E_Vcm_array, T, effective_mass) * self.A)
                return V_BTB, I_BTB
        elif material == 'InGaAs':
            V_InGaAs = np.asarray([V for V in self.RawIV[T].X if -self.Vf2[T] <= V <= -self.Vabs_net])
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
                    tau_p0 = 8e-9  # [s]
                    tau_n0 = 0.25e-9  # [s]
                    tau_p = tp * tau_p0 * (T / 300) ** alpha
                    tau_n = tn * tau_n0 * (T / 300) ** alpha
                    ND_abs = 7.5e14  # [cm-3]
                    d_InGaAs = 3e-4  # [cm]

                    ni = np.sqrt(Nc300 * Nv300) * (T / 300) ** 1.5 * np.exp(
                        - self.e * phys.Eg_InGaAs(T) / (2 * self.kB * T))
                    G_SRH = ni / (
                            2 * np.sqrt(tau_p * tau_n) * np.cosh(self.e * Eti / (kB * T) + 0.5 * np.log(tau_p / tau_n)))

                    dM = self.dm_InGaAs(Emax_Vcm, ND_abs, d_InGaAs)  # [cm]
                    return np.log10(self.A * self.e) + np.log10(prefactor * G_SRH) + np.log10(dM)
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
                    tau_p0 = 8e-9  # [s]
                    tau_n0 = 0.25e-9  # [s]
                    tau_p = tp * tau_p0 * (T / 300) ** alpha
                    tau_n = tn * tau_n0 * (T / 300) ** alpha
                    ND_abs = 7.5e14  # [cm-3]
                    d_InGaAs = 3e-4  # [cm]

                    ni = np.sqrt(Nc300 * Nv300) * (T / 300) ** 1.5 * np.exp(- self.e * phys.Eg_InGaAs(T) / (2 * kB * T))
                    G_SRH = ni / (2 * np.sqrt(tau_p * tau_n) * np.cosh(
                        self.e * Eti / (kB * T) + 0.5 * np.log(tau_p / tau_n)))

                    dM = self.dm_InGaAs(Emax_Vcm, ND_abs, d_InGaAs)  # [cm]
                    F_Gamma = np.sqrt(24 * (mt * me) * (kB * T) ** 3) / (self.e * self.hbar) / 100  # [V/cm]
                    E1 = Emax_Vcm
                    E2 = 0
                    d_Gamma = (np.sqrt(3 * np.pi) * self.eps_InGaAs * F_Gamma) / (self.e * ND_abs) * (
                            np.exp((E1 / F_Gamma) ** 2) - np.exp((E2 / F_Gamma) ** 2))  # [cm]
                    return np.log10(self.A * self.e) + np.log10(prefactor * G_SRH) + np.log10(dM + d_Gamma)
                return V_InGaAs, 10 ** TAT_InGaAs_IV((F_InGaAs, T), *parameter)
            elif type == 'BTB':
                V_range, effective_mass = parameter
                if V_range == 'All':
                    V_BTB = self.RawIV[T].X
                elif V_range == 'InP':
                    V_BTB = np.asarray([V for V in self.RawIV[T].X if -self.V0_TAT >= V > -self.Vmax_InP])
                elif V_range == 'InGaAs':
                    V_BTB = np.asarray([V for V in self.RawIV[T].X if -self.Vf2[T] <= V <= -self.Vabs_net])
                else:
                    raise BaseException("Wrong parameter: %s" % parameter)
                I_BTB = []
                for V in V_BTB:
                    Ej = self.Em_InGaAs(V)
                    E_Vcm_array = self.DopingProfile.FieldTransform(Ej, self.interface_um)
                    I_BTB.append(BTB.J_BTB_InGaAs(self.DopingProfile.X * 1e-4, E_Vcm_array, T, effective_mass) * self.A)
                return V_BTB, I_BTB
            else:
                raise BaseException("Wrong type: %s" % type)
        else:
            raise BaseException("Wrong material: %s" % material)

    def PlotEm(self, T, material):
        V_net = np.asarray([V for V in self.RawIV[T].X if -self.V0_TAT >= V > -self.Vf])
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
            V_InP = np.asarray([V for V in self.RawIV[T].X if -self.V0_TAT >= V > -self.Vmax_InP])
            F_InP = np.asarray([self.Em_InP(V) for V in self.RawIV[T].X if -self.V0_TAT >= V > -self.Vmax_InP])
            return V_InP, self.dm_InP(F_InP, 5e16, 7.8e16, 0.42e-4, 0.2e-4)
        elif material == 'InGaAs':
            V_InGaAs = np.asarray([V for V in self.RawIV[T].X if -self.Vf2[T] <= V <= -self.Vabs_net])
            F_InGaAs = np.asarray([self.Em_InGaAs(V) for V in V_InGaAs])
            return V_InGaAs, self.dm_InGaAs(F_InGaAs, 7.5e14, 3e-4)
        else:
            raise BaseException("Wrong material: %s" % material)


class TemperatureDependentIV(object):
    def __init__(self, Location, TemperatureList, Vpt, area, source):
        """
        輸入變溫 IV 數據 (Generation/Recombination rate)
        :param Location: TCAD 數據位置
        :param TemperatureList: 溫度範圍
        """
        self.__Location = Location
        self.__TemperatureList = np.asarray(TemperatureList)
        self.__Vpt = Vpt
        self.__Area = area
        self.__Curve = utils.parse_csv(Location, 2, 0, 0, 'Yes')
        if source == 'Simulation':
            self.__Mechanism = [self.__Curve[i].split("(")[0] for i in range(len(self.__Curve) - 1)
                                if i == 0 or (self.__Curve[i - 1].split("(")[0] != self.__Curve[i].split("(")[0])]
            self.__Node = [int(self.__Curve[i].split("_")[1].split("n")[1]) for i in range(len(self.__Curve))
                           if self.__Mechanism[0] in self.__Curve[i]]
        self.__IV = None
        self.__Rate = dict()
        self.__IofT = dict()
        self.__Eti = dict()

    def get_IV(self):
        if self.__IV is None:
            self.__IV = [Data.IV(self.__Location, curve, curve, self.__Area)
                         for curve in self.__Curve if 'TotalCurrent' in curve]
        return self.__IV

    def get_mechanism(self):
        for i, element in enumerate(self.__Mechanism):
            print('%2s, %25s, %25s, %50s' % (i, element.split(" ")[0], element.split(" ")[1], element))
        return self.__Mechanism

    def get_recombination_rate(self, rate_name, source):
        """
        以TCAD而言，這部分的邏輯是「只要該再生機制名稱在其數據名稱中，那麼就收集起來」。
        例如說，rate_name = 'IntegrSemiconductor TotalRecombination'，
        那麼就可以把不同溫度但有著此部分名稱的數據都放在一起。
        不過如果是針對實驗，那就不需要這麼做。只需要讀取不同溫度下的 IV 即可。
        :param rate_name: TCAD請輸入如'IntegrSemiconductor TotalRecombination'
                          但實驗數據則可空白，或是任何其他名稱。因為此時不再藉此名稱判斷是否存入數據。
        :param source: TCAD請輸入 'Simulation'，實驗數據則輸入 'Experiment'。
        :return: 若是 TCAD，則回傳屬於該再生機制的 rate 列表（記憶體位址）。
                 若是 Experiment，則回傳該實驗數據的不同溫度的 IV 列表（記憶體位址）。
        """
        if source == 'Simulation':
            if rate_name not in self.__Rate:
                self.__Rate[rate_name] = \
                    [Data.IV(self.__Location, curve, curve, self.__Area) for curve in self.__Curve if
                     rate_name in curve]
            return self.__Rate[rate_name]
        elif source == 'Experiment':
            if rate_name not in self.__Rate:
                self.__Rate[rate_name] = \
                    [Data.IV(self.__Location, 'V', curve, self.__Area) for curve in self.__Curve if
                     'V' not in curve]
            return self.__Rate[rate_name]
        else:
            raise BaseException("Wrong source: %s" % source)

    def get_I_of_T(self, voltage, rate_name, Tmin, Tmax, source):
        """
        :param voltage: [V]
        :param rate_name: 'IntegrSemiconductor TotalRecombination'
        :param Tmin: [K]
        :param Tmax: [K]
        :return: 輸出一個 list，依溫度順序排列且大於零的 I_gr (V) [A]
        """
        if rate_name not in self.__Rate:
            self.get_recombination_rate(rate_name, source)
        if rate_name not in self.__IofT:
            self.__IofT[rate_name] = dict()
            if source == 'Simulation':
                self.__IofT[rate_name][voltage] = \
                    [q * abs(rate.find_I(voltage)) for T, rate in enumerate(self.__Rate[rate_name])
                     if Tmax >= self.__TemperatureList[T] >= Tmin]
            elif source == 'Experiment':
                self.__IofT[rate_name][voltage] = \
                    [abs(rate.find_I(voltage)) for T, rate in enumerate(self.__Rate[rate_name])
                     if Tmax >= self.__TemperatureList[T] >= Tmin]
            else:
                raise BaseException('Wrong source: %s' % source)
        elif voltage not in self.__IofT[rate_name]:
            if source == 'Simulation':
                self.__IofT[rate_name][voltage] = \
                    [q * abs(rate.find_I(voltage)) for T, rate in enumerate(self.__Rate[rate_name])
                     if Tmax >= self.__TemperatureList[T] >= Tmin]
            elif source == 'Experiment':
                self.__IofT[rate_name][voltage] = \
                    [abs(rate.find_I(voltage)) for T, rate in enumerate(self.__Rate[rate_name])
                     if Tmax >= self.__TemperatureList[T] >= Tmin]
            else:
                raise BaseException('Wrong source: %s' % source)
        return np.asarray(self.__IofT[rate_name][voltage])

    def get_trap_level(self, type, I_of_T, voltage, Tmin, Tmax):
        """
        輸出在某電壓之變溫電流得到的 trap level
        :param type: 1, 2, 3（不同 fitting 方式，不同物理模型）
        :param I_of_T:
        :param voltage:
        :return: Eti = Et - Ei（應該是有絕對值，需要再從理論檢驗）
        """
        if voltage in self.__Eti:
            if type in self.__Eti[voltage]:
                return self.__Eti[voltage][type]
            else:
                self.__Eti[voltage][type] = self.__findeti(type, I_of_T, voltage, Tmin, Tmax)
                return self.__Eti[voltage][type]
        else:
            self.__Eti[voltage] = dict()
            self.__Eti[voltage][type] = self.__findeti(type, I_of_T, voltage, Tmin, Tmax)
            return self.__Eti[voltage][type]

    def __findeti(self, type, I_of_T, voltage, Tmin, Tmax):
        """
        這裡會有個問題，就是底下的 0 > voltage > self.__Vpt
        如果現在的 V 都是正的，那麼這條件似乎就永遠不會達成。
        :param type:
        :param I_of_T:
        :param voltage:
        :param Tmin:
        :param Tmax:
        :return:
        """
        inverseT = [1/T for T in self.__TemperatureList if Tmin <= T <= Tmax]

        def Eg_InP(T):
            alpha = 1.420
            beta = 0.000398
            gamma = 189
            return alpha - (beta * (T ** 2))/(T + gamma)

        def Eg_InGaAs(T):
            alpha = 0.8186754144533469
            beta = 0.000494420835956525
            gamma = 286.3332895176157
            return alpha - (beta * (T ** 2))/(T + gamma)

        def Eti(slope, V):
            if V < self.__Vpt:
                return -slope / np.log10(np.e) * kB / q - 0.5 * Eg300_InGaAs
            else:
                return -slope / np.log10(np.e) * kB / q - 0.5 * Eg300_InP

        if type == '1':
            def line(x, a, b):
                return a * x + b

            popt, pcov = curve_fit(line, inverseT, np.log10(I_of_T))
            return Eti(tuple(popt)[0], voltage)
        elif type == '2':
            def line(x, a, b):
                return a * x - 1.5 * np.log10(x) + b

            popt, pcov = curve_fit(line, inverseT, np.log10(I_of_T))
            return Eti(tuple(popt)[0], voltage)
        elif type == '3':
            if 0 > voltage > self.__Vpt:
                def line3_InP(x, a, b):
                    return a * x - 1.5 * np.log10(x) + b - np.log10(np.e) / (2 * kB) * q * Eg_InP(1/x) * x

                popt, pcov = curve_fit(line3_InP, inverseT, np.log10(I_of_T))
            else:
                def line3_InGaAs(x, a, b):
                    return a * x - 1.5 * np.log10(x) + b - np.log10(np.e) / (2 * kB) * q * Eg_InGaAs(1/x) * x

                popt, pcov = curve_fit(line3_InGaAs, inverseT, np.log10(I_of_T))
            return - tuple(popt)[0] / np.log10(np.e) * kB / q
        elif type == 'Hurkx':
            if abs(voltage) < abs(self.__Vpt):
                def line_TAT_InP(T, a, b, Eti):
                    return a * T + b / (T ** 2) - Eti * q / kB - q * Eg_InP(T) / (2 * kB)

                T_eff = np.asarray([1 / T for T in inverseT])
                guess = [5, 1e8, 0.33]
                popt, pcov = curve_fit(line_TAT_InP, T_eff, T_eff * np.log(I_of_T), p0=guess)
                print('V(InP): %s,   a: %.2f,   b: %.3e,   Eti: %.3f' % (voltage, tuple(popt)[0], tuple(popt)[1], tuple(popt)[2]))
                return tuple(popt)[2]
            else:
                def line_TAT_InGaAs(T, a, b, Eti):
                    return a * T + b / (T ** 2) - Eti * q / kB - q * Eg_InGaAs(T) / (2 * kB)

                T_eff = np.asarray([1 / T for T in inverseT])
                guess = [7, 1e7, 0.3]
                popt, pcov = curve_fit(line_TAT_InGaAs, T_eff, T_eff * np.log(I_of_T), p0=guess)
                print('V(InGaAs): %s,   a: %.2f,   b: %.3e,   Eti: %.3f' % (voltage, tuple(popt)[0], tuple(popt)[1], tuple(popt)[2]))
                return tuple(popt)[2]
        else:
            raise BaseException("Wrong type: %s" % type)


class QuantityDistribution(object):
    def __init__(self, Location):
        self.__Location = Location
        self.__Curve = utils.parse_csv(Location, 2, 0, 0, 'Yes')
        self.__DataList = [self.__Curve[i].rsplit("(", 1)[0] for i in range(len(self.__Curve)-1)
                            if i == 0 or (self.__Curve[i-1].split("(")[0] != self.__Curve[i].split("(")[0])]
        self.__Rate = dict()
        self.__Qb = dict()
        self.__Emax = dict()

    def get_DataList(self):
        for i, element in enumerate(self.__DataList):
            print('%2s, %25s, %25s, %50s' % (i, element.split(" ")[0], element.split(" ")[1], element))
        return self.__DataList

    def get_recombination_rate(self, rate_name):
        """
        得到一個 list，其內容為屬於 Data.CSV 類別的物件。
        :param rate_name: 例如 'IntegrMultiplication Band2BandGeneration'
        :return: RateList = [Rate(Thickness1), Rate(Thickness2), ... ]
        """
        if rate_name not in self.__Rate:
            self.__Rate[rate_name] = [Data.CSV(self.__Location, curve, curve)
                                      for curve in self.__Curve if rate_name in curve]
        return self.__Rate[rate_name]

    def get_QuantityNearBreakdown(self, num, breakdown_voltage, rate_name, factor=0.9):
        if rate_name not in self.__Rate:
            self.get_recombination_rate(rate_name)
        if rate_name not in self.__Qb:
            self.__Qb[rate_name] = dict()
            self.__Qb[rate_name][num] = utils.find(self.__Rate[rate_name][num].X, abs(self.__Rate[rate_name][num].Y),
                                                   - factor * abs(breakdown_voltage), 'log')
        elif num not in self.__Qb[rate_name]:
            self.__Qb[rate_name][num] = utils.find(self.__Rate[rate_name][num].X, abs(self.__Rate[rate_name][num].Y),
                                                   - factor * abs(breakdown_voltage), 'log')
        return self.__Qb[rate_name][num]