from Experiment import TemperatureDependentIV as TempIV
import numpy as np
import physics as phys
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import matplotlib.pylab as pylab
import utils
import matplotlib.patches as patches
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (20, 9.3),
          'axes.labelsize': 'x-large',
          'axes.titlesize':'x-large',
          'xtick.labelsize':'x-large',
          'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)
plt.rcParams.update({'font.size': 9})

q = 1.6e-19  # [C]
kB = 1.38e-23  # [J/K]
ColorList = ['maroon', 'r', 'crimson', 'orangered', 'darkorange', 'khaki', 'goldenrod', 'sienna',
             'olive', 'tan', 'salmon', 'silver', 'darkgreen', 'green', 'lime', 'lightblue', 'b',
             'darkblue', 'navy', 'cyan', 'turquoise', 'pink', 'm', 'plum', 'violet', 'indigo']


class TemperatureDependentIV(object):
    def __init__(self, Vpt, A, EtValue, EtY, Temperature, SRH_Field, Tmin, Tmax, Voltage):
        """
        CSV 必須是有著相同的物理模型，而只有溫度不同的元件電性模擬數據。基本上可分為有沒有開 trap-assisted tunneling。
        :param Vpt: 擊穿電壓，記得必須用負數。
        :param A: 元件面積 (cm2)
        :param EtValue: trap level，如 EtValue = {'InP': 0.336, 'InGaAs': 0.2}
        :param EtY: 畫 Eti vs. voltage 時的 y 軸上下限，如 EtY = [-0.5, 3]
        :param Temperature: 模擬的溫度範圍，如 Temperature = np.arange(450, 125, -25) (K)
        :param SRH_Field: ['On', 'Off']
        :param Tmin: 太低溫的模擬，收斂電壓可能很小，所以就不適合用來求 Et。必須用 dict 輸入，例如：
               Tmin = {'On': 225, 'Off': 275}
        :param Tmax: 同上
        :param Voltage: 想用來求 Eti 的電壓範圍，如 [-1, -3, -5, -10, -15, -17.5, -20, -23, -24, -25]
        """
        self.__Vpt = Vpt
        self.__A = A
        self.__EtValue = EtValue
        self.__EtY = EtY
        self.__Temperature = Temperature
        self.__SRH_Field = SRH_Field
        self.__Tmin = Tmin
        self.__Tmax = Tmax
        self.__Voltage = Voltage
        self.__ColorList = ColorList

        self.__figure = {i+1: plt.figure(i+1) for i in range(4)}
        self.__ax = {i+1: [self.__figure[i+1].add_subplot(1, 2, 1), self.__figure[i+1].add_subplot(1, 2, 2)] for i in range(4)}

    def addplot(self, num, key, x, y, colorset='k', name='AddPlot'):
        if key == 'On':
            self.__ax[num][0].plot(x, y, color=colorset, linestyle='-.', label='%s' % name)
        elif key == 'Off':
            self.__ax[num][1].plot(x, y, color=colorset, linestyle='-.', label='%s' % name)
        else:
            raise BaseException('Wrong key: %s' % key)

    def plot(self, source='Simulation'):
        """
        可以自動畫圖。同時他會讀取存放為 'data/Field_On.csv' 的檔案。
        :return:
        """
        Field = {key: TempIV('data/Field_' + key + '.csv', self.__Temperature, self.__Vpt, self.__A, source) for key in self.__SRH_Field}
        Te = {key: np.asarray([T for T in self.__Temperature if self.__Tmin[key] <= T <= self.__Tmax[key]]) for key in self.__SRH_Field}
        Mechanism = {key: Field[key].get_mechanism() for key in self.__SRH_Field}
        TotalGR = {key: Field[key].get_recombination_rate('IntegrSemiconductor TotalRecombination', source) for key in
                   self.__SRH_Field}
        SRH = {key: Field[key].get_recombination_rate('IntegrSemiconductor SRHRecombination', source) for key in self.__SRH_Field}
        IofT = {
            key: {V: Field[key].get_I_of_T(V, 'IntegrSemiconductor TotalRecombination', self.__Tmin[key], self.__Tmax[key], source) for V in
                  self.__Voltage} for key in self.__SRH_Field}
        Eti1 = {key: {V: Field[key].get_trap_level('1', IofT[key][V], V, self.__Tmin[key], self.__Tmax[key]) for V in self.__Voltage} for key
                in self.__SRH_Field}
        Eti2 = {key: {V: Field[key].get_trap_level('2', IofT[key][V], V, self.__Tmin[key], self.__Tmax[key]) for V in self.__Voltage} for key
                in self.__SRH_Field}
        Eti3 = {key: {V: Field[key].get_trap_level('3', IofT[key][V], V, self.__Tmin[key], self.__Tmax[key]) for V in self.__Voltage} for key
                in self.__SRH_Field}
        EtiHurkx = {key: {V: Field[key].get_trap_level('Hurkx', IofT[key][V], V, self.__Tmin[key], self.__Tmax[key]) for V in self.__Voltage}
                    for key in self.__SRH_Field}

        # 數據圖標題
        FigureTitle = r'E$_{ti}$(InP): %s; E$_{ti}$(InGaAs): %s' % (self.__EtValue['InP'], self.__EtValue['InGaAs'])

        for i, key in enumerate(self.__SRH_Field):
            for index, T in enumerate(self.__Temperature):
                if self.__Tmin[key] <= T <= self.__Tmax[key]:
                    self.__ax[1][i].plot(TotalGR[key][index].X, q * abs(TotalGR[key][index].Y), color=self.__ColorList[index],
                             label=T)
                    self.__ax[1][i].plot(SRH[key][index].X, q * abs(SRH[key][index].Y), color=self.__ColorList[index],
                             linestyle='-.')
            self.__ax[1][i].set_yscale('log')
            self.__ax[1][i].grid(True)
            self.__ax[1][i].legend(loc='best')
            self.__ax[1][i].set_xlabel('Voltage (V)')
            self.__ax[1][i].set_ylabel('Current (A)')
            self.__ax[1][i].set_title(r'SRH field enhancement: %s' % key)
            self.__ax[1][i].set_ylim((1e-22, 1e-2))
            self.__ax[1][i].set_xlim((-120, 0))
            for V in self.__Voltage:
                self.__ax[1][i].vlines(x=V, ymin=1e-22, ymax=1e-5, color='k', linestyles='-.')
        self.__figure[1].suptitle(FigureTitle)

        for i, key in enumerate(self.__SRH_Field):
            for index, V in enumerate(self.__Voltage):
                self.__ax[2][i].plot([1 / T for T in self.__Temperature if self.__Tmin[key] <= T <= self.__Tmax[key]],
                                     IofT[key][V], color=self.__ColorList[index], marker='o', markersize='9')
            self.__ax[2][i].grid(True)
            self.__ax[2][i].set_yscale('log')
            self.__ax[2][i].set_xlabel(r'1/T (K$^{-1}$)')
            self.__ax[2][i].set_ylabel('Current (A)')
            self.__ax[2][i].set_title('SRH field enhancement: %s\nlog(I) vs. T' % key)
            self.__ax[2][i].ticklabel_format(axis='x', scilimits=(0, 0), style='sci')
        self.__figure[2].suptitle(FigureTitle)

        for i, key in enumerate(self.__SRH_Field):
            for index, V in enumerate(self.__Voltage):
                self.__ax[3][i].plot(Te[key], Te[key] * np.log(IofT[key][V]),
                                     color=self.__ColorList[index], marker='o', markersize='9')
            self.__ax[3][i].grid(True)
            self.__ax[3][i].set_xlabel(r'T (K)')
            self.__ax[3][i].set_ylabel(r'T$\ln(IT)$')
            self.__ax[3][i].ticklabel_format(axis='y', scilimits=(0, 0), style='sci')
            self.__ax[3][i].set_title('SRH field enhancement: %s\nTln(I) vs. T' % key)
            self.__ax[3][i].legend(loc='best')
        self.__figure[3].suptitle(FigureTitle)

        for i, key in enumerate(self.__SRH_Field):
            self.__ax[4][i].plot(self.__Voltage, [Eti1[key][V] for V in self.__Voltage], '-bo', markersize='10', label='Line 1')
            self.__ax[4][i].plot(self.__Voltage, [Eti2[key][V] for V in self.__Voltage], '-ro', markersize='10', label='Line 2')
            self.__ax[4][i].plot(self.__Voltage, [Eti3[key][V] for V in self.__Voltage], '-go', markersize='10', label='Line 3')
            self.__ax[4][i].plot(self.__Voltage, [EtiHurkx[key][V] for V in self.__Voltage],
                                 '-mo', markersize='10', label='Line Hurkx')
            self.__ax[4][i].hlines(y=self.__EtValue['InP'], xmin=self.__Vpt, xmax=0, linestyles='-.', color='k')
            self.__ax[4][i].hlines(y=self.__EtValue['InGaAs'], xmin=self.__Voltage[-1], xmax=self.__Vpt, linestyles='-.', color='k')
            self.__ax[4][i].grid(True)
            self.__ax[4][i].set_ylim((self.__EtY[0], self.__EtY[1]))
            self.__ax[4][i].legend(loc='best')
            self.__ax[4][i].set_xlabel('Voltage (V)')
            self.__ax[4][i].set_ylabel(r'E$_{ti}$ (eV)')
            self.__ax[4][i].set_title('SRH field enhancement: %s\nTrap level distribution' % key)
        self.__figure[4].suptitle(FigureTitle)

        while True:
            try:
                plt.show()
            except UnicodeDecodeError:
                continue
            break


class DictData(object):
    def __init__(self, data_type, **kwargs):
        self.__dict = self.set_dict(dict(), data_type, kwargs)

    def get_dict(self):
        return self.__dict

    def set_dict(self, data, data_type, key_dict):
        if len(key_dict) == 1:
            for key, values in key_dict.items():
                for value in values:
                    if data_type == 'dict':
                        data[value] = dict()
                    elif data_type == 'list':
                        data[value] = []
                    else:
                        raise BaseException("Wrong data type: %s" % data_type)
            return data
        else:
            dict_temp = key_dict.copy()
            dict_temp.pop(next(iter(dict_temp)), None)
            for value in key_dict[next(iter(key_dict))]:
                data[value] = dict()
                data[value] = self.set_dict(data[value], data_type, dict_temp)
            return data

class IVfitting(object):
    def __init__(self, data):
        self.doping, self.T_analysis, self.RawIV, self.count, self.ColorSet10, self.Material, \
        self.BTB, self.FitType, self.FitIV, self.LineSet2, self.f_Vpt, self.AbruptOrConstant, \
        self.Lifetime, self.Lifetime_p, self.Lifetime_n, self.Voltage_IT, self.FitIT, self.V1, \
        self.V2, self.T_analysis_IT, self.V_InP, self.V_InGaAs, self.weight, self.FittingSwitch = data

    def plot_fig1(self, figtitle):
        fig_IV = plt.figure(1)
        plt.subplot(1, 2, 1)
        for i, T in enumerate(self.T_analysis):
            raw_V = [element for i, element in enumerate(self.RawIV[T].X) if i % self.count == 0]
            raw_I = [element for i, element in enumerate(self.RawIV[T].Y) if i % self.count == 0]
            plt.plot(raw_V, raw_I, linestyle='none', marker='o', fillstyle='none', color=self.ColorSet10[i], label=T)
            for material in self.Material:
                plt.plot(self.BTB[material][T][0], self.BTB[material][T][1], linestyle=':', color=self.ColorSet10[i])
            for j, model in enumerate(self.FitType):
                for material in self.Material:
                    plt.plot(self.FitIV[material][model][T][0], self.FitIV[material][model][T][1], color=self.ColorSet10[i],
                             linestyle=self.LineSet2[j])
        plt.vlines(x=-self.f_Vpt[self.doping], ymin=1e-13, ymax=1e-6, color='k', linestyle='-.', label=r'V$_{pt}$')
        plt.grid()
        plt.yscale('log')
        plt.xlabel('Voltage (V)')
        plt.ylabel('Current (A)')
        plt.title(r'IV modeling @ 300K @ d$_{InGaAsP,charge}=0$')
        plt.xlim((-45, 1))
        plt.ylim((1e-13, 1e-6))
        plt.legend(loc='best')
        plt.subplot(1, 2, 2)
        for i, material in enumerate(self.Material):
            plt.plot(*self.AbruptOrConstant.PlotEm(300, material), color=self.ColorSet10[i], label=material)
        plt.vlines(x=-self.f_Vpt[self.doping], ymin=0, ymax=7e5, color='k',
                   linestyle='-.', label=r'V$_{pt}=%s$' % -self.f_Vpt[self.doping])
        plt.grid()
        plt.xlabel('Voltage (V)')
        plt.ylabel('Electric field (V/cm)')
        plt.title('Maximum electric field distribution')
        plt.legend(loc='best')
        plt.ylim((-0.3e5, 8e5))
        plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
        fig_IV.suptitle(figtitle, fontsize=18)

    def plot_fig2(self, figtitle):
        fig_Eti_TAT = plt.figure(2)
        plt.subplot(1, 2, 1)
        for material in self.Material:
            for model in self.FitType:
                plt.plot(self.T_analysis, [self.FitIV[material][model][T][3] for T in self.T_analysis], marker='o', fillstyle='none',
                         label=r'[%s][%s]' % (material, model))
        plt.ylabel(r'E$_t$-E$_i$ (eV)')
        plt.legend(loc='best')
        plt.ylim((-0.1, 0.18))
        plt.grid()

        plt.subplot(1, 2, 2)
        for i, T in enumerate(self.T_analysis):
            for j, material in enumerate(self.Material):
                if material == 'InP':
                    plt.plot(*self.AbruptOrConstant.d_Gamma(T, material), color=self.ColorSet10[i], linestyle=self.LineSet2[j], label='%sK' % T)
                    plt.plot(*self.AbruptOrConstant.Plotdm(T, material), color=self.ColorSet10[i], linestyle=self.LineSet2[j])
                else:
                    plt.plot(*self.AbruptOrConstant.d_Gamma(T, material), color=self.ColorSet10[i], linestyle=self.LineSet2[j])
                    plt.plot(*self.AbruptOrConstant.Plotdm(T, material), color=self.ColorSet10[i], linestyle=self.LineSet2[j])
        plt.grid()
        plt.xlabel('Voltage (V)')
        plt.ylabel(r'$d_{TAT}$ (cm)')
        plt.title(r'$d_{TAT}\equiv\int\Gamma(x)dx$')
        plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
        plt.legend(loc='best', ncol=3)
        fig_Eti_TAT.suptitle(figtitle, fontsize=18)

    def plot_fig3(self, figtitle):
        fig_Lifetime = plt.figure(3)
        for i, material in enumerate(self.Material):
            plt.subplot(2, 2, i + 1)
            for model in self.FitType:
                plt.plot(self.T_analysis, [1e9 * self.Lifetime[material][model][T][0] for T in self.T_analysis], marker='o',
                         fillstyle='none',
                         linestyle='-', label=r'$\tau_p$[%s]' % model)
            plt.hlines(y=self.Lifetime_p[material], xmin=240, xmax=330, colors='k', linestyles='-.', label=material)
            plt.grid(True)
            plt.legend(loc='best')
            plt.title(material)
            plt.ylabel(r'$\tau_p$ (ns)')

            plt.subplot(2, 2, i + 3)
            for model in self.FitType:
                plt.plot(self.T_analysis, [1e9 * self.Lifetime[material][model][T][1] for T in self.T_analysis], marker='o',
                         fillstyle='none',
                         linestyle='-', label=r'$\tau_n$[%s]' % model)
            plt.hlines(y=self.Lifetime_n[material], xmin=240, xmax=330, colors='k', linestyles='-.', label=material)
            plt.grid(True)
            plt.legend(loc='best')
            plt.ylabel(r'$\tau_n$ (ns)')
            plt.xlabel('Temperature (K)')
        fig_Lifetime.suptitle(figtitle, fontsize=18)

    def plot_fig4(self, figtitle):
        ColorModel = {'SRH': 'r', 'TAT': 'b'}
        fig_IT = plt.figure(4)
        plt.subplot(1, 2, 1)
        for material in self.Material:
            if material == 'InP':
                for model in self.FitType:
                    plt.plot(self.Voltage_IT[material],
                             [self.FitIT[material][model][V][2] for V in self.Voltage_IT[material]],
                             color=ColorModel[model], label=model)
            else:
                for model in self.FitType:
                    plt.plot(self.Voltage_IT[material],
                             [self.FitIT[material][model][V][2] for V in self.Voltage_IT[material]],
                             color=ColorModel[model])
        plt.grid()
        plt.xlabel('Voltage (V)')
        plt.ylabel(r'E$_t$-E$_t$ (eV)')
        plt.title('Trap level distribution')
        plt.legend(loc='best')

        plt.subplot(1, 2, 2)
        IT_InP = np.asarray([utils.find(self.RawIV[T].X, self.RawIV[T].Y, self.V1, 'log') for T in self.T_analysis_IT])
        IT_InGaAs = np.asarray([utils.find(self.RawIV[T].X, self.RawIV[T].Y, self.V2, 'log') for T in self.T_analysis_IT])
        plt.plot(1 / self.T_analysis_IT, np.log10(IT_InP), linestyle='none', marker='o', markersize='15', fillstyle='none', label=r'$%s$ (V)' % self.V1)
        plt.plot(1 / self.T_analysis_IT, np.log10(IT_InGaAs), linestyle='none', marker='o', markersize='15', fillstyle='none', label=r'$%s$ (V)' % self.V2)
        for material in self.Material:
            for model in self.FitType:
                if material == 'InP':
                    plt.plot(1 / self.FitIT[material][model][self.V1][0],
                             np.log10(self.FitIT[material][model][self.V1][1]),
                             color=ColorModel[model], label=model)
                else:
                    plt.plot(1 / self.FitIT[material][model][self.V2][0],
                             np.log10(self.FitIT[material][model][self.V2][1]),
                             color=ColorModel[model])
        plt.grid()
        plt.xlabel(r'1/T (K$^{-1}$)')
        plt.ylabel('log10(I)')
        plt.legend(loc='best')
        plt.ticklabel_format(axis='x', style='sci', scilimits=(0, 0))
        fig_IT.suptitle(figtitle, fontsize=18)

    def plot_fig5(self, figtitle):
        tau_InP_p0 = 300  # [ns]
        tau_InP_n0 = 2.89  # [ns]
        tau_InGaAs_p0 = 8  # [ns]
        tau_InGaAs_n0 = 0.25  # [ns]
        ColorModel = {'SRH': 'r', 'TAT': 'b'}
        fig_IT_info = plt.figure(5)
        plt.subplot(2, 2, 1)
        for material in self.Material:
            for model in self.FitType:
                if material == 'InP':
                    plt.plot(self.V_InP, [tau_InP_p0 * self.FitIT[material][model][V][3][0] for V in self.V_InP],
                             color=ColorModel[model], label=model)
                else:
                    plt.plot(self.V_InGaAs,
                             [tau_InGaAs_p0 * self.FitIT[material][model][V][3][0] for V in self.V_InGaAs],
                             color=ColorModel[model])
        plt.grid()
        plt.yscale('log')
        plt.xlabel(r'Voltage (V)')
        plt.ylabel(r'$\tau_{p, 300}$ (ns)')
        plt.legend(loc='best')

        plt.subplot(2, 2, 2)
        for material in self.Material:
            for model in self.FitType:
                if material == 'InP':
                    plt.plot(self.V_InP, [self.FitIT[material][model][V][3][2] for V in self.V_InP],
                             color=ColorModel[model], label=model)
                else:
                    plt.plot(self.V_InGaAs, [self.FitIT[material][model][V][3][2] for V in self.V_InGaAs],
                             color=ColorModel[model])
        plt.grid()
        plt.yscale('log')
        plt.xlabel(r'Voltage (V)')
        plt.ylabel(r'$\alpha_p$ (-)')
        plt.legend(loc='best')

        plt.subplot(2, 2, 3)
        for material in self.Material:
            for model in self.FitType:
                if material == 'InP':
                    plt.plot(self.V_InP, [tau_InP_n0 * self.FitIT[material][model][V][3][1] for V in self.V_InP],
                             color=ColorModel[model], label=model)
                else:
                    plt.plot(self.V_InGaAs,
                             [tau_InGaAs_n0 * self.FitIT[material][model][V][3][1] for V in self.V_InGaAs],
                             color=ColorModel[model])
        plt.grid()
        plt.yscale('log')
        plt.xlabel(r'Voltage (V)')
        plt.ylabel(r'$\tau_{n,300}$ (ns)')
        plt.legend(loc='best')

        plt.subplot(2, 2, 4)
        for material in self.Material:
            for model in self.FitType:
                if material == 'InP':
                    plt.plot(self.V_InP, [self.FitIT[material][model][V][3][3] for V in self.V_InP],
                             color=ColorModel[model], label=model)
                else:
                    plt.plot(self.V_InGaAs, [self.FitIT[material][model][V][3][3] for V in self.V_InGaAs],
                             color=ColorModel[model])
        plt.grid()
        plt.yscale('log')
        plt.xlabel(r'Voltage (V)')
        plt.ylabel(r'$\alpha_n$ (-)')
        plt.legend(loc='best')
        fig_IT_info.suptitle(figtitle, fontsize=18)

    def plot_fig6(self, figtitle):
        def lifetime_p(T, tp, alpha_p, material):
            if material == 'InP':
                tau_p0 = 300e-9  # [s]
            elif material == 'InGaAs':
                tau_p0 = 8e-9  # [s]
            else:
                raise BaseException("Wrong material %s" % material)
            tau_p = tp * tau_p0 * (T / 300) ** alpha_p
            return tau_p

        def lifetime_n(T, tn, alpha_n, material):
            if material == 'InP':
                tau_n0 = 2.89e-9  # [s]
            elif material == 'InGaAs':
                tau_n0 = 0.25e-9  # [s]
            else:
                raise BaseException("Wrong material %s" % material)
            tau_n = tn * tau_n0 * (T / 300) ** alpha_n
            return tau_n

        def lifetime_approx(lifetime_array_at_T, weight_array):
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
        fig_lifetime = plt.figure(6)
        for i, material in enumerate(self.Material):
            plt.subplot(2, 2, i + 1)
            for model in self.FitType:
                LabelMarker1 = 0
                for V in self.Voltage_IT[material]:
                    tp = self.FitIT[material][model][V][3][0]
                    alpha_p = self.FitIT[material][model][V][3][2]
                    self.Lifetime_p_DB[material][model][V] = lifetime_p(self.T_analysis_IT, tp, alpha_p, material)  # [s]
                    '''
                    if LabelMarker1 == 0:
                        plt.plot(self.T_analysis_IT, 1e9 * Lifetime_p_DB[material][model][V],
                                 color=ColorModel[model], label=model)
                        LabelMarker1 = 1
                    else:
                        plt.plot(self.T_analysis_IT, 1e9 * Lifetime_p_DB[material][model][V],
                                 color=ColorModel[model])
                    '''
                Lifetime_avg_p = np.asarray([lifetime_approx([self.Lifetime_p_DB[material][model][V][i] for V
                                                              in self.Voltage_IT[material]], self.weight[material])
                                             for i in range(len(self.T_analysis_IT))])
                plt.plot(self.T_analysis_IT, 1e9 * Lifetime_avg_p, color=ColorModel[model], marker='o',
                         fillstyle='none', label=r'[IT][%s]' % model)
                if 'IV' in self.FittingSwitch:
                    plt.plot(self.T_analysis, [1e9 * self.Lifetime[material][model][T][0] for T in self.T_analysis],
                             marker='o', color=ColorModel[model], fillstyle='none', linestyle='--',
                             label=r'[IV][%s]' % model)
            plt.hlines(y=self.Lifetime_p[material], xmin=240, xmax=330, colors='k', linestyles='-.', label=material)
            plt.grid()
            plt.ylabel(r'$\tau_p$ (ns)')
            plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
            plt.title('%s lifetime' % material)
            plt.legend(loc='best')

            plt.subplot(2, 2, i + 3)
            for model in self.FitType:
                LabelMarker2 = 0
                for V in self.Voltage_IT[material]:
                    tn = self.FitIT[material][model][V][3][1]
                    alpha_n = self.FitIT[material][model][V][3][3]
                    self.Lifetime_n_DB[material][model][V] = lifetime_n(self.T_analysis_IT, tn, alpha_n, material)  # [s]
                    '''
                    if LabelMarker2 == 0:
                        plt.plot(self.T_analysis_IT, 1e9 * Lifetime_n_DB[material][model][V],
                                 color=ColorModel[model], label=model)
                        LabelMarker2 = 1
                    else:
                        plt.plot(self.T_analysis_IT, 1e9 * Lifetime_n_DB[material][model][V],
                                 color=ColorModel[model])
                    '''
                Lifetime_avg_n = np.asarray([lifetime_approx([self.Lifetime_n_DB[material][model][V][i] for V
                                                              in self.Voltage_IT[material]], self.weight[material])
                                             for i in range(len(self.T_analysis_IT))])
                plt.plot(self.T_analysis_IT, 1e9 * Lifetime_avg_n, color=ColorModel[model], marker='o',
                         fillstyle='none', label=r'[IT][%s]' % model)
                if 'IT' in self.FittingSwitch:
                    plt.plot(self.T_analysis, [1e9 * self.Lifetime[material][model][T][1] for T in self.T_analysis],
                             color=ColorModel[model], marker='o', fillstyle='none', linestyle='-.',
                             label=r'[IV][%s]' % model)
            plt.hlines(y=self.Lifetime_n[material], xmin=240, xmax=330, colors='k', linestyles='-.', label=material)
            plt.grid()
            plt.xlabel('Temperature (K)')
            plt.ylabel(r'$\tau_n$ (ns)')
            plt.ticklabel_format(axis='y', style='sci', scilimits=(0, 0))
            plt.legend(loc='best')
        fig_lifetime.suptitle(figtitle, fontsize=18)