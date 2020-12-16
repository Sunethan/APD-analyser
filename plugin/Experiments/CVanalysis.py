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
from scipy.optimize import curve_fit
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (20, 9.3),
          'axes.labelsize': 'x-large',
          'axes.titlesize':'x-large',
          'xtick.labelsize':'x-large',
          'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)
plt.rcParams.update({'font.size': 9})

ColorSet = ['dodgerblue', 'yellowgreen', 'goldenrod', 'darkviolet', 'b', 'brown', 'darkorange',
            'hotpink', 'fuchsia', 'g', 'royalblue', 'tomato', 'purple', 'olive', 'darkgreen']
LineStyleSet = ['-', '-.', '--', ':']

class CV(object):
    def __init__(self, cv_data_location, area):
        self._area = {device: area[i] for i, device in enumerate(cv_data_location.keys())}
        self._cv_data = {device: Data.Read(location, 'Standard', 'VHi', 'CpHi', show='No')
                         for device, location in cv_data_location.items()}
        self._doping_profile = None

    def save_doping(self, doping_profile_location, xlabel, ylabel):
        self._doping_profile = Data.Read(doping_profile_location, 'Standard', xlabel, ylabel, show='No')

    def doping_plot(self, fig_number, method, title, device, Xj):
        if self._doping_profile is None:
            raise BaseException("Please use CV.save_doping first")
        if not isinstance(method, list):
            method = [method]
        if not isinstance(device, list):
            device = [device]
        for element in method:
            if element not in ['direct', 'cylindrical', 'conical']:
                raise BaseException("Wrong method (direct/cylindrical/conical): %s" % method)
        Result = {'direct': {d: self.direct(self._cv_data[d].Y / self._area[d], self._cv_data[d].X)
                             for d in self._device_list('direct', device)},
                  'cylindrical': {list(list(d[0].items())[0])[0] + ' & ' + list(list(d[1].items())[0])[0]: DualCV(
                      d, self._area).cylindrical() for d in self._device_list('cylindrical', device)}}

        plt.figure(fig_number)
        plt.plot(self._doping_profile.X, abs(self._doping_profile.Y), '-b', label='Nominal profile')
        for i, m in enumerate(method):
            for j, d in enumerate(self._device_list(m, device)):
                if m != 'direct':
                    d = list(list(d[0].items())[0])[0] + ' & ' + list(list(d[1].items())[0])[0]
                plt.plot(Result[m][d]['Position'] + Xj, Result[m][d]['Concentration'],
                         color=ColorSet[j], linestyle=LineStyleSet[i - 1], label='[%s] %s' % (d, m))
        plt.yscale('log')
        plt.grid()
        plt.xlim((-4.5, 0))
        plt.xlabel(r'Position ($\mu$m)')
        plt.ylabel(r'Concentration (cm$^{-3}$)')
        plt.legend(loc='best')
        plt.title(title)

    def cv_plot(self, fig_number, title, *args):
        plt.figure(fig_number)
        for i, arg in enumerate(args):
            plt.plot(self._cv_data[arg].X, 1e12 * self._cv_data[arg].Y, color=ColorSet[i], label=arg)
        plt.grid()
        plt.xlabel('Voltage (V)')
        plt.ylabel('Capacitance (pF)')
        plt.legend(loc='best')
        plt.title(title)

    def _device_list(self, method, device):
        if method == 'direct':
            return device
        elif method == 'cylindrical' or method == 'conical':
            if len(set(device)) != len(device):
                raise BaseException("Check your device list. There must be double inputs.")
            pairs = []
            for i in range(len(device)):
                for j in range(i, len(device)):
                    pair = {device[i], device[j]}
                    if pair not in pairs and len(pair) == 2:
                        pairs.append(pair)
            pairs = [[{list(pair)[0]: self._cv_data[list(pair)[0]]},
                      {list(pair)[1]: self._cv_data[list(pair)[1]]}] for pair in pairs]
            return pairs
        else:
            raise BaseException("Wrong method (direct/cylindrical/conical): %s" % method)

    @staticmethod
    def direct(c, V):
        q = 1.6e-19
        eps = 13.9 * 8.85e-14
        position = np.array([eps / c_per_area for c_per_area in c]) * 1e4  # [um]
        dCdV = utils.dydx(V, c)  # [F/V/cm2]
        concentration = np.array(
            [abs(c_per_area ** 3 / (q * eps) / dCdV[i]) for i, c_per_area in enumerate(c)])  # [cm-3]
        return {'Position': position, 'Concentration': concentration}


class DualCV(object):
    def __init__(self, device_pair, area):
        self._data = {'V': {}, 'C': {}, 'r': {}}
        for i in [1, 2]:
            code, data = list(device_pair[i - 1].items())[0]
            self._data['V'][i] = data.X  # [V]
            self._data['C'][i] = data.Y  # [F]
            self._data['r'][i] = np.sqrt(area[code] / np.pi)  # [cm]
        if set(self._data['V'][1]) != set(self._data['V'][2]):
            self._data['V']['combination'] = list(set(list(self._data['V'][1]) + list(self._data['V'][2])))
            self._data['V']['combination'].sort()
            self._data['V']['combination'] = np.array(self._data['V']['combination'])
            for i in [1, 2]:
                self._data['C'][i] = utils.interpolation(self._data['V'][i], self._data['V']['combination'],
                                                         self._data['C'][i])
        else:
            self._data['V']['combination'] = self._data['V'][1]

    def cylindrical(self):
        r1 = self._data['r'][1]
        r2 = self._data['r'][2]
        c1 = self._data['C'][1]
        c2 = self._data['C'][2]

        cp = (r2 * c1 - r1 * c2) / (np.pi * r1 * r2 * (r1 - r2))
        # cl = (r1 ** 2 * c2 - r2 ** 2 * c1) / (2 * np.pi * r1 * r2 * (r1 - r2))
        return CV.direct(abs(cp), self._data['V']['combination'])

    def conical(self):
        pass