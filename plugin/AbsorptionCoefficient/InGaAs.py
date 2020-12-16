import utils
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import matplotlib.pylab as pylab
import DataAnalysis as Data
import inspect

params = {'legend.fontsize': 'x-large',
          'figure.figsize': (7.5, 7),
          'axes.labelsize': 'x-large',
          'axes.titlesize':'x-large',
          'xtick.labelsize':'x-large',
          'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)
plt.rcParams.update({'font.size': 13})

file_location = '/Users/Ethan/PycharmProjects/Python/Avalanche Photodiode/plugin/AbsorptionCoefficient/data/'


def Bacher(energy, color='b'):
    location = file_location + 'InGaAs_F.R.Bacher.csv'
    raw = Data.CSV(data=location, xlabel='energy', ylabel='absorption')
    if energy is not None:
        if energy < raw.X[0] or energy > raw.X[-1]:
            raise BaseException("Energy should be in [%.2f, %.2f]" % (raw.X[0], raw.X[-1]))
        return utils.find(raw.X, raw.Y, energy, 'log')
    else:
        plt.plot(raw.X, raw.Y, linestyle='-.', color=color, marker='s', label='F.P. Bacher (1988)')


def Adachi(energy, color='b'):
    location = file_location + 'InGaAs_Adachi.csv'
    raw = Data.CSV(data=location, xlabel='energy', ylabel='absorption')
    if energy is not None:
        if energy < raw.X[0] or energy > raw.X[-1]:
            raise BaseException("Energy should be in [%.2f, %.2f]" % (raw.X[0], raw.X[-1]))
        return utils.find(raw.X, raw.Y, energy, 'log')
    else:
        plt.plot(raw.X, raw.Y, linestyle='-.', color=color, marker='s', label='S. Adachi (1989)')


def Zielinski(energy, color='b'):
    location = file_location + 'InGaAs_E.Zielinski.csv'
    raw = Data.CSV(data=location, xlabel='energy', ylabel='absorption')
    if energy is not None:
        if energy < raw.X[0] or energy > raw.X[-1]:
            raise BaseException("Energy should be in [%.2f, %.2f]" % (raw.X[0], raw.X[-1]))
        return utils.find(raw.X, raw.Y, energy, 'log')
    else:
        plt.plot(raw.X, raw.Y, linestyle='-.', color=color, marker='s', label='E. Zielinski (1986)')