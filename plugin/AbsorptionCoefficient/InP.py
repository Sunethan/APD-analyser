import utils
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import matplotlib.pylab as pylab
import DataAnalysis as Data
import numpy as np
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


def Adachi(energy, color='b'):
    location = file_location + 'InP_Adachi.csv'
    raw = Data.CSV(data=location, xlabel='energy', ylabel='absorption')
    if energy is not None:
        if energy < raw.X[0] or energy > raw.X[-1]:
            raise BaseException("Energy should be in [%.2f, %.2f]" % (raw.X[0], raw.X[-1]))
        return utils.find(raw.X, raw.Y, energy, 'log')
    else:
        plt.plot(raw.X, raw.Y, linestyle='-.', color=color, marker='s', label='S. Adachi (1989)')


def Palik(wavelength, color='b'):
    """
    Input is wavelength (um), rather than energy(eV) or wavelength(nm).
    :param wavelength: (um)
    :param color: default is blue
    :return: absorption coefficient
    """
    location = file_location + 'InP_E.D.Palik.csv'
    raw = Data.CSV(data=location, xlabel='wavelength', ylabel='absorption')
    if wavelength is not None:
        if wavelength < np.amin(raw.X) or wavelength > np.amax(raw.X):
            raise BaseException("Wavelength should be in [%.2f, %.2f]" % (np.amin(raw.X), np.amax(raw.X)))
        return utils.find(raw.X, raw.Y, wavelength, 'log')
    else:
        plt.plot(raw.X, raw.Y, linestyle='-.', color=color, marker='s', label='E.D. Palik (1985)')