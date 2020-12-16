import numpy as np
import csv
import mpmath as mp
import math
import matplotlib
matplotlib.use('TkAgg')
import inspect
from scipy.optimize import curve_fit
import scipy.integrate as integrate
import os.path
import json
import utils
import time
start_time = time.time()

# Physics constants
e = 1.6e-19  # [Coul/e]
eps_InP = 12.5 * 8.85e-14  # [F/cm]
eps_InGaAs = 13.9 * 8.85e-14  # [F/cm] In 0.53 Ga 0.47 As
eps_InGaAsP = 13.436 * 8.85e-14  # [F/cm] Approximated by In 0.53 Ga 0.47 As 0.65 P 0.35
kB = 1.38e-23  # [J/K]
hbar = 1.054e-34  # [J-s]


def FermiIntegral_OneHalf(eta, num):
    def integrand(eta_F, energy):
        return 2 / np.sqrt(np.pi) * np.sqrt(energy) / (1 + np.exp(energy - eta_F))
    result = integrate.quad(lambda x: integrand(eta, x), 0, num)
    return result[0]


def HurkxTunnelingMass(temperature, current, vpt, voltage, guess, electricfield_Vcm, LifetimeSwitch):
    if LifetimeSwitch == 'On':
        def line(T, a, b, z, Eti):
            if abs(voltage) < abs(vpt):
                return a * T + b / (T ** 2) - Eti * e / kB - e * Eg_InP(T) / (2 * kB) - z * T * np.log(T / 300)
            else:
                return a * T + b / (T ** 2) - Eti * e / kB - e * Eg_InGaAs(T) / (2 * kB) - z * T * np.log(T / 300)

        popt, pcov = curve_fit(line, temperature, temperature * np.log(current), p0=guess)
        return (e * hbar * electricfield_Vcm * 100) ** 2 / (24 * kB ** 3 * tuple(popt)[1])
    else:
        def line(T, a, b, Eti):
            if abs(voltage) < abs(vpt):
                return a * T + b / (T ** 2) - Eti * e / kB - e * Eg_InP(T) / (2 * kB)
            else:
                return a * T + b / (T ** 2) - Eti * e / kB - e * Eg_InGaAs(T) / (2 * kB)

        popt, pcov = curve_fit(line, temperature, temperature * np.log(current), p0=guess)
        return (e * hbar * electricfield_Vcm * 100) ** 2 / (24 * kB ** 3 * tuple(popt)[1])


def get_Eti(temperature, current, vpt, voltage, formula, guess, switch='off'):
    if formula == 'Hurkx':
        def line(T, a, b, Eti):
            if abs(voltage) < abs(vpt):
                return a * T + b / (T ** 2) - Eti * e / kB - e * Eg_InP(T) / (2 * kB)
            else:
                return a * T + b / (T ** 2) - Eti * e / kB - e * Eg_InGaAs(T) / (2 * kB)

        popt, pcov = curve_fit(line, temperature, temperature * np.log(current), p0=guess)
        if switch == 'off':
            return tuple(popt)[2]
        else:
            return popt
    elif formula == 'FullHurkx':
        def line(T, a, b, Eti):
            if abs(voltage) < abs(vpt):
                return a * T + T * np.log((T/300) ** 1.5 + b * np.exp((300 * b) ** 2 / (12 * np.pi * T ** 3))) - \
                       Eti * e / kB - e * Eg_InP(T) / (2 * kB)
            else:
                return a * T + T * np.log((T/300) ** 1.5 + b * np.exp((300 * b) ** 2 / (12 * np.pi * T ** 3))) - \
                       Eti * e / kB - e * Eg_InGaAs(T) / (2 * kB)

        popt, pcov = curve_fit(line, temperature, temperature * np.log(current), p0=guess,
                               bounds=((5, 1e4, 0), (20, 1e5, 2)))
        if switch == 'off':
            return tuple(popt)[2]
        else:
            return popt
    elif formula == 'Hurkx+Lifetime':
        def line(T, a, b, Eti):
            z = -1.5
            if abs(voltage) < abs(vpt):
                return a * T + b / (T ** 2) - Eti * e / kB - e * Eg_InP(T) / (2 * kB) - z * T * np.log(T / 300)
            else:
                return a * T + b / (T ** 2) - Eti * e / kB - e * Eg_InGaAs(T) / (2 * kB) - z * T * np.log(T / 300)

        # bounds=((10, 1e10, 0, 0), (20, 1e11, 3, 2))
        popt, pcov = curve_fit(line, temperature, temperature * np.log(current), p0=guess)
        if switch == 'off':
            return tuple(popt)[2]
        else:
            return popt
    elif formula == 'SRH+Lifetime':
        def line(x, a, Eti, z):
            if abs(voltage) < abs(vpt):
                return a + z * np.log10(300) + (z - 1.5) * np.log10(x) - np.log10(np.e) * e * Eti * x / kB - \
                       np.log10(np.e) * e * Eg_InP(1 / x) / (2 * kB)
            else:
                return a + z * np.log10(300) + (z - 1.5) * np.log10(x) - np.log10(np.e) * e * Eti * x / kB - \
                       np.log10(np.e) * e * Eg_InGaAs(1 / x) / (2 * kB)

        popt, pcov = curve_fit(line, 1 / temperature, np.log10(current), p0=guess)
        if switch == 'off':
            return tuple(popt)[1]
        else:
            return popt
    else:
        def Eti(slope, V):
            if abs(V) < abs(vpt):
                return -slope / np.log10(np.e) * kB / e - 0.5 * Eg_InP(300)
            else:
                return -slope / np.log10(np.e) * kB / e - 0.5 * Eg_InGaAs(300)

        def line(x, a, b):
            if formula == '1':
                return a * x + b
            elif formula == '2':
                return a * x - 1.5 * np.log10(x) + b
            elif formula == '3':
                if abs(voltage) < abs(vpt):
                    return a * x - 1.5 * np.log10(x) + b - np.log10(np.e) / (2 * kB) * e * Eg_InP(1 / x) * x
                else:
                    return a * x - 1.5 * np.log10(x) + b - np.log10(np.e) / (2 * kB) * e * Eg_InGaAs(1 / x) * x
            else:
                raise BaseException("Wrong formula: %s" % formula)
        popt, pcov = curve_fit(line, 1 / temperature, np.log10(current))
        if formula == '1' or formula == '2':
            if switch == 'off':
                return Eti(tuple(popt)[0], voltage)
            else:
                return popt
        else:
            if switch == 'off':
                return - tuple(popt)[0] / np.log10(np.e) * kB / e
            else:
                return popt


def Eg_InP(T):
    """
    :param T: 絕對溫度 (K)
    :return: 能量 (eV)
    """
    alpha = 1.420
    beta = 0.000398
    gamma = 189
    return alpha - (beta * T ** 2) / (T + gamma)


def Eg_InGaAs(T):
    """
    :param T: 絕對溫度 (K)
    :return: 能量 (eV)
    """
    alpha = 0.8186754144533469
    beta = 0.000494420835956525
    gamma = 286.3332895176157
    return alpha - (beta * T ** 2) / (T + gamma)


def Eg_InGaAsP(x, y):
    # benzaquen1994alloy
    # Alloy broadening in photoluminescence spectra of GaxIn1−xAsyP1−y lattice matched to InP
    return 1.35 + 0.68 * x -1.068 * y + 0.758 * x ** 2 + 0.078 * \
           y ** 2 - 0.069 * x * y - 0.332 * x ** 2 * y + 0.3 * x * y ** 2


def Sze_Vb(ND, material):
    if material == 'InP':
        Eg = 1.42  # [eV]
    elif material == 'InGaAs':
        Eg = 0.752  # [eV]
    else:
        raise BaseException("No such material: %s" % material)
    return 60 * (Eg / 1.1) ** 1.5 * (ND / 1e16) ** (-3/4)


def dark_count_rate(area, position, prob_pair, g):
    return area * utils.ydx(position * 1e-4, prob_pair * g, 0, len(position) - 1)


def photon_detection_efficiency(position_um, absorption_coefficient, probability_pair, material):
    position_um, probability_pair = utils.Remesh(position_um, probability_pair, 0.001)
    index0 = utils.value_convert_index(-3.6, position_um)
    index1 = utils.value_convert_index(-3.5, position_um)
    index2 = utils.value_convert_index(-0.5, position_um)
    x = position_um
    x_wrt_top_um = np.asarray([x[i] - x[0] for i in range(len(x))])
    # 把 0.1um 的 InGaAsP 看作是完全透明！
    x_wrt_CG_um = np.asarray([x[i] - x[index0] for i in range(len(x))])
    absorb = absorption_coefficient
    Ppair = np.asarray(probability_pair)
    if material == 'InP':
        PDE_InP = utils.ydx(x * 1e-4, absorb * np.exp(- absorb * x_wrt_top_um * 1e-4) * Ppair, 0, index0) + \
                  utils.ydx(x * 1e-4, absorb * np.exp(- absorb * x_wrt_top_um * 1e-4) * Ppair, index2, len(x) - 1)
        return PDE_InP
    elif material == 'InGaAs':
        PDE_InGaAs = utils.ydx(x * 1e-4, absorb * np.exp(- absorb * x_wrt_CG_um * 1e-4) * Ppair, index1, index2)
        return PDE_InGaAs
    else:
        raise BaseException("Wrong material")


def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


def value_convert_index(value, _list):
    j = 0
    temp = value - _list[0]
    for i in range(len(_list)):
        if abs(value - _list[i]) < temp:
            j = i
            temp = abs(value - _list[i])
    return j


def savetmp(func, input, result):
    # [1] /Users/Ethan/PycharmProjects/Python/Practice/Cache_test.py
    # [2] /Users/Ethan/PycharmProjects/Python/Practice/
    # [3] tmp_Cache_test.csv
    # [4] /Users/Ethan/PycharmProjects/Python/Practice/data/
    # [5] /Users/Ethan/PycharmProjects/Python/Practice/data/tmp_Cache_test.csv
    location = inspect.stack()[1][1]  # [1]
    directory = os.getcwd() + '/'  # [2]
    filename = 'tmp_' + location.replace(directory, "").replace(".py", ".csv")  # [3]
    save_directory = directory + 'data/'  # [4]
    save_file = save_directory + filename  # [5]

    with open(save_file, newline='', encoding='utf-8-sig') as file2:
        rows = csv.DictReader(file2)
        for index, row in enumerate(rows):
            print(index, row['function_name'], row['result'])

    with open(save_file, 'w') as file:
        SaveData = csv.writer(file, delimiter=',')
        #SaveData.writerow(['function_name', 'result', 'input'])
        SaveData.writerow([func, input, result])


def orderfunction(x, y, type):
    if type == 'increasing':
        for i in range(len(x)):
            for j in range(i + 1, len(x)):
                if x[j] < x[i]:
                    temp1 = x[j]
                    temp2 = y[j]
                    x[j] = x[i]
                    y[j] = y[i]
                    x[i] = temp1
                    y[i] = temp2
    elif type == 'decreasing':
        for i in range(len(x)):
            for j in range(i + 1, len(x)):
                if x[j] > x[i]:
                    temp1 = x[j]
                    temp2 = y[j]
                    x[j] = x[i]
                    y[j] = y[i]
                    x[i] = temp1
                    y[i] = temp2
    else:
        raise BaseException("Wrong type : %s" % type)
    return x, y


def logBTB_JVcurve(x, F0, sigma, c):
    return -F0 * (x * 1e-6) - sigma * np.log(x * 1e-6) + np.log(c)


def logBTB_JVcurve_E2(x, F0, c):
    sigma = 1
    return -F0 * (x * 1e-6) - sigma * np.log(x * 1e-6) + np.log(c)


def BTBFitting(Emax, JV, lnJV_min, lnJV_max, boundset):
    if len(Emax) != len(JV):
        raise BaseException("Lengths of Emax and JV are not equal. len(Emax)=%s, len(JV)=%s" % (len(Emax), len(JV)))
    Emax_fit = []
    JV_BTB_fit = []
    for i in range(len(Emax)):
        if Emax[i] == 0:
            raise BaseException("Emax[%s] is zero." % i)
        elif JV[i] <= 0:
            raise BaseException("JV[%s]=%s is negative or zero." % (i, JV[i]))
        else:
            if lnJV_min < np.log(JV[i]) < lnJV_max:
                Emax_fit.append(Emax[i])  # [V/cm]
                JV_BTB_fit.append(JV[i])  # [A/(cm2-V)]

    if boundset[0][1] == 'E2':
        boundset[0].remove('E2')
        boundset[1].remove('E2')
        return curve_fit(logBTB_JVcurve_E2, 1 / np.asarray(Emax_fit) * 1e6, np.log(np.asarray(JV_BTB_fit)),
                         bounds=boundset)
    elif len(boundset[0]) == 3:
        return curve_fit(logBTB_JVcurve, 1 / np.asarray(Emax_fit) * 1e6, np.log(np.asarray(JV_BTB_fit)),
                         bounds=boundset)
    else:
        raise BaseException("Check boundset : %s" % boundset)


def save(location, xlabel, ylabel, xlist, ylist, form=',', title='Normal'):
    with open(location, 'w', newline='') as file:
        if form == ',':
            SaveData = csv.writer(file, delimiter=',')
            if title == 'Normal':
                SaveData.writerow([xlabel, ylabel])
            else:
                SaveData.writerow(title)
            for i in range(len(xlist)):
                SaveData.writerow([xlist[i], ylist[i]])
        elif form == '':
            SaveData = csv.writer(file, delimiter=' ')
            if title == 'Normal':
                SaveData.writerow([xlabel, ylabel])
            else:
                #SaveData.writerow(["""as"""])
                doubleQString = "{0}".format(title)
                json.dump(doubleQString, file)
                SaveData.writerow('')
            for i in range(len(xlist)):
                SaveData.writerow([xlist[i], ylist[i]])
    print('Save the file: %s' % location)


def interpolation(xi, xf, yi):
    for i in range(1, len(xi) - 1):
        if (xi[i + 1] - xi[i]) * (xi[i] - xi[i - 1]) < 0:
            raise BaseException("xi is not a monotonic function where "
                                "xi[%s], xi[%s], xi[%s] are %s, %s, %s, respectively."
                                % (i - 1, i, i + 1, xi[i - 1], xi[i], xi[i + 1]))
        elif xi[i] == xi[i + 1]:
            raise BaseException("xi[%s] = xi[%s] = %s and yi[%s]=%s, yi[%s]=%s" %
                                (i, i + 1, xi[i], i, yi[i], i + 1, yi[i + 1]))
    yf = np.empty([0])
    for f in range(len(xf)):
        if f != len(yf):
            raise BaseException("Error")
        if xf[f] < xi[0]:
            dx = xf[f] - xi[0]
            dx0 = xi[1] - xi[0]
            dy0 = yi[1] - yi[0]
            y_fit = yi[0] + dx * dy0 / dx0
            yf = np.append(yf, y_fit)
        elif xf[f] > xi[-1]:
            dx = xf[f] - xi[-1]
            dx_last = xi[-1] - xi[-2]
            dy_last = yi[-1] - yi[-2]
            y_fit = yi[-1] + dx * dy_last / dx_last
            '''
            # 在 89-4 中，由於這是 doping profile 的外插，所以要求緩衝層濃度最高為 1e18
            # 並且觀察 Emax(V) 與擬合結果是否隨之變化，後來都沒有發現明顯的影響。
            if y_fit > 1e18:
                y_fit = 1e18
            '''
            yf = np.append(yf, y_fit)
        else:
            for i in range(len(xi)):
                if xf[f] == xi[i]:
                    yf = np.append(yf, yi[i])
                else:
                    if 0 <= i < len(xi):
                        if xi[i] < xf[f] < xi[i + 1]:
                            a = abs(xf[f] - xi[i])
                            b = abs(xi[i + 1] - xf[f])
                            dyi = yi[i + 1] - yi[i]
                            y_fit = yi[i] + a * dyi / (a + b)
                            yf = np.append(yf, y_fit)
    return yf


def merge(v1, v2, form='np.array'):
    set_v1 = set(v1)
    set_v2 = set(v2)
    set_dif = set_v2 - set_v1
    v = list(v1) + list(set_dif)
    v.sort(reverse=False)
    if form == 'list':
        return v
    elif form == 'np.array':
        return np.asarray(v)
    else:
        raise BaseException("Wrong form.")


def multiple_merge(x):
    y = []
    for i in range(len(x)):
        if x[i] not in y:
            y.append(x[i])
    for i in range(len(y)):
        for j in range(i + 1, len(y)):
            if y[j] < y[i]:
                temp = y[j]
                y[j] = y[i]
                y[i] = temp
    return np.asarray(y)


def find(x, func, x_value, form):
    for i in range(len(x)):
        if x[i] < x[i + 1]:
            curve = 'increase'
            break
        elif x[i] > x[i + 1]:
            curve = 'decrease'
            break

    for i in range(len(x) - 1):
        #print('find x_value: ', x[i], x_value, x[i + 1], curve)
        if x[i] == x_value:
            return func[i]
        if (x_value > x[-1]) and (curve == 'increase'):
            goal = (func[-1] - func[-2]) / (x[-1] - x[-2]) * (x_value - x[-1])
            return goal
        if (x[i] < x_value < x[i + 1]) and (curve == 'increase'):
            a = x_value - x[i]
            b = x[i + 1] - x_value
            c1 = b / (a + b)
            c2 = a / (a + b)
            if form == 'linear':
                goal = c1 * func[i] + c2 * func[i + 1]
                return goal
            elif form == 'log':
                sign = func[i] / abs(func[i])
                goal = sign * 10 ** (c1 * np.log10(abs(func[i])) + c2 * np.log10(abs(func[i + 1])))
                return goal
            else:
                print('Check "goal" term.')
                return
        if (x[i + 1] < x_value < x[i]) and (curve == 'decrease'):
            a = x_value - x[i]
            b = x[i + 1] - x_value
            c1 = b / (a + b)
            c2 = a / (a + b)
            if form == 'linear':
                goal = c1 * func[i] + c2 * func[i + 1]
                return goal
            elif form == 'log':
                sign = func[i] / abs(func[i])
                goal = sign * 10 ** (c1 * np.log10(abs(func[i])) + c2 * np.log10(abs(func[i + 1])))
                return goal
            else:
                print('Check "goal" term.')
                return


def dydx(x, y):
    dydx = [(y[1] - y[0]) / (x[1] - x[0])]
    for i in range(1, len(y) - 1):
        dx1 = x[i] - x[i - 1]
        dy1 = y[i] - y[i - 1]
        dx2 = x[i + 1] - x[i]
        dy2 = y[i + 1] - y[i]
        dydx1 = dy1 / dx1
        dydx2 = dy2 / dx2
        dydx_avg = 0.5 * (dydx1 + dydx2)
        dydx.append(dydx_avg)
    dydx.append((y[-1] - y[-2]) / (x[-1] - x[-2]))
    return np.asarray(dydx)


def center_smooth(func, k):
    ysmooth = np.zeros(len(func))
    for i in range(0, len(func)):
        if k <= i <= len(func) - k - 1:
            for j in range(i - k, i + k + 1):
                ysmooth[i] = ysmooth[i] + func[j] / (2 * k + 1)
        elif i < k:
            for j in range(i + k + 1):
                ysmooth[i] = ysmooth[i] + func[j] / (i + k + 1)
        else:
            for j in range(i - k, len(func)):
                ysmooth[i] = ysmooth[i] + func[j] / (len(func) - i + k)
    return ysmooth


def partial_center_smooth(s1, s2, func, k, direction):
    if s2 == -1:
        s2 = len(func)
    ysmooth = np.zeros(len(func))
    if direction == 'f':
        for i in range(0, len(func)):
            if i < s1 or i > s2:
                ysmooth[i] = func[i]
            else:
                if k <= i <= len(func) - k - 1:
                    for j in range(i - k, i + k + 1):
                        ysmooth[i] = ysmooth[i] + func[j] / (2 * k + 1)
                elif i < k:
                    for j in range(i + k + 1):
                        ysmooth[i] = ysmooth[i] + func[j] / (i + k + 1)
                else:
                    for j in range(i - k, len(func)):
                        ysmooth[i] = ysmooth[i] + func[j] / (len(func) - i + k)
        return ysmooth
    elif direction == 'r':
        for i in range(0, len(func)):
            if i > len(func) - s1 or i < len(func) - s2:
                ysmooth[i] = func[i]
            else:
                if k <= i <= len(func) - k - 1:
                    for j in range(i - k, i + k + 1):
                        ysmooth[i] = ysmooth[i] + func[j] / (2 * k + 1)
                elif i < k:
                    for j in range(i + k + 1):
                        ysmooth[i] = ysmooth[i] + func[j] / (i + k + 1)
                else:
                    for j in range(i - k, len(func)):
                        ysmooth[i] = ysmooth[i] + func[j] / (len(func) - i + k)
        return ysmooth
    else:
        print('Wrong direction in partial_center_smooth function.')
        return 'error'


def parse_csv(csv_path, key, Cname, X, floatKey):
    result = dict()
    Curve_name = np.empty([0])
    with open(csv_path, newline='', encoding='utf-8-sig') as csv_file:
        rows = csv.DictReader(csv_file)
        fieldnames = rows.fieldnames
        # csvreader.fieldnames
        # If not passed as a parameter when creating the object,
        # this attribute is initialized upon first access or when
        # the first record is read from the file.
        if len(fieldnames[0].rsplit(' ', 1)) == 2:
            ## init  x,y list of curve of result
            for fieldname in fieldnames:
                curve_name, x_y = fieldname.rsplit(' ', 1)
                if curve_name not in result:
                    result[curve_name] = dict()
                    Curve_name = np.append(Curve_name, curve_name)
                result[curve_name][x_y] = []
            ## parse each field of each row
            for _index, row in enumerate(rows):
                ## format = "{curvename} [X|Y]"
                # print(_index)
                for fieldname in fieldnames:
                    if row[fieldname] == '-':
                        continue
                    curve_name, x_y = fieldname.rsplit(' ', 1)
                    if floatKey == 'Yes':
                        try:
                            result[curve_name][x_y].append(float(row[fieldname]))
                        except ValueError:
                            if len(row[fieldname]) == 0:
                                raise BaseException("row[%s][%s] is empty, that is, len(row[%s][%s])=0"
                                                    % (fieldname, _index, fieldname, _index))
                            else:
                                raise BaseException("ValueError: row[%s]=%s" % (fieldname, row[fieldname]))
                    else:
                        result[curve_name][x_y].append(row[fieldname])
        else:
            for fieldname in fieldnames:
                if fieldname not in result:
                    result[fieldname] = dict()
                    Curve_name = np.append(Curve_name, fieldname)
                result[fieldname] = []
            for _index, row in enumerate(rows):
                for fieldname in fieldnames:
                    if row[fieldname] == '-':
                        continue
                    if floatKey == 'Yes':
                        if isfloat(row[fieldname]) and len(row[fieldname]) > 0:
                            try:
                                result[fieldname].append(float(row[fieldname]))
                            except ValueError:
                                if len(row[fieldname]) == 0:
                                    raise BaseException("row[%s][%s] is empty, that is, len(row[%s][%s])=0"
                                                        % (fieldname, _index, fieldname, _index))
                                else:
                                    raise BaseException("ValueError: row[%s]=%s" % (fieldname, row[fieldname]))
                    else:
                        result[fieldname].append(row[fieldname])
    if key == 1:
        return result
    elif key == 2:
        return Curve_name
    elif key == 3:
        return result[Cname][X]
    else:
        return result, Curve_name


# Define functions of impact ionization coefficient
def alpha(electric_field, position, interface, form_InP='vanOverst', temperature=300):
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
                if abs(electric_field) < 3.8e5:
                    return 1.12e7 * np.exp(-3.11e6 / abs(electric_field))
                else:
                    return 2.93e6 * np.exp(-2.64e6 / abs(electric_field))
            else:
                raise
        else:
            return 0


def beta(electric_field, position, interface, form_InP='vanOverst', temperature=300):
    # InP / InGaAsP / InGaAs / InP
    # interface = [x_InP_InGaAsP, x_InGaAsP_InGaAs, x_InGaAs_InP]
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
                return a * (1 + c * (temperature - T0)) * abs(electric_field) ** gamma * np.exp(
                    - ((b * (1 + d * (temperature - T0))) / abs(electric_field)) ** delta)
            elif form_InP == 'vanOverst':
                if electric_field < 3.8e5:
                    return 4.79e6 * np.exp(-2.55e6 / abs(electric_field))
                else:
                    return 1.62e6 * np.exp(-2.11e6 / abs(electric_field))
            else:
                raise BaseException("Wrong form")
        else:
            return 0


def eps(position, interface):
    # InP / InGaAsP / InGaAs / InP
    # interface = [x_InP_InGaAsP, x_InGaAsP_InGaAs, x_InGaAs_InP]
    if position < interface[0] or position > interface[2]:
        return eps_InP
    elif interface[0] < position < interface[1]:
        return eps_InGaAsP
    else:
        return eps_InGaAs


def lifetime_Scharfetter(position, interface, mesh, total_doping):
    t_min_InP = [9e-9, 1.5e-9]
    t_min_InGaAs = [0.022e-9, 0.01e-9]
    t_min_InGaAsP = [0.5e-9, 0.85e-9]
    t_max_InP = [315e-9, 550e-9]
    t_max_InGaAs = [0.25e-9, 8e-9]
    t_max_InGaAsP = [1e-5, 1.5e-5]
    Nref_InP = [3.7e17, 4e15]
    Nref_InGaAs = [2.75e18, 1.7e18]
    Nref_InGaAsP = [1e17, 5.3e16]
    gamma_InP = [1.7, 2.2]
    gamma_InGaAs = [2.5, 1.9]
    gamma_InGaAsP = [2, 1.52]

    # InP / InGaAsP / InGaAs / InP
    # interface = [x_InP_InGaAsP, x_InGaAsP_InGaAs, x_InGaAs_InP]
    k, xk = position_converter(position, mesh, 0)
    if xk < interface[0] or xk > interface[2]:
        elifetime = t_min_InP[0] + (t_max_InP[0] - t_min_InP[0])/\
                    (1 + (total_doping[k]/Nref_InP[0]) ** gamma_InP[0])
        hlifetime = t_min_InP[1] + (t_max_InP[1] - t_min_InP[1])/\
                    (1 + (total_doping[k]/Nref_InP[1]) ** gamma_InP[1])
    elif interface[0] <= xk <= interface[1]:
        elifetime = t_min_InGaAsP[0] + (t_max_InGaAsP[0] - t_min_InGaAsP[0]) / \
                    (1 + (total_doping[k] / Nref_InGaAsP[0]) ** gamma_InGaAsP[0])
        hlifetime = t_min_InGaAsP[1] + (t_max_InGaAsP[1] - t_min_InGaAsP[1]) / \
                    (1 + (total_doping[k] / Nref_InGaAsP[1]) ** gamma_InGaAsP[1])
    else:
        elifetime = t_min_InGaAs[0] + (t_max_InGaAs[0] - t_min_InGaAs[0]) / \
                    (1 + (total_doping[k] / Nref_InGaAs[0]) ** gamma_InGaAs[0])
        hlifetime = t_min_InGaAs[1] + (t_max_InGaAs[1] - t_min_InGaAs[1]) / \
                    (1 + (total_doping[k] / Nref_InGaAs[1]) ** gamma_InGaAs[1])
    return elifetime, hlifetime


def mobility(position, interface, mesh, total_doping, E):
    # InP / InGaAsP / InGaAs / InP
    # interface = [x_InP_InGaAsP, x_InGaAsP_InGaAs, x_InGaAs_InP]
    # E is the electric field distribution for the specific terminal voltage
    k, xk = position_converter(position, mesh, 0)
    a_min_InP = [400, 10]
    a_d_InP = [4800, 160]
    a_N_InP = [3e17, 4.87e17]
    a_a_InP = [0.47, 0.62]
    a_min_InGaAs = [1600, 75]
    a_d_InGaAs = [1e4, 256]
    a_N_InGaAs = [2.1e17, 7.7e18]
    a_a_InGaAs = [0.76, 1.37]
    a_min_InGaAsP = [6e3, 5e3]
    a_d_InGaAsP = [0, 0]
    a_N_InGaAsP = [1e16, 1e16]
    a_a_InGaAsP = [1, 1]
    E0_InP = [4000, 4000]
    E0_InGaAsP = [4000, 4000]
    E0_InGaAs = [4000, 4000]
    vsat0_InP = [7.5e7, 8.0e6]
    vsat0_InGaAsP = [8e6, 6.25e6]
    vsat0_GaAs = [7.5e7, 8.0e6]
    vsat0_InAs = [1.636e7, 1.525e7]
    vsat0_InGaAs = [0.47 * vsat0_GaAs[0] + 0.53 * vsat0_InAs[0], 0.47 * vsat0_GaAs[1] + 0.53 * vsat0_InAs[1]]
    if xk < interface[0] or xk > interface[2]:
        miu_e = a_min_InP[0] + a_d_InP[0] / (1 + (total_doping[k] / a_N_InP[0]) ** a_a_InP[0])
        miu_h = a_min_InP[1] + a_d_InP[1] / (1 + (total_doping[k] / a_N_InP[1]) ** a_a_InP[1])
        if E[k] == 0:
            miu_hf_e = miu_e
            miu_hf_h = miu_h
        else:
            miu_hf_e = (miu_e + (vsat0_InP[0] / abs(E[k])) * (E[k] / E0_InP[0]) ** 4) / (1 + (E[k] / E0_InP[0]) ** 4)
            miu_hf_h = (miu_h + (vsat0_InP[1] / abs(E[k])) * (E[k] / E0_InP[1]) ** 4) / (1 + (E[k] / E0_InP[1]) ** 4)
        return miu_e, miu_h, miu_hf_e, miu_hf_h
    elif interface[0] <= xk <= interface[1]:
        miu_e = a_min_InGaAsP[0] + a_d_InGaAsP[0] / (1 + (total_doping[k] / a_N_InGaAsP[0]) ** a_a_InGaAsP[0])
        miu_h = a_min_InGaAsP[1] + a_d_InGaAsP[1] / (1 + (total_doping[k] / a_N_InGaAsP[1]) ** a_a_InGaAsP[1])
        if E[k] == 0:
            miu_hf_e = miu_e
            miu_hf_h = miu_h
        else:
            miu_hf_e = (miu_e + (vsat0_InGaAsP[0] / abs(E[k])) * (E[k] / E0_InGaAsP[0]) ** 4) / (1 + (E[k] / E0_InGaAsP[0]) ** 4)
            miu_hf_h = (miu_h + (vsat0_InGaAsP[1] / abs(E[k])) * (E[k] / E0_InGaAsP[1]) ** 4) / (1 + (E[k] / E0_InGaAsP[1]) ** 4)
        return miu_e, miu_h, miu_hf_e, miu_hf_h
    else:
        miu_e = a_min_InGaAs[0] + a_d_InGaAs[0] / (1 + (total_doping[k] / a_N_InGaAs[0]) ** a_a_InGaAs[0])
        miu_h = a_min_InGaAs[1] + a_d_InGaAs[1] / (1 + (total_doping[k] / a_N_InGaAs[1]) ** a_a_InGaAs[1])
        if E[k] == 0:
            miu_hf_e = miu_e
            miu_hf_h = miu_h
        else:
            miu_hf_e = (miu_e + (vsat0_InGaAs[0] / abs(E[k])) * (E[k] / E0_InGaAs[0]) ** 4) / (1 + (E[k] / E0_InGaAs[0]) ** 4)
            miu_hf_h = (miu_h + (vsat0_InGaAs[1] / abs(E[k])) * (E[k] / E0_InGaAs[1]) ** 4) / (1 + (E[k] / E0_InGaAs[1]) ** 4)
        return miu_e, miu_h, miu_hf_e, miu_hf_h

# Non-uniform mesh
# mesh_list : (P+)InP / M-Layer / InGaAsP / InGaAs / Buffer & Substrate
def mesh(n, line, x_raw, doping_raw):
    total_mesh = 0
    ix = np.zeros(7, dtype=int)  # temporary indices of position
    for j in range(len(n)):
        total_mesh = total_mesh + n[j]
    x = np.zeros(total_mesh)  # new mesh array
    doping = np.zeros(total_mesh)  # new doping profile in terms of new mesh array
    for position in range(0, total_mesh):
        # InP top bulk region  (p+ InP)
        if position < n[0]:
            x[position] = x_raw[0] + position * (line[0] - x_raw[0]) / n[0]
            ix[0] = position
        # Multiplication layer (p+ InP)
        elif n[0] <= position < n[0] + n[1]:
            x[position] = x[ix[0]] + (position - ix[0]) * (line[1] - x[ix[0]]) / n[1]
            ix[1] = position
        # Multiplication layer (n- InP)
        elif n[0] + n[1] <= position < n[0] + n[1] + n[2]:
            x[position] = x[ix[1]] + (position - ix[1]) * (line[2] - x[ix[1]]) / n[2]
            ix[2] = position
        # Charge layer (n+ InP)
        elif n[0] + n[1] + n[2] <= position < n[0] + n[1] + n[2] + n[3]:
            x[position] = x[ix[2]] + (position - ix[2]) * (line[3] - x[ix[2]]) / n[3]
            ix[3] = position
        # Grading layer (InGaAsP)
        elif n[0] + n[1] + n[2] + n[3] <= position < n[0] + n[1] + n[2] + n[3] + n[4]:
            x[position] = x[ix[3]] + (position - ix[3]) * (line[4] - x[ix[3]]) / n[4]
            ix[4] = position
        # Absorption layer (InGaAs)
        elif n[0] + n[1] + n[2] + n[3] + n[4] <= position < n[0] + n[1] + n[2] + n[3] + n[4] + n[5]:
            x[position] = x[ix[4]] + (position - ix[4]) * (line[5] - x[ix[4]]) / n[5]
            ix[5] = position
        # Buffer layer (n- InP)
        elif n[0] + n[1] + n[2] + n[3] + n[4] + n[5] \
                <= position < n[0] + n[1] + n[2] + n[3] + n[4] + n[5] + n[6]:
            x[position] = x[ix[5]] + (position - ix[5]) * (line[6] - x[ix[5]]) / n[6]
            ix[6] = position
        else:  # Substrate (n+ InP)
            x[position] = x[ix[6]] + (position - ix[6]) * (x_raw[len(x_raw) - 1] - x[ix[6]]) / n[7]
    doping_temp = 0
    for i in range(0, total_mesh):
        for k in range(0, len(x_raw) - 1):
            if x_raw[k] != x_raw[k+1]:
                if x_raw[k] <= x[i] <= x_raw[k + 1]:
                    if doping_raw[k] * doping_raw[k + 1] > 0:
                        doping_temp = math.log10(abs(float(doping_raw[k]))) + \
                                      (math.log10(abs(float(doping_raw[k + 1]))) -
                                       math.log10(abs(float(doping_raw[k])))) / float(x_raw[k + 1] -
                                                                           x_raw[k]) * float(x[i] - x_raw[k])
                        if doping_raw[k] > 0:
                            doping[i] = 10 ** doping_temp
                        if doping_raw[k] < 0:
                            doping[i] = - 10 ** doping_temp
                        doping_temp = 0
                    elif doping_raw[k] * doping_raw[k + 1] < 0:
                        doping[i] = doping_raw[k] + float(doping_raw[k + 1] - doping_raw[k]) / float(
                            x_raw[k + 1] - x_raw[k]) * float(x[i] - x_raw[k])
                    else:
                        d1 = doping_raw[k]
                        d2 = doping_raw[k + 1]
                        if doping_raw[k] == 0:
                            d1 = -20
                        if doping_raw[k + 1] == 0:
                            d2 = -20
                        doping_temp = d1 + (d2 - d1) / float(x_raw[k + 1] - x_raw[k]) * float(x[i] - x_raw[k])
                        doping[i] = 10 ** doping_temp
                        doping_temp = 0
                    # doping[i] = doping_raw[k] + float(doping_raw[k + 1] - doping_raw[k]) / float(
                    #   x_raw[k + 1] - x_raw[k]) * float(
                    #   x[i] - x_raw[k])
    return x, doping


# Construct depletion boundary data
def depletion_boundary(n_dep, xp_min, xp_max, x, doping, interface, v_goal, form, rj, xj):
    # InP / InGaAsP / InGaAs / InP
    # interface = [x_InP_InGaAsP, x_InGaAsP_InGaAs, x_InGaAs_InP]
    if xp_max != xp_min:
        xp = np.zeros(n_dep)  # Depletion boundary at p region
        for i in range(0, n_dep):
            xp[i] = xp_min + i * (xp_max - xp_min) / n_dep
    if xp_max == xp_min:
        xp = [xp_max]
        n_dep = 1

    # E(x), xn, xp
    ef_depbd = 0
    efcl_depbe = 0
    xn = np.zeros(n_dep)  # Depletion boundary at n region
    E = np.zeros((n_dep, len(x)))
    for i in range(0, len(xp)):
        for k in range(0, len(x) - 1):
            if x[k] > xp[i]:
                temp = ef_depbd
                temp_cl = efcl_depbe
                if form == 'planar':
                    l1 = doping[k] / eps(x[k], interface)
                    l2 = doping[k + 1] / eps(x[k + 1], interface)
                    ef_depbd = ef_depbd + e * (l1 + l2) * (x[k + 1] - x[k]) * 1e-4 / 2
                    E[i, k] = ef_depbd
                    if temp * ef_depbd < 0:
                        xn[i] = x[k]
                        temp = 0
                        ef_depbd = 0
                        E[i, k] = 0
                        l1 = 0
                        l2 = 0
                        break
                    if ef_depbd == 0:
                        temp = 0
                        ef_depbd = 0
                        l1 = 0
                        l2 = 0
                        break
                elif form == 'cylindrical':
                    l1 = (rj + x[k] - xj) * e * doping[k] / eps(x[k], interface) * 1e-4
                    l2 = (rj + x[k + 1] - xj) * e * doping[k + 1] / eps(x[k + 1], interface) * 1e-4
                    efcl_depbe = efcl_depbe + (l1 + l2) * (x[k + 1] - x[k]) * 1e-4 / 2
                    E[i, k] = efcl_depbe / ((rj + x[k] - xj) * 1e-4)
                    if temp_cl * efcl_depbe < 0:
                        xn[i] = x[k]
                        temp_cl = 0
                        efcl_depbe = 0
                        E[i, k] = 0
                        l1 = 0
                        l2 = 0
                        break
                    if efcl_depbe == 0:
                        temp_cl = 0
                        efcl_depbe = 0
                        l1 = 0
                        l2 = 0
                        break
                else:
                    return 'Check form', 'Check form', 'Check form', \
                           'Check form', 'Check form', 'Check form', \
                           'Check form', 'Check form'
        l1 = 0
        l2 = 0
    # "Measured" voltage
    xp_goal = 0
    xn_goal = 0
    v_fit_goal = 0
    ef_goal = 0
    voltage = np.zeros(n_dep)  # V(substrate) - V(anode)
    for i in range(0, len(voltage)):
        for k in range(0, len(x) - 1):
            l1 = - E[i, k]
            l2 = - E[i, k + 1]
            voltage[i] = voltage[i] + (l1 + l2) / 2 * (x[k + 1] - x[k]) * 1e-4
        l1 = 0
        l2 = 0
    if v_goal != 0:
        i_goal_temp = 0
        for i in range(0, len(voltage)):
            if abs(voltage[i_goal_temp] - v_goal) >= abs(voltage[i] - v_goal):
                i_goal_temp = i
                xp_goal = xp[i_goal_temp]
                xn_goal = xn[i_goal_temp]
                v_fit_goal = voltage[i_goal_temp]
        ef_goal = np.zeros(len(x))
        for i in range(len(x)):
            ef_goal[i] = E[i_goal_temp, i]
    if v_goal == 0:
        xp_goal = "Error: lacking of voltage goal"
        xn_goal = "Error: lacking of voltage goal"
        v_fit_goal = "Error: lacking of voltage goal"
        ef_goal = "Error: lacking of voltage goal"
    return xp, xn, voltage, E, xp_goal, xn_goal, v_fit_goal, ef_goal


def structure_analysis(x_raw, doping_raw, conc):
    # Definition of concentration
    # concentration = [ML, Charge_max, Grading, Absorption, Buffer]
    xj = 0
    xj_mc = 0
    xj_ga = 0
    xj_ab = 0
    xj_cg = 0
    i_count = 0
    cg_marker = 0
    for i in range(0, len(x_raw)-1):
        # Find the PN junction
        if (doping_raw[i + 1] * doping_raw[i] < 0 and
                i_count == 0):
            xj = x_raw[i+1]
            i_count = 1
        # Find the interface between multiplication and charge layer
        if conc[0] <= doping_raw[i] < conc[1] and i_count == 1:
            xj_mc = x_raw[i-1]
            i_count = 2
        # Find InGaAsP / InGaAs interface
        if doping_raw[i] < conc[3] <= doping_raw[i+1] and i_count == 2:
            xj_ga = (x_raw[i] + x_raw[i+1]) / 2
            i_count = 3
        # Find InGaAs / Buffer interface
        if doping_raw[i] < conc[4] <= doping_raw[i+1] and i_count == 3:
            xj_ab = (x_raw[i] + x_raw[i+1]) / 2
            i_count = 4
        # Find charge layer (InP) & Grading layer (InGaAsP) interface
        if doping_raw[i] > conc[0] and doping_raw[i+1] < conc[3] and cg_marker == 0:
            xj_cg = (x_raw[i] + x_raw[i+1]) / 2
            cg_marker = 1
    ml_thick = abs(xj - xj_mc)
    return xj, xj_mc, xj_cg, xj_ga, xj_ab, ml_thick


def multiplication_factor(x, ef, interface, model='vanOverst'):
    # Unit of input x in multiplication_factor is [um]
    # InP / InGaAsP / InGaAs / InP
    # interface = [x_InP_InGaAsP, x_InGaAsP_InGaAs, x_InGaAs_InP]
    un = np.zeros(len(x))
    up = np.zeros(len(x))
    temp_un = 0
    temp_up = 0
    temp_mn = 0
    temp_mp = 0
    side = np.zeros(2)
    for i in range(0, len(x)):
        for k in range(0, i):
            side[0] = - (alpha(ef[k + 1], x[k + 1], interface, form_InP=model) -
                         beta(ef[k + 1], x[k + 1], interface, form_InP=model))
            side[1] = - (alpha(ef[k], x[k], interface, form_InP=model) - beta(ef[k], x[k], interface, form_InP=model))
            temp_un = temp_un + (side[0] + side[1]) * (x[k + 1] - x[k]) / 2 * 1e-4
        side[0] = 0
        side[1] = 0
        for k in range(i, len(x) - 1):
            side[0] = (alpha(ef[k + 1], x[k + 1], interface, form_InP=model) -
                       beta(ef[k + 1], x[k + 1], interface, form_InP=model))
            side[1] = (alpha(ef[k], x[k], interface, form_InP=model) - beta(ef[k], x[k], interface, form_InP=model))
            temp_up = temp_up + (side[0] + side[1]) * (x[k + 1] - x[k]) / 2 * 1e-4
        un[i] = np.exp(temp_un)
        up[i] = np.exp(temp_up)
        temp_un = 0
        temp_up = 0
        side[0] = 0
        side[1] = 0
    for i in range(0, len(x) - 1):
        side1n = alpha(ef[i], x[i], interface, form_InP=model) * un[i]
        side2n = alpha(ef[i + 1], x[i + 1], interface, form_InP=model) * un[i + 1]
        temp_mn = temp_mn + (side1n + side2n) * (x[i + 1] - x[i]) * 1e-4 / 2
        side1p = beta(ef[i], x[i], interface, form_InP=model) * up[i]
        side2p = beta(ef[i + 1], x[i + 1], interface, form_InP=model) * up[i + 1]
        temp_mp = temp_mp + (side1p + side2p) * (x[i + 1] - x[i]) * 1e-4 / 2
    return temp_mn, temp_mp


def vb_solver(trial_number, step, error, x, e_goal, interface, s, v_number, doping, form, rj, xj, v0_goal):
    gate = 0
    v = v0_goal
    mn, mp = multiplication_factor(x, e_goal, interface)
    error_n = abs(mn - 1)
    error_p = abs(mp - 1)
    print('[0]', '  V: %.2f' % v0_goal, '  Mn: %.3f' % mn, '  Mp: %.3f' % mp,
          '  Error: %.3f' % max(error_n, error_p), '--- %.5s seconds ---' % (time.time() - start_time))
    for i in range(1, trial_number):
        end = trial_number - 1
        if mn > 1 + error or mp > 1 + error:
            if gate < 0:
                step = step / 2
            s = s + step
            xp1, xn1, v, ef, xp_goal, xn_goal, v_goal, e_goal = \
                depletion_boundary(v_number, s, s, x, doping, interface, 0, form, rj, xj)
            mn, mp = multiplication_factor(x, ef[0, :], interface)
            error_n = abs(mn - 1)
            error_p = abs(mp - 1)
            # print(type(i), type(v_goal), type(mn), type(mp))
            print('[%d]' % i, '  V: %.2f' % v[0], '  Mn: %.3f' % mn, '  Mp: %.3f' % mp,
                  '  Error: %.3f' % max(error_n, error_p), '--- %.5s seconds ---  (+s)' % (time.time() - start_time))
            gate = 1
            if i == end:
                error_n = abs(mn - 1)
                error_p = abs(mp - 1)
                print('Iterations larger than maximum %.3f' % v)
                print('[%d]' % i, '  V: %.2f' % v[0], '  Mn: %.3f' % mn, '  Mp: %.3f' % mp,
                      '  Error: %.3f' % max(error_n, error_p), '--- %.5s seconds ---' % (time.time() - start_time))
                return v, mn, mp, i
        elif mn < 1 - error or mp < 1 - error:
            if gate > 0:
                step = step / 2
            s = s - step
            xp1, xn1, v, ef, xp_goal, xn_goal, v_goal, e_goal = \
                depletion_boundary(v_number, s, s, x, doping, interface, 0, form, rj, xj)
            mn, mp = multiplication_factor(x, ef[0, :], interface)
            error_n = abs(mn - 1)
            error_p = abs(mp - 1)
            print('[%d]' % i, '  V: %.2f' % v[0], '  Mn: %.3f' % mn, '  Mp: %.3f' % mp,
                  '  Error: %.3f' % max(error_n, error_p), '--- %.5s seconds ---  (-s)' % (time.time() - start_time))
            gate = -1
            if i == end:
                error_n = abs(mn - 1)
                error_p = abs(mp - 1)
                print('Iterations larger than maximum %.3f' % v)
                print('[%d]' % i, '  V: %.2f' % v[0], '  Mn: %.3f' % mn, '  Mp: %.3f' % mp,
                      '  Error: %.3f' % max(error_n, error_p), '--- %.5s seconds ---' % (time.time() - start_time))
                return v, mn, mp, i
        else:
            return v, mn, mp, i-1


def position_converter(position, x_raw, s):
    if position >= x_raw[len(x_raw) - 1]:
        if s == 0:
            return len(x_raw) - 1, x_raw[len(x_raw) - 1]
        elif s == 1:
            return int(len(x_raw) - 1)
        else:
            return x_raw[len(x_raw) - 1]
    if position < x_raw[0]:
        if s == 0:
            return 0, x_raw[0]
        elif s == 1:
            return 0
        else:
            return x_raw[0]
    for i in range(len(x_raw)):
        if i+1 < len(x_raw):
            if x_raw[i] <= position < x_raw[i+1]:
                if s == 0:
                    return i, x_raw[i]
                elif s == 1:
                    return int(i)
                else:
                    return x_raw[i]


def ni(position, interface):
    # InP / InGaAsP / InGaAs / InP
    # interface = [x_InP_InGaAsP, x_InGaAsP_InGaAs, x_InGaAs_InP]
    if interface[0] <= position < interface[1]:  # InGaAsP
        # http://www.ioffe.ru/SVA/NSM/Semicond/GaInAsP/bandstr.html
        # Yamazoe et al. (1981) 看圖內插得來的數據 2e8，1000/T = 1000/300 ~ 3.3
        return 2e8  # (cm-3)
    elif interface[1] <= position < interface[2]:  # InGaAs
        return 6.3e11  # (cm-3) http://www.ioffe.ru/SVA/NSM/Semicond/GaInAs/bandstr.html
    else:  # InP
        return 1.3e7  # (cm-3) http://www.ioffe.ru/SVA/NSM/Semicond/InP/bandstr.html


def SRH_generation_rate(ni, tau_p0, tau_n0, et, ei):
    s1 = (et - ei) / 0.0259  # Et & Ev are in [eV]
    s2 = 0.5 * np.log(tau_p0 / tau_n0)
    return ni / (2 * (tau_p0 * tau_n0) ** 0.5 * np.cosh(s1 + s2))


def B2B_tunneling_rate(ni, tau_p0, tau_n0, et, ei, ef, A, p, B):
    s1 = (et - ei) / 0.0259  # Et & Ev are in [eV]
    s2 = 0.5 * np.log(tau_p0 / tau_n0)
    return ni / (2 * (tau_p0 * tau_n0) ** 0.5 * np.cosh(s1 + s2)) * A * ef ** p * np.exp(B/ef)  # ef in (V/cm)


def avalanche_generation(electric_field, position, mesh, interface, total_doping, E_v, n, p):
    # E_v : Earray at V (1D)
    an = alpha(electric_field, position, interface)
    ap = beta(electric_field, position, interface)
    miu_e, miu_h, miu_hf_n, miu_hf_p = mobility(position, interface, mesh, total_doping, E_v)
    k = position_converter(position, mesh, 1)
    Gava_n = an * n * miu_hf_n * abs(E_v[k])
    Gava_p = ap * p * miu_hf_p * abs(E_v[k])
    Gava = Gava_n + Gava_p
    return Gava


def dark_current(u, xp, xn, mesh, E, interface, area, Hurkx, D_InP,
                 D_InGaAsP, D_InGaAs, c, Q, avalanche, total_doping, electron_conc, hole_conc):
    # u[0, :] = InP
    # u[1, :] = InGaAsP
    # u[2, :] = InGaAs
    # u[3, :] = InP
    # u = [ InP / InGaAsP / InGaAs / .... ]
    # interface [ 1 / 2 / 3 / .... ]
    E = abs(E)
    J_data = np.zeros((len(xp), len(mesh)))
    a = [c * 1.72e20, c * 1.817e20, c * 1.785e20]
    b = [1.82e6, 10.14e6, 6.188e6]
    p = [2, 2, 2]
    q = 1.6e-19  # [C]
    current_density = np.zeros(len(xn))
    if Hurkx == 'Hurkx':
        for i in range(len(xn)):
            k1, xk1 = position_converter(xp[i], mesh, 0)
            k2, xk2 = position_converter(xn[i], mesh, 0)
            it1, xit1 = position_converter(interface[0], mesh, 0)
            it2, xit2 = position_converter(interface[1], mesh, 0)
            it3, xit3 = position_converter(interface[2], mesh, 0)
            #print(k1, k2, it1, it2)
            if xn[i] < interface[0]:
                temp = 0
                for j in range(k1, k2):
                    if avalanche == 'avalanche':
                        G_ava = avalanche_generation(E[i, j], mesh[j], mesh, interface,
                                                     total_doping, E[i, :], electron_conc[i, j], hole_conc[i, j])
                    else:
                        G_ava = 0
                    if E[i, j] == 0:  # u[0]
                        temp = temp + u[j] * (mesh[j] - mesh[j-1])
                    else:     # u[0]
                        temp = temp + (u[j] + G_ava + D_InP * a[0] * E[i, j] ** p[0] * np.exp(-b[0] / E[i, j])) * (mesh[j] - mesh[j-1])
                    J_data[i, j] = q * temp * 1e-4
                current_density[i] = q * temp * 1e-4
                for j in range(0, k2):
                    J_data[i, j] = J_data[i, j] - q * temp * 1e-4
            elif interface[0] <= xn[i] < interface[1]:
                temp = 0
                for j in range(k1, it1):
                    if avalanche == 'avalanche':
                        G_ava = avalanche_generation(E[i, j], mesh[j], mesh, interface,
                                                     total_doping, E[i, :], electron_conc[i, j], hole_conc[i, j])
                    else:
                        G_ava = 0
                    if E[i, j] == 0:  # u[0]
                        temp = temp + u[j] * (mesh[j] - mesh[j - 1])
                    else:  # u[0]
                        temp = temp + (u[j] + G_ava + D_InP * a[0] * E[i, j] ** p[0] * np.exp(-b[0] / E[i, j])) * (mesh[j] - mesh[j-1])
                    J_data[i, j] = q * temp * 1e-4
                for j in range(it1, k2):
                    if avalanche == 'avalanche':
                        G_ava = avalanche_generation(E[i, j], mesh[j], mesh, interface,
                                                     total_doping, E[i, :], electron_conc[i, j], hole_conc[i, j])
                    else:
                        G_ava = 0
                    if E[i, j] == 0:  # u[1]
                        temp = temp + u[j] * (mesh[j] - mesh[j - 1])
                    else:  # u[1]
                        temp = temp + (u[j] + G_ava + D_InGaAsP * a[1] * E[i, j] ** p[1] * np.exp(-b[1] / E[i, j])) * (mesh[j] - mesh[j-1])
                    J_data[i, j] = q * temp * 1e-4
                current_density[i] = q * temp * 1e-4
                for j in range(0, k2):
                    J_data[i, j] = J_data[i, j] - q * temp * 1e-4
            elif interface[1] <= xn[i] < interface[2]:
                temp = 0
                for j in range(k1, it1):
                    if avalanche == 'avalanche':
                        G_ava = avalanche_generation(E[i, j], mesh[j], mesh, interface,
                                                     total_doping, E[i, :], electron_conc[i, j], hole_conc[i, j])
                    else:
                        G_ava = 0
                    if E[i, j] == 0:  # u[0]
                        temp = temp + u[j] * (mesh[j] - mesh[j - 1])
                    else:  # u[0]
                        temp = temp + (u[j] + G_ava + D_InP * a[0] * E[i, j] ** p[0] * np.exp(-b[0] / E[i, j])) * (mesh[j] - mesh[j-1])
                    J_data[i, j] = q * temp * 1e-4
                for j in range(it1, it2):
                    if avalanche == 'avalanche':
                        G_ava = avalanche_generation(E[i, j], mesh[j], mesh, interface,
                                                     total_doping, E[i, :], electron_conc[i, j], hole_conc[i, j])
                    else:
                        G_ava = 0
                    if E[i, j] == 0:  # u[1]
                        temp = temp + u[j] * (mesh[j] - mesh[j - 1])
                    else:  # u[1]
                        temp = temp + (u[j] + G_ava + D_InGaAsP * a[1] * E[i, j] ** p[1] * np.exp(-b[1] / E[i, j])) * (mesh[j] - mesh[j-1])
                    J_data[i, j] = q * temp * 1e-4
                for j in range(it2, k2):
                    if avalanche == 'avalanche':
                        G_ava = avalanche_generation(E[i, j], mesh[j], mesh, interface,
                                                     total_doping, E[i, :], electron_conc[i, j], hole_conc[i, j])
                    else:
                        G_ava = 0
                    if E[i, j] == 0:  # u[2]
                        temp = temp + u[j] * (mesh[j] - mesh[j - 1])
                    else:  # u[2]
                        temp = temp + (u[j] + G_ava + D_InGaAs * a[2] * (Q * E[i, j]) ** p[2] * np.exp(-b[2] / (Q * E[i, j]))) * (mesh[j] - mesh[j-1])
                    J_data[i, j] = q * temp * 1e-4
                current_density[i] = q * temp * 1e-4
                for j in range(0, k2):
                    J_data[i, j] = J_data[i, j] - q * temp * 1e-4
            else:
                temp = 0
                for j in range(k1, it1):
                    if avalanche == 'avalanche':
                        G_ava = avalanche_generation(E[i, j], mesh[j], mesh, interface,
                                                     total_doping, E[i, :], electron_conc[i, j], hole_conc[i, j])
                    else:
                        G_ava = 0
                    if E[i, j] == 0:  # u[0]
                        temp = temp + u[j] * (mesh[j] - mesh[j - 1])
                    else:  # u[0]
                        temp = temp + (u[j] + G_ava + D_InP * a[0] * E[i, j] ** p[0] * np.exp(-b[0] / E[i, j])) * (
                                    mesh[j] - mesh[j - 1])
                    J_data[i, j] = q * temp * 1e-4
                for j in range(it1, it2):
                    if avalanche == 'avalanche':
                        G_ava = avalanche_generation(E[i, j], mesh[j], mesh, interface,
                                                     total_doping, E[i, :], electron_conc[i, j], hole_conc[i, j])
                    else:
                        G_ava = 0
                    if E[i, j] == 0:  # u[1]
                        temp = temp + u[j] * (mesh[j] - mesh[j - 1])
                    else:  # u[1]
                        temp = temp + (u[j] + G_ava + D_InGaAsP * a[1] * E[i, j] ** p[1] * np.exp(-b[1] / E[i, j])) * (
                                    mesh[j] - mesh[j - 1])
                    J_data[i, j] = q * temp * 1e-4
                for j in range(it2, it3):
                    if avalanche == 'avalanche':
                        G_ava = avalanche_generation(E[i, j], mesh[j], mesh, interface,
                                                     total_doping, E[i, :], electron_conc[i, j], hole_conc[i, j])
                    else:
                        G_ava = 0
                    if E[i, j] == 0:  # u[2]
                        temp = temp + u[j] + (mesh[j] - mesh[j - 1])
                    else:  # u[2]
                        temp = temp + (u[j] + G_ava + D_InGaAs * a[2] * (Q * E[i, j]) ** p[2] *
                                       np.exp(-b[2] / (Q * E[i, j]))) * (mesh[j] - mesh[j - 1])
                    J_data[i, j] = q * temp * 1e-4
                for j in range(it3, k2):
                    if avalanche == 'avalanche':
                        G_ava = avalanche_generation(E[i, j], mesh[j], mesh, interface,
                                                     total_doping, E[i, :], electron_conc[i, j], hole_conc[i, j])
                    else:
                        G_ava = 0
                    if E[i, j] == 0:  # u[0]
                        temp = temp + u[j] * (mesh[j] - mesh[j - 1])
                    else:  # u[0]
                        temp = temp + (u[j] + G_ava + D_InP * a[0] * E[i, j] ** p[0] *
                                       np.exp(-b[0] / E[i, j])) * (mesh[j] - mesh[j - 1])
                    J_data[i, j] = q * temp * 1e-4
                current_density[i] = q * temp * 1e-4
                for j in range(0, k2):
                    J_data[i, j] = J_data[i, j] - q * temp * 1e-4
        return current_density * area, J_data
    else:  # not yet edit its SRH generation rate
        for i in range(len(xn)):
            if xn[i] < interface[0]:
                current_density[i] = q * u[0] * (xn[i] - xp[i]) * 1e-4
            elif interface[0] <= xn[i] < interface[1]:
                current_density[i] = (q * u[0] * (interface[0] - xp[i]) + q * u[1] * (xn[i] - interface[0])) * 1e-4
            elif interface[1] <= xn[i] < interface[2]:
                current_density[i] = (q * u[0] * (interface[0] - xp[i])) * 1e-4 + \
                    (q * u[1] * (interface[1] - interface[0])) * 1e-4 + \
                    (q * u[2] * (xn[i] - interface[0])) * 1e-4
            else:
                current_density[i] = (q * u[0] * (interface[0] - xp[i])) * 1e-4 + \
                                     (q * u[1] * (interface[1] - interface[0])) * 1e-4 + \
                                     (q * u[2] * (interface[2] - interface[0])) * 1e-4 + \
                                     (q * u[0] * (xn[i] - interface[2])) * 1e-4
        return current_density * area


def potential(V_number, mesh, E):
    V = np.zeros((V_number, len(mesh)))
    for i in range(V_number):
        for j in range(1, len(mesh)):
            V[i, j] = V[i, j-1] - E[i, j-1] * (mesh[j] - mesh[j-1]) * 1e-4
        for j in range(len(mesh)):
            V[i, j] = V[i, j] - V[i, len(mesh)-1]
    return V


def remesh(mesh_initial, mesh_final, f):
    g = np.zeros(len(mesh_final))
    for i in range(len(g)):
        k = position_converter(mesh_final[i], mesh_initial, 1)
        if k < len(mesh_initial) - 1:
            g[i] = (f[k] + f[k+1]) / 2
        else:
            g[i] = f[k]
    return g


def qfl_auxiliary_function(mesh, position, J, Ei, interface, total_doping, E, f0, carrier):  # [-]
    k, xk = position_converter(position, mesh, 0)
    kT = 0.0259  # [eV]
    if carrier == 'h':
        if k == len(mesh) - 1:
            dx = (mesh[k] - mesh[k - 1]) * 1e-4  # [cm]
        else:
            dx = (mesh[k + 1] - mesh[k]) * 1e-4  # [cm]
        miu_e, miu_h, miu_hf_e, miu_hf_h = mobility(mesh[k], interface, mesh, total_doping, E)
        f = f0 + ((J[k] * mp.exp(- Ei[k] / kT)) / (miu_hf_h * ni(mesh[k], interface))) * dx / kT
        #print(((J[k] * mp.exp(- Ei[k] / kT)) / (miu_hf_h * ni(mesh[k], interface))) * dx / kT)
    # The following codes are useless
    elif carrier == 'e':
        if k == len(mesh) - 1:
            dx = (mesh[k] - mesh[k - 1]) * 1e-4  # [cm]
        else:
            dx = (mesh[k + 1] - mesh[k]) * 1e-4  # [cm]
        miu_e, miu_h, miu_hf_e, miu_hf_h = mobility(mesh[k], interface, mesh, total_doping, E)
        f = f0 + ((J[k] * mp.exp(Ei[k] / kT)) / (miu_hf_e * ni(mesh[k], interface))) * dx / kT
        # print('df: ', ((J[k] * mp.exp(Ei[k] / kT)) / (miu_hf_e * ni(mesh[k], interface))) * dx / kT)
    elif carrier == 'total':
        # for electron 2nd term
        f = 0
        df_plot = np.zeros(len(mesh))
        for i in range(len(mesh)):
            if i == len(mesh) - 1:
                dx = (mesh[i] - mesh[i - 1]) * 1e-4  # [cm]
            else:
                dx = (mesh[i + 1] - mesh[i]) * 1e-4  # [cm]
            miu_e, miu_h, miu_hf_e, miu_hf_h = mobility(mesh[i], interface, mesh, total_doping, E)
            f = f + ((J[i] * mp.exp(Ei[i] / kT)) / (miu_hf_e * ni(mesh[i], interface))) * dx / kT
            df_plot[i] = ((J[i] * mp.exp(Ei[i] / kT)) / (miu_hf_e * ni(mesh[i], interface))) * dx / kT
            print('i: ', i, 'df[i]: ', ((J[i] * mp.exp(Ei[i] / kT)) / (miu_hf_e * ni(mesh[i], interface))) * dx / kT)
            print('i: ', i, 'mesh[i]: ', mesh[i], 'ft % 10.f' % f)
    return f


def hole_QFL(V_number, mesh, Jp, Ei, Doping, interface, total_doping, E):  # [eV]
    p_qfl = np.zeros((V_number, len(mesh)))
    kT = 0.0259  # [eV]
    for i in range(V_number):
        f0 = 0
        for j in range(len(mesh)):
            f0 = qfl_auxiliary_function(mesh, mesh[j], Jp[i, :], Ei[i, :], interface, total_doping, E[i, :], f0, 'h')
            p_qfl[i, j] = \
                - kT * \
                mp.ln(abs(Doping[0])/ni(mesh[0], interface) *
                       mp.exp(-Ei[i, 0]/kT) - f0)
            #print('i: ', i, 'j: ', j, 'mesh[j]: ', mesh[j], 'f0: ', f0, 'p_QFL: ', p_qfl[i, j])
            #if -4.27 < mesh[j] < -4.25:
            #    slope = (p_qfl[i, j] - p_qfl[i, j-1]) / (mesh[j] - mesh[j - 1])
            #    print('slope: ', slope)
    return p_qfl


def electron_QFL(V_number, mesh, Jn, Ei, Doping, interface, total_doping, E):  # [eV]
    n_qfl = np.zeros((V_number, len(mesh)))
    kT = 0.0259  # [eV]
    for i in range(V_number):
        f0 = 0
        for j in range(len(mesh)):
            k = position_converter(mesh[j], mesh, 1)
            j2 = len(mesh) - 1 - k
            if k == len(mesh) - 1:
                dx = (mesh[1] - mesh[0]) * 1e-4  # [cm]
            else:
                dx = - (mesh[j2] - mesh[j2 - 1]) * 1e-4  # [cm] the minus sign comes from int_b^x, b < x
            miu_e, miu_h, miu_hf_e, miu_hf_h = mobility(mesh[j2], interface, mesh, total_doping, E[i, :])
            f0 = f0 + ((Jn[i, j2] * mp.exp(Ei[i, j2] / kT)) /
                       (miu_hf_e * ni(mesh[j2], interface))) * dx / kT
            s = (abs(Doping[1710])/ni(mesh[0], interface) * mp.exp(Ei[i, 1710]/kT)) * 1e50 + f0 * 1e50
            n_qfl[i, j2] = kT * (mp.ln(s) - 115.12925464970229)
            # print('--- %.5s seconds ---' % (time.time() - start_time), 'i: ', i, 'j: ', j2, 'N-QFL: ', n_qfl[i, j2])
    return n_qfl


def V_selection(V_meshline, V_meshlist, xp, xn, V, E):
    V_trans = np.zeros(len(V))
    xp_trans = np.zeros(len(V))
    xn_trans = np.zeros(len(V))
    E_trans = np.zeros((len(V), len(E[1])))
    for i in range(len(V)):
        V_trans[i] = V[len(V) - 1 - i]
        xp_trans[i] = xp[len(V) - 1 - i]
        xn_trans[i] = xn[len(V) - 1 - i]
        E_trans[i] = E[len(V) - 1 - i]
    N = 0
    for i in range(len(V_meshlist)):
        N = N + V_meshlist[i]
    V_mesh = np.zeros(N)
    xp_mesh = np.zeros(N)
    xn_mesh = np.zeros(N)
    E_mesh = np.zeros((N, len(E[1])))
    for i in range(len(V_meshline)):
        if i == 0:
            for j in range(V_meshlist[0]):
                V_mesh[j] = V_trans[0] + j * (V_meshline[0] - V_trans[0]) / V_meshlist[0]
        else:
            sum_mesh = 0
            for j in range(i):
                sum_mesh = sum_mesh + V_meshlist[j]
            for j in range(V_meshlist[i]):
                V_mesh[sum_mesh + j] = V_mesh[sum_mesh - 1] + (j + 1) * (V_meshline[i] - V_meshline[i-1]) / V_meshlist[i]
    for i in range(N):
        for j in range(len(V_trans)-1):
            if V_mesh[i] == V_trans[j]:
                xp_mesh[i] = xp_trans[j]
                xn_mesh[i] = xn_trans[j]
                E_mesh[i] = E_trans[j]
                break
            elif (V_mesh[i] - V_trans[j]) * (V_mesh[i] - V_trans[j+1]) < 0:
                if abs(V_mesh[i] - V_trans[j]) < abs(V_mesh[i] - V_trans[j+1]):
                    xp_mesh[i] = xp_trans[j]
                    xn_mesh[i] = xn_trans[j]
                    E_mesh[i] = E_trans[j]
                    break
                else:
                    xp_mesh[i] = xp_trans[j+1]
                    xn_mesh[i] = xn_trans[j+1]
                    E_mesh[i] = E_trans[j+1]
                    break
    xp_mesh[N - 1] = xp_trans[len(V_trans) - 1]
    xn_mesh[N - 1] = xn_trans[len(V_trans) - 1]
    E_mesh[N - 1] = E_trans[len(V_trans) - 1]
    return N, V_mesh, xp_mesh, xn_mesh, E_mesh

'''
def electron_quasi_fermi_level(V_number, mesh, Jn, Ei, Doping, interface, total_doping, E):  # [eV]
    n_qfl = np.zeros((V_number, len(mesh)))
    kT = 0.0259  # [eV]

    # plot
    df0 = np.zeros(len(mesh))
    # plot

    for i in range(25, 26):
        f0 = 0
        f_t = qfl_auxiliary_function(mesh, mesh[0], Jn[i, :], Ei[i, :], interface, total_doping, E[i, :], f0, 'total')
        for j in range(len(mesh)):
            k = position_converter(mesh[j], mesh, 1)
            if k == len(mesh) - 1:
                dx = (mesh[k] - mesh[k - 1]) * 1e-4  # [cm]
            else:
                dx = (mesh[k + 1] - mesh[k]) * 1e-4  # [cm]
            miu_e, miu_h, miu_hf_e, miu_hf_h = mobility(mesh[k], interface, mesh, total_doping, E[i, :])
            f0 = f0 + ((Jn[i, k] * mp.exp(Ei[i, k] / kT)) / (miu_hf_e * ni(mesh[k], interface))) * dx / kT
            # plot
            df0[j] = ((Jn[i, k] * mp.exp(Ei[i, k] / kT)) / (miu_hf_e * ni(mesh[k], interface))) * dx / kT
            # plot
            print('df0: ', ((Jn[i, k] * mp.exp(Ei[i, k] / kT)) / (miu_hf_e * ni(mesh[k], interface))) * dx / kT)
            s = (abs(Doping[1710])/ni(mesh[0], interface) * mp.exp(Ei[i, 1710]/kT)) * 1e50 + (- f_t + f0) * 1e50
            n_qfl[i, j] = kT * (mp.ln(s) - 115.12925464970229)
            print('--- %.5s seconds ---' % (time.time() - start_time), 'i: ', i, 'j: ', j, 'N-QFL: ', n_qfl[i, j])
            if -4 < n_qfl[i, j] < 9.86:
                print('ft: ', f_t)
                print('f0: ', f0)
                print('ft: % 10.f' % f_t)
                print('f0: % 10.f' % f0)
                print('f0-ft:', f0 - f_t)
                print('f0-ft: % 10.f' % (f0 - f_t))
            if n_qfl[i, j] < -4:
                print('ft: ', f_t)
                print('f0: ', f0)
                print('ft: % 10.f' % f_t)
                print('f0: % 10.f' % f0)
                print('f0-ft:', f0 - f_t)
                print('f0-ft: % 10.f' % (f0 - f_t))
                break
        # plot
        figure(num=10, figsize=(20, 9.3), dpi=80, facecolor='w', edgecolor='k')
        plt.plot(mesh, abs(df0), label='df0')
        plt.yscale('log')
        plt.legend(loc='upper right')
        plt.show()
        # plot
    return n_qfl



def electron_quasi_fermi_level(V_number, mesh, Jn, Ei, Doping, interface, total_doping, E):  # [eV]
    n_qfl = np.zeros((V_number, len(mesh)))
    kT = 0.0259  # [eV]
    for i in range(25, 26):
        f0 = 0
        df0_plot = np.zeros(len(mesh))
        f_t = qfl_auxiliary_function(mesh, mesh[0], Jn[i, :], Ei[i, :], interface, total_doping, E[i, :], f0,
                                         'total')
        for j in range(len(mesh)):
            temp = qfl_auxiliary_function(mesh, mesh[j], Jn[i, :], Ei[i, :], interface, total_doping, E[i, :], f0, 'electron')
            df0_plot[j] = temp - f0
            f0 = temp
            s = (abs(Doping[1710])/ni(mesh[0], interface) * mp.exp(Ei[i, 1710]/kT)) * 1e50 + (- f_t + f0) * 1e50
            n_qfl[i, j] = kT * (mp.ln(s) - 115.12925464970229)
            print('--- %.5s seconds ---' % (time.time() - start_time), 'i: ', i, 'j: ', j, 'N-QFL: ', n_qfl[i, j])
            if i >= 14:
                print('ND/ni * exp(Ei/kT): ', abs(Doping[1710]) / ni(mesh[0], interface) * mp.exp(Ei[i, 1710] / kT))
                print('f0: ', f0)
                print('ft: ', f_t)
                print('f0: % 10.f' % f0)
                print('ft: % 10.f' % f_t)
                print('f0 - ft: ', - f_t + f0)
                print('Total: ', abs(Doping[1710]) / ni(mesh[0], interface) * mp.exp(Ei[i, 1710] / kT) - f_t + f0)
                print('s: ', s, 'kT * mp.ln(s): ', kT * mp.ln(s))
                if n_qfl[i, j] < -400:
                    print('i: ', i, 'j: ', j, 'mesh[j]', mesh[j], 'Jn[i, j-1]', Jn[i, j-1], 'Jn[i, j]', Jn[i, j])
                    break
        print('--- %.5s seconds ---' % (time.time() - start_time), 'i: ', i, 'n-QFL[i, 1400]: ', n_qfl[i, 1400])
        figure(num=10, figsize=(20, 9.3), dpi=80, facecolor='w', edgecolor='k')
        plt.plot(mesh, abs(df0_plot), label='df0')
        plt.yscale('log')
        plt.show()
    return n_qfl
'''