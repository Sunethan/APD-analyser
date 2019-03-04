import numpy as np
import csv
import mpmath as mp
import math
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
import matplotlib.pylab as pylab
import time
start_time = time.time()

# Physics constants
e = 1.6e-19  # [Coul/e]
eps_InP = 12.5 * 8.85e-14  # [F/cm]
eps_InGaAs = 13.9 * 8.85e-14  # [F/cm] In 0.53 Ga 0.47 P
eps_InGaAsP = 13.436 * 8.85e-14  # [F/cm] Approximated by In 0.53 Ga 0.47 As 0.65 P 0.35


def ConvertToInt(list):
    for i in range(len(list)):
        temp = list[i]
        list[i] = int(temp)
    return list


def ConvertToFloat(list):
    for i in range(len(list)):
        temp = list[i]
        list[i] = float(temp)
    return list


def find(x, func, x_value, form):
    if x[0] < x[1]:
        curve = 'increase'
    else:
        curve = 'decrease'
    for i in range(len(x) - 1):
        #print(x[i], x_value, x[i + 1])
        if x[i] == x_value:
            return func[i]
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
        if (x[i] > x_value > x[i + 1]) and (curve == 'decrease'):
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
    x_dydx = np.zeros(len(y) - 1)
    dydx = np.zeros(len(y) - 1)
    for i in range(len(y) - 1):
        x_dydx[i] = x[i]
        dydx[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i])
    return x_dydx, dydx


def center_smooth(func, k):
    ysmooth = np.zeros(len(func))
    for i in range(k, len(func)):
        if k <= i <= len(func) - k - 1:
            for j in range(i - k, i + k + 1):
                ysmooth[i] = ysmooth[i] + func[j] / (2 * k + 1)
        elif i < k:
            for j in range(0, i + k + 1):
                ysmooth[i] = ysmooth[i] + func[j] / (i + k + 1)
        else:
            for j in range(i - k, len(func)):
                ysmooth[i] = ysmooth[i] + func[j] / (len(func) - i + k)
    return ysmooth


def parse_csv(csv_path, key, Cname, X, floatKey):
    result = dict()
    Curve_name = np.empty([0])
    with open(csv_path, newline='', encoding='utf-8-sig') as csv_file:
        rows = csv.DictReader(csv_file)
        fieldnames = rows.fieldnames
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
                    if row[fieldname] == '-' or row[fieldname] == '':
                        continue
                    curve_name, x_y = fieldname.rsplit(' ', 1)
                    if floatKey == 'Yes':
                        result[curve_name][x_y].append(float(row[fieldname]))
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
                    if row[fieldname] == '-' or row[fieldname] == '':
                        continue
                    if floatKey == 'Yes':
                        result[fieldname].append(float(row[fieldname]))
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
def alpha(electric_field, position, interface):
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
            if electric_field < 3.8e5:
                return 1.12e7 * np.exp(-3.11e6 / abs(electric_field))
            else:
                return 2.93e6 * np.exp(-2.64e6 / abs(electric_field))
        else:
            return 0


def beta(electric_field, position, interface):
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
            if electric_field < 3.8e5:
                return 4.79e6 * np.exp(-2.55e6 / abs(electric_field))
            else:
                return 1.62e6 * np.exp(-2.11e6 / abs(electric_field))
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
    ix = np.zeros(len(line), dtype=int)  # temporary indices of position
    for j in range(len(n)):
        total_mesh = total_mesh + n[j]
    x = np.zeros(total_mesh)  # new mesh array
    doping = np.zeros(total_mesh)  # new doping profile in terms of new mesh array
    for position in range(0, total_mesh):
        # InP top bulk region  (p+ InP)
        if position <= n[0]:
            #print(n)
            #print('x_raw[0]: ', x_raw[0])
            #print('position: ', position)
            #print('line[0]: ', line[0])
            #print('n[0]: ', n[0])
            x[position] = x_raw[0] + position * (line[0] - x_raw[0]) / n[0]
            ix[0] = position
        # Multiplication layer (p+ InP)
        elif n[0] < position < n[0] + n[1]:
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
        elif n[0] + n[1] + n[2] + n[3] + n[4] + n[5] + n[6] \
                <= position < n[0] + n[1] + n[2] + n[3] + n[4] + n[5] + n[6] + n[7]:
            x[position] = x[ix[6]] + (position - ix[6]) * (line[7] - x[ix[6]]) / n[7]
            ix[7] = position
        else:  # Substrate (n+ InP)
            x[position] = x[ix[7]] + (position - ix[7]) * (x_raw[len(x_raw) - 1] - x[ix[7]]) / n[8]
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
        if x_raw[i + 1] != x_raw[i]:
            m1 = (doping_raw[i + 1] - doping_raw[i]) / (x_raw[i + 1] - x_raw[i])
        else:
            m1 = 0
        '''
        print('-----------------')
        print('doping: ', doping_raw[i - 1], doping_raw[i], doping_raw[i + 1])
        print('position: ', x_raw[i - 1], x_raw[i], x_raw[i + 1])
        print('slope: ', m1, m2)
        print('-----------------')
        '''
        # Find the PN junction
        if (doping_raw[i + 1] * doping_raw[i] < 0 and
                i_count == 0):
            xj = x_raw[i+1]
            i_count = 1
        # Find the interface between multiplication and charge layer
        if conc[0] <= doping_raw[i] < conc[1] and i_count == 1 and m1 > 1e17 and doping_raw[i] > 1.1 * conc[0]:
            '''
            print(doping_raw[i - 1], doping_raw[i], doping_raw[i + 1])
            print(x_raw[i - 1], x_raw[i], x_raw[i + 1])
            print(m1, m2)
            '''
            xj_mc = x_raw[i - 1]
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


def multiplication_factor(x, ef, interface):
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
            side[0] = - (alpha(ef[k + 1], x[k + 1], interface) - beta(ef[k + 1], x[k + 1], interface))
            side[1] = - (alpha(ef[k], x[k], interface) - beta(ef[k], x[k], interface))
            temp_un = temp_un + (side[0] + side[1]) * (x[k + 1] - x[k]) / 2 * 1e-4
        side[0] = 0
        side[1] = 0
        for k in range(i, len(x) - 1):
            side[0] = (alpha(ef[k + 1], x[k + 1], interface) - beta(ef[k + 1], x[k + 1], interface))
            side[1] = (alpha(ef[k], x[k], interface) - beta(ef[k], x[k], interface))
            temp_up = temp_up + (side[0] + side[1]) * (x[k + 1] - x[k]) / 2 * 1e-4
        un[i] = np.exp(temp_un)
        up[i] = np.exp(temp_up)
        temp_un = 0
        temp_up = 0
        side[0] = 0
        side[1] = 0
    for i in range(0, len(x) - 1):
        side1n = alpha(ef[i], x[i], interface) * un[i]
        side2n = alpha(ef[i + 1], x[i + 1], interface) * un[i + 1]
        temp_mn = temp_mn + (side1n + side2n) * (x[i + 1] - x[i]) * 1e-4 / 2
        side1p = beta(ef[i], x[i], interface) * up[i]
        side2p = beta(ef[i + 1], x[i + 1], interface) * up[i + 1]
        temp_mp = temp_mp + (side1p + side2p) * (x[i + 1] - x[i]) * 1e-4 / 2
    return temp_mn, temp_mp


def vb_solver(trial_number, step, error, x, e_goal, interface, s, v_number, doping, form, rj, xj, v0_goal):
    gate = 0
    v = v0_goal
    mn, mp = multiplication_factor(x, e_goal, interface)
    error_n = abs(mn - 1)
    error_p = abs(mp - 1)
    print('[1]', '  V: %.2f' % v0_goal, '  Mn: %.3f' % mn, '  Mp: %.3f' % mp,
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
            print('[%d]' % (i + 1), '  V: %.2f' % v[0], '  Mn: %.3f' % mn, '  Mp: %.3f' % mp,
                  '  Error: %.3f' % max(error_n, error_p), '--- %.5s seconds ---  (+s: %s)' % ((time.time() - start_time), s))
            gate = 1
            if i == end:
                error_n = abs(mn - 1)
                error_p = abs(mp - 1)
                print('Iterations larger than maximum %.3f' % v)
                print('[%d]' % (i + 1), '  V: %.2f' % v[0], '  Mn: %.3f' % mn, '  Mp: %.3f' % mp,
                      '  Error: %.3f' % max(error_n, error_p), '--- %.5s seconds ---' % (time.time() - start_time))
                return v, mn, mp, i, xp1, xn1
        elif mn < 1 - error or mp < 1 - error:
            if gate > 0:
                step = step / 2
            s = s - step
            xp1, xn1, v, ef, xp_goal, xn_goal, v_goal, e_goal = \
                depletion_boundary(v_number, s, s, x, doping, interface, 0, form, rj, xj)
            mn, mp = multiplication_factor(x, ef[0, :], interface)
            error_n = abs(mn - 1)
            error_p = abs(mp - 1)
            print('[%d]' % (i + 1), '  V: %.2f' % v[0], '  Mn: %.3f' % mn, '  Mp: %.3f' % mp,
                  '  Error: %.3f' % max(error_n, error_p), '--- %.5s seconds ---  (-s: %s)' % ((time.time() - start_time), s))
            gate = -1
            if i == end:
                error_n = abs(mn - 1)
                error_p = abs(mp - 1)
                print('Iterations larger than maximum %.3f' % v)
                print('[%d]' % (i + 1), '  V: %.2f' % v[0], '  Mn: %.3f' % mn, '  Mp: %.3f' % mp,
                      '  Error: %.3f' % max(error_n, error_p), '--- %.5s seconds ---' % (time.time() - start_time))
                return v, mn, mp, i, xp1, xn1
        else:
            return v, mn, mp, i-1, xp1, xn1


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
    #xp_mesh[N - 1] = xp_trans[len(V_trans) - 1]
    #xn_mesh[N - 1] = xn_trans[len(V_trans) - 1]
    #E_mesh[N - 1] = E_trans[len(V_trans) - 1]

    NonZeroCount = 0
    for i in range(len(xp_mesh)):
        if xp_mesh[i] != 0:
            NonZeroCount += 1
    Effective_xp_mesh = np.zeros(NonZeroCount)
    Effective_xn_mesh = np.zeros(NonZeroCount)
    Effective_E_mesh = np.zeros((NonZeroCount, len(E[1])))
    Effective_V_mesh = np.zeros(NonZeroCount)
    j = 0
    for i in range(len(xp_mesh)):
        if xp_mesh[i] != 0:
            Effective_V_mesh[j] = V_mesh[i]
            Effective_xp_mesh[j] = xp_mesh[i]
            Effective_xn_mesh[j] = xn_mesh[i]
            Effective_E_mesh[j] = E_mesh[i]
            j += 1
    return N, Effective_V_mesh, Effective_xp_mesh, Effective_xn_mesh, Effective_E_mesh