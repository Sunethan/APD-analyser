import numpy as np
import csv
import matplotlib
matplotlib.use('TkAgg')
import inspect
import os.path
import json
import time
start_time = time.time()


def array_division(x1, y1, x2, y2):
    x = list(set(list(set(x1)) + list(set(x2))))
    x.sort(reverse=True)
    y1_prime = interpolation(x1, x, y1)
    y2_prime = interpolation(x2, x, y2)
    ans = y2_prime / y1_prime
    return x, ans


def CheckOverlap(numpy_x, numpy_y):
    x = numpy_x
    y = numpy_y
    length = len(x)
    i = 0
    while i < length - 1:
        if i < len(x) - 1 and x[i] == x[i + 1]:
            temp = (y[i] + y[i + 1]) / 2
            print('發現有重複點：x=%s, y[%s]=%s, y[%s]=%s, y(avg)=%s' % (x[i], i, y[i], i + 1, y[i + 1], temp))
            y[i] = temp
            x = np.delete(x, i + 1)
            y = np.delete(y, i + 1)
            i = i - 1
        i = i + 1
    '''
    for i in range(length-1):
        if i < len(x) - 1 and x[i] == x[i + 1]:
            temp = (y[i] + y[i + 1]) / 2
            print('發現有重複點：x=%s, y[%s]=%s, y[%s]=%s, y(avg)=%s' % (x[i], i, y[i], i + 1, y[i + 1], temp))
            y[i] = temp
            x = np.delete(x, i + 1)
            y = np.delete(y, i + 1)
    '''
    print('原本長度：%s，後來長度：%s' % (length, len(x)))
    return np.asarray(x), np.asarray(y)


def point_interpolation(x1, y1, x2, y2, goal_axis, goal_value, form):
    """
    目前還沒做好 goal_axis == 'x' 情況下的功能。
    point_interpolation()：由 (x1,y1)、(x2,y2) 與 x/y(goal) ，內插求 y/x(goal)
    find()：由 x(array)、y(array)、x(goal) 求 y(goal)。
    :param x1: 已知的 (x1,y1)
    :param y1: 已知的 (x1,y1)
    :param x2: 已知的 (x2,y2)
    :param y2: 已知的 (x2,y2)
    :param goal_axis: 欲求的已知值類型，例如已知 y(goal) 求 x(goal），則輸入'y'。
    :param goal_value: 欲求的已知值類型，例如已知 y(goal)=102 求 x(goal），則輸入102。
    :param form: 看是用 log 還是 linear 方式近似
    :return: 回傳欲求的目標的值
    """
    if goal_axis == 'x':
        if form == 'linear':
            yg = y1 + (y2 - y1) / (x2 - x1) * (goal_value - x1)
            return yg
        elif form == 'log':
            if y1 * y2 > 0:
                Y1 = np.log10(abs(y1))
                Y2 = np.log10(abs(y2))
                yg = 10 ** (Y1 + (Y2 - Y1) / (x2 - x1) * (goal_value - x1))
                return yg
            else:
                raise BaseException("y1=%.2e and y2=%.2e, they're not the same sign." % (y1, y2))
        else:
            raise BaseException("Wrong form: %s" % form)
    elif goal_axis == 'y':
        if form == 'linear':
            xg = x2 - (goal_value - y2) / (y1 - y2) * (x2 - x1)
            return xg
        elif form == 'log':
            if y1 * y2 > 0:
                Y1 = np.log10(abs(y1))
                Y2 = np.log10(abs(y2))
                Yg = np.log10(abs(goal_value))
                xg = x2 - (Yg - Y2) / (Y1 - Y2) * (x2 - x1)
                return xg
            else:
                raise BaseException("y1=%.2e and y2=%.2e, they're not the same sign." % (y1, y2))
        else:
            raise BaseException("Wrong form: %s" % form)


def Remesh(x, y, dx):
    x_new = np.arange(x[0], x[-1], dx)
    y_new = np.asarray([find(x, y, element, 'linear') for element in x_new])
    return x_new, y_new


def positive_array_filter(array):
    new_array = []
    for element in array:
        if element > 0:
            new_array.append(element)
        else:
            new_array.append(0)
    return np.asarray(new_array)


def find_intersection(func1, func2, x1, x2, tolerance):
    if x1 >= x2:
        raise BaseException("x1 should be smaller than x2.")

    def func(x):
        if callable(func2):
            return func1(x) - func2(x)
        else:
            return func1(x) - func2

    x1_final = x1
    x2_final = x2
    error = abs(x2 - x1)
    count = 0
    while error > tolerance and count < 100:
        #print('x1:%.3f, x2:%.3f, error(abs(x2-x1)):%.3e, func(x1):%.3e, df:%.3e' % (x1, x2, error, abs(func(x1)), abs(func((x1 + x2) / 2)) - abs(func(x1))))
        if func(x1) * func((x1 + x2) / 2) < 0:
            x2 = (x1 + x2) / 2
            if abs(func((x1 + x2) / 2)) < abs(func(x1)):
                x2_final = x2
        else:
            x1 = (x1 + x2) / 2
            if abs(func((x1 + x2) / 2)) < abs(func(x1)):
                x1_final = x1
        error = abs(x2 - x1)
        count += 1
    if abs(x1_final - x2_final) < error:
        print('[%s] error(abs(x2-x1)): %.3e, x: %.3e, df=%.3e, func1(x)=%.3e, func2(x)=%.3e' % (
            count, error, (x1_final + x2_final) / 2, func1((x1_final + x2_final) / 2) - func2((x1_final + x2_final) / 2), func1((x1_final + x2_final) / 2), func2((x1_final + x2_final) / 2)))
        return (x1_final + x2_final) / 2
    else:
        if callable(func2):
            print('[%s] *error(abs(x2-x1))*: %.3e, x: %.3e, df=%.3e, func1(x)=%.3e, func2(x)=%.3e' % (
                count, abs(x1_final - x2_final), x1_final, func1(x1_final) - func2(x1_final), func1(x1_final),
                func2(x1_final)))
            return x1_final
        else:
            print('[%s] *error(abs(x2-x1))*: %.3e, x: %.3e, df=%.3e, func1(x)=%.3e, func2(x)=%.3e' % (
                count, abs(x1_final - x2_final), x1_final, func1(x1_final) - func2, func1(x1_final),
                func2))
            return x1_final


def ydx(x, func, i1, i2):
    """
    :param x:
    :param func:
    :param i1: 起點的指標(index)
    :param i2: 終點的指標，最大為 len(x)-1
    :return:
    """
    ans = 0
    if i2 > i1:
        for i in range(i1, i2):
            ans += (func[i + 1] + func[i]) / 2 * (x[i + 1] - x[i])
        return ans
    elif i2 < i1:
        raise BaseException("i2(%s) < i1(%s)" % (i2, i1))
    else:
        return 0  # raise BaseException("積分範圍為空集合，i2(%s) = i1(%s)" % (i2, i1))


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


def saveData(name, label, x, y):
    # [1] /Users/Ethan/PycharmProjects/Python/Practice/Cache_test.py
    # [2] /Users/Ethan/PycharmProjects/Python/Practice/
    # [3] name.csv
    # [4] /Users/Ethan/PycharmProjects/Python/Practice/data/
    # [5] /Users/Ethan/PycharmProjects/Python/Practice/data/name.csv
    location = inspect.stack()[1][1]  # [1]
    directory = os.getcwd() + '/'  # [2]
    filename = name + '.csv'  # [3]
    save_directory = directory + 'data/'  # [4]
    save_file = save_directory + filename  # [5]

    if not os.path.isfile(save_file):
        with open(save_file, 'w') as file:
            SaveData = csv.writer(file, delimiter=',')
            SaveData.writerow(label.split(','))

    with open(save_file, 'a+') as file:
        SaveData = csv.writer(file, delimiter=',')
        SaveData.writerow([x, y])


def saveArray(filename, x, y, xlabel, ylabel):
    # [1] /Users/Ethan/PycharmProjects/Python/Practice/Cache_test.py
    # [2] /Users/Ethan/PycharmProjects/Python/Practice/
    # [3] /Users/Ethan/PycharmProjects/Python/Practice/data/
    # [4] /Users/Ethan/PycharmProjects/Python/Practice/data/filename.csv
    location = inspect.stack()[1][1]  # [1]
    directory = os.getcwd() + '/'  # [2]
    save_directory = directory + 'data/'  # [3]
    save_file = save_directory + filename  # [4]

    if not os.path.isfile(save_file):
        with open(save_file, 'w') as file:
            SaveData = csv.writer(file, delimiter=',')
            SaveData.writerow([xlabel, ylabel])

    with open(save_file, 'a+') as file:
        SaveData = csv.writer(file, delimiter=',')
        for i in range(len(x)):
            SaveData.writerow([x[i], y[i]])


def save_multiple_arrays(filename, data_dict, folder=None, overwrite='No'):
    #  directory = os.getcwd() + '/'  # [2]
    #  save_directory = directory + folder + '/'  # [3]
    if folder is None:
        save_file = '/Users/Ethan/PycharmProjects/Python/Avalanche Photodiode/Paper/FigData/' + filename
    else:
        save_file = folder + filename

    if not os.path.isfile(save_file):
        flag = 1
    else:
        if overwrite == 'Yes':
            flag = 1
        else:
            ans = input('[%s] has been created. Are you sure replace it? (Y/N): ' % filename)
            if ans == 'Y':
                flag = 1
            elif ans == 'N':
                flag = 0
            else:
                raise BaseException("Wrong answer (Y/N): %s" % ans)

    if flag == 1:
        def data(index):
            l = []
            for array in data_dict.values():
                if index > len(array) - 1:
                    l.append('')
                else:
                    l.append('%s' % array[index])
            return l

        with open(save_file, 'w', newline='') as file:
            SaveData = csv.writer(file, delimiter=',')
            SaveData.writerow([key for key in data_dict.keys()])
            max_length = max([len(value) for value in data_dict.values()])
            for i in range(max_length):
                SaveData.writerow(data(i))
    data_dict.clear()
    print('[Saved] %s' % filename)


def orderfunction(temp_x, temp_y, form):
    x = temp_x.copy()
    y = temp_y.copy()
    if form == 'increasing':
        for i in range(len(x)):
            for j in range(i + 1, len(x)):
                if x[j] < x[i]:
                    temp1 = x[j]
                    temp2 = y[j]
                    x[j] = x[i]
                    y[j] = y[i]
                    x[i] = temp1
                    y[i] = temp2
    elif form == 'decreasing':
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
        raise BaseException("Wrong type : %s" % form)
    return x, y


def interpolation(xi, xf, yi, extra='Yes'):
    for i in range(1, len(xi) - 1):
        if (xi[i + 1] - xi[i]) * (xi[i] - xi[i - 1]) < 0:
            raise BaseException("xi is not a monotonic function where "
                                "xi[%s], xi[%s], xi[%s] are %s, %s, %s, respectively."
                                % (i - 1, i, i + 1, xi[i - 1], xi[i], xi[i + 1]))
        elif xi[i] == xi[i + 1]:
            raise BaseException("有兩點的數值相同：xi[%s] = xi[%s] = %s and yi[%s]=%s, yi[%s]=%s" %
                                (i, i + 1, xi[i], i, yi[i], i + 1, yi[i + 1]))
    for i in range(len(xi) - 2):
        if xi[i] > xi[i + 1]:
            print('It is an decreasing function. Changing to increasing function...')
            xi, yi = orderfunction(xi, yi, 'increasing')
    yf = np.empty([0])
    for f in range(len(xf)):
        if f != len(yf):
            raise BaseException("Error")
        if xf[f] < xi[0]:
            if extra == 'Yes':
                dx = xf[f] - xi[0]
                dx0 = xi[1] - xi[0]
                dy0 = yi[1] - yi[0]
                y_fit = yi[0] + dx * dy0 / dx0
                yf = np.append(yf, y_fit)
            else:
                yf = np.append(yf, 0)
        elif xf[f] > xi[-1]:
            if extra == 'Yes':
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
                yf = np.append(yf, 0)
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
    """
    將兩個陣列合併，重複點會只取一個
    :param v1: 陣列 1
    :param v2: 陣列 2
    :param form: 預設為 'np.array'，也可以輸入 'list'。
    :return: 新陣列
    """
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


def find(x, func, x_value, form, extrapolation='No'):
    curve = None
    for i in range(len(x) - 1):
        if x[i] < x[i + 1]:
            if curve is None:
                curve = 'increase'
            elif curve == 'decrease':
                x, func = orderfunction(x, func, 'decreasing')
                print('Not monotonic. Automatically transform into a decreasing function.')
                break
        elif x[i] > x[i + 1]:
            if curve is None:
                curve = 'decrease'
            elif curve == 'increase':
                x, func = orderfunction(x, func, 'increasing')
                print('Not monotonic. Automatically transform into an increasing function.')
                break

    for i in range(len(x) - 1):
        '''
        這段是用來外插的，因為有些數據並不能這樣做，所以這段程式碼要再找時間修補。
                elif (x_value > x[-1]) and (curve == 'increase'):
                    goal = (func[-1] - func[-2]) / (x[-1] - x[-2]) * (x_value - x[-1])
                    return goal
                '''
        if (x_value - x[0]) * (x_value - x[-1]) > 0:
            if extrapolation == 'Yes':
                if x_value < x[0]:
                    m = (func[1] - func[0]) / (x[1] - x[0])
                    return m * (x_value - x[0]) + func[0]
                else:
                    m = (func[-1] - func[-2]) / (x[-1] - x[-2])
                    return m * (x_value - x[-1]) + func[-1]
            else:
                raise BaseException("x[0]=%s, x[-1]=%s, x_value=%s" % (x[0], x[-1], x_value))
        if x[-1] == x_value:
            return func[-1]
        elif x[i] == x_value:
            return func[i]
        elif (x[i] < x_value < x[i + 1]) and (curve == 'increase'):
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
        elif (x[i + 1] < x_value < x[i]) and (curve == 'decrease'):
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
        elif i == len(x) - 2:
            print(x)
            raise BaseException("Strange, curve type is '%s', x_value=%s, x[0]=%s, x[-1]=%s" % (curve, x_value, x[0], x[-1]))


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


def parse_text(text_path, xlabel, ylabel, show='No'):
    # Open file
    fp = open(text_path, "r",encoding='utf-8', errors='ignore')
    line = fp.readline()
    X_TurnKey = 'Off'
    Y_TurnKey = 'Off'
    X = []
    Y = []
    # 用 while 逐行讀取檔案內容，直至檔案結尾
    while line:
        if X_TurnKey == 'On':
            X.append(float(line.split(" ")[1]))
        if Y_TurnKey == 'On':
            Y.append(float(line.split(" ")[2].split("\n")[0]))
        if len(line.split(" ")) >= 3 and line.split(" ")[1] == xlabel:
            X_TurnKey = 'On'
        if len(line.split(" ")) >= 3 and line.split(" ")[2].split("\n")[0] == ylabel:
            Y_TurnKey = 'On'
        line = fp.readline()

    fp.close()
    if show == 'Yes':
        print('[utils.parse_text]  Read file: %s' % text_path)
    return np.asarray(X), np.asarray(Y)


def parse_csv(csv_path, key, Cname, X, floatKey, show='Yes'):
    """
    :param csv_path: *.csv 的位址
    :param key: 1-result
                2-Curvename
                3-result[Cname][X]
                4-result, Curve_name
    :param Cname: 如果 key == 3，那麼會輸出此曲線的數據（ x 座標或 y 座標）。
    :param X: 如果 key == 3，會輸出 Cname 曲線的 X 座標。
    :param floatKey: 如果要將所有數據都轉為 float type，那麼寫 'Yes'，否則為 'No'。
    :return: 取決於 key。
    """
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
    if show == 'Yes':
        print('[utils.parse_csv]  Read file: %s' % csv_path)
    if key == 1:
        return result
    elif key == 2:
        return Curve_name
    elif key == 3:
        return result[Cname][X]
    else:
        return result, Curve_name