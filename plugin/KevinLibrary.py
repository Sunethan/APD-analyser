import numpy as np


def find_index(x, value):
    """
    並不限定 x 的單調性與重複性，因為它僅僅是找這個元素的指標。
    因此，可能有多個元素，也可能沒有元素。
    :param x: [3,5,1,2,3,5,1,20,2,38]
    :param value: 5
    :return: [1, 5]
    """
    temp = []  # 給接下來多個元素位址的暫存空間
    for i in range(len(x)):
        if x[i] == value:
            temp.append(i)
    return temp  # [] or [1,5] or [8]...


def ydx(x, y, x_start_value=None, x_end_value=None):
    if not is_monotonic(x):
        x, y = sort_function(x, y, form='increasing')
    if is_duplicate(x):
        raise BaseException("幹")

    # 從現在起，保證 x, y 為單調且不重複的數據

    # 判斷積分範圍是否合理
    if x_start_value < x[0] or x_end_value > x[-1]:
        raise BaseException("積分範圍太大。")

    if x_start_value is None:
        x_start_index = 0
    else:
        x_start_index = find_index(x, x_start_value)  # [] or [number]

    if x_end_value is None:
        length = len(x) - 1
        x_end_index = length
    else:
        x_end_index = find_index(x, x_end_value)  # [] or [number]

    integral = 0
    if len(x_start_index) == 1 and len(x_end_index) == 1:
        for i in range(x_start_index[0], x_end_index[0]):
            y_avg = (y[i] + y[i + 1]) / 2
            dx = x[i + 1] - x[i]
            dA = y_avg * dx
            integral = integral + dA
    elif len(x_start_index) == 0 and len(x_end_index) == 1:
        x_start_index_replace = None
        for i in range(len(x)):
            # 因為 x 陣列已調為對指標 i 而言的單調遞增，所以遇到第一個比起始值（start value）
            # 還要大的數字時，此數字之指標 i 即為新積分的起始指標。
            # 若指標為 -1，則表示從未儲存過起始指標。
            # 反之，若已儲存過，那麼因為指標必大於零，不可能有 -1 的指標，所以就不會再進入 if 程式中。
            if x[i] > x_start_value:
                x_start_index_replace = i
                break

        # 將涵蓋在 (start, end) 內的「有指標起終點」範圍積分
        for i in range(x_start_index_replace, x_end_index[0]):
            y_avg = (y[i] + y[i + 1]) / 2
            dx = x[i + 1] - x[i]
            dA = y_avg * dx
            integral = integral + dA

        # 計算前段面積
        dx_front = x_start_value - x[x_start_index_replace - 1]
        slope_front = (y[x_start_index_replace] - y[x_start_index_replace - 1]) / \
                      (x[x_start_index_replace] - x[x_start_index_replace - 1])
        y_start_interpolate = y[x_start_index_replace - 1] + slope_front * dx_front
        dA_front = 0.5 * (y_start_interpolate + y[x_start_index_replace]) * \
                   (x[x_start_index_replace] - x_start_value)

        return integral + dA_front
    elif len(x_start_index) == 1 and len(x_end_index) == 0:
        x_end_index_replace = None
        for i in range(len(x)):
            # 因為 x 陣列已調為對指標 i 而言的單調遞增，所以遇到第一個比起始值（start value）
            # 還要大的數字時，此數字之指標 i 即為新積分的起始指標。
            # 若指標為 -1，則表示從未儲存過起始指標。
            # 反之，若已儲存過，那麼因為指標必大於零，不可能有 -1 的指標，所以就不會再進入 if 程式中。
            if x[i] > x_end_value:
                x_end_index_replace = i - 1
                break
        # 將涵蓋在 (start, end) 內的「有指標起終點」範圍積分
        for i in range(x_start_index[0], x_end_index_replace):
            y_avg = (y[i] + y[i + 1]) / 2
            dx = x[i + 1] - x[i]
            dA = y_avg * dx
            integral = integral + dA

        dx_back = x_end_value - x[x_end_index_replace]
        slope_back = (y[x_end_index_replace + 1] - y[x_end_index_replace]) / \
                     (x[x_end_index_replace + 1] - x[x_end_index_replace])
        y_end_interpolate = y[x_end_index_replace] + slope_back * dx_back
        dA_end = 0.5 * (y_end_interpolate + y[x_end_index_replace + 1]) * \
                 (x_end_value - x[x_end_index_replace])
        integral = integral + dA_end
        return integral
    else:
        x_start_index_replace = -1  # 用以防止重複儲存起始指標的標記（marker）
        x_end_index_replace = None
        for i in range(len(x)):
            # 因為 x 陣列已調為對指標 i 而言的單調遞增，所以遇到第一個比起始值（start value）
            # 還要大的數字時，此數字之指標 i 即為新積分的起始指標。
            # 若指標為 -1，則表示從未儲存過起始指標。
            # 反之，若已儲存過，那麼因為指標必大於零，不可能有 -1 的指標，所以就不會再進入 if 程式中。
            if x[i] > x_start_value and x_start_index_replace == -1:
                x_start_index_replace = i
            if x[i] > x_end_value:
                x_end_index_replace = i - 1
                break
        # 將涵蓋在 (start, end) 內的「有指標起終點」範圍積分
        for i in range(x_start_index_replace, x_end_index_replace):
            y_avg = (y[i] + y[i + 1]) / 2
            dx = x[i + 1] - x[i]
            dA = y_avg * dx
            integral = integral + dA

        dx_front = x_start_value - x[x_start_index_replace - 1]
        slope_front = (y[x_start_index_replace] - y[x_start_index_replace - 1]) / \
                      (x[x_start_index_replace] - x[x_start_index_replace - 1])
        y_start_interpolate = y[x_start_index_replace - 1] + slope_front * dx_front
        dA_front = 0.5 * (y_start_interpolate + y[x_start_index_replace]) * \
                         (x[x_start_index_replace] - x_start_value)

        dx_back = x_end_value - x[x_end_index_replace]
        slope_back = (y[x_end_index_replace + 1] - y[x_end_index_replace]) / \
                     (x[x_end_index_replace + 1] - x[x_end_index_replace])
        y_end_interpolate = y[x_end_index_replace] + slope_back * dx_back
        dA_end = 0.5 * (y_end_interpolate + y[x_end_index_replace + 1]) * \
                       (x_end_value - x[x_end_index_replace])
        integral = integral + dA_front + dA_end
    return integral


def dydx(x, y):
    if not is_monotonic(x):
        x, y = sort_function(x, y, form='increasing')
    if is_duplicate(x):
        raise BaseException("幹")

    # 從現在起，保證 x, y 為單調且不重複的數據

    temp = [(y[1] - y[0]) / (x[1] - x[0])]
    if len(y) >= 3:
        for i in range(1, len(y) - 1):
            dydx_front = (y[i] - y[i - 1]) / (x[i] - x[i - 1])
            dydx_back = (y[i + 1] - y[i]) / (x[i + 1] - x[i])
            temp.append((dydx_front + dydx_back) / 2)
    temp.append((y[-1] - y[-2]) / (x[-1] - x[-2]))
    return np.array(temp)


def find_y_from_x(x, y, x_given):
    if not is_monotonic(x):
        x, y = sort_function(x, y, form='increasing')
    if is_duplicate(x):
        raise BaseException("幹")

    # 從現在起，保證 x, y 為單調且不重複的數據

    for index, value in enumerate(x):
        if value == x_given:
            return y[index]
        elif index < len(x) - 1:
            if (x[index] - x_given) * (x[index + 1] - x_given) < 0:
                dx = x[index + 1] - x[index]
                dy = y[index + 1] - y[index]
                slope = dy / dx
                y_goal = y[index] + slope * (x_given - x[index])
                return y_goal
        else:
            if index == len(x) - 1:
                if x_given > x[-1]:
                    raise BaseException("你給的x值(%s)大於最大值(%s)" % (x_given, x[-1]))
                elif x_given < x[0]:
                    raise BaseException("你給的x值(%s)小於最小值(%s)" % (x_given, x[0]))
                else:
                    raise BaseException("程式碼的邏輯有問題！")


def find_x_from_y(x, y, y_given):
    if not is_monotonic(x):
        x, y = sort_function(x, y, form='increasing')
    if is_duplicate(x):
        raise BaseException("幹")

    # 從現在起，保證 x, y 為單調且不重複的數據

    x_goal = []
    final_index = len(y) - 1
    for index, value in enumerate(y):
        if value == y_given:
            x_goal.append(x[index])
        elif index < final_index:
            if (y[index] - y_given) * (y[index + 1] - y_given) < 0:
                dx = x[index + 1] - x[index]
                dy = y[index + 1] - y[index]
                slope = dy / dx
                x_find = x[index] + (y_given - y[index]) / slope
                x_goal.append(x_find)
        else:
            if index == final_index and len(x_goal) == 0:
                raise BaseException("沒找到%s哦！" % y_given)
    return x_goal


def is_monotonic(x):
    """
    邏輯：檢查 (x[i+1] - x[i]) * (x[i+2] - x[i+1]) 是否恆正？若為恆正，則具有單調性，否則需重新排序。
    i 是什麼？i 是從 0 開始， 但對任何 i 而言，程式也會讀取 x[i+2]，所以需特別留意 for i in range(...)
    的 i 列舉範圍。

    一般而言，倘若允許 i 遍歷陣列 x 中的所有元素，那麼 range 內應輸入 len(x)。
    因此，針對會同時讀取 x[i] 與 x[i+2] 的程式碼而言，我們應該於 range 內填入 len(x) - 2。
    :param x:
    :return:
    """
    if isinstance(x, (list, np.ndarray)):
        total_length_minus_two = len(x) - 2
        for i in range(0, total_length_minus_two):
            key = (x[i + 1] - x[i]) * (x[i + 2] - x[i + 1])
            if key < 0:
                return False
        return True
    else:
        raise BaseException("這不是 list or numpy.ndarray，而是 %s！" % type(x))


def is_duplicate(x):
    """
    思路：
    第 1 個元素(i=0)，有沒有跟剩餘的 N-1 個元素相同？ (N-1)：(0,1) & (0,2) & (0,3) & ....
    第 2 個元素(i=1)，有沒有跟剩餘的 N-2 個元素相同？ (N-2)
    ....
    第 N-1 個元素(i=N-2)，有沒有跟剩餘的 1 個元素相同？ (1)
    總共要 (N-1) + (N-2) + ... 1 = (N-1) * N/2 ~ N^2

    假設：10 個元素要 10 秒，那 100 個元素 10 * 100 = 1000 秒。
    :param x:
    :return:
    """
    for i in range(len(x) - 1):
        for j in range(i + 1, len(x)):
            if x[i] == x[j]:
                return True
    return False


def sort_function(x, y, form='increasing'):
    """
    考慮：
    x = [3,2,1,5]
    y = [5,1,6,8]
    若想將其改為 x 遞增，則為：
    x = [1,2,3,5]
    y = [6,1,5,8]
    若想將其改為 x 遞減：則為：
    x = [5,3,2,1]
    y = [8,5,1,6]

    作法，以遞增為例：
    x = [3,2,1,5] >> [1,2,3,5]
    比較第 1 個元素(A)與剩餘 N-1 元素(B)的大小。若 B<A，則能把最小元素放第 1 個。
    比較第 2 個元素(A)與剩餘 N-2 元素(B)的大小。若 B<A，則能把第二小元素放第 2 個。
    ...
    不過，現在除了要更改 x 以外，也要更改 y。因此，交換 x 時，也需要交換 y。
    :param x:
    :param y:
    :param form:
    :return:
    """
    x = x.copy()
    y = y.copy()
    if form == 'increasing':
        for i in range(0, len(x) - 1):
            for j in range(i + 1, len(x)):
                if x[j] < x[i]:
                    x_temp, y_temp = x[i], y[i]
                    x[i], y[i] = x[j], y[j]
                    x[j], y[j] = x_temp, y_temp
        return x, y
    elif form == 'decreasing':
        for i in range(0, len(x) - 1):
            for j in range(i + 1, len(x)):
                if x[j] > x[i]:
                    x_temp, y_temp = x[i], y[i]
                    x[i], y[i] = x[j], y[j]
                    x[j], y[j] = x_temp, y_temp
        return x, y
    else:
        raise BaseException("打錯了！(increasing/decreasing): %s" % form)
