import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import utils
import os
import time
import csv
import glob
from DataAnalysis import Read
from tempfile import NamedTemporaryFile
import shutil
import plotly.graph_objects as go
start_time = time.time()
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (7.5, 7),
          'axes.labelsize': 'x-large',
          'axes.titlesize': 'x-large',
          'xtick.labelsize': 'x-large',
          'ytick.labelsize': 'x-large'}
pylab.rcParams.update(params)
plt.rcParams.update({'font.size': 13})
import warnings
warnings.simplefilter("error")

ColorSet = ['dodgerblue', 'g', 'goldenrod', 'darkviolet', 'b', 'r', 'brown', 'darkorange',
            'hotpink', 'fuchsia', 'yellowgreen', 'royalblue', 'tomato', 'purple', 'olive', 'darkgreen']
LineStyleSet = ['-', '-.', '--', ':', (0, (3, 1, 1, 1)), (0, (5, 5)), (0, (3, 1, 1, 1, 1, 1)),
                (0, (3, 5, 1, 5)), (0, (5, 10))]


class IVData(object):
    def __init__(self, database, path):
        def set_name(m, d, filename):
            if 'D' in filename and 'F' not in filename:
                Code = int(filename.split("D")[0])
            elif 'L' in filename and 'F' not in filename:
                Code = int(filename.split("L")[0])
            elif 'F' in filename:
                Code = int(filename.split("F")[0])
            else:
                raise BaseException('"D" and "L" are not in the file name (%s). Please check the filename list.' % filename)
            Info = self.code_info(Code, show='No')
            return str(m) + '/' + str(d) + '-' + filename.split(".")[0] + '-%s' % [value for key, value in Info.items()]
        self._database = database
        self._path = path
        self._filename = set_name(path.split("/")[-3][:2], path.split("/")[-3][2:], path.split("/")[-1])
        info_name_list = self.key(1) + self.key(2)  # 這必須與 reverse_filename 對應
        self._info = {info_name_list[index]: element for index, element in enumerate(self.reverse_filename(self._filename))}

    @staticmethod
    def key(number):
        key1 = ['date', 'code', 'light', 'number']
        key2 = ['active', 'FGR', 'AGR', 'order', 'depth', 'Y']
        if number == 1:
            return key1
        elif number == 2:
            return key2
        elif number == 'all':
            return key1 + key2
        else:
            raise BaseException("Wrong key(1/2): %s" % number)

    @staticmethod
    def code_info(code, show='Yes'):
        """
        讀取編碼意義，例如 2613 為主動區直徑 400um，FGR spacing 4um、AGR直徑 450um、擴散順後先大後小、擴散深度(depth)不明。
        :param code: 元件代碼，如 2613。
        :param show: 印出結果，預設值為 'Yes'。
        :return: {'active': 400, 'FGR': 4, 'AGR': 450, 'order': 2, 'depth': 'null', 'Y': 3}

        >> from DataAnalysis import IVDataBase as IVDB
        >> DB = IVDB(data_path)：讀取所有數據
        >> DB.code_info(2613)
        >> Code [2613] is: {'active': 400, 'FGR': 4, 'AGR': 450, 'order': 2, 'depth': 'null', 'Y': 6}
        """
        if isinstance(code, int) and len(str(code)) == 4:
            X, Y, x, y = str(code)
            info = None
            # Active: 主動區直徑
            # FGR: 懸護環與側護環的距離
            # AGR: 側護環直徑
            # order: 1 為先小再大（主動區->淺護環），2 為先大再小（淺護環->主動區）
            NoneValue = 'null'
            if Y == '6' and y == '1':
                info = {'active': 1600, 'FGR': 4, 'AGR': 1650, 'order': 2, 'depth': NoneValue, 'Y': int(Y)}
            elif Y == '6' and y == '2':
                info = {'active': 800, 'FGR': 4, 'AGR': 850, 'order': 2, 'depth': NoneValue, 'Y': int(Y)}
            elif Y == '6' and y == '3':
                info = {'active': 400, 'FGR': 4, 'AGR': 450, 'order': 2, 'depth': NoneValue, 'Y': int(Y)}
            elif Y == '7' and y == '1':
                info = {'active': 200, 'FGR': 4, 'AGR': 250, 'order': 2, 'depth': NoneValue, 'Y': int(Y)}
            elif Y == '7' and y == '2':
                info = {'active': 100, 'FGR': 4, 'AGR': 150, 'order': 2, 'depth': NoneValue, 'Y': int(Y)}
            elif Y == '7' and y == '3':
                info = {'active': 50, 'FGR': 4, 'AGR': 100, 'order': 2, 'depth': NoneValue, 'Y': int(Y)}
            elif Y == '5' and y == '1':
                info = {'active': 200, 'FGR': NoneValue, 'AGR': 250, 'order': 2, 'depth': NoneValue, 'Y': int(Y)}
            elif Y == '5' and y == '2':
                info = {'active': 200, 'FGR': 5, 'AGR': 250, 'order': 2, 'depth': NoneValue, 'Y': int(Y)}
            elif Y == '5' and y == '3':
                info = {'active': 200, 'FGR': 4, 'AGR': 250, 'order': 2, 'depth': NoneValue, 'Y': int(Y)}
            elif Y == '5' and y == '4':
                info = {'active': 200, 'FGR': 3, 'AGR': 250, 'order': 2, 'depth': NoneValue, 'Y': int(Y)}
            elif Y == '4' and y == '1':
                info = {'active': 200, 'FGR': 4, 'AGR': 250, 'order': 1, 'depth': NoneValue, 'Y': int(Y)}
            elif Y == '4' and y == '2':
                info = {'active': 100, 'FGR': 4, 'AGR': 150, 'order': 1, 'depth': NoneValue, 'Y': int(Y)}
            elif Y == '4' and y == '3':
                info = {'active': 50, 'FGR': 4, 'AGR': 100, 'order': 1, 'depth': NoneValue, 'Y': int(Y)}
            elif Y == '3' and y == '3':
                info = {'active': 200, 'FGR': 3, 'AR': 250, 'order': 1, 'depth': NoneValue, 'Y': int(Y)}
            elif Y == '3' and y == '2':
                info = {'active': 200, 'FGR': 5, 'AGR': 250, 'order': 1, 'depth': NoneValue, 'Y': int(Y)}
            elif Y == '3' and y == '1':
                info = {'active': 200, 'FGR': NoneValue, 'AGR': 250, 'order': 1, 'depth': NoneValue, 'Y': int(Y)}
            elif Y == '2' and x == '1' and y == '3':
                info = {'active': 50, 'FGR': 4, 'AGR': 110, 'order': 2, 'depth': NoneValue, 'Y': int(Y)}
            elif Y == '2' and x == '2' and y == '3':
                info = {'active': 50, 'FGR': 4, 'AGR': 100, 'order': 2, 'depth': NoneValue, 'Y': int(Y)}
            elif Y == '2' and x == '3' and y == '3':
                info = {'active': 50, 'FGR': 4, 'AGR': 90, 'order': 2, 'depth': NoneValue, 'Y': int(Y)}
            elif Y == '2' and x == '1' and y == '2':
                info = {'active': 100, 'FGR': 4, 'AGR': 160, 'order': 2, 'depth': NoneValue, 'Y': int(Y)}
            elif Y == '2' and x == '2' and y == '2':
                info = {'active': 100, 'FGR': 4, 'AGR': 150, 'order': 2, 'depth': NoneValue, 'Y': int(Y)}
            elif Y == '2' and x == '3' and y == '2':
                info = {'active': 100, 'FGR': 4, 'AGR': 140, 'order': 2, 'depth': NoneValue, 'Y': int(Y)}
            elif Y == '2' and x == '1' and y == '1':
                info = {'active': 200, 'FGR': 4, 'AGR': 260, 'order': 2, 'depth': NoneValue, 'Y': int(Y)}
            elif Y == '2' and x == '2' and y == '1':
                info = {'active': 200, 'FGR': 4, 'AGR': 250, 'order': 2, 'depth': NoneValue, 'Y': int(Y)}
            elif Y == '2' and x == '3' and y == '1':
                info = {'active': 200, 'FGR': 4, 'AGR': 240, 'order': 2, 'depth': NoneValue, 'Y': int(Y)}
            elif Y == '8' or Y == '9':
                info = {'active': 240, 'FGR': NoneValue, 'AGR': NoneValue, 'order': NoneValue, 'depth': int(x),
                        'Y': int(Y)}
            elif Y == '1':
                if y == '1':
                    info = {'active': 200, 'FGR': 5, 'AGR': 250, 'order': 2, 'depth': NoneValue, 'Y': int(Y)}
                elif y == '2':
                    info = {'active': 200, 'FGR': 5, 'AGR': 250, 'order': 2, 'depth': NoneValue, 'Y': int(Y)}
                elif y == '3':
                    info = {'active': 200, 'FGR': 5, 'AGR': 250, 'order': 2, 'depth': NoneValue, 'Y': int(Y)}
                elif y == '4':
                    info = {'active': 200, 'FGR': 5, 'AGR': 250, 'order': 2, 'depth': NoneValue, 'Y': int(Y)}
                else:
                    raise BaseException("Wrong code: %s" % code)
            else:
                raise BaseException("%s isn't prepared" % code)
            if show == 'Yes':
                print('Code [%s] is: %s' % (code, info))
            return info
        else:
            raise BaseException("Wrong code format: %s" % code)

    @staticmethod
    def reverse_filename(filename):
        date = filename.split("-")[0].split("/")[0] + filename.split("-")[0].split("/")[1]
        code = filename.split("-")[1][:-1]
        light = filename.split("-")[1][-1]
        number = filename.split("-")[2]
        info = [filename.split("[")[1].split("]")[0].split(",")[0]] + \
               [value.split(" ")[1] for i, value in enumerate(filename.split("[")[1].split("]")[0].split(",")) if i > 0]

        return [date, code, light, number] + info  # 滿足 key1 + key2 的格式

    def name(self):
        return self._filename

    def path(self, relative='off'):
        if relative == 'off':
            return self._path
        elif relative == 'on':
            return self._path[len(self._database) + 1:]

    def info(self):
        self._info['X'] = self._info['code'][0]
        self._info['x'] = self._info['code'][2]
        self._info['y'] = self._info['code'][3]

        Xmin = 1
        Ymin = 1
        Xmax = 13
        Ymax = 9
        X_List = np.arange(Xmin, Xmax + 1, 1)
        xmax = {posX: {1: 2, 2: 3, 3: 2, 4: 3, 5: 2, 6: 1, 7: 3, 8: 3, 9: 3} for posX in X_List}
        ymax = {posX: {1: 4, 2: 3, 3: 3, 4: 3, 5: 4, 6: 3, 7: 3, 8: 3, 9: 3} for posX in X_List}
        dx = 1 / (xmax[1][int(self._info['Y'])] + 1)
        dy = 1 / (ymax[1][int(self._info['Y'])] + 1)
        d = {'X': dx, 'Y': dy}
        self._info['Xpos'] = "{0:.2f}".format(float(self._info['X']) - 1 + d['X'] * float(self._info['x']))
        self._info['Ypos'] = "{0:.2f}".format(float(self._info['Y']) - 1 + d['Y'] * float(self._info['y']))
        return self._info


class IVDataBase(object):
    def __init__(self, location, date='*', code='*', light='*', number='*', **kwargs):
        """
        只要資料庫符合 date/code/codeD-1.TXT 格式，就能夠藉已知參數快速讀取彙整資料。
        :param location:資料庫目錄。例如 /user/..../IV。而 IV 資料夾中需有實驗日期目錄，如 /IV/1120 為11月20日的量測結果。
        :param kwargs:可輸入條件為 date, code, light, number, Y, active, FGR, AGR, order, depth
            date: 日期，如1120
            code: 元件座標代碼，如2613
            light: 照光為'L'，沒照光則為'D'
            number: 實驗編號，例如 2613D-3.TXT 的編號為 3
            Y: 元件Y座標，用以指定特定比較組。例如Y=3表示不同FGR間隔的組別，Y=6為不同元件大小的組別。
            active: 主動區直徑（不包含與主動區接觸的側護環），單位為 um
            FGR: 浮護環與側護環的間隔，單位為 um
            AGR: 與主動區相接觸的側護環直徑，單位為 um
            order: 擴散次序。先擴散主動區再擴散護環為 1，先擴散護環再擴散主動區為 2
            depth: 擴散深度，僅限於 test key 使用。1：淺，2：中，3：深

        >> from DataAnalysis import IVDataBase as IVDB
        >> DB = IVDB(data_path)：讀取所有數據
        >> DB = IVDB(data_path, date=1120, active=800, light='D')：選擇主動區直徑為 800um 的11月20號暗電流量側數據
        >> DB = IVDB(data_path, code=3612, light='L')：選擇編號為 3612 的所有日期光電流量側數據
        >> DB = IVDB(data_path, Y=6, active=400）：在不同元件尺寸組別(Y=6)中，選擇主動區直徑為400um的數據。
        """

        def data_path():
            return location + '/' + str(date) + '/' + str(code) + '/' + \
                   str(code) + str(light) + '-' + str(number) + '.TXT'

        def coordinate(data):
            flag = 0
            for key in kwargs:
                if key in ['X', 'Y', 'x', 'y']:
                    if not isinstance(kwargs[key], list):
                        kwargs[key] = [kwargs[key]]
                    if int(data.info()[key]) not in kwargs[key]:
                        flag = 1
            if flag == 0:
                return True
            else:
                return False

        def path_check(path):
            filename = path.split("/")[-1]
            code = path.split("/")[-2]
            if code in filename:
                if '-' in filename and '.TXT' in filename and ('L' in filename or 'D' in filename):
                    if filename.split("-")[0][:-1].isdigit() and (filename.split("-")[0][-1] == 'L'
                                                                  or filename.split("-")[0][-1] == 'D'):
                        if filename.split("-")[1].split(".")[0].isdigit():
                            return True
                        else:
                            return False
                    else:
                        return False
                else:
                    return False
            else:
                print('Check %s folder: %s' % (code, path))
                return False

        # 輸入資料
        self._date = date
        self._code = code
        self._light = light
        self._number =number

        # 預設電性
        self._EC = ['Vpt', 'Vb', 'Vpt2', 'Vr']

        # 檢查 key 是否衝突
        self._location = location
        self._key_contradiction(2, kwargs)

        # 開始讀檔
        all_data = [IVData(location, path) for path in glob.glob(data_path()) if path_check(path)]
        self._data = []

        for data in all_data:
            if coordinate(data):
                flag = 0
                kwargs_left = list(set(kwargs.keys()) - {'X', 'Y', 'x', 'y'})
                for key in kwargs_left:
                    if not isinstance(kwargs[key], list):
                        kwargs[key] = [kwargs[key]]
                    for element in kwargs[key]:
                        if data.info()[key] == str(element):
                            flag = 1
                            break
                if flag == 1 or (flag == 0 and len(kwargs_left) == 0):
                    self._data.append(data)

        if len(self._data) == 0:
            raise BaseException("No such curve!")

    def _key_contradiction(self, number, kwargs):
        key = IVData.key(number) + ['X', 'x', 'y']
        for set_key in kwargs:
            if set_key not in key:
                raise BaseException('"%s" is not the default key %s' % (set_key, key))
            if number == 2:
                if set_key in key and (self._code != '*' or self._number != '*'):
                    raise BaseException('(code, number) and (X, Y, x, y, active, FGR, AGR, order) '
                                        'cannot be simultaneously indicated.')
        if self._date == '*' and self._number != '*':
            raise BaseException("Indication of number while date is arbitrary is meaningless.")
        if self._number != '*' and self._date == '*' and self._code == '*' and self._light == '*':
            raise BaseException('Only number information is not enough.')

    def print_key(self):
        """
        印出允許的 IVDB key 值。
        :return: ['date', 'code', 'light', 'number', 'active', 'FGR', 'AGR', 'order', 'depth', 'Y']

        >> IVDB.print_key()
        >> Acceptable keys: ['date', 'code', 'light', 'number', 'active', 'FGR', 'AGR', 'order', 'depth', 'Y']
        """
        key_list = IVData.key(1) + IVData.key(2)
        print('Acceptable keys: %s' % key_list)
        return key_list

    def read(self):
        """
        讀取IV數據，其資料結構為 [IV1, IV2, ...]，例如：

        >> IV = IVDB.read()  # 其資料結構為列表（list），所以是由 IV[0]、IV[1] 讀取數據。
        >> Voltage = IV[0].X  # 第一筆資料的電壓 numpy 陣列
        >> print(Voltage)
        >> [0, -0.1, -0.2, -0.3, ......]
        >> Current = IV[0].Y  # 第一筆資料的電流 numpy 陣列
        >> print(Current)
        >> [1e-12, 1.5e-12, 1.7e-12, ....]
        >> plt.plot(IV[0].X, IV[0].Y)  # 繪製第一條電性數據

        其中，numpy 擴充程式庫的陣列(Array) 與 python 內建的列表（list）不同，numpy 陣列允許數學運算。

        >> Area = np.pi * (5e-4) ** 2  # [cm2] 假設元件半徑為 5um
        >> CurrentDensity = IV[0].Y / Area  # [A/cm2] 計算電流密度

        :return: [IV1, IV2, IV3, ....]
        """

        return [Read(data.path(), 'HP4156', show='No') for data in self._data]

    def filename(self, show='Yes'):
        """
        回傳讀取到的數據名稱清單，例如：

        >> Filename = IVDB.filename()
        >> Filename list: ["11/20-1712D-4-[100, 4, 150, 2, 'null', 7]", "11/20-1712D-1-[100, 4, 150, 2, 'null', 7]"]

        其中，
        第一組數據名稱為 11/20-1712D-4-[100, 4, 150, 2, 'null', 7]
        第二組數據名稱為 11/20-1712D-1-[100, 4, 150, 2, 'null', 7]

        :return: [IV1_name, IV2_name, ....]
        """
        filename_list = []
        if show == 'Yes':
            print('[index]     mm/dd-code-number-[active, FGR, AGR, order, depth, Y]')
            for i, data in enumerate(self._data):
                filename_list.append(data.name())
                number = '[' + str(i) + ']'
                print('{:>5} {:>55}'.format(number, data.name()))
        return filename_list

    def find_index(self, name=None, **kwargs):
        """
        尋找符合設定條件的元件在 IV.read() 輸出的列表索引值。
        倘若有個名為 "11/20-1712D-4-[100, 4, 150, 2, 'null', 7]" 的數據，那麼可以藉此找到其索引。

        >> from DataAnalysis import IVDataBase as IVDB
        >> DB = IVDB(data_path)：讀取所有數據
        >> IV = DB.read()
        >> j = DB.find_index(name="11/20-1712D-4-[100, 4, 150, 2, 'null', 7]")

        倘若想在已挑選出來的數據中，再另行挑選符合某些條件的元件，那麼也可藉此挑選。

        >> j = DB.find_index(active=200)  # 挑選出主動區直徑200um的元件數據
        >> print(j)
        >> [8, 9, 10, 11]

        :param name: 預設值為 None，即可接受空白。
        :param kwargs: [date, code, light, number, Y, active, FGR, AGR, order, depth]
        :return: [index1, index2, index3, ...]
        """
        if name is not None:
            if len(kwargs) > 0:
                raise BaseException('Cannot search indicated specifications %s while using the specific file name "%s"'
                                    % (kwargs, name))
            for index, data in enumerate(self._data):
                if name == data.name():
                    return [index]
            raise BaseException('"%s" not found!' % name)
        else:
            self._key_contradiction(2, kwargs)
            index_list = []
            for i, data in enumerate(self._data):

                flag = 1
                for key, value in kwargs.items():
                    if str(value) != data.info()[key]:
                        flag = flag * 0

                if flag == 1:
                    index_list.append(i)
            if len(index_list) == 0:
                raise BaseException('Index not found: %s' % kwargs)
            else:
                return index_list

    def path(self, relative='on', show='on'):
        """
        輸出檔案路徑清單
        :param relative: 有 'on' 與 'off' 兩種選項，預設值為 'on'，無填寫則使用預設值。
                         開啟則使用相對資料庫之檔案路徑，否則為絕對路徑。
        :param show: 有 'on' 與 'off' 兩種選項，預設值為 'on'，無填寫則使用預設值。 開啟則印出各檔案之路徑。
        :return: [path1, path2, path3, ....]

        >> from DataAnalysis import IVDataBase as IVDB
        >> # 取出 11/18 測量的 1713D-1.TXT 絕對路徑。由於此為單一檔案，所以其路徑會存在列表（list）的第一個元素，最後用[0]取出。
        >> location = IVDB(data_path, date=1118, code=1713, light='D', number=1).path(relative='off', show='on')[0]
        >> [index]
        >>   [0]     /Users/Ethan/PycharmProjects/Python/Avalanche Photodiode/Experiment/IV/1118/1713/1713D-1.TXT
        """
        location = []
        if show == 'on':
            print('[index]')
        for i, data in enumerate(self._data):
            if relative == 'on' or relative == 'off':
                location.append(data.path(relative=relative))
            else:
                raise BaseException("Data-path relative parameter: on/off only.")
            if show == 'on':
                str_index = '[' + str(i) + ']'
                print('{:>5}     {:>40}'.format(str_index, location[-1]))
        return location

    def table(self, datalimit='off', Xmin=1, Xmax=6, Ymin=1, Ymax=9, **kwargs):
        """
        元件測量記錄簿，自動開啟網頁呈現「實驗目的」、「編號」、「量測結果」、「主動區直徑」等資訊。
        :param datalimit: 是否僅呈現已量測之元件資訊列表。預設值為 'off'，同時顯示已量測與未量測之元件資訊列表。
        :param Xmin: 元件X座標之最小值，預設值為 1。
        :param Xmax: 元件X座標之最大值，預設值為 6。
        :param Ymin: 元件Y座標之最小值，預設值為 1。
        :param Ymax: 元件Y座標之最大值，預設值為 9。
        :param kwargs: 可設定欲顯示之元件條件，例如主動區直徑、FGR spacing等，例如 active=200 為僅顯示直徑200um之元件數據。
                       可設定之條件：Y, active、FGR、AGR、order
        :return: 開啟記錄簿

        >> IVDB = Data.IVDataBase(data_path)
        >> IVDB.table()  # 顯示所有元件資料（已測量與未測量），且元件座標落於 1<X<6，1<Y<9。
        >> IVDB.table(Y=3)  # 只顯示 Y=3 的資料，其對應至觀察 FGR spacing 效應的測驗組（擴散順序為先小再大）
        >> IVDB.table(datalimit='on', Ymin=3, Ymax=7, active=100)  # 只顯示已測量之元件資料，且限定3<Y<7、主動區直徑為100um。

        表格欄位說明：
        Y：元件Y座標，用以快速對照位置與其設計目的
        Test：設計目的
        Code：元件座標代碼
        Measurement：量測日期與光暗電流註記
        Vpt：擊穿電壓
        Vb：崩潰電壓
        Active：主動區直徑，單位為 um
        FGR：Floating guard ring 與主動區側護環的間隔，單位為 um
        AGR：Attached guard ring，與主動區相接觸的淺護環直徑，單位為 um
        Order：主動區擴散（小實心圓）與護環擴散與（大空心環）之擴散順序
        """
        # 元件座標定義

        def day_DL(measurement, Code):
            if Code in measurement:
                return measurement[Code]
            else:
                return ''

        def code_generator(location, Xmin, Xmax, Ymin, Ymax, kwargs):
            def range_list(value):
                return range(1, value + 1)
            X_List = np.arange(Xmin, Xmax + 1, 1)
            Y_List = np.arange(Ymin, Ymax + 1, 1)
            xmax = {posX: {1: 2, 2: 3, 3: 2, 4: 3, 5: 2, 6: 1, 7: 3, 8: 3, 9: 3} for posX in X_List}
            ymax = {posX: {1: 4, 2: 3, 3: 3, 4: 3, 5: 4, 6: 3, 7: 3, 8: 3, 9: 3} for posX in X_List}
            #ymax[13] = {1: 4, 2: 3, 3: 3, 4: 3, 5: 4, 6: 2, 7: 3, 8: 3, 9: 3}
            TestName = {1: 'Pad pattern', 2: 'GR size', 3: 'FGR spacing (1)', 4: 'Active size (1)',
                        5: 'FGR spacing (2)', 6: 'Active size(2)', 7: 'Active size (2)',
                        8: 'Depth (大Pad)', 9: 'Depth (小Pad)'}
            ActiveSize = {1: {x: {y: 200 for y in range_list(ymax[Xmin][1])} for x in range_list(xmax[Xmin][1])},
                          2: {x: {1: 200, 2: 100, 3: 50} for x in range_list(xmax[Xmin][2])},
                          3: {x: {y: 200 for y in range_list(ymax[Xmin][3])} for x in range_list(xmax[Xmin][3])},
                          4: {x: {1: 200, 2: 100, 3: 50} for x in range_list(xmax[Xmin][4])},
                          5: {x: {y: 200 for y in range_list(ymax[Xmin][5])} for x in range_list(xmax[Xmin][5])},
                          6: {x: {1: 1600, 2: 800, 3: 400} for x in range_list(xmax[Xmin][6])},
                          7: {x: {1: 200, 2: 100, 3: 50} for x in range_list(xmax[Xmin][7])},
                          8: {x: {y: 240 for y in range_list(ymax[Xmin][8])} for x in range_list(xmax[Xmin][8])},
                          9: {x: {y: 240 for y in range_list(ymax[Xmin][9])} for x in range_list(xmax[Xmin][9])}}
            AGRSize = {1: {x: {y: 250 for y in range_list(ymax[Xmin][1])} for x in range_list(xmax[Xmin][1])},
                       2: {1: {1: 260, 2: 160, 3: 110},
                                   2: {1: 250, 2: 150, 3: 100},
                                   3: {1: 240, 2: 140, 3: 90}},
                       3: {x: {y: 250 for y in range_list(ymax[Xmin][3])} for x in range_list(xmax[Xmin][3])},
                       4: {x: {1: 250, 2: 150, 3: 100} for x in range_list(xmax[Xmin][4])},
                       5: {x: {y: 250 for y in range_list(ymax[Xmin][5])} for x in range_list(xmax[Xmin][5])},
                       6: {x: {1: 1650, 2: 850, 3: 450} for x in range_list(xmax[Xmin][6])},
                       7: {x: {1: 250, 2: 150, 3: 100} for x in range_list(xmax[Xmin][7])},
                       8: {x: {y: 0 for y in range_list(ymax[Xmin][8])} for x in range_list(xmax[Xmin][8])},
                       9: {x: {y: 0 for y in range_list(ymax[Xmin][9])} for x in range_list(xmax[Xmin][9])}}
            FGRSpacing = {1: {x: {y: 5 for y in range_list(ymax[Xmin][1])} for x in range_list(xmax[Xmin][1])},
                          2: {x: {y: 4 for y in range_list(ymax[Xmin][2])} for x in range_list(xmax[Xmin][2])},
                          3: {x: {1: np.inf, 2: 5, 3: 3} for x in range_list(xmax[Xmin][3])},
                          4: {x: {y: 4 for y in range_list(ymax[Xmin][4])} for x in range_list(xmax[Xmin][4])},
                          5: {x: {1: np.inf, 2: 5, 3: 4, 4: 3} for x in range_list(xmax[Xmin][5])},
                          6: {x: {y: 4 for y in range_list(ymax[Xmin][6])} for x in range_list(xmax[Xmin][6])},
                          7: {x: {y: 4 for y in range_list(ymax[Xmin][7])} for x in range_list(xmax[Xmin][7])},
                          8: {x: {y: np.inf for y in range_list(ymax[Xmin][8])} for x in range_list(xmax[Xmin][8])},
                          9: {x: {y: np.inf for y in range_list(ymax[Xmin][9])} for x in range_list(xmax[Xmin][9])}}
            Order = {1: {x: {y: '先大再小(1)' for y in range_list(ymax[Xmin][1])} for x in range_list(xmax[Xmin][1])},
                     2: {x: {y: '先大再小(2)' for y in range_list(ymax[Xmin][2])} for x in range_list(xmax[Xmin][2])},
                     3: {x: {y: '先小再大(1)' for y in range_list(ymax[Xmin][3])} for x in range_list(xmax[Xmin][3])},
                     4: {x: {y: '先小再大(1)' for y in range_list(ymax[Xmin][4])} for x in range_list(xmax[Xmin][4])},
                     5: {x: {y: '先大再小(2)' for y in range_list(ymax[Xmin][5])} for x in range_list(xmax[Xmin][5])},
                     6: {x: {y: '先大再小(2)' for y in range_list(ymax[Xmin][6])} for x in range_list(xmax[Xmin][6])},
                     7: {x: {y: '先大再小(2)' for y in range_list(ymax[Xmin][7])} for x in range_list(xmax[Xmin][7])},
                     8: {x: {y: '無護環' for y in range_list(ymax[Xmin][8])} for x in range_list(xmax[Xmin][8])},
                     9: {x: {y: '無護環' for y in range_list(ymax[Xmin][9])} for x in range_list(xmax[Xmin][9])}}

            Day_DL = {}
            for data in self._data:
                filename = data.name()
                code_tmp = int(filename.split("-")[1][:-1])
                day_tmp = filename.split("-")[0]
                light_tmp = filename.split("-")[1][-1]
                if code_tmp not in Day_DL:
                    Day_DL[code_tmp] = [day_tmp + light_tmp]
                else:
                    if (day_tmp + light_tmp) not in Day_DL[code_tmp]:
                        Day_DL[code_tmp].append(day_tmp + light_tmp)

            CodeList = dict()
            for X in X_List:
                for Y in Y_List:
                    for x in range(1, xmax[X][Y] + 1):
                        for y in range(1, ymax[X][Y] + 1):
                            if Y < 10:
                                code = X * 1000 + Y * 100 + x * 10 + y
                            else:
                                code = X * 10000 + Y * 100 + x * 10 + y
                            CodeList[code] = {'X': X, 'Y': Y, 'x': x, 'y': y,
                                              'Test': TestName[Y],
                                              'active': ActiveSize[Y][x][y],
                                              'AGR': AGRSize[Y][x][y],
                                              'FGR': FGRSpacing[Y][x][y],
                                              'order': Order[Y][x][y],
                                              'day-DL': day_DL(Day_DL, code)}
                            for key in self._EC:
                                CodeList[code][key] = ''

                            if datalimit == 'on':
                                if code not in Day_DL:
                                    del CodeList[code]
                            if code in CodeList:
                                for key in kwargs:
                                    if not isinstance(kwargs[key], list):
                                        kwargs[key] = [kwargs[key]]
                                    if CodeList[code][key] not in kwargs[key]:
                                        del CodeList[code]
                                        break

            statistics_path = location + '/statistics.csv'
            if os.path.isfile(statistics_path):
                with open(statistics_path, newline='') as csvfile:
                    reader = csv.DictReader(csvfile)
                    headers = reader.fieldnames
                    start = len(IVData.key(1)) + len(IVData.key(2)) + 4
                    tempfile = {element: {} for element in headers[start:]}
                    for row in reader:
                        if int(row['code']) in CodeList:
                            for key in tempfile:
                                if len(row[key]) > 0:
                                    if row[key].replace('.', '', 1).replace('-', '', 1).isdigit():
                                        if row['code'] in tempfile[key]:
                                            tempfile[key][row['code']].append(float(row[key]))
                                        else:
                                            tempfile[key][row['code']] = [float(row[key])]
                                    else:
                                        CodeList[int(row['code'])][key] = row[key]
                for key in list(tempfile.items()).copy():
                    tempfile[key[0] + '_avg'] = {}
                    tempfile[key[0] + '_std'] = {}
                    for code in tempfile[key[0]]:
                        tempfile[key[0] + '_avg'][code] = sum(tempfile[key[0]][code]) / len(tempfile[key[0]][code])
                        tempfile[key[0] + '_std'][code] = np.sqrt(sum([(element - tempfile[key[0] + '_avg'][code]) ** 2
                                                                        for element in tempfile[key[0]][code]]) /
                                                                   len(tempfile[key[0]][code]))
                for key in tempfile:
                    for code in tempfile[key]:
                        if key in CodeList[int(code)]:
                            std = round(tempfile[key + '_std'][code], 1)
                            if std == 0:
                                CodeList[int(code)][key] = '%.1f' % tempfile[key + '_avg'][code]
                            else:
                                CodeList[int(code)][key] = '%.1f&#177;%.1f' % \
                                                           (tempfile[key + '_avg'][code],
                                                            round(tempfile[key + '_std'][code], 1))
            return CodeList

        for set_key in kwargs:
            if set_key not in ['X', 'Y', 'x', 'y', 'active', 'FGR', 'AGR', 'order']:
                raise BaseException('"%s" is not in the permissible data list %s' %
                                    (set_key, ['X', 'Y', 'x', 'y', 'active', 'FGR', 'AGR', 'order']))

        DataTable = {'X': [], 'Y': [], 'x': [], 'y': [], 'code': [], 'day-DL': [],
                     'Test': [], 'active': [], 'AGR': [], 'FGR': [], 'order': []}
        for key in self._EC:
            DataTable[key] = []
        for code, info in code_generator(self._location, Xmin, Xmax, Ymin, Ymax, kwargs).items():
            DataTable['X'].append(info['X'])
            DataTable['Y'].append(info['Y'])
            DataTable['x'].append(int(str(code)[2]))
            DataTable['y'].append(int(str(code)[3]))
            DataTable['code'].append(code)
            DataTable['Test'].append(info['Test'])
            DataTable['active'].append(info['active'])
            DataTable['AGR'].append(info['AGR'])
            DataTable['FGR'].append(info['FGR'])
            DataTable['order'].append(info['order'])
            DataTable['day-DL'].append(info['day-DL'])
            for key in self._EC:
                DataTable[key].append(info[key])

        headerColor = 'grey'
        rowEvenColor = 'lightgrey'
        rowOddColor = 'white'
        HeadTitle = ['<b>Y</b>', '<b>Test</b>', '<b>Code</b>', '<b>Measurement</b>', '<b>V<sub>pt</sub></b>',
                     '<b>V<sub>pt2</sub></b>', '<b>V<sub>b</sub></b>', '<b>V<sub>r</sub></b>', '<b>Active</b>',
                     '<b>FGR</b>', '<b>AGR</b>', '<b>Order</b>']

        for key in DataTable:
            if key != 'Y':
                utils.orderfunction(DataTable['Y'].copy(), DataTable[key], 'decreasing')
        DataTable['Y'].sort(reverse=True)

        fig = go.Figure(data=
                        [go.Table(columnorder=[i for i in range(len(HeadTitle))],
                                  columnwidth=[3, 12, 6, 38, 9, 9, 9, 9, 7, 5, 5, 10],
                                  header=dict(values=HeadTitle, line_color='darkslategray', fill_color=headerColor,
                                              align=['center'], font=dict(color='white', size=15), height=30),
                                  cells=dict(values=[DataTable['Y'], DataTable['Test'], DataTable['code'],
                                                     DataTable['day-DL'], DataTable['Vpt'], DataTable['Vpt2'],
                                                     DataTable['Vb'], DataTable['Vr'], DataTable['active'], DataTable['FGR'],
                                                     DataTable['AGR'], DataTable['order']],
                                             line_color='darkslategray',
                                             fill_color=[[rowOddColor, rowEvenColor] * int(len(DataTable['Test']) / 2)],
                                             align=['center'],
                                             font=dict(color='darkslategray', size=14), height=30))])
        fig.show()

    def iv_plot(self, fig_number, label_flag='on', label='partial', title=None,
                draw=None, remove=None, colorparm=None, label_dict=None, show_number='Yes',
                hlines=None, show_index='Yes', show_date='Yes', show_gain_points='No', **kwargs):

        if 'Failure' in kwargs:
            temp_kwargs = {key: value for key, value in kwargs.items() if key != 'Failure'}
            self._key_contradiction('all', temp_kwargs)
        else:
            self._key_contradiction('all', kwargs)
        cll_set = dict()

        def col_number(CurveList):
            if colorparm is None:
                if int(len(CurveList) / 5) == 0:
                    return 1
                else:
                    return int(len(CurveList) / 5) + 1
            else:
                temp = cll_set.keys()  # list(Color_Storage.keys())
                if int(len(temp) / 5) == 0:
                    return 1
                else:
                    return int(len(temp) / 5) + 1

        def check_key(name_info, kwargs):
            '''
            info_name_list = self._key1 + self._key2
            name_info = {info_name_list[index]: element for index, element in enumerate(info)}
            '''
            flag = 0
            for key in kwargs:
                if isinstance(kwargs[key], list):
                    for element in kwargs[key]:
                        if str(element) == name_info[key]:
                            flag = 2
                            break
                    if flag == 2:
                        flag = 0
                    else:
                        flag = 1
                else:
                    if str(kwargs[key]) != name_info[key]:
                        flag = 1
                        break
            if flag == 0:
                return True
            else:
                return False

        def cll_index(j):
            if colorparm is None:
                color_index = j % len(ColorSet)
                linestyle_index = int(j / len(ColorSet))
                label_marker = 'default'
            else:
                i = CurveList[j]
                if self._data[i].info()[colorparm] not in cll_set:
                    cll_set[self._data[i].info()[colorparm]] = [len(cll_set.keys())]
                else:
                    cll_set[self._data[i].info()[colorparm]].append('more-than-one')
                j = cll_set[self._data[i].info()[colorparm]][0]
                color_index = j % len(ColorSet)
                linestyle_index = int(j / len(ColorSet))
                label_marker = 'colorparm'
            return color_index, linestyle_index, label_marker

        def label_func(i, label_marker):
            curve_number = CurveList[i]
            name_info = self._data[curve_number].info()  # {'date': 1130, 'code': 4513, ...}
            if label_marker == 'default':
                if label == 'partial':
                    if show_date == 'Yes':
                        show_info = name_info['date'] + '-' + name_info['code'] + \
                                    name_info['light'] + '-' + name_info['number']
                    elif show_date == 'No':
                        show_info = name_info['code'] + name_info['light'] + '-' + name_info['number']
                    else:
                        raise BaseException("Wrong show_date (Yes/No): %s" % show_date)

                    if show_index == 'Yes':
                        return '[%s] %s' % (curve_number, show_info)
                    elif show_index == 'No':
                        return show_info
                    else:
                        raise BaseException("Wrong show_index (Yes/No): %s" % show_index)
                elif label == 'full':
                    return self._data[curve_number].name()
                else:
                    raise BaseException("Wrong name settings(partial/full): %s" % label)
            else:
                if len(cll_set[self._data[curve_number].info()[colorparm]]) == 1:  # if len(Color_Storage[name_info[colorparm]]) <= 1:
                    if label_dict is None:
                        return '%s' % (colorparm + ': ' + str(name_info[colorparm]))
                    else:
                        return label_dict[name_info[colorparm]]
                else:
                    return None

        def CurveNumber():
            if colorparm is None:
                return len(CurveList)
            else:
                return len(cll_set.keys())

        plt.figure(fig_number)
        ColorSet = ['dodgerblue', 'yellowgreen', 'goldenrod', 'darkviolet',
                    'darkorange', 'brown', 'b', 'hotpink', 'fuchsia', 'g', 'royalblue', 'tomato',
                    'purple', 'olive', 'darkgreen']
        LineStyleSet = ['-', '-.', '--', ':', (0, (3, 1, 1, 1)), (0, (5, 5)),
                        (0, (3, 1, 1, 1, 1, 1)), (0, (3, 5, 1, 5)), (0, (5, 10))]

        if draw is not None and remove is not None:
            raise BaseException("Cannot use draw_list and except_list simultaneously")

        CurveList = []
        IV = self.read()
        for i in range(len(IV)):
            if draw is None:
                if remove is None:
                    if check_key(self._data[i].info(), kwargs):
                        CurveList.append(i)
                else:
                    if i not in remove:
                        if check_key(self._data[i].info(), kwargs):
                            CurveList.append(i)
            else:
                if i in draw:
                    if check_key(self._data[i].info(), kwargs):
                        CurveList.append(i)

        if len(CurveList) == 0:
            raise BaseException("No such curve!")
        if len(CurveList) >= len(ColorSet) * len(LineStyleSet):
            raise BaseException("Too many curves (%s) (Max:%s)" % (len(IV), len(ColorSet) * len(LineStyleSet)))
        for j in range(len(CurveList)):
            color_index, linestyle_index, label_marker = cll_index(j)
            if color_index >= len(ColorSet) or linestyle_index >= len(LineStyleSet):
                print('Too many curves (max:%s): %s' % (len(ColorSet) * len(LineStyleSet), len(CurveList)))
            plt.plot(IV[CurveList[j]].X, abs(IV[CurveList[j]].Y),
                     color=ColorSet[color_index],
                     linestyle=LineStyleSet[linestyle_index],
                     label=label_func(j, label_marker))
        if label_flag == 'on':
            plt.legend(loc='best', ncol=col_number(CurveList))
        elif label_flag != 'off':
            raise BaseException("Wrong label_flag (on/off): %s" % label_flag)
        if title is not None:
            if show_number == 'Yes':
                plt.title(title + ' (%s)' % CurveNumber(), fontsize=15)
            elif show_number == 'No':
                plt.title(title, fontsize=15)
            else:
                raise BaseException("Wrong show_number (Yes/No): %s" % show_number)
        else:
            plt.title('%s (%s)' % (kwargs, CurveNumber()), fontsize=15)
        if hlines is not None:
            if not isinstance(hlines, list):
                hlines = [hlines]
            for y in hlines:
                plt.hlines(y=y, xmin=-60, xmax=0, color='k', linestyle='-.')
        plt.grid()
        plt.yscale('log')
        plt.xlabel('Voltage (V)')
        plt.ylabel('Current (A)')


class Statistics(IVDataBase):
    def __init__(self, location, failure_check='No'):
        """
        1. 目的：這是用以統計製程與電性數據的工具，可以在此新增、移除、修改與繪製統計資料，例如擊穿電壓與元件座標、FGR spacing等。
        2. 輸入：因為統計表單是針對所有數據統計，而不是針對特定幾種元件的數據，所以參數僅有資料庫路徑。
        3. 輸出：執行後會在資料庫的最上層目錄中自動生成統計數據 statistics.csv，其數據列表可由 IVDataBase 中的 self._EC 修改。
        :param location: 資料庫路徑，例如 '/Users/Ethan/PycharmProjects/Python/Avalanche Photodiode/Experiment/IV'。
        :param show: 顯示目錄與檔案之成功讀取訊息，用以排除故障。
        """
        super().__init__(location)
        self._location = location
        self._statistics = location + '/statistics.csv'
        self._characteristics = ['X', 'Y', 'x', 'y'] + self._EC
        self._state = ['Failure', 'RelativeGain', 'Vbm', 'Vm1', 'GainReference']
        self._data_list = IVData.key('all') + self._characteristics + self._state

        # 檢查儲存目錄是否已有統計資料，若有則檢查元件清單是否一致，否則新建統計資料
        if not os.path.isfile(self._statistics):
            with open(self._statistics, 'w', newline='') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=self._data_list)
                writer.writeheader()
                for data in self._data:
                    # 初始化數據
                    for key in self._characteristics:
                        if key not in data.info():
                            data.info()[key] = ''
                    # 儲存數據
                    writer.writerow(data.info())
        else:
            # 檢查抬頭是否相同
            flag = 0
            with open(self._statistics, newline='') as csvfile:
                reader = csv.DictReader(csvfile)
                headers = reader.fieldnames
                if headers != self._data_list:
                    flag = 1
            if flag == 1:
                ans = input("Are you sure you want to modify the data label? It's irreversible. (Y/N):  ")
                if ans == 'Y':
                    tempfile = NamedTemporaryFile(mode='wt', delete=False)
                    with open(self._statistics, newline='') as csvfile, tempfile:
                        reader = csv.DictReader(csvfile)
                        writer = csv.DictWriter(tempfile, fieldnames=self._data_list)
                        writer.writeheader()
                        for row in reader:
                            for key in self._data_list:
                                if key not in row:
                                    row[key] = ''
                            writer.writerow(row)
                    shutil.move(tempfile.name, self._statistics)

            # 逐一確認每個量測數據是否都已被存在記錄簿中
            for data in self._data:

                flag = 0  # 檢查元件資料是否不同於所有已存資料。
                with open(self._statistics, newline='') as csvfile:
                    reader = csv.DictReader(csvfile)
                    row_count = sum(1 for row in reader)

                with open(self._statistics, newline='') as csvfile:
                    reader = csv.DictReader(csvfile)
                    for row in reader:
                        for key in IVData.key(1):
                            if row[key] != data.info()[key]:
                                flag += 1
                                break

                if flag == row_count:
                    for key in self._characteristics:
                        if key not in data.info():
                            data.info()[key] = ''
                    with open(self._statistics, 'a', newline='') as csvfile:
                        writer = csv.DictWriter(csvfile, fieldnames=self._data_list)
                        writer.writerow(data.info())

        # 更新 failure 狀態
        if failure_check == 'Yes':
            tempfile = NamedTemporaryFile(mode='wt', delete=False)
            with open(self._statistics, newline='') as csvfile, tempfile:
                reader = csv.DictReader(csvfile)
                writer = csv.DictWriter(tempfile, fieldnames=self._data_list)
                writer.writeheader()
                for row in reader:
                    if row['Failure'] not in ['Yes', 'No']:
                        if row['Vpt'] == 'null' and row['Vb'] == 'null' \
                                and row['Vpt2'] == 'null' and row['Vr'] == 'null':
                            row['Failure'] = 'Yes'
                    writer.writerow(row)
            shutil.move(tempfile.name, self._statistics)

    def add(self, **kwargs):
        """
        用以輸入元件參數，例如針對 3612 元件在 11/24 測量的第 2 組光電流數據，新增其擊穿電壓與崩潰電壓。
        :param kwargs: (date, code, light, number) 與 (Vpt, Vpt2, Vr, Vb)，其中 Vr 為崩潰前電流突然異常減小(reduced)的電壓
        :return: 在 statistics.csv 更新數據

        from IVanalysis import Statistics as STA
        DB = STA(data_path)

        # 針對 11/24-2732L-1 元件數據加上 Vpt=-19 參數。
        DB.add(date=1124, code=2732, light='L', number=1, Vpt=-19)

        # 針對 2732 元件的所有量測數據，加上 Vr='null' 參數（null表示沒有這個參數）。
        DB.add(code=2732, Vr='null')

        # 針對在 11/24 測量的所有 1521 元件數據，加上相同的Vpt、Vb、Vpt2與Vr。
        DB.add(date=1124, code=1521, Vpt=-19.43, Vb=-50, Vpt2='null', Vr='null')
        """
        if 'Failure' in kwargs and ('Vpt' in kwargs or 'Vpt2' in kwargs or 'Vb' in kwargs or 'Vr' in kwargs):
            if kwargs['Failure'] == 'Yes' and (kwargs['Vpt'] != 'null' and kwargs['Vpt2'] != 'null'
                                               and kwargs['Vb'] != 'null' and kwargs['Vr'] != 'null'):
                ans = input("Are you sure [%s] is failure? (Yes/No): " % kwargs)
                if ans not in ['Yes', 'No']:
                    raise BaseException("Wrong answer (Yes/No): %s" % ans)
                else:
                    if ans == 'No':
                        raise BaseException("Check failure state.")

            flag = 0
            if kwargs['Failure'] == 'No':
                for key in kwargs:
                    if key in ['Vpt', 'Vpt2', 'Vb', 'Vr'] and kwargs[key] != 'null':
                        flag = 1
                        break
                if flag == 0:
                    raise BaseException("Check failure state.")

        if 'Failure' in kwargs:
            if kwargs['Failure'] == 'Yes':
                for key in kwargs:
                    if key in ['Vpt', 'Vpt2', 'Vb', 'Vr'] and kwargs[key] != 'null':
                        answer = input('Are you sure device failure occurred (Yes/No)? : ')
                        if answer == 'Yes' or answer == 'No':
                            kwargs['Failure'] = answer
                            break
                        else:
                            raise BaseException("Wrong answer (Yes/No): %s" % answer)

        for key in kwargs:
            if key not in self._characteristics[4:] + IVData.key(1) + self._state + ['Y'] and key != 'Forced':
                raise BaseException("%s is not in the permissible list: %s" %
                                    (key, IVData.key(1) + self._characteristics[4:] + self._state + ['Y']))

        tempkey = []
        for key in kwargs:
            if key in IVData.key(1):
                tempkey.append(key)

        def gate(row, kwargs, tempkey):
            flag = 0
            for key in tempkey:
                if row[key] != str(kwargs[key]):
                    flag = 1

            # 檢查 Y 座標
            if 'Y' in kwargs:
                if not isinstance(kwargs['Y'], list):
                    kwargs['Y'] = [kwargs['Y']]
                if int(row['code'][1]) not in kwargs['Y']:
                    flag = 1
            if flag == 0:
                return True
            else:
                return False

        tempfile = NamedTemporaryFile(mode='wt', delete=False)
        with open(self._statistics, 'r', newline='') as csvfile, tempfile:
            reader = csv.DictReader(csvfile)
            writer = csv.DictWriter(tempfile, fieldnames=self._data_list)
            writer.writeheader()
            flag = 0
            for row in reader:
                if gate(row, kwargs, tempkey):
                    flag = 1
                    for key, value in kwargs.items():
                        if key not in ['Y', 'Forced']:
                            if row[key] != '' and row[key] != str(value):
                                if 'Forced' in kwargs:
                                    if kwargs['Forced'] == 'Yes':
                                        row[key] = value
                                else:
                                    ans = input('[%s-%s%s-%s] Change %s = %s to %s? (Y/N):  '
                                                % (row['date'], row['code'], row['light'], row['number'], key, row[key],
                                                   value))
                                    if ans == 'Y':
                                        row[key] = value
                                    elif ans not in ['Y', 'N']:
                                        raise BaseException("Wrong answer (Y/N): %s" % ans)
                            else:
                                row[key] = value
                    writer.writerow(row)
                else:
                    writer.writerow(row)
            if flag == 0:
                if len(tempkey) == 0 and 'Y' in kwargs:
                    raise BaseException('Y=%s not found' % kwargs['Y'])
                else:
                    raise BaseException('%s not found' % tempkey)
        shutil.move(tempfile.name, self._statistics)

    def delete(self, date, code, light, number):
        """
        清空特定元件量測數據的所有參數
        :param date: 量測日期，格式為 mmdd，如1124
        :param code: 元件座標代碼，如1712
        :param light: 照光條件（'L'/'D'）
        :param number: 實驗編號
        :return: 刪除 statistics.csv 資料庫中的某次元件量測的參數（Vpt,Vpt2,Vr,Vb)
        """
        tempfile = NamedTemporaryFile(mode='wt', delete=False)
        with open(self._statistics, 'r', newline='') as csvfile, tempfile:
            reader = csv.DictReader(csvfile)
            writer = csv.DictWriter(tempfile, fieldnames=self._data_list)
            writer.writeheader()
            for row in reader:
                if row['date'] == str(date) and row['code'] == str(code) \
                        and row['light'] == light and row['number'] == str(number):
                    for key in row:
                        if key in self._characteristics[4:]:
                            row[key] = ''
                writer.writerow(row)
        shutil.move(tempfile.name, self._statistics)
        print('[%s-%s%s-%s] is deleted' % (date, code, light, number))

    def print(self, **kwargs):
        """
        印出特定元件的現有參數列表。
        :param kwargs: 為 (date, code, light, number) 的組合。
        :return: 印出其 Vpt, Vpt2, Vr, Vb 參數值。

        DB.print(code=3712)  # 印出 3712 元件的所有參數
        DB.print(date=1124, code=1712)  # 印出在 1124 量測的所有 3712 元件參數。
        DB.print(date=1124, code=1712, light='L', number=1)  # 印出 11/24 測量的 1712L-1.TXT 數據之參數。
        """
        for key in kwargs:
            if key not in IVData.key(1):
                raise BaseException("%s is not in permissible list %s" % (key, IVData.key(1)))
        if 'number' in kwargs and ('date' not in kwargs or 'code' not in kwargs or 'light' not in kwargs):
            raise BaseException("Only number indication is not enough.")
        with open(self._statistics, 'r', newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                temprow = {}
                flag = 0
                for key in kwargs:
                    if row[key] != str(kwargs[key]):
                        flag = 1
                if flag == 0:
                    for key, value in row.items():
                        if key in self._characteristics[4:] + ['Failure']:
                            temprow[key] = value
                    print('[%s-%s%s-%s]  %s' % (row['date'], row['code'], row['light'], row['number'],
                                                ['%s: %s' % (key, value) for key, value in temprow.items()]))

    def plot(self, x, y, avg='No', std='No', xnull='No', ynull='No', absvalue='Yes', str_number=7, **kwargs):
        """
        挑選兩參數互相做圖，觀察關係。
        :param x: 橫軸參數，如元件主動區直徑（Ex: x='active'）。
        :param y: 縱軸參數，如元件擊穿電壓（Ex: y='Vpt'）。
        :param avg: 取出該在相同x值下，各種y值的平均值。
        :param std: 取出該在相同x值下，各種y值的標準差。
        :param null: 是否允許讀取 'null'（屬於字串，string）作為y值，還在開發中，預設為'no'。遇到參數為 'null' 之數據自動跳過。
        :param absvalue: 是否自動將y值取絕對值（倘若y值為數字，而非字串）。
        :param kwargs: 元件範圍條件，例如針對 11/24 測量的 1512 元件，或是針對 active=800 的元件，此外也能加上 Ymin 與 Ymax 的條件。
        :return: 回傳字典 dictionary 資料結構，如 {'X': [1,2,3,4], 'Y': ['1e-10, 1.2e-10, 3e-10, 4e-10]}

        範例一：針對 6 <= Y <= 7 的元件數據，匯出擊穿電壓對主動區直徑的關係圖，並附上其平均值。
        Data = DB.plot(x='active', y='Vpt', avg='Yes', Ymin=6, Ymax=7)
        plt.figure(1)
        plt.plot(Data['X'], Data['Y'], 'bo', linestyle='none', label='Vpt')
        plt.plot(Data['X_avg'], Data['Y_avg'], '-.b')
        plt.grid()
        plt.show()

        範例二：針對 11/24 量測的主動區直徑 100um 元件數據，繪製元件的 X 晶圓座標與崩潰電壓的關係，並附上平均值與標準差。
        DB.plot(x='X', y='Vb', avg='Yes', std='Yes', active=100, date=1124)
        plt.figure(1)
        plt.plot(Data['X'], Data['Y'], 'bo', linestyle='none', label='Vb')
        plt.errorbar(Data['X_avg'], Data['Y_avg'], Data['Y_std'], '-.b')
        plt.grid()
        plt.show()
        """
        coordinate = ['Xpos', 'Ypos']
        OtherKeys = ['Ymin', 'Ymax', 'contact', 'VbmRatio']
        Xmin = 1
        Ymin = 1
        Xmax = 13
        Ymax = 9
        X_List = np.arange(Xmin, Xmax + 1, 1)
        Y_List = np.arange(Ymin, Ymax + 1, 1)
        xmax = {posX: {1: 2, 2: 3, 3: 2, 4: 3, 5: 2, 6: 1, 7: 3, 8: 3, 9: 3} for posX in X_List}
        ymax = {posX: {1: 4, 2: 3, 3: 3, 4: 3, 5: 4, 6: 3, 7: 3, 8: 3, 9: 3} for posX in X_List}

        if x not in self._data_list and x not in coordinate:
            raise BaseException('xlabel %s is not in the permissible list %s' % (x, self._data_list + coordinate))
        elif y not in self._data_list and y not in coordinate:
            raise BaseException('ylabel %s is not in the permissible list %s' % (y, self._data_list + coordinate))
        for key in kwargs:
            if key not in self._data_list + OtherKeys:
                temp = self._data_list + OtherKeys
                raise BaseException('key %s is not in the permissible list %s' % (key, temp))

        def gate(row, kwargs):
            flag = 0
            contact_list = {'design1': '4', 'design2': '3', 'design3': '2'}
            if 'contact' in kwargs and int(row['Y']) != 1:
                return False
            for key in kwargs:
                if key not in OtherKeys:
                    if "'" not in row[key]:
                        if row[key] != str(kwargs[key]):
                            flag = 1
                    else:
                        if row[key].split("'")[1] != str(kwargs[key]):
                            flag = 1
                elif key == 'Ymin' and int(row['Y']) < kwargs[key]:
                    flag = 1
                elif key == 'Ymax' and int(row['Y']) > kwargs[key]:
                    flag = 1
                elif key == 'contact':
                    if kwargs['contact'] not in list(contact_list.keys()):
                        raise BaseException("Contact value is wrong (design1/design2/design3): %s" %
                                            kwargs['contact'])
                    else:
                        if int(row['y']) != int(contact_list[kwargs['contact']]):
                            flag = 1
            if y == 'RelativeGain' and row['Vbm'] != '':
                if row['Vr'] not in ['', 'null']:
                    if 'VbmRatio' in kwargs:
                        if abs(float(row['Vr'])) < kwargs['VbmRatio'] * abs(float(row['Vbm'])) or \
                                kwargs['VbmRatio'] * abs(float(row['Vbm'])) < abs(float(row['Vm1'])):
                            flag = 1
                    else:
                        raise BaseException("'VbmRatio' is needed when plotting 'RelativeGain'")
                else:
                    if 'VbmRatio' in kwargs:
                        if kwargs['VbmRatio'] * abs(float(row['Vbm'])) < abs(float(row['Vm1'])):
                            flag = 1
                    else:
                        raise BaseException("'VbmRatio' is needed when plotting 'RelativeGain'")
            if flag == 0:
                return True
            else:
                return False

        X = []
        Y = []
        str_flag = 0
        xnull_flag = 0
        with open(self._statistics, 'r', newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                if ynull == 'No':
                    if row[y] != 'null':
                        if x not in coordinate:
                            if str_flag == 0 and not row[x].replace('.', '', 1).replace('-', '', 1).isdigit() \
                                    and xnull == 'Yes':
                                str_flag = 1
                        if gate(row, kwargs) and row[y] != '':
                            if x not in coordinate:
                                if row[x].replace('.', '', 1).replace('-', '', 1).isdigit():
                                    X.append(float(row[x]))
                                    xnull_flag = 0
                                else:
                                    if xnull == 'Yes':
                                        X.append(str_number)
                                    else:
                                        xnull_flag = 1
                            else:
                                dx = 1 / (xmax[1][int(row['Y'])] + 1)
                                dy = 1 / (ymax[1][int(row['Y'])] + 1)
                                d = {'X': dx, 'Y': dy}
                                X.append(float(row[x[0]]) - 1 + d[x[0]] * float(row[x[0].lower()]))  # 'Xpos'[0] = 'X'
                            if xnull_flag == 0:
                                if row[y].replace('.', '', 1).replace('-', '', 1).isdigit():
                                    if absvalue == 'Yes':
                                        Y.append(abs(float(row[y])))
                                    else:
                                        Y.append(float(row[y]))
                                else:
                                    Y.append(row[y])
                else:
                    if str_flag == 0 and not row[x].replace('.', '', 1).replace('-', '', 1).isdigit():
                        str_flag = 1
                    if gate(row, kwargs) and row[y] != '':
                        if row[x].replace('.', '', 1).replace('-', '', 1).isdigit():
                            X.append(float(row[x]))
                        else:
                            if xnull == 'No':
                                break
                            else:
                                X.append(str_number)
                        if row[y].replace('.', '', 1).replace('-', '', 1).isdigit():
                            if absvalue == 'Yes':
                                Y.append(abs(float(row[y])))
                            else:
                                Y.append(float(row[y]))
                        else:
                            Y.append(str_number)
        if str_flag == 0:
            X, Y = utils.orderfunction(X.copy(), Y.copy(), 'increasing')
        Data = {'X': np.array(X), 'Y': np.array(Y)}
        if avg == 'Yes':
            temp = {}
            for i, x in enumerate(X):
                if x not in temp:
                    temp[x] = [Y[i]]
                else:
                    temp[x].append(Y[i])
            Data['X_avg'] = np.array([key for key in temp])
            Data['Y_avg'] = np.array([sum(temp[key]) / len(temp[key]) for key in temp])
        else:
            if avg != 'No':
                raise BaseException('Wrong avg value (Yes/No):  %s' % avg)
        if std == 'Yes':
            temp = {}
            for i, x in enumerate(X):
                if x not in temp:
                    temp[x] = [Y[i]]
                else:
                    temp[x].append(Y[i])
            Data['X_std'] = np.array([key for key in temp])
            temp_std = {key: [(value - sum(temp[key]) / len(temp[key])) ** 2 for value in temp[key]] for key in temp}
            Data['Y_std'] = np.array([np.sqrt(sum(temp_std[key]) / len(temp_std[key])) for key in temp])
        return Data

    def check(self, **kwargs):
        """
        給定 (date, code, light, number)，針對該元件量測結果，檢查其電性參數是否已登錄齊全。
        :param kwargs: 可選 (date, code, light, number, Y, active, FGR, AGR, order, depth)。
        :return: 印出結果。

        DB.check(code=1521)  # 檢查代碼為 1512 的所有元件的參數登錄狀況。
        DB.check(date=1124, code=1721, light='L', number=1)  # 檢查此量測數據的參數登錄狀況
        """
        Remaining = dict()
        found_flag = 0
        for key in kwargs:
            if key not in IVData.key('all'):
                raise BaseException("%s is not in permissible list %s" % (key, IVData.key('all')))
        if 'number' in kwargs and ('date' not in kwargs or 'code' not in kwargs or 'light' not in kwargs):
            raise BaseException("Only number indication is not enough.")
        with open(self._statistics, 'r', newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                flag = 0
                for key in kwargs:
                    if row[key] != str(kwargs[key]):
                        flag = 1
                if flag == 0:
                    found_flag = 1
                    for key in row:
                        if key in self._characteristics[4:] + ['Failure'] and row[key] == '':
                            name = row['date'] + '-' + row['code'] + row['light'] + '-' + row['number']
                            if key not in Remaining:
                                Remaining[key] = [name]
                            else:
                                Remaining[key].append(name)
        if found_flag == 0:
            raise BaseException("%s not found in the IV database" % kwargs)
        else:
            if len(list(Remaining.items())) == 0:
                print('%s is completed!' % kwargs)
            else:
                for key in Remaining:
                    for item in Remaining[key]:
                        print('[%s]   %s' % (key, item))

    def failure(self, fig_number, fig_title, draw, draw_color=None):
        def coordinate(row):
            Xmin = 1
            Ymin = 1
            Xmax = 13
            Ymax = 9
            X_List = np.arange(Xmin, Xmax + 1, 1)
            xmax = {posX: {1: 2, 2: 3, 3: 2, 4: 3, 5: 2, 6: 1, 7: 3, 8: 3, 9: 3} for posX in X_List}
            ymax = {posX: {1: 4, 2: 3, 3: 3, 4: 3, 5: 4, 6: 3, 7: 3, 8: 3, 9: 3} for posX in X_List}
            dx = 1 / (xmax[1][int(row['Y'])] + 1)
            dy = 1 / (ymax[1][int(row['Y'])] + 1)
            d = {'X': dx, 'Y': dy}
            Xpos = float(row['X']) - 1 + d['X'] * float(row['x'])
            Ypos = float(row['Y']) - 1 + d['Y'] * float(row['y'])
            return Xpos, Ypos

        Dataset = {'Failure': [], 'Vpt2': [], 'Vr': [], 'Measured': [], 'Leakage': []}
        Measured_Code = []
        with open(self._statistics, 'r', newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                if row['Failure'] == 'Yes':
                    Dataset['Failure'].append(coordinate(row))
                if row['Vpt2'] not in ['', 'null']:
                    Dataset['Vpt2'].append(coordinate(row))
                if row['Vr'] not in ['', 'null']:
                    Dataset['Vr'].append(coordinate(row))
                if row['code'] not in Measured_Code:
                    Measured_Code.append(row['code'])
                    Dataset['Measured'].append(coordinate(row))
                if row['Failure'] == 'Leakage':
                    Dataset['Leakage'].append(coordinate(row))

        for key, value in Dataset.items():
            Dataset[key] = list(set(value))
        x = {key: [xy[0] for xy in Dataset[key]] for key in Dataset.keys()}
        y = {key: [xy[1] for xy in Dataset[key]] for key in Dataset.keys()}

        area = np.pi * 3
        Color = {'Failure': 'red', 'Vpt2': 'green', 'Vr': 'b', 'Measured': 'royalblue', 'Leakage': 'red'}
        Marker = {'Failure': 'x', 'Vpt2': 's', 'Vr': 'v', 'Measured': 'o', 'Leakage': '^'}
        size = {'Failure': 9 * area, 'Vpt2': 3 * area, 'Vr': 3 * area, 'Measured': 9 * area, 'Leakage': 3 * area}
        Transparency = {'Failure': 1, 'Vpt2': 1, 'Vr': 1, 'Measured': 0.3, 'Leakage': 1}
        Label = {'Failure': 'S/C', 'Vpt2': r'V$_{pt2}$', 'Vr': r'V$_r$', 'Measured': 'Measured', 'Leakage': 'Leakage'}

        # Plot
        plt.figure(fig_number, figsize=(9, 8))
        plt.scatter(x['Measured'], y['Measured'], marker=Marker['Measured'], c=Color['Measured'],
                    alpha=Transparency['Measured'], s=size['Measured'], label='Measured')
        for element in draw:
            if draw_color is not None:
                Color[element] = draw_color
            plt.scatter(x[element], y[element], marker=Marker[element], c=Color[element],
                        alpha=Transparency[element], s=size[element], label=Label[element])
        plt.title(fig_title)
        for i in range(7):
            plt.vlines(x=i, ymin=0, ymax=9, color='grey')
        for i in range(9):
            plt.hlines(y=i, xmin=0, xmax=6, color='grey')
        number = len(draw) + 1
        plt.legend(loc='center', ncol=number, bbox_to_anchor=(0.5, 1.1))
        plt.xlim((0, 6))
        plt.ylim((0, 9))
        plt.xticks([0.5 + i for i in range(6)], [i for i in range(1, 7)])
        plt.yticks([0.5 + i for i in range(9)], [i for i in range(1, 10)])
        plt.xlabel('X')
        plt.ylabel('Y')

    def gain(self, date, code, number, vbm, vm1, ratio=0.9, forced='No'):
        if vbm >= 0 or vm1 >= 0:
            raise BaseException("Wrong vb or vm1 (negative number): %s, %s" % (vbm, vm1))
        IV = IVDataBase(self._location, date=date, code=code, light='L', number=number).read()[0]
        Im1 = utils.find(IV.X, abs(IV.Y), vm1, 'log')
        Im = utils.find(IV.X, abs(IV.Y), ratio * vbm, 'log')
        relative_gain = Im / Im1
        self.add(date=date, code=code, light='L', number=number, GainReference='Yes')
        self.add(code=code, Vbm=vbm, Vm1=vm1, RelativeGain=relative_gain, Forced=forced)
        print('[%s]: %.2f' % (code, relative_gain))

    def new_gain(self, ratio):
        IVList = dict()
        with open(self._statistics, 'r', newline='') as IVfile:
            IVreader = csv.DictReader(IVfile)
            for row in IVreader:
                if row['GainReference'] == 'Yes' and row['code'] not in IVList:
                    IVList[row['code']] = IVDataBase(self._location, date=row['date'], code=row['code'],
                                                     light='L', number=row['number']).read()[0]

        tempfile = NamedTemporaryFile(mode='wt', delete=False)
        with open(self._statistics, 'r', newline='') as csvfile, tempfile:
            reader = csv.DictReader(csvfile)
            writer = csv.DictWriter(tempfile, fieldnames=self._data_list)
            writer.writeheader()
            for row in reader:
                if row['RelativeGain'] != '':
                    vm1 = float(row['Vm1'])
                    vbm = float(row['Vbm'])
                    if row['code'] in IVList:
                        IV = IVList[row['code']]
                        Im1 = utils.find(IV.X, abs(IV.Y), vm1, 'log')
                        Im = utils.find(IV.X, abs(IV.Y), ratio * vbm, 'log')
                        print('[%s] %.2f  >>>  %.2f' % (row['code'], float(row['RelativeGain']), Im / Im1))
                        row['RelativeGain'] = Im / Im1
                        writer.writerow(row)
                    else:
                        raise BaseException("[%s] GainReference is unknown" % row['code'])
                else:
                    writer.writerow(row)
        shutil.move(tempfile.name, self._statistics)