import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.pylab as pylab
import utils
import os
import time
import csv
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


class IVDataBase(object):
    def __init__(self, location, show='no', **kwargs):
        """
        只要資料庫符合 date/code/codeD-1.TXT 格式，就能夠藉已知參數快速讀取彙整資料。
        :param location:資料庫目錄。例如 /user/..../IV。而 IV 資料夾中需有實驗日期目錄，如 /IV/1120 為11月20日的量測結果。
        :param show: 顯示讀取資料夾與檔案的成功訊息，用以排除故障。
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
        def set_location(code_folder, filename):
            return code_folder + '/' + filename

        def set_name(m, d, filename):
            if 'D' in filename:
                Code = int(filename.split("D")[0])
            elif 'L' in filename:
                Code = int(filename.split("L")[0])
            else:
                raise BaseException('"D" and "L" are not in the file name. Please check the filename list.')
            Info = self.code_info(Code, show='No')
            return str(m) + '/' + str(d) + '-' + filename.split(".")[0] + '-%s' % [value for key, value in Info.items()]

        def filename_check(filename, code, path):
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
                print('Check %s folder: %s' % (code, path + '/' + code + '/' + filename))
                return False

        def number_format():
            return '-' + str(kwargs['number']) + '.TXT'

        def check_key(default_keys):
            key2_list = []
            for key in kwargs:
                if key in default_keys:
                    key2_list.append(key)
            return key2_list

        # 預設電性
        self._EC = ['Vpt', 'Vb', 'Vpt2', 'Vr']

        # 檢查 key 是否衝突
        self._location = location
        self._key1 = ['date', 'code', 'light', 'number']
        self._key2 = ['active', 'FGR', 'AGR', 'order', 'depth', 'Y']
        self._key_contradiction(kwargs)

        # 開始讀檔
        os.chdir(self._location)
        self._file_location = []  # 單個/多個檔案位址
        self._filename_list = []  # 單個/多個檔名
        self._key2_foundfolder = []  # 符合 key2 的資料夾
        self._DateDirectory = [name for name in os.listdir(".") if os.path.isdir(name) and len(name) == 4]
        if 'date' in kwargs:
            self._month = int(kwargs['date'] / 100)
            self._day = kwargs['date'] % 100
            if str(kwargs['date']) in self._DateDirectory:
                self._DateFolder = location + '/' + str(kwargs['date'])
                if show == 'Yes':
                    print('Successfully open %s/%s folder' % (self._month, self._day))
                os.chdir(self._DateFolder)
                self._CodeDirectory = [name for name in os.listdir(".") if os.path.isdir(name) and name.isdigit()]
                if 'code' in kwargs:
                    if str(kwargs['code']) in self._CodeDirectory:
                        self._CodeFolder = self._DateFolder + '/' + str(kwargs['code'])
                        if show == 'Yes':
                            print('Successfully read %s in %s/%s folder' % (kwargs['code'], self._month, self._day))
                        Filename = [name for name in os.listdir(self._CodeFolder)
                                    if filename_check(name, str(kwargs['code']), self._DateFolder)]
                        if 'light' in kwargs:
                            on_off_file = {'L': [], 'D': []}
                            if kwargs['light'] in on_off_file:
                                for file in Filename:
                                    if kwargs['light'] in file:
                                        on_off_file[kwargs['light']].append(file)
                                if len(on_off_file[kwargs['light']]) == 0:
                                    raise BaseException("No photo-IV data in (%s/%s, %s) folder" %
                                                        (self._month, self._day, kwargs['code']))
                                else:
                                    if 'number' in kwargs:
                                        for file in on_off_file[kwargs['light']]:
                                            if number_format() in file:
                                                self._file_location.append(set_location(self._CodeFolder, file))
                                                self._filename_list.append(set_name(self._month, self._day, file))
                                                if show == 'Yes':
                                                    print('Successfully load %s' % file)
                                    else:
                                        self._file_location = [set_location(self._CodeFolder, file)
                                                               for file in on_off_file[kwargs['light']]]
                                        self._filename_list = [set_name(self._month, self._day, file) for file in on_off_file[kwargs['light']]]
                                        if show == 'Yes':
                                            print('Successfully load all files in %s at %s/%s' %
                                                  (kwargs['code'], self._month, self._day))
                                            print(self._filename_list)
                            else:
                                raise BaseException('Wrong light status: "%s" (L/D)' % kwargs['light'])
                        else:
                            self._file_location = [set_location(self._CodeFolder, file) for file in Filename]
                            self._filename_list = [set_name(self._month, self._day, file) for file in Filename]
                    else:
                        raise BaseException("%s files not found at %s/%s folder!" %
                                            (kwargs['code'], self._month, self._day))
                else:
                    key2_list = check_key(self._key2)
                    if len(key2_list) == 0:
                        for code_subdir in self._CodeDirectory:
                            self._CodeFolder = self._DateFolder + '/' + code_subdir
                            Filename = [name for name in os.listdir(self._CodeFolder)
                                        if filename_check(name, str(code_subdir), self._DateFolder)]
                            if 'light' in kwargs:
                                on_off_file = {'L': [], 'D': []}
                                if kwargs['light'] in on_off_file:
                                    for file in Filename:
                                        if kwargs['light'] in file:
                                            on_off_file[kwargs['light']].append(file)
                                    if len(on_off_file[kwargs['light']]) == 0:
                                        raise BaseException("No %s-IV data in (%s/%s, %s) folder" %
                                                            (kwargs['light'], self._month, self._day, code_subdir))
                                    else:
                                        if 'number' in kwargs:
                                            for file in on_off_file[kwargs['light']]:
                                                if number_format() in file:
                                                    self._file_location.append(set_location(self._CodeFolder, file))
                                                    self._filename_list.append(set_name(self._month, self._day, file))
                                                    if show == 'Yes':
                                                        print('Successfully load %s' % file)
                                        else:
                                            for file in on_off_file[kwargs['light']]:
                                                self._file_location.append(set_location(self._CodeFolder, file))
                                                self._filename_list.append(set_name(self._month, self._day, file))
                                else:
                                    raise BaseException('Wrong light status: "%s" (L/D)' % kwargs['light'])
                            else:
                                for file in Filename:
                                    self._file_location.append(set_location(self._CodeFolder, file))
                                    self._filename_list.append(set_name(self._month, self._day, file))
                    else:
                        code_info = {codeFolder: self.code_info(int(codeFolder), show='No')
                                     for codeFolder in self._CodeDirectory}
                        pass_codeFolder = {key: [] for key in key2_list}
                        for folder in code_info:
                            for key in key2_list:
                                if code_info[folder][key] == kwargs[key]:
                                    pass_codeFolder[key].append(folder)

                        # 檢查是否有完全沒找到的 key2
                        NotFound = dict()
                        for key in pass_codeFolder:
                            if len(pass_codeFolder[key]) == 0:
                                NotFound[key] = kwargs[key]
                        if len(NotFound.keys()) > 0:
                            raise BaseException('%s not found at %s/%s folder.' % (NotFound, self._month, self._day))
                        if show == 'Yes':
                            print('Search result: %s' % pass_codeFolder)

                        # 找出交集
                        FoundFolder = set(pass_codeFolder[list(pass_codeFolder.keys())[0]])
                        if len(pass_codeFolder) > 1:
                            for element in list(pass_codeFolder.keys())[1:]:
                                FoundFolder.intersection_update(set(pass_codeFolder[element]))
                        self._key2_foundfolder = list(FoundFolder)
                        self._key2_foundfolder.sort()  # 固定 code 資料夾順序，否則每次開啟都會不同。
                        if len(self._key2_foundfolder) > 0:
                            if show == 'Yes':
                                print('Search result in %s/%s folder for %s : %s ' %
                                      (self._month, self._day,
                                       {key: kwargs[key] for key in key2_list}, self._key2_foundfolder))
                        else:
                            raise BaseException('Set intersection is null.')

                        for code in self._key2_foundfolder:
                            self._CodeFolder = self._DateFolder + '/' + code
                            if show == 'Yes':
                                print('Successfully read %s in %s/%s folder' % (code, self._month, self._day))
                            Filename = [name for name in os.listdir(self._CodeFolder)
                                        if filename_check(name, str(code), self._DateFolder)]
                            if 'light' in kwargs:
                                on_off_file = {'L': [], 'D': []}
                                if kwargs['light'] in on_off_file:
                                    for file in Filename:
                                        if kwargs['light'] in file:
                                            on_off_file[kwargs['light']].append(file)
                                    if len(on_off_file[kwargs['light']]) == 0:
                                        raise BaseException("No photo-IV data in (%s/%s, %s) folder" %
                                                            (self._month, self._day, kwargs['code']))
                                    else:
                                        if 'number' in kwargs:
                                            for file in on_off_file[kwargs['light']]:
                                                if number_format() in file:
                                                    self._file_location.append(set_location(self._CodeFolder, file))
                                                    self._filename_list.append(set_name(self._month, self._day, file))
                                                    if show == 'Yes':
                                                        print('Successfully load %s' % file)
                                        else:
                                            for file in on_off_file[kwargs['light']]:
                                                self._file_location.append(set_location(self._CodeFolder, file))
                                                self._filename_list.append(set_name(self._month, self._day, file))
                                else:
                                    raise BaseException('Wrong light status: "%s" (L/D)' % kwargs['light'])
                            else:
                                for file in Filename:
                                    self._file_location.append(set_location(self._CodeFolder, file))
                                    self._filename_list.append(set_name(self._month, self._day, file))
            else:
                raise BaseException("%s/%s files not found!" % (self._month, self._day))
        else:
            for date in self._DateDirectory:
                self._month = int(int(date) / 100)
                self._day = int(date) % 100
                self._DateFolder = location + '/' + date
                os.chdir(self._DateFolder)
                self._CodeDirectory = [name for name in os.listdir(".") if os.path.isdir(name) and name.isdigit()]
                if 'code' in kwargs:
                    if str(kwargs['code']) in self._CodeDirectory:
                        self._CodeFolder = self._DateFolder + '/' + str(kwargs['code'])
                        if show == 'Yes':
                            print('Successfully read %s in %s/%s folder' % (kwargs['code'], self._month, self._day))
                        Filename = [name for name in os.listdir(self._CodeFolder)
                                    if filename_check(name, str(kwargs['code']), self._DateFolder)]
                        if 'light' in kwargs:
                            on_off_file = {'L': [], 'D': []}
                            if kwargs['light'] in on_off_file:
                                for file in Filename:
                                    if kwargs['light'] in file:
                                        on_off_file[kwargs['light']].append(file)
                                for file in on_off_file[kwargs['light']]:
                                    self._file_location.append(set_location(self._CodeFolder, file))
                                    self._filename_list.append(set_name(self._month, self._day, file))
                                if show == 'Yes':
                                    print('Successfully load all files in %s at %s/%s' %
                                          (kwargs['code'], self._month, self._day))
                                    print(self._filename_list)
                            else:
                                raise BaseException('Wrong light status: "%s" (L/D)' % kwargs['light'])
                        else:
                            for file in Filename:
                                self._file_location.append(set_location(self._CodeFolder, file))
                                self._filename_list.append(set_name(self._month, self._day, file))
                else:
                    key2_list = check_key(self._key2)
                    if len(key2_list) == 0:
                        for code_subdir in self._CodeDirectory:
                            self._CodeFolder = self._DateFolder + '/' + code_subdir
                            Filename = [name for name in os.listdir(self._CodeFolder)
                                        if filename_check(name, str(code_subdir), self._DateFolder)]
                            if 'light' in kwargs:
                                on_off_file = {'L': [], 'D': []}
                                if kwargs['light'] in on_off_file:
                                    for file in Filename:
                                        if kwargs['light'] in file:
                                            on_off_file[kwargs['light']].append(file)
                                    else:
                                        if 'number' in kwargs:
                                            for file in on_off_file[kwargs['light']]:
                                                if number_format() in file:
                                                    self._file_location.append(set_location(self._CodeFolder, file))
                                                    self._filename_list.append(set_name(self._month, self._day, file))
                                                    if show == 'Yes':
                                                        print('Successfully load %s' % file)
                                        else:
                                            for file in on_off_file[kwargs['light']]:
                                                self._file_location.append(set_location(self._CodeFolder, file))
                                                self._filename_list.append(set_name(self._month, self._day, file))
                                else:
                                    raise BaseException('Wrong light status: "%s" (L/D)' % kwargs['light'])
                            else:
                                for file in Filename:
                                    self._file_location.append(set_location(self._CodeFolder, file))
                                    self._filename_list.append(set_name(self._month, self._day, file))
                    else:
                        code_info = {codeFolder: self.code_info(int(codeFolder), show='No')
                                     for codeFolder in self._CodeDirectory}
                        pass_codeFolder = {key: [] for key in key2_list}
                        for folder in code_info:
                            for key in key2_list:
                                if code_info[folder][key] == kwargs[key]:
                                    pass_codeFolder[key].append(folder)

                        # 檢查是否有完全沒找到的 key2
                        NotFound = dict()
                        for key in pass_codeFolder:
                            if len(pass_codeFolder[key]) == 0:
                                NotFound[key] = kwargs[key]
                        if len(NotFound.keys()) == 0:
                            # 找出交集
                            FoundFolder = set(pass_codeFolder[list(pass_codeFolder.keys())[0]])
                            if len(pass_codeFolder) > 1:
                                for element in list(pass_codeFolder.keys())[1:]:
                                    FoundFolder.intersection_update(set(pass_codeFolder[element]))
                            self._key2_foundfolder = list(FoundFolder)
                            self._key2_foundfolder.sort()  # 固定 code 資料夾順序，否則每次開啟都會不同。
                            if len(self._key2_foundfolder) > 0:
                                if show == 'Yes':
                                    print('Search result in %s/%s folder for %s : %s ' %
                                          (self._month, self._day,
                                           {key: kwargs[key] for key in key2_list}, self._key2_foundfolder))
                                for code in self._key2_foundfolder:
                                    self._CodeFolder = self._DateFolder + '/' + code
                                    if show == 'Yes':
                                        print('Successfully read %s in %s/%s folder' % (code, self._month, self._day))
                                    Filename = [name for name in os.listdir(self._CodeFolder)
                                                if filename_check(name, str(code), self._DateFolder)]
                                    if 'light' in kwargs:
                                        on_off_file = {'L': [], 'D': []}
                                        if kwargs['light'] in on_off_file:
                                            for file in Filename:
                                                if kwargs['light'] in file:
                                                    on_off_file[kwargs['light']].append(file)
                                            if len(on_off_file[kwargs['light']]) > 0:
                                                if 'number' in kwargs:
                                                    for file in on_off_file[kwargs['light']]:
                                                        if number_format() in file:
                                                            self._file_location.append(
                                                                set_location(self._CodeFolder, file))
                                                            self._filename_list.append(
                                                                set_name(self._month, self._day, file))
                                                            print('Successfully load %s' % file)
                                                else:
                                                    for file in on_off_file[kwargs['light']]:
                                                        self._file_location.append(set_location(self._CodeFolder, file))
                                                        self._filename_list.append(
                                                            set_name(self._month, self._day, file))
                                        else:
                                            raise BaseException('Wrong light status: "%s" (L/D)' % kwargs['light'])
                                    else:
                                        for file in Filename:
                                            self._file_location.append(set_location(self._CodeFolder, file))
                                            self._filename_list.append(set_name(self._month, self._day, file))

            # 顯示錯誤訊息
            if 'code' in kwargs:
                if 'light' in kwargs:
                    if kwargs['light'] == 'L':
                        if len(self._filename_list) == 0:
                            raise BaseException("No %s photo-IV data is recorded." % kwargs['code'])
                    elif kwargs['light'] == 'D':
                        if len(self._filename_list) == 0:
                            raise BaseException("No %s dark-IV data is recorded." % kwargs['code'])
                else:
                    if len(self._filename_list) == 0:
                        raise BaseException("%s files not found in every folder!" % kwargs['code'])
            else:
                if 'light' in kwargs:
                    if kwargs['light'] == 'L':
                        if len(self._filename_list) == 0:
                            raise BaseException("No photo-IV data is recorded.")
                    elif kwargs['light'] == 'D':
                        if len(self._filename_list) == 0:
                            raise BaseException("No dark-IV data is recorded.")
                else:
                    if len(self._filename_list) == 0:
                        raise BaseException("No IV data is recorded.")

    def _key_contradiction(self, kwargs):
        for set_key in kwargs:
            if set_key not in self._key1 and set_key not in self._key2:
                raise BaseException('"%s" is not the default key %s and %s' % (set_key, self._key1, self._key2))
            if set_key in self._key2 and ('code' in kwargs or 'number' in kwargs):
                raise BaseException('(code, number) and (Y, active, FGR, AGR, order) '
                                    'cannot be simultaneously indicated.')
        if 'date' not in kwargs and 'number' in kwargs and len(kwargs) > 1:
            raise BaseException("Indication of number while date is arbitrary is meaningless.")
        if 'number' in kwargs and len(kwargs) == 1:
            raise BaseException('Only "%s" information is not enough.' % list(kwargs.keys())[0])

    @staticmethod
    def _reverse_filename(filename):
        date = filename.split("-")[0].split("/")[0] + filename.split("-")[0].split("/")[1]
        code = filename.split("-")[1][:-1]
        light = filename.split("-")[1][-1]
        number = filename.split("-")[2]
        info = [filename.split("[")[1].split("]")[0].split(",")[0]] + \
               [value.split(" ")[1] for i, value in enumerate(filename.split("[")[1].split("]")[0].split(",")) if i > 0]

        return [date, code, light, number] + info

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
                info = {'active': 240, 'FGR': NoneValue, 'AGR': NoneValue, 'order': NoneValue, 'depth': int(x), 'Y': int(Y)}
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

    def print_key(self):
        """
        印出允許的 IVDB key 值。
        :return: ['date', 'code', 'light', 'number', 'active', 'FGR', 'AGR', 'order', 'depth', 'Y']

        >> IVDB.print_key()
        >> Acceptable keys: ['date', 'code', 'light', 'number', 'active', 'FGR', 'AGR', 'order', 'depth', 'Y']
        """
        key_list = self._key1 + self._key2
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
        self._file_list = []
        for file in self._file_location:
            self._file_list.append(Read(file, 'HP4156'))
        return self._file_list

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
        if show == 'Yes':
            print('[index]     mm/dd-code-number-[active, FGR, AGR, order, depth, Y]')
            for i, filename in enumerate(self._filename_list):
                number = '[' + str(i) + ']'
                print('{:>5} {:>55}'.format(number, filename))
        return self._filename_list

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
            for index, filename in enumerate(self._filename_list):
                if name == filename:
                    return [index]
            raise BaseException('"%s" not found!' % name)
        else:
            self._key_contradiction(kwargs)
            index_list = []
            for index, filename in enumerate(self._filename_list):
                info_name_list = self._key1 + self._key2
                info_list = self._reverse_filename(filename)
                name_info = {info_name_list[index]: element for index, element in enumerate(info_list)}

                flag = 1
                for key, value in kwargs.items():
                    if str(value) != name_info[key]:
                        flag = flag * 0

                if flag == 1:
                    index_list.append(index)
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
        database = self._location
        location = []

        if show == 'on':
            print('[index]')
        for i, filename in enumerate(self._filename_list):
            filename_split = filename.split("-")
            datefolder = filename_split[0].split("/")[0] + filename_split[0].split("/")[1]
            codefolder = filename_split[1][:-1]
            txt_filename = filename_split[1] + '-' + filename_split[2] + '.TXT'
            relative_path = datefolder + '/' + codefolder + '/' + txt_filename
            if relative == 'on':
                location.append(relative_path)
            elif relative == 'off':
                location.append(database + '/' + relative_path)
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
        for set_key in kwargs:
            if set_key not in ['Y', 'active', 'FGR', 'AGR', 'order']:
                raise BaseException('"%s" is not in the permissible data list %s' %
                                    (set_key, ['Y', 'active', 'FGR', 'AGR', 'order']))

        def day_DL(measurement, Code):
            if Code in measurement:
                return measurement[Code]
            else:
                return ''

        def range_list(value):
            return range(1, value + 1)

        def code_generator():
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
                          3: {x: {1: 3, 2: 5, 3: np.inf} for x in range_list(xmax[Xmin][3])},
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
            for index, filename in enumerate(self._filename_list):
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
                            CodeList[code] = {'X': X, 'Y': Y,
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
                                    if kwargs[key] != CodeList[code][key]:
                                        del CodeList[code]
                                        break

            statistics_path = self._location + '/statistics.csv'
            if os.path.isfile(statistics_path):
                with open(statistics_path, newline='') as csvfile:
                    reader = csv.DictReader(csvfile)
                    headers = reader.fieldnames
                    start = len(self._key1) + len(self._key2) + 4
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
        DataTable = {'X': [], 'Y': [], 'x': [], 'y': [], 'code': [], 'day-DL': [],
                     'Test': [], 'active': [], 'AGR': [], 'FGR': [], 'order': []}
        for key in self._EC:
            DataTable[key] = []
        for code, info in code_generator().items():
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

    def iv_plot(self, number, label_flag='on', label='partial', title=None,
                draw=None, remove=None, colorparm=None, **kwargs):

        self._key_contradiction(kwargs)
        Color_Storage = dict()

        def col_number(CurveList):
            if int(len(CurveList) / 5) == 0:
                return 1
            else:
                return int(len(CurveList) / 5)

        def check_key(info, kwargs):
            info_name_list = self._key1 + self._key2
            name_info = {info_name_list[index]: element for index, element in enumerate(info)}
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

        def color_func(info, colorparm):
            if colorparm is None:
                return ColorSet[i % len(ColorSet)]
            else:
                if colorparm not in self._key1 and colorparm not in self._key2:
                    raise BaseException("Wrong color parameter: %s" % colorparm)
                else:
                    info_name_list = self._key1 + self._key2
                    name_info = {info_name_list[index]: element for index, element in enumerate(info)}
                    if name_info[colorparm] not in Color_Storage:
                        Color_Storage[name_info[colorparm]] = [ColorSet[len(list(Color_Storage.items())) - 1]]
                    else:
                        Color_Storage[name_info[colorparm]].append(1)
                    return Color_Storage[name_info[colorparm]][0]

        def label_func(i):
            if colorparm is None:
                if label == 'partial':
                    date, code, light, number = curve_info[i][:4]
                    return '[%s] %s' % (i, date + '-' + code + light + '-' + number)
                elif label == 'full':
                    return curve_info[i]
                else:
                    raise BaseException("Wrong name settings(partial/full): %s" % label)
            else:
                name_info = {info_name_list[index]: element for index, element in enumerate(curve_info[i])}
                i_name_info = {info_name_list[index]: element for index, element in enumerate(curve_info[i])}
                if len(Color_Storage[i_name_info[colorparm]]) <= 1:
                    return '%s' % (colorparm + ': ' + name_info[colorparm])
                else:
                    return None

        def linestyle_func(i):
            if colorparm is None:
                return LineStyleSet[int(i / len(ColorSet))]
            else:
                return '-'

        plt.figure(number)
        ColorSet = ['dodgerblue', 'yellowgreen', 'goldenrod', 'darkviolet',
                    'darkorange', 'brown', 'b', 'hotpink', 'fuchsia', 'g', 'royalblue', 'tomato',
                    'purple', 'olive', 'darkgreen']
        LineStyleSet = ['-', '-.', '--', ':', (0, (3, 1, 1, 1)), (0, (3, 1, 1, 1, 1, 1)), (0, (3, 5, 1, 5))]

        info_name_list = self._key1 + self._key2
        curve_info = [self._reverse_filename(filename) for i, filename in enumerate(self.filename(show='No'))]

        if draw is not None and remove is not None:
            raise BaseException("Cannot use draw_list and except_list simultaneously")

        CurveList = []
        IV = self.read()
        for i in range(len(IV)):
            if draw is None:
                if remove is None:
                    if check_key(curve_info[i], kwargs):
                        CurveList.append(i)
                else:
                    if i not in remove:
                        if check_key(curve_info[i], kwargs):
                            CurveList.append(i)
            else:
                if i in draw:
                    if check_key(curve_info[i], kwargs):
                        CurveList.append(i)

        if len(CurveList) == 0:
            raise BaseException("No such curve!")
        if len(CurveList) >= len(ColorSet) * len(LineStyleSet):
            raise BaseException("Too many curves (%s) (Max:%s)" % (len(IV), len(ColorSet) * len(LineStyleSet)))
        for i in CurveList:
            plt.plot(IV[i].X, abs(IV[i].Y), color=color_func(curve_info[i], colorparm),
                     linestyle=linestyle_func(i),
                     label=label_func(i))
        if label_flag == 'on':
            plt.legend(loc='best', ncol=col_number(CurveList))
        elif label_flag != 'off':
            raise BaseException("Wrong label_flag (on/off): %s" % label_flag)
        if title is not None:
            plt.title(title, fontsize=15)
        plt.grid()
        plt.yscale('log')
        plt.xlabel('Voltage (V)')
        plt.ylabel('Current (A)')


class Statistics(IVDataBase):
    def __init__(self, location, show='no'):
        """
        1. 目的：這是用以統計製程與電性數據的工具，可以在此新增、移除、修改與繪製統計資料，例如擊穿電壓與元件座標、FGR spacing等。
        2. 輸入：因為統計表單是針對所有數據統計，而不是針對特定幾種元件的數據，所以參數僅有資料庫路徑。
        3. 輸出：執行後會在資料庫的最上層目錄中自動生成統計數據 statistics.csv，其數據列表可由 IVDataBase 中的 self._EC 修改。
        :param location: 資料庫路徑，例如 '/Users/Ethan/PycharmProjects/Python/Avalanche Photodiode/Experiment/IV'。
        :param show: 顯示目錄與檔案之成功讀取訊息，用以排除故障。
        """
        super().__init__(location, show)
        self._statistics = location + '/statistics.csv'
        self._characteristics = ['X', 'Y', 'x', 'y'] + self._EC
        self._data_list = self._key1 + self._key2 + self._characteristics

        # 檢查儲存目錄是否已有統計資料，若有則檢查元件清單是否一致，否則新建統計資料
        if not os.path.isfile(self._statistics):
            with open(self._statistics, 'w', newline='') as csvfile:
                writer = csv.DictWriter(csvfile, fieldnames=self._data_list)
                writer.writeheader()
                for filename in self._filename_list:
                    # 初始化數據
                    info_name_list = self._key1 + self._key2
                    info_list = self._reverse_filename(filename)
                    name_info = {info_name_list[index]: element for index, element in enumerate(info_list)}
                    name_info['X'], name_info['Y'], name_info['x'], name_info['y'] = name_info['code']
                    for key in self._characteristics:
                        if key not in name_info:
                            name_info[key] = ''
                    # 儲存數據
                    writer.writerow(name_info)
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
            for filename in self._filename_list:
                info_name_list = self._key1 + self._key2
                info_list = self._reverse_filename(filename)
                name_info = {info_name_list[index]: element for index, element in enumerate(info_list)}

                flag = 0  # 檢查元件資料是否不同於所有已存資料。
                with open(self._statistics, newline='') as csvfile:
                    reader = csv.DictReader(csvfile)
                    row_count = sum(1 for row in reader)

                with open(self._statistics, newline='') as csvfile:
                    reader = csv.DictReader(csvfile)
                    for row in reader:
                        for key in self._key1:
                            if row[key] != name_info[key]:
                                flag += 1
                                break

                if flag == row_count:
                    name_info['X'], name_info['Y'], name_info['x'], name_info['y'] = name_info['code']
                    for key in self._characteristics:
                        if key not in name_info:
                            name_info[key] = ''
                    with open(self._statistics, 'a', newline='') as csvfile:
                        writer = csv.DictWriter(csvfile, fieldnames=self._data_list)
                        writer.writerow(name_info)

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
        for key in kwargs:
            if key not in self._characteristics[4:] and key not in self._key1:
                raise BaseException("%s is not in the permissible list: %s" %
                                    (key, self._key1 + self._characteristics[4:]))
        tempkey = []
        for key in kwargs:
            if key in self._key1:
                tempkey.append(key)

        def gate(row, kwargs, tempkey):
            flag = 0
            for key in tempkey:
                if row[key] != str(kwargs[key]):
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
                        if row[key] != '' and row[key] != str(value):
                            ans = input('[%s-%s%s-%s] Change %s = %s to %s? (Y/N):  '
                                        % (row['date'], row['code'], row['light'], row['number'], key, row[key], value))
                            if ans == 'Y':
                                row[key] = value
                        else:
                            row[key] = value
                    writer.writerow(row)
                else:
                    writer.writerow(row)
            if flag == 0:
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
            if key not in self._key1:
                raise BaseException("%s is not in permissible list %s" % (key, self._key1))
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
                        if key in self._characteristics[4:]:
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
        if x not in self._data_list:
            raise BaseException('xlabel %s is not in the permissible list %s' % (x, self._data_list))
        elif y not in self._data_list:
            raise BaseException('ylabel %s is not in the permissible list %s' % (y, self._data_list))

        for key in kwargs:
            if key not in self._data_list and key != 'Ymin' and key != 'Ymax':
                temp = self._data_list + ['Ymin', 'Ymax']
                raise BaseException('key %s is not in the permissible list %s' % (key, temp))

        def gate(row, kwargs):
            flag = 0
            for key in kwargs:
                if key != 'Ymin' and key != 'Ymax':
                    if row[key] != str(kwargs[key]):
                        flag = 1
                elif key == 'Ymin' and int(row['Y']) < kwargs[key]:
                    flag = 1
                elif key == 'Ymax' and int(row['Y']) > kwargs[key]:
                    flag = 1
            if flag == 0:
                return True
            else:
                return False

        X = []
        Y = []
        str_flag = 0
        with open(self._statistics, 'r', newline='') as csvfile:
            reader = csv.DictReader(csvfile)
            for row in reader:
                if row['code'] == str(6613) and 'active' in kwargs:
                    print(kwargs['active'], row['active'], str(kwargs['active']) == row['active'], row['Vpt'], row['Vb'])
                if x == 'Vpt' and 'active' in kwargs:
                    if kwargs['active'] == 400:
                        print('[Out] Code: %s      Active: %s    Vpt: %s    Vb: %s' % (row['code'], row['active'], row['Vpt'], row['Vb']))
                        if row['active'] == str(400):
                            print('[------] Active: %s     Vpt: %s    Vb: %s' % (row['active'], row['Vpt'], row['Vb']))
                if ynull == 'No':
                    if row[y] != 'null':
                        if str_flag == 0 and not row[x].replace('.', '', 1).replace('-', '', 1).isdigit() and xnull == 'Yes':
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
            if x == 'Vpt' and 'active' in kwargs:
                if kwargs['active'] == 400:
                    print(X, Y)
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
            if key not in self._key1 and key not in self._key2:
                raise BaseException("%s is not in permissible list %s" % (key, self._key1 + self._key2))
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
                        if key in self._characteristics[4:] and row[key] == '':
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
                print('Uncompleted parameters:')
                for key in Remaining:
                    for item in Remaining[key]:
                        print('[%s]   %s' % (key, item))
