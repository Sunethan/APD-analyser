import numpy as np
import csv
import tkinter as tk
import tkinter.simpledialog as ts
from tkinter import ttk
from PIL import ImageTk, Image
import tkinter.filedialog as tf
import tkinter.messagebox as tM
import physfunc as phys
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.backends.backend_tkagg as tkagg
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import matplotlib.pylab as pylab
import datetime
import time
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (7.5, 7),
          'axes.labelsize': 'x-large',
          'axes.titlesize':'x-large',
          'xtick.labelsize':'x-large',
          'ytick.labelsize':'x-large'}
pylab.rcParams.update(params)
plt.rcParams.update({'font.size': 10})

defaultsettings, defaultList = phys.parse_csv('defaults.csv', 4, 0, 0, 'No')
mesh_list = phys.ConvertToInt(defaultsettings[defaultList[0]])
regions = defaultsettings[defaultList[1]]
mesh_line_set = defaultsettings[defaultList[2]]
default_meshlineposition = phys.ConvertToFloat(defaultsettings[defaultList[3]])
concentration_data = defaultsettings[defaultList[4]]
default_concs = phys.ConvertToFloat(defaultsettings[defaultList[5]])
default_interface = defaultsettings[defaultList[6]]
default_IntPosition = phys.ConvertToFloat(defaultsettings[defaultList[7]])
default_Vnumbers = phys.ConvertToInt(defaultsettings[defaultList[8]])
default_Vmeshlines = [round(float(i)) for i in defaultsettings[defaultList[9]]]
default_others = phys.ConvertToFloat(defaultsettings[defaultList[10]])
default_info = defaultsettings[defaultList[11]]
default_infoValues = defaultsettings[defaultList[12]]
default_Plot = defaultsettings[defaultList[13]]
default_PlotValue = defaultsettings[defaultList[14]]


# temp data
vlines = [0] * len(default_meshlineposition)
line_temp = [0] * 2
temp_Vpt = 0


# IO function
def SetDefault():
    with open('defaults.csv', 'w') as default_file:
        default_file.write('%s' % defaultList[0])
        for i in range(1, len(defaultList)):
            default_file.write(',%s' % defaultList[i])

        lengthlist = [len(mesh_list), len(regions), len(mesh_line_set),
                      len(default_meshlineposition), len(concentration_data),
                      len(default_concs), len(default_interface),
                      len(default_IntPosition), len(default_Vnumbers),
                      len(default_Vmeshlines), len(default_others),
                      len(default_info), len(default_infoValues),
                      len(default_Plot), len(default_PlotValue)]

        for i in range(max(lengthlist)):
            datalist = str(mesh_values[i].get())
            writelist = [0] * len(lengthlist)
            for j in range(len(lengthlist)):
                if i >= lengthlist[j]:
                    writelist[j] = ''
                else:
                    if j == 0:
                        writelist[j] = str(mesh_values[i].get())
                    elif j == 1:
                        writelist[j] = regions[i]
                    elif j == 2:
                        writelist[j] = mesh_line_set[i]
                    elif j == 3:
                        writelist[j] = str(mesh_line_positions[i].get())
                    elif j == 4:
                        writelist[j] = concentration_data[i]
                    elif j == 5:
                        writelist[j] = str(concentrations[i].get())
                    elif j == 6:
                        writelist[j] = default_interface[i]
                    elif j == 7:
                        writelist[j] = str(interface_list[i].get())
                    elif j == 8:
                        writelist[j] = str(Vmesh_number[i].get())
                    elif j == 9:
                        if i == 0:
                            writelist[j] = str(0)
                        else:
                            writelist[j] = str(Vmesh_lineEnt[i - 1].get())
                    elif j == 10:
                        writelist[j] = str(othersettings[i].get())
                    elif j == 11:
                        writelist[j] = default_info[i]
                    elif j == 12:
                        if len(InfoArray[i].get()) == 0:
                            writelist[j] = 'None'
                        else:
                            writelist[j] = str(InfoArray[i].get())
                    elif j == 13:
                        writelist[j] = default_Plot[i]
                    elif j == 14:
                        if len(PlotDataArray[i].get()) == 0:
                            writelist[j] = 'None'
                        else:
                            writelist[j] = str(PlotDataArray[i].get())
                if j > 0:
                    datalist = datalist + ',' + writelist[j]
            default_file.write('\n%s' % datalist)


def readcurveDoubleClick(event):
    global x_raw
    global doping_raw
    try:
        data
    except NameError:
        tM.showerror("Error", "You haven't load data.")
    else:
        try:
            load_index
        except NameError:
            tM.showerror("Error", "You need to select a data first.")
        else:
            x_raw = np.asarray(data[curvename[load_index]]['X'])
            doping_raw = np.asarray(data[curvename[load_index]]['Y'])
            StoreVariables('x_raw(length)', len(x_raw))
            StoreVariables('doping_raw(length)', len(doping_raw))
            log_region.insert(tk.END, "\nSuccefully read raw doping data (x, conc.)")
            log_region.see(tk.END)
            ax1.set_title(r'Doping profile @ %s' % curvename[load_index])
            bar1.draw()
            ax2.set_title(
                r'Depletion region boundary position relative to the P/N junction @ %s' % curvename[load_index])
            bar2.draw()


def readcurve():
    global x_raw
    global doping_raw
    try:
        data
    except NameError:
        tM.showerror("Error", "You haven't load data.")
    else:
        try:
            load_index
        except NameError:
            tM.showerror("Error", "You need to select a data first.")
        else:
            x_raw = np.asarray(data[curvename[load_index]]['X'])
            doping_raw = np.asarray(data[curvename[load_index]]['Y'])
            StoreVariables('x_raw(length)', len(x_raw))
            StoreVariables('doping_raw(length)', len(doping_raw))
            log_region.insert(tk.END, "\nSuccefully read raw doping data (x, conc.)")
            log_region.see(tk.END)
            ax1.set_title(r'Doping profile @ %s' % curvename[load_index])
            bar1.draw()
            ax2.set_title(
                r'Depletion region boundary position relative to the P/N junction @ %s' % curvename[load_index])
            bar2.draw()


def openfile():
    global datapath
    datapath = tf.askopenfilename()
    csv_path.delete('1.0', tk.END)
    csv_path.insert(tk.END, datapath)
    tv.yview_moveto(1)


def importsettingsDoubleClick(event):

    try:
        conc
    except NameError:
        tM.showerror("Error", "You should find layer structure first, and then import "
                              "data to replace those values.")
    else:
        global importdata, importlist, import_data, import_curvename, import_x_raw, import_doping_raw

        SavedLog, SavedPath = phys.parse_csv('log/AllSaved.csv', 4, 0, 0, 'No')
        importdata, importlist = phys.parse_csv(SavedLog[SavedPath[1]][data_index], 4, 0, 0, 'No')

        import_mesh_list = phys.ConvertToInt(importdata[importlist[0]])
        import_regions = importdata[importlist[1]]
        import_mesh_line_set = importdata[importlist[2]]
        import_meshlineposition = phys.ConvertToFloat(importdata[importlist[3]])
        import_concentration_data = importdata[importlist[4]]
        import_concs = phys.ConvertToFloat(importdata[importlist[5]])
        import_interface = importdata[importlist[6]]
        import_IntPosition = phys.ConvertToFloat(importdata[importlist[7]])
        import_Vnumbers = phys.ConvertToInt(importdata[importlist[8]])
        import_Vmeshlines = [round(float(i)) for i in importdata[importlist[9]]]
        import_others = phys.ConvertToFloat(importdata[importlist[10]])
        import_info = importdata[importlist[11]]
        import_infoValues = importdata[importlist[12]]

        for i in range(len(mesh_values)):
            mesh_values[i].delete(0, 'end')
            mesh_values[i].insert(0, import_mesh_list[i])
        for i in range(len(concentrations)):
            concentrations[i].delete(0, 'end')
            concentrations[i].insert(0, '%.2e' % import_concs[i])
        for i in range(len(mesh_line_positions)):
            mesh_line_positions[i].delete(0, 'end')
            mesh_line_positions[i].insert(0, import_meshlineposition[i])
        for i in range(len(interface_list)):
            interface_list[i].delete(0, 'end')
            interface_list[i].insert(0, import_IntPosition[i])
        for i in range(len(Vmesh_number)):
            Vmesh_number[i].delete(0, 'end')
            Vmesh_number[i].insert(0, import_Vnumbers[i])
        for i in range(len(Vmesh_lineEnt)):
            Vmesh_lineEnt[i].delete(0, 'end')
            Vmesh_lineEnt[i].insert(0, import_Vmeshlines[i + 1])
        for i in range(len(othersettings)):
            othersettings[i].delete(0, 'end')
            othersettings[i].insert(0, import_others[i])
        for i in range(len(InfoArray)):
            InfoArray[i].delete(0, 'end')
            InfoArray[i].insert(0, import_infoValues[i])

        StatisticsData, StatisticsDataList = phys.parse_csv('log/Statistics.csv', 4, 0, 0, 'No')
        XArrayEnt.insert(0, StatisticsData[StatisticsDataList[len(StatisticsDataList) - 2]][data_index])
        YArrayEnt.insert(0, StatisticsData[StatisticsDataList[-1]][data_index])


def importsettings():

    try:
        conc
    except NameError:
        tM.showerror("Error", "You should find layer structure first, and then import "
                              "data to replace those values.")
    else:
        global importdata, importlist, import_data, import_curvename, import_x_raw, import_doping_raw

        SavedLog, SavedPath = phys.parse_csv('log/AllSaved.csv', 4, 0, 0, 'No')
        importdata, importlist = phys.parse_csv(SavedLog[SavedPath[1]][data_index], 4, 0, 0, 'No')

        import_mesh_list = phys.ConvertToInt(importdata[importlist[0]])
        import_regions = importdata[importlist[1]]
        import_mesh_line_set = importdata[importlist[2]]
        import_meshlineposition = phys.ConvertToFloat(importdata[importlist[3]])
        import_concentration_data = importdata[importlist[4]]
        import_concs = phys.ConvertToFloat(importdata[importlist[5]])
        import_interface = importdata[importlist[6]]
        import_IntPosition = phys.ConvertToFloat(importdata[importlist[7]])
        import_Vnumbers = phys.ConvertToInt(importdata[importlist[8]])
        import_Vmeshlines = [round(float(i)) for i in importdata[importlist[9]]]
        import_others = phys.ConvertToFloat(importdata[importlist[10]])
        import_info = importdata[importlist[11]]
        import_infoValues = importdata[importlist[12]]

        for i in range(len(mesh_values)):
            mesh_values[i].delete(0, 'end')
            mesh_values[i].insert(0, import_mesh_list[i])
        for i in range(len(concentrations)):
            concentrations[i].delete(0, 'end')
            concentrations[i].insert(0, '%.2e' % import_concs[i])
        for i in range(len(mesh_line_positions)):
            mesh_line_positions[i].delete(0, 'end')
            mesh_line_positions[i].insert(0, import_meshlineposition[i])
        for i in range(len(interface_list)):
            interface_list[i].delete(0, 'end')
            interface_list[i].insert(0, import_IntPosition[i])
        for i in range(len(Vmesh_number)):
            Vmesh_number[i].delete(0, 'end')
            Vmesh_number[i].insert(0, import_Vnumbers[i])
        for i in range(len(Vmesh_lineEnt)):
            Vmesh_lineEnt[i].delete(0, 'end')
            Vmesh_lineEnt[i].insert(0, import_Vmeshlines[i + 1])
        for i in range(len(othersettings)):
            othersettings[i].delete(0, 'end')
            othersettings[i].insert(0, import_others[i])
        for i in range(len(InfoArray)):
            InfoArray[i].delete(0, 'end')
            InfoArray[i].insert(0, import_infoValues[i])

        StatisticsData, StatisticsDataList = phys.parse_csv('log/Statistics.csv', 4, 0, 0, 'No')
        XArrayEnt.insert(0, StatisticsData[StatisticsDataList[19]][data_index])
        YArrayEnt.insert(0, StatisticsData[StatisticsDataList[20]][data_index])

def StoreVariables(name, element):
    if abs(element) >= 1e3:
        VarListBox.insert(tk.END, "%s :     %.3e" % (name, element))
        VarListBox.see(tk.END)
    else:
        VarListBox.insert(tk.END, "%s :     %.5f" % (name, element))
        VarListBox.see(tk.END)


def LoadTreeviewData():
    if len(tv.get_children()) > 0:
        tv.delete(*tv.get_children())
    HistoryData, HisDataList = phys.parse_csv('log/Statistics.csv', 4, 0, 0, 'No')
    for i in range(len(HistoryData[HisDataList[0]])):
        datatree = [HistoryData[HisDataList[0]][i], HistoryData[HisDataList[1]][i], HistoryData[HisDataList[2]][i], HistoryData[HisDataList[3]][i],
                    HistoryData[HisDataList[4]][i], HistoryData[HisDataList[5]][i], HistoryData[HisDataList[6]][i],
                    HistoryData[HisDataList[7]][i], HistoryData[HisDataList[8]][i], HistoryData[HisDataList[9]][i],
                    HistoryData[HisDataList[10]][i], HistoryData[HisDataList[11]][i], HistoryData[HisDataList[12]][i], HistoryData[HisDataList[13]][i],
                    HistoryData[HisDataList[14]][i], HistoryData[HisDataList[17]][i], HistoryData[HisDataList[18]][i]]
        tv.insert('', 'end', text=(i + 1), values=datatree)
    tv.yview_moveto(1)

def Save():
    SetDefault()
    Settings = datetime.datetime.now().strftime("%Y-%m-%d-%H-%M-%S-Settings.csv")
    default = csv.reader(open('defaults.csv', 'r'))
    SettingsSavePath = 'log/' + Settings
    savefile = csv.writer(open(SettingsSavePath, 'w'))
    for row in default:
        savefile.writerow(row)

    NoEntryValue = 0
    NodeEntryValue = 0
    RingEntryValue = 0
    FormEntryValue = 0
    OtherInfoValue = 0
    formvalue = 0
    NdValue = 0
    NcValue = 0
    XArray = 0
    YArray = 0
    if len(NoEntry.get()) == 0:
        NoEntryValue = 'None'
    else:
        NoEntryValue = str(NoEntry.get())
    if len(NodeEntry.get()) == 0:
        NodeEntryValue = 'None'
    else:
        NodeEntryValue = str(NodeEntry.get())
    if len(RingEntry.get()) == 0:
        RingEntryValue = 'None'
    else:
        RingEntryValue = str(RingEntry.get())
    if len(FormEntry.get()) == 0:
        FormEntryValue = 'None'
    else:
        FormEntryValue = str(FormEntry.get())
    if len(NDEntry.get()) == 0:
        NdValue = 'None'
    else:
        NdValue = str(NDEntry.get())
    if len(NCEntry.get()) == 0:
        NcValue = 'None'
    else:
        NcValue = str(NCEntry.get())
    if len(VptEntry.get()) == 0:
        VptValue = 'None'
    else:
        VptValue = str(VptEntry.get())
    if len(VbEntry.get()) == 0:
        VbValue = 'None'
    else:
        VbValue = str(VbEntry.get())
    if len(OtherInfo.get()) == 0:
        OtherInfoValue = 'None'
    else:
        OtherInfoValue = str(OtherInfo.get())
    if f.get() == 1:
        formvalue = 'Planar'
    else:
        formvalue = 'Cylindrical'

    if len(XArrayEnt.get()) == 0:
        XArray = 'None'
    else:
        XArray = str(XArrayEnt.get())

    if len(YArrayEnt.get()) == 0:
        YArray = 'None'
    else:
        YArray = str(YArrayEnt.get())

    if temp_VptValue[0] == 0:
        Vpt = 'None'
    else:
        Vpt = "{:.2f}".format(temp_VptValue[0])

    try:
        datapath
    except NameError:
        tM.showerror("Error", "Sorry, file input is needed.")
    else:
        try:
            curvename
        except NameError:
            tM.showerror("Error", "Sorry, you should choose a curve from datalist first.")
        else:
            # save the saved data

            # 這是 import 所需的程式碼
            ithData = len(tv.get_children()) + 1
            with open('log/AllSaved.csv', 'a') as allsaved:
                allsaved.write('\n%s,%s,%s,%s' % (ithData, SettingsSavePath, datapath, curvename[load_index]))

            Vb_value = 0
            Mn_value = 0
            Mp_value = 0
            maxerror_value = 0

            if Vb_data[0] == 0:
                Vb_value = 'None'
                Mn_value = 'None'
                Mp_value = 'None'
                maxerror_value = 'None'
                with open('log/Statistics.csv', 'a') as result:
                    Date = datetime.datetime.now().strftime("%Y-%m-%d")
                    Time = datetime.datetime.now().strftime("%H:%M:%S")
                    result.write('\n%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s' %
                                 (Date, Time, NoEntryValue, NodeEntryValue, NdValue, NcValue,
                                  Vpt, VptValue, Vb_value, VbValue, maxerror_value, RingEntryValue,
                                  FormEntryValue, formvalue, OtherInfoValue, Mn_value, Mp_value,
                                  curvename[load_index], datapath, XArray, YArray))
            else:
                Vb_value = str(round(Vb, 2))
                Mn_value = str(round(Mn, 2))
                Mp_value = str(round(Mp, 2))
                maxerror_value = str(round(maxerror, 3))
                with open('log/Statistics.csv', 'a') as result:
                    Date = datetime.datetime.now().strftime("%Y-%m-%d")
                    Time = datetime.datetime.now().strftime("%H:%M:%S")
                    result.write('\n%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s' %
                                 (Date, Time, NoEntryValue, NodeEntryValue, NdValue, NcValue, Vpt, VptValue,
                                  Vb_value, VbValue, maxerror_value, RingEntryValue, FormEntryValue,
                                  formvalue, OtherInfoValue, Mn_value, Mp_value, curvename[load_index],
                                  datapath, XArray, YArray))
            '''
            try:
                Vb
            except NameError:
                Vb_value = 'None'
                Mn_value = 'None'
                Mp_value = 'None'
                maxerror_value = 'None'
                with open('log/Statistics.csv', 'a') as result:
                    Date = datetime.datetime.now().strftime("%Y-%m-%d")
                    Time = datetime.datetime.now().strftime("%H:%M:%S")
                    result.write('\n%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s' %
                                 (Date, Time, NoEntryValue, NodeEntryValue, NdValue, NcValue,
                                  Vpt, VptValue, Vb_value, VbValue, maxerror_value, RingEntryValue,
                                  FormEntryValue, formvalue, OtherInfoValue, Mn_value, Mp_value,
                                  curvename[load_index], datapath, XArray, YArray))
            else:
                Vb_value = str(round(Vb, 2))
                Mn_value = str(round(Mn, 2))
                Mp_value = str(round(Mp, 2))
                maxerror_value = str(round(maxerror, 3))
                with open('log/Statistics.csv', 'a') as result:
                    Date = datetime.datetime.now().strftime("%Y-%m-%d")
                    Time = datetime.datetime.now().strftime("%H:%M:%S")
                    result.write('\n%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s' %
                                 (Date, Time, NoEntryValue, NodeEntryValue, NdValue, NcValue, Vpt, VptValue,
                                  Vb_value, VbValue, maxerror_value, RingEntryValue, FormEntryValue,
                                  formvalue, OtherInfoValue, Mn_value, Mp_value, curvename[load_index],
                                  datapath, XArray, YArray))
            '''

            LoadTreeviewData()
            tv.yview_moveto(1)
            log_region.insert(tk.END, "\nSuccessfully save & reload!")
            log_region.see(tk.END)


def retrieve_input():
    global data
    global curvename
    data, curvename = phys.parse_csv(csv_path.get("1.0", 'end-1c'), 4, 0, 0, 'Yes')
    log_region.insert(tk.END, "\n %s" % curvename)
    log_region.see(tk.END)
    if len(LoadListBox.get(0, tk.END)) > 0:
        LoadListBox.delete(0, tk.END)
    for index, element in enumerate(curvename):
        LoadListBox.insert(index + 1, curvename[index])


# Tool functions
def plotfitting():
    wrongrelativevalue = 0
    for i in range(len(mesh_line_set) - 1):
        if float(mesh_line_positions[i].get()) > float(mesh_line_positions[i + 1].get()):
            wrongrelativevalue = 1
            print(i, float(mesh_line_positions[i].get()))
            print(i + 1, float(mesh_line_positions[i + 1].get()))
            break
    if wrongrelativevalue == 1:
        tM.showerror("Error", "Some meshlines have wrong relative value.")
    else:
        if line_temp[1] == 0:
            global line2
            readcurve()
            FitDoping()
            meshlineplot()
            line2, = ax1.plot(new_mesh, abs(fit_doping),
                              marker='o', markersize=4,
                              linestyle='-.', color='red',
                              label='fitted profile')
            line_temp[1] = line2
            ax1.legend(loc='upper right')
            bar1.draw()
            log_region.insert(tk.END, "\nSuccessfully plot fitting curve!")
            log_region.see(tk.END)
        else:
            line_temp[1].remove()
            readcurve()
            FitDoping()
            meshlineplot()
            line2, = ax1.plot(new_mesh, abs(fit_doping),
                              marker='o', markersize=4,
                              linestyle='-.', color='red',
                              label='fitted profile')
            line_temp[1] = line2
            ax1.legend(loc='upper right')
            bar1.draw()
            log_region.insert(tk.END, "\nSuccessfully plot fitting curve!")
            log_region.see(tk.END)


def ShowChoice():
    print('ShowChoice', load_index)


def FormChoice():
    print('FormChoice', f.get())


def onselect(event):
    w = event.widget
    global load_index
    load_index = w.curselection()[0]
    value = w.get(load_index)

'''
def DataSelect(event):
    w = event.widget
    global data_index
    data_index = w.curselection()[0]
    value = w.get(data_index)
'''

def selectTreeItem(event):
    global data_index
    curItem = tv.focus()
    data_index = int(round(tv.item(curItem)['text']) - 1)
    #print(tv.item(curItem))
    #print(tv.item(curItem)['text'])


def VarSelect(event):
    w = event.widget
    global var_index
    var_index = w.curselection()[0]
    value = w.get(var_index)


def LayerAnalysis():
    try:
        x_raw
    except NameError:
        tM.showerror("Error", "Data is not loaded.")
    else:
        global conc, xj, xj_mc, xj_cg, xj_ga, xj_ab, MLThick, xj_raw

        conc = np.zeros(len(concentrations))
        for i in range(len(conc)):
            conc[i] = concentrations[i].get()
        xj, xj_mc, xj_cg, xj_ga, xj_ab, MLThick = phys.structure_analysis(x_raw, doping_raw, conc)

        # find fit xj
        for i in range(len(doping_raw) - 1):
            if doping_raw[i] * doping_raw[i + 1] < 0:
                xj_raw = x_raw[i + 1]
                break

        for i in range(len(doping_raw) - 1):
            if x_raw[i] > xj_raw \
                    and doping_raw[i + 1] < doping_raw[i]:
                fit_Chconc = doping_raw[i]
                xj_cg = (x_raw[i] + x_raw[i + 1]) / 2
                break
        StoreVariables('Chconc(fit)', fit_Chconc)
        NCEntry.delete(0, 'end')
        NCEntry.insert(0, '%.2e' % fit_Chconc)
        concentrations[1].delete(0, 'end')
        concentrations[1].insert(0, '%.2e' % fit_Chconc)
        mesh_line_positions[4].delete(0, 'end')
        mesh_line_positions[4].insert(0, xj_cg)

        Nd_x = (xj_raw + float(mesh_line_positions[3].get())) / 2
        for i in range(len(doping_raw) - 1):
            if x_raw[i] > Nd_x:
                fit_ND = doping_raw[i]
                break
        StoreVariables('ND(fit)', fit_ND)
        NDEntry.delete(0, 'end')
        NDEntry.insert(0, '%.2e' % fit_ND)
        concentrations[0].delete(0, 'end')
        concentrations[0].insert(0, '%.2e' % fit_ND)

        StoreVariables('xj', xj)
        StoreVariables('xj_mc', xj_mc)
        StoreVariables('xj_cg', xj_cg)
        StoreVariables('xj_ga', xj_ga)
        StoreVariables('xj_ab', xj_ab)
        StoreVariables('MLThick', MLThick)
        # delete values
        if len(mesh_line_positions[1].get()) > 0:
            for i in range(2, 7):
                if (i == 4 and xj_cg == 0) or (i == 5 and xj_ga == 0) or (i == 6 and xj_ab == 0):
                    continue
                mesh_line_positions[i].delete(0, 'end')
        mesh_line_positions[2].insert(0, xj)
        mesh_line_positions[3].insert(0, xj_mc)
        if xj_cg != 0:
            mesh_line_positions[4].insert(0, xj_cg)
        if xj_ga != 0:
            mesh_line_positions[5].insert(0, xj_ga)
        if xj_ab != 0:
            mesh_line_positions[6].insert(0, xj_ab)
        log_region.insert(tk.END, "\nStructure analysis:"
                                  "\nxj: %s"
                                  "\nxj_mc: %s"
                                  "\nxj_cg: %s"
                                  "\nxj_ga: %s"
                                  "\nxj_ab: %s"
                                  "\nMLThick: %s" % (xj, xj_mc, xj_cg, xj_ga, xj_ab, MLThick))
        log_region.insert(tk.END, "\nInterface info:"
                                  "\nInP/InGaAsP : %s"
                                  "\nInGaAsP/InGaAs : %s"
                                  "\nInGaAs/InP : %s" % (xj_cg, xj_ga, xj_ab))


def FitDoping():
    global temp_import_a, temp_import_b, mesh_line, new_mesh, fit_doping, \
        xj_fit, InP_InGaAsP_interface, InGaAsP_InGaAs_interface, InGaAs_Buffer_interface, interface_lines
    try:
        x_raw
    except NameError:
        tM.showerror("Error", "you haven't read / create data")
    else:
        log_region.insert(tk.END, "\n------Start computing------")

        for i in range(len(mesh_list)):
            mesh_list[i] = int(mesh_values[i].get())
        log_region.insert(tk.END, "\nSuccessfully load mesh points!")
        mesh_line = []
        for i in range(len(mesh_line_set)):
            mesh_line.append(float(mesh_line_positions[i].get()))
            if i in range(2, 7):
                if i == 2:
                    StoreVariables('xj', float(mesh_line_positions[i].get()))
                elif i == 3:
                    StoreVariables('xj_mc', float(mesh_line_positions[i].get()))
                elif i == 4:
                    StoreVariables('xj_cg', float(mesh_line_positions[i].get()))
                elif i == 5:
                    StoreVariables('xj_ga', float(mesh_line_positions[i].get()))
                elif i == 6:
                    StoreVariables('xj_ab', float(mesh_line_positions[i].get()))
        log_region.insert(tk.END, "\nSuccessfully load mesh line positions!")
        new_mesh, fit_doping = phys.mesh(mesh_list, mesh_line, x_raw, doping_raw)


        line_density = []
        for i in range(len(mesh_list)):
            if i == 0:
                dx = float(mesh_line_positions[i].get()) - min(new_mesh)
            elif i == len(mesh_list) - 1:
                dx = max(new_mesh) - float(mesh_line_positions[-1].get())
            else:
                dx = float(mesh_line_positions[i].get()) - float(mesh_line_positions[i - 1].get())
            dn = mesh_list[i]
            line_density.append(dx / dn)
            StoreVariables('[%d] Density (um/grid)' % i, line_density[i])

        # find fit xj
        for i in range(len(fit_doping) - 1):
            if fit_doping[i] * fit_doping[i + 1] < 0:
                xj_fit = (new_mesh[i] + new_mesh[i + 1]) / 2
                break
        StoreVariables('xj(raw)', xj_raw)
        StoreVariables('xj(fit)', xj_fit)

        # set interface line
        InP_InGaAsP_interface = float(interface_list[0].get()) - xj_fit
        InGaAsP_InGaAs_interface = float(interface_list[1].get()) - xj_fit
        InGaAs_Buffer_interface = float(interface_list[2].get()) - xj_fit
        interface_lines = [InP_InGaAsP_interface, InGaAsP_InGaAs_interface, InGaAs_Buffer_interface]

        log_region.insert(tk.END, "\nSuccessfully fit doping profile!")
        log_region.insert(tk.END, "\n--Computing finished!--")
        log_region.see(tk.END)


def meshlineplot():
    wrongrelativevalue = 0
    for i in range(len(default_meshlineposition) - 1):
        if float(mesh_line_positions[i].get()) >= float(mesh_line_positions[i + 1].get()):
            wrongrelativevalue = 1
            break
    if wrongrelativevalue == 1:
        tM.showerror("Error", "Some meshline has wrong relative value.")
    else:
        log_region.insert(tk.END, "\nDrawing mesh lines...")
        log_region.see(tk.END)
        if vlines[0] == 0:
            for i in range(len(default_meshlineposition)):
                temp = ax1.vlines(x=float(mesh_line_positions[i].get()), ymin=5e13, ymax=1e19, linestyle='-.',
                                  color='k')
                vlines[i] = temp
        else:
            for i in range(len(default_meshlineposition)):
                vlines[i].remove()
            for i in range(len(default_meshlineposition)):
                temp = ax1.vlines(x=float(mesh_line_positions[i].get()), ymin=5e13, ymax=1e19, linestyle='-.',
                                  color='k')
                vlines[i] = temp
        ax1.set_ylim((5e13, 1e19))
        bar1.draw()
        log_region.insert(tk.END, "\nSuccessfully plot mesh lines!")
        log_region.see(tk.END)


def AuToFindVb(trial_number, step, error, x, e_goal, interface, s, sn, v_number, doping, form, rj, xj, v0_goal):
    gate = 0
    v = [v0_goal, 0]
    start_time = time.time()
    mn, mp = phys.multiplication_factor(x, e_goal, interface)
    error_n = abs(mn - 1)
    error_p = abs(mp - 1)
    print('Total mesh number: %s' % len(x))

    turn_marker = 0
    iteration = 1
    step_size1 = 2
    step_size2 = 2

    print('[1]', '  V: %.2f' % v0_goal, '  Mn: %.3f' % mn, '  Mp: %.3f' % mp,
          '  Error: %.3f' % max(error_n, error_p), '--- %.5s seconds ---' % (time.time() - start_time))
    with open('VbLog.csv', 'a') as VbLog:
        VbLog.write('\n[1] V: %.2f  Mn: %.3f  Mp: %.3f  (%.5ss)' % (v0_goal, mn, mp, (time.time() - start_time)))
    for i in range(1, trial_number):
        end = trial_number - 1
        if mn > 1 + error or mp > 1 + error:
            if gate == -1:
                if turn_marker == 0:
                    turn_marker = 1
                step = step / step_size1
                step_size2 = 2
                iteration = 1
            if gate == 1 and turn_marker == 1:
                iteration += 1
                if iteration == step_size1:
                    step_size1 += 1
                    step = step / step_size1
                    iteration = 1
            s = s + step
            gate = 1

            xp1, xn1, v, ef, xp_goal, xn_goal, v_goal, e_goal = \
                phys.depletion_boundary(v_number, s, s, x, doping, interface, 0, form, rj, xj)
            mn, mp = phys.multiplication_factor(x, ef[0, :], interface)
            error_n = abs(mn - 1)
            error_p = abs(mp - 1)
            print('[%d]' % (i + 1), '  V: %.2f' % v[0], '  Mn: %.3f' % mn, '  Mp: %.3f' % mp,
                  '  Error: %.3f' % max(error_n, error_p), '--- %.5s seconds ---  (+s: %s)' % ((time.time() - start_time), s))
            with open('VbLog.csv', 'a') as VbLog:
                VbLog.write('\n[%d] V: %.2f  Mn: %.3f  Mp: %.3f  (%.5ss)' % ((i + 1), v[0], mn, mp, (time.time() - start_time)))
            if i == end:
                error_n = abs(mn - 1)
                error_p = abs(mp - 1)
                print('Iterations larger than maximum %.3f' % v)
                print('[%d]' % (i + 1), '  V: %.2f' % v[0], '  Mn: %.3f' % mn, '  Mp: %.3f' % mp,
                      '  Error: %.3f' % max(error_n, error_p), '--- %.5s seconds ---' % (time.time() - start_time))
                with open('VbLog.csv', 'a') as VbLog:
                    VbLog.write('\n[%d] V: %.2f  Mn: %.3f  Mp: %.3f  (%.5ss)' % (
                    (i + 1), v[0], mn, mp, (time.time() - start_time)))
                return v, mn, mp, i, xp1[0], xn1[0]
        elif mn < 1 - error or mp < 1 - error:
            if gate == 1:
                if turn_marker == 0:
                    turn_marker = 1
                step = step / step_size2
                step_size1 = 2
                iteration = 1
            if gate == -1 and turn_marker == 1:
                iteration += 1
                if iteration == step_size2:
                    step_size2 += 1
                    step = step / step_size2
                    iteration = 1
            s = s - step
            gate = -1

            xp1, xn1, v, ef, xp_goal, xn_goal, v_goal, e_goal = \
                phys.depletion_boundary(v_number, s, s, x, doping, interface, 0, form, rj, xj)
            mn, mp = phys.multiplication_factor(x, ef[0, :], interface)
            error_n = abs(mn - 1)
            error_p = abs(mp - 1)
            print('[%d]' % (i + 1), '  V: %.2f' % v[0], '  Mn: %.3f' % mn, '  Mp: %.3f' % mp,
                  '  Error: %.3f' % max(error_n, error_p), '--- %.5s seconds ---  (-s: %s)' % ((time.time() - start_time), s))
            with open('VbLog.csv', 'a') as VbLog:
                VbLog.write('\n[%d] V: %.2f  Mn: %.3f  Mp: %.3f  (%.5ss)' % (
                (i + 1), v[0], mn, mp, (time.time() - start_time)))
            if i == end:
                error_n = abs(mn - 1)
                error_p = abs(mp - 1)
                print('Iterations larger than maximum %.3f' % v)
                print('[%d]' % (i + 1), '  V: %.2f' % v[0], '  Mn: %.3f' % mn, '  Mp: %.3f' % mp,
                      '  Error: %.3f' % max(error_n, error_p), '--- %.5s seconds ---' % (time.time() - start_time))
                with open('VbLog.csv', 'a') as VbLog:
                    VbLog.write('\n[%d] V: %.2f  Mn: %.3f  Mp: %.3f  (%.5ss)' % (
                    (i + 1), v[0], mn, mp, (time.time() - start_time)))
                return v, mn, mp, i, xp1[0], xn1[0]
        else:
            return v[0], mn, mp, i-1, s, sn


class Page(tk.Frame):
    def __init__(self, *args, **kwargs):
        tk.Frame.__init__(self, *args, **kwargs)
    def show(self):
        self.lift()


class Page1(Page):
    def __init__(self, *args, **kwargs):
        Page.__init__(self, *args, **kwargs)
        frame = tk.Frame(self)
        frame.pack(side="top", fill="both")
        Userguide = "1.開啟 *.CSV 資料檔\n" \
                    "2.Load file.\n" \
                    "3.\n"
        Read_button = "Read 按鈕：\n讀取 *.CSV 資料檔，儲存為 x_raw 與 doping_raw"
        Compute_button = "Compute 按鈕：\n1. 讀取濃度列表，計算結構參數（各層位置）\n" \
                         "2. 讀取分配網格區域的切割線（mesh lines），設定新的網格陣列。\n" \
                         "3. 儲存各材料之間的介面位置（interface）。\n" \
                         "4. 以新的網格，以對數方式擬合濃度分佈曲線。"
        tk.Message(frame, width=300, text=Userguide).pack(side="top")
        tk.Message(frame, width=300, text=Read_button).pack()
        tk.Message(frame, width=300, text=Compute_button).pack(side="bottom")


class Page2(Page):
    def __init__(self, *args, **kwargs):
        Page.__init__(self, *args, **kwargs)
        frame = tk.Frame(self)
        frame.pack()
        frame_0 = tk.Frame(frame)
        frame_1 = tk.Frame(frame_0)
        frame_4 = tk.Frame(frame_0)
        frame_2 = tk.Frame(frame)
        frame_3 = tk.Frame(frame)
        frame_0.pack(side="top", fill='x', expand=True, pady=5)
        frame_4.pack(side="left", fill='y', ipadx=8)  # file info frame
        frame_1.pack(side="left", fill='y')  # Settings frame
        frame_2.pack(fill='x')
        frame_3.pack(side="bottom", fill='x')

        global ax1, bar1, fig1
        fig1 = Figure(figsize=(18, 5), dpi=80, facecolor='w', edgecolor='k')
        ax1 = fig1.add_subplot(111)
        ax1.set_yscale('log')
        ax1.set_xlabel(r'Position ($\mu m$)')
        ax1.set_ylabel(r'Concentration ($cm^{-3}$)')
        ax1.set_title(r'Doping profile')
        ax1.set_ylim((5e13, 1e19))
        ax1.grid(True)
        bar1 = FigureCanvasTkAgg(fig1, frame_2)
        bar1.get_tk_widget().pack(side='right')
        tkagg.NavigationToolbar2Tk(bar1, frame_3)
        bar1.draw()


        def plot():
            global line1
            if line_temp[0] == 0:
                try:
                    conc
                except NameError:
                    tM.showerror("Error", "You should analyze the layer structure.")
                    #LayerAnalysis()
                    #readcurve()
                    #FitDoping()
                else:
                    line1, = ax1.plot(x_raw, abs(doping_raw), '-b', label='Raw data')
                    line_temp[0] = line1
                    ax1.legend(loc='upper right')
                    bar1.draw()
                    log_region.insert(tk.END, "\nSuccessfully plot!")
                    log_region.see(tk.END)
            else:
                line_temp[0].remove()
                readcurve()
                FitDoping()
                line1, = ax1.plot(x_raw, abs(doping_raw), '-b', label='Raw data')
                line_temp[0] = line1
                ax1.legend(loc='upper right')
                bar1.draw()
                log_region.insert(tk.END, "\nSuccessfully plot!")
                log_region.see(tk.END)


        def refreshFigure():
            if line_temp[0] == 0:
                tM.showerror("Error", "You haven't plotted any figure so there's no need to refresh.")
            else:
                readcurve()
                LayerAnalysis()
                meshlineplot()
                line1.set_xdata(x_raw)
                line1.set_ydata(abs(doping_raw))
                line1.set_label('Raw data')
                ax1.legend(loc='upper right')
                bar1.draw()
                log_region.insert(tk.END, "\nSuccessfully refresh selected curve and its legend!")
                log_region.see(tk.END)


        def deletelines():
            for i in range(len(default_meshlineposition)):
                vlines[i].remove()
                vlines[i] = 0
            bar1.draw()
            log_region.insert(tk.END, "\nSuccessfully remove vlines!")
            log_region.see(tk.END)


        def deletefittingcurve():
            if line_temp[1] == 0:
                tM.showerror("Error", "There's no fitted curve.")
            else:
                log_region.insert(tk.END, "\nDeleting fitting curve...")
                line_temp[1].remove()
                line_temp[1] = 0
                ax1.legend(loc='upper right')
                bar1.draw()
                log_region.insert(tk.END, "\nSuccessfully delete fitting curve!")
                log_region.see(tk.END)


        # data info
        tk.Label(frame_4, text="File Info.", font=('Times', 20, 'bold')).grid(row=0, column=0, columnspan=2)
        tk.Label(frame_4, text="No.", font=('Times', 16)).grid(row=1, column=0, sticky='w')
        tk.Label(frame_4, text="Node", font=('Times', 16)).grid(row=2, column=0, sticky='w')
        tk.Label(frame_4, text="Ring(Y/N)", font=('Times', 16)).grid(row=3, column=0, sticky='w')
        tk.Label(frame_4, text="Form(M/P)", font=('Times', 16)).grid(row=4, column=0, sticky='w')
        tk.Label(frame_4, text="Nd(cm-3)", font=('Times', 16)).grid(row=5, column=0, sticky='w')
        tk.Label(frame_4, text="Nc(cm-3)", font=('Times', 16)).grid(row=6, column=0, sticky='w')
        tk.Label(frame_4, text="Vpt", font=('Times', 16)).grid(row=7, column=0, sticky='w')
        tk.Label(frame_4, text="Vb", font=('Times', 16)).grid(row=8, column=0, sticky='w')
        tk.Label(frame_4, text="Other info.", font=('Times', 16)).grid(row=9, column=0, sticky='nw')
        global NoEntry, NodeEntry, RingEntry, FormEntry, NDEntry, NCEntry, VptEntry, VbEntry, OtherInfo, InfoArray
        NoEntry = tk.Entry(frame_4, width=12)
        NodeEntry = tk.Entry(frame_4, width=12)
        RingEntry = tk.Entry(frame_4, width=12)
        FormEntry = tk.Entry(frame_4, width=12)
        NDEntry = tk.Entry(frame_4, width=12)
        NCEntry = tk.Entry(frame_4, width=12)
        VptEntry = tk.Entry(frame_4, width=12)
        VbEntry = tk.Entry(frame_4, width=12)
        OtherInfo = tk.Entry(frame_4, width=12)
        NoEntry.grid(row=1, column=1, sticky='w')
        NodeEntry.grid(row=2, column=1, sticky='w')
        RingEntry.grid(row=3, column=1, sticky='w')
        FormEntry.grid(row=4, column=1, sticky='w')
        NDEntry.grid(row=5, column=1, sticky='w')
        NCEntry.grid(row=6, column=1, sticky='w')
        VptEntry.grid(row=7, column=1, sticky='w')
        VbEntry.grid(row=8, column=1, sticky='w')
        OtherInfo.grid(row=9, column=1)
        NoEntry.insert(0, default_infoValues[0])
        NodeEntry.insert(0, default_infoValues[1])
        RingEntry.insert(0, default_infoValues[2])
        FormEntry.insert(0, default_infoValues[3])
        NDEntry.insert(0, default_infoValues[4])
        NCEntry.insert(0, default_infoValues[5])
        VptEntry.insert(0, default_infoValues[6])
        VbEntry.insert(0, default_infoValues[7])
        OtherInfo.insert(0, default_infoValues[8])
        InfoArray = [NoEntry, NodeEntry, RingEntry, FormEntry, NDEntry, NCEntry, VptEntry, VbEntry, OtherInfo]


        # load/read data
        tk.Button(frame_1, text="Read", command=readcurve, width=10).grid(row=0, column=0, sticky='w')
        tk.Button(frame_1, text="F-Layer", command=LayerAnalysis, width=10).grid(row=1, column=0, sticky='e')
        tk.Button(frame_1, text="Plot", command=plot, width=10).grid(row=2, column=0)
        tk.Button(frame_1, text="Refresh", command=refreshFigure, width=10).grid(row=3, column=0)
        tk.Button(frame_1, text="Fitting!", command=plotfitting, width=10).grid(row=4, column=0)
        tk.Button(frame_1, text="Delete fitcurve", command=deletefittingcurve, width=10).grid(row=5, column=0)
        # mesh points
        tk.Label(frame_1, text="Region", font=('times', 20, 'bold')).grid(row=0, column=1)
        tk.Label(frame_1, text="Mesh number", font=('times', 20, 'bold')).grid(row=0, column=2)
        global mesh_values
        mesh_values = []
        for i in range(len(mesh_list)):
            j = i + 1
            tk.Label(frame_1, text=regions[i]).grid(row=j, column=1)
            mesh_points = tk.Entry(frame_1, width=4)
            mesh_points.grid(row=j, column=2)
            mesh_values.append(mesh_points)
            mesh_points.insert(0, mesh_list[i])

        # concentration settings
        tk.Label(frame_1, text="Region", font=('times', 20, 'bold')).grid(row=0, column=3)
        tk.Label(frame_1, text="Concentration", font=('times', 20, 'bold')).grid(row=0, column=4)
        global concentrations
        concentrations = []
        for i in range(len(concentration_data)):
            j = i + 1
            tk.Label(frame_1, text=concentration_data[i]).grid(row=j, column=3)
            concs = tk.Entry(frame_1, width=8)
            concs.grid(row=j, column=4)
            concentrations.append(concs)
            concs.insert(0, '%.2e' % default_concs[i])

        # mesh line positions
        tk.Label(frame_1, text="Mesh line", font=('times', 20, 'bold')).grid(row=0, column=5)
        tk.Label(frame_1, text="Position", font=('times', 20, 'bold')).grid(row=0, column=6)
        tk.Button(frame_1, text="Draw / Refresh lines", command=meshlineplot).grid(row=9, column=5)
        tk.Button(frame_1, text="Delete", command=deletelines).grid(row=9, column=6)
        global mesh_line_positions
        mesh_line_positions = []
        for i in range(len(mesh_line_set)):
            j = i + 1
            tk.Label(frame_1, text=mesh_line_set[i], width=15).grid(row=j, column=5)
            mesh_lines = tk.Entry(frame_1, width=7)
            mesh_lines.grid(row=j, column=6)
            mesh_lines.insert(0, default_meshlineposition[i])
            mesh_line_positions.append(mesh_lines)

class Page3(Page):
    def __init__(self, *args, **kwargs):
        Page.__init__(self, *args, **kwargs)

        frame = tk.Frame(self)
        frame.pack()
        frame_0 = tk.Frame(frame, pady=5)
        frame_1 = tk.Frame(frame_0)
        frame_2 = tk.Frame(frame)
        frame_3 = tk.Frame(frame)
        frame_4 = tk.Frame(frame_0)
        frame_41 = tk.Frame(frame_4)
        frame_42 = tk.Frame(frame_4)
        frame_43 = tk.Frame(frame_4)
        frame_44 = tk.Frame(frame_4)
        frame_0.pack(side="top", fill='x', expand=True)
        frame_1.pack(side="left")
        frame_4.pack(side="left")
        frame_2.pack(fill='x')
        frame_3.pack(side="bottom", fill='x')
        frame_41.pack(side="top")
        frame_42.pack(side="top", pady=5)
        frame_43.pack(side="top")
        frame_44.pack(side="top", pady=5)

        global ax2, bar2
        fig2 = Figure(figsize=(18, 5.8), dpi=80, facecolor='w', edgecolor='k')
        ax2 = fig2.add_subplot(111)
        ax2.set_xlabel(r'Voltage ($V$)')
        ax2.set_ylabel(r'Position ($\mu m$)')
        ax2.set_title(r'Depletion region boundary position relative to the P/N junction')
        ax2.grid(True)
        bar2 = FigureCanvasTkAgg(fig2, frame_2)
        bar2.get_tk_widget().pack(side='right')
        tkagg.NavigationToolbar2Tk(bar2, frame_3)

        # temp data
        global temp_VptValue, Vb_Marker, Vb_data
        temp_curvename = [0]
        temp_hline = [0] * len(default_interface)
        temp_P_hline = [0] * 2
        temp_N_hline = [0] * 2
        line_temp = [0] * 2
        vmeshlines = [0] * (len(default_Vnumbers) - 1)
        temp_VptVb_Lines = [0] * 3
        temp_VptValue = [0] * 2
        Vb_Marker = [0]
        Vb_data = [0] * 4

        def depletionwidth():
            start_time = time.time()

            # delete V mesh lines
            if vmeshlines[0] != 0:
                for i in range(len(default_Vnumbers) - 1):
                    vmeshlines[i].remove()
                    vmeshlines[i] = 0
                bar2.draw()
                log_region.insert(tk.END, "\nSuccessfully remove voltage mesh lines!")
                log_region.see(tk.END)

            # delete interface lines
            if temp_hline[0] != 0:
                for i in range(len(default_interface)):
                    temp_hline[i].remove()
                    temp_hline[i] = 0
                bar2.draw()
                log_region.insert(tk.END, "\nSuccessfully remove %s interface lines!" % len(default_interface))
                log_region.see(tk.END)

            # delete Vpt & Vb lines
            if temp_VptVb_Lines[0] != 0:
                temp_VptVb_Lines[0].remove()
                temp_VptVb_Lines[0] = 0
                temp_VptValue[0] = 0
                log_region.insert(tk.END, "\nSuccessfully remove Vpt lines!")
                log_region.see(tk.END)

            if temp_VptVb_Lines[1] != 0:
                temp_VptVb_Lines[1].remove()
                temp_VptVb_Lines[1] = 0
                temp_VptValue[1] = 0
                log_region.insert(tk.END, "\nSuccessfully remove Vpt2 lines!")
                log_region.see(tk.END)

            if temp_VptVb_Lines[2] != 0:
                temp_VptVb_Lines[2].remove()
                temp_VptVb_Lines[2] = 0
                temp_VptValue[2] = 0
                log_region.insert(tk.END, "\nSuccessfully remove Vb lines!")
                log_region.see(tk.END)

            # delete P-Line 1 & 2
            if temp_P_hline[0] != 0:
                for i in range(2):
                    temp_P_hline[i].remove()
                    temp_P_hline[i] = 0
                bar2.draw()
                log_region.insert(tk.END, "\nSuccessfully remove 2 P lines!")
                log_region.see(tk.END)

            # delete P-Line 1 & 2
            if temp_N_hline[0] != 0:
                for i in range(2):
                    temp_N_hline[i].remove()
                    temp_N_hline[i] = 0
                bar2.draw()
                log_region.insert(tk.END, "\nSuccessfully remove 2 N lines!")
                log_region.see(tk.END)

            # reset plot
            bar2.draw()


            try:
                xj_fit
            except NameError:
                tM.showerror("Error", "You haven't fitted the data.")
            else:
                for i in range(1, 5):
                    V1_meshlineEnt[i].config(text=(Vmesh_lineEnt[i - 1].get()))
                if float(XmaxEnt.get()) > xj_fit:
                    tM.showerror("Warning", "Xmax should be smaller than xj = %.5f." % xj_fit)
                else:
                    global interface, xp, xn, V, E, xp_goal, xn_goal, V_goal, E_goal, temp_Vgoal
                    Vnumber = round(float(VnumberEnt.get()))
                    xpmin = float(XminEnt.get())
                    xpmax = float(XmaxEnt.get())
                    Vgoal = float(VgoalEnt.get())
                    interface = [float(interface_list[0].get()),
                                 float(interface_list[1].get()),
                                 float(interface_list[2].get())]
                    if f.get() == 1:
                        form = 'planar'
                    else:
                        form = 'cylindrical'

                    xp, xn, V, E, xp_goal, xn_goal, V_goal, E_goal = \
                        phys.depletion_boundary(Vnumber, xpmin, xpmax, new_mesh, fit_doping,
                                                interface, Vgoal, form, 'none', 'none')
                    temp_Vgoal = V_goal
                    V_meshline = []
                    V_meshlist = []
                    WhileCount = 1
                    Vmaxtemp1 = 0
                    Vmaxtemp2 = 0
                    Vstep = abs(xpmin - xj_fit) / 5
                    Vgate = 0
                    temp_xpmin = 0
                    GenerateMesh = 0

                    if float(Vmesh_lineEnt[-1].get()) > max(V):
                        Vmaxtemp1 = max(V)
                    else:
                        Vmaxtemp2 = max(V)
                    print('[%s] Vmax: %.5f  at  xp = %.5f  or  xpj = %.5f --- %.5s seconds ---' % (
                        WhileCount, max(V), xpmin, xpmin - xj_fit, (time.time() - start_time)))

                    with open('VptLog.csv', 'w') as VptLog:
                        VptLog.write('[%s] Vmax: %.2f when xp = %.5f (%.5ss)' % (
                            WhileCount, max(V), xpmin, (time.time() - start_time)))

                    for i in range(1, len(default_Vmeshlines)):
                        V_meshline.append(float(Vmesh_lineEnt[i - 1].get()))
                        V_meshlist.append(int(Vmesh_number[i - 1].get()))

                    turn_marker = 0
                    iteration = 1
                    step_size1 = 2
                    step_size2 = 2
                    while abs(float(Vmesh_lineEnt[-1].get()) - max(V)) > float(VmaxErrorEntry.get()):
                        if float(Vmesh_lineEnt[-1].get()) > max(V):
                            temp_xpmin = xpmin

                            if Vgate == -1:
                                if turn_marker == 0:
                                    turn_marker = 1
                                Vstep = Vstep / step_size1
                                step_size2 = 2
                                iteration = 1
                            if Vgate == 1 and turn_marker == 1:
                                iteration += 1
                                if iteration == step_size1:
                                    step_size1 += 1
                                    Vstep = Vstep / step_size1
                                    iteration = 1
                            xpmin = xpmin - Vstep
                            Vgate = 1

                            #if Vgate == -1:
                            #    Vstep = Vstep / 2

                            xp, xn, V, E, xp_goal, xn_goal, V_goal, E_goal = \
                                phys.depletion_boundary(Vnumber, xpmin, xpmax, new_mesh, fit_doping,
                                                        interface, Vgoal, form, 'none', 'none')
                            WhileCount += 1
                            print('[%s] Vmax: %.5f  at  xp = %.5f  or  xpj = %.5f --- %.5s seconds ---' % (
                                WhileCount, max(V), xpmin, xpmin - xj_fit, (time.time() - start_time)))
                            with open('VptLog.csv', 'a') as VptLog:
                                VptLog.write('\n[%s] Vmax: %.2f when xp = %.5f (%.5ss)' % (
                                    WhileCount, max(V), xpmin, (time.time() - start_time)))
                            if (float(Vmesh_lineEnt[-1].get()) > max(V) and Vmaxtemp1 == max(V)) \
                                    or (float(Vmesh_lineEnt[-1].get()) < max(V) and Vmaxtemp2 == max(V)):
                                if GenerateMesh == 1:
                                    if abs(float(Vmesh_lineEnt[-1].get()) - Vmaxtemp1) < abs(
                                            float(Vmesh_lineEnt[-1].get()) - Vmaxtemp2):
                                        break
                                    else:
                                        xp, xn, V, E, xp_goal, xn_goal, V_goal, E_goal = \
                                            phys.depletion_boundary(Vnumber, xpmin, xpmax, new_mesh, fit_doping,
                                                                    interface, Vgoal, form, 'none', 'none')
                                        break
                                print('Generating new mesh...')
                                with open('VptLog.csv', 'a') as VptLog:
                                    VptLog.write('\nGenerating new mesh...')
                                mesh_line[0] = temp_xpmin
                                mesh_line[1] = (mesh_line[0] + mesh_line[2]) / 2
                                mesh_line_positions[0].delete(0, 'end')
                                mesh_line_positions[0].insert(0, temp_xpmin)
                                mesh_line_positions[1].delete(0, 'end')
                                mesh_line_positions[1].insert(0, mesh_line[1])
                                plotfitting()
                                print('New mesh generation finished!')
                                with open('VptLog.csv', 'a') as VptLog:
                                    VptLog.write('\nNew mesh generation finished!')
                                xp, xn, V, E, xp_goal, xn_goal, V_goal, E_goal = \
                                    phys.depletion_boundary(Vnumber, xpmin, xpmax, new_mesh, fit_doping,
                                                            interface, Vgoal, form, 'none', 'none')
                                WhileCount += 1
                                print('[%s] Vmax: %.5f  at  xp = %.5f  or  xpj = %.5f --- %.5s seconds ---' % (
                                    WhileCount, max(V), xpmin, xpmin - xj_fit, (time.time() - start_time)))
                                with open('VptLog.csv', 'a') as VptLog:
                                    VptLog.write('\n[%s] Vmax: %.2f when xp = %.5f (%.5ss)' % (
                                        WhileCount, max(V), xpmin, (time.time() - start_time)))
                                GenerateMesh = 1
                            if float(Vmesh_lineEnt[-1].get()) > max(V):
                                Vmaxtemp1 = max(V)
                            else:
                                Vmaxtemp2 = max(V)
                        else:

                            if Vgate == 1:
                                if turn_marker == 0:
                                    turn_marker = 1
                                Vstep = Vstep / step_size2
                                step_size1 = 2
                                iteration = 1
                            if Vgate == -1 and turn_marker == 1:
                                iteration += 1
                                if iteration == step_size2:
                                    step_size2 += 1
                                    Vstep = Vstep / step_size2
                                    iteration = 1
                            xpmin = xpmin + Vstep
                            Vgate = -1

                            #if Vgate == 1:
                            #    Vstep = Vstep / 2
                            #xpmin = xpmin + Vstep
                            #Vgate = -1

                            xp, xn, V, E, xp_goal, xn_goal, V_goal, E_goal = \
                                phys.depletion_boundary(Vnumber, xpmin, xpmax, new_mesh, fit_doping,
                                                        interface, Vgoal, form, 'none', 'none')
                            WhileCount += 1
                            print('[%s] Vmax: %.5f  at  xp = %.5f  or  xpj = %.5f --- %.5s seconds ---' % (
                            WhileCount, max(V), xpmin, xpmin - xj_fit, (time.time() - start_time)))
                            with open('VptLog.csv', 'a') as VptLog:
                                VptLog.write('\n[%s] Vmax: %.2f when xp = %.5f (%.5ss)' % (
                                    WhileCount, max(V), xpmin, (time.time() - start_time)))
                            if (float(Vmesh_lineEnt[-1].get()) > max(V) and Vmaxtemp1 == max(V)) \
                                    or (float(Vmesh_lineEnt[-1].get()) < max(V) and Vmaxtemp2 == max(V)):
                                if GenerateMesh == 1:
                                    if abs(float(Vmesh_lineEnt[-1].get()) - Vmaxtemp1) < abs(
                                            float(Vmesh_lineEnt[-1].get()) - Vmaxtemp2):
                                        xp, xn, V, E, xp_goal, xn_goal, V_goal, E_goal = \
                                            phys.depletion_boundary(Vnumber, xpmin, xpmax, new_mesh, fit_doping,
                                                                    interface, Vgoal, form, 'none', 'none')
                                        break
                                    else:
                                        break
                                print('Generating new mesh...')
                                with open('VptLog.csv', 'a') as VptLog:
                                    VptLog.write('\nGenerating new mesh...')
                                mesh_line[0] = temp_xpmin
                                mesh_line[1] = (mesh_line[0] + mesh_line[2]) / 2
                                mesh_line_positions[0].delete(0, 'end')
                                mesh_line_positions[0].insert(0, temp_xpmin)
                                mesh_line_positions[1].delete(0, 'end')
                                mesh_line_positions[1].insert(0, mesh_line[1])
                                plotfitting()
                                print('New mesh generation finished!')
                                with open('VptLog.csv', 'a') as VptLog:
                                    VptLog.write('\nNew mesh generation finished!')
                                xp, xn, V, E, xp_goal, xn_goal, V_goal, E_goal = \
                                    phys.depletion_boundary(Vnumber, xpmin, xpmax, new_mesh, fit_doping,
                                                            interface, Vgoal, form, 'none', 'none')
                                WhileCount += 1
                                print('[%s] Vmax: %.5f  at  xp = %.5f  or  xpj = %.5f --- %.5s seconds ---' % (
                                    WhileCount, max(V), xpmin, xpmin - xj_fit, (time.time() - start_time)))
                                with open('VptLog.csv', 'a') as VptLog:
                                    VptLog.write('\n[%s] Vmax: %.2f when xp = %.5f (%.5ss)' % (
                                        WhileCount, max(V), xpmin, (time.time() - start_time)))
                                GenerateMesh = 1
                            if float(Vmesh_lineEnt[-1].get()) > max(V):
                                Vmaxtemp1 = max(V)
                            else:
                                Vmaxtemp2 = max(V)

                    if GenerateMesh == 1:
                        mesh_line[0] = temp_xpmin
                        mesh_line[1] = (mesh_line[0] + mesh_line[2]) / 2
                        mesh_line_positions[0].delete(0, 'end')
                        mesh_line_positions[0].insert(0, temp_xpmin)
                        mesh_line_positions[1].delete(0, 'end')
                        mesh_line_positions[1].insert(0, mesh_line[1])
                        plotfitting()
                    XminEnt.delete(0, 'end')
                    XminEnt.insert(0, xpmin)

                    global TotalVmesh, V_mesh, xp_mesh, xn_mesh, E_mesh, Wn, Wp

                    TotalVmesh, V_mesh, xp_mesh, xn_mesh, E_mesh = \
                        phys.V_selection(V_meshline, V_meshlist, xp, xn, V, E)
                    Wn = xn_mesh - xj_fit
                    Wp = xj_fit - xp_mesh

                    # punch-through voltage
                    global Vpt, Vpt2
                    for i in range(len(Wn)):
                        if Wn[i] > InP_InGaAsP_interface:
                            if temp_VptValue[0] == 0:
                                Vpt = V_mesh[i] - (V_mesh[i] - V_mesh[i - 1]) * (Wn[i] - InP_InGaAsP_interface) / (
                                            Wn[i] - Wn[i - 1])
                                StoreVariables('Vpt', Vpt)
                                temp_VptValue[0] = Vpt
                            else:
                                if Wn[i] > InGaAsP_InGaAs_interface:
                                    Vpt2 = V_mesh[i] - (V_mesh[i] - V_mesh[i - 1]) * (Wn[i] - InGaAsP_InGaAs_interface) / (
                                            Wn[i] - Wn[i - 1])
                                    StoreVariables('Vpt2', Vpt2)
                                    temp_VptValue[1] = Vpt2
                                    break

                    StoreVariables('Vmax', max(V))
                    StoreVariables('V(Goal)', V_goal)

                    # reset Vb
                    if Vb_Marker[0] != 0:
                        Vb_data[0] = 0

                    global WnV_line, WpV_line
                    log_region.insert(tk.END, "\nPlotting depletion width vs voltage ...")
                    log_region.see(tk.END)

                    VptLogTemp = csv.reader(open('VptLog.csv', 'r'))
                    for row in VptLogTemp:
                        Vpt_Log.insert(tk.END, '%s\n' % str(row[0]))
                    Vpt_Log.insert(tk.END, '---------------------------------------------------')
                    Vpt_Log.see(tk.END)

                    if line_temp[0] == 0:
                        WnV_line, = ax2.plot(V_mesh, Wn, '-bo', markersize=3, label=r'$x_n$')
                        WpV_line, = ax2.plot(V_mesh, -Wp, '-ro', markersize=3, label=r'$x_p$')
                        line_temp[0] = WnV_line
                        line_temp[1] = WpV_line
                        ax2.legend(loc='upper left')
                        bar2.draw()
                        log_region.insert(tk.END, "\nSuccessfully plot!")
                        log_region.see(tk.END)
                    else:
                        log_region.insert(tk.END, "\nRemoving plot...")
                        log_region.see(tk.END)
                        line_temp[0].remove()
                        line_temp[1].remove()
                        log_region.insert(tk.END, "\nSuccessfully remove plot!")
                        log_region.insert(tk.END, "\nCreate new plot...!")
                        log_region.see(tk.END)
                        WnV_line, = ax2.plot(V_mesh, Wn, '-bo', markersize=3,
                                             label=r'$x_n$')
                        WpV_line, = ax2.plot(V_mesh, -Wp, '-ro', markersize=3,
                                             label=r'$x_p$')
                        line_temp[0] = WnV_line
                        line_temp[1] = WpV_line
                        ax2.legend(loc='upper left')
                        bar2.draw()
                        log_region.insert(tk.END, "\nSuccessfully plot!")
                        log_region.see(tk.END)

                    LogStatistics, LogStatisticsList = phys.parse_csv('log/Statistics.csv', 4, 0, 0, 'No')
                    if temp_curvename[0] != curvename[load_index] or len(LogStatistics[LogStatisticsList[0]]) == 0:
                        Save()
                    else:
                        SetDefault()
                        with open('log/TempStatistics.csv', 'w', newline='') as file:
                            tempfile = csv.writer(file, delimiter=',')
                            tempfile.writerow(LogStatisticsList)
                            if len(LogStatistics[LogStatisticsList[0]]) == 1:
                                # print('len = 1')
                                LogData = []
                                for j in range(len(LogStatisticsList)):
                                    LogData.append(str(LogStatistics[LogStatisticsList[j]][0]))
                                LogData[1] = str(datetime.datetime.now().strftime("%H:%M:%S"))
                                LogData[6] = str("{:.2f}".format(Vpt))
                                tempfile.writerow(LogData)
                                # print(LogData)
                            else:
                                # print('len > 1')
                                for i in range(len(LogStatistics[LogStatisticsList[0]])):
                                    if i < len(LogStatistics[LogStatisticsList[0]]) - 1:
                                        LogData = []
                                        for j in range(len(LogStatisticsList)):
                                            LogData.append(str(LogStatistics[LogStatisticsList[j]][i]))
                                        tempfile.writerow(LogData)
                                        #print(i, LogData)
                                    else:
                                        LogData = []
                                        for j in range(len(LogStatisticsList)):
                                            LogData.append(str(LogStatistics[LogStatisticsList[j]][i]))
                                        LogData[1] = str(datetime.datetime.now().strftime("%H:%M:%S"))
                                        LogData[6] = str("{:.2f}".format(Vpt))
                                        tempfile.writerow(LogData)

                        tempLog, tempLogList = phys.parse_csv('log/TempStatistics.csv', 4, 0, 0, 'No')
                        with open('log/Statistics.csv', 'w', newline='') as file2:
                            ReplaceStatistics = csv.writer(file2, delimiter=',')
                            ReplaceStatistics.writerow(LogStatisticsList)
                            # Date,Time,No.,Node,Nd,Nc,Vpt,Vpt(TCAD),Vb,Vb(TCAD),
                            # maxerror,Ring,Form(N/P),Coordinate(C/P),OtherInfo,Mn,
                            # Mp,curvename,datapath,XArray,YArray
                            for i in range(len(tempLog[tempLogList[0]])):
                                LogData = []
                                for j in range(len(tempLogList)):
                                    LogData.append(str(tempLog[tempLogList[j]][i]))
                                ReplaceStatistics.writerow(LogData)
                                # print(LogData)

                        LoadTreeviewData()
                    temp_curvename[0] = curvename[load_index]
                    print('-------------------------------------------------------------------------------')


        def vbsolver():
            if float(V2.get()) < float(VgoalEnt.get()):
                tM.showerror("Warning", "To check the mesh condition, V5(max) should be larger than goal voltage.")
            try:
                E_goal
            except NameError:
                tM.showerror("Error", "You need to find depletion width distribution first because it would "
                                      "generate a the initial electric fields, which correspondes to the"
                                      "voltage desired, and use it to find the multiplication factor to"
                                      "be the initial starting iteration point.")
            else:
                start_time_remesh = time.time()
                print('SPeedUpRemesh: ', SpeedUpRemesh.get())

                global Vb, Mn, Mp, count, maxerror, xp_b, xn_b

                if SpeedUpRemesh.get() != 1:
                    print('--- Use breakdown voltage solver ---')
                    log_region.insert(tk.END, "\n--- Use breakdown voltage solver --")
                    log_region.see(tk.END)

                    if f.get() == 1:
                        form = 'planar'
                    else:
                        form = 'cylindrical'

                    trial_number = round(float(TrialsEnt.get()))
                    step = float((StepEnt.get()))
                    error = float(TolEnt.get())
                    V_number = round(float((VnumberEnt.get())))
                    start_time_vbsolver = time.time()
                    Vb, Mn, Mp, count, xp_b, xn_b = \
                        AuToFindVb(trial_number, step, error, new_mesh,
                                       E_goal, interface, xp_goal, xn_goal, V_number, fit_doping, form, 'none', xj, V_goal)
                else:
                    print('Generating new mesh...')
                    print('Initial mesh number: %s' % len(new_mesh))
                    log_region.see(tk.END)

                    if f.get() == 1:
                        form = 'planar'
                    else:
                        form = 'cylindrical'

                    Modified_Vgoal = V_goal + 10
                    xpmin = float(mesh_line_positions[0].get())
                    xpmax = float(mesh_line_positions[1].get())
                    Vnumber = round(float(VnumberEnt.get()))
                    xp, xn, V, E, xpg_modified, xng_modified, Vg_modified, Eg_modified = \
                        phys.depletion_boundary(Vnumber, xpmin, xpmax, new_mesh, fit_doping,
                                                interface, Modified_Vgoal, form, 'none', 'none')

                    temp_mesh_list = mesh_list
                    temp_mesh_line = mesh_line
                    line_n_marker = 0
                    for i in range(len(mesh_line_positions)):
                        if float(mesh_line_positions[i].get()) < xpg_modified:
                            temp_mesh_list[i] = 5
                            line_p_marker = i
                        elif float(mesh_line_positions[i].get()) >= xng_modified:
                            temp_mesh_list[i + 1] = 5
                            if line_n_marker == 0:
                                line_n_marker = i
                    p_mesh_number = temp_mesh_list[line_p_marker + 1] * \
                                    abs(xpg_modified - float(mesh_line_positions[line_p_marker + 1].get())) / \
                                    abs(float(mesh_line_positions[line_p_marker].get()) - float(
                                        mesh_line_positions[line_p_marker + 1].get()))
                    temp_mesh_list[line_p_marker + 1] = int(round(p_mesh_number))
                    n_mesh_number = temp_mesh_list[line_n_marker] * \
                                    abs(xng_modified - float(mesh_line_positions[line_n_marker - 1].get())) / \
                                    abs(float(mesh_line_positions[line_n_marker].get()) - float(
                                        mesh_line_positions[line_n_marker - 1].get()))
                    if n_mesh_number != 0:
                        temp_mesh_list[line_n_marker] = int(round(n_mesh_number))
                    else:
                        temp_mesh_list[-1] = 5
                    temp_mesh_line[line_p_marker] = xpg_modified
                    temp_mesh_line[line_n_marker] = xng_modified
                    modified_mesh, modified_doping = phys.mesh(temp_mesh_list, temp_mesh_line, x_raw, doping_raw)

                    xp_2, xn_2, V_2, E_2, xp_goal_2, xn_goal_2, V_goal_2, E_goal_2 = \
                        phys.depletion_boundary(Vnumber, xpmin, xpmax, modified_mesh, modified_doping,
                                                interface, float(VgoalEnt.get()), form, 'none', 'none')
                    remesh_dt = time.time() - start_time_remesh

                    print('New mesh generation finished! Use %.5s seconds' % remesh_dt)
                    print('Final mesh number: %s' % len(modified_mesh))
                    log_region.see(tk.END)
                    print('--- Use breakdown voltage solver ---')
                    log_region.insert(tk.END, "\n--- Use breakdown voltage solver ---")
                    log_region.see(tk.END)
                    trial_number = round(float(TrialsEnt.get()))
                    step = float((StepEnt.get()))
                    error = float(TolEnt.get())
                    V_number = round(float((VnumberEnt.get())))

                    if f.get() == 1:
                        form = 'planar'
                    else:
                        form = 'cylindrical'

                    start_time_vbsolver = time.time()
                    '''
                    Vb, Mn, Mp, count, xp_b, xn_b = \
                        AuToFindVb(trial_number, step, error, new_mesh,
                                       E_goal, interface, xp_goal, xn_goal, V_number, fit_doping, form, 'none', xj, V_goal)
                    '''
                    Vb, Mn, Mp, count, xp_b, xn_b = \
                        AuToFindVb(trial_number, step, error, modified_mesh,
                                   E_goal_2, interface, xp_goal_2, xn_goal_2, V_number, modified_doping, form, 'none',
                                   xj, V_goal_2)

                error_n = abs(Mn - 1)
                error_p = abs(Mp - 1)
                maxerror = max(error_p, error_n)

                Vb_data[0] = Vb
                Vb_data[1] = Mn
                Vb_data[2] = Mp
                Vb_data[3] = maxerror

                dt = time.time() - start_time_vbsolver
                print('Count: %d  Vb: %.2f Done!  Mn: %.3f  Mp: %.3f  Error: %.3f  Use %.5s seconds' % (count, Vb, Mn, Mp, maxerror, dt))
                log_region.insert(tk.END, "\nCount:[%d]  Vb: %.2f Done!"
                                          "\nMn: %.3f   Mp: %.3f"
                                          "\nError: %.3f"
                                          "\nUse %.5s seconds" % (count + 1, Vb, Mn, Mp, maxerror, dt))
                log_region.insert(tk.END, "\nSuccessfully find breakdown voltage!")
                log_region.see(tk.END)
                Vb_Marker[0] = 1
                #Save()
                tv.yview_moveto(1)
                StoreVariables('Vb', Vb)
                StoreVariables('Mn', Mn)
                StoreVariables('Mp', Mp)
                StoreVariables('count', int(count))
                StoreVariables('maxerror', maxerror)
                StoreVariables('Xpb', xp_b)
                StoreVariables('Xnb', xn_b)

                VbLogTemp = csv.reader(open('VbLog.csv', 'r'))
                for row in VbLogTemp:
                    Vb_Log.insert(tk.END, '%s\n' % str(row[0]))
                Vb_Log.insert(tk.END, '---------------------------------------------------')
                Vb_Log.see(tk.END)

                log_region.insert(tk.END, "\nSuccessfully save & reload!")
                log_region.see(tk.END)

                LogStatistics, LogStatisticsList = phys.parse_csv('log/Statistics.csv', 4, 0, 0, 'No')

                with open('log/TempStatistics.csv', 'w', newline='') as file:
                    tempfile = csv.writer(file, delimiter=',')
                    tempfile.writerow(LogStatisticsList)
                    if len(LogStatistics[LogStatisticsList[0]]) == 1:
                        #print('len = 1')
                        LogData = []
                        for j in range(len(LogStatisticsList)):
                            LogData.append(str(LogStatistics[LogStatisticsList[j]][0]))
                        LogData[1] = str(datetime.datetime.now().strftime("%H:%M:%S"))
                        LogData[8] = str("{:.2f}".format(Vb))
                        LogData[10] = str("{:.4f}".format(maxerror))
                        tempfile.writerow(LogData)
                    else:
                        for i in range(len(LogStatistics[LogStatisticsList[0]])):
                            if i < len(LogStatistics[LogStatisticsList[0]]) - 1:
                                LogData = []
                                for j in range(len(LogStatisticsList)):
                                    LogData.append(str(LogStatistics[LogStatisticsList[j]][i]))
                                tempfile.writerow(LogData)
                            else:
                                LogData = []
                                for j in range(len(LogStatisticsList)):
                                    LogData.append(str(LogStatistics[LogStatisticsList[j]][i]))
                                LogData[1] = str(datetime.datetime.now().strftime("%H:%M:%S"))
                                LogData[8] = str("{:.2f}".format(Vb))
                                LogData[10] = str("{:.4f}".format(maxerror))
                                tempfile.writerow(LogData)

                tempLog, tempLogList = phys.parse_csv('log/TempStatistics.csv', 4, 0, 0, 'No')
                with open('log/Statistics.csv', 'w', newline='') as file2:
                    ReplaceStatistics = csv.writer(file2, delimiter=',')
                    ReplaceStatistics.writerow(LogStatisticsList)
                    # Date,Time,No.,Node,Nd,Nc,Vpt,Vpt(TCAD),Vb,Vb(TCAD),
                    # maxerror,Ring,Form(N/P),Coordinate(C/P),OtherInfo,Mn,
                    # Mp,curvename,datapath,XArray,YArray
                    for i in range(len(tempLog[tempLogList[0]])):
                        LogData = []
                        for j in range(len(tempLogList)):
                            LogData.append(str(tempLog[tempLogList[j]][i]))
                        ReplaceStatistics.writerow(LogData)

                LoadTreeviewData()
                print('-------------------------------------------------------------------------------')


        def ResizePlot():
            try:
                V_mesh
            except NameError:
                tM.showerror("Error", "You should make a plot first.")
            else:
                XaxisMin = min(V_mesh) - (max(V_mesh) - min(V_mesh)) / 6
                XaxisMax = max(V_mesh) + (max(V_mesh) - min(V_mesh)) / 6
                YaxisMin = -max(Wp) - (interface_lines[2] + max(Wp)) / 6
                YaxisMax = interface_lines[2] + (interface_lines[2] + min(Wp)) / 6
                ax2.set_xlim((XaxisMin, XaxisMax))
                ax2.set_ylim((YaxisMin, YaxisMax))
                bar2.draw()


        def Vmeshlineplot():
            for i in range(1, 5):
                V1_meshlineEnt[i].config(text=(Vmesh_lineEnt[i - 1].get()))
            try:
                Wn
            except NameError:
                tM.showerror("Error", "You should plot D-width first.")
            else:
                wrongrelativevalue = 0
                for i in range(len(Vmesh_lineEnt) - 1):
                    if float(Vmesh_lineEnt[i].get()) >= float(Vmesh_lineEnt[i + 1].get()):
                        wrongrelativevalue = 1
                        break
                if wrongrelativevalue == 1:
                    tM.showerror("Error", "Some voltage meshlines have wrong relative value.")
                else:
                    log_region.insert(tk.END, "\nDraw voltage mesh lines...")
                    log_region.see(tk.END)
                    if vmeshlines[0] == 0:
                        for i in range(len(Vmesh_lineEnt) - 1):
                            temp = ax2.vlines(x=float(Vmesh_lineEnt[i].get()), ymin=-max(Wp), ymax=max(Wn), linestyle=':',
                                              color='k')
                            vmeshlines[i] = temp
                    else:
                        for i in range(len(Vmesh_lineEnt) - 1):
                            vmeshlines[i].remove()
                        for i in range(len(Vmesh_lineEnt) - 1):
                            temp = ax2.vlines(x=float(Vmesh_lineEnt[i].get()), ymin=-max(Wp), ymax=max(Wn), linestyle=':',
                                              color='k')
                            vmeshlines[i] = temp
                    XaxisMin = min(V_mesh) - (max(V_mesh) - min(V_mesh)) / 6
                    XaxisMax = max(V_mesh) + (max(V_mesh) - min(V_mesh)) / 6
                    YaxisMin = -max(Wp) - (max(Wn) + max(Wp)) / 6
                    YaxisMax = max(Wn) + (max(Wn) + min(Wp)) / 6
                    ax2.set_xlim((XaxisMin, XaxisMax))
                    ax2.set_ylim((YaxisMin, YaxisMax))
                    bar2.draw()
                    log_region.insert(tk.END, "\nSuccessfully plot voltage mesh lines!")
                    log_region.see(tk.END)


        def deleteVmeshlines():
            if vmeshlines[0] == 0:
                tM.showerror("Error", "You haven't plotted voltage mesh lines.")
            else:
                for i in range(len(default_Vnumbers) - 1):
                    vmeshlines[i].remove()
                    vmeshlines[i] = 0
                bar2.draw()
                log_region.insert(tk.END, "\nSuccessfully remove voltage mesh lines!")
                log_region.see(tk.END)


        def DrawInterfaceLine():
            try:
                xj_fit
            except NameError:
                tM.showerror("Error", "You haven't fitted the data.")
            else:
                try:
                    V
                except NameError:
                    tM.showerror("Error", "You haven't plotted the depletion width distribution.")
                else:
                    log_region.insert(tk.END, "\nDrawing interface lines...")
                    log_region.see(tk.END)

                    if temp_hline[0] == 0:
                        for i in range(len(default_interface)):
                            temp_hline[i] = ax2.hlines(y=interface_lines[i], xmin=0, xmax=max(V_mesh),
                                                       linestyle=':', color='black')
                    else:
                        for i in range(len(default_interface)):
                            temp_hline[i].remove()
                        for i in range(len(default_interface)):
                            temp_hline[i] = ax2.hlines(y=interface_lines[i], xmin=0, xmax=max(V_mesh),
                                                       linestyle=':', color='black')

                    XaxisMin = min(V_mesh) - (max(V_mesh) - min(V_mesh)) / 6
                    XaxisMax = max(V_mesh) + (max(V_mesh) - min(V_mesh)) / 6
                    YaxisMin = -max(Wp) - (interface_lines[2] + max(Wp)) / 6
                    YaxisMax = interface_lines[2] + (interface_lines[2] + min(Wp)) / 6
                    ax2.set_xlim((XaxisMin, XaxisMax))
                    ax2.set_ylim((YaxisMin, YaxisMax))
                    bar2.draw()
                    log_region.insert(tk.END, "\nSuccessfully plot %s interface lines!" % len(default_interface))
                    log_region.see(tk.END)


        def DeletePunchLine():
            if temp_hline[0] == 0:
                tM.showerror("Error", "You haven't plotted interface lines.")
            else:
                for i in range(len(default_interface)):
                    temp_hline[i].remove()
                    temp_hline[i] = 0
                bar2.draw()
                log_region.insert(tk.END, "\nSuccessfully remove %s interface lines!" % len(default_interface))
                log_region.see(tk.END)


        def DrawVptVb():
            try:
                V_mesh
            except NameError:
                tM.showerror("Error", "You need to plot D-width first.")
            else:
                Vb_state = 0
                if Vb_data[0] == 0:
                    Vb_state = 1

                log_region.insert(tk.END, "\nDrawing Vpt line...")
                log_region.see(tk.END)
                if temp_VptVb_Lines[0] == 0:
                    try:
                        Vpt
                    except NameError:
                        tM.showerror("Warning", "you haven't found the punch-through voltage, please try larger maximum voltage.")
                    else:
                        temp1 = ax2.vlines(x=float(Vpt), ymin=-max(Wp), ymax=interface_lines[2], linestyle='-.',
                                           color='orange', label=r'$V_{pt}$')
                        temp2 = ax2.vlines(x=float(Vpt2), ymin=-max(Wp), ymax=interface_lines[2], linestyle='-.',
                                           color='lawngreen', label=r'$V_{pt2}$')
                        temp_VptVb_Lines[0] = temp1
                        temp_VptVb_Lines[1] = temp2
                        log_region.insert(tk.END, "\nVpt line finished")
                        log_region.see(tk.END)

                if Vb_state == 0:
                    log_region.insert(tk.END, "\nDrawing Vb line...")
                    log_region.see(tk.END)
                    if temp_VptVb_Lines[2] == 0:
                        temp = ax2.vlines(x=float(Vb), ymin=-max(Wp), ymax=interface_lines[2], linestyle='-.',
                                          color='m', label=r'$V_{b}$')
                        temp_VptVb_Lines[2] = temp
                    log_region.insert(tk.END, "\nVb line finished")
                    log_region.see(tk.END)

                ax2.legend(loc='upper left')
                bar2.draw()


        def ClearVptVb():
            if temp_VptVb_Lines[0] == 0:
                tM.showerror("Error", "You haven't plotted Vptmesh lines.")
            else:
                temp_VptVb_Lines[0].remove()
                temp_VptVb_Lines[0] = 0
                log_region.insert(tk.END, "\nSuccessfully remove Vpt lines!")
                log_region.see(tk.END)

            if temp_VptVb_Lines[1] == 0:
                tM.showerror("Error", "You haven't plotted Vb mesh lines.")
            else:
                temp_VptVb_Lines[1].remove()
                temp_VptVb_Lines[1] = 0
                log_region.insert(tk.END, "\nSuccessfully remove Vpt2 lines!")
                log_region.see(tk.END)

            if temp_VptVb_Lines[2] == 0:
                tM.showerror("Error", "You haven't plotted Vb mesh lines.")
            else:
                temp_VptVb_Lines[2].remove()
                temp_VptVb_Lines[2] = 0
                log_region.insert(tk.END, "\nSuccessfully remove Vb lines!")
                log_region.see(tk.END)
            ax2.legend(loc='upper left')
            bar2.draw()


        def DrawLine12():
            try:
                xj_fit
            except NameError:
                tM.showerror("Error", "You haven't fitted the data.")
            else:
                try:
                    V
                except NameError:
                    tM.showerror("Error", "You haven't plotted the depletion width distribution.")
                else:
                    log_region.insert(tk.END, "\nDrawing P-Line1 and P-Line2...")
                    log_region.see(tk.END)

                    x_Pline1 = float(mesh_line_positions[0].get()) - xj_fit
                    x_Pline2 = float(mesh_line_positions[1].get()) - xj_fit
                    x_PlineSet = [x_Pline1, x_Pline2]

                    if temp_P_hline[0] == 0:
                        for i in range(2):
                            temp_P_hline[i] = ax2.hlines(y=x_PlineSet[i], xmin=0, xmax=max(V_mesh),
                                                       linestyle='-.', color='green')
                    else:
                        for i in range(2):
                            temp_P_hline[i].remove()
                        for i in range(2):
                            temp_P_hline[i] = ax2.hlines(y=x_PlineSet[i], xmin=0, xmax=max(V_mesh),
                                                       linestyle='-.', color='green')

                    XaxisMin = min(V_mesh) - (max(V_mesh) - min(V_mesh)) / 6
                    XaxisMax = max(V_mesh) + (max(V_mesh) - min(V_mesh)) / 6
                    YaxisMin = -max(Wp) - (interface_lines[2] + max(Wp)) / 6
                    YaxisMax = interface_lines[2] + (interface_lines[2] + min(Wp)) / 6
                    ax2.set_xlim((XaxisMin, XaxisMax))
                    ax2.set_ylim((YaxisMin, YaxisMax))
                    bar2.draw()
                    log_region.insert(tk.END, "\nSuccessfully plot P-Line1 and P-Line2!")
                    log_region.see(tk.END)


        def Draw_NLine12():
            try:
                xj_fit
            except NameError:
                tM.showerror("Error", "You haven't fitted the data.")
            else:
                try:
                    V
                except NameError:
                    tM.showerror("Error", "You haven't plotted the depletion width distribution.")
                else:
                    log_region.insert(tk.END, "\nDrawing P-Line1 and P-Line2...")
                    log_region.see(tk.END)

                    x_Nline1 = float(mesh_line_positions[6].get()) - xj_fit
                    x_Nline2 = float(mesh_line_positions[-1].get()) - xj_fit
                    x_NlineSet = [x_Nline1, x_Nline2]

                    if temp_N_hline[0] == 0:
                        for i in range(2):
                            temp_N_hline[i] = ax2.hlines(y=x_NlineSet[i], xmin=0, xmax=max(V_mesh),
                                                       linestyle='-.', color='green')
                    else:
                        for i in range(2):
                            temp_N_hline[i].remove()
                        for i in range(2):
                            temp_N_hline[i] = ax2.hlines(y=x_NlineSet[i], xmin=0, xmax=max(V_mesh),
                                                       linestyle='-.', color='green')

                    XaxisMin = min(V_mesh) - (max(V_mesh) - min(V_mesh)) / 6
                    XaxisMax = max(V_mesh) + (max(V_mesh) - min(V_mesh)) / 6
                    YaxisMin = -max(Wp) - (interface_lines[2] + max(Wp)) / 6
                    YaxisMax = interface_lines[2] + (interface_lines[2] + min(Wp)) / 6
                    ax2.set_xlim((XaxisMin, XaxisMax))
                    ax2.set_ylim((YaxisMin, YaxisMax))
                    bar2.draw()
                    log_region.insert(tk.END, "\nSuccessfully plot P-Line1 and P-Line2!")
                    log_region.see(tk.END)


        def DeleteDrawLine12():
            if temp_P_hline[0] == 0:
                tM.showerror("Error", "You haven't plotted voltage mesh lines.")
            else:
                for i in range(2):
                    temp_P_hline[i].remove()
                    temp_P_hline[i] = 0
                bar2.draw()
                log_region.insert(tk.END, "\nSuccessfully remove %s interface lines!" % len(default_interface))
                log_region.see(tk.END)

        def DeleteDraw_NLine12():
            if temp_N_hline[0] == 0:
                tM.showerror("Error", "You haven't plotted voltage mesh lines.")
            else:
                for i in range(2):
                    temp_N_hline[i].remove()
                    temp_N_hline[i] = 0
                bar2.draw()
                log_region.insert(tk.END, "\nSuccessfully remove %s interface lines!" % len(default_interface))
                log_region.see(tk.END)


        tk.Button(frame_1, text="Plot D-width", width=10, command=depletionwidth).grid(row=0, column=0)
        tk.Button(frame_1, text="Find Vb", width=10, command=vbsolver).grid(row=1, column=0)
        tk.Button(frame_1, text="Resize", width=10, command=ResizePlot).grid(row=2, column=0)
        global SpeedUpRemesh, f
        f = tk.IntVar()
        f.set(1)
        SpeedUpRemesh = tk.IntVar()
        SpeedUpRemesh.set(1)
        tk.Button(frame_1, text="Draw Vpt & Vb", width=10, command=DrawVptVb).grid(row=3, column=0)
        tk.Button(frame_1, text="Clear Vpt & Vb", width=10, command=ClearVptVb).grid(row=4, column=0)
        tk.Checkbutton(frame_1, text="Vb-Remsh", variable=SpeedUpRemesh).grid(row=5, column=0, sticky='w')
        tk.Radiobutton(frame_1, indicatoron=0, text="Planar", justify=tk.LEFT, command=FormChoice, variable=f, value=1).grid(row=6, column=0, sticky='w')
        tk.Radiobutton(frame_1, indicatoron=0, text="Cylindrical", justify=tk.LEFT, command=FormChoice, variable=f, value=2).grid(row=7, column=0, sticky='w')

        # Interface position
        tk.Label(frame_1, text="Interface", font=('times', 20, 'bold')).grid(row=0, column=7)
        tk.Label(frame_1, text="Position", font=('times', 20, 'bold'), padx=10).grid(row=0, column=8)
        tk.Button(frame_1, text="Draw interface", command=DrawInterfaceLine).grid(row=4, column=7)
        tk.Button(frame_1, text="Delete", command=DeletePunchLine).grid(row=4, column=8)
        tk.Button(frame_1, text="Draw P-Line 1 & 2", command=DrawLine12).grid(row=5, column=7)
        tk.Button(frame_1, text="Delete", command=DeleteDrawLine12).grid(row=5, column=8)
        tk.Button(frame_1, text="Draw N-Line 1 & 2", command=Draw_NLine12).grid(row=6, column=7)
        tk.Button(frame_1, text="Delete", command=DeleteDraw_NLine12).grid(row=6, column=8)

        global interface_list
        interface_list = []
        for i in range(len(default_interface)):
            j = i + 1
            tk.Label(frame_1, text=default_interface[i]).grid(row=j, column=7)
            interfaces = tk.Entry(frame_1, width=4)
            interfaces.grid(row=j, column=8)
            interfaces.insert(0, default_IntPosition[i])
            interface_list.append(interfaces)

        # Voltage distribution
        tk.Label(frame_1, text="V region", font=('times', 20, 'bold')).grid(row=0, column=4)
        tk.Label(frame_1, text="Numbers", font=('times', 20, 'bold')).grid(row=0, column=6)
        global Vmesh_number, Vmesh_lineEnt
        Vmesh_number = []
        Vmesh_lineEnt = []
        V1_meshlineEnt = []
        for i in range(1, 6):
            V1 = tk.Label(frame_1, text=str(default_Vmeshlines[i - 1]), width=2)
            V2 = tk.Entry(frame_1, width=4)
            V1.grid(row=i, column=3)
            V2.grid(row=i, column=5, sticky='w')
            V2.insert(0, default_Vmeshlines[i])
            tk.Label(frame_1, text="<  V%s  <" % i).grid(row=i, column=4)
            Vmesh_temp = tk.Entry(frame_1, width=3)
            Vmesh_temp.grid(row=i, column=6)
            Vmesh_temp.insert(0, default_Vnumbers[i - 1])
            Vmesh_number.append(Vmesh_temp)
            V1_meshlineEnt.append(V1)
            Vmesh_lineEnt.append(V2)
        tk.Button(frame_1, text="Draw", command=Vmeshlineplot).grid(row=6, column=4)
        tk.Button(frame_1, text="Delete", command=deleteVmeshlines).grid(row=7, column=4)

        tk.Label(frame_1, text="Xmin", font=('times', 16)).grid(row=0, column=1, sticky='e')
        tk.Label(frame_1, text="Xmax(<Xj)", font=('times', 16)).grid(row=1, column=1, sticky='e')
        tk.Label(frame_1, text="VmaxError", font=('times', 16)).grid(row=2, column=1, sticky='e')
        tk.Label(frame_1, text="V(Goal)", font=('times', 16)).grid(row=3, column=1, sticky='e')
        tk.Label(frame_1, text="Tolerance", font=('times', 16)).grid(row=4, column=1, sticky='e')
        tk.Label(frame_1, text="Trials", font=('times', 16)).grid(row=5, column=1, sticky='e')
        tk.Label(frame_1, text="Vb-Step", font=('times', 16)).grid(row=6, column=1, sticky='e')
        tk.Label(frame_1, text="V-number", font=('times', 16)).grid(row=7, column=1, sticky='e')
        XminEnt = tk.Entry(frame_1, width=7)
        XmaxEnt = tk.Entry(frame_1, width=7)
        VmaxErrorEntry = tk.Entry(frame_1, width=7)
        VgoalEnt = tk.Entry(frame_1, width=7)
        TolEnt = tk.Entry(frame_1, width=7)
        TrialsEnt = tk.Entry(frame_1, width=7)
        StepEnt = tk.Entry(frame_1, width=7)
        VnumberEnt = tk.Entry(frame_1, width=7)

        global othersettings
        othersettings = [XminEnt, XmaxEnt, VmaxErrorEntry, VgoalEnt, TolEnt, TrialsEnt, StepEnt, VnumberEnt]

        XminEnt.grid(row=0, column=2, padx=10)
        XmaxEnt.grid(row=1, column=2)
        VmaxErrorEntry.grid(row=2, column=2)
        VgoalEnt.grid(row=3, column=2)
        TolEnt.grid(row=4, column=2)
        TrialsEnt.grid(row=5, column=2)
        StepEnt.grid(row=6, column=2)
        VnumberEnt.grid(row=7, column=2)
        XminEnt.insert(0, default_others[0])
        XmaxEnt.insert(0, default_others[1])
        VmaxErrorEntry.insert(0, default_others[2])
        VgoalEnt.insert(0, default_others[3])
        TolEnt.insert(0, default_others[4])
        TrialsEnt.insert(0, round(default_others[5]))
        StepEnt.insert(0, default_others[6])
        VnumberEnt.insert(0, round(default_others[7]))

        s = ttk.Style()
        s.theme_use('clam')
        s.configure("red.Horizontal.TProgressbar", foreground='red', background='green')

        Vpt_Log = tk.Text(frame_42, height=5.5, width=51, highlightbackground='grey')
        Vpt_Log.pack(side='left', fill='y')
        VptScroll = tk.Scrollbar(frame_42)
        VptScroll.pack(side='right', fill='y')
        VptScroll.config(command=Vpt_Log.yview)
        Vpt_Log.config(yscrollcommand=VptScroll.set)
        Vpt_Log.insert(tk.END, "------------Punch-through voltage result-----------\n")

        Vb_Log = tk.Text(frame_44, height=5.5, width=51, highlightbackground='grey')
        Vb_Log.pack(side='left', fill='y')
        VptScroll = tk.Scrollbar(frame_44)
        VptScroll.pack(side='right', fill='y')
        VptScroll.config(command=Vb_Log.yview)
        Vb_Log.config(yscrollcommand=VptScroll.set)
        Vb_Log.insert(tk.END, "--------------Breakdown voltage result-------------\n")

class Page4(Page):
   def __init__(self, *args, **kwargs):
       Page.__init__(self, *args, **kwargs)
       frame = tk.Frame(self)
       frame.pack()
       frame_1 = tk.Frame(frame, pady=20)
       frame_2 = tk.Frame(frame)
       frame_3 = tk.Frame(frame)
       frame_4 = tk.Frame(frame)
       frame_1.pack(side="top", fill='x')
       frame_3.pack(side="top", fill='x')
       frame_4.pack(side="bottom", fill='x')
       frame_2.pack(side="bottom", fill='x')

       def DrawAny():
           SetDefault()

           temp_xdata = XArrayEnt.get()
           temp_ydata = YArrayEnt.get()
           temp_xdata = temp_xdata.rsplit(',')
           temp_ydata = temp_ydata.rsplit('||')
           temp_Vpt_fit = temp_ydata[0].rsplit(',')
           temp_Vpt_TCAD = temp_ydata[1].rsplit(',')
           temp_Vb_fit = temp_ydata[2].rsplit(',')
           temp_Vb_TCAD = temp_ydata[3].rsplit(',')

           if len(temp_Vpt_fit) == len(temp_Vpt_TCAD) \
                   and len(temp_Vpt_TCAD) == len(temp_Vb_fit) \
                   and len(temp_Vb_fit) == len(temp_Vb_TCAD) \
                   and len(temp_Vb_TCAD) == len(temp_xdata):

               AnyXData = [float(i) for i in temp_xdata]
               temp_Vpt_fit = [float(i) for i in temp_Vpt_fit]
               temp_Vpt_TCAD = [float(i) for i in temp_Vpt_TCAD]
               temp_Vb_fit = [float(i) for i in temp_Vb_fit]
               temp_Vb_TCAD = [float(i) for i in temp_Vb_TCAD]

               if temp_DrawAny[0] == 0:
                   VptTCAD_line, = ax3.plot(AnyXData, temp_Vpt_TCAD, '-.bo', markersize=3, label=r'$V_{pt}$ (TCAD)')
                   VptFit_line, = ax3.plot(AnyXData, temp_Vpt_fit, '-bo', markersize=3, label=r'$V_{pt}$ (fit)')
                   VbTCAD_line, = ax3.plot(AnyXData, temp_Vb_TCAD, '-.ro', markersize=3, label=r'$V_{b}$ (TCAD)')
                   VbFit_line, = ax3.plot(AnyXData, temp_Vb_fit, '-ro', markersize=3, label=r'$V_{b}$ (fit)')
                   temp_DrawAny[0] = VptTCAD_line
                   temp_DrawAny[1] = VptFit_line
                   temp_DrawAny[2] = VbTCAD_line
                   temp_DrawAny[3] = VbFit_line
                   XaxisMin3 = min(AnyXData) - (max(AnyXData) - min(AnyXData)) / 6
                   XaxisMax3 = max(AnyXData) + (max(AnyXData) - min(AnyXData)) / 6
                   YaxisMin3 = 0
                   YaxisMax3 = max(temp_Vb_fit) + max(temp_Vb_fit) / 6
                   ax3.set_xlim((XaxisMin3, XaxisMax3))
                   ax3.set_ylim((YaxisMin3, YaxisMax3))
                   ax3.legend(loc='center left')
                   bar3.draw()
                   log_region.insert(tk.END, "\nSuccessfully plot!")
                   log_region.see(tk.END)
               else:
                   log_region.insert(tk.END, "\nRemoving plot...")
                   log_region.see(tk.END)
                   temp_DrawAny[0].remove()
                   temp_DrawAny[1].remove()
                   temp_DrawAny[2].remove()
                   temp_DrawAny[3].remove()
                   log_region.insert(tk.END, "\nSuccessfully remove plot!")
                   log_region.insert(tk.END, "\nCreate new plot...!")
                   log_region.see(tk.END)
                   VptTCAD_line, = ax3.plot(AnyXData, temp_Vpt_TCAD, '-.bo', markersize=3, label=r'$V_{pt}$ (TCAD)')
                   VptFit_line, = ax3.plot(AnyXData, temp_Vpt_fit, '-bo', markersize=3, label=r'$V_{pt}$ (fit)')
                   VbTCAD_line, = ax3.plot(AnyXData, temp_Vb_TCAD, '-.ro', markersize=3, label=r'$V_{b}$ (TCAD)')
                   VbFit_line, = ax3.plot(AnyXData, temp_Vb_fit, '-ro', markersize=3, label=r'$V_{b}$ (fit)')
                   temp_DrawAny[0] = VptTCAD_line
                   temp_DrawAny[1] = VptFit_line
                   temp_DrawAny[2] = VbTCAD_line
                   temp_DrawAny[3] = VbFit_line
                   XaxisMin3 = min(AnyXData) - (max(AnyXData) - min(AnyXData)) / 6
                   XaxisMax3 = max(AnyXData) + (max(AnyXData) - min(AnyXData)) / 6
                   YaxisMin3 = 0
                   YaxisMax3 = max(temp_Vb_fit) + max(temp_Vb_fit) / 6
                   ax3.set_xlim((XaxisMin3, XaxisMax3))
                   ax3.set_ylim((YaxisMin3, YaxisMax3))
                   ax3.legend(loc='center left')
                   bar3.draw()
                   log_region.insert(tk.END, "\nSuccessfully plot!")
                   log_region.see(tk.END)
           else:
               tM.showerror("Error", "Each X and Y data must have different length.")


       def ClearAny():
           if temp_DrawAny[0] == 0:
               tM.showerror("Error", "You haven't plotted any figure.")
           else:
               for i in range(4):
                   temp_DrawAny[i].remove()
                   temp_DrawAny[i] = 0
               bar3.draw()
               log_region.insert(tk.END, "\nSuccessfully remove %s interface lines!" % len(default_interface))
               log_region.see(tk.END)


       global tv
       tv = ttk.Treeview(frame_1, height=5, selectmode='browse')
       StatisticsDataList = ['Date', 'Time', 'APD', 'Node', 'ND',
                        'Nc', 'Vpt', 'Vpt(T)', 'Vb', 'Vb(T)', '|Error|',
                        'Ring', 'Form', 'Coordinate',
                        'OtherInfo.', 'CurveName', 'FilePath']
       StatisticsDataListWidth = [105, 80, 40, 50, 75, 75, 73, 73, 77, 77, 60, 45, 45, 70, 65, 45, 45]
       tv['columns'] = StatisticsDataList
       tv.heading("#0", text='No.', anchor='w')
       tv.column("#0", anchor="w", width=40)
       for index, element in enumerate(StatisticsDataList):
           tv.heading(element, text=element)
           tv.column(element, anchor='center', width=StatisticsDataListWidth[index])
       tv.pack(side='left')
       LoadTreeviewData()
       tv.bind('<ButtonRelease-1>', selectTreeItem)
       tv.bind('<Double-Button>', importsettingsDoubleClick)

       vsb = ttk.Scrollbar(frame_1, orient="vertical", command=tv.yview)
       vsb.pack(side='right', fill='y')
       tv.configure(yscrollcommand=vsb.set)

       tk.Label(frame_3, text="VptFit   ||   VptTCAD   ||   VbFit   ||   VbTCAD", font=('Times', 20, 'bold')).grid(row=0, column=0, columnspan=3)
       tk.Label(frame_3, text="Concentration: ", font=('Times', 16)).grid(row=1, column=0, sticky='e')

       global XArrayEnt, YArrayEnt, PlotDataArray

       XArrayEnt = tk.Entry(frame_3, width=70)
       XArrayEnt.grid(row=1, column=1)
       XArrayEnt.insert(0, default_PlotValue[0])
       tk.Label(frame_3, text="Voltage: ", font=('Times', 16)).grid(row=2, column=0, sticky='e')
       YArrayEnt = tk.Entry(frame_3, width=70)
       YArrayEnt.grid(row=2, column=1)
       YArrayEnt.insert(0, default_PlotValue[1])
       PlotDataArray = [XArrayEnt, YArrayEnt]
       DrawButton = tk.Button(frame_3, text="Draw", command=DrawAny)
       DrawButton.grid(row=1, column=2, padx=20)
       ClearDrawButton = tk.Button(frame_3, text="Clear", command=ClearAny)
       ClearDrawButton.grid(row=2, column=2, padx=20)

       global ax3, bar3, temp_DrawAny

       temp_DrawAny = [0] * 4
       fig3 = Figure(figsize=(18, 6.1), dpi=80, facecolor='w', edgecolor='k')
       ax3 = fig3.add_subplot(111)
       ax3.set_xlabel(r'X')
       ax3.set_ylabel(r'Y')
       ax3.grid(True)
       bar3 = FigureCanvasTkAgg(fig3, frame_2)
       bar3.get_tk_widget().pack(side='right')
       tkagg.NavigationToolbar2Tk(bar3, frame_4)


class MainView(tk.Frame):
    def __init__(self, *args, **kwargs):
        tk.Frame.__init__(self, *args, **kwargs)
        p1 = Page1(self)
        p2 = Page2(self)
        p3 = Page3(self)
        p4 = Page4(self)

        buttonframe = tk.Frame(self)
        container = tk.Frame(self)
        buttonframe.pack(side="top", fill="x", expand=False)
        container.pack(side="top", fill="both", expand=True)

        p1.place(in_=container, x=0, y=0, relwidth=1, relheight=1)
        p2.place(in_=container, x=0, y=0, relwidth=1, relheight=1)
        p3.place(in_=container, x=0, y=0, relwidth=1, relheight=1)
        p4.place(in_=container, x=0, y=0, relwidth=1, relheight=1)

        b1 = tk.Button(buttonframe, text="User Guide", command=p1.lift)
        b2 = tk.Button(buttonframe, text="Doping profile", command=p2.lift)
        b3 = tk.Button(buttonframe, text="Breakdown Voltage", command=p3.lift)
        b4 = tk.Button(buttonframe, text="Statistics", command=p4.lift)
        b5 = tk.Button(buttonframe, text="Default", command=SetDefault)
        b6 = tk.Button(buttonframe, text="Save & Default", command=Save)
        b7 = tk.Button(buttonframe, text="Import", command=importsettings)

        b1.pack(side="left")
        b2.pack(side="left")
        b3.pack(side="left")
        b4.pack(side="left")
        b5.pack(side="right")
        b6.pack(side="right")
        b7.pack(side="right")

        p2.show()

if __name__ == "__main__":
    root = tk.Tk()
    # left frame
    LeftFrame = tk.Frame(root, width=150, height=300, borderwidth=5)
    LeftFrame.pack(side="left", fill="both", expand=True)
    LeftFrame.pack_propagate(0)
    # left frame top
    LeftFrame_top = tk.Frame(LeftFrame, height=50)
    LeftFrame_top.pack(side="top", fill="x")
    # left frame log region
    LeftFrame_LOG = tk.Frame(LeftFrame, height=100)
    LeftFrame_LOG.pack(side="top", fill="x")
    # left frame Read/Load Buttons
    LeftFrame_Load_title = tk.Frame(LeftFrame, height=20)
    LeftFrame_Load_title.pack(side="top", fill="x")
    LeftFrame_Load = tk.Frame(LeftFrame, height=50)
    LeftFrame_Load.pack(side="top", fill="x")
    # left frame variables window
    LeftFrame_variables_title = tk.Frame(LeftFrame, height=20)
    LeftFrame_variables_title.pack(side="top", fill="x")
    LeftFrame_variables = tk.Frame(LeftFrame, height=50)
    LeftFrame_variables.pack(side="top", fill="x")

    # right frame
    RightFrame = tk.Frame(root, width=1050, height=500, borderwidth=5)
    RightFrame.pack(side="right", fill="both", expand=True)
    RightFrame.pack_propagate(0)

    # Load right frame
    main = MainView(RightFrame)
    main.pack(side="top", fill="both", expand=True)

    # Main root size
    root.wm_geometry("1600x900")
    root.title("Avalanche Photodiode analysis tool")

    # Left frame content
    tk.Label(LeftFrame_top, text="DopingProfile.csv", padx=20, font=('times', 25, 'bold')).grid(row=0, column=0, columnspan=2)
    csv_path = tk.Text(LeftFrame_top, height=5, width=35, highlightbackground='grey')
    csv_path.grid(row=1, column=0, columnspan=2)
    tk.Button(LeftFrame_top, text='File Open', command=openfile).grid(row=2, column=0)
    tk.Button(LeftFrame_top, text='Load file', command=retrieve_input).grid(row=2, column=1)

    tk.Label(LeftFrame_LOG, text="Log", padx=20, font=('times', 25, 'bold')).pack()
    log_region = tk.Text(LeftFrame_LOG, height=15, width=34, highlightbackground='grey')
    log_region.pack(side='left', fill='y')
    Scroll = tk.Scrollbar(LeftFrame_LOG)
    Scroll.pack(side='right', fill='y')
    Scroll.config(command=log_region.yview)
    log_region.config(yscrollcommand=Scroll.set)
    log_region.insert(tk.END, "--------Python Log region--------\n")

    # Load data list
    tk.Label(LeftFrame_Load_title, text="Data List", padx=20, font=('times', 25, 'bold')).pack()
    LoadListBox = tk.Listbox(LeftFrame_Load, height=5, width=27)
    LoadListBox.pack(side="left", fill="y")
    LoadScrollbar = tk.Scrollbar(LeftFrame_Load, orient="vertical")
    LoadScrollbar.config(command=LoadListBox.yview)
    LoadScrollbar.pack(side="right", fill="y")
    LoadListBox.config(yscrollcommand=LoadScrollbar.set)
    LoadListBox.bind("<<ListboxSelect>>", onselect)
    LoadListBox.bind("<Double-Button>", readcurveDoubleClick)

    # Variables
    varTitle = tk.Label(LeftFrame_variables_title, text="Variable History", padx=20, font=('times', 25, 'bold'))
    varTitle.pack()
    VarListBox = tk.Listbox(LeftFrame_variables, height=10, width=27)
    VarListBox.pack(side="left", fill="y")
    VarScrollbar = tk.Scrollbar(LeftFrame_variables, orient="vertical")
    VarScrollbar.config(command=VarListBox.yview)
    VarScrollbar.pack(side="right", fill="y")
    VarListBox.config(yscrollcommand=VarScrollbar.set)
    VarListBox.bind("<<ListboxSelect>>", VarSelect)

    while True:
        try:
            root.mainloop()
        except UnicodeDecodeError:
            continue
        break