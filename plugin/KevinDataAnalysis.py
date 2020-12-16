import numpy as np
import KevinLibrary as elite
import matplotlib.pyplot as plt  # 畫圖必備
import matplotlib.pylab as pylab
params = {'legend.fontsize': 'x-large',
          'figure.figsize': (20, 10),
          'axes.labelsize': 'x-large',
          'axes.titlesize': 'x-large',
          'xtick.labelsize': 'x-large',
          'ytick.labelsize': 'x-large'}
pylab.rcParams.update(params)
plt.rcParams.update({'font.size': 9})


class Crate(object):
    def __init__(self, ti, tf, x, temperature, width, height, mass, t_interval=0.1):
        """
        物件初始化
        :param ti: 初時刻
        :param tf: 末時刻
        :param x: 位置陣列
        :param temperature: 溫度
        :param width: 寬度
        :param height: 高度
        :param mass: 質量
        :param t_interval: 時距
        """
        self.mass = mass

        # 只要給定ti, tf, t_interval，就利用 np.arange() 製作時間陣列，再把此陣列賦予 self.t
        self.t = np.arange(ti, tf, t_interval)

        self.x = x
        self.T = temperature
        self.width = width
        self.height = height

    def velocity(self, t=None):
        if t is None:
            return elite.dydx(self.t, self.x)
        else:
            v = elite.dydx(self.t, self.x)
            return elite.find_y_from_x(self.t, v, t)

    def acceleration(self, t=None):
        """
        由給定的 x(t) 計算 a(t)
        :param t: 可以計算特定時間 t 下的加速度。但也可以不給特定時間，而輸出所有時刻的 a(t)
        :return: 回傳 a(t) 陣列，如 [3,5,1,2,6,8,1....]
        """
        if t is None:
            v = self.velocity()
            return elite.dydx(self.t, v)
        else:
            v = self.velocity()
            a = elite.dydx(self.t, v)
            return elite.find_y_from_x(self.t, a, t)

    def displacement(self, ti, tf):
        if ti < self.t[0]:
            raise BaseException("時間(%s)太早了！必須大於 %s" % (ti, self.t[0]))
        if tf > self.t[-1]:
            raise BaseException("時間(%s)太晚了！必須小於 %s" % (tf, self.t[-1]))

        xi = elite.find_y_from_x(self.t, self.x, ti)
        xf = elite.find_y_from_x(self.t, self.x, tf)
        return xf - xi

    def average_velocity(self, ti, tf):
        dx = self.displacement(ti, tf)
        dt = tf - ti
        return dx / dt

    def force(self, t=None):
        if t is None:
            return self.mass * self.acceleration()
        else:
            a_t = elite.find_y_from_x(self.t, self.acceleration(), t)
            return self.mass * a_t

    def plot_xt(self, fig_num, color, label, linestyle='-'):
        plt.figure(fig_num)
        plt.plot(self.t, self.x, color=color, linestyle=linestyle, label=label)
        plt.grid()
        plt.xlabel('t (s)')
        plt.ylabel('x (cm)')
        plt.legend(loc='best')

