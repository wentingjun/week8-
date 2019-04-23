# 第一题（拟合巨磁实验）
# -*- coding: utf-8 -*-
import numpy as np
from scipy.optimize import leastsq
import pylab as pl

def func(x, p):
   # 数据拟合所用的函数: A*sin(2*pi*k*x + theta)
    A, k, theta = p
    return A*np.sin(2*np.pi*k*x+theta)   

def residuals(p, y, x):
    #实验数据x, y和拟合函数之间的差，p为拟合需要找到的系数
    return y - func(x, p)

x = np.linspace(0, 48*np.pi/180, 17)  #转为弧度制
A, k, theta = 10, 0.34, np.pi/6 # 真实数据的函数参数
y = [48.5,52.6,27.0,-13.8,-38.0,-29.5,-4.9,25.2,48.6,53.2,26.7,-16.1,-39.4,-29.9,-3.5,25.2,48.5] # 真实数据

p0 = [55, 0.2, 3/np.pi] # 第一次猜测的函数拟合参数

# 调用leastsq进行数据拟合
# residuals为计算误差的函数
# p0为拟合参数的初始值
# args为需要拟合的实验数据
plsq = leastsq(residuals, p0, args=(y0, x))

print("真实参数:", [A, k, theta] )
print("拟合参数", plsq[0]) # 实验数据拟合后的参数
pl.title('《巨磁》实验中的实验五曲线函数') 
pl.xlabel('转动角度/弧度') 
pl.ylabel('输出电压/mV')
pl.rcParams['font.family'] = ['simHei']  #中文正常显示
pl.rcParams['font.sans-serif'] = ['simHei'] # 步骤一（替换sans-serif字体）
pl.rcParams['axes.unicode_minus'] = False   # 步骤二（解决坐标轴负数的负号显示问题）
pl.plot(x, y, label=u"真实数据")
pl.plot(x, func(x, plsq[0]), label=u"拟合数据")
pl.legend()
pl.show()

# -*- coding: utf-8 -*-
import numpy as np
from scipy.optimize import leastsq
import pylab as pl
from scipy.optimize import curve_fit
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

pl.rcParams['font.family'] = ['simHei']
pl.rcParams['font.sans-serif'] = ['simHei'] # 步骤一（替换sans-serif字体）
pl.rcParams['axes.unicode_minus'] = False   # 步骤二（解决坐标轴负数的负号显示问题）
pl.title("温度数据的曲线拟合")
pl.xlabel("月份")
pl.ylabel("温度/摄氏度")
pl.grid()  # 生成网格

x=np.arange(1,13,1)
xmajorLocator  = MultipleLocator(1) #将x主刻度标签设置为1的倍数
ymax=np.array([17, 19, 21, 28, 33, 38, 37, 37, 31, 23, 19, 18 ])  #每月温度最高值真实数据
ymin=np.array([-62, -59, -56, -46, -32, -18, -9, -13, -25, -46, -52, -58]) #每月温度最低值真实数据

def fmax(x,a,b,c):
    #数据拟合所用的函数: a*np.sin(x*np.pi/6+b)+c
    return a*np.sin(x*np.pi/6+b)+c


def fmin(x,a,b,c):
    #数据拟合所用的函数: a*np.sin(x*np.pi/6+b)+c
    return a*np.sin(x*np.pi/6+b)+c

fita,fitb=curve_fit(fmax,x,ymax)
fitc,fitd=curve_fit(fmin,x,ymin)
print(fita,fitb)
print(fitc,fitd)

pl.plot(x,ymax,label=u"最大值真实数据",color='m',linestyle='-',marker='')#绘制每月温度最大值真实数据
pl.plot(x,fmax(x,fita[0],fita[1],fita[2]),label=u"最大值拟合数据",color='g',linestyle='-',marker='')#绘制每月温度最大值拟合曲线

pl.plot(x,ymin,label=u"最小值真实数据",color='b',linestyle='-',marker='')#绘制每月温度最小值真实数据
pl.plot(x,fmin(x,fitc[0],fitc[1],fitc[2]),label=u"最小值拟合数据",color='r',linestyle='-',marker='')#绘制每月温度最小值拟合曲线

pl.legend()
pl.show()