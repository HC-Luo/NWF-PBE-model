#coding=utf-8
#正压原始方程模式 primitive barotropic equation model
#2017 Apr 19 leonidas

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import subroutine as sub

NX  = 20       #为x方向格点数
NY  = 16       #为y方向格点数
d   = 300000.0 #网格距
cla = 51.0     #区域中心纬度
clo = 118.0    #区域中心经度
dt  = 600.0    #时间步长
zo  = 2500.0   #是为了减小重力惯性外波的波速，增加差分格式的稳定性而引入的位势高度
s   = 0.5      #平滑系数
nt2 = 72       #用于判别是否积分12小时，是否该做内点平滑
nt4 = 6        #用于判定是否该做边界平滑
nt5 = 36       #用于判定是否该做时间平滑
c1  = dt/2.0
c2  = dt*2.0
ur = np.zeros((NX,NY),float) #初始的x方向风速
vr = np.zeros((NX,NY),float) #初始的y方向风速
zr = np.zeros((NX,NY),float) #初始的位势高度
ua = np.zeros((NX,NY),float) #n-1时间层的x方向风速
va = np.zeros((NX,NY),float) #n-1时间层的y方向风速
za = np.zeros((NX,NY),float) #n-1时间层的位势高度
ub = np.zeros((NX,NY),float) #n时间层的x方向风速
vb = np.zeros((NX,NY),float) #n时间层的y方向风速
zb = np.zeros((NX,NY),float) #n时间层的位势高度
uc = np.zeros((NX,NY),float) #n+1时间层的x方向风速
vc = np.zeros((NX,NY),float) #n+1时间层的y方向风速
zc = np.zeros((NX,NY),float) #n+1时间层的位势高度
rm = np.zeros((NX,NY),float) #放大系数
f  = np.zeros((NX,NY),float) #地转参数
w  = np.zeros((NX,NY),float) #工作数组

#===============================读入初始资料场====================================
filename1 = '/Users/leonidas/AnacondaProjects/sx/NWPsx/za.dat'
fileRead1 = file(filename1,'r')
for i in range(NX):
    for j in range(NY):
        za[i,j] = float(fileRead1.readline()[6*i:5+6*i])
        zr[i,j] = za[i,j]
    fileRead1.seek(0)

#=============================计算放大系数和地转参数================================
rm,f = sub.cmf(rm,f,d,cla,NX,NY)

#================================计算地转风初值===================================
ua,va = sub.cgw(ua,va,za,rm,f,d,NX,NY)
for i in range(NX):
    for j in range(NY):
        ur[i,j] = ua[i,j]
        vr[i,j] = va[i,j]
'''
#this may not be used

filename2 = '/Users/leonidas/AnacondaProjects/sx/NWPsx/ua.dat'
filename3 = '/Users/leonidas/AnacondaProjects/sx/NWPsx/va.dat'
fileRead2 = file(filename2,'r')
fileRead3 = file(filename3,'r')
for i in range(NX):
    for j in range(NY):
        ua[i,j] = float(fileRead2.readline()[6*i:5+6*i])
    fileRead2.seek(0)
'''

#==================================边值传送======================================
ub,vb,zb = sub.tbv(ub,vb,zb,ua,va,za,NX,NY)
uc,vc,zc = sub.tbv(uc,vc,zc,ua,va,za,NX,NY)

#==================================开始预报======================================
for na in [1,2]:
    nb = 0
    #欧拉后差积分1小时
    for nn in range(6):
        ub,vb,zb = sub.ti(ua,va,za,ua,va,za,ub,vb,zb,rm,f,d,dt,zo,NX,NY)
        ua,va,za = sub.ti(ua,va,za,ub,vb,zb,ua,va,za,rm,f,d,dt,zo,NX,NY)
        nb = nb+1
        continue

    #边界平滑子程序
    za = sub.ssbp(za,w,s,NX,NY)
    ua = sub.ssbp(ua,w,s,NX,NY)
    va = sub.ssbp(va,w,s,NX,NY)

    #前差积分半步
    ub,vb,zb = sub.ti(ua,va,za,ua,va,za,ub,vb,zb,rm,f,d,c1,zo,NX,NY)
    #中央差积分半步
    uc,vc,zc = sub.ti(ua,va,za,ub,vb,zb,uc,vc,zc,rm,f,d,dt,zo,NX,NY)
    #数组传送子程序
    ub,vb,zb = sub.ta(ub,vb,zb,uc,vc,zc,NX,NY)
    #中央差积分一步,共积分11小时
    for nn in range(66):
        uc,vc,zc = sub.ti(ua,va,za,ub,vb,zb,uc,vc,zc,rm,f,d,c2,zo,NX,NY)
        nb = nb+1
        #打印积分步数,na大循环步,nb小循环步
        print na,nb
        #判断是否积分12小时
        if nb==nt2:
            break
        #判断是否做边界平滑
        if nb/nt4*nt4==nb:
            zc = sub.ssbp(zc,w,s,NX,NY)
            uc = sub.ssbp(uc,w,s,NX,NY)
            vc = sub.ssbp(vc,w,s,NX,NY)
        #判断是否做时间平滑
        if nb==nt5:
            ub,vb,zb = sub.ts(ua,ub,uc,va,vb,vc,za,zb,zc,s,NX,NY)
        if nb==nt5+1:
            ub,vb,zb = sub.ts(ua,ub,uc,va,vb,vc,za,zb,zc,s,NX,NY)
        #数组传送,为下一轮积分做准备
        ua,va,za = sub.ta(ua,va,za,ub,vb,zb,NX,NY)
        ub,vb,zb = sub.ta(ub,vb,zb,uc,vc,zc,NX,NY)

        continue

    #区域内点平滑
    zc = sub.ssip(zc,w,s,NX,NY,2)
    uc = sub.ssip(uc,w,s,NX,NY,2)
    vc = sub.ssip(vc,w,s,NX,NY,2)
    #打印积分步数
    print na,nb
    #数组传送,为后12小时的积分做准备
    ua,va,za = sub.ta(ua,va,za,uc,vc,zc,NX,NY)

#===================================作 图=======================================
plt.figure(figsize = [18,10])
scale = '50m'
lim=[85,151.5,32.5,70]
xstep = 5
ystep = 5
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_extent(lim,crs=ccrs.PlateCarree())
land = cfeature.NaturalEarthFeature('physical', 'land', scale,edgecolor='face',facecolor=cfeature.COLORS['land'])
ax.add_feature(land, facecolor='0.75')
ax.coastlines(scale)
lon_formatter = LongitudeFormatter(zero_direction_label=False)
lat_formatter = LatitudeFormatter()
ax.xaxis.set_major_formatter(lon_formatter)
ax.yaxis.set_major_formatter(lat_formatter)
ax.set_xticks(np.arange(lim[0],lim[1]+0.5*xstep,xstep), crs=ccrs.PlateCarree())
ax.set_yticks(np.arange(lim[2]+2.5,lim[3]+0.5*ystep+2,ystep), crs=ccrs.PlateCarree())

dd = 50
lev = np.arange(5200,5800+dd,dd)
x = np.arange(85,155,3.5)
y = np.arange(32.5,72.5,2.5)
lons, lats = np.meshgrid(x, y)
pic = ax.contour(x, y, zc.T,levels = lev,cmap='rainbow',linewidths = 3)
plt.colorbar(pic,shrink = 0.8)
#plt.savefig('/Users/leonidas/Desktop/After_NWP.pdf')
plt.show()
