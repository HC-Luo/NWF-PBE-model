#coding=utf-8
#computing map factors and coriolis parameter
'''
rk为圆锥常数,rlq为兰勃特投影映像平面上赤道到北极点的距离,a为地球半径
sita为标准余纬,psx为区域中心余纬,r为模式中心到北极的距离
'''
import numpy as np
import math

def cmf(rm,f,d,cla,NX,NY):
    m = NX
    n = NY
    rk=0.7156
    rlq=11423370.0
    a=6371000.0
    conv=57.29578
    w1=2.0/rk
    sita=30.0/conv
    psx=(90.0-cla)/conv

    #计算模式中心到北极的距离r
    cel0=a*math.sin(sita)/rk
    cel=(math.tan(psx/2.0))/(math.tan(sita/2.0))
    r=cel0*cel**rk

    #确定网格坐标原点在地图坐标系中的位置
    xi0=-(m-1)/2.0
    yj0=r/d+(n-1)/2.0

    #各格点至北极点的距离rl,(xj,yi)为模式各格点在地图坐标系中的位置
    for i in range(m):
        for j in range(n):
            xi=xi0+(i-1)
            yj=yj0-(j-1)
            rl=math.sqrt(xi**2+yj**2)*d
            #求放大系数rm和柯氏参数f
            w2=(rl/rlq)**w1
            sinl=(1.0-w2)/(1.0+w2)
            rm[i,j]=rk*rl/(a*math.sqrt(1.0-sinl**2))
            f[i,j] = (1.4584e-4)*sinl
    print 'End computing map factors and coriolis parameter'
    return rm,f
