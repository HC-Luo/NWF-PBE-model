#coding=utf-8
import numpy as np
import math
#===============computing map factors and coriolis parameter====================
'''
rk为圆锥常数,rlq为兰勃特投影映像平面上赤道到北极点的距离,a为地球半径
sita为标准余纬,psx为区域中心余纬,r为模式中心到北极的距离
'''
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

#=======================calcualting geostrophic winds===========================
def cgw(ua,va,za,rm,f,d,NX,NY):
    m = NX
    n = NY
    g = 9.8
    for i in np.arange(0,m):
        for j in np.arange(1,n-1):
            ua[i,j] = -1.0*rm[i,j]*g/f[i,j]*((za[i,j+1]-za[i,j-1])/(2*d))
        ua[i,0] = -1.0*rm[i,0]*g/f[i,0]*((za[i,1]-za[i,0])/d)
        ua[i,n-1] = -1.0*rm[i,n-1]*g/f[i,n-1]*((za[i,n-1]-za[i,n-2])/d)

    for j in np.arange(0,n):
        for i in np.arange(1,m-1):
            va[i,j] =  1.0*rm[i,j]*g/f[i,j]*((za[i+1,j]-za[i-1,j])/(2*d))
        va[0,j] =  1.0*rm[0,j]*g/f[0,j]*((za[1,j]-za[0,j])/d)
        va[m-1,j] =  1.0*rm[m-1,j]*g/f[m-1,j]*((za[m-1,j]-za[m-2,j])/d)
    print 'End calcualting ua, va'
    return ua,va

#================space smoothing for boundary points 边界九点平滑=================
def ssbp(a,w,s,NX,NY):
    m  =  NX
    n  =  NY
    m1 = m-1
    m2 = m-2
    n1 = n-1
    n2 = n-2
    for i in np.arange(1,m1):
        for j in [1,n2]:
            w[i,j] = a[i,j]+0.5*s*(1.0-s)\
                    *(a[i-1,j]+a[i+1,j]+a[i,j-1]+a[i,j+1]-4.0*a[i,j])+0.25*s*s\
                    *(a[i-1,j-1]+a[i-1,j+1]+a[i+1,j-1]+a[i+1,j+1]-4.0*a[i,j])

    for i in [1,m2]:
        for j in np.arange(2,n2):
            w[i,j] = a[i,j]+0.5*s*(1.0-s)\
                    *(a[i-1,j]+a[i+1,j]+a[i,j-1]+a[i,j+1]-4.0*a[i,j])+0.25*s*s\
                    *(a[i-1,j-1]+a[i-1,j+1]+a[i+1,j-1]+a[i+1,j+1]-4.0*a[i,j])

    for i in np.arange(1,m1):
        for j in [1,n2]:
            a[i,j] = w[i,j]

    for i in [1,m2]:
        for j in np.arange(2,n2):
            a[i,j] = w[i,j]

    return a

#=============space smoothing for internal points 区域内5点平滑(正逆平滑)==========
#l=1为只执行正平滑，l=2为执行正逆平滑.
def ssip(a,w,s,NX,NY,l):
    m = NX
    n = NY
    m2 = m-2
    n2 = n-2
    if l==1:
        for i in np.arange(2,m2):
            for j in np.arange(2,n2):
                w[i,j] = a[i,j]+s/4.0*(a[i-1,j]+a[i+1,j]+a[i,j-1]+a[i,j+1]-4.0*a[i,j])

        for i in np.arange(2,m2):
            for j in np.arange(2,n2):
                a[i,j] = w[i,j]
    elif l==2:
        for i in np.arange(2,m2):
            for j in np.arange(2,n2):
                w[i,j] = a[i,j]+s/4.0*(a[i-1,j]+a[i+1,j]+a[i,j-1]+a[i,j+1]-4.0*a[i,j])

        for i in np.arange(2,m2):
            for j in np.arange(2,n2):
                a[i,j] = w[i,j]

        s = -0.5

        for i in np.arange(2,m2):
            for j in np.arange(2,n2):
                w[i,j] = a[i,j]+s/4.0*(a[i-1,j]+a[i+1,j]+a[i,j-1]+a[i,j+1]-4.0*a[i,j])

        for i in np.arange(2,m2):
            for j in np.arange(2,n2):
                a[i,j] = w[i,j]

    else:
        print 'Space smoothing ERROR!'

    return a

#===========================transmiting arrays  数组传送=========================
def ta(ua,va,za,ub,vb,zb,NX,NY):
    m = NX
    n = NY
    for i in range(m):
        for j in range(n):
            ua[i,j] = ub[i,j]
            va[i,j] = vb[i,j]
            za[i,j] = zb[i,j]

    return ua,va,za

#====================transmiting boundary valaus  赋固定边界值====================
def tbv(ua,va,za,ub,vb,zb,NX,NY):
    m = NX
    n = NY
    m1 = m-1
    n1 = n-1
    for i in range(m):
        for j in [0,n1]:
            ua[i,j] = ub[i,j]
            va[i,j] = vb[i,j]
            za[i,j] = zb[i,j]
    for i in [0,m1]:
        for j in range(n):
            ua[i,j] = ub[i,j]
            va[i,j] = vb[i,j]
            za[i,j] = zb[i,j]

    return ua,va,za

#================================time integrations==============================
def ti(ua,va,za,ub,vb,zb,uc,vc,zc,rm,f,d,dt,zo,NX,NY):
    m = NX
    n = NY
    c = 0.25/d
    m1 = m-1
    n1 = n-1
    for i in np.arange(1,m1):
        for j in np.arange(1,n1):
            e = -c*rm[i,j]*((ub[i+1,j]+ub[i,j])*(ub[i+1,j]-ub[i,j])\
                +(ub[i,j]+ub[i-1,j])*(ub[i,j]-ub[i-1,j])\
                +(vb[i,j-1]+vb[i,j])*(ub[i,j]-ub[i,j-1])\
                +(vb[i,j]+vb[i,j+1])*(ub[i,j+1]-ub[i,j])\
                +19.6*(zb[i+1,j]-zb[i-1,j]))+f[i,j]*vb[i,j]

            uc[i,j] = ua[i,j]+e*dt

            g = -c*rm[i,j]*((ub[i+1,j]+ub[i,j])*(vb[i+1,j]-vb[i,j])\
                +(ub[i,j]+ub[i-1,j])*(vb[i,j]-vb[i-1,j])\
                +(vb[i,j-1]+vb[i,j])*(vb[i,j]-vb[i,j-1])\
                +(vb[i,j]+vb[i,j+1])*(vb[i,j+1]-vb[i,j])\
                +19.6*(zb[i,j+1]-zb[i,j-1]))-f[i,j]*ub[i,j]

            vc[i,j] = va[i,j]+g*dt

    for i in np.arange(1,m1):
        for j in np.arange(1,n1):
            h = -c*rm[i,j]*((ub[i+1,j]+ub[i,j])*(zb[i+1,j]-zb[i,j])\
                +(ub[i,j]+ub[i-1,j])*(zb[i,j]-zb[i-1,j])\
                +(vb[i,j-1]+vb[i,j])*(zb[i,j]-zb[i,j-1])\
                +(vb[i,j]+vb[i,j+1])*(zb[i,j+1]-zb[i,j])\
                +2.0*(zb[i,j]-zo)*(ub[i+1,j]-ub[i-1,j]+vb[i,j+1]-vb[i,j-1]))

            zc[i,j] = za[i,j]+h*dt

    return uc,vc,zc

#================================time smoothimg=================================
def ts(ua,ub,uc,va,vb,vc,za,zb,zc,s,NX,NY):
    m = NX
    n = NY
    m1 = m-1
    n1 = n-1
    for i in np.arange(1,m1):
        for j in np.arange(1,n1):
            ub[i,j]=ub[i,j]+s*(ua[i,j]+uc[i,j]-2.0*ub[i,j])/2.0
            vb[i,j]=vb[i,j]+s*(va[i,j]+vc[i,j]-2.0*vb[i,j])/2.0
            zb[i,j]=zb[i,j]+s*(za[i,j]+zc[i,j]-2.0*zb[i,j])/2.0

    return ub,vb,zb
