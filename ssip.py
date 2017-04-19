#coding=utf-8
#space smoothing for internal points 区域内5点平滑(正逆平滑)
#l=1为只执行正平滑，l=2为执行正逆平滑.

import numpy as np

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
