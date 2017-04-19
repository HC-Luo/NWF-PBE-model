#coding=utf-8
#space smoothing for boundary points 边界九点平滑
import numpy as np

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
