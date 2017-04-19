#coding=utf-8
import numpy as np

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
