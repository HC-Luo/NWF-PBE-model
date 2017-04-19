#coding=utf-8
#time integrations
import numpy as np

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
