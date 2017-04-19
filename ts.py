#coding=utf-8
#time smoothimg
import numpy as np

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
