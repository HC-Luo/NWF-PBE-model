#coding=utf-8
#transmiting arrays  数组传送

def ta(ua,va,za,ub,vb,zb,NX,NY):
    m = NX
    n = NY
    for i in range(m):
        for j in range(n):
            ua[i,j] = ub[i,j]
            va[i,j] = vb[i,j]
            za[i,j] = zb[i,j]

    return ua,va,za
