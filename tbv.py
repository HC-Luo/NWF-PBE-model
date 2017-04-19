#coding=utf-8
#transmiting boundary valaus  赋固定边界值

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
