import numpy as np
def _ss(Y,g):
    gm=Y.mean(0); s=0.0
    for lv in np.unique(g):
        m=Y[g==lv]; s+=len(m)*((m.mean(0)-gm)**2).sum()
    return s
def _resid(Y,g):
    R=Y.copy()
    for lv in np.unique(g):
        R[g==lv]-=Y[g==lv].mean(0)
    return R
def adonis2(Y,f1,f2,nperm=9999,seed=1):
    """Euclidea, terminos secuenciales (equivale a vegan::adonis2 by='terms')."""
    Y=np.asarray(Y,float); n=len(Y)
    SST=((Y-Y.mean(0))**2).sum()
    def stats(f1,f2):
        SS1=_ss(Y,f1)
        R1=_resid(Y,f1); SS2=_ss(R1,f2)
        SSR=SST-SS1-SS2
        d1=len(np.unique(f1))-1; d2=len(np.unique(f2))-1; dr=n-1-d1-d2
        return SS1,SS2,SSR,d1,d2,dr
    SS1,SS2,SSR,d1,d2,dr=stats(f1,f2)
    F1=(SS1/d1)/(SSR/dr); F2=(SS2/d2)/(SSR/dr)
    rng=np.random.default_rng(seed); c1=c2=0
    for _ in range(nperm):
        p=rng.permutation(n)
        a,b,r,_,_,_=stats(f1[p],f2[p])
        if (a/d1)/(r/dr)>=F1: c1+=1
        if (b/d2)/(r/dr)>=F2: c2+=1
    return dict(R2_f1=SS1/SST,F_f1=F1,p_f1=(c1+1)/(nperm+1),
                R2_f2=SS2/SST,F_f2=F2,p_f2=(c2+1)/(nperm+1))
def betadisper(Y,g,nperm=9999,seed=1):
    """Distancia al centroide por grupo, ANOVA con permutaciones."""
    Y=np.asarray(Y,float); d=np.zeros(len(Y))
    for lv in np.unique(g):
        m=Y[g==lv]; d[g==lv]=np.sqrt(((m-m.mean(0))**2).sum(1))
    gm=d.mean(); dfb=len(np.unique(g))-1; dfw=len(Y)-1-dfb
    def F(d,g):
        b=sum((g==lv).sum()*(d[g==lv].mean()-d.mean())**2 for lv in np.unique(g))
        w=sum(((d[g==lv]-d[g==lv].mean())**2).sum() for lv in np.unique(g))
        return (b/dfb)/(w/dfw)
    F0=F(d,g); rng=np.random.default_rng(seed); c=0
    for _ in range(nperm):
        if F(d,rng.permutation(g))>=F0: c+=1
    return dict(F=F0,p=(c+1)/(nperm+1))
