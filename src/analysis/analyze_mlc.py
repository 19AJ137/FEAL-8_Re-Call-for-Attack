import sys, numpy as np
def load(fn):
    rows=[]
    for line in open(fn):
        if line.startswith('#'): continue
        a,b,c,d=line.split('|')
        x=a.split(); r=c.split()
        rows.append([float(v) for v in x[3:9]]+[float(v) for v in r]+[d.strip()=='ok'])
    return np.array(rows)
for fn in sys.argv[1:]:
    D=load(fn); n=len(D)
    print(f"{fn}: {n} keys, consistent {int(D[:,10].sum())}/{n}")
    names=['MLC (paper, 8 approx)','fam0 only (4 approx)','single approx 1','quadratic 8 approx']
    print("  stat                     [1,2^8) [2^8,2^19) [2^19,2^24) [2^24,2^37]   <=2^24  <=2^26  <=2^28  <=2^30  median log2 rank")
    for k,nm in enumerate(names):
        r=D[:,6+k]
        b=[np.sum(r<8),np.sum((r>=8)&(r<19)),np.sum((r>=19)&(r<24)),np.sum(r>=24)]
        print(f"  {nm:24s} {b[0]:6d} {b[1]:9d} {b[2]:10d} {b[3]:11d}   {np.mean(r<24):6.0%} {np.mean(r<26):6.0%} {np.mean(r<28):6.0%} {np.mean(r<30):6.0%}   {np.median(r):.1f}")
    zA0,zA1,zB0,zB1,p0,p1=[D[:,i] for i in range(6)]
    cur=(abs(zB0)>=3.5)&(abs(p0)>=2)&(abs(zA0)>=7.3)
    new=(np.hypot(zB0,zB1)>=3.9175)&(np.hypot(p0,p1)>=2.486)&(np.hypot(zA0,zA1)>=8.5)
    mlc24=D[:,6]<24
    print(f"  same keys, fwd only: 方式C current (>=7.3) {cur.mean():.0%}, 2-family {new.mean():.0%}, paper MLC R=2^24 {mlc24.mean():.0%}")
    print(f"  paper MLC succeeds where current 方式C fails: {np.sum(mlc24&~cur)}/{np.sum(~cur)};  where 2-family fails: {np.sum(mlc24&~new)}/{np.sum(~new)}")
