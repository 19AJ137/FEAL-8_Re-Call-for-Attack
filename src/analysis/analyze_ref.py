"""方式C：現行（第1族のみ）と第2族追加の通過率を比べる。
ref_real.txt（refstat の出力）と、段階2の誤り候補数モデル（attackC 実測で較正）を使う。"""
import numpy as np
from scipy import stats

d = np.loadtxt("ref_real.txt", usecols=range(1, 13))
F, B = d[:, 0:6], d[:, 6:12]  # zA0 zA1 zB0 zB1 pre0 pre1
N = len(d)
c = 0.2004  # 8ビット分割の容量。誤り仮説では corr(zA, zB) = sqrt(c)
r1_35 = np.sqrt(-2 * np.log(2 * stats.norm.sf(3.5)))   # |zB0|>=3.5 と同じ生存数
r1_30 = np.sqrt(-2 * np.log(2 * stats.norm.sf(3.0)))
rpre = np.sqrt(-2 * np.log(2 * stats.norm.sf(2.0)))    # 先頭256組の 2σ と同じ通過率
print(f"keys {N}; r1(3.5)={r1_35:.4f} r1(3.0)={r1_30:.4f} rpre={rpre:.4f}")

def corr(a, b):
    return np.corrcoef(a, b)[0, 1]
print("correct key, fwd: mean|zA0| %.3f mean|zA1| %.3f mean|zB0| %.3f mean|zB1| %.3f" %
      tuple(np.abs(F[:, :4]).mean(0)))
print("correct key, fwd: corr(zA0,zA1)=%.3f corr(|zB0|,|zB1|)=%.3f corr(zB0,zB1)=%.3f" %
      (corr(F[:, 0] * np.sign(F[:, 1]), np.abs(F[:, 1])), corr(np.abs(F[:, 2]), np.abs(F[:, 3])), corr(F[:, 2], F[:, 3])))
print("correct key: |zA0| fwd vs rev identical in %.1f%% of keys; |zA1| fwd = |zA0| rev in %.1f%%; |zA1| fwd = |zA1| rev in %.1f%%" %
      (100 * np.mean(np.isclose(np.abs(F[:, 0]), np.abs(B[:, 0]))),
       100 * np.mean(np.isclose(np.abs(F[:, 1]), np.abs(B[:, 0]))),
       100 * np.mean(np.isclose(np.abs(F[:, 1]), np.abs(B[:, 1])))))

# --- 誤り候補数のモデル ---
NS = 2.0**45 * 2 * stats.norm.sf(3.5)  # 段階1の生存仮説（全空間、1方向）
COMP = 2.0**15                          # 補完（DK7 残り8bit × DK6 残り7bit）

def false_cur(thr, thr1=3.5, m=400000):
    # 生存仮説の zB0 ~ N(0,1) | |zB0|>=thr1、zA0 = sqrt(c) zB0 + sqrt(1-c) ξ
    rng = np.random.default_rng(1)
    u = rng.uniform(0, 1, m)
    zb = stats.norm.isf(u * stats.norm.sf(thr1))
    mu = np.sqrt(c) * zb
    s = np.sqrt(1 - c)
    p = stats.norm.sf((thr - mu) / s) + stats.norm.cdf((-thr - mu) / s)
    ns = 2.0**45 * 2 * stats.norm.sf(thr1)
    return ns * COMP * p.mean()

def false_new(r2, r1=r1_35, m=400000):
    # 生存仮説の |zB| ~ Rayleigh | >= r1、|zA|^2/(1-c) ~ ncx2(2, c|zB|^2/(1-c))
    rng = np.random.default_rng(2)
    rb2 = r1**2 + rng.exponential(2.0, m)
    p = stats.ncx2.sf(r2**2 / (1 - c), 2, c * rb2 / (1 - c))
    ns = 2.0**45 * np.exp(-r1**2 / 2)
    return ns * COMP * p.mean()

print("\nfalse candidates per direction (model):")
for t in (6.5, 7.0, 7.3, 7.5):
    print(f"  current |zA0|>={t}: {false_cur(t):.3g}   (attackC measured: 6.5→2e7, 7.0→7.5e5, 7.3→1e5, 7.5→2.5e4)")
for r2 in (7.5, 8.0, 8.5, 9.0, 9.5):
    print(f"  new sqrt(zA0^2+zA1^2)>={r2}: {false_new(r2):.3g}")

def cur(X, thr1=3.5, thrA=6.5):
    return (np.abs(X[:, 2]) >= thr1) & (np.abs(X[:, 4]) >= 2.0) & (np.abs(X[:, 0]) >= thrA)

def new(X, r1=r1_35, r2=8.5):
    s1 = np.hypot(X[:, 2], X[:, 3]) >= r1
    pre = np.hypot(X[:, 4], X[:, 5]) >= rpre
    s2 = np.hypot(X[:, 0], X[:, 1]) >= r2
    return s1 & pre & s2

def show(name, f, b):
    print(f"  {name:44s} fwd {f.mean():6.1%}  rev {b.mean():6.1%}  either {(f | b).mean():6.1%}")

print("\nsuccess of stages 1-2 (correct key; stage 3 not included):")
show("current thr1=3.5, |zA0|>=6.5", cur(F), cur(B))
show("current thr1=3.5, |zA0|>=7.3 (sequential set)", cur(F, thrA=7.3), cur(B, thrA=7.3))
show("current thr1=3.0, |zA0|>=6.5", cur(F, 3.0), cur(B, 3.0))
show("step1 only: |zB0|>=3.5", np.abs(F[:, 2]) >= 3.5, np.abs(B[:, 2]) >= 3.5)
show("step1 only: hypot(zB0,zB1)>=3.917", np.hypot(F[:, 2], F[:, 3]) >= r1_35, np.hypot(B[:, 2], B[:, 3]) >= r1_35)
for r2 in (8.0, 8.5, 9.0):
    show(f"new r1=3.917 (same survivors), r2={r2}", new(F, r2=r2), new(B, r2=r2))
show("new r1=3.440 (same as thr 3.0), r2=8.5", new(F, r1_30, 8.5), new(B, r1_30, 8.5))
