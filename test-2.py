import numpy as np

d = np.load("pile_1000_1.npz")
C, U, L, R = d["centers"], d["directions"], d["lengths"], d["radii"]
n = len(L)
if not np.allclose(np.linalg.norm(U, axis=1), 1.0, atol=1e-6):
    print("ATTENTION : directions non unitaires, zmin fausse")

P = C - 0.5 * L[:, None] * U
Q = C + 0.5 * L[:, None] * U
D1 = Q - P
A = (D1 * D1).sum(1)

def min_gap(i):
    p1, d1 = P[i], D1[i]
    rv = p1 - P
    f = (D1 * rv).sum(1)
    cc = D1 @ d1
    b = rv @ d1
    denom = A[i] * A - cc * cc
    s = np.zeros(n)
    ok = denom > 1e-18
    s[ok] = np.clip((b[ok] * A[ok] - cc[ok] * f[ok]) / denom[ok], 0.0, 1.0)
    t = (s * cc - f) / A
    lo, hi = t < 0, t > 1
    t[lo] = 0.0; s[lo] = np.clip(-b[lo] / A[i], 0.0, 1.0)
    t[hi] = 1.0; s[hi] = np.clip((cc[hi] - b[hi]) / A[i], 0.0, 1.0)
    dist = np.linalg.norm(p1 + s[:, None] * d1 - (P + t[:, None] * D1), axis=1)
    gap = dist - (R[i] + R)
    gap[i] = np.inf
    return gap.min()

gaps = np.array([min_gap(i) for i in range(n)])
zmin = C[:, 2] - 0.5 * L * np.abs(U[:, 2])
SEUIL = 0.04
fl = np.where((gaps > SEUIL) & (zmin > SEUIL))[0]

print(f"tas       : zmin min {zmin.min()*1e3:7.1f} µm | médiane {np.median(zmin)*1e3:7.1f} | max {zmin.max()*1e3:7.1f}")
print(f"témoin    : gap médian des fibres en contact = {np.median(gaps[gaps <= 0])*1e3:.1f} µm")
print(f"flottants : {len(fl)} / {n}")
if len(fl):
    q = np.percentile(zmin[fl], [0, 25, 50, 75, 100]) * 1e3
    print(f"  altitude (µm)      : min {q[0]:.0f} / q25 {q[1]:.0f} / méd {q[2]:.0f} / q75 {q[3]:.0f} / max {q[4]:.0f}")
    q = np.percentile(gaps[fl], [0, 50, 100]) * 1e3
    print(f"  gap au voisin (µm): min {q[0]:.0f} / méd {q[1]:.0f} / max {q[2]:.0f}")
    hist, _ = np.histogram(fl, bins=np.arange(0, n + 100, 100))
    print(f"  flottants par tranche de 100 indices : {hist.tolist()}")
