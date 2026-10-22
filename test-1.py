import numpy as np

d = np.load("pile_1000_1.npz")
C, U, L, R = d["centers"], d["directions"], d["lengths"], d["radii"]
P = C - 0.5 * L[:, None] * U
Q = C + 0.5 * L[:, None] * U
D1 = Q - P
A = (D1 * D1).sum(1)
n = len(L)

def min_gap(i):
    # distance axe-axe min à toutes les autres fibres, vectorisé (RTCD 5.1.9)
    p1, d1 = P[i], D1[i]
    rv = p1 - P
    e = (D1 * D1).sum(1)
    f = (D1 * rv).sum(1)
    cc = D1 @ d1
    b = rv @ d1
    denom = A[i] * e - cc * cc
    s = np.zeros(n)
    ok = denom > 1e-18
    s[ok] = np.clip((b[ok] * e[ok] - cc[ok] * f[ok]) / denom[ok], 0.0, 1.0)
    t = (s * cc - f) / e
    lo, hi = t < 0, t > 1
    t[lo] = 0.0; s[lo] = np.clip(-b[lo] / A[i], 0.0, 1.0)
    t[hi] = 1.0; s[hi] = np.clip((cc[hi] - b[hi]) / A[i], 0.0, 1.0)
    dist = np.linalg.norm(p1 + s[:, None] * d1 - (P + t[:, None] * D1), axis=1)
    gap = dist - (R[i] + R)
    gap[i] = np.inf
    return gap.min()

zmin = C[:, 2] - 0.5 * L * np.abs(U[:, 2])
SEUIL = 0.04  # mm ; contact Bullet = pénétration -20/+11 µm, donc 40 µm = sûr
floaters = [i for i in range(n) if min_gap(i) > SEUIL and zmin[i] > SEUIL]
print(f"{len(floaters)} fibres sans contact et au-dessus du sol (seuil 40 µm)")
for i in floaters[:20]:
    print(f"  fibre {i:4d} : x={C[i,0]:6.2f} y={C[i,1]:6.2f} zmin={zmin[i]*1e3:7.1f} µm")
