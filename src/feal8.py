"""Reference FEAL-8 implementation (Shimizu & Miyaguchi, 1987) and the ISO/DIS 9797-style CBC-MAC
used to check the transcribed challenge data."""


def rot2(x):
    return ((x << 2) | (x >> 6)) & 0xFF


def S(x, y, d):
    return rot2((x + y + d) & 0xFF)


def f(a, b):
    """f function: a = 4 bytes, b = 2 bytes"""
    a0, a1, a2, a3 = a
    b0, b1 = b
    f1 = a1 ^ b0 ^ a0
    f2 = a2 ^ b1 ^ a3
    f1 = S(f1, f2, 1)
    f2 = S(f2, f1, 0)
    f0 = S(a0, f1, 0)
    f3 = S(a3, f2, 1)
    return [f0, f1, f2, f3]


def fK(a, b):
    """fK function of the key schedule: a, b = 4 bytes"""
    a0, a1, a2, a3 = a
    b0, b1, b2, b3 = b
    k1 = a1 ^ a0
    k2 = a2 ^ a3
    k1 = S(k1, k2 ^ b0, 1)
    k2 = S(k2, k1 ^ b1, 0)
    k0 = S(a0, k1 ^ b2, 0)
    k3 = S(a3, k2 ^ b3, 1)
    return [k0, k1, k2, k3]


def xor(x, y):
    return [p ^ q for p, q in zip(x, y)]


def key_schedule(key):
    """64-bit key -> sixteen 16-bit extended keys K0..K15 (2 bytes each)"""
    A, B, D = list(key[:4]), list(key[4:]), [0, 0, 0, 0]
    sub = []
    for _ in range(8):
        Bn = fK(A, xor(B, D))
        D, A, B = A, B, Bn
        sub.append(Bn[0:2])
        sub.append(Bn[2:4])
    return sub


def encrypt(key, pt, sk=None):
    K = sk or key_schedule(key)
    L, R = list(pt[:4]), list(pt[4:])
    L = xor(L, K[8] + K[9])
    R = xor(R, K[10] + K[11])
    R = xor(R, L)
    for r in range(8):
        L, R = R, xor(L, f(R, K[r]))
    # output transformation: (R8, L8)
    L, R = R, L
    R = xor(R, L)
    L = xor(L, K[12] + K[13])
    R = xor(R, K[14] + K[15])
    return bytes(L + R)


def cbc_mac(key, data, iv=bytes(8)):
    sk = key_schedule(key)
    x = iv
    for i in range(0, len(data), 8):
        x = encrypt(key, bytes(a ^ b for a, b in zip(data[i:i + 8], x)), sk)
    return x
