"""Check data/plaintext.bin and data/ciphertext.bin against the CBC-MAC values
printed with the challenge (data/mac_fig8.txt and data/mac_fig9.txt).

Usage: python3 check_mac.py [data_dir]

MAC rule (as determined by matching all 64 printed values):
  K = 0, IV = 0, x_i = FEAL-8_K(D_i xor x_{i-1}), chained over all 1024 blocks.
  After every 32 blocks (x_32, x_64, ..., x_1024) the left 32 bits of x are printed.
"""
import os
import sys

from feal8 import encrypt, key_schedule


def macs(data):
    K = bytes(8)
    sk = key_schedule(K)
    x = bytes(8)
    out = []
    for i in range(0, len(data), 8):
        x = encrypt(K, bytes(a ^ b for a, b in zip(data[i:i + 8], x)), sk)
        if (i // 8 + 1) % 32 == 0:
            out.append(x.hex()[:8])
    return out


def main():
    d = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(__file__), "..", "data")
    assert encrypt(bytes.fromhex("0123456789abcdef"), bytes(8)).hex() == "ceef2c86f2490752"
    ok_all = True
    for name, fig in (("plaintext.bin", "mac_fig8.txt"), ("ciphertext.bin", "mac_fig9.txt")):
        data = open(os.path.join(d, name), "rb").read()
        exp = open(os.path.join(d, fig)).read().split()
        got = macs(data)
        ok = sum(g == e for g, e in zip(got, exp))
        print(f"{name}: {len(data)} bytes, CBC-MAC {ok}/{len(exp)} match {fig}")
        for k, (g, e) in enumerate(zip(got, exp), 1):
            if g != e:
                print(f"  #{k:02d}: computed {g}, printed {e}")
        ok_all &= ok == len(exp) == 32
    sys.exit(0 if ok_all else 1)


if __name__ == "__main__":
    main()
