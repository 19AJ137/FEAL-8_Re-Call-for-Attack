"""Verify a FEAL-8 key against the challenge data (written independently of feal8.h).
Usage: python3 verify_key.py <key hex> plaintext.bin ciphertext.bin
Checks the test vector (key 0123456789ABCDEF, plaintext 0 -> CEEF2C86F2490752), then encrypts
all 1024 blocks in ECB mode and compares them with the ciphertext file."""
import sys
def rot2(x): return ((x<<2)|(x>>6))&0xff
def S0(a,b): return rot2((a+b)&0xff)
def S1(a,b): return rot2((a+b+1)&0xff)
def f(a,b):
    f1=a[1]^b[0]^a[0]; f2=a[2]^b[1]^a[3]
    f1=S1(f1,f2); f2=S0(f2,f1)
    return [S0(a[0],f1),f1,f2,S1(a[3],f2)]
def fK(a,b):
    fk1=a[1]^a[0]; fk2=a[2]^a[3]
    fk1=S1(fk1,fk2^b[0]); fk2=S0(fk2,fk1^b[1])
    return [S0(a[0],fk1^b[2]),fk1,fk2,S1(a[3],fk2^b[3])]
def sched(key):
    A=list(key[0:4]); B=list(key[4:8]); D=[0,0,0,0]; K=[]
    for r in range(8):
        Bn=fK(A,[B[i]^D[i] for i in range(4)]); K+=[Bn[0:2],Bn[2:4]]; D,A,B=A,B,Bn
    return K
def enc(K,p):
    L=[p[i]^K[8+i//2][i%2] for i in range(4)]; R=[p[4+i]^K[10+i//2][i%2] for i in range(4)]
    R=[R[i]^L[i] for i in range(4)]
    for r in range(8):
        t=f(R,K[r]); L,R=R,[L[i]^t[i] for i in range(4)]
    L=[L[i]^R[i] for i in range(4)]; o=R+L
    return bytes([o[i]^K[12+i//2][i%2] for i in range(4)]+[o[4+i]^K[14+i//2][i%2] for i in range(4)])
assert enc(sched(bytes.fromhex("0123456789abcdef")),bytes(8)).hex()=="ceef2c86f2490752"
key=bytes.fromhex(sys.argv[1]); K=sched(key)
P=open(sys.argv[2],"rb").read(); C=open(sys.argv[3],"rb").read()
ok=sum(enc(K,P[8*i:8*i+8])==C[8*i:8*i+8] for i in range(1024))
print(f"key {key.hex()}: {ok}/1024 blocks match")
print("block0", P[0:8].hex(), "->", enc(K,P[0:8]).hex(), "file", C[0:8].hex())
