/* 崎小山ら（IEICE 2016）付録の 6段近似 8本
 * l = fam*4 + x。fam 0: ΓY2=04031104, ΓY6=04031004 / fam 1: ΓY2=04011904, ΓY6=04011804
 * x: ΓX6 = 10117d7c, 18197574, 18195d5c, 10115554（ε = 97, 65, 97, 65 × 2^-12）
 * (L1,R1) のマスク = (00010504, ΓY2^00000100)、(L7,R7) のマスク = (ΓX6^00010504, ΓY6) */
#ifndef APPROX8_H
#define APPROX8_H
static const u32 A8_Y2[2] = {0x04031104u, 0x04011904u};
static const u32 A8_Y6[2] = {0x04031004u, 0x04011804u};
static const u32 A8_X6[4] = {0x10117d7cu, 0x18197574u, 0x18195d5cu, 0x10115554u};
static const int A8_EPS[4] = {97, 65, 97, 65};
#define A8_L1 0x00010504u
static inline u32 a8_R1(int fam) { return A8_Y2[fam] ^ 0x00000100u; }
static inline u32 a8_L7(int x) { return A8_X6[x] ^ 0x00010504u; }
static inline u32 a8_R7(int fam) { return A8_Y6[fam]; }
/* 真の中間値 (L1,R1,L7,R7) から近似 l のビット */
static inline int a8_bit(int l, u32 L1, u32 R1, u32 L7, u32 R7) {
    int fam = l >> 2, x = l & 3;
    return par(A8_L1 & L1) ^ par(a8_R1(fam) & R1) ^ par(a8_L7(x) & L7) ^ par(a8_R7(fam) & R7);
}
#endif
