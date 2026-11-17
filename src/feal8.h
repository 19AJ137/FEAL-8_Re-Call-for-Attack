/* FEAL-8 と等価表現・線形近似の共通部分
 * 32 ビット語は左端（バイト0）を最上位に置く。 */
#ifndef FEAL8_H
#define FEAL8_H
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

typedef uint8_t u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;

static inline u8 rol2(u8 x) { return (u8)((x << 2) | (x >> 6)); }
static inline u8 S0(u8 x, u8 y) { return rol2((u8)(x + y)); }
static inline u8 S1(u8 x, u8 y) { return rol2((u8)(x + y + 1)); }
static inline u8 BY(u32 x, int i) { return (u8)(x >> (24 - 8 * i)); }
static inline u32 PK(u8 a, u8 b, u8 c, u8 d) {
    return ((u32)a << 24) | ((u32)b << 16) | ((u32)c << 8) | d;
}

/* f(α, β)：β は 16 ビット */
static inline u32 feal_f(u32 a, u16 b) {
    u8 a0 = BY(a, 0), a1 = BY(a, 1), a2 = BY(a, 2), a3 = BY(a, 3);
    u8 f1 = a1 ^ (u8)(b >> 8) ^ a0, f2 = a2 ^ (u8)b ^ a3;
    f1 = S1(f1, f2);
    f2 = S0(f2, f1);
    return PK(S0(a0, f1), f1, f2, S1(a3, f2));
}

/* N：S 関数4個の部分。F(A,K) = N(M(A) ^ (0,K,0)) */
static inline u32 feal_N(u32 x) {
    u8 x0 = BY(x, 0), x1 = BY(x, 1), x2 = BY(x, 2), x3 = BY(x, 3);
    u8 y1 = S1(x1, x2), y2 = S0(x2, y1);
    return PK(S0(x0, y1), y1, y2, S1(x3, y2));
}
static inline u32 feal_M(u32 a) {
    u8 a0 = BY(a, 0), a1 = BY(a, 1), a2 = BY(a, 2), a3 = BY(a, 3);
    return PK(a0, a0 ^ a1, a2 ^ a3, a3);
}

static inline u32 fK(u32 a, u32 b) {
    u8 a0 = BY(a, 0), a1 = BY(a, 1), a2 = BY(a, 2), a3 = BY(a, 3);
    u8 k1 = a1 ^ a0, k2 = a2 ^ a3;
    k1 = S1(k1, k2 ^ BY(b, 0));
    k2 = S0(k2, k1 ^ BY(b, 1));
    return PK(S0(a0, k1 ^ BY(b, 2)), k1, k2, S1(a3, k2 ^ BY(b, 3)));
}

typedef struct { u16 k[16]; } sched_t;

static void key_schedule(u64 key, sched_t *s) {
    u32 A = (u32)(key >> 32), B = (u32)key, D = 0;
    for (int r = 0; r < 8; r++) {
        u32 Bn = fK(A, B ^ D);
        D = A; A = B; B = Bn;
        s->k[2 * r] = (u16)(Bn >> 16);
        s->k[2 * r + 1] = (u16)Bn;
    }
}

/* 暗号化。st が非 NULL なら st[2*i], st[2*i+1] に (L_i, R_i)（i=0..8）を残す */
static u64 feal8_enc(const sched_t *s, u64 p, u32 *st) {
    u32 L = (u32)(p >> 32), R = (u32)p;
    L ^= ((u32)s->k[8] << 16) | s->k[9];
    R ^= ((u32)s->k[10] << 16) | s->k[11];
    R ^= L;
    if (st) { st[0] = L; st[1] = R; }
    for (int r = 0; r < 8; r++) {
        u32 t = L ^ feal_f(R, s->k[r]);
        L = R; R = t;
        if (st) { st[2 * r + 2] = L; st[2 * r + 3] = R; }
    }
    u32 CL = R, CR = L ^ R;
    CL ^= ((u32)s->k[12] << 16) | s->k[13];
    CR ^= ((u32)s->k[14] << 16) | s->k[15];
    return ((u64)CL << 32) | CR;
}

static inline int par(u32 x) { return __builtin_parity(x); }

/* 乱数（xoshiro256**） */
static u64 rs[4];
static inline u64 rotl64(u64 x, int k) { return (x << k) | (x >> (64 - k)); }
static inline u64 rnd(void) {
    u64 r = rotl64(rs[1] * 5, 7) * 9, t = rs[1] << 17;
    rs[2] ^= rs[0]; rs[3] ^= rs[1]; rs[1] ^= rs[2]; rs[0] ^= rs[3];
    rs[2] ^= t; rs[3] = rotl64(rs[3], 45);
    return r;
}
static void seed(u64 x) {
    for (int i = 0; i < 4; i++) {
        x += 0x9e3779b97f4a7c15ULL;
        u64 z = x;
        z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
        rs[i] = z ^ (z >> 31);
    }
}

/* 線形近似のマスク（第1〜5段、偏り 2^-3）と第6段 */
#define G_L1 0x00010504u
#define G_R1 0x04031004u
#define G_L6 0x04031004u
#define G_R6 0x00010504u
#define G6_Y 0x04031004u /* 第6段 F 出力マスク */
#define G_L7 0x10105050u /* 6段近似 Approximation 1 の (L7,R7) マスク */
#define G_R7 0x04031004u

/* 8ビット分割：第6段 N の入力 w1 のビット2〜4、w2 のビット2〜6 */
static inline int part_idx(u8 w1, u8 w2) { return ((w1 >> 2) & 7) | (((w2 >> 2) & 31) << 3); }
static double e_tab[256];
static void build_etab(void) {
    double acc[256] = {0};
    int cnt[256] = {0};
    for (int w1 = 0; w1 < 256; w1++)
        for (int w2 = 0; w2 < 256; w2++) {
            int g = part_idx((u8)w1, (u8)w2);
            /* w0_0, w3_0 は単独では分からないので、既知の A1_0 = w0_0^w1_0、A2_0 = w3_0^w2_0 に吸収させる */
            acc[g] += (par(G6_Y & feal_N(PK(0, (u8)w1, (u8)w2, 0))) ^ (w1 & 1) ^ (w2 & 1)) ? -1 : 1;
            cnt[g]++;
        }
    for (int g = 0; g < 256; g++) e_tab[g] = acc[g] / cnt[g];
}

static u64 *load_unique(const char *fn, int *n) {
    FILE *fp = fopen(fn, "rb");
    if (!fp) { perror(fn); exit(1); }
    u8 buf[8192];
    size_t got = fread(buf, 1, sizeof buf, fp);
    fclose(fp);
    u64 *v = (u64 *)malloc(sizeof(u64) * 1024);
    int m = 0;
    for (size_t i = 0; i + 8 <= got; i += 8) {
        u64 x = 0;
        for (int j = 0; j < 8; j++) x = (x << 8) | buf[i + j];
        int dup = 0;
        for (int j = 0; j < m; j++) if (v[j] == x) { dup = 1; break; }
        if (!dup) v[m++] = x;
    }
    *n = m;
    return v;
}
#endif
