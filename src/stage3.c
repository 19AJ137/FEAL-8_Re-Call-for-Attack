/* 方式C 段階3：段階2の候補 (EK0 15bit, DK7 30bit, DK6 15bit) から 64 ビット鍵を復元する
 *
 *   ./stage3 self plaintext.bin seed [cands.txt]   乱数鍵で検証（候補ファイルが無ければ真の候補だけ）
 *   ./stage3 real plaintext.bin ciphertext.bin cands.txt
 *
 * 手順（いずれも 5 段近似（相関 1/4）を1段ずらして使う）
 *   (a) 第2〜6段の近似：EK0 のバイト0・3（14bit）と正規化ビット、EK1 のバイト1・2（15bit）
 *   (c) 第0〜4段の近似：DK6 のバイト0・3（14bit）と正規化ビット、DK5 のバイト1・2（15bit）
 *   (t) 第2〜4段の 3 段近似（相関 1/2、鍵の推測なし）で (a)(c) の組を判定
 *   (b) 第3〜7段の近似：EK1 のバイト0・3（14bit）と正規化ビット、EK2 のバイト1・2
 *   (e) EK0, EK1, DK6, DK7 と鍵スケジュールから 64 ビット鍵（未知 32bit × 最上位ビットの曖昧さ 8bit を全数）
 * (a)(b)(c) は2段構え：分割統計（8bit 分割、WHT で κ 全通り）→ 生存した仮説を補完して全統計量
 */
#include "feal8.h"
#include <omp.h>
#include <sys/time.h>

#define W 16
static int n;
static u32 PLa[1024], PRa[1024], CLa[1024], CRa[1024];
static float Ehat[256], E2hat[256];
static u64 ROW[65536][2];
static double now(void) { struct timeval tv; gettimeofday(&tv, 0); return tv.tv_sec + tv.tv_usec * 1e-6; }
static inline int hprime(int w1, int w2) {
    return par(G6_Y & feal_N(PK(0, (u8)w1, (u8)w2, 0))) ^ (w1 & 1) ^ (w2 & 1);
}
static inline int partx(u32 x) { return part_idx((u8)(x >> 8), (u8)x); }
static inline u32 emb_part(int kp) { return ((u32)(kp & 7) << 10) | ((u32)((kp >> 3) & 31) << 2); }
static inline u32 emb_rest(int r) {
    u32 w1 = (r & 3) | (((r >> 2) & 7) << 5), w2 = (r >> 5) & 3;
    return (w1 << 8) | w2;
}
static inline u32 dk7_22(u32 t) { return PK((u8)((t >> 16) & 3), (u8)t, (u8)(t >> 8), (u8)((t >> 18) & 15)); }
static inline u32 dk7_ext(int e8) { return ((u32)(e8 & 31) << 26) | ((u32)((e8 >> 5) & 7) << 4); }
static inline u8 ror2(u8 x) { return (u8)((x >> 2) | (x << 6)); }
static void wht(float *a) {
    for (int h = 1; h < 256; h <<= 1)
        for (int i = 0; i < 256; i += 2 * h)
            for (int j = i; j < i + h; j++) { float x = a[j], y = a[j + h]; a[j] = x + y; a[j + h] = x - y; }
}
static void transpose64(u64 *a) {
    int j = 32; u64 m = 0x00000000FFFFFFFFULL;
    for (; j != 0; j >>= 1, m ^= (m << j))
        for (int k = 0; k < 64; k = ((k | j) + 1) & ~j) {
            u64 t = ((a[k] >> j) ^ a[k | j]) & m;
            a[k] ^= t << j; a[k | j] ^= t;
        }
}
static void setup_tables(void) {
    build_etab();
    float e2[256];
    for (int g = 0; g < 256; g++) { Ehat[g] = (float)e_tab[g]; e2[g] = (float)(e_tab[g] * e_tab[g]); }
    wht(Ehat); memcpy(E2hat, e2, sizeof e2); wht(E2hat);
    static u8 HP[65536];
    for (u32 x = 0; x < 65536; x++) HP[x] = (u8)hprime(x >> 8, x & 255);
    for (u32 x = 0; x < 65536; x++) {
        u64 r0 = 0, r1 = 0;
        for (int r = 0; r < 128; r++) {
            int b = HP[x ^ emb_rest(r)];
            if (r < 64) r0 |= (u64)b << r; else r1 |= (u64)b << (r - 64);
        }
        ROW[x][0] = r0; ROW[x][1] = r1;
    }
}
static inline u32 L7of(int i, u32 DK7) { u32 R7 = CLa[i] ^ CRa[i]; return CLa[i] ^ feal_N(feal_M(R7) ^ DK7); }

/* ---------- 拡張ステップ（汎用） ---------- */
/* 標本ごとの既知量：推測する段の F 出力のバイト0 = Z0 ^ S0(X0^g0, Y1)、バイト3 = Z3 ^ S1(X3^g3, Y2)、
 * 次段の F 入力（鍵前）a = ((byte0^m1)<<8) | (m2^byte3)、既知パリティ O */
typedef struct { u8 X0[1024], Y1[1024], Z0[1024], X3[1024], Y2[1024], Z3[1024], m1[1024], m2[1024]; u64 O[W]; } sd_t;
typedef struct { u8 nrm, g0, g3; u16 k; float z; } hyp_t;

static int cmp_hyp(const void *a, const void *b) {
    float x = fabsf(((const hyp_t *)a)->z), y = fabsf(((const hyp_t *)b)->z);
    return x < y ? 1 : x > y ? -1 : 0;
}
/* 仮説はいったん大きな作業領域に集め、最後に |z| の大きい順に maxout 件を残す（生成順で打ち切らない） */
#define HBIG (1 << 18)
static int ext_step(const sd_t *sd, hyp_t *out_final, int maxout, double t1, double t2, hyp_t *best, long long *work) {
    static __thread hyp_t *big = NULL;
    if (!big) big = malloc(sizeof(hyp_t) * HBIG);
    hyp_t *out = big; int nbig = 0;
    int nout = 0; best->z = 0;
    long long wk = 0;
    for (int nr = 0; nr < 2; nr++) {
        const sd_t *s = &sd[nr];
        for (int gl = 0; gl < 64; gl++) {
            int g0l = gl & 3, g3l = gl >> 2;
            u8 u[1024];
            float T[256] = {0}, Ng[256] = {0};
            for (int i = 0; i < n; i++) {
                u8 b0 = s->Z0[i] ^ S0(s->X0[i] ^ (u8)g0l, s->Y1[i]);
                u8 b3 = s->Z3[i] ^ S1(s->X3[i] ^ (u8)g3l, s->Y2[i]);
                u[i] = (u8)part_idx(b0 ^ s->m1[i], s->m2[i] ^ b3);
                T[u[i]] += ((s->O[i >> 6] >> (i & 63)) & 1) ? -1.f : 1.f;
                Ng[u[i]] += 1;
            }
            wht(T); wht(Ng);
            for (int g = 0; g < 256; g++) { T[g] *= Ehat[g]; Ng[g] *= E2hat[g]; }
            wht(T); wht(Ng);
            int kl[256], nk = 0;
            for (int k = 0; k < 256; k++)
                if (fabsf(T[k]) >= t1 * sqrtf(Ng[k] * 256.f)) kl[nk++] = k;
            if (!nk) continue;
            /* 補完：g0 のビット2〜6、g3 のビット4〜6 */
            for (int ext = 0; ext < 256; ext++) {
                u8 g0 = (u8)(g0l | ((ext & 31) << 2)), g3 = (u8)(g3l | (((ext >> 5) & 7) << 4));
                u32 a[1024];
                for (int i = 0; i < n; i++) {
                    u8 b0 = s->Z0[i] ^ S0(s->X0[i] ^ g0, s->Y1[i]);
                    u8 b3 = s->Z3[i] ^ S1(s->X3[i] ^ g3, s->Y2[i]);
                    a[i] = ((u32)(b0 ^ s->m1[i]) << 8) | (u8)(s->m2[i] ^ b3);
                }
                int delta = partx(a[0]) ^ u[0];
                for (int q = 0; q < nk; q++) {
                    u32 base = emb_part(kl[q] ^ delta);
                    u64 V[128][W], blk[2][64];
                    for (int b = 0; b < W; b++) {
                        for (int j = 0; j < 64; j++) {
                            int i = b * 64 + j;
                            if (i < n) { const u64 *rw = ROW[a[i] ^ base]; blk[0][j] = rw[0]; blk[1][j] = rw[1]; }
                            else blk[0][j] = blk[1][j] = 0;
                        }
                        transpose64(blk[0]); transpose64(blk[1]);
                        for (int r = 0; r < 64; r++) { V[r][b] = blk[0][r]; V[64 + r][b] = blk[1][r]; }
                    }
                    for (int r = 0; r < 128; r++) {
                        int d = 0;
                        for (int w = 0; w < W; w++) d += __builtin_popcountll(V[r][w] ^ s->O[w]);
                        float z = (float)(n - 2 * d) / sqrtf((float)n);
                        wk++;
                        hyp_t h = { (u8)nr, g0, g3, (u16)(base | emb_rest(r)), z };
                        if (fabsf(z) > fabsf(best->z)) *best = h;
                        if (fabsf(z) >= t2 && nout < HBIG) out[nout++] = h;
                    }
                }
            }
        }
    }
    if (work) *work += wk;
    (void)nbig;
    if (nout > maxout) qsort(out, nout, sizeof(hyp_t), cmp_hyp);
    if (nout > maxout) nout = maxout;
    memcpy(out_final, out, sizeof(hyp_t) * nout);
    return nout;
}

static inline void setO(sd_t *s, int i, int bit) { if (bit) s->O[i >> 6] |= 1ULL << (i & 63); }

/* (a)：EK0 の中央2バイト e1,e2（e2 の最上位 0）、DK7 30bit */
static void build_a(u8 e1, u8 e2, u32 DK7, sd_t sd[2]) {
    for (int nr = 0; nr < 2; nr++) {
        sd_t *s = &sd[nr]; memset(s->O, 0, sizeof s->O);
        u8 k1 = e1 ^ (nr ? 0x80 : 0), k2 = e2 ^ (nr ? 0x80 : 0);
        for (int i = 0; i < n; i++) {
            u32 PL = PLa[i], PR = PRa[i], L1 = PL ^ PR, R7 = CLa[i] ^ CRa[i], L7 = L7of(i, DK7);
            u8 x1 = BY(L1, 0) ^ BY(L1, 1) ^ k1, x2 = BY(L1, 2) ^ BY(L1, 3) ^ k2;
            u8 Y1 = S1(x1, x2), Y2 = S0(x2, Y1);
            s->X0[i] = BY(L1, 0); s->Y1[i] = Y1; s->Z0[i] = BY(PL, 0);
            s->X3[i] = BY(L1, 3); s->Y2[i] = Y2; s->Z3[i] = BY(PL, 3);
            s->m1[i] = BY(PL, 1) ^ Y1; s->m2[i] = BY(PL, 2) ^ Y2;
            int o = ((s->m2[i] >> 2) & 1) ^ ((BY(PL, 3) >> 2) & 1) ^ (BY(L1, 3) & 1) ^ (Y2 & 1)
                  ^ par(0x04031004u & L1) ^ par(0x04031004u & L7) ^ par(0x00010504u & R7);
            setO(s, i, o);
        }
    }
}
/* (c)：DK7 30bit、DK6 の中央2バイト d1,d2（d2 の最上位 0） */
static void build_c(u32 DK7, u8 d1, u8 d2, sd_t sd[2]) {
    for (int nr = 0; nr < 2; nr++) {
        sd_t *s = &sd[nr]; memset(s->O, 0, sizeof s->O);
        u8 k1 = d1 ^ (nr ? 0x80 : 0), k2 = d2 ^ (nr ? 0x80 : 0);
        for (int i = 0; i < n; i++) {
            u32 PL = PLa[i], PR = PRa[i], R7 = CLa[i] ^ CRa[i], L7 = L7of(i, DK7);
            u8 w1 = BY(L7, 0) ^ BY(L7, 1) ^ k1, w2 = BY(L7, 2) ^ BY(L7, 3) ^ k2;
            u8 Y1 = S1(w1, w2), Y2 = S0(w2, Y1);
            s->X0[i] = BY(L7, 0); s->Y1[i] = Y1; s->Z0[i] = BY(R7, 0);
            s->X3[i] = BY(L7, 3); s->Y2[i] = Y2; s->Z3[i] = BY(R7, 3);
            s->m1[i] = BY(R7, 1) ^ Y1; s->m2[i] = BY(R7, 2) ^ Y2;
            int o = par(0x00010504u & PL) ^ par(0x04031004u & (PL ^ PR)) ^ par(0x04031004u & L7)
                  ^ ((s->m2[i] >> 2) & 1) ^ ((BY(R7, 3) >> 2) & 1) ^ (BY(L7, 3) & 1) ^ (Y2 & 1);
            setO(s, i, o);
        }
    }
}
/* (b)：EK0 全体、EK1 の中央2バイト k1,k2（k2 の最上位 0） */
static void build_b(u32 EK0, u8 k1n, u8 k2n, sd_t sd[2]) {
    for (int nr = 0; nr < 2; nr++) {
        sd_t *s = &sd[nr]; memset(s->O, 0, sizeof s->O);
        u8 k1 = k1n ^ (nr ? 0x80 : 0), k2 = k2n ^ (nr ? 0x80 : 0);
        for (int i = 0; i < n; i++) {
            u32 PL = PLa[i], PR = PRa[i], L1 = PL ^ PR, R1 = PL ^ feal_N(feal_M(L1) ^ EK0);
            u32 L8 = CLa[i] ^ CRa[i], R8 = CLa[i];
            u8 w1 = BY(R1, 0) ^ BY(R1, 1) ^ k1, w2 = BY(R1, 2) ^ BY(R1, 3) ^ k2;
            u8 Y1 = S1(w1, w2), Y2 = S0(w2, Y1);
            s->X0[i] = BY(R1, 0); s->Y1[i] = Y1; s->Z0[i] = BY(L1, 0);
            s->X3[i] = BY(R1, 3); s->Y2[i] = Y2; s->Z3[i] = BY(L1, 3);
            s->m1[i] = BY(L1, 1) ^ Y1; s->m2[i] = BY(L1, 2) ^ Y2;
            int o = ((s->m2[i] >> 2) & 1) ^ ((BY(L1, 3) >> 2) & 1) ^ (BY(R1, 3) & 1) ^ (Y2 & 1)
                  ^ par(0x04031004u & R1) ^ par(0x04031004u & L8) ^ par(0x00010504u & R8);
            setO(s, i, o);
        }
    }
}
/* (t)：第2〜4段の 3 段近似。EK0 全体、EK1 中央、DK7、DK6 全体、DK5 中央 */
static double test3(u32 EK0, u8 ek1_1, u8 ek1_2, u32 DK7, u32 DK6, u8 dk5_1, u8 dk5_2) {
    int s = 0;
    for (int i = 0; i < n; i++) {
        u32 PL = PLa[i], PR = PRa[i], L1 = PL ^ PR, R1 = PL ^ feal_N(feal_M(L1) ^ EK0);
        u32 R7 = CLa[i] ^ CRa[i], L7 = L7of(i, DK7), L6 = R7 ^ feal_N(feal_M(L7) ^ DK6);
        u32 a = feal_N((feal_M(R1) & 0x00FFFF00u) ^ PK(0, ek1_1, ek1_2, 0));
        u32 c = feal_N((feal_M(L6) & 0x00FFFF00u) ^ PK(0, dk5_1, dk5_2, 0));
        u8 y1 = BY(a, 1), y2 = BY(a, 2), v1 = BY(c, 1), v2 = BY(c, 2);
        int p = par(0x00010504u & R1) ^ par(0x04031004u & L1) ^ (BY(R1, 0) & 1) ^ ((y1 >> 1) & 1) ^ ((y2 >> 4) & 1)
              ^ (BY(R1, 3) & 1) ^ (y2 & 1)
              ^ par(0x04031104u & L7) ^ (BY(L6, 0) & 1) ^ ((v1 >> 1) & 1) ^ ((v2 >> 4) & 1) ^ (BY(L6, 3) & 1);
        s += p ? -1 : 1;
    }
    return s / sqrt(n);
}

/* (t) をビット列で：平文側の寄与（EK0 全体, EK1 中央）と暗号文側の寄与（DK7, DK6 全体, DK5 中央）に分ける */
static void p3vec(u32 EK0, u8 k1, u8 k2, u64 *v) {
    memset(v, 0, sizeof(u64) * W);
    for (int i = 0; i < n; i++) {
        u32 PL = PLa[i], PR = PRa[i], L1 = PL ^ PR, R1 = PL ^ feal_N(feal_M(L1) ^ EK0);
        u32 a = feal_N((feal_M(R1) & 0x00FFFF00u) ^ PK(0, k1, k2, 0));
        u8 y1 = BY(a, 1), y2 = BY(a, 2);
        int p = par(0x00010504u & R1) ^ par(0x04031004u & L1) ^ (BY(R1, 0) & 1) ^ ((y1 >> 1) & 1) ^ ((y2 >> 4) & 1)
              ^ (BY(R1, 3) & 1) ^ (y2 & 1);
        if (p) v[i >> 6] |= 1ULL << (i & 63);
    }
}
static void c3vec(u32 DK7, u32 DK6, u8 k1, u8 k2, u64 *v) {
    memset(v, 0, sizeof(u64) * W);
    for (int i = 0; i < n; i++) {
        u32 R7 = CLa[i] ^ CRa[i], L7 = L7of(i, DK7), L6 = R7 ^ feal_N(feal_M(L7) ^ DK6);
        u32 c = feal_N((feal_M(L6) & 0x00FFFF00u) ^ PK(0, k1, k2, 0));
        u8 v1 = BY(c, 1), v2 = BY(c, 2);
        int p = par(0x04031104u & L7) ^ (BY(L6, 0) & 1) ^ ((v1 >> 1) & 1) ^ ((v2 >> 4) & 1) ^ (BY(L6, 3) & 1);
        if (p) v[i >> 6] |= 1ULL << (i & 63);
    }
}

/* ---------- (e) 鍵スケジュールの逆算 ---------- */
static inline u32 fKinv(u32 a, u32 o) {
    u8 a0 = BY(a, 0), a3 = BY(a, 3), t1 = BY(a, 0) ^ BY(a, 1), t2 = BY(a, 2) ^ BY(a, 3);
    u8 o0 = BY(o, 0), o1 = BY(o, 1), o2 = BY(o, 2), o3 = BY(o, 3);
    u8 b0 = (u8)(ror2(o1) - t1 - 1) ^ t2, b1 = (u8)(ror2(o2) - t2) ^ o1;
    u8 b2 = (u8)(ror2(o0) - a0) ^ o1, b3 = (u8)(ror2(o3) - a3 - 1) ^ o2;
    return PK(b0, b1, b2, b3);
}
static u64 pt_chk[2], ct_chk[2];   /* 順写像用：読み込んだ平文 X → 暗号文 Y */
static u64 rpt_chk[2], rct_chk[2]; /* 逆写像用：Y → X */
static int BOTHMAP = 1;            /* (e) で両方の写像（Y=E(X) と X=E(Y)）を試す */
static int amb_lo = 0, amb_hi = 128, q0_lo = 0, q0_hi = 256;
/* (e) 入力：EK0（バイト1,2 正確、バイト0,3 はビット0〜6）、EK1 中央（(a) の正規化済み k）、
 *          DK7（30bit）、DK6（バイト1,2 正確、バイト0,3 はビット0〜6）
 * 曖昧さ 7bit：EK0 の最上位 m0,m3、DK7 の最上位 d0,d3、DK6 の最上位 f0,f3、EK1 中央の正規化 nr1
 * EK1 のバイト0（q0）は全数、バイト3（q3）は 8bit の整合条件でほぼ一意。残り X2, Q2 の 2^16 を全数 */
static int recover_key(u32 EK0f, u16 ek1mid, u32 DK7f, u32 DK6f, u64 *keys, int maxk) {
    int nk = 0;
    #pragma omp parallel for schedule(dynamic, 1) collapse(2)
    for (int amb = amb_lo; amb < amb_hi; amb++)
        for (int q0i = q0_lo; q0i < q0_hi; q0i++) {
            int m0 = amb & 1, m3 = (amb >> 1) & 1, d0 = (amb >> 2) & 1, d3 = (amb >> 3) & 1;
            int f0 = (amb >> 4) & 1, f3 = (amb >> 5) & 1, nr1 = (amb >> 6) & 1;
            u32 EK0 = EK0f ^ PK(m0 << 7, 0, 0, m3 << 7);
            u8 E11 = (u8)(ek1mid >> 8) ^ (nr1 ? 0x80 : 0) ^ (m0 << 1), E12 = (u8)ek1mid ^ (nr1 ? 0x80 : 0) ^ (m3 << 1);
            u32 DK7 = DK7f ^ PK(d0 << 7, 0, 0, d3 << 7);
            u32 DK6 = DK6f ^ feal_M(PK(d0 << 1, 0, 0, d3 << 1)) ^ PK(f0 << 7, 0, 0, f3 << 7);
            u8 c00 = BY(EK0, 0), c03 = BY(EK0, 3), k0 = BY(DK6, 0), k3 = BY(DK6, 3), D0 = BY(DK7, 0), D3 = BY(DK7, 3);
            u8 DK61 = BY(DK6, 1), DK62 = BY(DK6, 2), DK71 = BY(DK7, 1), DK72 = BY(DK7, 2);
            u8 q0 = (u8)q0i, B60 = q0 ^ c00, B80 = k0 ^ D0, B83 = k3 ^ D3;
            u8 cA = (u8)(ror2(k0) - q0) ^ DK71 ^ D0, cU = (u8)(ror2(B80) - B60);
            for (int q3i = 0; q3i < 256; q3i++) {
                u8 q3 = (u8)q3i, B63 = q3 ^ c03;
                u8 o2 = (u8)(ror2(B83) - B63 - 1) ^ k3 ^ q3;
                if ((u8)(o2 ^ B63 ^ DK72 ^ D3) != (u8)(ror2(k3) - q3 - 1)) continue;
                u8 r2 = ror2(o2);
                for (int X2i = 0; X2i < 256; X2i++) {
                    u8 X2 = (u8)X2i, o1 = cA ^ X2, zz = X2 ^ B63 ^ k0 ^ q0;
                    u8 X1 = (u8)(ror2(o1) - zz - 1) ^ B60, t2 = X2 ^ B63;
                    u8 cq = (u8)(r2 - t2) ^ o1;
                    for (int Q2i = 0; Q2i < 256; Q2i++) {
                        u8 Q2 = (u8)Q2i, U2 = cU ^ o1 ^ Q2;
                        u8 U1 = (u8)(ror2(U2) - (Q2 ^ q3)) ^ X1 ^ DK62 ^ U2 ^ k3;
                        u8 Q1 = U1 ^ cq;
                        if (S1(q0 ^ Q1, Q2 ^ q3 ^ B60 ^ DK61 ^ k0 ^ U1) != U1) continue;
                        u8 V1 = o1 ^ U1, V2 = o2 ^ U2;
                        u32 B5 = PK(q0, Q1, Q2, q3), B6 = PK(B60, X1, X2, B63), B7 = PK(k0, U1, U2, k3);
                        u32 B8 = B7 ^ PK(D0, V1, V2, D3);
                        u32 B4 = PK(DK61 ^ k0 ^ U1, DK62 ^ U2 ^ k3, DK71 ^ D0 ^ V1, DK72 ^ V2 ^ D3);
                        if (fK(B5, B6 ^ B4) != B7 || fK(B6, B7 ^ B5) != B8) continue;
                        u8 P1 = Q1 ^ X1, P2 = Q2 ^ X2;
                        u32 B3 = B5 ^ fKinv(B4, B6), B2 = B4 ^ fKinv(B3, B5), B1 = B3 ^ fKinv(B2, B4);
                        u32 Kexp = PK(BY(EK0, 1) ^ c00 ^ P1, BY(EK0, 2) ^ P2 ^ c03, E11 ^ q0 ^ Q1, E12 ^ Q2 ^ q3);
                        if (B1 != Kexp) continue;
                        u32 B0 = B2 ^ fKinv(B1, B3), A0 = B1 ^ fKinv(B0, B2);
                        if (fK(A0, B0) != B1) continue;
                        u64 key = ((u64)A0 << 32) | B0;
                        sched_t s; key_schedule(key, &s);
                        if (feal8_enc(&s, pt_chk[0], NULL) != ct_chk[0] || feal8_enc(&s, pt_chk[1], NULL) != ct_chk[1]) continue;
                        #pragma omp critical
                        { if (nk < maxk) keys[nk++] = key; }
                    }
                }
            }
        }
    return nk;
}

/* 逆向き（平文と暗号文を入れ替えて段階1〜3を回した場合）の (e)。
 * 入れ替えた暗号の実効鍵は、本来の鍵で EK*0 = M(D)^(0,K7,0)、EK*1 = M(Kcd)^(0,K6,0)、
 * DK*7 = EK0 = M(c0)^(0,K0,0)、DK*6 = EK1 = M(K89)^(0,K1,0)（D = Kcd^Kef、c0 = K89^Kab）。
 * 未知：Kcd（U0..U3）、K89 中央（Q1,Q2）、Kab 中央（X1,X2）、D 中央（V1,V2）。U3 は 8bit 条件でほぼ一意、U0,X2,Q2 を全数 */
static int REV = 0;
static int recover_key_rev(u32 EK0f, u16 ek1mid, u32 DK7f, u32 DK6f, u64 *keys, int maxk) {
    int nk = 0;
    #pragma omp parallel for schedule(dynamic, 1) collapse(2)
    for (int amb = amb_lo; amb < amb_hi; amb++)
        for (int U0i = 0; U0i < 256; U0i++) {
            int m0 = amb & 1, m3 = (amb >> 1) & 1, d0 = (amb >> 2) & 1, d3 = (amb >> 3) & 1;
            int f0 = (amb >> 4) & 1, f3 = (amb >> 5) & 1, nr1 = (amb >> 6) & 1;
            u32 E0 = EK0f ^ PK(m0 << 7, 0, 0, m3 << 7);
            u8 E11 = (u8)(ek1mid >> 8) ^ (nr1 ? 0x80 : 0) ^ (m0 << 1), E12 = (u8)ek1mid ^ (nr1 ? 0x80 : 0) ^ (m3 << 1);
            u32 R0 = DK7f ^ PK(d0 << 7, 0, 0, d3 << 7);                                     /* 本来の EK0 */
            u32 R1 = DK6f ^ feal_M(PK(d0 << 1, 0, 0, d3 << 1)) ^ PK(f0 << 7, 0, 0, f3 << 7); /* 本来の EK1 */
            u8 D0 = BY(E0, 0), D3 = BY(E0, 3), E01 = BY(E0, 1), E02 = BY(E0, 2);
            u8 c00 = BY(R0, 0), c03 = BY(R0, 3), q0 = BY(R1, 0), q3 = BY(R1, 3);
            u8 B60 = q0 ^ c00, B63 = q3 ^ c03, U0 = (u8)U0i;
            for (int U3i = 0; U3i < 256; U3i++) {
                u8 U3 = (u8)U3i;
                u8 o2 = (u8)(ror2(U3) - q3 - 1) ^ B63 ^ E02 ^ D3;
                if (S1(B63, o2 ^ U3 ^ q3) != (u8)(U3 ^ D3)) continue;
                u8 r2 = ror2(o2);
                for (int X2i = 0; X2i < 256; X2i++) {
                    u8 X2 = (u8)X2i;
                    u8 o1 = (u8)(ror2(U0) - q0) ^ X2 ^ E01 ^ D0;
                    u8 X1 = (u8)(ror2(o1) - (X2 ^ B63 ^ U0 ^ q0) - 1) ^ B60;
                    u8 t2 = X2 ^ B63, cU2 = (u8)(ror2(U0 ^ D0) - B60) ^ o1, cq = (u8)(r2 - t2) ^ o1;
                    for (int Q2i = 0; Q2i < 256; Q2i++) {
                        u8 Q2 = (u8)Q2i, U2 = cU2 ^ Q2;
                        u8 U1 = (u8)(ror2(U2) - (Q2 ^ q3)) ^ X1 ^ E12 ^ U2 ^ U3;
                        u8 Q1 = cq ^ U1;
                        if (S1(q0 ^ Q1, Q2 ^ q3 ^ B60 ^ E11 ^ U0 ^ U1) != U1) continue;
                        u8 V1 = o1 ^ U1, V2 = o2 ^ U2;
                        u32 B5 = PK(q0, Q1, Q2, q3), B6 = PK(B60, X1, X2, B63), B7 = PK(U0, U1, U2, U3);
                        u32 B8 = B7 ^ PK(D0, V1, V2, D3);
                        u32 B4 = PK(E11 ^ U0 ^ U1, E12 ^ U2 ^ U3, E01 ^ D0 ^ V1, E02 ^ V2 ^ D3);
                        if (fK(B5, B6 ^ B4) != B7 || fK(B6, B7 ^ B5) != B8) continue;
                        u8 P1 = Q1 ^ X1, P2 = Q2 ^ X2;
                        u32 B3 = B5 ^ fKinv(B4, B6), B2 = B4 ^ fKinv(B3, B5), B1 = B3 ^ fKinv(B2, B4);
                        u32 Kexp = PK(BY(R0, 1) ^ c00 ^ P1, BY(R0, 2) ^ P2 ^ c03, BY(R1, 1) ^ q0 ^ Q1, BY(R1, 2) ^ Q2 ^ q3);
                        if (B1 != Kexp) continue;
                        u32 B0 = B2 ^ fKinv(B1, B3), A0 = B1 ^ fKinv(B0, B2);
                        if (fK(A0, B0) != B1) continue;
                        u64 key = ((u64)A0 << 32) | B0;
                        sched_t s; key_schedule(key, &s);
                        if (feal8_enc(&s, rpt_chk[0], NULL) != rct_chk[0] || feal8_enc(&s, rpt_chk[1], NULL) != rct_chk[1]) continue;
                        #pragma omp critical
                        { if (nk < maxk) keys[nk++] = key; }
                    }
                }
            }
        }
    return nk;
}

static void load_pairs(const char *pf, const char *cf, int self, u64 key) {
    FILE *fp = fopen(pf, "rb"); u8 pb[8192]; if (fread(pb, 1, 8192, fp) != 8192) exit(1); fclose(fp);
    u8 cb[8192];
    if (!self) { fp = fopen(cf, "rb"); if (fread(cb, 1, 8192, fp) != 8192) exit(1); fclose(fp); }
    sched_t s; if (self) key_schedule(key, &s);
    u64 seen[1024]; n = 0;
    for (int b = 0; b < 1024; b++) {
        u64 p = 0, c = 0;
        for (int j = 0; j < 8; j++) { p = (p << 8) | pb[8 * b + j]; if (!self) c = (c << 8) | cb[8 * b + j]; }
        int dup = 0;
        for (int j = 0; j < n; j++) if (seen[j] == p) { dup = 1; break; }
        if (dup) continue;
        if (self) c = feal8_enc(&s, p, NULL);
        seen[n] = p; PLa[n] = p >> 32; PRa[n] = (u32)p; CLa[n] = c >> 32; CRa[n] = (u32)c; n++;
    }
    if (REV && self) { /* 自己試験：作った (P,C) を入れ替える */
        for (int i = 0; i < n; i++) { u32 t = PLa[i]; PLa[i] = CLa[i]; CLa[i] = t; t = PRa[i]; PRa[i] = CRa[i]; CRa[i] = t; }
    }
    /* 照合：順写像は X → Y、逆写像は Y → X（X = 読み込んだ平文側、Y = 暗号文側） */
    pt_chk[0] = ((u64)PLa[0] << 32) | PRa[0]; ct_chk[0] = ((u64)CLa[0] << 32) | CRa[0];
    pt_chk[1] = ((u64)PLa[1] << 32) | PRa[1]; ct_chk[1] = ((u64)CLa[1] << 32) | CRa[1];
    rpt_chk[0] = ct_chk[0]; rct_chk[0] = pt_chk[0]; rpt_chk[1] = ct_chk[1]; rct_chk[1] = pt_chk[1];
}

/* 段階2の候補 (t,e,κ,e8,r) → (EK0 中央, DK7 30bit, DK6 中央) */
typedef struct { u8 e1, e2; u32 DK7; u8 d1, d2; float z; } cand_t;
static cand_t conv_cand(u32 t, int e, int kap, int e8, int r, float z) {
    cand_t c; c.e1 = (u8)e; c.e2 = (u8)(e >> 8) & 0x7F; c.z = z;
    u32 DKs = dk7_22(t); c.DK7 = DKs | dk7_ext(e8);
    u32 L7s = L7of(0, DKs), L7f = L7of(0, c.DK7);
    int u0 = part_idx(BY(L7s, 0) ^ BY(L7s, 1), BY(L7s, 2) ^ BY(L7s, 3));
    u32 a0 = ((u32)(BY(L7f, 0) ^ BY(L7f, 1)) << 8) | (BY(L7f, 2) ^ BY(L7f, 3));
    u32 k16 = emb_part(kap ^ partx(a0) ^ u0) | emb_rest(r);
    c.d1 = (u8)(k16 >> 8); c.d2 = (u8)k16;
    return c;
}

typedef struct { int ok; u64 key; double ta, tc, tt, te; float za, zc, z3; int npairs, ntried; } res_t;
typedef struct { int i, j; double sc, z; } pair_t;
static int cmp_pair(const void *a, const void *b) { double x = ((const pair_t *)a)->sc, y = ((const pair_t *)b)->sc; return x < y ? 1 : x > y ? -1 : 0; }

static double T1 = 2.0, T2 = 5.0, T3 = 8.0; static int MAXE = 64; static int NO_E = 0;
static res_t run_cand(const cand_t *c, int verbose) {
    res_t R = {0};
    static __thread sd_t sd[2];
    static __thread hyp_t ha[4096], hc[4096];
    hyp_t best;
    double t0 = now();
    build_a(c->e1, c->e2, c->DK7, sd);
    int na = ext_step(sd, ha, 4096, T1, T2, &best, NULL);
    R.za = best.z; R.ta = now() - t0;
    if (verbose) printf("  (a) %d hyps >= %.1f, best z=%.2f nrm=%d g0=%02x g3=%02x k=%04x  (%.2fs)\n", na, T2, best.z, best.nrm, best.g0, best.g3, best.k, R.ta);
    if (!na) return R;
    t0 = now();
    build_c(c->DK7, c->d1, c->d2, sd);
    int nc = ext_step(sd, hc, 4096, T1, T2, &best, NULL);
    R.zc = best.z; R.tc = now() - t0;
    if (verbose) printf("  (c) %d hyps >= %.1f, best z=%.2f nrm=%d g0=%02x g3=%02x k=%04x  (%.2fs)\n", nc, T2, best.z, best.nrm, best.g0, best.g3, best.k, R.tc);
    if (!nc) return R;
    t0 = now();
    const int PPMAX = 1 << 16;
    pair_t *pp = malloc(sizeof(pair_t) * PPMAX); int np = 0; double bz = 0;
    u64 (*PV)[W] = malloc(sizeof(*PV) * na), (*CV)[W] = malloc(sizeof(*CV) * nc);
    for (int i = 0; i < na; i++) {
        u8 n0 = ha[i].nrm ? 0x80 : 0;
        p3vec(PK(ha[i].g0, c->e1 ^ n0, c->e2 ^ n0, ha[i].g3), (u8)(ha[i].k >> 8), (u8)ha[i].k, PV[i]);
    }
    for (int j = 0; j < nc; j++) {
        u8 n6 = hc[j].nrm ? 0x80 : 0;
        c3vec(c->DK7, PK(hc[j].g0, c->d1 ^ n6, c->d2 ^ n6, hc[j].g3), (u8)(hc[j].k >> 8), (u8)hc[j].k, CV[j]);
    }
    for (int i = 0; i < na; i++)
        for (int j = 0; j < nc; j++) {
            int d = 0;
            for (int w = 0; w < W; w++) d += __builtin_popcountll(PV[i][w] ^ CV[j][w]);
            double z = (n - 2.0 * d) / sqrt(n);
            if (fabs(z) > fabs(bz)) bz = z;
            if (fabs(z) >= T3 && np < PPMAX) { pp[np].i = i; pp[np].j = j; pp[np].z = z; pp[np].sc = z * z + (double)ha[i].z * ha[i].z + (double)hc[j].z * hc[j].z; np++; }
        }
    free(PV); free(CV);
    qsort(pp, np, sizeof(pair_t), cmp_pair);
    R.z3 = (float)bz; R.tt = now() - t0; R.npairs = np;
    if (verbose) printf("  (t) best 3-round z=%.2f; %d of %d x %d pairs pass %.1f (%.2fs)\n", bz, np, na, nc, T3, R.tt);
    t0 = now();
    for (int q = 0; q < np && q < MAXE && !R.ok && !NO_E; q++) {
        const hyp_t *A = &ha[pp[q].i], *C = &hc[pp[q].j];
        u8 n0 = A->nrm ? 0x80 : 0, n6 = C->nrm ? 0x80 : 0;
        u32 EK0 = PK(A->g0, c->e1 ^ n0, c->e2 ^ n0, A->g3);
        u32 DK6 = PK(C->g0, c->d1 ^ n6, c->d2 ^ n6, C->g3);
        u64 keys[8];
        int nk = REV ? recover_key_rev(EK0, A->k, c->DK7, DK6, keys, 8) : recover_key(EK0, A->k, c->DK7, DK6, keys, 8);
        int alt = 0;
        if (!nk && BOTHMAP) { nk = REV ? recover_key(EK0, A->k, c->DK7, DK6, keys, 8) : recover_key_rev(EK0, A->k, c->DK7, DK6, keys, 8); alt = nk > 0; }
        if (alt && verbose) printf("  (e) key found with the OTHER mapping (%s)\n", REV ? "data satisfies Y = E(X)" : "data satisfies X = E(Y), i.e. Y = D(X)");
        R.ntried++;
        if (nk) { R.ok = 1; R.key = keys[0];
            if (verbose) printf("  (e) key found at pair %d (EK0=%08x EK1mid=%04x DK6=%08x, z3=%.2f)\n", q, EK0, A->k, DK6, pp[q].z); }
    }
    R.te = now() - t0;
    if (verbose) printf("  (e) tried %d pair(s), %s (%.1fs)\n", R.ntried, R.ok ? "key recovered" : "no key", R.te);
    free(pp);
    return R;
}

int main(int argc, char **argv) {
    int self = !strcmp(argv[1], "self");
    if (getenv("REV")) REV = 1;
    if (getenv("ONEMAP")) BOTHMAP = 0;
    if (getenv("T1")) T1 = atof(getenv("T1"));
    if (getenv("T2")) T2 = atof(getenv("T2"));
    if (getenv("T3")) T3 = atof(getenv("T3"));
    if (getenv("MAXE")) MAXE = atoi(getenv("MAXE"));
    setup_tables();
    const char *cf = NULL;
    u64 key = 0;
    if (self) {
        seed(strtoull(argv[3], 0, 0)); key = rnd();
        load_pairs(argv[2], NULL, 1, key);
        if (argc > 4) cf = argv[4];
    } else { load_pairs(argv[2], argv[3], 0, 0); cf = argv[4]; }
    printf("N=%d T1=%.2f T2=%.2f T3=%.2f threads=%d%s\n", n, T1, T2, T3, omp_get_max_threads(), REV ? " REVERSE" : "");
    if (self) {
        sched_t s; key_schedule(key, &s);
        u32 K[16]; for (int i = 0; i < 16; i++) K[i] = s.k[i];
        u32 K89 = (K[8] << 16) | K[9], Kab = (K[10] << 16) | K[11], Kcd = (K[12] << 16) | K[13], Kef = (K[14] << 16) | K[15];
        if (REV) { /* 入れ替えた暗号の副鍵：段鍵は逆順、白色化鍵は K89<->Kcd, Kab<->Kef */
            u32 Kr[8]; for (int i = 0; i < 8; i++) Kr[i] = K[7 - i];
            for (int i = 0; i < 8; i++) K[i] = Kr[i];
            u32 a = K89, b = Kab; K89 = Kcd; Kab = Kef; Kcd = a; Kef = b;
        }
        u32 EK[8], DK[8];
        for (int i = 0; i < 8; i++) {
            EK[i] = feal_M((i & 1) ? K89 : (K89 ^ Kab)) ^ (K[i] << 8);
            DK[i] = feal_M((i & 1) ? (Kcd ^ Kef) : Kcd) ^ (K[i] << 8);
        }
        printf("key %016llx\ntrue EK0=%08x EK1=%08x EK2=%08x  DK7=%08x DK6=%08x DK5=%08x\n", (unsigned long long)key,
               EK[0], EK[1], EK[2], DK[7], DK[6], DK[5]);
        /* 段階2の出力に相当する真の候補 */
        cand_t c; u8 e1 = BY(EK[0], 1), e2 = BY(EK[0], 2);
        if (e2 & 0x80) { e1 ^= 0x80; e2 ^= 0x80; }
        c.e1 = e1; c.e2 = e2; c.DK7 = DK[7] & 0x7FFFFF7Fu;
        u32 d0 = (DK[7] >> 31) & 1, d3 = (DK[7] >> 7) & 1;
        u32 DK6f = DK[6] ^ feal_M(PK(d0 << 1, 0, 0, d3 << 1));
        u8 d1 = BY(DK6f, 1), d2 = BY(DK6f, 2);
        if (d2 & 0x80) { d1 ^= 0x80; d2 ^= 0x80; }
        c.d1 = d1; c.d2 = d2; c.z = 0;
        printf("true candidate: e1=%02x e2=%02x DK7=%08x d1=%02x d2=%02x\n", c.e1, c.e2, c.DK7, c.d1, c.d2);
        /* 期待される (a)(c) の答え */
        u32 m0 = (EK[0] >> 31) & 1, m3 = (EK[0] >> 7) & 1;
        u32 EK1f = EK[1] ^ feal_M(PK(m0 << 1, 0, 0, m3 << 1));
        u32 f0 = (DK6f >> 31) & 1, f3 = (DK6f >> 7) & 1;
        u32 L6off = PK(f0 << 1, 0, 0, f3 << 1);
        u32 DK5f = DK[5] ^ feal_M(L6off);
        printf("expect (a): g0=%02x g3=%02x EK1 mid=%02x%02x ; (c): g0=%02x g3=%02x DK5 mid=%02x%02x (top bits of mid only via XOR)\n",
               BY(EK[0], 0) & 0x7F, BY(EK[0], 3) & 0x7F, BY(EK1f, 1), BY(EK1f, 2),
               BY(DK6f, 0) & 0x7F, BY(DK6f, 3) & 0x7F, BY(DK5f, 1), BY(DK5f, 2));
        u8 k1n = BY(EK1f, 1), k2n = BY(EK1f, 2); int nr1 = 0;
        if (k2n & 0x80) { k1n ^= 0x80; k2n ^= 0x80; nr1 = 1; }
        if (getenv("EFAST")) {
            int amb = m0 | (m3 << 1) | (d0 << 2) | (d3 << 3) | (f0 << 4) | (f3 << 5) | (nr1 << 6);
            int q0 = BY(EK[1], 0);
            amb_lo = amb; amb_hi = amb + 1; if (!REV) { q0_lo = q0; q0_hi = q0 + 1; }
            printf("EFAST: restricting (e) to amb=%02x q0=%02x\n", amb, q0);
        }
        if (getenv("ETEST")) {
            u32 EK0f = EK[0] & 0x7FFFFF7Fu, DK6x = DK6f & 0x7FFFFF7Fu;
            u64 keys[8]; double t0 = now();
            int nk = recover_key(EK0f, (u16)((k1n << 8) | k2n), c.DK7, DK6x, keys, 8);
            printf("ETEST: (e) on expected inputs -> %d key(s) %016llx (%s) %.1fs\n", nk, nk ? (unsigned long long)keys[0] : 0ULL,
                   nk && keys[0] == key ? "matches" : "NO", now() - t0);
        }
        res_t R = run_cand(&c, 1);
        printf("TRUE candidate result: %s key=%016llx (%s)\n", R.ok ? "RECOVERED" : "failed", (unsigned long long)R.key,
               R.ok && R.key == key ? "matches" : "no match");
    }
    if (cf) {
        FILE *fp = fopen(cf, "r");
        cand_t *cs = malloc(sizeof(cand_t) * 4000000); int nc = 0;
        u32 t; int e, kap, e8, r; float z; char line[256];
        while (nc < 4000000 && fgets(line, sizeof line, fp)) {
            if (sscanf(line, "%x %x %x %x %x %f", &t, &e, &kap, &e8, &r, &z) != 6) continue;
            cs[nc++] = conv_cand(t, e, kap, e8, r, z);
        }
        fclose(fp);
        int lim = getenv("LIMIT") ? atoi(getenv("LIMIT")) : nc;
        if (lim < nc) nc = lim;
        printf("processing %d candidates from %s\n", nc, cf);
        double W0 = now();
        int pa = 0, pc = 0, pt = 0; double sa = 0;
        NO_E = 1;
        res_t *RR = calloc(nc, sizeof(res_t));
        #pragma omp parallel for schedule(dynamic, 1) reduction(+:pa,pc,pt,sa)
        for (int i = 0; i < nc; i++) {
            res_t R = run_cand(&cs[i], 0);
            RR[i] = R;
            sa += R.ta + R.tc + R.tt;
            if (R.za != 0 && fabsf(R.za) >= T2) pa++;
            if (R.zc != 0 && fabsf(R.zc) >= T2) pc++;
            if (R.npairs > 0) { pt++;
                #pragma omp critical
                { printf("  cand %d (z2=%.2f) passed 3-round test: z3=%.2f, pairs %d\n", i, cs[i].z, R.z3, R.npairs); fflush(stdout); } }
        }
        double W1 = now();
        printf("filter: %d cands, pass(a)=%d (%.1f%%), pass(c)=%d (%.1f%% of a), pass(t)=%d; wall %.1fs (%.3f s/cand wall), thread time %.2f s/cand\n",
               nc, pa, 100.0 * pa / nc, pc, pa ? 100.0 * pc / pa : 0, pt, W1 - W0, (W1 - W0) / nc, sa / nc);
        if (getenv("RECOVER")) { NO_E = 0;
            for (int i = 0; i < nc; i++) if (RR[i].npairs > 0) {
                printf("recovering from cand %d:\n", i); res_t R = run_cand(&cs[i], 1);
                if (R.ok) { printf("KEY %016llx\n", (unsigned long long)R.key); break; } } }
    }
    return 0;
}
