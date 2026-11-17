/* 崎小山ら（IEICE 2016）4.3 の第1段階（6段近似 8本の MLC）を FEAL-8 で再現し、
 * 正解 (mK1 15bit, mK8 22bit) の 2^37 中の順位を見積もる。
 *
 *   ./mlc real plaintext.bin ntrials logS seed
 *   ./mlc rand N             ntrials logS seed
 *
 * 平文側は EK0 の 2^15 通りを全部、暗号文側は真の DK7 と無作為な DK7 を 2^logS 個。
 * 真の DK7 の行は厳密に数え、無作為な行は 2^22 / 2^logS 倍して全空間の順位に換算する。
 * 統計量：
 *   s_mlc  = Σ ε_l |dev_l|（論文の式(10)、ε = 97/65）
 *   s_fam0 = 第1族（ΓY2=04031104）4本だけの Σ ε|dev|
 *   s_one  = 近似1本（Approximation 1）の |dev|
 *   s_quad = Σ (ε_l dev_l)^2
 * あわせて同じ鍵で方式C の参照統計量（真の中間値、attackC.c と同じ定義）を出す。 */
#include "feal8.h"
#include "approx8.h"
#include <omp.h>

#define W 16
#define NS 4
static int n;
static u32 PL[1024], PR[1024], CL[1024], CR[1024];
static u64 (*P)[2][W];
static u64 (*C)[8][W];
static double etab[2][256];

static inline u32 dk7_22(u32 t) {
    return PK((u8)((t >> 16) & 3), (u8)t, (u8)(t >> 8), (u8)((t >> 18) & 15));
}
static void build_etab2(void) {
    for (int f = 0; f < 2; f++) {
        double acc[256] = {0}; int cnt[256] = {0};
        for (int w1 = 0; w1 < 256; w1++)
            for (int w2 = 0; w2 < 256; w2++) {
                int g = part_idx((u8)w1, (u8)w2);
                acc[g] += (par(A8_Y6[f] & feal_N(PK(0, (u8)w1, (u8)w2, 0))) ^ (w1 & 1) ^ (w2 & 1)) ? -1 : 1;
                cnt[g]++;
            }
        double cap = 0;
        for (int g = 0; g < 256; g++) { etab[f][g] = acc[g] / cnt[g]; cap += etab[f][g] * etab[f][g] * cnt[g] / 65536.0; }
        fprintf(stderr, "partition capacity fam%d = %.4f\n", f, cap);
    }
}

static inline void scores(const u64 *p0, const u64 *p1, u64 c[8][W], long long s[NS]) {
    int d[8];
    for (int l = 0; l < 8; l++) {
        const u64 *p = (l < 4) ? p0 : p1;
        int pc = 0;
        for (int w = 0; w < W; w++) pc += __builtin_popcountll(p[w] ^ c[l][w]);
        d[l] = n - 2 * pc;
        if (d[l] < 0) d[l] = -d[l];
    }
    s[0] = 97LL * (d[0] + d[2] + d[4] + d[6]) + 65LL * (d[1] + d[3] + d[5] + d[7]);
    s[1] = 97LL * (d[0] + d[2]) + 65LL * (d[1] + d[3]);
    s[2] = d[0];
    long long q = 0;
    for (int l = 0; l < 8; l++) { long long v = (long long)A8_EPS[l & 3] * d[l]; q += v * v; }
    s[3] = q;
}

int main(int argc, char **argv) {
    if (argc < 6) return 1;
    int real = !strcmp(argv[1], "real");
    int ntr = atoi(argv[3]), logS = atoi(argv[4]);
    seed(strtoull(argv[5], 0, 0));
    int S = 1 << logS;
    u64 *pre = NULL; int npre = 0;
    if (real) pre = load_unique(argv[2], &npre);
    int Nr = real ? npre : atoi(argv[2]);
    build_etab2();
    P = aligned_alloc(64, sizeof(*P) * 32768);
    C = aligned_alloc(64, sizeof(*C) * (S + 1));
    printf("# mode=%s N=%d trials=%d S=2^%d threads=%d\n", argv[1], Nr, ntr, logS, omp_get_max_threads());
    printf("# trial key N zA0 zA1 zB0 zB1 pre0 pre1 | true s_mlc | rank_mlc rank_fam0 rank_one rank_quad (log2) | consistent\n");
    for (int tr = 0; tr < ntr; tr++) {
        u64 key = rnd();
        sched_t s; key_schedule(key, &s);
        n = Nr;
        if (real) for (int i = 0; i < n; i++) { PL[i] = pre[i] >> 32; PR[i] = (u32)pre[i]; }
        else {
            for (int i = 0; i < n; i++) {
                u64 p; int dup;
                do { p = rnd(); dup = 0; for (int j = 0; j < i; j++) if (((u64)PL[j] << 32 | PR[j]) == p) { dup = 1; break; } } while (dup);
                PL[i] = p >> 32; PR[i] = (u32)p;
            }
        }
        /* 真値と参照統計量 */
        static u32 T[1024][4];
        double sA[2] = {0, 0}, s4[2] = {0, 0}, Sg[2][256] = {{0}}; int Ng[256] = {0};
        for (int i = 0; i < n; i++) {
            u32 st[18]; u64 c = feal8_enc(&s, ((u64)PL[i] << 32) | PR[i], st);
            CL[i] = c >> 32; CR[i] = (u32)c;
            T[i][0] = st[2]; T[i][1] = st[3]; T[i][2] = st[14]; T[i][3] = st[15];
            u32 w = feal_M(st[13]) ^ ((u32)s.k[6] << 8);
            int g = part_idx(BY(w, 1), BY(w, 2)); Ng[g]++;
            for (int f = 0; f < 2; f++) {
                int q = par(A8_L1 & st[2]) ^ par(a8_R1(f) & st[3]);
                int a = q ^ par(A8_Y6[f] & st[12]) ^ par(0x00010504u & st[13]);
                sA[f] += a ? -1 : 1;
                if (i < 256) s4[f] += a ? -1 : 1;
                int qp = q ^ par(A8_Y6[f] & st[15]) ^ par(0x00010504u & st[14]) ^ ((BY(w, 0) ^ BY(w, 1)) & 1) ^ ((BY(w, 3) ^ BY(w, 2)) & 1);
                Sg[f][g] += qp ? -1 : 1;
            }
        }
        double zA[2], zB[2], pz[2];
        for (int f = 0; f < 2; f++) {
            double num = 0, den = 0;
            for (int g = 0; g < 256; g++) { num += etab[f][g] * Sg[f][g]; den += etab[f][g] * etab[f][g] * Ng[g]; }
            zA[f] = sA[f] / sqrt(n); zB[f] = num / sqrt(den); pz[f] = s4[f] / 16.0;
        }
        u32 K89 = ((u32)s.k[8] << 16) | s.k[9], Kab = ((u32)s.k[10] << 16) | s.k[11];
        u32 Kcd = ((u32)s.k[12] << 16) | s.k[13], Kef = ((u32)s.k[14] << 16) | s.k[15];
        u32 EK0 = feal_M(K89 ^ Kab) ^ ((u32)s.k[0] << 8);
        u32 DK7 = feal_M(Kcd ^ Kef) ^ ((u32)s.k[7] << 8);
        int e1 = BY(EK0, 1), e2 = BY(EK0, 2);
        if (e2 & 0x80) { e1 ^= 0x80; e2 ^= 0x80; }
        int e_true = e1 | (e2 << 8);
        u32 t_true = BY(DK7, 1) | (BY(DK7, 2) << 8) | ((BY(DK7, 0) & 3) << 16) | ((BY(DK7, 3) & 15) << 18);
        /* 平文側 2^15 */
        #pragma omp parallel for schedule(static)
        for (int e = 0; e < 32768; e++) {
            u32 ek = PK(0, (u8)e, (u8)(e >> 8), 0);
            u64 v0[W] = {0}, v1[W] = {0};
            for (int i = 0; i < n; i++) {
                u32 L1 = PL[i] ^ PR[i], R1 = PL[i] ^ feal_N(feal_M(L1) ^ ek);
                int b = par(A8_L1 & L1);
                if (b ^ par(a8_R1(0) & R1)) v0[i >> 6] |= 1ULL << (i & 63);
                if (b ^ par(a8_R1(1) & R1)) v1[i >> 6] |= 1ULL << (i & 63);
            }
            memcpy(P[e][0], v0, sizeof v0); memcpy(P[e][1], v1, sizeof v1);
        }
        /* 暗号文側：0 番が真、1..S が無作為（真と同じものは除く） */
        u32 *tl = malloc(sizeof(u32) * (S + 1));
        tl[0] = t_true;
        for (int j = 1; j <= S; j++) { u32 t; do t = (u32)(rnd() & 0x3FFFFF); while (t == t_true); tl[j] = t; }
        #pragma omp parallel for schedule(static)
        for (int j = 0; j <= S; j++) {
            u32 dk = dk7_22(tl[j]);
            u64 v[8][W]; memset(v, 0, sizeof v);
            for (int i = 0; i < n; i++) {
                u32 R7 = CL[i] ^ CR[i], L7 = CL[i] ^ feal_N(feal_M(R7) ^ dk);
                for (int l = 0; l < 8; l++)
                    if (par(a8_L7(l & 3) & L7) ^ par(a8_R7(l >> 2) & R7)) v[l][i >> 6] |= 1ULL << (i & 63);
            }
            memcpy(C[j], v, sizeof v);
        }
        /* 等価表現の確認：真の (e,t) で計算したビットと真の中間値のビットの差が一定 */
        int consistent = 1;
        for (int l = 0; l < 8; l++) {
            int d0 = -1;
            for (int i = 0; i < n; i++) {
                int bc = (int)((P[e_true][l >> 2][i >> 6] ^ C[0][l][i >> 6]) >> (i & 63)) & 1;
                int bt = a8_bit(l, T[i][0], T[i][1], T[i][2], T[i][3]);
                if (d0 < 0) d0 = bc ^ bt; else if ((bc ^ bt) != d0) consistent = 0;
            }
        }
        long long strue[NS];
        scores(P[e_true][0], P[e_true][1], C[0], strue);
        long long gt_ex[NS] = {0}, eq_ex[NS] = {0}, gt_s[NS] = {0}, eq_s[NS] = {0};
        #pragma omp parallel
        {
            long long g1[NS] = {0}, q1[NS] = {0}, g2[NS] = {0}, q2[NS] = {0};
            #pragma omp for schedule(dynamic, 64)
            for (int e = 0; e < 32768; e++) {
                for (int j = 0; j <= S; j++) {
                    if (j == 0 && e == e_true) continue;
                    long long sc[NS];
                    scores(P[e][0], P[e][1], C[j], sc);
                    for (int k = 0; k < NS; k++) {
                        if (j == 0) { g1[k] += sc[k] > strue[k]; q1[k] += sc[k] == strue[k]; }
                        else { g2[k] += sc[k] > strue[k]; q2[k] += sc[k] == strue[k]; }
                    }
                }
            }
            #pragma omp critical
            for (int k = 0; k < NS; k++) { gt_ex[k] += g1[k]; eq_ex[k] += q1[k]; gt_s[k] += g2[k]; eq_s[k] += q2[k]; }
        }
        double scale = 4194304.0 / S;
        printf("%d %016llx %d %.3f %.3f %.3f %.3f %.2f %.2f | %lld |", tr, (unsigned long long)key, n,
               zA[0], zA[1], zB[0], zB[1], pz[0], pz[1], strue[0]);
        for (int k = 0; k < NS; k++) {
            double rank = 1.0 + gt_ex[k] + 0.5 * eq_ex[k] + scale * (gt_s[k] + 0.5 * eq_s[k]);
            printf(" %.2f", log2(rank));
        }
        printf(" | %s\n", consistent ? "ok" : "NG");
        fflush(stdout);
        free(tl);
    }
    return 0;
}
