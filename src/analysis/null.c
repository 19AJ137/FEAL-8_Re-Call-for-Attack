/* 誤り鍵の統計量の分布（等価表現で計算）
 *   ./null real plaintext.bin keys seed ncipher
 *   ./null rand N             keys seed ncipher
 * 平文側は EK0 の中央2バイト（2^16 通り、ほかは真値）を全数。暗号文側は
 *   (i)  真の DK7/DK6
 *   (ii) 乱数の DK7/DK6 を ncipher 個（完全な誤り仮説）
 * 方式A：z = Σ(-1)^{p^c}/√N。分割：z = Σ e_g(-1)^{p^c}/√(Σe_g^2)
 */
#include "feal8.h"

#define W 16 /* 1024 ビット */

static void eqkeys(const sched_t *s, u32 *EK0, u32 *DK7, u32 *DK6) {
    u32 K89 = ((u32)s->k[8] << 16) | s->k[9], Kab = ((u32)s->k[10] << 16) | s->k[11];
    u32 Kcd = ((u32)s->k[12] << 16) | s->k[13], Kef = ((u32)s->k[14] << 16) | s->k[15];
    *EK0 = feal_M(K89 ^ Kab) ^ ((u32)s->k[0] << 8);
    *DK7 = feal_M(Kcd ^ Kef) ^ ((u32)s->k[7] << 8);
    *DK6 = feal_M(Kcd) ^ ((u32)s->k[6] << 8);
}

static int n;
static u64 *pt, *ct;
static u64 (*pv)[W];   /* 平文側ビット列（EK0 の中央2バイトごと） */
static float (*ps)[1024]; /* 同じものを ±1 で */

static void cipher_side(u32 DK7, u32 DK6, u64 *cv, float *wv) {
    memset(cv, 0, sizeof(u64) * W);
    double den = 0;
    for (int i = 0; i < n; i++) {
        u32 CL = ct[i] >> 32, CR = (u32)ct[i];
        u32 R7 = CL ^ CR, L7 = CL ^ feal_N(feal_M(R7) ^ DK7);
        u32 w = feal_M(L7) ^ DK6;
        u32 L6 = R7 ^ feal_N(w);
        int c = par(G_L6 & L6) ^ par(G_R6 & L7);
        if (c) cv[i >> 6] |= 1ULL << (i & 63);
        /* 分割用：Q' の暗号文側 = <G_L6,R7>^<G_R6,L7>^(w0^w1)_0^(w3^w2)_0 */
        int cq = par(G_L6 & R7) ^ par(G_R6 & L7) ^ ((BY(w, 0) ^ BY(w, 1)) & 1) ^ ((BY(w, 3) ^ BY(w, 2)) & 1);
        double e = e_tab[part_idx(BY(w, 1), BY(w, 2))];
        wv[i] = (float)(cq ? -e : e);
        den += e * e;
    }
    float sc = (float)(1.0 / sqrt(den));
    for (int i = 0; i < n; i++) wv[i] *= sc;
}

static inline double zA(const u64 *a, const u64 *b) {
    int d = 0;
    for (int j = 0; j < W; j++) d += __builtin_popcountll(a[j] ^ b[j]);
    return (n - 2.0 * d) / sqrt(n);
}
static inline double zB(const float *a, const float *b) {
    float s = 0;
    for (int i = 0; i < n; i++) s += a[i] * b[i];
    return s;
}

int main(int argc, char **argv) {
    build_etab();
    int real = !strcmp(argv[1], "real");
    int keys = atoi(argv[3]);
    seed(strtoull(argv[4], 0, 0));
    int ncip = atoi(argv[5]);
    if (real) pt = load_unique(argv[2], &n);
    else { n = atoi(argv[2]); pt = malloc(8 * n); }
    ct = malloc(8 * n);
    pv = malloc(sizeof(*pv) * 65536);
    ps = malloc(sizeof(*ps) * 65536);
    double thr[] = {3.0, 4.0, 4.5, 5.0, 5.5, 6.0};
    long long hA[6] = {0}, hB[6] = {0}, totA = 0, totB = 0;
    double sA = 0, sA2 = 0, sB = 0, sB2 = 0, maxA = 0, maxB = 0;
    for (int t = 0; t < keys; t++) {
        sched_t s; key_schedule(rnd(), &s);
        if (!real) for (int i = 0; i < n; i++) pt[i] = rnd();
        for (int i = 0; i < n; i++) ct[i] = feal8_enc(&s, pt[i], NULL);
        u32 EK0, DK7, DK6; eqkeys(&s, &EK0, &DK7, &DK6);
        for (int e = 0; e < 65536; e++) {
            u32 ek = (EK0 & 0xFF0000FFu) | ((u32)e << 8);
            memset(pv[e], 0, sizeof pv[e]);
            for (int i = 0; i < n; i++) {
                u32 PL = pt[i] >> 32, PR = (u32)pt[i];
                u32 L1 = PL ^ PR, R1 = PL ^ feal_N(feal_M(L1) ^ ek);
                int p = par(G_L1 & L1) ^ par(G_R1 & R1);
                if (p) pv[e][i >> 6] |= 1ULL << (i & 63);
                ps[e][i] = p ? -1.f : 1.f;
            }
        }
        int etrue = (EK0 >> 8) & 0xFFFF;
        /* (i) 正しい暗号文側 × 全 EK0 */
        u64 cv[W]; float wv[1024];
        cipher_side(DK7, DK6, cv, wv);
        double ztrue = zA(pv[etrue], cv), zbtrue = zB(ps[etrue], wv);
        int same = 0, hi65 = 0, hiB3 = 0; double maxw = 0, maxwB = 0;
        for (int e = 0; e < 65536; e++) {
            int d = 0;
            for (int j = 0; j < W; j++) d += __builtin_popcountll(pv[e][j] ^ pv[etrue][j]);
            int eq = d == 0, cp = d == n;
            double z = fabs(zA(pv[e], cv)), zb = fabs(zB(ps[e], wv));
            if (eq || cp) { same++; continue; }
            if (z >= 6.5) hi65++;
            if (zb >= 3.0) hiB3++;
            if (z > maxw) maxw = z;
            if (zb > maxwB) maxwB = zb;
        }
        printf("key %d: true zA=%.2f zB=%.2f | EK0 values giving identical/complement p-vector: %d of 65536 | "
               "wrong EK0 with right cipher side: max zA %.2f, #>=6.5: %d ; max zB %.2f, #zB>=3: %d\n",
               t, fabs(ztrue), fabs(zbtrue), same, maxw, hi65, maxwB, hiB3);
        /* (ii) 完全に誤った仮説 */
        for (int c = 0; c < ncip; c++) {
            cipher_side((u32)rnd(), (u32)rnd(), cv, wv);
            for (int e = 0; e < 65536; e++) {
                double z = fabs(zA(pv[e], cv));
                sA += z * z; totA++;
                if (z > maxA) maxA = z;
                for (int k = 0; k < 6; k++) hA[k] += z >= thr[k];
            }
            if (c < ncip / 16 + 1) {
                for (int e = 0; e < 65536; e++) {
                    double z = fabs(zB(ps[e], wv));
                    sB += z * z; totB++;
                    if (z > maxB) maxB = z;
                    for (int k = 0; k < 6; k++) hB[k] += z >= thr[k];
                }
            }
        }
    }
    printf("\nnull (fully wrong) %s N=%d\n", argv[1], n);
    printf("  A: %lld hyps, E[z^2]=%.4f max %.2f\n", totA, sA / totA, maxA);
    printf("  B: %lld hyps, E[z^2]=%.4f max %.2f\n", totB, sB / totB, maxB);
    printf("  thr   P_A(obs)    P_B(obs)    P_normal\n");
    for (int k = 0; k < 6; k++)
        printf("  %.1f  %.3e  %.3e  %.3e\n", thr[k], (double)hA[k] / totA, (double)hB[k] / totB, erfc(thr[k] / sqrt(2)));
    return 0;
}
