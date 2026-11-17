/* 正解鍵での統計量（方式A・6段近似1本・8ビット分割）を測る
 *   ./corr check                      テストベクタ、等価表現、近似の偏り、分割の容量
 *   ./corr rand  N    keys seed       一様乱数の平文 N 個
 *   ./corr real  file keys seed [a b] 実データの平文（重複除去）。a..b でユニーク番号の範囲を指定
 */
#include "feal8.h"

static void eqkeys(const sched_t *s, u32 *EK0, u32 *DK7, u32 *DK6) {
    u32 K89 = ((u32)s->k[8] << 16) | s->k[9], Kab = ((u32)s->k[10] << 16) | s->k[11];
    u32 Kcd = ((u32)s->k[12] << 16) | s->k[13], Kef = ((u32)s->k[14] << 16) | s->k[15];
    *EK0 = feal_M(K89 ^ Kab) ^ ((u32)s->k[0] << 8);
    *DK7 = feal_M(Kcd ^ Kef) ^ ((u32)s->k[7] << 8);
    *DK6 = feal_M(Kcd) ^ ((u32)s->k[6] << 8);
}

static void check(void) {
    sched_t s;
    key_schedule(0x0123456789ABCDEFULL, &s);
    u64 c = feal8_enc(&s, 0, NULL);
    printf("test vector: %016llx (expect ceef2c86f2490752) %s\n", (unsigned long long)c,
           c == 0xceef2c86f2490752ULL ? "OK" : "NG");

    /* 等価表現：真の中間値との差が平文によらず一定か */
    seed(1);
    int ok = 1;
    for (int t = 0; t < 200 && ok; t++) {
        key_schedule(rnd(), &s);
        u32 EK0, DK7, DK6; eqkeys(&s, &EK0, &DK7, &DK6);
        u32 d[4] = {0};
        for (int i = 0; i < 64; i++) {
            u64 p = rnd(); u32 st[18];
            u64 cc = feal8_enc(&s, p, st);
            u32 PL = p >> 32, PR = (u32)p, CL = cc >> 32, CR = (u32)cc;
            u32 L1 = PL ^ PR, R1 = PL ^ feal_N(feal_M(L1) ^ EK0);
            u32 R7 = CL ^ CR, L7 = CL ^ feal_N(feal_M(R7) ^ DK7);
            u32 L6 = R7 ^ feal_N(feal_M(L7) ^ DK6);
            u32 e[4] = {L1 ^ st[2], R1 ^ st[3], L7 ^ st[14], L6 ^ st[12]};
            for (int j = 0; j < 4; j++) {
                if (i == 0) d[j] = e[j]; else if (d[j] != e[j]) ok = 0;
            }
        }
    }
    printf("equivalent representation (EK0/DK7/DK6, 200 keys): %s\n", ok ? "OK" : "NG");

    /* 近似の相関（一様乱数 2^24 個、1 鍵） */
    key_schedule(rnd(), &s);
    long long a5 = 0, a6 = 0, r6 = 0;
    const long long NN = 1LL << 24;
    for (long long i = 0; i < NN; i++) {
        u32 st[18];
        feal8_enc(&s, rnd(), st);
        int q = par(G_L1 & st[2]) ^ par(G_R1 & st[3]);
        int q5 = q ^ par(G_L6 & st[12]) ^ par(G_R6 & st[13]);
        int q6 = q ^ par(G_L7 & st[14]) ^ par(G_R7 & st[15]);
        a5 += q5 ? -1 : 1; a6 += q6 ? -1 : 1;
        (void)r6;
    }
    printf("5-round approx corr %.4f (expect +-0.25), 6-round approx corr %.5f (expect +-1/32=0.03125)\n",
           (double)a5 / NN, (double)a6 / NN);

    build_etab();
    double cap = 0;
    for (int g = 0; g < 256; g++) cap += e_tab[g] * e_tab[g];
    printf("8-bit partition capacity (mean e_g^2) = %.4f (survey: 0.200)\n", cap / 256);
}

typedef struct { double zA, z6, zB, zA256, rA, rB; } zs_t;

static zs_t stats(const sched_t *s, const u64 *pt, int n) {
    long long sA = 0, s6 = 0, sA256 = 0;
    double Sg[256] = {0}; int Ng[256] = {0};
    double Tg[256] = {0}; int Mg[256] = {0}; long long rA = 0;
    for (int i = 0; i < n; i++) {
        u32 st[18];
        feal8_enc(s, pt[i], st);
        int q = par(G_L1 & st[2]) ^ par(G_R1 & st[3]);
        int qA = q ^ par(G_L6 & st[12]) ^ par(G_R6 & st[13]);
        int q6 = q ^ par(G_L7 & st[14]) ^ par(G_R7 & st[15]);
        sA += qA ? -1 : 1; s6 += q6 ? -1 : 1;
        if (i == 255) sA256 = sA;
        /* 分割：w = M(R6) ^ (0,K6,0)、Q' = A ^ w0_0 ^ w3_0 */
        u32 R6 = st[13], R7 = st[15], L7 = st[14];
        u32 w = feal_M(R6) ^ ((u32)s->k[6] << 8);
        int qp = q ^ par(G_L6 & R7) ^ par(G_R6 & L7) ^ ((BY(w, 0) ^ BY(w, 1)) & 1) ^ ((BY(w, 3) ^ BY(w, 2)) & 1);
        int g = part_idx(BY(w, 1), BY(w, 2));
        Sg[g] += qp ? -1 : 1; Ng[g]++;
        /* 逆向き：5段近似を第2〜6段、分割を第1段（平文側）に置く */
        {
            u32 L1 = st[2], R1 = st[3], L2 = st[4], R2 = st[5];
            int qr = par(G_L1 & L2) ^ par(G_R1 & R2) ^ par(G_L6 & L7) ^ par(G_R6 & R7);
            rA += qr ? -1 : 1;
            u32 v = feal_M(R1) ^ ((u32)s->k[1] << 8);
            int qq = par(G_L6 & L7) ^ par(G_R6 & R7) ^ par(G_L1 & R1) ^ par(G_R1 & L1)
                   ^ ((BY(v, 0) ^ BY(v, 1)) & 1) ^ ((BY(v, 3) ^ BY(v, 2)) & 1);
            int h = part_idx(BY(v, 1), BY(v, 2));
            Tg[h] += qq ? -1 : 1; Mg[h]++;
        }
    }
    double num = 0, den = 0;
    for (int g = 0; g < 256; g++) { num += e_tab[g] * Sg[g]; den += e_tab[g] * e_tab[g] * Ng[g]; }
    double num2 = 0, den2 = 0;
    for (int g = 0; g < 256; g++) { num2 += e_tab[g] * Tg[g]; den2 += e_tab[g] * e_tab[g] * Mg[g]; }
    zs_t z = { fabs((double)sA / sqrt(n)), fabs((double)s6 / sqrt(n)), fabs(num / sqrt(den)), fabs((double)sA256 / 16.0),
               fabs((double)rA / sqrt(n)), fabs(num2 / sqrt(den2)) };
    return z;
}

static void report(const char *name, double *v, int k, const double *th, int nth) {
    double m = 0, m2 = 0;
    for (int i = 0; i < k; i++) { m += v[i]; m2 += v[i] * v[i]; }
    m /= k; m2 = sqrt(m2 / k - m * m);
    printf("  %-4s |z| mean %.2f sd %.2f", name, m, m2);
    for (int t = 0; t < nth; t++) {
        int c = 0;
        for (int i = 0; i < k; i++) c += v[i] >= th[t];
        printf("  P(|z|>=%.1f)=%.3f", th[t], (double)c / k);
    }
    printf("\n");
}

int main(int argc, char **argv) {
    if (argc < 2) return 1;
    if (!strcmp(argv[1], "check")) { check(); return 0; }
    build_etab();
    int keys = atoi(argv[3]);
    seed(strtoull(argv[4], 0, 0));
    int n; u64 *pt;
    int real = !strcmp(argv[1], "real");
    if (real) {
        pt = load_unique(argv[2], &n);
        if (argc >= 7) { int a = atoi(argv[5]), b = atoi(argv[6]); pt += a; n = b - a; }
    } else {
        n = atoi(argv[2]);
        pt = malloc(sizeof(u64) * n);
    }
    int jB3A65 = 0, jB3e2A65 = 0, jB35A65 = 0, rev = 0, either = 0;
    double rAm = 0, rBm = 0;
    double *A = malloc(sizeof(double) * keys), *B = malloc(sizeof(double) * keys), *C6 = malloc(sizeof(double) * keys);
    for (int t = 0; t < keys; t++) {
        sched_t s; key_schedule(rnd(), &s);
        if (!real) for (int i = 0; i < n; i++) pt[i] = rnd();
        zs_t z = stats(&s, pt, n);
        A[t] = z.zA; B[t] = z.zB; C6[t] = z.z6;
        jB3A65 += z.zB >= 3.0 && z.zA >= 6.5;
        jB3e2A65 += z.zB >= 3.0 && z.zA256 >= 2.0 && z.zA >= 6.5;
        jB35A65 += z.zB >= 3.5 && z.zA >= 6.5;
        int fw = z.zB >= 3.0 && z.zA >= 6.5, rv = z.rB >= 3.0 && z.rA >= 6.5;
        rev += rv; either += fw || rv; rAm += z.rA; rBm += z.rB;
    }
    printf("%s N=%d keys=%d  (theory: A=sqrt(N)/4=%.2f, B=%.2f)\n", argv[1], n, keys, sqrt(n) / 4,
           0.25 * sqrt(n * 0.2));
    double thA[] = {6.5, 7.0}, thB[] = {3.0, 3.5}, th6[] = {2.0};
    report("A", A, keys, thA, 2);
    report("B", B, keys, thB, 2);
    report("6r", C6, keys, th6, 1);
    printf("  joint: P(B>=3.0 & A>=6.5)=%.3f  P(B>=3.0 & A256>=2 & A>=6.5)=%.3f  P(B>=3.5 & A>=6.5)=%.3f\n",
           (double)jB3A65 / keys, (double)jB3e2A65 / keys, (double)jB35A65 / keys);
    printf("  reverse: mean |zA| %.2f |zB| %.2f  P(rev success)=%.3f  P(fwd or rev)=%.3f\n",
           rAm / keys, rBm / keys, (double)rev / keys, (double)either / keys);
    return 0;
}
