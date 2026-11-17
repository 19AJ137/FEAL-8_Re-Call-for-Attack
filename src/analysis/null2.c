/* 方式C 段階1（8ビット分割、κ 256 通りを WHT で一括）の誤り仮説の分布を、第1族・第2族で同時に測る。
 * 実データの平文、乱数鍵、無作為な DK7(22bit)、EK0 全 2^15、κ 全 256。
 *   ./null2 plaintext.bin nkeys nd seed */
#include "feal8.h"
#include "approx8.h"

static double etab[2][256];
static float Eh[2][256], E2h[2][256];
static void wht(float *a) {
    for (int h = 1; h < 256; h <<= 1)
        for (int i = 0; i < 256; i += 2 * h)
            for (int j = i; j < i + h; j++) { float x = a[j], y = a[j + h]; a[j] = x + y; a[j + h] = x - y; }
}
static inline u32 dk7_22(u32 t) {
    return PK((u8)((t >> 16) & 3), (u8)t, (u8)(t >> 8), (u8)((t >> 18) & 15));
}

int main(int argc, char **argv) {
    int n; u64 *pt = load_unique(argv[1], &n);
    int nk = atoi(argv[2]), nd = atoi(argv[3]); seed(strtoull(argv[4], 0, 0));
    for (int f = 0; f < 2; f++) {
        double acc[256] = {0}; int cnt[256] = {0};
        for (int w1 = 0; w1 < 256; w1++)
            for (int w2 = 0; w2 < 256; w2++) {
                int g = part_idx((u8)w1, (u8)w2);
                acc[g] += (par(A8_Y6[f] & feal_N(PK(0, (u8)w1, (u8)w2, 0))) ^ (w1 & 1) ^ (w2 & 1)) ? -1 : 1;
                cnt[g]++;
            }
        for (int g = 0; g < 256; g++) { etab[f][g] = acc[g] / cnt[g]; Eh[f][g] = (float)etab[f][g]; E2h[f][g] = (float)(etab[f][g] * etab[f][g]); }
        wht(Eh[f]); wht(E2h[f]);
    }
    static u32 PLa[1024], PRa[1024], CLa[1024], CRa[1024];
    static u8 pb[32768][2][1024];
    const double R[] = {3.0, 3.4395, 3.9175, 4.5, 5.0};
    double cnt_n = 0, s0 = 0, s1 = 0, s00 = 0, s11 = 0, s01 = 0, t35 = 0, chi[5] = {0};
    for (int k = 0; k < nk; k++) {
        sched_t s; key_schedule(rnd(), &s);
        for (int i = 0; i < n; i++) {
            u64 c = feal8_enc(&s, pt[i], NULL);
            PLa[i] = pt[i] >> 32; PRa[i] = (u32)pt[i]; CLa[i] = c >> 32; CRa[i] = (u32)c;
        }
        for (int e = 0; e < 32768; e++) {
            u32 ek = PK(0, (u8)e, (u8)(e >> 8), 0);
            for (int i = 0; i < n; i++) {
                u32 L1 = PLa[i] ^ PRa[i], R1 = PLa[i] ^ feal_N(feal_M(L1) ^ ek);
                int b = par(A8_L1 & L1);
                pb[e][0][i] = b ^ par(a8_R1(0) & R1);
                pb[e][1][i] = b ^ par(a8_R1(1) & R1);
            }
        }
        for (int d = 0; d < nd; d++) {
            u32 DK7 = dk7_22((u32)(rnd() & 0x3FFFFF));
            u8 u[1024], c[2][1024]; float Ng[256] = {0}, sd[2][256];
            for (int i = 0; i < n; i++) {
                u32 R7 = CLa[i] ^ CRa[i], L7 = CLa[i] ^ feal_N(feal_M(R7) ^ DK7);
                for (int f = 0; f < 2; f++) c[f][i] = par(A8_Y6[f] & R7) ^ par(0x00000404u & L7);
                u[i] = (u8)part_idx(BY(L7, 0) ^ BY(L7, 1), BY(L7, 2) ^ BY(L7, 3));
                Ng[u[i]] += 1;
            }
            wht(Ng);
            for (int f = 0; f < 2; f++) {
                float t[256];
                for (int g = 0; g < 256; g++) t[g] = Ng[g] * E2h[f][g];
                wht(t);
                for (int q = 0; q < 256; q++) sd[f][q] = sqrtf(t[q] / 256.f) * 256.f;
            }
            for (int e = 0; e < 32768; e++) {
                float T[2][256] = {{0}};
                for (int i = 0; i < n; i++)
                    for (int f = 0; f < 2; f++) T[f][u[i]] += (pb[e][f][i] ^ c[f][i]) ? -1.f : 1.f;
                for (int f = 0; f < 2; f++) {
                    wht(T[f]);
                    for (int g = 0; g < 256; g++) T[f][g] *= Eh[f][g];
                    wht(T[f]);
                }
                for (int q = 0; q < 256; q++) {
                    double z0 = T[0][q] / sd[0][q], z1 = T[1][q] / sd[1][q];
                    cnt_n++; s0 += z0; s1 += z1; s00 += z0 * z0; s11 += z1 * z1; s01 += z0 * z1;
                    t35 += fabs(z0) >= 3.5;
                    double r2 = z0 * z0 + z1 * z1;
                    for (int j = 0; j < 5; j++) chi[j] += r2 >= R[j] * R[j];
                }
            }
        }
    }
    printf("samples %.3g (keys %d x DK7 %d x EK0 2^15 x kappa 256)\n", cnt_n, nk, nd);
    printf("fam0: mean %.4f var %.4f   fam1: mean %.4f var %.4f   corr(z0,z1) %.4f\n",
           s0 / cnt_n, s00 / cnt_n, s1 / cnt_n, s11 / cnt_n, s01 / cnt_n);
    printf("P(|z0|>=3.5) %.3e (normal %.3e)\n", t35 / cnt_n, erfc(3.5 / sqrt(2)));
    for (int j = 0; j < 5; j++)
        printf("P(z0^2+z1^2 >= %.3f^2) %.3e (chi2_2 %.3e)\n", R[j], chi[j] / cnt_n, exp(-R[j] * R[j] / 2));
    return 0;
}
