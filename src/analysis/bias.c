/* 8本の 6段近似の相関を真の中間値で測る。あわせて
 *  - 5段部分（第1〜5段）の fam0 / fam1 の相関
 *  - 8本の指示ビットどうしの相関（鍵ごとの符号をそろえたもの）
 *   ./bias nkeys log2texts seed */
#include "feal8.h"
#include "approx8.h"

int main(int argc, char **argv) {
    int nk = atoi(argv[1]), lt = atoi(argv[2]);
    seed(strtoull(argv[3], 0, 0));
    long long NT = 1LL << lt;
    double sumabs[8] = {0}, sum5[2] = {0}, cross[8][8] = {{0}}, x01 = 0;
    for (int k = 0; k < nk; k++) {
        sched_t s; key_schedule(rnd(), &s);
        long long c[8] = {0}, c5[2] = {0}, cc[8][8] = {{0}}, cx = 0;
        for (long long t = 0; t < NT; t++) {
            u32 st[18]; feal8_enc(&s, rnd(), st);
            u32 L1 = st[2], R1 = st[3], L6 = st[12], R6 = st[13], L7 = st[14], R7 = st[15];
            int b[8];
            for (int l = 0; l < 8; l++) { b[l] = a8_bit(l, L1, R1, L7, R7); c[l] += b[l] ? -1 : 1; }
            for (int f = 0; f < 2; f++) {
                int q = par(A8_L1 & L1) ^ par(a8_R1(f) & R1) ^ par(A8_Y6[f] & L6) ^ par(0x00010504u & R6);
                c5[f] += q ? -1 : 1;
            }
            int q0 = par(A8_L1 & L1) ^ par(a8_R1(0) & R1) ^ par(A8_Y6[0] & L6) ^ par(0x00010504u & R6);
            int q1 = par(A8_L1 & L1) ^ par(a8_R1(1) & R1) ^ par(A8_Y6[1] & L6) ^ par(0x00010504u & R6);
            cx += (q0 ^ q1) ? -1 : 1;
            for (int i = 0; i < 8; i++) for (int j = i + 1; j < 8; j++) cc[i][j] += (b[i] ^ b[j]) ? -1 : 1;
        }
        double cr[8];
        for (int l = 0; l < 8; l++) { cr[l] = (double)c[l] / NT; sumabs[l] += fabs(cr[l]); }
        for (int f = 0; f < 2; f++) sum5[f] += fabs((double)c5[f] / NT);
        x01 += fabs((double)cx / NT);
        /* 符号をそろえた相関：E[(-1)^(bi^bj)] から sign(ci)sign(cj) ci cj を引いた残り＝指示ビットの共分散 */
        for (int i = 0; i < 8; i++) for (int j = i + 1; j < 8; j++) {
            double e = (double)cc[i][j] / NT;
            double sg = (cr[i] < 0) != (cr[j] < 0) ? -1 : 1;
            cross[i][j] += sg * e;
        }
    }
    printf("keys=%d texts/key=2^%d\n", nk, lt);
    printf("approx  ΓY2      ΓX6       |corr|     bias(x2^-12)  paper ε\n");
    for (int l = 0; l < 8; l++)
        printf("  %d     %08x %08x  %.5f   %6.1f         %d\n", l + 1, A8_Y2[l >> 2], A8_X6[l & 3],
               sumabs[l] / nk, sumabs[l] / nk / 2 * 4096, A8_EPS[l & 3]);
    printf("5-round part |corr|: fam0 %.4f  fam1 %.4f ; corr of (q0 xor q1) %.4f (independent would be %.4f)\n",
           sum5[0] / nk, sum5[1] / nk, x01 / nk, (sum5[0] / nk) * (sum5[1] / nk));
    printf("sign-aligned E[(-1)^(bi^bj)] (independent would be |ci||cj| ~ 0.001-0.002):\n");
    for (int i = 0; i < 8; i++) {
        printf("  %d:", i + 1);
        for (int j = 0; j < 8; j++) printf(j > i ? " %7.4f" : "        ", j > i ? cross[i][j] / nk : 0.0);
        printf("\n");
    }
    return 0;
}
