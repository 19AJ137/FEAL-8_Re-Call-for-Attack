/* 方式C の正解仮説の統計量を、第1族（現行）と第2族（論文の ΓY2=04011904）で同時に出す。
 * 真の中間値を使う（attackC.c の reference scores と同じ定義）。順方向と逆方向（平文と暗号文を入れ替え）。
 *   ./refstat plaintext.bin nkeys seed  > out
 *   KEY=a5ea6ec2fb0cc741 ./refstat plaintext.bin 1 0   （指定した鍵 1 本だけ）
 * 出力：key  fwd: zA0 zA1 zB0 zB1 pre0 pre1   rev: 同じ 6 個 */
#include "feal8.h"
#include "approx8.h"

static double etab[2][256];
static void build_etab2(void) {
    for (int f = 0; f < 2; f++) {
        double acc[256] = {0}; int cnt[256] = {0};
        for (int w1 = 0; w1 < 256; w1++)
            for (int w2 = 0; w2 < 256; w2++) {
                int g = part_idx((u8)w1, (u8)w2);
                acc[g] += (par(A8_Y6[f] & feal_N(PK(0, (u8)w1, (u8)w2, 0))) ^ (w1 & 1) ^ (w2 & 1)) ? -1 : 1;
                cnt[g]++;
            }
        for (int g = 0; g < 256; g++) etab[f][g] = acc[g] / cnt[g];
    }
}

static void stats(const sched_t *s, const u64 *pt, int n, double out[6]) {
    double sA[2] = {0, 0}, s4[2] = {0, 0}, Sg[2][256] = {{0}}; int Ng[256] = {0};
    for (int i = 0; i < n; i++) {
        u32 st[18]; feal8_enc(s, pt[i], st);
        u32 w = feal_M(st[13]) ^ ((u32)s->k[6] << 8);
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
    for (int f = 0; f < 2; f++) {
        double num = 0, den = 0;
        for (int g = 0; g < 256; g++) { num += etab[f][g] * Sg[f][g]; den += etab[f][g] * etab[f][g] * Ng[g]; }
        out[f] = sA[f] / sqrt(n); out[2 + f] = num / sqrt(den); out[4 + f] = s4[f] / 16.0;
    }
}

int main(int argc, char **argv) {
    int n; u64 *pt = load_unique(argv[1], &n);
    int nk = atoi(argv[2]); seed(strtoull(argv[3], 0, 0));
    const char *fixed = getenv("KEY");   /* KEY=<16 hex digits>: report this key only (nkeys, seed ignored) */
    if (fixed) nk = 1;
    build_etab2();
    u64 ct[1024];
    for (int k = 0; k < nk; k++) {
        u64 key = fixed ? strtoull(fixed, 0, 16) : rnd();
        sched_t s, r; key_schedule(key, &s);
        /* 入れ替えた暗号：段鍵を逆順、白色化鍵 K89<->Kcd、Kab<->Kef */
        for (int i = 0; i < 8; i++) r.k[i] = s.k[7 - i];
        r.k[8] = s.k[12]; r.k[9] = s.k[13]; r.k[10] = s.k[14]; r.k[11] = s.k[15];
        r.k[12] = s.k[8]; r.k[13] = s.k[9]; r.k[14] = s.k[10]; r.k[15] = s.k[11];
        for (int i = 0; i < n; i++) {
            ct[i] = feal8_enc(&s, pt[i], NULL);
            if (k < 4 && feal8_enc(&r, ct[i], NULL) != pt[i]) { fprintf(stderr, "reverse cipher mismatch\n"); return 1; }
        }
        double f[6], b[6];
        stats(&s, pt, n, f); stats(&r, ct, n, b);
        printf("%016llx", (unsigned long long)key);
        for (int j = 0; j < 6; j++) printf(" %.4f", f[j]);
        for (int j = 0; j < 6; j++) printf(" %.4f", b[j]);
        printf("\n");
    }
    return 0;
}
