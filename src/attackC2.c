/* 方式C（2族版）：attackC.c に先小山ら（IEICE 2016）の第2族の近似（ΓY2=04011904、ΓY6=04011804）を足したもの
 *
 *   ./attackC2 self  plaintext.bin seed nd [r1 rpre r2]
 *   ./attackC2 real  plaintext.bin ciphertext.bin start count out.txt [r1 rpre r2]
 *
 * 仮説 (DK7 22bit, EK0 15bit, κ 8bit) と補完 (e8, r) は attackC.c と同じ。違いは統計量だけ：
 *   第 f 族（f=0: 現行、f=1: 第2族）
 *     平文側   p_f = <00010504,L1> ^ <a1_f,R1>      a1_0=04031004, a1_1=04011804
 *     暗号文側 c_f = <a6_f,R7> ^ <00000404,L7>      a6_0=04031004, a6_1=04011804
 *     分割表   e_f(g)、h'_f(w1,w2) は第6段の出力マスク a6_f で作る
 *   段階1：z0^2 + z1^2 >= r1^2（既定 3.9175 = 現行の |z|>=3.5 と同じ生存数）
 *   段階2：先頭256組で z0^2+z1^2 >= rpre^2（既定 2.486 = 現行の 2σ と同じ通過率）、
 *          全体で z0^2+z1^2 >= r2^2 を候補として出力
 * 2族の符号の関係は鍵によるかもしれないので、二乗和で合わせる。 */
#include "feal8.h"
#include "approx8.h"
#include <omp.h>
#include <sys/time.h>

#define W 16
static int n;
static u32 PLa[1024], PRa[1024], CLa[1024], CRa[1024];
static u64 (*Pv[2])[W];
static float Ehat[2][256], E2hat[2][256];
static u64 ROW[2][65536][2];
static double etab[2][256];
static double r1 = 3.9175, rpre = 2.486, r2 = 6.0;

static double now(void) { struct timeval tv; gettimeofday(&tv, 0); return tv.tv_sec + tv.tv_usec * 1e-6; }
static inline int hprime(int f, int w1, int w2) {
    return par(A8_Y6[f] & feal_N(PK(0, (u8)w1, (u8)w2, 0))) ^ (w1 & 1) ^ (w2 & 1);
}
static inline int partx(u32 x) { return part_idx((u8)(x >> 8), (u8)x); }
static inline u32 emb_part(int kp) { return ((u32)(kp & 7) << 10) | ((u32)((kp >> 3) & 31) << 2); }
static inline u32 emb_rest(int r) {
    u32 w1 = (r & 3) | (((r >> 2) & 7) << 5), w2 = (r >> 5) & 3;
    return (w1 << 8) | w2;
}
static inline u32 dk7_22(u32 t) {
    return PK((u8)((t >> 16) & 3), (u8)t, (u8)(t >> 8), (u8)((t >> 18) & 15));
}
static inline u32 dk7_ext(int e8) { return ((u32)(e8 & 31) << 26) | ((u32)((e8 >> 5) & 7) << 4); }

static void wht(float *a) {
    for (int h = 1; h < 256; h <<= 1)
        for (int i = 0; i < 256; i += 2 * h)
            for (int j = i; j < i + h; j++) { float x = a[j], y = a[j + h]; a[j] = x + y; a[j + h] = x - y; }
}

static void setup_tables(void) {
    for (int f = 0; f < 2; f++) {
        double acc[256] = {0}; int cnt[256] = {0};
        for (int w1 = 0; w1 < 256; w1++)
            for (int w2 = 0; w2 < 256; w2++) {
                int g = part_idx((u8)w1, (u8)w2);
                acc[g] += hprime(f, w1, w2) ? -1 : 1; cnt[g]++;
            }
        for (int g = 0; g < 256; g++) {
            etab[f][g] = acc[g] / cnt[g];
            Ehat[f][g] = (float)etab[f][g]; E2hat[f][g] = (float)(etab[f][g] * etab[f][g]);
        }
        wht(Ehat[f]); wht(E2hat[f]);
        static u8 HP[65536];
        for (u32 x = 0; x < 65536; x++) HP[x] = (u8)hprime(f, x >> 8, x & 255);
        for (u32 x = 0; x < 65536; x++) {
            u64 a = 0, b = 0;
            for (int r = 0; r < 128; r++) {
                int bit = HP[x ^ emb_rest(r)];
                if (r < 64) a |= (u64)bit << r; else b |= (u64)bit << (r - 64);
            }
            ROW[f][x][0] = a; ROW[f][x][1] = b;
        }
    }
}

static void build_plain_side(void) {
    for (int f = 0; f < 2; f++) Pv[f] = aligned_alloc(64, sizeof(*Pv[f]) * 32768);
    #pragma omp parallel for schedule(static)
    for (int e = 0; e < 32768; e++) {
        u32 ek = PK(0, (u8)e, (u8)(e >> 8), 0);
        u64 v0[W] = {0}, v1[W] = {0};
        for (int i = 0; i < n; i++) {
            u32 L1 = PLa[i] ^ PRa[i], R1 = PLa[i] ^ feal_N(feal_M(L1) ^ ek);
            int b = par(A8_L1 & L1);
            if (b ^ par(a8_R1(0) & R1)) v0[i >> 6] |= 1ULL << (i & 63);
            if (b ^ par(a8_R1(1) & R1)) v1[i >> 6] |= 1ULL << (i & 63);
        }
        memcpy(Pv[0][e], v0, sizeof v0); memcpy(Pv[1][e], v1, sizeof v1);
    }
}

static void transpose64(u64 *a) {
    int j = 32; u64 m = 0x00000000FFFFFFFFULL;
    for (; j != 0; j >>= 1, m ^= (m << j))
        for (int k = 0; k < 64; k = ((k | j) + 1) & ~j) {
            u64 t = ((a[k] >> j) ^ a[k | j]) & m;
            a[k] ^= t << j; a[k | j] ^= t;
        }
}

typedef struct { u32 t; u16 e; u8 kap, e8, r; float z, z0, z1; } cand_t;
typedef struct { double t1, t2; long long surv, cmp_pre, cmp_full; int ncand; } dstat_t;

static int process_d(u32 t, cand_t *cand, int maxc, dstat_t *st, u32 *survbuf,
                     int want_e, int want_k, float *zb_out) {
    double t0 = now();
    u32 DK7 = dk7_22(t);
    u8 u[1024]; int ci[2][1024]; float Ng[256] = {0}; u64 Cw[2][W] = {{0}};
    for (int i = 0; i < n; i++) {
        u32 R7 = CLa[i] ^ CRa[i], L7 = CLa[i] ^ feal_N(feal_M(R7) ^ DK7);
        for (int f = 0; f < 2; f++) {
            ci[f][i] = par(A8_Y6[f] & R7) ^ par(0x00000404u & L7);
            if (ci[f][i]) Cw[f][i >> 6] |= 1ULL << (i & 63);
        }
        u[i] = (u8)part_idx(BY(L7, 0) ^ BY(L7, 1), BY(L7, 2) ^ BY(L7, 3));
        Ng[u[i]] += 1;
    }
    float isd[2][256];
    wht(Ng);
    for (int f = 0; f < 2; f++) {
        float tt[256];
        for (int g = 0; g < 256; g++) tt[g] = Ng[g] * E2hat[f][g];
        wht(tt);
        for (int k = 0; k < 256; k++) isd[f][k] = 1.f / (sqrtf(tt[k] / 256.f) * 256.f);
    }
    const float r1sq = (float)(r1 * r1);
    int ns = 0;
    for (int e = 0; e < 32768; e++) {
        float T[2][256] = {{0}};
        for (int f = 0; f < 2; f++) {
            const u64 *pv = Pv[f][e];
            for (int w = 0; w < W; w++) {
                u64 sx = pv[w] ^ Cw[f][w];
                int base = w * 64, lim = n - base < 64 ? n - base : 64;
                for (int b = 0; b < lim; b++) T[f][u[base + b]] += 1.f - 2.f * (float)((sx >> b) & 1);
            }
            wht(T[f]);
            for (int g = 0; g < 256; g++) T[f][g] *= Ehat[f][g];
            wht(T[f]);
        }
        if (e == want_e && zb_out) { zb_out[0] = T[0][want_k] * isd[0][want_k]; zb_out[1] = T[1][want_k] * isd[1][want_k]; }
        for (int k = 0; k < 256; k++) {
            float a = T[0][k] * isd[0][k], b = T[1][k] * isd[1][k];
            if (a * a + b * b >= r1sq) survbuf[ns++] = (u32)e | ((u32)k << 15);
        }
    }
    double t1 = now();
    st->t1 = t1 - t0; st->surv = ns;
    int cnt[257] = {0};
    for (int s = 0; s < ns; s++) cnt[(survbuf[s] >> 15) + 1]++;
    for (int k = 0; k < 256; k++) cnt[k + 1] += cnt[k];
    u16 *byk = malloc(sizeof(u16) * (ns + 1));
    { int pos[256]; memcpy(pos, cnt, sizeof pos);
      for (int s = 0; s < ns; s++) byk[pos[survbuf[s] >> 15]++] = (u16)(survbuf[s] & 0x7FFF); }
    int nc = 0; long long cpre = 0, cfull = 0;
    static __thread u64 V[2][128][W] __attribute__((aligned(64)));
    u64 blk[2][64];
    u32 a[1024];
    const float isq = 1.f / sqrtf((float)n);
    const int prelim = (int)ceil(rpre * rpre * 256.0);   /* z4_0^2+z4_1^2 >= rpre^2 * 256（z4 は ±256 の整数） */
    const float r2sq = (float)(r2 * r2);
    for (int e8 = 0; e8 < 256; e8++) {
        u32 DKf = DK7 | dk7_ext(e8);
        for (int i = 0; i < n; i++) {
            u32 R7 = CLa[i] ^ CRa[i], L7 = CLa[i] ^ feal_N(feal_M(R7) ^ DKf);
            a[i] = ((u32)(BY(L7, 0) ^ BY(L7, 1)) << 8) | (BY(L7, 2) ^ BY(L7, 3));
        }
        int delta = partx(a[0]) ^ u[0];
        for (int k = 0; k < 256; k++) {
            int lo = cnt[k], hi = cnt[k + 1];
            if (lo == hi) continue;
            u32 base = emb_part(k ^ delta);
            for (int f = 0; f < 2; f++)
                for (int b = 0; b < W; b++) {
                    for (int j = 0; j < 64; j++) {
                        int i = b * 64 + j;
                        if (i < n) {
                            u64 m = ci[f][i] ? ~0ULL : 0ULL;
                            const u64 *rw = ROW[f][a[i] ^ base];
                            blk[0][j] = rw[0] ^ m; blk[1][j] = rw[1] ^ m;
                        } else { blk[0][j] = 0; blk[1][j] = 0; }
                    }
                    transpose64(blk[0]); transpose64(blk[1]);
                    for (int r = 0; r < 64; r++) { V[f][r][b] = blk[0][r]; V[f][64 + r][b] = blk[1][r]; }
                }
            for (int s = lo; s < hi; s++) {
                int e = byk[s];
                const u64 *p0 = Pv[0][e], *p1 = Pv[1][e];
                for (int r = 0; r < 128; r++) {
                    const u64 *v0 = V[0][r], *v1 = V[1][r];
                    int d0 = __builtin_popcountll(p0[0] ^ v0[0]) + __builtin_popcountll(p0[1] ^ v0[1]) +
                             __builtin_popcountll(p0[2] ^ v0[2]) + __builtin_popcountll(p0[3] ^ v0[3]);
                    int d1 = __builtin_popcountll(p1[0] ^ v1[0]) + __builtin_popcountll(p1[1] ^ v1[1]) +
                             __builtin_popcountll(p1[2] ^ v1[2]) + __builtin_popcountll(p1[3] ^ v1[3]);
                    int z40 = 256 - 2 * d0, z41 = 256 - 2 * d1;
                    cpre++;
                    if (z40 * z40 + z41 * z41 < prelim) continue;
                    cfull++;
                    for (int w = 4; w < W; w++) { d0 += __builtin_popcountll(p0[w] ^ v0[w]); d1 += __builtin_popcountll(p1[w] ^ v1[w]); }
                    float z0 = (float)(n - 2 * d0) * isq, z1 = (float)(n - 2 * d1) * isq;
                    float rr = z0 * z0 + z1 * z1;
                    if (rr >= r2sq && nc < maxc) {
                        cand[nc].t = t; cand[nc].e = (u16)e; cand[nc].kap = (u8)k; cand[nc].e8 = (u8)e8; cand[nc].r = (u8)r;
                        cand[nc].z = sqrtf(rr); cand[nc].z0 = z0; cand[nc].z1 = z1; nc++;
                    }
                }
            }
        }
    }
    free(byk);
    st->t2 = now() - t1; st->cmp_pre = cpre; st->cmp_full = cfull; st->ncand = nc;
    return nc;
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
        seen[n] = p;
        PLa[n] = p >> 32; PRa[n] = (u32)p; CLa[n] = c >> 32; CRa[n] = (u32)c; n++;
    }
}

static int cmpz(const void *a, const void *b) {
    float x = ((const cand_t *)a)->z, y = ((const cand_t *)b)->z;
    return x < y ? 1 : x > y ? -1 : 0;
}

int main(int argc, char **argv) {
    if (argc < 2) return 1;
    int self = !strcmp(argv[1], "self");
    int ai = self ? 5 : 7;
    if (argc > ai) r1 = atof(argv[ai]);
    if (argc > ai + 1) rpre = atof(argv[ai + 1]);
    if (argc > ai + 2) r2 = atof(argv[ai + 2]);
    double T0 = now();
    setup_tables();
    u64 key = 0; u32 tlist_n; u32 *tlist;
    u32 t_true = 0; int e_true = -1, k_true = -1, e8_true = -1, r_true = -1;
    if (self) {
        seed(strtoull(argv[3], 0, 0));
        key = rnd();
        load_pairs(argv[2], NULL, 1, key);
        sched_t s; key_schedule(key, &s);
        u32 K89 = ((u32)s.k[8] << 16) | s.k[9], Kab = ((u32)s.k[10] << 16) | s.k[11];
        u32 Kcd = ((u32)s.k[12] << 16) | s.k[13], Kef = ((u32)s.k[14] << 16) | s.k[15];
        u32 EK0 = feal_M(K89 ^ Kab) ^ ((u32)s.k[0] << 8);
        u32 DK7 = feal_M(Kcd ^ Kef) ^ ((u32)s.k[7] << 8);
        u32 DK6 = feal_M(Kcd) ^ ((u32)s.k[6] << 8);
        int e1 = BY(EK0, 1), e2 = BY(EK0, 2);
        if (e2 & 0x80) { e1 ^= 0x80; e2 ^= 0x80; }
        e_true = e1 | (e2 << 8);
        t_true = BY(DK7, 1) | (BY(DK7, 2) << 8) | ((BY(DK7, 0) & 3) << 16) | ((BY(DK7, 3) & 15) << 18);
        e8_true = ((BY(DK7, 0) >> 2) & 31) | (((BY(DK7, 3) >> 4) & 7) << 5);
        u32 DKf = dk7_22(t_true) | dk7_ext(e8_true), Kx = 0; int u0 = 0;
        {
            u32 R7 = CLa[0] ^ CRa[0];
            u32 L7t = CLa[0] ^ feal_N(feal_M(R7) ^ DK7);
            u32 w = feal_M(L7t) ^ DK6;
            u32 L7e = CLa[0] ^ feal_N(feal_M(R7) ^ DKf);
            u32 ax = ((u32)(BY(L7e, 0) ^ BY(L7e, 1)) << 8) | (BY(L7e, 2) ^ BY(L7e, 3));
            Kx = (((u32)BY(w, 1) << 8) | BY(w, 2)) ^ ax;
            u32 L7s = CLa[0] ^ feal_N(feal_M(R7) ^ dk7_22(t_true));
            u0 = part_idx(BY(L7s, 0) ^ BY(L7s, 1), BY(L7s, 2) ^ BY(L7s, 3));
            if (Kx & 0x80) Kx ^= 0x8080;
            int kp = partx(Kx);
            for (int r = 0; r < 128; r++) if ((emb_rest(r) | emb_part(kp)) == Kx) r_true = r;
            k_true = kp ^ partx(ax) ^ u0;
        }
        /* 参照値（真の中間値） */
        double sA[2] = {0, 0}, Sg[2][256] = {{0}}; int Ngi[256] = {0};
        for (int i = 0; i < n; i++) {
            u32 st[18]; feal8_enc(&s, ((u64)PLa[i] << 32) | PRa[i], st);
            u32 w = feal_M(st[13]) ^ ((u32)s.k[6] << 8);
            int g = part_idx(BY(w, 1), BY(w, 2)); Ngi[g]++;
            for (int f = 0; f < 2; f++) {
                int q = par(A8_L1 & st[2]) ^ par(a8_R1(f) & st[3]);
                sA[f] += (q ^ par(A8_Y6[f] & st[12]) ^ par(0x00010504u & st[13])) ? -1 : 1;
                int qp = q ^ par(A8_Y6[f] & st[15]) ^ par(0x00010504u & st[14]) ^ ((BY(w, 0) ^ BY(w, 1)) & 1) ^ ((BY(w, 3) ^ BY(w, 2)) & 1);
                Sg[f][g] += qp ? -1 : 1;
            }
        }
        double zB[2];
        for (int f = 0; f < 2; f++) {
            double num = 0, den = 0;
            for (int g = 0; g < 256; g++) { num += etab[f][g] * Sg[f][g]; den += etab[f][g] * etab[f][g] * Ngi[g]; }
            zB[f] = num / sqrt(den);
        }
        printf("key %016llx N=%d truth: t=%06x e=%04x kappa=%02x e8=%02x r=%02x\n",
               (unsigned long long)key, n, t_true, e_true, k_true, e8_true, r_true);
        printf("reference (true intermediates): zA0=%.3f zA1=%.3f zB0=%.3f zB1=%.3f\n",
               sA[0] / sqrt(n), sA[1] / sqrt(n), zB[0], zB[1]);
        int nd = atoi(argv[4]);
        tlist_n = nd; tlist = malloc(sizeof(u32) * nd);
        tlist[0] = t_true;
        for (int i = 1; i < nd; i++) { u32 t; do t = (u32)(rnd() & 0x3FFFFF); while (t == t_true); tlist[i] = t; }
    } else {
        load_pairs(argv[2], argv[3], 0, 0);
        u32 start = strtoul(argv[4], 0, 0), count = strtoul(argv[5], 0, 0);
        tlist_n = count; tlist = malloc(sizeof(u32) * count);
        for (u32 i = 0; i < count; i++) tlist[i] = start + i;
    }
    build_plain_side();
    double T1 = now();
    printf("setup %.2fs, threads %d, r1=%.4f rpre=%.3f r2=%.3f\n", T1 - T0, omp_get_max_threads(), r1, rpre, r2);
    fflush(stdout);

    const int MAXC = 1 << 18;
    cand_t *all = malloc(sizeof(cand_t) * MAXC * 4); int nall = 0;
    double sum_t1 = 0, sum_t2 = 0; long long sum_surv = 0, sum_pre = 0, sum_full = 0, sum_cand = 0;
    float zb_true[2] = {0, 0}; int true_surv = 0;
    double W0 = now();
    #pragma omp parallel
    {
        u32 *sb = malloc(sizeof(u32) * (1 << 23));
        cand_t *cb = malloc(sizeof(cand_t) * MAXC);
        #pragma omp for schedule(dynamic, 1) reduction(+:sum_t1,sum_t2,sum_surv,sum_pre,sum_full,sum_cand)
        for (u32 j = 0; j < tlist_n; j++) {
            dstat_t st; float zb[2] = {0, 0};
            int want = self && j == 0;
            int nc = process_d(tlist[j], cb, MAXC, &st, sb, want ? e_true : -1, want ? k_true : -1, want ? zb : NULL);
            sum_t1 += st.t1; sum_t2 += st.t2; sum_surv += st.surv; sum_pre += st.cmp_pre; sum_full += st.cmp_full; sum_cand += nc;
            if (want) {
                zb_true[0] = zb[0]; zb_true[1] = zb[1];
                for (int s = 0; s < st.surv; s++) if (sb[s] == ((u32)e_true | ((u32)k_true << 15))) true_surv = 1;
            }
            #pragma omp critical
            { for (int c = 0; c < nc && nall < MAXC * 4; c++) all[nall++] = cb[c]; }
        }
        free(sb); free(cb);
    }
    double W1 = now();
    printf("processed %u DK7 values in %.1fs wall (%.3f s/DK7 wall)\n", tlist_n, W1 - W0, (W1 - W0) / tlist_n);
    printf("per DK7 (thread time): stage1 %.3fs, stage2 %.3fs; survivors %.0f (expect %.0f); pre %.3g, full %.3g (%.2f%%); cands(r>=%.2f) %.2f\n",
           sum_t1 / tlist_n, sum_t2 / tlist_n, (double)sum_surv / tlist_n, 8388608.0 * exp(-r1 * r1 / 2),
           (double)sum_pre / tlist_n, (double)sum_full / tlist_n, 100.0 * sum_full / sum_pre, r2, (double)sum_cand / tlist_n);
    printf("extrapolated full 2^22 DK7: %.1f core-days (thread time)\n", (sum_t1 + sum_t2) / tlist_n * 4194304.0 / 86400);
    qsort(all, nall, sizeof(cand_t), cmpz);
    if (self) {
        printf("true hypothesis: stage1 z0=%.3f z1=%.3f r=%.3f survived=%s\n", zb_true[0], zb_true[1],
               hypot(zb_true[0], zb_true[1]), true_surv ? "yes" : "no");
        int rank = -1;
        for (int c = 0; c < nall; c++)
            if (all[c].t == t_true && all[c].e == e_true && all[c].kap == k_true && all[c].e8 == e8_true && all[c].r == r_true) { rank = c; break; }
        if (rank >= 0) printf("true key in candidates: z0=%.3f z1=%.3f r=%.3f\n", all[rank].z0, all[rank].z1, all[rank].z);
        else printf("true key in candidates: no\n");
        const double th[] = {6.0, 6.5, 7.0, 7.5, 8.0, 8.5};
        long long cw[6] = {0}, ct[6] = {0};
        for (int c = 0; c < nall; c++)
            for (int j = 0; j < 6; j++) { if (all[c].z >= th[j]) { if (all[c].t == t_true) ct[j]++; else cw[j]++; } }
        double sc = tlist_n > 1 ? 4194304.0 / (tlist_n - 1) : 0;
        printf("candidates by r threshold (wrong DK7 -> full-space estimate | from the true DK7 row):\n");
        for (int j = 0; j < 6; j++) printf("  r>=%.1f: wrong %lld -> %.3g | true-row %lld\n", th[j], cw[j], cw[j] * sc, ct[j]);
        printf("top candidates:\n");
        for (int c = 0; c < nall && c < 6; c++)
            printf("  t=%06x e=%04x kap=%02x e8=%02x r=%02x z0=%.3f z1=%.3f r=%.3f%s\n", all[c].t, all[c].e, all[c].kap, all[c].e8, all[c].r,
                   all[c].z0, all[c].z1, all[c].z,
                   (all[c].t == t_true && all[c].e == e_true && all[c].kap == k_true && all[c].e8 == e8_true && all[c].r == r_true) ? "  <== TRUE" : "");
    } else {
        FILE *fo = fopen(argv[6], "w");
        for (int c = 0; c < nall; c++)
            fprintf(fo, "%06x %04x %02x %02x %02x %.4f %.4f %.4f\n", all[c].t, all[c].e, all[c].kap, all[c].e8, all[c].r, all[c].z, all[c].z0, all[c].z1);
        fclose(fo);
        printf("wrote %d candidates\n", nall);
    }
    return 0;
}
