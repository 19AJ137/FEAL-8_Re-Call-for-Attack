/* 方式C：8ビット分割の統計量でふるい（段階1）、生存仮説を補完して方式Aの統計量で確定（段階2）
 *
 *   ./attackC self  plaintext.bin seed nd [thr1 thr2pre thrrec]
 *       乱数鍵で暗号文を作り、真の DK7(22bit) と乱数の DK7 を nd-1 個処理する。正解の通過と時間を報告
 *   ./attackC real  plaintext.bin ciphertext.bin start count out.txt [thr1 thr2pre thrrec]
 *       DK7(22bit) の番号 start..start+count-1 を処理し、候補を out.txt に書く
 *
 * 記号（等価表現）
 *   平文側  L1 = PL^PR, R1 = PL ^ N(M(L1)^EK0)、p = <00010504,L1>^<04031004,R1>
 *           EK0 はバイト1・2の 15 ビット（e2 の最上位は 0 に正規化）
 *   暗号文側 R7 = CL^CR, L7 = CL ^ N(M(R7)^DK7)
 *           段階1：DK7 は 03FFFF0F の 22 ビット。c = <04031004,R7>^<00000404,L7>、
 *                 群 u = part((L7_0^L7_1), (L7_2^L7_3))、κ は DK6 の対応する 8 ビット（XOR）
 *           段階2：DK7 の残り 8 ビット（バイト0のビット2〜6、バイト3のビット4〜6）と
 *                 DK6 の残り 7 ビット（w1 のビット0,1,5,6,7、w2 のビット0,1）を補う。
 *                 方式A の暗号文側ビット = c ^ h'(w1,w2)、(w1,w2) = (L7_0^L7_1, L7_2^L7_3) ^ k
 */
#include "feal8.h"
#include <omp.h>
#include <sys/time.h>

#define W 16
static int n;
static u32 PLa[1024], PRa[1024], CLa[1024], CRa[1024];
static u64 (*Pv)[W];        /* 平文側ビット列 2^15 本 */
static float Ehat[256];     /* WHT(e) */
static float E2hat[256];    /* WHT(e^2) */
static u64 ROW[65536][2];   /* ROW[x] のビット r = h'(x ^ emb(r)) */
static u8 HP[65536];        /* h'(x)、x = (w1<<8)|w2 */
static double thr1 = 3.0, thr2pre = 2.0, thrrec = 6.0;

static double now(void) { struct timeval tv; gettimeofday(&tv, 0); return tv.tv_sec + tv.tv_usec * 1e-6; }

static inline int hprime(int w1, int w2) {
    return par(G6_Y & feal_N(PK(0, (u8)w1, (u8)w2, 0))) ^ (w1 & 1) ^ (w2 & 1);
}
/* 16 ビット x = (w1<<8)|w2 の分割ビット（8 ビット） */
static inline int partx(u32 x) { return part_idx((u8)(x >> 8), (u8)x); }
/* 分割ビット kp を x の位置へ */
static inline u32 emb_part(int kp) { return ((u32)(kp & 7) << 10) | ((u32)((kp >> 3) & 31) << 2); }
/* 残り 7 ビット r を x の位置へ：w1 のビット 0,1,5,6,7 と w2 のビット 0,1 */
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
            for (int j = i; j < i + h; j++) {
                float x = a[j], y = a[j + h];
                a[j] = x + y; a[j + h] = x - y;
            }
}

static void setup_tables(void) {
    build_etab();
    float e2[256];
    for (int g = 0; g < 256; g++) { Ehat[g] = (float)e_tab[g]; e2[g] = (float)(e_tab[g] * e_tab[g]); }
    wht(Ehat);
    memcpy(E2hat, e2, sizeof e2); wht(E2hat);
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

static void build_plain_side(void) {
    Pv = aligned_alloc(64, sizeof(*Pv) * 32768);
    #pragma omp parallel for schedule(static)
    for (int e = 0; e < 32768; e++) {
        u32 ek = PK(0, (u8)e, (u8)(e >> 8), 0);
        u64 v[W] = {0};
        for (int i = 0; i < n; i++) {
            u32 L1 = PLa[i] ^ PRa[i], R1 = PLa[i] ^ feal_N(feal_M(L1) ^ ek);
            if (par(G_L1 & L1) ^ par(G_R1 & R1)) v[i >> 6] |= 1ULL << (i & 63);
        }
        memcpy(Pv[e], v, sizeof v);
    }
}

/* 64x64 ビット行列の転置（a[i] のビット j → 結果 a[j] のビット i） */
static void transpose64(u64 *a) {
    int j = 32; u64 m = 0x00000000FFFFFFFFULL;
    for (; j != 0; j >>= 1, m ^= (m << j)) {
        for (int k = 0; k < 64; k = ((k | j) + 1) & ~j) {
            u64 t = ((a[k] >> j) ^ a[k | j]) & m;
            a[k] ^= t << j; a[k | j] ^= t;
        }
    }
}

typedef struct { u32 t; u16 e; u8 kap, e8, r; float z; } cand_t;
typedef struct {
    double t1, t2; long long surv, cmp_pre, cmp_full; int ncand;
} dstat_t;

/* 1 つの DK7(22bit) を処理。truth が非 NULL なら真の (e,κ) の段階1スコアを返す */
static int process_d(u32 t, cand_t *cand, int maxc, dstat_t *st, u32 *survbuf,
                     int want_e, int want_k, float *zb_out) {
    double t0 = now();
    u32 DK7 = dk7_22(t);
    u8 u[1024]; u64 C[W] = {0}; float Ng[256] = {0};
    int ci[1024];
    for (int i = 0; i < n; i++) {
        u32 R7 = CLa[i] ^ CRa[i], L7 = CLa[i] ^ feal_N(feal_M(R7) ^ DK7);
        int c = par(G_L6 & R7) ^ par(0x00000404u & L7);
        ci[i] = c;
        if (c) C[i >> 6] |= 1ULL << (i & 63);
        u[i] = (u8)part_idx(BY(L7, 0) ^ BY(L7, 1), BY(L7, 2) ^ BY(L7, 3));
        Ng[u[i]] += 1;
    }
    /* D^2(κ)*256 = WHT(WHT(N)·WHT(e^2)) */
    float thr[256];
    wht(Ng);
    for (int g = 0; g < 256; g++) Ng[g] *= E2hat[g];
    wht(Ng);
    for (int k = 0; k < 256; k++) thr[k] = (float)(thr1 * sqrt(Ng[k] / 256.0) * 256.0);
    /* 段階1 */
    int ns = 0;
    for (int e = 0; e < 32768; e++) {
        float T[256] = {0};
        const u64 *pv = Pv[e];
        for (int w = 0; w < W; w++) {
            u64 s = pv[w] ^ C[w];
            int base = w * 64, lim = n - base < 64 ? n - base : 64;
            for (int b = 0; b < lim; b++) T[u[base + b]] += 1.f - 2.f * (float)((s >> b) & 1);
        }
        wht(T);
        for (int g = 0; g < 256; g++) T[g] *= Ehat[g];
        wht(T);
        if (e == want_e && zb_out) *zb_out = T[want_k] / (thr[want_k] / (float)thr1);
        for (int k = 0; k < 256; k++)
            if (fabsf(T[k]) >= thr[k]) survbuf[ns++] = (u32)e | ((u32)k << 15);
    }
    double t1 = now();
    st->t1 = t1 - t0; st->surv = ns;
    /* 段階2：κ ごとにまとめる */
    int cnt[257] = {0};
    for (int s = 0; s < ns; s++) cnt[(survbuf[s] >> 15) + 1]++;
    for (int k = 0; k < 256; k++) cnt[k + 1] += cnt[k];
    u16 *byk = malloc(sizeof(u16) * (ns + 1));
    { int pos[256]; memcpy(pos, cnt, sizeof pos);
      for (int s = 0; s < ns; s++) byk[pos[survbuf[s] >> 15]++] = (u16)(survbuf[s] & 0x7FFF); }
    int nc = 0; long long cpre = 0, cfull = 0;
    u64 V[128][W] __attribute__((aligned(64)));
    u64 blk[2][64];
    u32 a[1024];
    const float isq = 1.f / sqrtf((float)n);
    const int pre_lim = (int)ceil(thr2pre * 16.0); /* |256-2d| >= thr2pre*sqrt(256) */
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
            /* 128 本の暗号文側ビット列を転置で作る */
            for (int b = 0; b < W; b++) {
                for (int j = 0; j < 64; j++) {
                    int i = b * 64 + j;
                    if (i < n) {
                        u64 m = ci[i] ? ~0ULL : 0ULL;
                        const u64 *rw = ROW[a[i] ^ base];
                        blk[0][j] = rw[0] ^ m; blk[1][j] = rw[1] ^ m;
                    } else { blk[0][j] = 0; blk[1][j] = 0; }
                }
                transpose64(blk[0]); transpose64(blk[1]);
                for (int r = 0; r < 64; r++) { V[r][b] = blk[0][r]; V[64 + r][b] = blk[1][r]; }
            }
            for (int s = lo; s < hi; s++) {
                int e = byk[s];
                const u64 *pv = Pv[e];
                for (int r = 0; r < 128; r++) {
                    const u64 *v = V[r];
                    int d4 = __builtin_popcountll(pv[0] ^ v[0]) + __builtin_popcountll(pv[1] ^ v[1]) +
                             __builtin_popcountll(pv[2] ^ v[2]) + __builtin_popcountll(pv[3] ^ v[3]);
                    int z4 = 256 - 2 * d4;
                    cpre++;
                    if (z4 < pre_lim && z4 > -pre_lim) continue;
                    cfull++;
                    int d = d4;
                    for (int w = 4; w < W; w++) d += __builtin_popcountll(pv[w] ^ v[w]);
                    float z = (float)(n - 2 * d) * isq;
                    if (fabsf(z) >= thrrec && nc < maxc) {
                        cand[nc].t = t; cand[nc].e = (u16)e; cand[nc].kap = (u8)k;
                        cand[nc].e8 = (u8)e8; cand[nc].r = (u8)r; cand[nc].z = z; nc++;
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
    FILE *fp = fopen(pf, "rb"); u8 pb[8192]; fread(pb, 1, 8192, fp); fclose(fp);
    u8 cb[8192];
    if (!self) { fp = fopen(cf, "rb"); fread(cb, 1, 8192, fp); fclose(fp); }
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
    float x = fabsf(((const cand_t *)a)->z), y = fabsf(((const cand_t *)b)->z);
    return x < y ? 1 : x > y ? -1 : 0;
}

int main(int argc, char **argv) {
    if (argc < 2) return 1;
    int self = !strcmp(argv[1], "self");
    int ai = self ? 5 : 7;
    if (argc > ai) thr1 = atof(argv[ai]);
    if (argc > ai + 1) thr2pre = atof(argv[ai + 1]);
    if (argc > ai + 2) thrrec = atof(argv[ai + 2]);
    double T0 = now();
    setup_tables();
    u64 key = 0;
    u32 tlist_n; u32 *tlist;
    u32 t_true = 0; int e_true = -1, k_true = -1, e8_true = -1, r_true = -1;
    if (self) {
        seed(strtoull(argv[3], 0, 0));
        key = rnd();
        load_pairs(argv[2], NULL, 1, key);
        /* 真値 */
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
        /* k の実効値：真の w と、列挙した DK7 から作る a の差 */
        u32 DKf = dk7_22(t_true) | dk7_ext(e8_true), Kx = 0; int ok = 1, u0 = 0;
        for (int i = 0; i < n; i++) {
            u32 R7 = CLa[i] ^ CRa[i];
            u32 L7t = CLa[i] ^ feal_N(feal_M(R7) ^ DK7);
            u32 w = feal_M(L7t) ^ DK6;
            u32 L7e = CLa[i] ^ feal_N(feal_M(R7) ^ DKf);
            u32 ax = ((u32)(BY(L7e, 0) ^ BY(L7e, 1)) << 8) | (BY(L7e, 2) ^ BY(L7e, 3));
            u32 wx = ((u32)BY(w, 1) << 8) | BY(w, 2);
            /* 真の w の分割ビットと残りビットのうち、h' に効く部分だけ一致すればよい */
            u32 kx = (wx ^ ax);
            if (i == 0) {
                Kx = kx;
                u32 L7s = CLa[i] ^ feal_N(feal_M(R7) ^ dk7_22(t_true));
                u0 = part_idx(BY(L7s, 0) ^ BY(L7s, 1), BY(L7s, 2) ^ BY(L7s, 3));
            } else if ((kx ^ Kx) & 0x7F7F) ok = 0; /* ビット7・15 は XOR だけが効く */
            else if (((kx ^ Kx) >> 15 ^ (kx ^ Kx) >> 7) & 1) ok = 0;
        }
        if (Kx & 0x80) Kx ^= 0x8080;
        int kp = partx(Kx);
        r_true = 0;
        for (int r = 0; r < 128; r++) if ((emb_rest(r) | emb_part(kp)) == Kx) r_true = r;
        if ((emb_rest(r_true) | emb_part(kp)) != Kx) ok = 0;
        u32 a0; { u32 R7 = CLa[0] ^ CRa[0], L7e = CLa[0] ^ feal_N(feal_M(R7) ^ DKf);
                  a0 = ((u32)(BY(L7e, 0) ^ BY(L7e, 1)) << 8) | (BY(L7e, 2) ^ BY(L7e, 3)); }
        k_true = kp ^ partx(a0) ^ u0;
        /* corr.c と同じ定義での正解鍵スコア */
        long long sA = 0; double num = 0, den = 0, Sg[256] = {0}; int Ngi[256] = {0};
        for (int i = 0; i < n; i++) {
            u32 st[18]; feal8_enc(&s, ((u64)PLa[i] << 32) | PRa[i], st);
            int q = par(G_L1 & st[2]) ^ par(G_R1 & st[3]);
            sA += (q ^ par(G_L6 & st[12]) ^ par(G_R6 & st[13])) ? -1 : 1;
            u32 w = feal_M(st[13]) ^ ((u32)s.k[6] << 8);
            int qp = q ^ par(G_L6 & st[15]) ^ par(G_R6 & st[14]) ^ ((BY(w, 0) ^ BY(w, 1)) & 1) ^ ((BY(w, 3) ^ BY(w, 2)) & 1);
            int g = part_idx(BY(w, 1), BY(w, 2));
            Sg[g] += qp ? -1 : 1; Ngi[g]++;
        }
        for (int g = 0; g < 256; g++) { num += e_tab[g] * Sg[g]; den += e_tab[g] * e_tab[g] * Ngi[g]; }
        printf("key %016llx N=%d  truth: t=%06x e=%04x kappa=%02x e8=%02x r=%02x consistent=%s\n",
               (unsigned long long)key, n, t_true, e_true, k_true, e8_true, r_true, ok ? "yes" : "NO");
        printf("reference scores (true intermediates): zA=%.3f zB=%.3f\n", (double)sA / sqrt(n), num / sqrt(den));
        int nd = atoi(argv[4]);
        tlist_n = nd; tlist = malloc(sizeof(u32) * nd);
        tlist[0] = t_true;
        for (int i = 1; i < nd; i++) tlist[i] = (u32)(rnd() & 0x3FFFFF);
    } else {
        load_pairs(argv[2], argv[3], 0, 0);
        u32 start = strtoul(argv[4], 0, 0), count = strtoul(argv[5], 0, 0);
        tlist_n = count; tlist = malloc(sizeof(u32) * count);
        for (u32 i = 0; i < count; i++) tlist[i] = start + i;
    }
    build_plain_side();
    double T1 = now();
    printf("setup %.2fs, threads %d, thr1=%.2f thr2pre=%.2f thrrec=%.2f\n", T1 - T0, omp_get_max_threads(), thr1, thr2pre, thrrec);
    fflush(stdout);

    const int MAXC = 1 << 16;
    cand_t *all = malloc(sizeof(cand_t) * MAXC * 4); int nall = 0;
    double sum_t1 = 0, sum_t2 = 0; long long sum_surv = 0, sum_pre = 0, sum_full = 0, sum_cand = 0;
    float zb_true = 0; int true_surv = 0;
    double W0 = now();
    #pragma omp parallel
    {
        u32 *sb = malloc(sizeof(u32) * (1 << 23));
        cand_t *cb = malloc(sizeof(cand_t) * MAXC);
        #pragma omp for schedule(dynamic, 1) reduction(+:sum_t1,sum_t2,sum_surv,sum_pre,sum_full,sum_cand)
        for (u32 j = 0; j < tlist_n; j++) {
            dstat_t st; float zb = 0;
            int want = self && j == 0;
            int nc = process_d(tlist[j], cb, MAXC, &st, sb, want ? e_true : -1, want ? k_true : -1, want ? &zb : NULL);
            sum_t1 += st.t1; sum_t2 += st.t2; sum_surv += st.surv; sum_pre += st.cmp_pre; sum_full += st.cmp_full; sum_cand += nc;
            if (want) {
                zb_true = zb;
                for (int s = 0; s < st.surv; s++) if (sb[s] == ((u32)e_true | ((u32)k_true << 15))) true_surv = 1;
            }
            #pragma omp critical
            {
                for (int c = 0; c < nc && nall < MAXC * 4; c++) all[nall++] = cb[c];
            }
        }
        free(sb); free(cb);
    }
    double W1 = now();
    printf("processed %u DK7 values in %.1fs wall (%.3f s/DK7 wall)\n", tlist_n, W1 - W0, (W1 - W0) / tlist_n);
    printf("per DK7 (thread time): stage1 %.3fs, stage2 %.3fs; survivors %.0f (expect %.0f); "
           "pre-compares %.3g, full %.3g (%.2f%%); candidates(|z|>=%.1f) %.2f\n",
           sum_t1 / tlist_n, sum_t2 / tlist_n, (double)sum_surv / tlist_n, 8388608.0 * erfc(thr1 / sqrt(2)),
           (double)sum_pre / tlist_n, (double)sum_full / tlist_n, 100.0 * sum_full / sum_pre, thrrec, (double)sum_cand / tlist_n);
    double full_wall = (W1 - W0) / tlist_n * 4194304.0;
    printf("extrapolated full 2^22 DK7: %.1f h wall on this machine (%d threads), %.1f core-days (thread time)\n",
           full_wall / 3600, omp_get_max_threads(), (sum_t1 + sum_t2) / tlist_n * 4194304.0 / 86400);
    qsort(all, nall, sizeof(cand_t), cmpz);
    if (self && getenv("CANDOUT")) {
        FILE *fo = fopen(getenv("CANDOUT"), "w");
        for (int c = 0; c < nall; c++)
            fprintf(fo, "%06x %04x %02x %02x %02x %.4f %d\n", all[c].t, all[c].e, all[c].kap, all[c].e8, all[c].r, all[c].z,
                    all[c].t == t_true ? 1 : 0);
        fclose(fo);
    }
    if (self) {
        printf("true hypothesis: stage1 zB=%.3f survived=%s\n", zb_true, true_surv ? "yes" : "no");
        int rank = -1; float zt = 0;
        for (int c = 0; c < nall; c++)
            if (all[c].t == t_true && all[c].e == e_true && all[c].kap == k_true && all[c].e8 == e8_true && all[c].r == r_true) { rank = c; zt = all[c].z; }
        int above = 0;
        for (int c = 0; c < nall; c++) if (all[c].t != t_true && fabsf(all[c].z) >= fabsf(zt)) above++;
        printf("true key in candidates: %s zA=%.3f; wrong candidates with |z|>=|z_true| in processed set: %d\n",
               rank >= 0 ? "yes" : "no", zt, above);
        long long c65 = 0, c70 = 0, c75 = 0;
        for (int c = 0; c < nall; c++) {
            if (all[c].t == t_true) continue;
            float z = fabsf(all[c].z); c65 += z >= 6.5; c70 += z >= 7.0; c75 += z >= 7.5;
        }
        double sc = 4194304.0 / (tlist_n - 1);
        printf("wrong candidates from random DK7: >=6.5: %lld, >=7.0: %lld, >=7.5: %lld  -> full space x%.0f: %.3g, %.3g, %.3g\n",
               c65, c70, c75, sc, c65 * sc, c70 * sc, c75 * sc);
        printf("top candidates:\n");
        for (int c = 0; c < nall && c < 8; c++)
            printf("  t=%06x e=%04x kap=%02x e8=%02x r=%02x z=%.3f%s\n", all[c].t, all[c].e, all[c].kap, all[c].e8, all[c].r, all[c].z,
                   (all[c].t == t_true && all[c].e == e_true && all[c].kap == k_true && all[c].e8 == e8_true && all[c].r == r_true) ? "  <== TRUE" : "");
    } else {
        FILE *fo = fopen(argv[6], "w");
        for (int c = 0; c < nall; c++)
            fprintf(fo, "%06x %04x %02x %02x %02x %.4f\n", all[c].t, all[c].e, all[c].kap, all[c].e8, all[c].r, all[c].z);
        fclose(fo);
        printf("wrote %d candidates\n", nall);
    }
    return 0;
}
