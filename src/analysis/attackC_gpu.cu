// 方式C 段階1・2 の GPU 版（attackC.c と同じ統計量）
//   ./attackC_gpu self plaintext.bin seed nd [thr1 thrrec]       乱数鍵、真の DK7 + 乱数 DK7 (nd-1 個)
//   ./attackC_gpu real plaintext.bin ciphertext.bin start count out.txt [thr1 thrrec]
// 段階1：warp 1 本が (DK7, EK0) 1 組。群別符号和 T[256] を共有メモリに作り、WHT で κ 全通り
// 段階2：block 1 個が DK7 1 個。κ ごとに生存 EK0 のビット列を共有メモリへ、e8 ごとに L7 を計算し、
//        128 本の暗号文側ビット列を warp ballot で作って XOR＋popcount（先頭 256 組で 2σ 打ち切り）
#include "feal8.h"
#include <cuda_runtime.h>
#include <chrono>
#include <vector>
#include <algorithm>

#define NW 32           // 1024 ビット = 32 語
#define BD 2048         // 1 回に処理する DK7 の数
#define CAP 32768       // DK7 1 個あたりの生存仮説の上限（しきい値 3.0 で平均 22,600）
#define CCAP (1 << 22)  // 候補の上限（1 バッチ）
#define SMAX 192        // κ 1 個あたりの生存数の上限（超えた分は分割）。共有メモリ 48KB に収める

static double nowsec() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

__constant__ float cEhat[256], cE2hat[256];
__constant__ u32 cCL[1024], cCR[1024];
__constant__ int cN;

struct Cand { u32 t; u16 e; u8 kap, e8, r; float z; };

__device__ __forceinline__ u8 dS0(u8 x, u8 y) { u8 s = x + y; return (u8)((s << 2) | (s >> 6)); }
__device__ __forceinline__ u8 dS1(u8 x, u8 y) { u8 s = x + y + 1; return (u8)((s << 2) | (s >> 6)); }
__device__ __forceinline__ u8 dBY(u32 x, int i) { return (u8)(x >> (24 - 8 * i)); }
__device__ __forceinline__ u32 dPK(u8 a, u8 b, u8 c, u8 d) { return ((u32)a << 24) | ((u32)b << 16) | ((u32)c << 8) | d; }
__device__ __forceinline__ u32 dN(u32 x) {
    u8 x0 = dBY(x, 0), x1 = dBY(x, 1), x2 = dBY(x, 2), x3 = dBY(x, 3);
    u8 y1 = dS1(x1, x2), y2 = dS0(x2, y1);
    return dPK(dS0(x0, y1), y1, y2, dS1(x3, y2));
}
__device__ __forceinline__ u32 dM(u32 a) { u8 a0 = dBY(a, 0), a1 = dBY(a, 1), a2 = dBY(a, 2), a3 = dBY(a, 3); return dPK(a0, a0 ^ a1, a2 ^ a3, a3); }
__device__ __forceinline__ int dpart(u8 w1, u8 w2) { return ((w1 >> 2) & 7) | (((w2 >> 2) & 31) << 3); }
__device__ __forceinline__ u32 ddk7_22(u32 t) { return dPK((u8)((t >> 16) & 3), (u8)t, (u8)(t >> 8), (u8)((t >> 18) & 15)); }
__device__ __forceinline__ u32 ddk7_ext(int e8) { return ((u32)(e8 & 31) << 26) | ((u32)((e8 >> 5) & 7) << 4); }
__device__ __forceinline__ u32 demb_part(int kp) { return ((u32)(kp & 7) << 10) | ((u32)((kp >> 3) & 31) << 2); }

// ---- K0：DK7 ごとの前処理（群 u、暗号文側ビット c、しきい値） ----
__global__ void k_prep(const u32 *tl, int nd, float thr1, u8 *U, u32 *C, float *THR) {
    int d = blockIdx.x; if (d >= nd) return;
    u32 DK7 = ddk7_22(tl[d]);
    __shared__ float Ng[256];
    int i = threadIdx.x;
    if (i < 256) Ng[i] = 0;
    __syncthreads();
    int c = 0;
    if (i < cN) {
        u32 R7 = cCL[i] ^ cCR[i], L7 = cCL[i] ^ dN(dM(R7) ^ DK7);
        c = __popc(0x04031004u & R7) & 1; c ^= __popc(0x00000404u & L7) & 1;
        u8 u = (u8)dpart(dBY(L7, 0) ^ dBY(L7, 1), dBY(L7, 2) ^ dBY(L7, 3));
        U[(size_t)d * 1024 + i] = u;
        atomicAdd(&Ng[u], 1.f);
    }
    unsigned b = __ballot_sync(0xffffffffu, c);
    if ((i & 31) == 0) C[(size_t)d * NW + (i >> 5)] = b;
    __syncthreads();
    // WHT(N)·Ê2 → WHT
    for (int h = 1; h < 256; h <<= 1) {
        if (i < 128) { int j = (i / h) * 2 * h + (i % h); float x = Ng[j], y = Ng[j + h]; Ng[j] = x + y; Ng[j + h] = x - y; }
        __syncthreads();
    }
    if (i < 256) Ng[i] *= cE2hat[i];
    __syncthreads();
    for (int h = 1; h < 256; h <<= 1) {
        if (i < 128) { int j = (i / h) * 2 * h + (i % h); float x = Ng[j], y = Ng[j + h]; Ng[j] = x + y; Ng[j + h] = x - y; }
        __syncthreads();
    }
    if (i < 256) THR[(size_t)d * 256 + i] = thr1 * sqrtf(fmaxf(Ng[i], 0.f) * 256.f);
}

// ---- K1：段階1。block = (d, EK0 64 個)、warp 8 本が 8 個ずつ ----
__global__ void __launch_bounds__(256) k_stage1(const u32 *P, const u8 *U, const u32 *C, const float *THR, int nd,
                                                u32 *SURV, int *SCNT, int want_d, int want_e, float *zb_out) {
    int d = blockIdx.y; if (d >= nd) return;
    __shared__ u8 su[1024];
    __shared__ u32 sc[NW];
    __shared__ float sthr[256];
    __shared__ float T[8][256];
    int tid = threadIdx.x, w = tid >> 5, lane = tid & 31;
    for (int i = tid; i < 1024; i += 256) su[i] = U[(size_t)d * 1024 + i];
    if (tid < NW) sc[tid] = C[(size_t)d * NW + tid];
    sthr[tid] = THR[(size_t)d * 256 + tid];
    __syncthreads();
    float *t = T[w];
    for (int q = 0; q < 8; q++) {
        int e = blockIdx.x * 64 + w * 8 + q;
        for (int j = lane; j < 256; j += 32) t[j] = 0;
        __syncwarp();
        const u32 *pe = P + (size_t)e * NW;
        for (int ww = 0; ww < NW; ww++) {
            u32 s = pe[ww] ^ sc[ww];
            int i = ww * 32 + lane;
            if (i < cN) atomicAdd(&t[su[i]], ((s >> lane) & 1) ? -1.f : 1.f);
        }
        __syncwarp();
        for (int h = 1; h < 256; h <<= 1) {
            for (int k = lane; k < 128; k += 32) { int j = (k / h) * 2 * h + (k % h); float x = t[j], y = t[j + h]; t[j] = x + y; t[j + h] = x - y; }
            __syncwarp();
        }
        for (int j = lane; j < 256; j += 32) t[j] *= cEhat[j];
        __syncwarp();
        for (int h = 1; h < 256; h <<= 1) {
            for (int k = lane; k < 128; k += 32) { int j = (k / h) * 2 * h + (k % h); float x = t[j], y = t[j + h]; t[j] = x + y; t[j + h] = x - y; }
            __syncwarp();
        }
        for (int k = lane; k < 256; k += 32) {
            if (d == want_d && e == want_e && zb_out) zb_out[k] = t[k] / (sthr[k] > 0 ? sthr[k] : 1.f);
            if (fabsf(t[k]) >= sthr[k]) {
                int pos = atomicAdd(&SCNT[d], 1);
                if (pos < CAP) SURV[(size_t)d * CAP + pos] = (u32)e | ((u32)k << 15);
            }
        }
        __syncwarp();
    }
}

// ---- K2：DK7 ごとに κ で並べ替え ----
__global__ void k_sort(const u32 *SURV, const int *SCNT, int nd, u16 *SORTED, int *OFF) {
    int d = blockIdx.x; if (d >= nd) return;
    __shared__ int cnt[257], pos[256];
    int tid = threadIdx.x, n = min(SCNT[d], CAP);
    if (tid < 257) cnt[tid] = 0;
    __syncthreads();
    for (int s = tid; s < n; s += blockDim.x) atomicAdd(&cnt[(SURV[(size_t)d * CAP + s] >> 15) + 1], 1);
    __syncthreads();
    if (tid == 0) { for (int k = 0; k < 256; k++) cnt[k + 1] += cnt[k]; }
    __syncthreads();
    if (tid < 256) { pos[tid] = cnt[tid]; OFF[(size_t)d * 257 + tid] = cnt[tid]; }
    if (tid == 0) OFF[(size_t)d * 257 + 256] = cnt[256];
    __syncthreads();
    for (int s = tid; s < n; s += blockDim.x) {
        u32 v = SURV[(size_t)d * CAP + s];
        int p = atomicAdd(&pos[v >> 15], 1);
        SORTED[(size_t)d * CAP + p] = (u16)(v & 0x7FFF);
    }
}

// ---- K3：段階2。block = DK7 1 個、1024 スレッド（スレッド = 標本） ----
__device__ __forceinline__ u32 wtrans_step(u32 x, int j, u32 lo, int lane) {
    u32 y = __shfl_xor_sync(0xffffffffu, x, j);
    return (lane & j) ? ((x & ~lo) | ((y & ~lo) >> j)) : ((x & lo) | ((y & lo) << j));
}
// warp 内 32x32 ビット転置：lane l の語のビット c → lane c の語のビット l
__device__ __forceinline__ u32 wtrans(u32 x, int lane) {
    x = wtrans_step(x, 16, 0x0000FFFFu, lane);
    x = wtrans_step(x, 8, 0x00FF00FFu, lane);
    x = wtrans_step(x, 4, 0x0F0F0F0Fu, lane);
    x = wtrans_step(x, 2, 0x33333333u, lane);
    x = wtrans_step(x, 1, 0x55555555u, lane);
    return x;
}
__global__ void __launch_bounds__(1024) k_stage2(const u32 *tl, const u32 *P, const u32 *ROW, const u8 *U, const u32 *C,
                                                 const u16 *SORTED, const int *OFF, int nd, float thrrec, Cand *CB, int *CCNT,
                                                 unsigned long long *stats) {
    int d = blockIdx.x; if (d >= nd) return;
    __shared__ u32 V[128][NW + 1];
    __shared__ u32 Ps[SMAX][NW];
    __shared__ u16 Es[SMAX];
    __shared__ int soff[257];
    __shared__ int sdelta;
    int tid = threadIdx.x, w = tid >> 5, lane = tid & 31;
    u32 t = tl[d], DK7 = ddk7_22(t);
    int cbit = (C[(size_t)d * NW + w] >> lane) & 1;
    u32 cmask = cbit ? 0xffffffffu : 0u;
    int u0 = U[(size_t)d * 1024];
    if (tid < 257) soff[tid] = OFF[(size_t)d * 257 + tid];
    const float isq = rsqrtf((float)cN);
    unsigned long long npre = 0, nfull = 0;
    u32 R7 = 0, CLi = 0;
    if (tid < cN) { R7 = cCL[tid] ^ cCR[tid]; CLi = cCL[tid]; }
    u32 MR7 = dM(R7);
    for (int e8 = 0; e8 < 256; e8++) {
        u32 a = 0;
        __syncthreads();
        if (tid < cN) {
            u32 L7 = CLi ^ dN(MR7 ^ (DK7 | ddk7_ext(e8)));
            a = ((u32)(dBY(L7, 0) ^ dBY(L7, 1)) << 8) | (u8)(dBY(L7, 2) ^ dBY(L7, 3));
            if (tid == 0) sdelta = dpart((u8)(a >> 8), (u8)a) ^ u0;
        }
        __syncthreads();
        int delta = sdelta;
        for (int k = 0; k < 256; k++) {
            int lo = soff[k], hi = soff[k + 1];
            for (int c0 = lo; c0 < hi; c0 += SMAX) {
                int s = min(SMAX, hi - c0);
                u32 base = demb_part(k ^ delta);
                u32 r0 = 0, r1 = 0, r2 = 0, r3 = 0;
                if (tid < cN) {
                    const uint4 rw = ((const uint4 *)ROW)[a ^ base];
                    r0 = rw.x ^ cmask; r1 = rw.y ^ cmask; r2 = rw.z ^ cmask; r3 = rw.w ^ cmask;
                }
                r0 = wtrans(r0, lane); r1 = wtrans(r1, lane); r2 = wtrans(r2, lane); r3 = wtrans(r3, lane);
                __syncthreads(); // 前の比較が終わるまで V, Ps を上書きしない
                V[lane][w] = r0; V[32 + lane][w] = r1; V[64 + lane][w] = r2; V[96 + lane][w] = r3;
                for (int x = tid; x < s * NW; x += 1024) {
                    int si = x >> 5, ww = x & 31;
                    u16 e = SORTED[(size_t)d * CAP + c0 + si];
                    Ps[si][ww] = P[(size_t)e * NW + ww];
                    if (ww == 0) Es[si] = e;
                }
                __syncthreads();
                for (int pr = tid; pr < s * 128; pr += 1024) {
                    int si = pr >> 7, r = pr & 127;
                    int dd = 0;
#pragma unroll
                    for (int ww = 0; ww < 8; ww++) dd += __popc(Ps[si][ww] ^ V[r][ww]);
                    int z8 = 256 - 2 * dd;
                    npre++;
                    if (z8 < 32 && z8 > -32) continue;
                    nfull++;
#pragma unroll
                    for (int ww = 8; ww < NW; ww++) dd += __popc(Ps[si][ww] ^ V[r][ww]);
                    float z = (float)(cN - 2 * dd) * isq;
                    if (fabsf(z) >= thrrec) {
                        int p = atomicAdd(CCNT, 1);
                        if (p < CCAP) { Cand cc; cc.t = t; cc.e = Es[si]; cc.kap = (u8)k; cc.e8 = (u8)e8; cc.r = (u8)r; cc.z = z; CB[p] = cc; }
                    }
                }
            }
        }
    }
    atomicAdd(&stats[0], npre); atomicAdd(&stats[1], nfull);
}

// ---------------- ホスト側 ----------------
static int n;
static u32 PLa[1024], PRa[1024], CLa[1024], CRa[1024];
static inline int hprime(int w1, int w2) { return par(G6_Y & feal_N(PK(0, (u8)w1, (u8)w2, 0))) ^ (w1 & 1) ^ (w2 & 1); }
static inline u32 emb_rest(int r) { u32 w1 = (r & 3) | (((r >> 2) & 7) << 5), w2 = (r >> 5) & 3; return (w1 << 8) | w2; }
static inline u32 emb_part(int kp) { return ((u32)(kp & 7) << 10) | ((u32)((kp >> 3) & 31) << 2); }
static inline int partx(u32 x) { return part_idx((u8)(x >> 8), (u8)x); }
static inline u32 dk7_22(u32 t) { return PK((u8)((t >> 16) & 3), (u8)t, (u8)(t >> 8), (u8)((t >> 18) & 15)); }
static inline u32 dk7_ext(int e8) { return ((u32)(e8 & 31) << 26) | ((u32)((e8 >> 5) & 7) << 4); }
static void hwht(float *a) { for (int h = 1; h < 256; h <<= 1) for (int i = 0; i < 256; i += 2 * h) for (int j = i; j < i + h; j++) { float x = a[j], y = a[j + h]; a[j] = x + y; a[j + h] = x - y; } }

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
}

#define CK(x) do { cudaError_t e_ = (x); if (e_ != cudaSuccess) { printf("CUDA %s at %d\n", cudaGetErrorString(e_), __LINE__); exit(1); } } while (0)

int main(int argc, char **argv) {
    int self = !strcmp(argv[1], "self");
    float thr1 = 3.5f, thrrec = self ? 6.0f : 6.5f;
    int ai = self ? 5 : 7;
    if (argc > ai) thr1 = atof(argv[ai]);
    if (argc > ai + 1) thrrec = atof(argv[ai + 1]);
    double T0 = nowsec();
    build_etab();
    float Eh[256], E2[256];
    for (int g = 0; g < 256; g++) { Eh[g] = (float)e_tab[g]; E2[g] = (float)(e_tab[g] * e_tab[g]); }
    hwht(Eh); hwht(E2);
    std::vector<u32> tlist;
    u64 key = 0; u32 t_true = 0; int e_true = -1, k_true = -1, e8_true = -1, r_true = -1;
    FILE *fo = NULL;
    if (self) {
        seed(strtoull(argv[3], 0, 0)); key = rnd();
        load_pairs(argv[2], NULL, 1, key);
        sched_t s; key_schedule(key, &s);
        u32 Kcd = ((u32)s.k[12] << 16) | s.k[13], Kef = ((u32)s.k[14] << 16) | s.k[15];
        u32 K89 = ((u32)s.k[8] << 16) | s.k[9], Kab = ((u32)s.k[10] << 16) | s.k[11];
        u32 EK0 = feal_M(K89 ^ Kab) ^ ((u32)s.k[0] << 8), DK7 = feal_M(Kcd ^ Kef) ^ ((u32)s.k[7] << 8), DK6 = feal_M(Kcd) ^ ((u32)s.k[6] << 8);
        int e1 = BY(EK0, 1), e2 = BY(EK0, 2); if (e2 & 0x80) { e1 ^= 0x80; e2 ^= 0x80; }
        e_true = e1 | (e2 << 8);
        t_true = BY(DK7, 1) | (BY(DK7, 2) << 8) | ((BY(DK7, 0) & 3) << 16) | ((BY(DK7, 3) & 15) << 18);
        e8_true = ((BY(DK7, 0) >> 2) & 31) | (((BY(DK7, 3) >> 4) & 7) << 5);
        u32 DKf = dk7_22(t_true) | dk7_ext(e8_true);
        u32 R7 = CLa[0] ^ CRa[0];
        u32 L7t = CLa[0] ^ feal_N(feal_M(R7) ^ DK7), w = feal_M(L7t) ^ DK6;
        u32 L7e = CLa[0] ^ feal_N(feal_M(R7) ^ DKf), L7s = CLa[0] ^ feal_N(feal_M(R7) ^ dk7_22(t_true));
        u32 ax = ((u32)(BY(L7e, 0) ^ BY(L7e, 1)) << 8) | (BY(L7e, 2) ^ BY(L7e, 3));
        u32 Kx = ((((u32)BY(w, 1) << 8) | BY(w, 2)) ^ ax);
        if (Kx & 0x80) Kx ^= 0x8080;
        int kp = partx(Kx);
        for (int r = 0; r < 128; r++) if ((emb_rest(r) | emb_part(kp)) == Kx) r_true = r;
        int u0 = part_idx(BY(L7s, 0) ^ BY(L7s, 1), BY(L7s, 2) ^ BY(L7s, 3));
        k_true = kp ^ partx(ax) ^ u0;
        printf("key %016llx truth t=%06x e=%04x kappa=%02x e8=%02x r=%02x\n", (unsigned long long)key, t_true, e_true, k_true, e8_true, r_true);
        int nd = atoi(argv[4]);
        tlist.push_back(t_true);
        for (int i = 1; i < nd; i++) tlist.push_back((u32)(rnd() & 0x3FFFFF));
    } else {
        load_pairs(argv[2], argv[3], 0, 0);
        u32 start = strtoul(argv[4], 0, 0), count = strtoul(argv[5], 0, 0);
        for (u32 i = 0; i < count; i++) tlist.push_back(start + i);
        fo = fopen(argv[6], "a");
    }
    // 平文側ビット列（32 語）、ROW 表
    std::vector<u32> Ph((size_t)32768 * NW, 0);
    for (int e = 0; e < 32768; e++) {
        u32 ek = PK(0, (u8)e, (u8)(e >> 8), 0);
        for (int i = 0; i < n; i++) {
            u32 L1 = PLa[i] ^ PRa[i], R1 = PLa[i] ^ feal_N(feal_M(L1) ^ ek);
            if (par(G_L1 & L1) ^ par(G_R1 & R1)) Ph[(size_t)e * NW + (i >> 5)] |= 1u << (i & 31);
        }
    }
    std::vector<u32> ROWh((size_t)65536 * 4);
    static u8 HP[65536];
    for (u32 x = 0; x < 65536; x++) HP[x] = (u8)hprime(x >> 8, x & 255);
    for (u32 x = 0; x < 65536; x++) {
        u32 r4[4] = {0, 0, 0, 0};
        for (int r = 0; r < 128; r++) if (HP[x ^ emb_rest(r)]) r4[r >> 5] |= 1u << (r & 31);
        for (int j = 0; j < 4; j++) ROWh[(size_t)x * 4 + j] = r4[j];
    }
    CK(cudaMemcpyToSymbol(cEhat, Eh, sizeof Eh)); CK(cudaMemcpyToSymbol(cE2hat, E2, sizeof E2));
    CK(cudaMemcpyToSymbol(cCL, CLa, sizeof CLa)); CK(cudaMemcpyToSymbol(cCR, CRa, sizeof CRa)); CK(cudaMemcpyToSymbol(cN, &n, sizeof n));
    u32 *dP, *dROW, *dT, *dC, *dSURV; u8 *dU; float *dTHR, *dZB; int *dSCNT, *dOFF, *dCCNT; u16 *dSORT; Cand *dCB; unsigned long long *dST;
    CK(cudaMalloc(&dP, Ph.size() * 4)); CK(cudaMemcpy(dP, Ph.data(), Ph.size() * 4, cudaMemcpyHostToDevice));
    CK(cudaMalloc(&dROW, ROWh.size() * 4)); CK(cudaMemcpy(dROW, ROWh.data(), ROWh.size() * 4, cudaMemcpyHostToDevice));
    CK(cudaMalloc(&dT, BD * 4)); CK(cudaMalloc(&dU, (size_t)BD * 1024)); CK(cudaMalloc(&dC, (size_t)BD * NW * 4));
    CK(cudaMalloc(&dTHR, (size_t)BD * 256 * 4)); CK(cudaMalloc(&dSURV, (size_t)BD * CAP * 4)); CK(cudaMalloc(&dSCNT, BD * 4));
    CK(cudaMalloc(&dSORT, (size_t)BD * CAP * 2)); CK(cudaMalloc(&dOFF, (size_t)BD * 257 * 4)); CK(cudaMalloc(&dCB, (size_t)CCAP * sizeof(Cand)));
    CK(cudaMalloc(&dCCNT, 4)); CK(cudaMalloc(&dST, 16)); CK(cudaMalloc(&dZB, 256 * 4));
    printf("N=%d thr1=%.2f thrrec=%.2f setup %.1fs, %zu DK7 values\n", n, thr1, thrrec, nowsec() - T0, tlist.size());
    fflush(stdout);
    std::vector<Cand> all;
    std::vector<int> scnt(BD);
    long long tot_surv = 0, overflow = 0; unsigned long long st[2] = {0, 0};
    float zb_true[256] = {0}; int true_surv = 0;
    double t1s = 0, t2s = 0, W0 = nowsec();
    for (size_t b0 = 0; b0 < tlist.size(); b0 += BD) {
        int nd = (int)std::min((size_t)BD, tlist.size() - b0);
        CK(cudaMemcpy(dT, tlist.data() + b0, nd * 4, cudaMemcpyHostToDevice));
        CK(cudaMemset(dSCNT, 0, BD * 4)); CK(cudaMemset(dCCNT, 0, 4)); CK(cudaMemset(dST, 0, 16));
        double a = nowsec();
        k_prep<<<nd, 1024>>>(dT, nd, thr1, dU, dC, dTHR);
        dim3 g1(32768 / 64, nd);
        k_stage1<<<g1, 256>>>(dP, dU, dC, dTHR, nd, dSURV, dSCNT, (self && b0 == 0) ? 0 : -1, e_true, dZB);
        k_sort<<<nd, 1024>>>(dSURV, dSCNT, nd, dSORT, dOFF);
        CK(cudaDeviceSynchronize());
        double bq = nowsec(); t1s += bq - a;
        k_stage2<<<nd, 1024>>>(dT, dP, dROW, dU, dC, dSORT, dOFF, nd, thrrec, dCB, dCCNT, dST);
        CK(cudaDeviceSynchronize());
        t2s += nowsec() - bq;
        int ncand; CK(cudaMemcpy(&ncand, dCCNT, 4, cudaMemcpyDeviceToHost));
        if (ncand > CCAP) { printf("candidate buffer overflow %d\n", ncand); ncand = CCAP; }
        std::vector<Cand> cb(ncand);
        CK(cudaMemcpy(cb.data(), dCB, (size_t)ncand * sizeof(Cand), cudaMemcpyDeviceToHost));
        CK(cudaMemcpy(scnt.data(), dSCNT, nd * 4, cudaMemcpyDeviceToHost));
        unsigned long long s2[2]; CK(cudaMemcpy(s2, dST, 16, cudaMemcpyDeviceToHost)); st[0] += s2[0]; st[1] += s2[1];
        for (int i = 0; i < nd; i++) { tot_surv += std::min(scnt[i], CAP); if (scnt[i] > CAP) overflow++; }
        if (self && b0 == 0) {
            CK(cudaMemcpy(zb_true, dZB, sizeof zb_true, cudaMemcpyDeviceToHost));
            std::vector<u32> sv(std::min(scnt[0], CAP));
            CK(cudaMemcpy(sv.data(), dSURV, sv.size() * 4, cudaMemcpyDeviceToHost));
            for (u32 v : sv) if (v == ((u32)e_true | ((u32)k_true << 15))) true_surv = 1;
        }
        if (fo) { for (auto &c : cb) fprintf(fo, "%06x %04x %02x %02x %02x %.4f\n", c.t, c.e, c.kap, c.e8, c.r, c.z); fflush(fo); }
        else all.insert(all.end(), cb.begin(), cb.end());
        if (!self) { printf("batch %zu-%zu: %d cands, %.1fs elapsed\n", (size_t)tlist[b0], (size_t)tlist[b0 + nd - 1], ncand, nowsec() - W0); fflush(stdout); }
    }
    double W1 = nowsec(), nd = (double)tlist.size();
    printf("processed %.0f DK7 in %.1fs (stage1 %.1fs, stage2 %.1fs) -> %.4f s/DK7; full 2^22: %.2f h\n", nd, W1 - W0, t1s, t2s,
           (W1 - W0) / nd, (W1 - W0) / nd * 4194304.0 / 3600);
    printf("survivors/DK7 %.0f (overflow DK7s %lld); pre-compares/DK7 %.3g, full %.3g (%.2f%%); candidates/DK7 %.2f\n",
           tot_surv / nd, overflow, st[0] / nd, st[1] / nd, 100.0 * st[1] / std::max(1ULL, st[0]), all.size() / nd);
    if (self) {
        printf("true hypothesis: stage1 zB=%.3f survived=%s\n", zb_true[k_true] * thr1, true_surv ? "yes" : "no");
        std::sort(all.begin(), all.end(), [](const Cand &x, const Cand &y) { return fabsf(x.z) > fabsf(y.z); });
        int rank = -1;
        for (size_t c = 0; c < all.size(); c++)
            if (all[c].t == t_true && all[c].e == e_true && all[c].kap == k_true && all[c].e8 == e8_true && all[c].r == r_true) { rank = (int)c; break; }
        long long c65 = 0, c70 = 0;
        for (auto &c : all) if (c.t != t_true) { c65 += fabsf(c.z) >= 6.5f; c70 += fabsf(c.z) >= 7.0f; }
        printf("true key in candidates: %s (rank %d, z=%.3f); wrong >=6.5: %lld, >=7.0: %lld (full-space x%.0f)\n",
               rank >= 0 ? "yes" : "no", rank, rank >= 0 ? all[rank].z : 0.f, c65, c70, 4194304.0 / (nd - 1));
        if (getenv("CANDOUT")) {
            FILE *f = fopen(getenv("CANDOUT"), "w");
            for (auto &c : all) fprintf(f, "%06x %04x %02x %02x %02x %.4f %d\n", c.t, c.e, c.kap, c.e8, c.r, c.z, c.t == t_true ? 1 : 0);
            fclose(f);
        }
    } else fclose(fo);
    return 0;
}
