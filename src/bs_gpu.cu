// FEAL-8 ビットスライス総当たり（CUDA）。32 鍵を u32 の各ビットに載せる。
//   ./bs_gpu test   : 既知鍵を範囲に埋めて探索し、見つかるかと CL 一致数を確認
//   ./bs_gpu bench  [launches] : 速度測定
// 鍵 64 ビット = A0(上位32) | B0(下位32)。鍵ビット 0..4 = レーン、5..24 = スレッド番号、
// 25..35 = スレッド内ループ、36..63 = 起動ごとの基準値。
// 鍵スケジュールは 2 回に分けて計算する（レジスタ節約）：
//   1 回目：B1..B7 を計算し、B5,B6 で前処理、B7 を照合値に畳み込む
//   2 回目：B1..B4 を再計算しながら 8 段を回す。CL 一致時だけ B5..B8 を再計算して CR も照合
#include "feal8.h"
#include <cuda_runtime.h>
#include <chrono>

#define ONES 0xFFFFFFFFu
#define TPB 256
#define NBLK 4096
#define ITER_BITS 11

__constant__ u32 cPL[32], cPR[32], cCL[32], cCR[32]; // 平文・暗号文の各ビット（0 か ~0）
__device__ unsigned long long d_found[64];
__device__ unsigned int d_nfound, d_nclmatch;

// 1 バイトのビット列：b[0] が最下位
__device__ __forceinline__ void add8(const u32 *x, const u32 *y, u32 c, u32 *s) {
#pragma unroll
    for (int i = 0; i < 8; i++) {
        u32 xi = x[i], yi = y[i];
        s[i] = xi ^ yi ^ c;
        if (i < 7) c = (xi & yi) | (c & (xi ^ yi));
    }
}
// S_d(x,y) = ROL2(x+y+d)：出力ビット k = 和のビット (k-2) mod 8
__device__ __forceinline__ void sbox(const u32 *x, const u32 *y, u32 d, u32 *o) {
    u32 s[8];
    add8(x, y, d, s);
#pragma unroll
    for (int k = 0; k < 8; k++) o[k] = s[(k + 6) & 7];
}
// 32 ビット語：w[j*8+b] = バイト j（0 が最上位）のビット b
__device__ __forceinline__ void fK(const u32 *a, const u32 *b, u32 *o) {
    u32 k1[8], k2[8], t[8];
#pragma unroll
    for (int i = 0; i < 8; i++) { k1[i] = a[8 + i] ^ a[i]; k2[i] = a[16 + i] ^ a[24 + i]; t[i] = k2[i] ^ b[i]; }
    sbox(k1, t, ONES, o + 8);
#pragma unroll
    for (int i = 0; i < 8; i++) t[i] = o[8 + i] ^ b[8 + i];
    sbox(k2, t, 0, o + 16);
#pragma unroll
    for (int i = 0; i < 8; i++) t[i] = o[8 + i] ^ b[16 + i];
    sbox(a, t, 0, o);
#pragma unroll
    for (int i = 0; i < 8; i++) t[i] = o[16 + i] ^ b[24 + i];
    sbox(a + 24, t, ONES, o + 24);
}
// L ^= f(R, (k0,k1))、k0,k1 は 8 ビットずつ
__device__ __forceinline__ void round_xor(u32 *L, const u32 *a, const u32 *k0, const u32 *k1) {
    u32 f1[8], f2[8], o1[8], o2[8], o[8];
#pragma unroll
    for (int i = 0; i < 8; i++) { f1[i] = a[8 + i] ^ k0[i] ^ a[i]; f2[i] = a[16 + i] ^ k1[i] ^ a[24 + i]; }
    sbox(f1, f2, ONES, o1);
    sbox(f2, o1, 0, o2);
#pragma unroll
    for (int i = 0; i < 8; i++) { L[8 + i] ^= o1[i]; L[16 + i] ^= o2[i]; }
    sbox(a, o1, 0, o);
#pragma unroll
    for (int i = 0; i < 8; i++) L[i] ^= o[i];
    sbox(a + 24, o2, ONES, o);
#pragma unroll
    for (int i = 0; i < 8; i++) L[24 + i] ^= o[i];
}
// 鍵ビット（語内ビット位置 p、最下位=0）を w[] の添字へ
__device__ __forceinline__ int idx(int p) { return (3 - (p >> 3)) * 8 + (p & 7); }

__device__ __forceinline__ void key_slices(u64 kbase, u32 *A0, u32 *B0) {
    const u32 lane[5] = {0xAAAAAAAAu, 0xCCCCCCCCu, 0xF0F0F0F0u, 0xFF00FF00u, 0xFFFF0000u};
#pragma unroll
    for (int p = 0; p < 32; p++) {
        A0[idx(p)] = (u32)0 - (u32)((kbase >> (32 + p)) & 1);
        B0[idx(p)] = p < 5 ? lane[p] : (u32)0 - (u32)((kbase >> p) & 1);
    }
}

__global__ void __launch_bounds__(TPB) search(u64 base) {
    u64 tid = (u64)blockIdx.x * TPB + threadIdx.x;
    for (u32 it = 0; it < (1u << ITER_BITS); it++) {
        u64 kb = base | ((u64)it << 25) | (tid << 5);
        u32 A0[32], B0[32];
        key_slices(kb, A0, B0);
        u32 x[32], y[32], z[32];       // 窓（鍵スケジュールの直近 3 語）
        u32 L[32], R[32], tg[32];
        // 1 回目：B1..B7
        fK(A0, B0, x);                                  // B1
#pragma unroll
        for (int i = 0; i < 32; i++) z[i] = x[i] ^ A0[i];
        fK(B0, z, y);                                   // B2
#pragma unroll
        for (int i = 0; i < 32; i++) z[i] = y[i] ^ B0[i];
        fK(x, z, R);                                    // B3（R を一時置き場に）
#pragma unroll
        for (int i = 0; i < 32; i++) z[i] = R[i] ^ x[i];
        fK(y, z, L);                                    // B4
#pragma unroll
        for (int i = 0; i < 32; i++) z[i] = L[i] ^ y[i];
        fK(R, z, x);                                    // B5
#pragma unroll
        for (int i = 0; i < 32; i++) z[i] = x[i] ^ R[i];
        fK(L, z, y);                                    // B6
#pragma unroll
        for (int i = 0; i < 32; i++) z[i] = y[i] ^ L[i];
        fK(x, z, tg);                                   // B7
#pragma unroll
        for (int i = 0; i < 32; i++) {
            tg[i] ^= cCL[i];                            // CL ^ B7 = 目標の R8
            L[i] = cPL[i] ^ x[i];                       // L0 = PL ^ (K8,K9)
        }
#pragma unroll
        for (int i = 0; i < 32; i++) R[i] = cPR[i] ^ y[i] ^ L[i]; // R0 = PR ^ (K10,K11) ^ L0
        // 2 回目：B1..B4 を作りながら 8 段
        fK(A0, B0, x);                                  // B1 -> 段 0,1
        round_xor(L, R, x, x + 8);                      // 段0：L1=R0, R1=L0^f → (R, L) の役割交換
        round_xor(R, L, x + 16, x + 24);                // 段1
#pragma unroll
        for (int i = 0; i < 32; i++) z[i] = x[i] ^ A0[i];
        fK(B0, z, y);                                   // B2 -> 段 2,3
        round_xor(L, R, y, y + 8);
        round_xor(R, L, y + 16, y + 24);
#pragma unroll
        for (int i = 0; i < 32; i++) z[i] = y[i] ^ B0[i];
        u32 w3[32];
        fK(x, z, w3);                                   // B3 -> 段 4,5
        round_xor(L, R, w3, w3 + 8);
        round_xor(R, L, w3 + 16, w3 + 24);
#pragma unroll
        for (int i = 0; i < 32; i++) z[i] = w3[i] ^ x[i];
        fK(y, z, x);                                    // B4 -> 段 6,7
        round_xor(L, R, x, x + 8);
        round_xor(R, L, x + 16, x + 24);
        // 8 段後：状態 (L8,R8) = (L, R)。CL = R8 ^ B7 を照合
        u32 diff = 0;
#pragma unroll
        for (int i = 0; i < 32; i++) diff |= R[i] ^ tg[i];
        u32 hit = ~diff;
        if (hit) {
            atomicAdd(&d_nclmatch, __popc(hit));
            // CR = L8 ^ R8 ^ B8 を確認するため B5..B8 をもう一度
            u32 b1[32], b2[32], b3[32], b4[32], b5[32], b6[32], b7[32], b8[32];
            fK(A0, B0, b1);
            for (int i = 0; i < 32; i++) z[i] = b1[i] ^ A0[i];
            fK(B0, z, b2);
            for (int i = 0; i < 32; i++) z[i] = b2[i] ^ B0[i];
            fK(b1, z, b3);
            for (int i = 0; i < 32; i++) z[i] = b3[i] ^ b1[i];
            fK(b2, z, b4);
            for (int i = 0; i < 32; i++) z[i] = b4[i] ^ b2[i];
            fK(b3, z, b5);
            for (int i = 0; i < 32; i++) z[i] = b5[i] ^ b3[i];
            fK(b4, z, b6);
            for (int i = 0; i < 32; i++) z[i] = b6[i] ^ b4[i];
            fK(b5, z, b7);
            for (int i = 0; i < 32; i++) z[i] = b7[i] ^ b5[i];
            fK(b6, z, b8);
            u32 d2 = 0;
            for (int i = 0; i < 32; i++) d2 |= L[i] ^ R[i] ^ b8[i] ^ cCR[i];
            u32 full = hit & ~d2;
            while (full) {
                int l = __ffs(full) - 1; full &= full - 1;
                unsigned int k = atomicAdd(&d_nfound, 1);
                if (k < 64) d_found[k] = kb | (u64)l;
            }
        }
    }
}

static void set_consts(u64 p, u64 c) {
    u32 PLs[32], PRs[32], CLs[32], CRs[32];
    for (int pidx = 0; pidx < 32; pidx++) {
        int j = (3 - (pidx >> 3)) * 8 + (pidx & 7);
        PLs[j] = ((p >> (32 + pidx)) & 1) ? ONES : 0; PRs[j] = ((p >> pidx) & 1) ? ONES : 0;
        CLs[j] = ((c >> (32 + pidx)) & 1) ? ONES : 0; CRs[j] = ((c >> pidx) & 1) ? ONES : 0;
    }
    cudaMemcpyToSymbol(cPL, PLs, sizeof PLs); cudaMemcpyToSymbol(cPR, PRs, sizeof PRs);
    cudaMemcpyToSymbol(cCL, CLs, sizeof CLs); cudaMemcpyToSymbol(cCR, CRs, sizeof CRs);
}

int main(int argc, char **argv) {
    const u64 per_launch = (u64)NBLK * TPB * 32 << ITER_BITS; // 2^36
    cudaDeviceProp prop; cudaGetDeviceProperties(&prop, 0);
    int clk = 0; cudaDeviceGetAttribute(&clk, cudaDevAttrClockRate, 0);
    printf("GPU %s, SMs %d, clock %d MHz; keys per launch 2^%d\n", prop.name, prop.multiProcessorCount,
           clk / 1000, 63 - __builtin_clzll(per_launch));
    seed(12345);
    if (argc > 1 && !strcmp(argv[1], "test")) {
        for (int t = 0; t < 3; t++) {
            u64 key = rnd(), p = rnd();
            sched_t s; key_schedule(key, &s);
            u64 c = feal8_enc(&s, p, NULL);
            set_consts(p, c);
            u64 base = key & ~(per_launch - 1);
            unsigned int z = 0;
            cudaMemcpyToSymbol(d_nfound, &z, 4); cudaMemcpyToSymbol(d_nclmatch, &z, 4);
            search<<<NBLK, TPB>>>(base);
            cudaError_t e = cudaDeviceSynchronize();
            unsigned int nf, ncl; u64 f[64];
            cudaMemcpyFromSymbol(&nf, d_nfound, 4); cudaMemcpyFromSymbol(&ncl, d_nclmatch, 4);
            cudaMemcpyFromSymbol(f, d_found, sizeof f);
            int ok = 0;
            for (unsigned i = 0; i < nf && i < 64; i++) if (f[i] == key) ok = 1;
            printf("test %d: key %016llx  %s  found %u (CL matches %u, expect ~16)  -> %s\n", t,
                   (unsigned long long)key, cudaGetErrorString(e), nf, ncl, ok ? "FOUND" : "MISSING");
        }
        return 0;
    }
    int launches = argc > 2 ? atoi(argv[2]) : 8;
    u64 key = rnd(), p = rnd();
    sched_t s; key_schedule(key, &s);
    set_consts(p, feal8_enc(&s, p, NULL));
    search<<<NBLK, TPB>>>(0); cudaDeviceSynchronize(); // ウォームアップ
    auto t0 = std::chrono::steady_clock::now();
    for (int l = 0; l < launches; l++) search<<<NBLK, TPB>>>((u64)(l + 1) * per_launch);
    cudaDeviceSynchronize();
    double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
    double rate = (double)per_launch * launches / sec;
    printf("bench: %d launches, %.3f s, %.3e keys/s; 2^63 keys -> %.2f years\n", launches, sec, rate,
           9.223372036854776e18 / rate / 3.156e7);
    return 0;
}
