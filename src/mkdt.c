/* 既知鍵で plaintext.bin を FEAL-8 で「復号」して書き出す（C = D_K(P) の試験データ）。D(E(x)) = x も確認する
 *   ./mkdt plaintext.bin seed out.bin   （鍵は self と同じ seed → rnd() の最初の値） */
#include "feal8.h"
static u64 feal8_dec(const sched_t *s, u64 c) {
    u32 CL = (u32)(c >> 32), CR = (u32)c;
    u32 R = CL ^ (((u32)s->k[12] << 16) | s->k[13]);
    u32 L = CR ^ (((u32)s->k[14] << 16) | s->k[15]) ^ R;
    for (int r = 7; r >= 0; r--) { u32 Lp = R ^ feal_f(L, s->k[r]); R = L; L = Lp; }
    R ^= L;
    R ^= ((u32)s->k[10] << 16) | s->k[11];
    L ^= ((u32)s->k[8] << 16) | s->k[9];
    return ((u64)L << 32) | R;
}
int main(int argc, char **argv) {
    seed(strtoull(argv[2], 0, 0)); u64 key = rnd();
    sched_t s; key_schedule(key, &s);
    for (int i = 0; i < 100000; i++) { u64 x = rnd(); if (feal8_dec(&s, feal8_enc(&s, x, NULL)) != x) { printf("D(E(x)) != x\n"); return 1; } }
    FILE *fp = fopen(argv[1], "rb"); u8 b[8192]; if (fread(b, 1, 8192, fp) != 8192) return 1; fclose(fp);
    for (int i = 0; i < 1024; i++) {
        u64 p = 0; for (int j = 0; j < 8; j++) p = (p << 8) | b[8 * i + j];
        u64 c = feal8_dec(&s, p);
        for (int j = 0; j < 8; j++) b[8 * i + j] = (u8)(c >> (56 - 8 * j));
    }
    fp = fopen(argv[3], "wb"); fwrite(b, 1, 8192, fp); fclose(fp);
    printf("key %016llx, D(E(x))=x OK (1e5)\n", (unsigned long long)key);
    return 0;
}
