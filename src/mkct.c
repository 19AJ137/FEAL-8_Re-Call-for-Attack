/* 既知鍵で plaintext.bin を ECB 暗号化して書き出す（実データ経路の試験用） */
#include "feal8.h"
int main(int argc, char **argv) {
    u64 key = strtoull(argv[1], 0, 16);
    FILE *fp = fopen(argv[2], "rb"); u8 b[8192]; if (fread(b, 1, 8192, fp) != 8192) return 1; fclose(fp);
    sched_t s; key_schedule(key, &s);
    for (int i = 0; i < 1024; i++) {
        u64 p = 0; for (int j = 0; j < 8; j++) p = (p << 8) | b[8 * i + j];
        u64 c = feal8_enc(&s, p, NULL);
        for (int j = 0; j < 8; j++) b[8 * i + j] = (u8)(c >> (56 - 8 * j));
    }
    fp = fopen(argv[3], "wb"); fwrite(b, 1, 8192, fp); fclose(fp);
    return 0;
}
