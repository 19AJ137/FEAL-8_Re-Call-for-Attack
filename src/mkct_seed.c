/* attackC2/stage3 の self と同じ鍵（seed → rnd() の最初の値）で plaintext.bin を ECB 暗号化する
 *   ./mkct_seed plaintext.bin seed out.bin */
#include "feal8.h"
int main(int argc, char **argv) {
    seed(strtoull(argv[2], 0, 0));
    u64 key = rnd();
    sched_t s; key_schedule(key, &s);
    FILE *fp = fopen(argv[1], "rb"); u8 b[8192];
    if (fread(b, 1, 8192, fp) != 8192) return 1;
    fclose(fp);
    for (int i = 0; i < 1024; i++) {
        u64 p = 0;
        for (int j = 0; j < 8; j++) p = (p << 8) | b[8 * i + j];
        u64 c = feal8_enc(&s, p, NULL);
        for (int j = 0; j < 8; j++) b[8 * i + j] = (u8)(c >> (56 - 8 * j));
    }
    fp = fopen(argv[3], "wb"); fwrite(b, 1, 8192, fp); fclose(fp);
    printf("key %016llx\n", (unsigned long long)key);
    return 0;
}
