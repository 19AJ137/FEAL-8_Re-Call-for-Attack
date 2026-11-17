# References / 参考文献

## The challenge and the data / 懸賞問題とデータ

1. S. Miyaguchi, "The FEAL-8 Cryptosystem and a Call for Attack", Advances in Cryptology – CRYPTO '89,
   LNCS 435, pp. 624–627, Springer, 1990. https://doi.org/10.1007/0-387-34805-0_59
   — Rules of the call: preliminary problem (8192 bytes of plaintext and ciphertext, ECB, one fixed key),
   final problems, prize, application period (1989-08-20 to 1991-08-31).
2. 松井甲子雄『コンピュータによる暗号解読法入門』森北出版, 1990（第1版第1刷 1990-02-24）, ISBN 4-627-83350-4.
   — Source of the plaintext, ciphertext and CBC-MAC values in [`../data`](../data) (appendix figures 6–9,
   around pp. 197–206; transcribed from the 2nd printing of 1997).
   本リポジトリの平文・暗号文・MAC 値の出典（付図6〜9、pp. 197–206 付近。1997 年の第2刷から書き起こし）。
3. S. Miyaguchi, A. Shiraishi, A. Shimizu, "Fast Data Encipherment Algorithm FEAL-8",
   Review of the Electrical Communication Laboratories, Vol. 36, No. 4, 1988.
   — The text in the plaintext is the beginning of this paper (Japanese and English).
   平文の文章部分はこの論文の冒頭（日本語と英語）。
4. A. Shimizu, S. Miyaguchi, "Fast Data Encipherment Algorithm FEAL", Advances in Cryptology –
   EUROCRYPT '87, LNCS 304, Springer, 1988.
5. M. Matsui, "Celebrating the 25th year of FEAL – A New Prize Problem", CRYPTO 2012 rump session.
   https://crypto.2012.rump.cr.yp.to/19997d5a295baee62c05ba73534745ef.pdf
   — Notes that the 1989 problem remained unsolved.

## Attacks used / 使った攻撃法

6. E. Biham, Y. Carmeli, "An Improvement of Linear Cryptanalysis with Addition Operations with
   Applications to FEAL-8X", Selected Areas in Cryptography – SAC 2014, LNCS 8781, Springer, 2014.
   https://doi.org/10.1007/978-3-319-13051-4_4
   — Equivalent representation, the 5- and 6-round approximations, partitioning, and the 2^10-pair
   attack of Section 5.
7. Y. Sakikoyama, Y. Todo, K. Aoki, M. Morii, "Efficient Implementations for Practical Linear
   Cryptanalysis and Its Application to FEAL-8X", IEICE Trans. Fundamentals, Vol. E99-A, No. 1,
   pp. 31–38, 2016. https://doi.org/10.1587/transfun.E99.A.31
   — Multiple linear cryptanalysis with 2^10 pairs; the second family of approximations
   (ΓY = 04011904) used here comes from its appendix.
8. G. Leurent, "Improved Differential-Linear Cryptanalysis of 7-round Chaskey with Partitioning",
   Advances in Cryptology – EUROCRYPT 2016 (Appendix A: generalised partitioning of FEAL-8X).
   https://who.rocq.inria.fr/Gaetan.Leurent/files/Chaskey_EC16.pdf

## Background / その他

9. M. Matsui, A. Yamagishi, "A New Method for Known Plaintext Attack of FEAL Cipher",
   Advances in Cryptology – EUROCRYPT '92, LNCS 658, Springer, 1993.
   https://doi.org/10.1007/3-540-47555-9_7
10. K. Ohta, K. Aoki, "Linear Cryptanalysis of the Fast Data Encipherment Algorithm",
    Advances in Cryptology – CRYPTO '94, LNCS 839, Springer, 1994. https://doi.org/10.1007/3-540-48658-5_2
11. E. Biham, A. Shamir, "Differential Cryptanalysis of Feal and N-Hash", Advances in Cryptology –
    EUROCRYPT '91, LNCS 547, Springer, 1991. https://doi.org/10.1007/3-540-46416-6_1
