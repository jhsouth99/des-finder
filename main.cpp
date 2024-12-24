#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>

#define WIN32
#if defined(WIN32)
#include "os_windows.h"
#elif defined(__linux__)
#include "os_linux.h"
#elif defined(__APPLE__)
#include "os_osx.h"
#else
#error Unknown OS!
#endif

#define USE_BIT_PERMUTATION

#ifdef USE_BIT_PERMUTATION
#include <immintrin.h>

#define pext(source, mask) _pext_u32(source, mask)

int bitCount_u32_PeterWegner(uint32_t n) {
	uint32_t c;
	for (c = 0; n; c++)
		n &= (n - 1);
	return c;
}

int bitCount_u32_HammingWeight(uint32_t n) {
	const uint32_t m1 = 0x55555555; //01 01 01 01...
	const uint32_t m2 = 0x33333333; //00 11 00 11 ...
	const uint32_t m3 = 0x0f0f0f0f; //00 00 11 11
	const uint32_t h1 = 0x01010101; //00 00 00 01
	n = n - ((n >> 1) & m1);
	n = (n & m2) + ((n >> 2) & m2);
	return (((n + (n >> 4)) & m3) * h1) >> 24;
}

#define bitCount_u32(n) bitCount_u32_HammingWeight(n)

#define USE_DELTA_SWAP

#ifdef USE_DELTA_SWAP
#define DELTA_SWAP(a, b, delta, mask)   \
    b = (a ^ (a >> delta)) & mask;      \
    a = a ^ b ^ (b << delta);
#endif

#endif

uint64_t ArrToUll(const uint8_t _In_(&_Arr)[8]) {
	uint64_t _Output = 0ULL;
	for (size_t i = 0; i < 8; ++i) {
		_Output |= (uint64_t)_Arr[i] << (8 * (7 - i));
	}
	return _Output;
}

void UllToArr(
	uint64_t _In_ _Input,
	uint8_t _Out_(&_Output)[8]
) {
	for (size_t i = 0; i < 8; ++i) {
		_Output[7 - i] = (_Input >> (8 * i)) & 0xFF;
	}
}

static const uint32_t des_skb[8][64] = {
	{
		/* for C bits (numbered as per FIPS 46) 1 2 3 4 5 6 */
		// 0010 0000 0000 1001 0000 1000 0011 0000
		0x00000000L, 0x00000010L, 0x20000000L, 0x20000010L,
		0x00010000L, 0x00010010L, 0x20010000L, 0x20010010L,
		0x00000800L, 0x00000810L, 0x20000800L, 0x20000810L,
		0x00010800L, 0x00010810L, 0x20010800L, 0x20010810L,
		0x00000020L, 0x00000030L, 0x20000020L, 0x20000030L,
		0x00010020L, 0x00010030L, 0x20010020L, 0x20010030L,
		0x00000820L, 0x00000830L, 0x20000820L, 0x20000830L,
		0x00010820L, 0x00010830L, 0x20010820L, 0x20010830L,
		0x00080000L, 0x00080010L, 0x20080000L, 0x20080010L,
		0x00090000L, 0x00090010L, 0x20090000L, 0x20090010L,
		0x00080800L, 0x00080810L, 0x20080800L, 0x20080810L,
		0x00090800L, 0x00090810L, 0x20090800L, 0x20090810L,
		0x00080020L, 0x00080030L, 0x20080020L, 0x20080030L,
		0x00090020L, 0x00090030L, 0x20090020L, 0x20090030L,
		0x00080820L, 0x00080830L, 0x20080820L, 0x20080830L,
		0x00090820L, 0x00090830L, 0x20090820L, 0x20090830L,
	},
	{
		/* for C bits (numbered as per FIPS 46) 7 8 10 11 12 13 */
		// 0001 0010 0010 0000 0010 0100 0000 0100
		0x00000000L, 0x02000000L, 0x00002000L, 0x02002000L,
		0x00200000L, 0x02200000L, 0x00202000L, 0x02202000L,
		0x00000004L, 0x02000004L, 0x00002004L, 0x02002004L,
		0x00200004L, 0x02200004L, 0x00202004L, 0x02202004L,
		0x00000400L, 0x02000400L, 0x00002400L, 0x02002400L,
		0x00200400L, 0x02200400L, 0x00202400L, 0x02202400L,
		0x00000404L, 0x02000404L, 0x00002404L, 0x02002404L,
		0x00200404L, 0x02200404L, 0x00202404L, 0x02202404L,
		0x10000000L, 0x12000000L, 0x10002000L, 0x12002000L,
		0x10200000L, 0x12200000L, 0x10202000L, 0x12202000L,
		0x10000004L, 0x12000004L, 0x10002004L, 0x12002004L,
		0x10200004L, 0x12200004L, 0x10202004L, 0x12202004L,
		0x10000400L, 0x12000400L, 0x10002400L, 0x12002400L,
		0x10200400L, 0x12200400L, 0x10202400L, 0x12202400L,
		0x10000404L, 0x12000404L, 0x10002404L, 0x12002404L,
		0x10200404L, 0x12200404L, 0x10202404L, 0x12202404L,
	},
	{
		/* for C bits (numbered as per FIPS 46) 14 15 16 17 19 20 */
		// 0000 1001 0000 0100 0000 0010 0000 0011
		0x00000000L, 0x00000001L, 0x00040000L, 0x00040001L,
		0x01000000L, 0x01000001L, 0x01040000L, 0x01040001L,
		0x00000002L, 0x00000003L, 0x00040002L, 0x00040003L,
		0x01000002L, 0x01000003L, 0x01040002L, 0x01040003L,
		0x00000200L, 0x00000201L, 0x00040200L, 0x00040201L,
		0x01000200L, 0x01000201L, 0x01040200L, 0x01040201L,
		0x00000202L, 0x00000203L, 0x00040202L, 0x00040203L,
		0x01000202L, 0x01000203L, 0x01040202L, 0x01040203L,
		0x08000000L, 0x08000001L, 0x08040000L, 0x08040001L,
		0x09000000L, 0x09000001L, 0x09040000L, 0x09040001L,
		0x08000002L, 0x08000003L, 0x08040002L, 0x08040003L,
		0x09000002L, 0x09000003L, 0x09040002L, 0x09040003L,
		0x08000200L, 0x08000201L, 0x08040200L, 0x08040201L,
		0x09000200L, 0x09000201L, 0x09040200L, 0x09040201L,
		0x08000202L, 0x08000203L, 0x08040202L, 0x08040203L,
		0x09000202L, 0x09000203L, 0x09040202L, 0x09040203L,
	},
	{
		/* for C bits (numbered as per FIPS 46) 21 23 24 26 27 28 */
		// 0000 0100 0001 0010 0001 0001 0000 1000
		0x00000000L, 0x00100000L, 0x00000100L, 0x00100100L,
		0x00000008L, 0x00100008L, 0x00000108L, 0x00100108L,
		0x00001000L, 0x00101000L, 0x00001100L, 0x00101100L,
		0x00001008L, 0x00101008L, 0x00001108L, 0x00101108L,
		0x04000000L, 0x04100000L, 0x04000100L, 0x04100100L,
		0x04000008L, 0x04100008L, 0x04000108L, 0x04100108L,
		0x04001000L, 0x04101000L, 0x04001100L, 0x04101100L,
		0x04001008L, 0x04101008L, 0x04001108L, 0x04101108L,
		0x00020000L, 0x00120000L, 0x00020100L, 0x00120100L,
		0x00020008L, 0x00120008L, 0x00020108L, 0x00120108L,
		0x00021000L, 0x00121000L, 0x00021100L, 0x00121100L,
		0x00021008L, 0x00121008L, 0x00021108L, 0x00121108L,
		0x04020000L, 0x04120000L, 0x04020100L, 0x04120100L,
		0x04020008L, 0x04120008L, 0x04020108L, 0x04120108L,
		0x04021000L, 0x04121000L, 0x04021100L, 0x04121100L,
		0x04021008L, 0x04121008L, 0x04021108L, 0x04121108L,
	},
	{
		/* for D bits (numbered as per FIPS 46) 1 2 3 4 5 6 */
		// 0011 0000 0001 0001 0001 0000 0000 0100
		0x00000000L, 0x10000000L, 0x00010000L, 0x10010000L,
		0x00000004L, 0x10000004L, 0x00010004L, 0x10010004L,
		0x20000000L, 0x30000000L, 0x20010000L, 0x30010000L,
		0x20000004L, 0x30000004L, 0x20010004L, 0x30010004L,
		0x00100000L, 0x10100000L, 0x00110000L, 0x10110000L,
		0x00100004L, 0x10100004L, 0x00110004L, 0x10110004L,
		0x20100000L, 0x30100000L, 0x20110000L, 0x30110000L,
		0x20100004L, 0x30100004L, 0x20110004L, 0x30110004L,
		0x00001000L, 0x10001000L, 0x00011000L, 0x10011000L,
		0x00001004L, 0x10001004L, 0x00011004L, 0x10011004L,
		0x20001000L, 0x30001000L, 0x20011000L, 0x30011000L,
		0x20001004L, 0x30001004L, 0x20011004L, 0x30011004L,
		0x00101000L, 0x10101000L, 0x00111000L, 0x10111000L,
		0x00101004L, 0x10101004L, 0x00111004L, 0x10111004L,
		0x20101000L, 0x30101000L, 0x20111000L, 0x30111000L,
		0x20101004L, 0x30101004L, 0x20111004L, 0x30111004L,
	},
	{
		/* for D bits (numbered as per FIPS 46) 8 9 11 12 13 14 */
		// 0000 1010 0000 0010 0000 0100 0000 1001
		0x00000000L, 0x08000000L, 0x00000008L, 0x08000008L,
		0x00000400L, 0x08000400L, 0x00000408L, 0x08000408L,
		0x00020000L, 0x08020000L, 0x00020008L, 0x08020008L,
		0x00020400L, 0x08020400L, 0x00020408L, 0x08020408L,
		0x00000001L, 0x08000001L, 0x00000009L, 0x08000009L,
		0x00000401L, 0x08000401L, 0x00000409L, 0x08000409L,
		0x00020001L, 0x08020001L, 0x00020009L, 0x08020009L,
		0x00020401L, 0x08020401L, 0x00020409L, 0x08020409L,
		0x02000000L, 0x0A000000L, 0x02000008L, 0x0A000008L,
		0x02000400L, 0x0A000400L, 0x02000408L, 0x0A000408L,
		0x02020000L, 0x0A020000L, 0x02020008L, 0x0A020008L,
		0x02020400L, 0x0A020400L, 0x02020408L, 0x0A020408L,
		0x02000001L, 0x0A000001L, 0x02000009L, 0x0A000009L,
		0x02000401L, 0x0A000401L, 0x02000409L, 0x0A000409L,
		0x02020001L, 0x0A020001L, 0x02020009L, 0x0A020009L,
		0x02020401L, 0x0A020401L, 0x02020409L, 0x0A020409L,
	},
	{
		/* for D bits (numbered as per FIPS 46) 16 17 18 19 20 21 */
		// 0000 0001 0010 1000 0000 0011 0001 0000
		0x00000000L, 0x00000100L, 0x00080000L, 0x00080100L,
		0x01000000L, 0x01000100L, 0x01080000L, 0x01080100L,
		0x00000010L, 0x00000110L, 0x00080010L, 0x00080110L,
		0x01000010L, 0x01000110L, 0x01080010L, 0x01080110L,
		0x00200000L, 0x00200100L, 0x00280000L, 0x00280100L,
		0x01200000L, 0x01200100L, 0x01280000L, 0x01280100L,
		0x00200010L, 0x00200110L, 0x00280010L, 0x00280110L,
		0x01200010L, 0x01200110L, 0x01280010L, 0x01280110L,
		0x00000200L, 0x00000300L, 0x00080200L, 0x00080300L,
		0x01000200L, 0x01000300L, 0x01080200L, 0x01080300L,
		0x00000210L, 0x00000310L, 0x00080210L, 0x00080310L,
		0x01000210L, 0x01000310L, 0x01080210L, 0x01080310L,
		0x00200200L, 0x00200300L, 0x00280200L, 0x00280300L,
		0x01200200L, 0x01200300L, 0x01280200L, 0x01280300L,
		0x00200210L, 0x00200310L, 0x00280210L, 0x00280310L,
		0x01200210L, 0x01200310L, 0x01280210L, 0x01280310L,
	},
	{
		/* for D bits (numbered as per FIPS 46) 22 23 24 25 27 28 */
		// 0000 0100 0000 0100 0010 1000 0010 0010
		0x00000000L, 0x04000000L, 0x00040000L, 0x04040000L,
		0x00000002L, 0x04000002L, 0x00040002L, 0x04040002L,
		0x00002000L, 0x04002000L, 0x00042000L, 0x04042000L,
		0x00002002L, 0x04002002L, 0x00042002L, 0x04042002L,
		0x00000020L, 0x04000020L, 0x00040020L, 0x04040020L,
		0x00000022L, 0x04000022L, 0x00040022L, 0x04040022L,
		0x00002020L, 0x04002020L, 0x00042020L, 0x04042020L,
		0x00002022L, 0x04002022L, 0x00042022L, 0x04042022L,
		0x00000800L, 0x04000800L, 0x00040800L, 0x04040800L,
		0x00000802L, 0x04000802L, 0x00040802L, 0x04040802L,
		0x00002800L, 0x04002800L, 0x00042800L, 0x04042800L,
		0x00002802L, 0x04002802L, 0x00042802L, 0x04042802L,
		0x00000820L, 0x04000820L, 0x00040820L, 0x04040820L,
		0x00000822L, 0x04000822L, 0x00040822L, 0x04040822L,
		0x00002820L, 0x04002820L, 0x00042820L, 0x04042820L,
		0x00002822L, 0x04002822L, 0x00042822L, 0x04042822L,
	},
};

#define PERM_OP(a, b, t, n, m)          \
  do {                                  \
    (t) = ((((a) >> (n)) ^ (b)) & (m)); \
    (b) ^= (t);                         \
    (a) ^= ((t) << (n));                \
  } while (0)

#define HPERM_OP(a,t,n,m) ((t)=((((a)<<(16-(n)))^(a))&(m)),\
        (a)=(a)^(t)^(t>>(16-(n))))

static uint32_t find_c(uint32_t s) {
	uint32_t c = 0;
	uint32_t temp = s & 0x20090830L;
	for (int i = 0; i < 64; i++) {
		if (des_skb[0][i] == temp) {
			c |= i;
			break;
		}
	}
	temp = s & 0x12202404L;
	for (int i = 0; i < 64; i++) {
		if (des_skb[1][i] == temp) {
			c |= ((i & 0x03) | ((i & 0x3c) << 1)) << 6;
			break;
		}
	}
	temp = s & 0x09040203L;
	for (int i = 0; i < 64; i++) {
		if (des_skb[2][i] == temp) {
			c |= ((i & 0x0f) | ((i & 0x30) << 1)) << 13;
			break;
		}
	}
	temp = s & 0x04121108L;
	for (int i = 0; i < 64; i++) {
		if (des_skb[3][i] == temp) {
			c |= ((i & 0x01) | ((i & 0x06) << 1) | ((i & 0x38) << 2)) << 20;
			break;
		}
	}
	return c;
}

static uint32_t find_d(uint32_t t) {
	uint32_t d = 0;
	uint32_t temp = t & 0x30111004L;
	for (int i = 0; i < 64; i++) {
		if (des_skb[4][i] == temp) {
			d |= i;
			break;
		}
	}
	temp = t & 0x0A020409L;
	for (int i = 0; i < 64; i++) {
		if (des_skb[5][i] == temp) {
			d |= ((i & 0x03) | ((i & 0x3c) << 1)) << 7;
			break;
		}
	}
	temp = t & 0x01280310L;
	for (int i = 0; i < 64; i++) {
		if (des_skb[6][i] == temp) {
			d |= i << 15;
			break;
		}
	}
	temp = t & 0x04042822L;
	for (int i = 0; i < 64; i++) {
		if (des_skb[7][i] == temp) {
			d |= ((i & 0x0f) | ((i & 0x30) << 1)) << 21;
			break;
		}
	}
	return d;
}

static bool check_entropy(const uint32_t* rk) {
	for (int i = 0; i < 6; i++) {
		if (!rk[i] || !~rk[i])
			return false;
	}
	for (int i = 0; i < 6; i++) {
		if ((rk[i] & 0xff00ff00) == 0 || (rk[i] & 0xff00ff00) == 0xff00ff00)
			return false;
		if ((rk[i] & 0xffff0000) == 0 || (rk[i] & 0xffff0000) == 0xffff0000)
			return false;
		if ((rk[i] & 0x0000ffff) == 0 || (rk[i] & 0x0000ffff) == 0x0000ffff)
			return false;
		if ((rk[i] & 0x00ff00ff) == 0 || (rk[i] & 0x00ff00ff) == 0x00ff00ff)
			return false;
	}
	return true;
}

static uint64_t des_key_detect(const uint32_t* rk) {
	const static uint32_t __shift[] = { 1, 2, 4, 6, 8, 10, 12, 14, 15, 17, 19, 21, 23, 25, 27, 28 };

	if (!check_entropy(rk))
		return 0;
	uint32_t orig_c = 0, orig_d = 0;
	for (int i = 0; i < 2; i++) {
		int t2_1 = _lrotr(rk[i * 2], 2);
		int t2_2 = _lrotr(rk[i * 2 + 1], 6);
		uint32_t t = ((t2_1 >> 16L) & 0x0000ffff) | (t2_2 & 0xffff0000);
		uint32_t s = (t2_2 << 16L) | (t2_1 & 0x0000ffff);
		uint32_t c = find_c(s);
		uint32_t d = find_d(t);
		c = ((c << __shift[i]) | (c >> (28 - __shift[i])));
		d = ((d << __shift[i]) | (d >> (28 - __shift[i])));
		c &= 0x0fffffffU;
		d &= 0x0fffffffU;
		orig_c |= c;
		orig_d |= d;
	}
	//if (!check_entropy(orig_c)) return 0;

	//if (!check_entropy(orig_d)) return 0;
	
	for (int i = 0; i < 16; i++) {
		uint32_t c = ((orig_c >> __shift[i]) | (orig_c << (28 - __shift[i])));
		uint32_t d = ((orig_d >> __shift[i]) | (orig_d << (28 - __shift[i])));
		c &= 0x0fffffffL;
		d &= 0x0fffffffL;

		uint32_t s = des_skb[0][(c) & 0x3f] |
			des_skb[1][((c >> 6L) & 0x03) | ((c >> 7L) & 0x3c)] |
			des_skb[2][((c >> 13L) & 0x0f) | ((c >> 14L) & 0x30)] |
			des_skb[3][((c >> 20L) & 0x01) | ((c >> 21L) & 0x06) |
			((c >> 22L) & 0x38)];
		uint32_t t = des_skb[4][(d) & 0x3f] |
			des_skb[5][((d >> 7L) & 0x03) | ((d >> 8L) & 0x3c)] |
			des_skb[6][(d >> 15L) & 0x3f] |
			des_skb[7][((d >> 21L) & 0x0f) | ((d >> 22L) & 0x30)];

		/* table contained 0213 4657 */
		uint32_t t2 = ((t << 16L) | (s & 0x0000ffffL)) & 0xffffffffL;
		if (rk[i * 2] != _lrotr(t2, 30))
			return 0;

		t2 = ((s >> 16L) | (t & 0xffff0000L));
		if (rk[i * 2 + 1] != _lrotr(t2, 26))
			return 0;
	}

	orig_c |= (orig_d & 0x0f000000L) << 4L;
	orig_d = ((orig_d & 0x000000ffL) << 16L) | (orig_d & 0x0000ff00L) |
		((orig_d & 0x00ff0000L) >> 16L);
	uint32_t t;
	PERM_OP(orig_d, orig_c, t, 1, 0x55555555L);
	PERM_OP(orig_c, orig_d, t, 8, 0x00ff00ffL);
	PERM_OP(orig_d, orig_c, t, 1, 0x55555555L);
	HPERM_OP(orig_d, t, -2, 0xcccc0000L);
	HPERM_OP(orig_c, t, -2, 0xcccc0000L);
	PERM_OP(orig_d, orig_c, t, 4, 0x0f0f0f0fL);
	return ((uint64_t)orig_d << 32) | (uint64_t)orig_c;
}

static void find_keys(uint32_t pid)
{
	printf("Searching PID %u ...\n", pid);

	if (!os_process_begin(pid))
	{
		printf("Failed to open process\n");
		return;
	}

	uint8_t buffer[8192] = { 0, };
	uint32_t avail = 16;

	uint64_t region = 0;
	uint64_t region_size = 0;
	uint64_t size = 0;

	uint64_t addr = 0;

	clock_t t0 = clock();
	uint64_t total = 0;

	uint64_t next_size;
	uint64_t next = os_process_next(&next_size);
	addr = region = next;
	size = region_size = next_size;

	for (;;)
	{
		if (size == 0)
		{
			uint64_t next_size;
			uint64_t next = os_process_next(&next_size);
			if (next == 0)
			{
				break;
			}

			if (region + region_size != next)
			{
				avail = 0;
			}

			addr = region = next;
			size = region_size = next_size;
		}

		uint32_t read = sizeof(buffer) - avail;
		if (read > size) read = (uint32_t)size;

		read = os_process_read(region, buffer + avail, read);
		if (read == 0)
		{
			avail = 0;
			size = 0;
			continue;
		}
		total += read;
		region += read;
		size -= read;
		avail += read;

		uint32_t offset = 0;
		if (avail >= 32 * sizeof(uint32_t))
		{
			while (offset <= avail - 32 * sizeof(uint32_t))
			{
				uint64_t key;
				if (key = des_key_detect((const uint32_t*)&buffer[offset])) {
					printf("[%p] Found DES key: %016llx\n", (void*)addr, _byteswap_uint64(key));

					offset += 32 * sizeof(uint32_t);
					addr += 32 * sizeof(uint32_t);
				}
				else {
					offset += 4;
					addr += 4;
				}
			}

			avail -= offset;
		}

		memmove(buffer, buffer + offset, avail);
	}

	clock_t t1 = clock();
	double time = (double)(t1 - t0) / CLOCKS_PER_SEC;
	const double MB = 1024.0 * 1024.0;
	printf("Processed %.2f MB, speed = %.2f MB/s\n", total / MB, total / MB / time);

	os_process_end();
}

static void find_keys(char* arg) {
	FILE* fp;
	size_t read_len = 0;
	size_t remain = 0;
	size_t i = 0;
	size_t total = 0;
	uint8_t buffer[8192] = { 0, };
	if ((fp = fopen(arg, "rb")) == NULL)
	{
		printf("파일 오픈 에러 '%s'\n", arg);
		return;
	}
	while (!feof(fp)) {
		read_len = fread(buffer + remain, sizeof(char), sizeof(buffer) - remain, fp);
		for (i = 0; i <= read_len + remain - 32 * sizeof(uint32_t);) {
			uint64_t key;
			if (key = des_key_detect((const uint32_t*)&buffer[i])) {
				printf("[%p] Found DES key: %016llx\n", (void*)(total + i), _byteswap_uint64(key));

				i += 32 * sizeof(uint32_t);
			}
			else {
				i += sizeof(uint32_t);
			}
		}
		total += read_len;
		remain = read_len + remain - i;
		memmove(buffer, buffer + i, remain);
	}
	fclose(fp);
	printf("Processed %d Byte\n", (int)total);
}

static void printHelp() {
	printf("Usage: aria-finder [ {-i/--pid} [PID1 PID2 ...]]\n");
	printf("                   [ {-p/--process-name} [PNAME1 PNAME2 ...]]\n");
	printf("                   [ {-f/--file-name} [FNAME1 FNAME2 ...]]\n");
}

int __cdecl main(int argc, char* argv[])
{
	enum ArgType
	{
		PID, PNAME, FNAME, NONE
	} argType = PID;
	int pid = 0;
	if (argc < 2)
	{
		printHelp();
		return EXIT_FAILURE;
	}

	os_startup();

	for (int i = 1; i < argc; i++) {
		if (argv[i][0] == '-') {
			switch (argv[i][1]) {
			case '-':
				if (!strcmp(argv[i] + 2, "pid"))
					argType = PID;
				else if (!strcmp(argv[i] + 2, "process-name"))
					argType = PNAME;
				else if (!strcmp(argv[i] + 2, "file-name"))
					argType = FNAME;
				else {
					printHelp();
					return EXIT_FAILURE;
				}
				break;
			case 'i':
				argType = PID;
				break;
			case 'p':
				argType = PNAME;
				break;
			case 'f':
				argType = FNAME;
				break;
			default:
				printHelp();
				return EXIT_FAILURE;
			}
		}
		else {
			switch (argType) {
			case PID:
				printf("search key in (PID) %s\n", argv[i]);
				pid = atoi(argv[i]);
				find_keys(pid);
				break;
			case PNAME:
				printf("search key in (process) %s\n", argv[i]);
				if (os_enum_start()) {
					while (pid = os_enum_next(argv[i])) {
						find_keys(pid);
					}
					os_enum_end();
				}
				break;
			case FNAME:
				printf("search key in %s\n", argv[i]);
				find_keys(argv[i]);
				break;
			default:
				printHelp();
				return EXIT_FAILURE;
			}
		}
	}

	printf("Done!\n\n");
	return EXIT_SUCCESS;
}
