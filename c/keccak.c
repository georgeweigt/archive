/* Keccak-256 coded to look like the FIPS standard

Regarding bit indexing of strings in FIPS PUB 202

In the standard, X[i] indicates the i'th bit of string X (see page 5)

Given a 200 byte string X, we have

X[0] is the LS bit of the first byte (lowest memory address)

X[1599] is the MS bit of the last byte (highest memory address)

Regarding Keccak-256 (see Table 3 on page 22 of FIPS PUB 202)

Rate r = 1088 bits (136 bytes)

Capacity c = 512 bits ( 64 bytes)

*/

#include <stdio.h>
#include <stdint.h>
#include <string.h>

#define RATE 136

#define A(x,y,z) A[320 * (x) + 64 * (y) + (z)]
#define Aprime(x,y,z) Aprime[320 * (x) + 64 * (y) + (z)]

uint8_t RC[64][24]; // round constants

uint8_t *
theta(uint8_t *A)
{
	int x, y, z;
	static uint8_t Aprime[1600], C[5][64], D[5][64];

	for (x = 0; x < 5; x++)
		for (z = 0; z < 64; z++)
			C[x][z] = A(x,0,z) ^ A(x,1,z) ^ A(x,2,z) ^ A(x,3,z) ^ A(x,4,z);

	for (x = 0; x < 5; x++)
		for (z = 0; z < 64; z++)
			D[x][z] = C[(5 + x - 1) % 5][z] ^ C[(x + 1) % 5][(64 + z - 1) % 64];

	for (x = 0; x < 5; x++)
		for (y = 0; y < 5; y++)
			for (z = 0; z < 64; z++)
				Aprime(x,y,z) = A(x,y,z) ^ D[x][z];

	return Aprime;
}

uint8_t *
rho(uint8_t *A)
{
	int t, u, x, y, z;
	static uint8_t Aprime[1600];

	for (z = 0; z < 64; z++)
		Aprime(0,0,z) = A(0,0,z);

	x = 1;
	y = 0;

	for (t = 0; t < 24; t++) {
		for (z = 0; z < 64; z++)
			Aprime(x,y,z) = A(x,y,(320 + z - (t + 1) * (t + 2) / 2) % 64);
		u = y;
		y = (2 * x + 3 * y) % 5;
		x = u;
	}

	return Aprime;
}

uint8_t *
pi(uint8_t *A)
{
	int x, y, z;
	static uint8_t Aprime[1600];

	for (x = 0; x < 5; x++)
		for (y = 0; y < 5; y++)
			for (z = 0; z < 64; z++)
				Aprime(x,y,z) = A((x + 3 * y) % 5,x,z);

	return Aprime;
}

uint8_t *
chi(uint8_t *A)
{
	int x, y, z;
	static uint8_t Aprime[1600];

	for (x = 0; x < 5; x++)
		for (y = 0; y < 5; y++)
			for (z = 0; z < 64; z++)
				Aprime(x,y,z) = A(x,y,z) ^ ((A((x + 1) % 5,y,z) ^ 1) & A((x + 2) % 5,y,z));

	return Aprime;
}

uint8_t
rc(int t)
{
	int i, R;

	if (t % 255 == 0)
		return 1;

	R = 1;

	for (i = 1; i <= t % 255; i++) {
		R <<= 1;
		if (R & 0x100)
			R ^= 0x171;
	}

	return R & 1;
}

uint8_t *
iota(uint8_t *A, int ir)
{
	int z;

	for (z = 0; z < 64; z++)
		A(0,0,z) ^= RC[z][ir];

	return A;
}

uint8_t *
Rnd(uint8_t *A, int ir)
{
	return iota(chi(pi(rho(theta(A)))), ir);
}

uint8_t mask[8] = {1,2,4,8,0x10,0x20,0x40,0x80};

void
Keccak(uint8_t *S)
{
	int ir, k, x, y, z;
	static uint8_t a[1600], *A = a;

	// convert S to A

	memset(A, 0, 1600);

	for (x = 0; x < 5; x++)
		for (y = 0; y < 5; y++)
			for (z = 0; z < 64; z++) {
				k = 64 * (5 * y + x) + z;
				if (S[k / 8] & mask[k % 8])
					A(x,y,z) = 1;
			}

	for (ir = 0; ir < 24; ir++)
		A = Rnd(A, ir);

	// convert A to S

	memset(S, 0, 200);

	for (x = 0; x < 5; x++)
		for (y = 0; y < 5; y++)
			for (z = 0; z < 64; z++)
				if (A(x,y,z)) {
					k = 64 * (5 * y + x) + z;
					S[k / 8] |= mask[k % 8];
				}
}

uint8_t *
sponge(uint8_t *N, int len)
{
	int i, j, k, n;
	static uint8_t S[200]; // 1600 bits

	memset(S, 0, 200);

	n = len / RATE; // number of full blocks

	for (i = 0; i < n; i++) {
		for (j = 0; j < RATE; j++)
			S[j] ^= N[RATE * i + j];
		Keccak(S);
	}

	// pad last block

	k = len % RATE;

	for (i = 0; i < k; i++)
		S[i] ^= N[RATE * n + i];

	S[k] ^= 0x01;
	S[RATE - 1] ^= 0x80;

	Keccak(S);

	return S;
}

void
keccak256(uint8_t *outbuf, uint8_t *inbuf, int inbuflen)
{
	uint8_t *S = sponge(inbuf, inbuflen);
	memcpy(outbuf, S, 32);
}

char *
keccak256str(uint8_t *buf, int len)
{
	int i;
	uint8_t *S;
	static char Z[65];

	S = sponge(buf, len);

	for (i = 0; i < 32; i++)
		sprintf(Z + 2 * i, "%02x", S[i]);

	return Z;
}

void
test_keccak256(void)
{
	int err;
	char *Z;
	static uint8_t buf[RATE + 1];

	printf("Testing keccak256\n");

	memset(buf, 'a', sizeof buf);

	Z = keccak256str(NULL, 0);
	err = strcmp(Z, "c5d2460186f7233c927e7db2dcc703c0e500b653ca82273b7bfad8045d85a470");
	printf("%s %s\n", Z, err ? "err" : "ok");

	Z = keccak256str((uint8_t *) "hello", 5);
	err = strcmp(Z, "1c8aff950685c2ed4bc3174f3472287b56d9517b9c948127319a09a7a36deac8");
	printf("%s %s\n", Z, err ? "err" : "ok");

	Z = keccak256str(buf, RATE - 1);
	err = strcmp(Z, "34367dc248bbd832f4e3e69dfaac2f92638bd0bbd18f2912ba4ef454919cf446");
	printf("%s %s\n", Z, err ? "err" : "ok");

	Z = keccak256str(buf, RATE);
	err = strcmp(Z, "a6c4d403279fe3e0af03729caada8374b5ca54d8065329a3ebcaeb4b60aa386e");
	printf("%s %s\n", Z, err ? "err" : "ok");

	Z = keccak256str(buf, RATE + 1);
	err = strcmp(Z, "d869f639c7046b4929fc92a4d988a8b22c55fbadb802c0c66ebcd484f1915f39");
	printf("%s %s\n", Z, err ? "err" : "ok");
}

int
main(void)
{
	int i, j;

	for (i = 0; i < 24; i++)
		for (j = 0; j < 7; j++)
			RC[(1 << j) - 1][i] = rc(j + 7 * i);

	test_keccak256();
}
