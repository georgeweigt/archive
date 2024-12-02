/* This is a plain and simple unsigned big number library.

The "main" function runs a self test.

To compile and run:

	gcc bignum.c
	./a.out
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define MLENGTH(p) (p)[-1]
#define MZERO(p) (MLENGTH(p) == 1 && (p)[0] == 0)
#define MEQUAL(p, n) (MLENGTH(p) == 1 && (p)[0] == (n))

void verify(uint32_t *v);
uint32_t * madd(uint32_t *u, uint32_t *v);
uint32_t * msub(uint32_t *u, uint32_t *v);
uint32_t * mmul(uint32_t *u, uint32_t *v);
uint32_t * mdiv(uint32_t *u, uint32_t *v);
void mmod(uint32_t *u, uint32_t *v);
uint32_t * mpow(uint32_t *u, uint32_t *v);
void mshr(uint32_t *u);
int mcmp(uint32_t *u, uint32_t *v);
uint32_t * mint(uint32_t n);
uint32_t * mnew(int n);
void mfree(uint32_t *u);
uint32_t * mcopy(uint32_t *u);
void mnorm(uint32_t *u);
uint32_t * mscan(char *s);
char * mstr(uint32_t *u);
int mdivby1billion(uint32_t *u);
uint32_t * mgcd(uint32_t *u, uint32_t *v);
uint32_t * mroot(uint32_t *a, uint32_t *n);
void msetbit(uint32_t *x, uint32_t k);
void mclrbit(uint32_t *x, uint32_t k);
uint32_t * mmodinv(uint32_t *a, uint32_t *p);

uint32_t w[7] = {
	0x00000000,
	0x00000001,
	0x7fffffff,
	0x80000000,
	0x80000001,
	0xfffffffe,
	0xffffffff,
};

int mcount;

int
main(int argc, char *argv[])
{
	int a, b, c, d, e, f, g, h;
	uint32_t v[8];

	for (a = 0; a < 7; a++) { v[0] = w[a];
	for (b = 0; b < 7; b++) { v[1] = w[b];
	for (c = 0; c < 7; c++) { v[2] = w[c]; putchar('.'); fflush(stdout); // progress indicator
	for (d = 0; d < 7; d++) { v[3] = w[d];
	for (e = 0; e < 7; e++) { v[4] = w[e];
	for (f = 0; f < 7; f++) { v[5] = w[f];
	for (g = 0; g < 7; g++) { v[6] = w[g];
	for (h = 0; h < 7; h++) { v[7] = w[h];
		verify(v);
	}}}}}}}}

	if (mcount) {
		printf("memory leak err\n");
		exit(1);
	}

	printf("ok\n");
}

void
verify(uint32_t *v)
{
	uint32_t *a, *b, *q, *r, *s, *t;

	a = mnew(4);
	b = mnew(4);

	a[0] = v[0];
	a[1] = v[1];
	a[2] = v[2];
	a[3] = v[3];

	b[0] = v[4];
	b[1] = v[5];
	b[2] = v[6];
	b[3] = v[7];

	mnorm(a);
	mnorm(b);

	q = mdiv(a, b);

	if (q == NULL) {
		if (MLENGTH(b) == 1 && b[0] == 0) {
			mfree(a);
			mfree(b);
			return; // divide by zero
		}
		printf("err line %d\n", __LINE__);
		exit(1);
	}

	r = mcopy(a);

	mmod(r, b);

	// check r < b

	if (mcmp(r, b) >= 0) {
		printf("err line %d\n", __LINE__);
		exit(1);
	}

	// check a = b * q + r

	t = mmul(b, q);
	s = madd(t, r);

	if (mcmp(a, s) != 0) {
		printf("err line %d\n", __LINE__);
		exit(1);
	}

	// check a - r = b * q

	mfree(s);
	s = msub(a, r);

	if (mcmp(s, t) != 0) {
		printf("err line %d\n", __LINE__);
		exit(1);
	}

	mfree(a);
	mfree(b);
	mfree(q);
	mfree(r);
	mfree(s);
	mfree(t);
}

// returns u + v

uint32_t *
madd(uint32_t *u, uint32_t *v)
{
	int i, nu, nv, nw;
	uint64_t t;
	uint32_t *w;
	nu = MLENGTH(u);
	nv = MLENGTH(v);
	if (nu > nv)
		nw = nu + 1;
	else
		nw = nv + 1;
	w = mnew(nw);
	for (i = 0; i < nu; i++)
		w[i] = u[i];
	for (i = nu; i < nw; i++)
		w[i] = 0;
	t = 0;
	for (i = 0; i < nv; i++) {
		t += (uint64_t) w[i] + v[i];
		w[i] = t;
		t >>= 32;
	}
	for (i = nv; i < nw; i++) {
		t += w[i];
		w[i] = t;
		t >>= 32;
	}
	mnorm(w);
	return w;
}

// returns u - v

uint32_t *
msub(uint32_t *u, uint32_t *v)
{
	int i, nu, nv, nw;
	uint64_t t;
	uint32_t *w;
	nu = MLENGTH(u);
	nv = MLENGTH(v);
	if (nu > nv)
		nw = nu;
	else
		nw = nv;
	w = mnew(nw);
	for (i = 0; i < nu; i++)
		w[i] = u[i];
	for (i = nu; i < nw; i++)
		w[i] = 0;
	t = 0;
	for (i = 0; i < nv; i++) {
		t += (uint64_t) w[i] - v[i];
		w[i] = t;
		t = (long long) t >> 32; // cast to extend sign
	}
	for (i = nv; i < nw; i++) {
		t += w[i];
		w[i] = t;
		t = (long long) t >> 32; // cast to extend sign
	}
	mnorm(w);
	return w;
}

// returns u * v

uint32_t *
mmul(uint32_t *u, uint32_t *v)
{
	int i, j, nu, nv, nw;
	uint64_t t;
	uint32_t *w;
	nu = MLENGTH(u);
	nv = MLENGTH(v);
	nw = nu + nv;
	w = mnew(nw);
	for (i = 0; i < nu; i++)
		w[i] = 0;
	for (j = 0; j < nv; j++) {
		t = 0;
		for (i = 0; i < nu; i++) {
			t += (uint64_t) u[i] * v[j] + w[i + j];
			w[i + j] = t;
			t >>= 32;
		}
		w[i + j] = t;
	}
	mnorm(w);
	return w;
}

// returns floor(u / v)

uint32_t *
mdiv(uint32_t *u, uint32_t *v)
{
	int i, k, nu, nv;
	uint32_t *q, qhat, *w;
	uint64_t a, b, t;
	mnorm(u);
	mnorm(v);
	if (MLENGTH(v) == 1 && v[0] == 0)
		return NULL; // v = 0
	nu = MLENGTH(u);
	nv = MLENGTH(v);
	k = nu - nv;
	if (k < 0) {
		q = mnew(1);
		q[0] = 0;
		return q; // u < v, return zero
	}
	u = mcopy(u);
	q = mnew(k + 1);
	w = mnew(nv + 1);
	b = v[nv - 1];
	do {
		q[k] = 0;
		while (nu >= nv + k) {
			// estimate 32-bit partial quotient
			a = u[nu - 1];
			if (nu > nv + k)
				a = a << 32 | u[nu - 2];
			if (a < b)
				break;
			qhat = a / (b + 1);
			if (qhat == 0)
				qhat = 1;
			// w = qhat * v
			t = 0;
			for (i = 0; i < nv; i++) {
				t += (uint64_t) qhat * v[i];
				w[i] = t;
				t >>= 32;
			}
			w[nv] = t;
			// u = u - w
			t = 0;
			for (i = k; i < nu; i++) {
				t += (uint64_t) u[i] - w[i - k];
				u[i] = t;
				t = (long long) t >> 32; // cast to extend sign
			}
			if (t) {
				// u is negative, restore u
				t = 0;
				for (i = k; i < nu; i++) {
					t += (uint64_t) u[i] + w[i - k];
					u[i] = t;
					t >>= 32;
				}
				break;
			}
			q[k] += qhat;
			mnorm(u);
			nu = MLENGTH(u);
		}
	} while (--k >= 0);
	mnorm(q);
	mfree(u);
	mfree(w);
	return q;
}

// u = u mod v

void
mmod(uint32_t *u, uint32_t *v)
{
	int i, k, nu, nv;
	uint32_t qhat, *w;
	uint64_t a, b, t;
	mnorm(u);
	mnorm(v);
	if (MLENGTH(v) == 1 && v[0] == 0)
		return; // v = 0
	nu = MLENGTH(u);
	nv = MLENGTH(v);
	k = nu - nv;
	if (k < 0)
		return; // u < v
	w = mnew(nv + 1);
	b = v[nv - 1];
	do {
		while (nu >= nv + k) {
			// estimate 32-bit partial quotient
			a = u[nu - 1];
			if (nu > nv + k)
				a = a << 32 | u[nu - 2];
			if (a < b)
				break;
			qhat = a / (b + 1);
			if (qhat == 0)
				qhat = 1;
			// w = qhat * v
			t = 0;
			for (i = 0; i < nv; i++) {
				t += (uint64_t) qhat * v[i];
				w[i] = t;
				t >>= 32;
			}
			w[nv] = t;
			// u = u - w
			t = 0;
			for (i = k; i < nu; i++) {
				t += (uint64_t) u[i] - w[i - k];
				u[i] = t;
				t = (long long) t >> 32; // cast to extend sign
			}
			if (t) {
				// u is negative, restore u
				t = 0;
				for (i = k; i < nu; i++) {
					t += (uint64_t) u[i] + w[i - k];
					u[i] = t;
					t >>= 32;
				}
				break;
			}
			mnorm(u);
			nu = MLENGTH(u);
		}
	} while (--k >= 0);
	mfree(w);
}

// returns u ** v

uint32_t *
mpow(uint32_t *u, uint32_t *v)
{
	uint32_t *t, *w;
	u = mcopy(u);
	v = mcopy(v);
	// w = 1
	w = mnew(1);
	w[0] = 1;
	for (;;) {
		if (v[0] & 1) {
			// w = w * u
			t = mmul(w, u);
			mfree(w);
			w = t;
		}
		// v = v >> 1
		mshr(v);
		// v = 0?
		if (MLENGTH(v) == 1 && v[0] == 0)
			break;
		// u = u * u
		t = mmul(u, u);
		mfree(u);
		u = t;
	}
	mfree(u);
	mfree(v);
	return w;
}

// u = u >> 1

void
mshr(uint32_t *u)
{
	int i;
	for (i = 0; i < MLENGTH(u) - 1; i++) {
		u[i] >>= 1;
		if (u[i + 1] & 1)
			u[i] |= 0x80000000;
	}
	u[i] >>= 1;
	mnorm(u);
}

// compare u and v

int
mcmp(uint32_t *u, uint32_t *v)
{
	int i;
	mnorm(u);
	mnorm(v);
	if (MLENGTH(u) < MLENGTH(v))
		return -1;
	if (MLENGTH(u) > MLENGTH(v))
		return 1;
	for (i = MLENGTH(u) - 1; i >= 0; i--) {
		if (u[i] < v[i])
			return -1;
		if (u[i] > v[i])
			return 1;
	}
	return 0; // u = v
}

// convert unsigned to bignum

uint32_t *
mint(uint32_t n)
{
	uint32_t *p;
	p = mnew(1);
	p[0] = n;
	return p;
}

uint32_t *
mnew(int n)
{
	uint32_t *u;
	u = (uint32_t *) malloc((n + 1) * sizeof (uint32_t));
	if (u == NULL) {
		printf("malloc kaput\n");
		exit(1);
	}
	mcount++;
	*u = n;
	return u + 1;
}

void
mfree(uint32_t *u)
{
	free(u - 1);
	mcount--;
}

uint32_t *
mcopy(uint32_t *u)
{
	int i;
	uint32_t *v;
	v = mnew(MLENGTH(u));
	for (i = 0; i < MLENGTH(u); i++)
		v[i] = u[i];
	return v;
}

// remove leading zeroes

void
mnorm(uint32_t *u)
{
	while (MLENGTH(u) > 1 && u[MLENGTH(u) - 1] == 0)
		MLENGTH(u)--;
}

// convert string to bignum (9 decimal digits fits in 32 bits)

uint32_t *
mscan(char *s)
{
	int i, k, len;
	uint32_t *a, *b, *t;
	a = mint(0);
	len = (int) strlen(s);
	if (len == 0)
		return a;
	k = len % 9;
	if (k == 0)
		k = 9;
	for (i = 0; i < k; i++)
		a[0] = 10 * a[0] + s[i] - '0';
	if (k == len)
		return a;
	b = mint(0);
	while (k < len) {
		b[0] = 1000000000; // 10^9
		t = mmul(a, b);
		mfree(a);
		a = t;
		b[0] = 0;
		for (i = 0; i < 9; i++)
			b[0] = 10 * b[0] + s[k++] - '0';
		t = madd(a, b);
		mfree(a);
		a = t;
	}
	mfree(b);
	return a;
}

// convert bignum to string (returned value points to static buffer)

char *
mstr(uint32_t *u)
{
	int i, k, n, r;
	static char *buf;
	static int len;

	// estimate string length

	// note that 0xffffffff -> 000000004 294967295

	// hence space for 8 leading zeroes is required

	n = 10 * MLENGTH(u) + 10;

	n = 1000 * (n / 1000 + 1);

	if (n > len) {
		if (buf)
			free(buf);
		buf = malloc(n);
		if (buf == NULL)
			exit(1);
		len = n;
	}

	u = mcopy(u);

	k = len - 1;
	buf[k] = '\0'; // string terminator

	for (;;) {
		r = mdivby1billion(u);
		for (i = 0; i < 9; i++) {
			buf[--k] = r % 10 + '0';
			r /= 10;
		}
		if (MZERO(u))
			break;
	}

	mfree(u);

	// remove leading zeroes

	while (buf[k] == '0' && buf[k + 1])
		k++;

	return buf + k;
}

// returns remainder, quotient returned in u

int
mdivby1billion(uint32_t *u)
{
	int i;
	uint64_t r = 0;
	for (i = MLENGTH(u) - 1; i >= 0; i--) {
		r = r << 32 | u[i];
		u[i] = (uint32_t) (r / 1000000000);
		r -= (uint64_t) 1000000000 * u[i];
	}
	mnorm(u);
	return (int) r;
}

uint32_t *
mgcd(uint32_t *u, uint32_t *v)
{
	int i, k, n, sign;
	uint32_t *t;

	if (MZERO(u)) {
		t = mcopy(v);
		return t;
	}

	if (MZERO(v)) {
		t = mcopy(u);
		return t;
	}

	u = mcopy(u);
	v = mcopy(v);

	k = 0;

	while ((u[0] & 1) == 0 && (v[0] & 1) == 0) {
		mshr(u);
		mshr(v);
		k++;
	}

	if (u[0] & 1) {
		t = mcopy(v);
		sign = -1;
	} else {
		t = mcopy(u);
		sign = 1;
	}

	for (;;) {

		while ((t[0] & 1) == 0)
			mshr(t);

		if (sign == 1) {
			mfree(u);
			u = mcopy(t);
		} else {
			mfree(v);
			v = mcopy(t);
		}

		mfree(t);

		if (mcmp(u, v) < 0) {
			t = msub(v, u);
			sign = -1;
		} else {
			t = msub(u, v);
			sign = 1;
		}

		if (MZERO(t)) {
			mfree(t);
			mfree(v);
			n = (k / 32) + 1;
			v = mnew(n);
			for (i = 0; i < n; i++)
				v[i] = 0;
			msetbit(v, k);
			t = mmul(u, v);
			mfree(u);
			mfree(v);
			return t;
		}
	}
}

// returns NULL if not perfect root, otherwise returns a^(1/n)

uint32_t *
mroot(uint32_t *a, uint32_t *n)
{
	int i, j, k;
	uint32_t *b, *c, m;

	if (MLENGTH(n) > 1 || n[0] == 0)
		return NULL;

	// k is bit length of a

	k = 32 * (MLENGTH(a) - 1);

	m = a[MLENGTH(a) - 1];

	while (m) {
		m >>= 1;
		k++;
	}

	if (k == 0)
		return mint(0);

	// initial guess of index of ms bit in result

	k = (k - 1) / n[0];

	j = k / 32 + 1; // k is bit index, not number of bits

	b = mnew(j);

	for (i = 0; i < j; i++)
		b[i] = 0;

	while (k >= 0) {
		msetbit(b, k);
		mnorm(b);
		c = mpow(b, n);
		switch (mcmp(c, a)) {
		case -1:
			break;
		case 0:
			mfree(c);
			return b;
		case 1:
			mclrbit(b, k);
			break;
		}
		mfree(c);
		k--;
	}

	mfree(b);

	return NULL;
}

void
msetbit(uint32_t *x, uint32_t k)
{
	x[k / 32] |= 1 << (k % 32);
}

void
mclrbit(uint32_t *x, uint32_t k)
{
	x[k / 32] &= ~(1 << (k % 32));
}

// See 'Mathematical routines for the NIST prime elliptic curves'

// Returns (1 / a) mod p

uint32_t *
mmodinv(uint32_t *a, uint32_t *p)
{
	uint32_t *k, *r, *u, *v, *t, *x1, *x2;
	u = mcopy(a);
	v = mcopy(p);
	x1 = mint(1);
	x2 = mint(0);
	while (!MEQUAL(u, 1) && !MEQUAL(v, 1)) {
		while ((u[0] & 1) == 0) {
			mshr(u);
			if (x1[0] & 1) {
				t = madd(x1, p);
				mfree(x1);
				x1 = t;
			}
			mshr(x1);
		}
		while ((v[0] & 1) == 0) {
			mshr(v);
			if (x2[0] & 1) {
				t = madd(x2, p);
				mfree(x2);
				x2 = t;
			}
			mshr(x2);
		}
		if (mcmp(u, v) >= 0) {
			t = msub(u, v);
			mfree(u);
			u = t;
			// x1 = x1 - x2
			k = msub(p, x2);
			t = madd(x1, k);
			mfree(x1);
			x1 = t;
			mmod(x1, p);
			mfree(k);
		} else {
			t = msub(v, u);
			mfree(v);
			v = t;
			// x2 = x2 - x1
			k = msub(p, x1);
			t = madd(x2, k);
			mfree(x2);
			x2 = t;
			mmod(x2, p);
			mfree(k);
		}
	}
	if (MEQUAL(u, 1)) {
		r = x1;
		mfree(x2);
	} else {
		r = x2;
		mfree(x1);
	}
	mfree(u);
	mfree(v);
	return r;
}
