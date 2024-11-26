// Compute diagonal of Hamiltonian matrix H from
// "Accurate energies of the He atom with undergraduate quantum mechanics"
// by Robert C. Masse and Thad G. Walker

#include <stdio.h>
#include <math.h>
#include <string.h>

double E0(int n1, int n2);
double integral(int n1, int n2, int n3, int n4);
double P(double r, int n, int l);
double laguerre(double x, int n, int a);
double factorial(int n);

char str[100];

int
main(int argc, char *argv[])
{
	int err;
	double h;

	h = integral(1, 1, 1, 1) + E0(1, 1);
	sprintf(str, "%0.3f", h);
	err = strcmp(str, "-2.750");
	printf("1s1s %g %s\n", h, err ? "err" : "ok");

	h = integral(1, 2, 1, 2) + E0(1, 2);
	sprintf(str, "%0.3f", h);
	err = strcmp(str, "-2.080");
	printf("1s2s %g %s\n", h, err ? "err" : "ok");

	h = integral(1, 3, 1, 3) + E0(1, 3);
	sprintf(str, "%0.3f", h);
	err = strcmp(str, "-2.023");
	printf("1s3s %g %s\n", h, err ? "err" : "ok");

	h = integral(1, 4, 1, 4) + E0(1, 4);
	sprintf(str, "%0.3f", h);
	err = strcmp(str, "-2.010");
	printf("1s4s %g %s\n", h, err ? "err" : "ok");
}

// E0(n1, n2) = -2 / n1^2 - 2 / n2^2

double
E0(int n1, int n2)
{
	return -2.0 / (double) n1 / (double) n1 - 2.0 / (double) n2 / (double) n2;
}

#define R 25.0
#define N 1000

double p13[N], p24[N];

double
integral(int n1, int n2, int n3, int n4)
{
	int i, j;
	double dr, r, sum;

	dr = R / N;

	for (i = 0; i < N; i++) {
		r = (i + 0.5) * dr;
		p13[i] = P(r, n1, 0) * P(r, n3, 0);
		p24[i] = P(r, n2, 0) * P(r, n4, 0);
	}

	sum = 0.0;

	for (i = 0; i < N; i++)
		for (j = 0; j < N; j++) {
			r = ((i > j ? i : j) + 0.5) * dr;
			sum += p13[i] * p24[j] / r;
		}

	return sum * dr * dr;
}

// P(r,n,l) = sqrt(2 (n - l - 1)! / n^2 / (n + l)!) *
//            (4 r / n)^(l + 1) *
//            exp(-2 r / n) *
//            L(4 r / n, n - l - 1, 2 l + 1)

double
P(double r, int n, int l)
{
	double a, b, c, d;
	a = sqrt(2.0 * factorial(n - l - 1) / n / n / factorial(n + l));
	b = pow(4.0 * r / (double) n, (double) (l + 1));
	c = exp(-2.0 * r / (double) n);
	d = laguerre(4.0 * r / (double) n, n - l - 1, 2 * l + 1);
	return a * b * c * d;
}

// L(x,n,a) = (n + a)! sum(k,0,n,(-x)^k / ((n - k)! (a + k)! k!))

double
laguerre(double x, int n, int a)
{
	int k;
	double sum = 0.0;
	for (k = 0; k <= n; k++)
		sum += pow(-x, (double) k) / (factorial(n - k) * factorial(a + k) * factorial(k));
	return factorial(n + a) * sum;
}

double
factorial(int n)
{
	int i;
	double f = 1.0;
	for (i = 1; i <= n; i++)
		f *= (double) i;
	return f;
}
