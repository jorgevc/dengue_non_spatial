/*
Copyright 2012 Jorge Velazquez
*/
#include <stdlib.h>
#include <fcntl.h>
#include <math.h>

static unsigned int xj=123456789,yj=987654321,zj=43219876,cj=6543217;
static int llamadas=0;

#define PHI 0x9e3779b9
#define PI 3.141592654
static unsigned int Q[4096];
static unsigned int cw = 362436;

#pragma omp threadprivate(xj,yj,zj,cj)

unsigned int JKISS(void)
{
	unsigned long long t;

	xj = 314527869 * xj + 1234567;
	yj ^= (yj << 5); yj ^= (yj >> 7); yj ^= (yj << 22);
	t = 4294584393ULL * zj + cj; cj = (t >> 32); zj = t;

	return xj + yj + zj;
}

unsigned int devrand(void)
{
	int fn;
	unsigned int r;

	fn = open("/dev/urandom", O_RDONLY);
	if (fn == -1){exit(-1);}
	if (read(fn, &r, 4) != 4){ exit(-1);}

	close(fn);

	return r;
}

void init_JKISS(void)
{
xj = devrand();
while (!(yj = devrand()));
zj = devrand();
cj = devrand() % 698769068 + 1;
}

void Seed_JKISS(unsigned int x, unsigned int y, unsigned int z, unsigned int c)
{
xj = x;
if(y!=0){yj = y;}else{yj=987654321;}
zj = z;
cj = c;
}

double FS_JKISS(void)
{
return (JKISS() / 4294967296.0);
}

double F_JKISS(void)
{
double d;
unsigned int a, b;

a = JKISS() >> 6;
b = JKISS() >> 5;
d = (a * 134217728.0 + b) / 9007199254740992.0;

return d;
}

int I_JKISS(int l, int u)
{
int i;
int delta;

delta = (u-l) + 1;
	i = (int)(delta * FS_JKISS());
	i = i + l;
return i;
}

void init_CMWC(unsigned int x)
{
int i;

Q[0] = x;
Q[1] = x + PHI;
Q[2] = x + PHI + PHI;

	for (i = 3; i < 4096; i++)
		Q[i] = Q[i-3] ^ Q[i-2] ^ PHI ^ i;
}

unsigned int CMWC(void)
{
unsigned long long t, a = 18782ULL;
static unsigned int i = 4095;
unsigned int x, r = 0xfffffffe;
i = (i + 1) & 4095;
t = a * Q[i] + cw;
cw = (t >> 32);
x = t + cw;
	if (x < cw) {
		x++;
		cw++;
	}
return (Q[i] = r - x);
}

double gaussM()
{
	static double V1, V2, S;
	static int phase = 0;
	double X;

	if(phase == 0) {
		do {
			double U1 = F_JKISS();
			double U2 = F_JKISS();

			V1 = 2 * U1 - 1;
			V2 = 2 * U2 - 1;
			S = V1 * V1 + V2 * V2;
			} while(S >= 1 || S == 0);

		X = V1 * sqrt(-2 * log(S) / S);
	} else
		X = V2 * sqrt(-2 * log(S) / S);

	phase = 1 - phase;

	return X;
}



double gaussBM()
{
	static double U, V;
	static int phase = 0;
	double Z;

	if(phase == 0) {
		U = F_JKISS() + 1./(9007199254740992.0);
		V = F_JKISS();
		Z = sqrt(-2 * log(U)) * sin(2 * PI * V);
	} else
		Z = sqrt(-2 * log(U)) * cos(2 * PI * V);

	phase = 1 - phase;

	return Z;
}


