FTT实现速度快，精度差
```cpp
#include <cstdio>
#include <complex>
#define debug(...) fprintf(stderr, __VA_ARGS__)
typedef long long lolong;
typedef std::complex<double> complex;
inline int input() { int x; scanf("%d", &x); return x; }
inline lolong linput() { lolong x; scanf("%lld", &x); return x; }
const int maxn = 400005, maxk = 20;
const complex I(0, 1);
int R[maxn];
complex Wn[maxn];
void FFT(complex *A, int n, int t) {
	if(t == -1)
		for(int i = 1; i < n; i ++)
			if(i < (n - i))
				std::swap(A[i], A[n - i]);
	for(int i = 0; i < n; i ++)
		if(i < R[i])
			std::swap(A[i], A[R[i]]);

	for(int m = 1, l = 0; m < n; m <<= 1, l ++) {
		/* complex Wn(cos(M_PI / m), sin(M_PI / m) * t); */
		for(int i = 0; i < n; i += m << 1) {
			/* complex W = 1; */
			for(int k = i; k < i + m; k ++) {
				/* complex W(cos(M_PI / m * (k - i)), sin(M_PI / m * (k - i)) * t); */
				complex W = Wn[1ll * (k - i) * n / m];
				/* if(t == -1) W = std::conj(W); */
				complex a0 = A[k], a1 = A[k + m] * W;
				A[k] = a0 + a1;
				A[k + m] = a0 - a1;
				/* W *= Wn; */
			}
		}
	}

	if(t == -1)
		for(int i = 0; i < n; i ++)
			A[i] /= n;
}

int mod;
inline lolong num(complex x) {
	double d = x.real();
	return d < 0 ? lolong(d - 0.5) % mod : lolong(d + 0.5) % mod;
}

inline void FFTFFT(complex *a, complex *b, int len, int t) {
	for(int i = 0; i < len; i ++)
		a[i] = a[i] + I * b[i];
	FFT(a, len, t);
	for(int i = 0; i < len; i ++)
		b[i] = std::conj(a[i ? len - i : 0]);
	for(int i = 0; i < len; i ++) {
		complex p = a[i], q = b[i];
		a[i] = (p + q) * 0.5;
		b[i] = (q - p) * 0.5 * I;
	}
}

complex a0[maxn], a1[maxn], b0[maxn], b1[maxn];
/* complex a0b0[maxn], a1b0[maxn], a0b1[maxn], a1b1[maxn]; */
complex p[maxn], q[maxn];

int main() {
	int n = input(), m = input();
	mod = input();
	int M = int(sqrt(mod) + 1);

	for(int i = 0; i <= n; i ++) {
		int x = input() % mod;
		a0[i] = x / M;
		a1[i] = x % M;
	}
	for(int i = 0; i <= m; i ++) {
		int x = input() % mod;
		b0[i] = x / M;
		b1[i] = x % M;
	}
	int len = 1;
	while(len < n + m + 1)
		len <<= 1;

	for(int i = 1; i < len; i ++)
		R[i] = R[i >> 1] >> 1 | ((i & 1) * (len >> 1));

	for(int i = 0; i < len; i ++)
		Wn[i] = complex(cos(M_PI / len * i), sin(M_PI / len * i));

	FFTFFT(a0, a1, len, 1);
	FFTFFT(b0, b1, len, 1);

	for(int i = 0; i < len; i ++) {
		p[i] = a0[i] * b0[i] + I * a1[i] * b0[i];
		q[i] = a0[i] * b1[i] + I * a1[i] * b1[i];
	}

	FFT(p, len, -1);
	FFT(q, len, -1);

	for(int i = 0; i <= n + m; i ++)
		printf("%lld ", (M * M * num(p[i].real()) % mod +
				M * (num(p[i].imag()) + num(q[i].real())) % mod +
				num(q[i].imag())) % mod);
	puts("");
}

```
NTT实现，速度慢精度高
```cpp
#include <bits/stdc++.h>
//#include <bits/extc++.h>
#define N 800005 
#define __BEGIN_MULTITEST__ \
	signed T; \
	scanf("%d",&T); \
	while(T--) \
	{
#define __END_MULTITEST__ }
using namespace std;
//using namespace __gnu_cxx;
//using namespace __gnu_pbds;
namespace FNTT
{
	const long long modA=998244353,modB=1004535809,modC=469762049,G=3,M=1ll*modA*modB;
	long long mod;
	int rev[N];
	long long quick_pow(long long a,long long b,const long long mod)
	{
		int ret=1;
		while(b)
		{
			if(b&1)
				ret=1ll*ret*a%mod;
			a=1ll*a*a%mod;
			b>>=1;
		}
		return ret;
	}
	const long long MA=quick_pow(modB%modA,modA-2,modA),MB=quick_pow(modA%modB,modB-2,modB),MC=quick_pow(M%modC,modC-2,modC);
	long long slow_mul(long long a,long long b,const long long mod)
	{
		long long ret=0;
		while(b)
		{
			if(b&1)
				ret=(ret+a)%mod;
			a=(a<<1)%mod;
			b>>=1;
		}
		return ret;
	}
	void NTTinit(int &lim,int n)
	{
		lim=1;
		int l=-1;
		while(lim<=(n<<1))
		{
			lim<<=1;
			l++;
		}
		for(int i=0;i<lim;i++)	
			rev[i]=(rev[i>>1]>>1)|((i&1)<<l);
	}
	void NTT(int f[],int type,int lim,int mod)
	{
		for(int i=0;i<lim;i++)
			if(i<rev[i])
				swap(f[i],f[rev[i]]);
		for(int k=1;k<lim;k<<=1)
			for(int i=0,w=quick_pow(G,(mod-1)/(k<<1),mod);i<lim;i+=(k<<1))
				for(int j=0,prodw=1;j<k;j++,prodw=1ll*prodw*w%mod)
				{
					int x=f[i|j],y=1ll*f[i|j|k]*prodw%mod;
					f[i|j]=(x+y)%mod;
					f[i|j|k]=(x-y+mod)%mod;
				}
		if(type==-1)
		{
			int tmpinv=quick_pow(lim,mod-2,mod);
			f[0]=1ll*f[0]*tmpinv%mod;
			for(int i=1;i<=(lim>>1);i++)
			{
				f[i]=1ll*f[i]*tmpinv%mod;
				if(i!=lim-i)
					f[lim-i]=1ll*f[lim-i]*tmpinv%mod;
				swap(f[i],f[lim-i]);
			}
		}
	}
	long long CRT(int a,int b,int c)
	{
		long long A=(slow_mul(1ll*modB*a%M,MA,M)+slow_mul(1ll*modA*b%M,MB,M))%M;
		long long k=(c-A%modC+modC)%modC*1ll*MC%modC;
		return (1ll*k%mod*1ll*(M%mod)%mod+A%mod)%mod;
	}
	void mulNTT(const int f[],const int g[],int h[],int mod,int n,int lim)
	{
		static int a[N],b[N];
		memcpy(a,f,sizeof(int)*n);
		memcpy(b,g,sizeof(int)*n);
		memset(a+n,0,sizeof(int)*(lim-n));
		memset(b+n,0,sizeof(int)*(lim-n));
		NTT(a,1,lim,mod);
		NTT(b,1,lim,mod);
		for(int i=0;i<lim;i++)
			h[i]=1ll*a[i]*b[i]%mod;
		NTT(h,-1,lim,mod);
		memset(h+n,0,sizeof(int)*(lim-n));
	}
	void MTT(int f[],int g[],int h[],int n,int m)
	{
		int lim;
		NTTinit(lim,n+m);
		static int tmpa[N],tmpb[N],tmpc[N];
		mulNTT(f,g,tmpa,modA,n+m+1,lim);
		mulNTT(f,g,tmpb,modB,n+m+1,lim);
		mulNTT(f,g,tmpc,modC,n+m+1,lim);
		for(int i=0;i<=n+m;i++)
			h[i]=(CRT(tmpa[i],tmpb[i],tmpc[i])%mod+mod)%mod;
	}
}
using namespace FNTT;
int n,m;
int f[N],g[N],h[N];
signed main()
{
//	__BEGIN_MULTITEST__
	scanf("%d%d%lld",&n,&m,&mod);
	for(int i=0;i<=n;i++)
		scanf("%d",&f[i]);
	for(int i=0;i<=m;i++)
		scanf("%d",&g[i]);
	MTT(f,g,h,n,m);
	for(int i=0;i<=n+m;i++)
		printf("%d ",h[i]);
	printf("\n");
//	__END_MULTITEST__
	return 0;
}

```