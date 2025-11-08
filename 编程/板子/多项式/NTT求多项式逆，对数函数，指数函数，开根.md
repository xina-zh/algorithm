已知 F(x)，要求 G(x) 令 G(x)≡lnF(x)。
令函数 f(x)=ln(x)，则原式可以化作：
G(x)≡f(F(x)) (mod  xn)
两边求导，发现 f(F(x)) 是个复合函数，复合函数求导公式为 f(g(x))′=f′(g(x))g′(x)，所以左右求导之后为：
G′(x)=f′(F(x))F′(x) (mod  xn)
此时我们可以想一想 ln 的求导公式，ln′(x)=x1​，所以接着可以化为：
G′(x)=F(x)F′(x)​ (mod  xn)

这个时候你可能要问了：诶多项式除法？你在逗我吗？
当然不是辣……我们刚刚才学多项式的逆元，反正我们不要求商只要求余数，为什么不拿出来用呢？
所以我们只需要将读入的 F 求导作为 a，求逆作为 b，计算出 a×b (mod  998244353)，此时求出的是 G′，对它求积分就可以得出我们要求的 G 了。
```cpp
#include<bits/stdc++.h>
using namespace std;
const int P=998244353;
long long qpow(long long a,long long k){
	long long mulv=1;
	while(k)
	{
		if(k&1)mulv=mulv*a%P;
		a=a*a%P;
		k>>=1;
	}
	return mulv;
}
long long inv(long long x){
	return qpow(x%P,P-2);
}
struct Poly:vector<long long>
{
	Poly(int _n):vector<long long>(_n){}
	void change(){
		vector<int> rev(size());
		for(int i=1;i<size();i++)
			rev[i]=rev[i/2]/2+(i&1)*size()/2;
		for(int i=1;i<size();i++)
			if(rev[i]<i)std::swap((*this)[i],(*this)[rev[i]]);
	}
	//ntt 前数组长度需要变为 2^k
	void NTT(int op){
		change();
		for(int l=2;l<=size();l<<=1){
			int gn=qpow(3,(P-1)/l);
			if(op==-1)gn=inv(gn);
			for(int s=0;s<size();s+=l){
				long long g=1;
				for(int i=0;i<l/2;i++,g=g*gn%P){
					long long L=(*this)[s+i],R=(*this)[s+i+l/2];
					(*this)[s+i]=(L+g*R)%P;
					(*this)[s+i+l/2]=(L+P-g*R%P)%P;
				}
			}
		}
		if(op==-1){
			int tmp=inv(size());
			for(auto &ai:*this)ai=ai*tmp%P;
		}
	}
	void display(string name){
		cout<<name<<endl;
		for(auto ai:*this)
			cout<<ai<<" ";
		cout<<endl;
	}
};
Poly polyInv(Poly &f){
	Poly f_inv(1);
	f_inv[0]=inv(f[0]);
	for(int i=2;i<=f.size();i<<=1){
		f_inv.resize(i);
		int L=i<<1;
		Poly g(L),inv_t(L);
		for(int j=0;j<i;j++){
			g[j]=f[j];
			inv_t[j]=f_inv[j];
		}
		g.NTT(1),inv_t.NTT(1);
		for(int j=0;j<L;j++)
			g[j]=g[j]*inv_t[j]%P*inv_t[j]%P;
		g.NTT(-1);
		for(int j=0;j<i;j++)
			f_inv[j]=(f_inv[j]+f_inv[j]+P-g[j])%P;
	}
	return f_inv;
}
Poly polyDeriv(Poly &f){
	Poly ret(f.size());
	for(int i=0;i<f.size()-1;i++)
		ret[i]=f[i+1]*(i+1)%P;
	return ret;
}
Poly polyInt(Poly &f){
	Poly ret(f.size());
	for(int i=1;i<f.size();i++)
		ret[i]=f[i-1]*inv(i)%P;
	return ret;
}
Poly polyLn(Poly &f){
	auto df=polyDeriv(f);
	auto f_inv=polyInv(f);
	df.resize(f.size()*2);
	f_inv.resize(f.size()*2);
	df.NTT(1);
	f_inv.NTT(1);
	for(int i=0;i<df.size();i++)
		df[i]=df[i]*f_inv[i]%P;
	df.NTT(-1);
	auto ret_int=polyInt(df);
	ret_int.resize(f.size());
	return ret_int;
}
Poly polyExp(Poly &f){
	Poly exp_f(1);
	exp_f[0]=1;
	for(int i=2;i<=f.size()*2;i<<=1){
		auto temp=polyLn(exp_f);
		for(int j=0;j<i/2;j++)
			temp[j]=(f[j]-temp[j]+P)%P;
		temp[0]+=1;
		exp_f.resize(i),temp.resize(i);
		temp.NTT(1),exp_f.NTT(1);
		for(int j=0;j<i;j++)
			exp_f[j]=exp_f[j]*temp[j]%P;
		exp_f.NTT(-1);
	}
	exp_f.resize(f.size());
	return exp_f;
}
Poly polySqrt(Poly &f) {
    Poly g(1);
    g[0] = 1;
    long long inv2 = inv(2);
    for (int len = 2; len <= f.size(); len <<= 1) {
        Poly g0 = g;
        g0.resize(len);
        Poly h = polyInv(g0);
        h.resize(len);
        Poly f_len(len);
        for (int i = 0; i < min(len, (int)f.size()); i++) {
            f_len[i] = f[i];
        }
        f_len.resize(2 * len);
        h.resize(2 * len);
        f_len.NTT(1);
        h.NTT(1);
        for (int i = 0; i < 2 * len; i++) {
            f_len[i] = f_len[i] * h[i] % P;
        }
        f_len.NTT(-1);
        f_len.resize(len);
        for (int i = 0; i < len; i++) {
            g0[i] = (g0[i] + f_len[i]) * inv2 % P;
        }
        g = g0;
    }
    return g;
}

int main()
{
    ios::sync_with_stdio(false);
    int n;
    cin>>n;
    int len=1<<(__lg(n)+1);
    Poly a(len);
    for(int i=0;i<n;i++)cin>>a[i];
    auto exp_a=polySqrt(a);
    for(int i=0;i<n;i++)
        cout<<exp_a[i]<<" ";
    return 0;
}

```