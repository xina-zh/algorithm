```cpp
#include<iostream>
#include<cmath>
#define int long long
using namespace std;
const double pi=acos(-1);
struct complex{
    double x,y;
    complex(double x=0,double y=0) :x(x),y(y){}

};
complex operator+ (complex a,complex b){
    return {a.x+b.x,a.y+b.y};
}
complex operator- (complex a,complex b){
    return {a.x-b.x,a.y-b.y};
}
complex operator* (complex a,complex b){
    return {a.x*b.x-a.y*b.y,a.x*b.y+a.y*b.x};
}
const int N=3e5+10;
int limit,L=1;
complex a[N],b[N];
int r[N];
void fft(complex *a,int tmp){
    for(int i=0;i<limit;i++){
        if(r[i]>i){
            swap(a[i],a[r[i]]);
        }
    }
    for(int mid=1;mid<limit;mid<<=1){
        complex wn(cos(pi/mid),tmp*sin(pi/mid));
        for(int len=mid<<1,pos=0;pos<limit;pos+=len){
            complex w(1,0);
            for(int k=0;k<mid;k++,w=w*wn){
                complex x=a[k+pos];
                complex y=w*a[k+pos+mid];
                a[pos+k]=x+y;
                a[pos+mid+k]=x-y;
            }
        }
    }
}
signed main(){
    ios::sync_with_stdio(false);
    cin.tie(0);cout.tie(0);
    int n,m;cin>>n>>m;
    for(int i=0;i<=n;i++){
        cin>>a[i].x;
    }
    for(int i=0;i<=m;i++){
        cin>>b[i].x;
    }
    while((n+m)>=(1<<L)){
        L++;
    }
    limit=(1<<L);
    for(int i=0;i<limit;i++){
        r[i]=((r[i>>1]>>1)|((i&1)<<(L-1)));
    }
    fft(a,1);
    fft(b,1);
    for(int i=0;i<limit;i++){
        a[i]=a[i]*b[i];
    }
    fft(a,-1);
    for(int i=0;i<=n+m;i++){
        cout<<(int)(a[i].x/limit+0.5)<<' ';
    }
    return 0;
}
```