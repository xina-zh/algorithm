
```c++
#include<iostream>
#include<queue>
using namespace std;

#define int long long
int n,m,S,T;
const int N=205,M=5005,inf=1e9+7;
int q[N],e[M*2],ne[M*2],w[M*2],idx;
int depth[N];
int now[N];

void add(int a,int b,int c){
    e[idx]=b;w[idx]=c;ne[idx]=q[a];q[a]=idx++;
    e[idx]=a;w[idx]=0;ne[idx]=q[b];q[b]=idx++;
}

bool bfs(){
    for(int i=0;i<=n;i++) depth[i]=inf;
    now[S]=q[S];
    queue<int> que;
    que.push(S);
    depth[S]=1;
    while(que.size()){
        int to=que.front();que.pop();
        for(int i=q[to];i!=-1;i=ne[i]){
            if(!w[i]) continue;
            int j=e[i];
            if(depth[j]!=inf) continue;
            now[j]=q[j];
            que.push(j);
            depth[j]=depth[to]+1;
            if(j==T) return true;
        }
    }
    return false;
}

int dfs(int u,int sum){
    if(u==T) return sum;
    int res=0;
    for(int i=now[u];i!=-1&&res<sum;i=ne[i]){
        now[u]=i;
        if(!w[i]) continue;
        int j=e[i];
        if(depth[j]!=depth[u]+1) continue;
        int k=dfs(j,min(w[i],sum-res));
        if(k==0) depth[j]=inf;
        res+=k;
        w[i]-=k;w[i^1]+=k;
    }
    return res;
}

signed main(){
    ios::sync_with_stdio(false);
    cin.tie(0);cout.tie(0);
    cin>>n>>m>>S>>T;
    for(int i=0;i<=n;i++) q[i]=-1;
    for(int i=1;i<=m;i++){
        int a,b,c;
        cin>>a>>b>>c;
        add(a,b,c);
    }
    int ans=0;
    while(bfs()){
        ans+=dfs(S,inf);
    }
    cout<<ans<<'\n';
    return 0;
}


```