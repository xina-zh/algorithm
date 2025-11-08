## D 分牛

通分 $\frac{1}{a}$, $\frac{1}{b}$,  $\frac{1}{c}$ 后观察, 答案就是 $lcm(a,b,c)-n$.

```c++
#define int long long

int lcm(int a, int b) {
    return a * b / gcd(a, b);
}
int solve(int _) {
    int n, a, b, c;
    cin >> n >> a >> b >> c;

    int Lcm = lcm(lcm(a, b), c);

    return Lcm - n;
}
signed main()
{
    ios::sync_with_stdio(false); cin.tie(nullptr); cout.tie(nullptr);
    int _ = 1;
    cin >> _;
    while (_--) {
        cout << solve(_) << endl;
    }
}
```

## E 加训

考虑 DP.

### DP_1

`dp[i][j]`表示 lanhuo 恰好在第 i 秒解出第 j 题的期望.

故可直接枚举当前解决第 j 题是第几次尝试, 注意要额外处理一下:

1.   尝试满 n 次后放弃.
2.   还没尝试满 n 次, 未放弃, 但时间 t 结束了, 在最后计算贡献时处理.

代码:

```c++
//------取模机------//
using i64 = long long;
template<class T>
constexpr T power(T a, i64 b) {
    T res {1};
    for (; b; b /= 2, a *= a) {
        if (b % 2) {
            res *= a;
        }
    }
    return res;
} // 快速幂

constexpr i64 mul(i64 a, i64 b, i64 p) {
    i64 res = a * b - i64(1.L * a * b / p) * p;
    res %= p;
    if (res < 0) {
        res += p;
    }
    return res;
} // 取模乘

template<i64 P>
struct MInt {
    i64 x;
    constexpr MInt() : x {0} {}
    constexpr MInt(i64 x) : x {norm(x % getMod())} {}

    static i64 Mod;
    constexpr static i64 getMod() {
        if (P > 0) {
            return P;
        } else {
            return Mod;
        }
    }
    constexpr static void setMod(i64 Mod_) {
        Mod = Mod_;
    }//只有P<=0, setMod才生效
    constexpr i64 norm(i64 x) const {
        if (x < 0) {
            x += getMod();
        }
        if (x >= getMod()) {
            x -= getMod();
        }
        return x;
    }
    constexpr i64 val() const {
        return x;
    }
    constexpr MInt operator-() const {
        MInt res;
        res.x = norm(getMod() - x);
        return res;
    }
    constexpr MInt inv() const {
        return power(*this, getMod() - 2);
    }
    constexpr MInt &operator*=(MInt rhs) & {
        if (getMod() < (1ULL << 31)) {
            x = x * rhs.x % int(getMod());
        } else {
            x = mul(x, rhs.x, getMod());
        }
        return *this;
    }
    constexpr MInt &operator+=(MInt rhs) & {
        x = norm(x + rhs.x);
        return *this;
    }
    constexpr MInt &operator-=(MInt rhs) & {
        x = norm(x - rhs.x);
        return *this;
    }
    constexpr MInt &operator/=(MInt rhs) & {
        return *this *= rhs.inv();
    }
    friend constexpr MInt operator*(MInt lhs, MInt rhs) {
        MInt res = lhs;
        res *= rhs;
        return res;
    }
    friend constexpr MInt operator+(MInt lhs, MInt rhs) {
        MInt res = lhs;
        res += rhs;
        return res;
    }
    friend constexpr MInt operator-(MInt lhs, MInt rhs) {
        MInt res = lhs;
        res -= rhs;
        return res;
    }
    friend constexpr MInt operator/(MInt lhs, MInt rhs) {
        MInt res = lhs;
        res /= rhs;
        return res;
    }
    friend constexpr std::istream &operator>>(std::istream &is, MInt &a) {
        i64 v;
        is >> v;
        a = MInt(v);
        return is;
    }
    friend constexpr std::ostream &operator<<(std::ostream &os, const MInt &a) {
        return os << a.val();
    }
    friend constexpr bool operator==(MInt lhs, MInt rhs) {
        return lhs.val() == rhs.val();
    }
    friend constexpr bool operator!=(MInt lhs, MInt rhs) {
        return lhs.val() != rhs.val();
    }
    friend constexpr bool operator<(MInt lhs, MInt rhs) {
        return lhs.val() < rhs.val();
    }
};

template<>
i64 MInt<0>::Mod = 998244353; //只有P<=0, Mod才生效

constexpr int P = 1e9 + 7; //在这设置要用的模数
using Z = MInt<P>;
//------取模机------//

int n, t;
Z a, b, dp[2005][2005], fail[25], win[25];
signed main()
{
    ios::sync_with_stdio(false); cin.tie(nullptr); cout.tie(nullptr);
    cin >> n >> t;

    fail[0] = 1;
    for (int i = 1; i <= n; ++i) {
        cin >> a >> b;
        win[i] = a / b;
        fail[i] = fail[i - 1] * (1 - win[i]); // 连败 i 场的概率
    }

    // dp[i][j] 表示恰好在第 i 秒 AC 了第 j 道题
    dp[0][0] = 1;
    for (int i = 1; i <= t; ++i) { // 在第 i 秒,
        for (int j = 0; j <= i; ++j) { // AC 第 j 道题
            for (int y = 1; j && y <= min(i, n); ++y) { // 尝试第 y 次 AC
                dp[i][j] += dp[i - y][j - 1] * fail[y - 1] * win[y];
            }
            if (i >= n)
                dp[i][j] += dp[i - n][j] * fail[n]; // 在第 j 秒尝试了 n 次然后放弃
        }
    }

    Z ans = 0;
    for (int i = 0; i < n; ++i) { // 前面有 i 次无用的尝试, 但未满 n 次
        for (int j = 1; j <= t; ++j) { // 过题数
            ans += dp[t - i][j] * fail[i] * j;
        }
    }

    cout << ans << endl;
}
```

### DP_2

`dp[i][j]`表示当前是第 i 次尝试, 已经解出了 j 题.

```c++
ll n, t, a[25], b[25], c[25], d[25];
ll dp[25][2005];
ll ans, sum;
void solve() {
    cin >> n >> t;
    for (ll i = 1; i <= n; ++i)cin >> a[i] >> b[i];
    for (ll i = 1; i <= n; ++i) {
        c[i] = a[i] * ksm(b[i], mod - 2) % mod;
        d[i] = (1 + mod - c[i]) % mod;
    }
    //c[i]成功概率 d[i]失败概率

    for (ll i = 0; i < 25; ++i) {
        for (ll j = 0; j < 2005; ++j) {
            dp[i][j] = 0;
        }
    }
    dp[0][0] = 1;
    //初始化

    for (ll i = 0; i < t; ++i) {
        for (ll j = i + 1; j >= 0; --j) {
            dp[0][j + 1] = (dp[0][j + 1] + dp[n - 1][j] * c[n] % mod) % mod; //第n次成功
            dp[n][j] = (dp[n][j] + dp[n - 1][j] * d[n] % mod) % mod; //失败，过题数不变，失败次数清零
            dp[n - 1][j] = 0; //第n次，没过直接放弃

            for (ll k = n - 2; k >= 0; --k) { //计算前n-1次
                dp[0][j + 1] = (dp[0][j + 1] + dp[k][j] * c[k + 1] % mod) % mod; //成功，题数+1，且失败次数更新为0
                dp[k + 1][j] = (dp[k + 1][j] + dp[k][j] * d[k + 1] % mod) % mod; //失败，题数不变，失败次数+1
                dp[k][j] = 0;
            }
            dp[0][j] = (dp[0][j] + dp[n][j]) % mod;
            dp[n][j] = 0;
        }
    }

    ans = 0;
    for (ll i = 1; i <= t; ++i) {
        sum = 0;
        for (ll j = 0; j < n; ++j) {
            sum = (sum + dp[j][i]) % mod;
        }
        ans = (ans + sum * i % mod) % mod;
    }
    cout << ans << endl;
}
```

## G 我要遍历这棵树 (1)

1. 对于固定根的树: 对于每个节点, **以它为根节点的子树**的答案是每个**以它的儿子节点为根节点的子树**的答案的乘积, 再乘上 $A_{num}^{num}$ (num是它的儿子数量). 最后 dfs 处理出来根节点的答案就好.

2. 不固定根时: 考虑换根 DP. 

    对于以 x 为根的树, 转移方程为: 

    ```c++
    ans[x] = (ans[fa] / (fa所有连边的全排列) / ans[x] * (fa所有连边减掉x的全排列) ) * (x 所有儿子子树的答案乘积) * (x所有连边的全排列);
    ```

实际上就是在 x 的父亲的答案 `ans[fa]`中去掉子树 `ans[x]`贡献. 再把 fa 去掉 son 后的部分作为 son 的子树重新计算以 x 为根的整颗树的答案.

```c++
#include<bits/stdc++.h>
using namespace std;
#ifdef LOCAL
#include "D:\Users\TauLee\OneDrive\Mine\c++\debug.h"
#else
#define debug(...) 42;
#define endl '\n'
#endif
// #pragma comment(linker, "/STACK:102400000,102400000")

//------取模机------//
using i64 = long long;
template<class T>
constexpr T power(T a, i64 b) {
    T res {1};
    for (; b; b /= 2, a *= a) {
        if (b % 2) {
            res *= a;
        }
    }
    return res;
} // 快速幂

constexpr i64 mul(i64 a, i64 b, i64 p) {
    i64 res = a * b - i64(1.L * a * b / p) * p;
    res %= p;
    if (res < 0) {
        res += p;
    }
    return res;
} // 取模乘

template<i64 P>
struct MInt {
    i64 x;
    constexpr MInt() : x {0} {}
    constexpr MInt(i64 x) : x {norm(x % getMod())} {}

    static i64 Mod;
    constexpr static i64 getMod() {
        if (P > 0) {
            return P;
        } else {
            return Mod;
        }
    }
    constexpr static void setMod(i64 Mod_) {
        Mod = Mod_;
    }//只有P<=0, setMod才生效
    constexpr i64 norm(i64 x) const {
        if (x < 0) {
            x += getMod();
        }
        if (x >= getMod()) {
            x -= getMod();
        }
        return x;
    }
    constexpr i64 val() const {
        return x;
    }
    constexpr MInt operator-() const {
        MInt res;
        res.x = norm(getMod() - x);
        return res;
    }
    constexpr MInt inv() const {
        return power(*this, getMod() - 2);
    }
    constexpr MInt &operator*=(MInt rhs) & {
        if (getMod() < (1ULL << 31)) {
            x = x * rhs.x % int(getMod());
        } else {
            x = mul(x, rhs.x, getMod());
        }
        return *this;
    }
    constexpr MInt &operator+=(MInt rhs) & {
        x = norm(x + rhs.x);
        return *this;
    }
    constexpr MInt &operator-=(MInt rhs) & {
        x = norm(x - rhs.x);
        return *this;
    }
    constexpr MInt &operator/=(MInt rhs) & {
        return *this *= rhs.inv();
    }
    friend constexpr MInt operator*(MInt lhs, MInt rhs) {
        MInt res = lhs;
        res *= rhs;
        return res;
    }
    friend constexpr MInt operator+(MInt lhs, MInt rhs) {
        MInt res = lhs;
        res += rhs;
        return res;
    }
    friend constexpr MInt operator-(MInt lhs, MInt rhs) {
        MInt res = lhs;
        res -= rhs;
        return res;
    }
    friend constexpr MInt operator/(MInt lhs, MInt rhs) {
        MInt res = lhs;
        res /= rhs;
        return res;
    }
    friend constexpr std::istream &operator>>(std::istream &is, MInt &a) {
        i64 v;
        is >> v;
        a = MInt(v);
        return is;
    }
    friend constexpr std::ostream &operator<<(std::ostream &os, const MInt &a) {
        return os << a.val();
    }
    friend constexpr bool operator==(MInt lhs, MInt rhs) {
        return lhs.val() == rhs.val();
    }
    friend constexpr bool operator!=(MInt lhs, MInt rhs) {
        return lhs.val() != rhs.val();
    }
    friend constexpr bool operator<(MInt lhs, MInt rhs) {
        return lhs.val() < rhs.val();
    }
};

template<>
i64 MInt<0>::Mod = 998244353; //只有P<=0, Mod才生效

constexpr int P = 1e9 + 7; //在这设置要用的模数
using Z = MInt<P>;
//------取模机------//

//----计算组合数----//
struct Comb {
    int n;
    std::vector<Z> _fac; //阶乘
    std::vector<Z> _invfac; //阶乘的逆元
    std::vector<Z> _inv; //数字的逆元

    Comb() : n{0}, _fac{1}, _invfac{1}, _inv{0} {}
    Comb(int n) : Comb() {
        init(n);
    }

    void init(int m) {
        m = std::min<i64>(m, Z::getMod() - 1);
        if (m <= n) return;
        _fac.resize(m + 1);
        _invfac.resize(m + 1);
        _inv.resize(m + 1);

        for (int i = n + 1; i <= m; i++) {
            _fac[i] = _fac[i - 1] * i;
        }
        _invfac[m] = _fac[m].inv();
        for (int i = m; i > n; i--) {
            _invfac[i - 1] = _invfac[i] * i;
            _inv[i] = _invfac[i] * _fac[i - 1];
        }
        n = m;
    }

    Z fac(int m) {
        if (m > n) init(2 * m);
        return _fac[m];
    }
    Z invfac(int m) {
        if (m > n) init(2 * m);
        return _invfac[m];
    }
    Z inv(int m) {
        if (m > n) init(2 * m);
        return _inv[m];
    }
    Z C(int n, int m) {
        if (n < m || m < 0) return 0;
        return fac(n) * invfac(m) * invfac(n - m);
    }
    Z A(int n, int m) {
        if (n < m || m < 0 ) return 0;
        return fac(n) * invfac(n - m);
    }
} comb; // 调用时会自动扩展2倍
//----计算组合数----//

const int Maxn = 1e5 + 10;
int n, v, u, sonum[Maxn];
vector<int>eg[Maxn];
Z ans[Maxn];

Z dfs(int x, int fa) {
    Z tp = 1;
    int num = 0;
    for (auto c : eg[x]) {
        if (c == fa)continue;
        num++;
        tp *= dfs(c, x);
    }
    sonum[x] = num;
    return ans[x] = tp * comb.fac(num);
}

void dfs1(int x, int fa) {
    Z tp = 1;
    for (auto c : eg[x]) {
        if (c == fa)continue;
        tp *= ans[c];
    }

    sonum[x]++;
    ans[x] = (ans[fa] / comb.fac(sonum[fa]) / ans[x] * comb.fac(sonum[fa] - 1)) * \
             tp * comb.fac(sonum[x]);

    for (auto c : eg[x]) {
        if (c == fa)continue;
        dfs1(c, x);
    }
}
signed main()
{
    ios::sync_with_stdio(false); cin.tie(nullptr); cout.tie(nullptr);
    freopen("6.in", "r", stdin);
    freopen("6.out", "w", stdout);
    cin >> n;

    for (int i = 1; i < n; ++i) {
        cin >> v >> u;
        eg[u].emplace_back(v);
        eg[v].emplace_back(u);
    }

    ans[1] = dfs(1, 1);
    for (auto c : eg[1]) {
        dfs1(c, 1);
    }

    Z sum = 0;
    for (int i = 1; i <= n; ++i) {
        sum += ans[i];
    }
    cout << sum << endl;
}
```

## I 来打瓦

从 $(1,1)$ 到 $(i,j)$ 所有路径数量为 $C_{i+j-2}^{i-1}$, 即总共 $i+j-2$ 步, 选择其中 $i-1$ 步向右. 

用 `f[i][j]`表示从 $(1,1)$ 到 $(i,j)$ 的合法路径数量. 枚举每个地雷作为 $(i,j)$, 更新 f 数组.

更新 `f[i][j]`时, 遍历所有在 **以(1,1)为左下角, (i, j)为右上角的矩形** 内的地雷. 对于每个地雷 (a,b), 对 `f[i][j]`的贡献为负的 `f[a][b]*[从(a,b)到(i,j)的路径数量]`.

代码:

```c++
//------取模机------//
using i64 = long long;
template<class T>
constexpr T power(T a, i64 b) {
    T res {1};
    for (; b; b /= 2, a *= a) {
        if (b % 2) {
            res *= a;
        }
    }
    return res;
} // 快速幂

constexpr i64 mul(i64 a, i64 b, i64 p) {
    i64 res = a * b - i64(1.L * a * b / p) * p;
    res %= p;
    if (res < 0) {
        res += p;
    }
    return res;
} // 取模乘

template<i64 P>
struct MInt {
    i64 x;
    constexpr MInt() : x {0} {}
    constexpr MInt(i64 x) : x {norm(x % getMod())} {}

    static i64 Mod;
    constexpr static i64 getMod() {
        if (P > 0) {
            return P;
        } else {
            return Mod;
        }
    }
    constexpr static void setMod(i64 Mod_) {
        Mod = Mod_;
    }//只有P<=0, setMod才生效
    constexpr i64 norm(i64 x) const {
        if (x < 0) {
            x += getMod();
        }
        if (x >= getMod()) {
            x -= getMod();
        }
        return x;
    }
    constexpr i64 val() const {
        return x;
    }
    constexpr MInt operator-() const {
        MInt res;
        res.x = norm(getMod() - x);
        return res;
    }
    constexpr MInt inv() const {
        return power(*this, getMod() - 2);
    }
    constexpr MInt &operator*=(MInt rhs) & {
        if (getMod() < (1ULL << 31)) {
            x = x * rhs.x % int(getMod());
        } else {
            x = mul(x, rhs.x, getMod());
        }
        return *this;
    }
    constexpr MInt &operator+=(MInt rhs) & {
        x = norm(x + rhs.x);
        return *this;
    }
    constexpr MInt &operator-=(MInt rhs) & {
        x = norm(x - rhs.x);
        return *this;
    }
    constexpr MInt &operator/=(MInt rhs) & {
        return *this *= rhs.inv();
    }
    friend constexpr MInt operator*(MInt lhs, MInt rhs) {
        MInt res = lhs;
        res *= rhs;
        return res;
    }
    friend constexpr MInt operator+(MInt lhs, MInt rhs) {
        MInt res = lhs;
        res += rhs;
        return res;
    }
    friend constexpr MInt operator-(MInt lhs, MInt rhs) {
        MInt res = lhs;
        res -= rhs;
        return res;
    }
    friend constexpr MInt operator/(MInt lhs, MInt rhs) {
        MInt res = lhs;
        res /= rhs;
        return res;
    }
    friend constexpr std::istream &operator>>(std::istream &is, MInt &a) {
        i64 v;
        is >> v;
        a = MInt(v);
        return is;
    }
    friend constexpr std::ostream &operator<<(std::ostream &os, const MInt &a) {
        return os << a.val();
    }
    friend constexpr bool operator==(MInt lhs, MInt rhs) {
        return lhs.val() == rhs.val();
    }
    friend constexpr bool operator!=(MInt lhs, MInt rhs) {
        return lhs.val() != rhs.val();
    }
    friend constexpr bool operator<(MInt lhs, MInt rhs) {
        return lhs.val() < rhs.val();
    }
};

template<>
i64 MInt<0>::Mod = 998244353; //只有P<=0, Mod才生效

constexpr int P = 1e9 + 7; //在这设置要用的模数
using Z = MInt<P>;
//------取模机------//

//----计算组合数----//
struct Comb {
    int n;
    std::vector<Z> _fac; //阶乘
    std::vector<Z> _invfac; //阶乘的逆元
    std::vector<Z> _inv; //数字的逆元

    Comb() : n{0}, _fac{1}, _invfac{1}, _inv{0} {}
    Comb(int n) : Comb() {
        init(n);
    }

    void init(int m) {
        m = std::min<i64>(m, Z::getMod() - 1);
        if (m <= n) return;
        _fac.resize(m + 1);
        _invfac.resize(m + 1);
        _inv.resize(m + 1);

        for (int i = n + 1; i <= m; i++) {
            _fac[i] = _fac[i - 1] * i;
        }
        _invfac[m] = _fac[m].inv();
        for (int i = m; i > n; i--) {
            _invfac[i - 1] = _invfac[i] * i;
            _inv[i] = _invfac[i] * _fac[i - 1];
        }
        n = m;
    }

    Z fac(int m) {
        if (m > n) init(2 * m);
        return _fac[m];
    }
    Z invfac(int m) {
        if (m > n) init(2 * m);
        return _invfac[m];
    }
    Z inv(int m) {
        if (m > n) init(2 * m);
        return _inv[m];
    }
    Z C(int n, int m) {
        if (n < m || m < 0) return 0;
        return fac(n) * invfac(m) * invfac(n - m);
    }
    Z A(int n, int m) {
        if (n < m || m < 0 ) return 0;
        return fac(n) * invfac(m);
    }
} comb; // 调用时会自动扩展2倍
//----计算组合数----//
signed main()
{
    ios::sync_with_stdio(false); cin.tie(nullptr); cout.tie(nullptr);
    int n, m, x;
    cin >> n >> m >> x;
    vector<pair<int, int>>v(x);
    for (auto &[a, b] : v)cin >> a >> b;

    sort(v.begin(), v.end());
    Z AllRoad = comb.C(n + m - 2, n - 1);
    vector<Z>Road(x, 0);
    for (int i = 0; i < x; ++i) {
        int x = v[i].first, y = v[i].second;
        Road[i] = comb.C(x + y - 2, x - 1);

        for (int j = 0; j < i; ++j) {
            int xx = v[j].first, yy = v[j].second;
            if (xx <= x && yy <= y) {
                Road[i] -= Road[j] * comb.C(x - xx + y - yy, x - xx);
            }
        }
    }
    Z sum = 0;
    for (int i = 0; i < x; ++i) {
        Road[i] *= comb.C(n - v[i].first + m - v[i].second, n - v[i].first);
        sum += Road[i];
    }

    Z ans = 1 - (sum / AllRoad);
    cout << ans << endl;
}
```

