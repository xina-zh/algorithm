```cpp
#include<iostream>
#include<vector>
#include<algorithm>
#define int long long
using namespace std;
const int MOD = 998244353;
const int G = 3; // 原根
int pow_mod(int a, int b) {
    int res = 1;
    while (b) {
        if (b & 1) res = res * a % MOD;
        a = a * a % MOD;
        b >>= 1;
    }
    return res;
}
const int N = 3e5 + 10;
int limit, L = 1;
int a[N], b[N];
int r[N];
void ntt(int *a, int type) {
    for (int i = 0; i < limit; i++) {
        if (r[i] < i) {
            swap(a[i], a[r[i]]);
        }
    }
    for (int mid = 1; mid < limit; mid <<= 1) {
        int wn = pow_mod(G, (MOD - 1) / (mid << 1));
        if (type == -1) wn = pow_mod(wn, MOD - 2);
        for (int len = mid << 1, pos = 0; pos < limit; pos += len) {
            int w = 1;
            for (int k = 0; k < mid; k++, w = w * wn % MOD) {
                int x = a[pos + k];
                int y = w * a[pos + k + mid] % MOD;
                a[pos + k] = (x + y) % MOD;
                a[pos + k + mid] = (x - y + MOD) % MOD;
            }
        }
    }
    if (type == -1) {
        int inv_limit = pow_mod(limit, MOD - 2);
        for (int i = 0; i < limit; i++) {
            a[i] = a[i] * inv_limit % MOD;
        }
    }
}
signed main() {
    ios::sync_with_stdio(false);
    cin.tie(0); cout.tie(0);
    int n, m; cin >> n >> m;
    for (int i = 0; i <= n; i++) {
        cin >> a[i];
        a[i] %= MOD;
    }
    for (int i = 0; i <= m; i++) {
        cin >> b[i];
        b[i] %= MOD;
    }
    while ((1 << L) < n + m + 1) {
        L++;
    }
    limit = (1 << L);
    for (int i = 0; i < limit; i++) {
        r[i] = (r[i >> 1] >> 1) | ((i & 1) << (L - 1));
    }
    ntt(a, 1);
    ntt(b, 1);
    for (int i = 0; i < limit; i++) {
        a[i] = a[i] * b[i] % MOD;
    }
    ntt(a, -1);
    for (int i = 0; i <= n + m; i++) {
        cout << a[i] << ' ';
    }
    return 0;
}
```