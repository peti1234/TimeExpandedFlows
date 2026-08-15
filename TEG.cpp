#include <bits/stdc++.h>

using namespace std;
const int c=1000005; // az uj csucsok szama tobb lehet
const int c2=1005;
const int inf=1e9, eps=1e-6;

int n;
double t;


int node_cnt;


int szint[c], pos[c], kezd, veg;
bool v[c];
vector<int> sz[c], ki[c];
vector<double> s[c];
queue<int> q;
void add_edge(int a, int b, double s1, double s2=0) {
    //cout << "add " << a << " " << b << " " << s1 << " " << s2 << "\n";
    ki[a].push_back(sz[b].size()), ki[b].push_back(sz[a].size());
    sz[a].push_back(b), sz[b].push_back(a);
    s[a].push_back(s1), s[b].push_back(s2);
}
bool bfs() {
    for (int i=1; i<=node_cnt; i++) {
        szint[i]=0;
    }
    szint[kezd]=1;
    q.push(kezd);
    while(q.size()>0) {
        int a=q.front();
        q.pop();
        for (int i=0; i<sz[a].size(); i++) {
            int x=sz[a][i];
            double y=s[a][i];
            if (y>eps && !szint[x]) {
                szint[x]=szint[a]+1;
                q.push(x);
            }
        }
    }
    return szint[veg];
}
double dfs(int a, double d) {
    if (a==veg) {
        return d;
    }
    double sum=0;
    for (; pos[a]<sz[a].size(); pos[a]++) {
        int x=sz[a][pos[a]];
        double y=s[a][pos[a]];
        if (szint[a]+1==szint[x] && y>eps) {
            double z=dfs(x, min(y, d-sum));
            if (z>eps) {
                s[a][pos[a]]-=z;
                s[x][ki[a][pos[a]]]+=z;
                sum+=z;
            }
        }
    }
    return sum;
}
double dinic() {
    //cout << "folyam: " << node_cnt << "\n";
    if (node_cnt>5000) {
        cout << "tul sok csucs\n";
        exit(0);
    }
    kezd=1, veg=node_cnt;
    double ans=0;
    while(bfs()) {
        for (int i=1; i<=node_cnt; i++) {
            pos[i]=0;
        }

        double x=dfs(kezd, inf);
        while(x>eps) {
            ans+=x;
            x=dfs(kezd, inf);
        }
    }
    return ans;
}










double connection_len=0.1;
bool local_test=1;
vector<pair<double, double> > meeting_points[c2][c2];

void generate_moving_objects() {
    for (double time=0; time<=t; time+=connection_len) {
        vector<int> pos(n+1, 0);
        for (int i=1; i<=n; i++) {
            pos[i]=rand()%100;
        }
        for (int i=1; i<=n; i++) {
            for (int j=1; j<=n; j++) {
                if (i!=j && abs(pos[i]-pos[j])<=5) {
                    meeting_points[i][j].push_back({time, abs(pos[i]-pos[j])});
                }
            }
        }
    }
}

double get_integral(int x, int y, double l, double r) {
    double ans=0;
    for (auto p:meeting_points[x][y]) {
        double t=p.first, dist=p.second;
        double t_kezd=t-connection_len/2, t_veg=t+connection_len/2;
        double val=1/(dist+1);

        double kozos=min(r, t_veg)-max(l, t_kezd);
        if (kozos>0) {
            ans+=kozos*val;
        }
    }
    return ans;
}




double find_integral(int x, int y, double l, double r) {
    if (local_test) {
        return get_integral(x, y, l, r);
    }
    //cout << "integrate " << x << " " << y << " " << l << " " << r << endl;
    double res;
    cin >> res;
    return res;
}



vector<double> belso[c2];
vector<double> kozos[c2][c2];
map<double, int> m[c];
void rec(int a, int b, double l, double r, double lim) {
    double val=find_integral(a, b, l, r);
    if (val>lim) {
        double m=(l+r)/2;
        belso[a].push_back(m), belso[b].push_back(m);
        kozos[a][b].push_back(m);
        rec(a, b, l, m, lim);
        rec(a, b, m, r, lim);
    }
}
void generate_dynamic_graph(double lim) {
    for (int i=1; i<=n; i++) {
        belso[i].push_back(0), belso[i].push_back(t);
    }
    for (int i=1; i<=n; i++) {
        for (int j=1; j<=n; j++) {
            if (i!=j) {
                rec(i, j, 0, t, lim);
            }
        }
    }
}

void generate_fixed_gap_graph(double fix_gap) {
    for (int i=1; i<=n; i++) {
        belso[i].push_back(0), belso[i].push_back(t);
    }
    for (int i=1; i<=n; i++) {
        for (int j=1; j<=n; j++) {
            if (i!=j) {
                for (double time=fix_gap; time<t; time+=fix_gap) {
                    kozos[i][j].push_back(time);
                }
            }
        }
        for (double time=fix_gap; time<t; time+=fix_gap) {
            belso[i].push_back(time);
        }
    }
}
void construct_graph(bool is_upper_estimate=1) {
    for (int i=1; i<=n; i++) {
        sort(belso[i].begin(), belso[i].end());
    }
    for (int i=1; i<=n; i++) {
        for (auto x:belso[i]) {
            m[i][x]=++node_cnt;
        }
        int si=belso[i].size();
        for (int i=node_cnt-si+1; i<node_cnt; i++) {
            add_edge(i, i+1, inf);
        }
    }

    for (int i=1; i<=n; i++) {
        for (int j=1; j<=n; j++) {
            if (i==j) continue;
            kozos[i][j].push_back(0), kozos[i][j].push_back(t);
            sort(kozos[i][j].begin(), kozos[i][j].end());

            int si=kozos[i][j].size();

            for (int k=1; k<si; k++) {
                double a=kozos[i][j][k-1], b=kozos[i][j][k];
                double val=find_integral(i, j, a, b);
                int kezd=m[i][a];
                int veg=(is_upper_estimate ? m[j][a] : m[j][b]);
                add_edge(kezd, veg, val);
            }
        }
    }
}

void calc_flow(bool is_fixed_gap, double precision, bool is_upper_estimate) {
    if (is_fixed_gap) {
        generate_fixed_gap_graph(precision);
    } else {
        generate_dynamic_graph(precision);
    }
    construct_graph(is_upper_estimate);

    double ans=dinic();
    cout << "size of maxflow: " << ans << "\n";

    kezd=0, veg=0;
    for (int i=0; i<=node_cnt; i++) {
        szint[i]=0, pos[i]=0, v[i]=0;
        sz[i].clear(), ki[i].clear(), s[i].clear();
    }
    node_cnt=0;
    for (int i=0; i<=n; i++) {
        belso[i].clear();
        for (int j=0; j<=n; j++) {
            kozos[i][j].clear();
        }
    }
}

int main()
{
    /*int x;
    cin >> x;
    cout << x << "\n";
    return 0;*/

    //cout << "csucsok szama: ";
    cin >> n;
    //n=10;
    //cout << "\n";

    //cout << "ido: ";
    cin >> t;
    //t=20;
    //cout << "\n";

    //double lim=1;
    //cout << "elvart pontossag: ";
    //cin >> lim;
    //cout << "\n";

    if (local_test) {
        generate_moving_objects();
    }

    calc_flow(1, 1, 0);
    calc_flow(1, 1, 1);
    return 0;
    calc_flow(1, 0.8, 0);
    calc_flow(1, 0.8, 1);
    calc_flow(1, 0.2, 0);
    calc_flow(1, 0.2, 1);
    calc_flow(0, 1, 0);
    calc_flow(0, 1, 1);
    calc_flow(0, 0.5, 0);
    calc_flow(0, 0.5, 1);
    calc_flow(0, 0.2, 0);
    calc_flow(0, 0.2, 1);
    calc_flow(0, 0.1, 0);
    calc_flow(0, 0.1, 1);
    calc_flow(0, 0.05, 0);
    calc_flow(0, 0.05, 1);

    //generate_fixed_gap_graph(lim);
    //generate_dynamic_graph(lim);
    //construct_graph(0);
    //double ans=dinic();
    //cout << ans << "\n";
    return 0;
}
/*
2
1 1
1
1
1
1
*/

/*

hogyan lehet bizonyitani, hogy a valasz legalabb valamennyi?


tobbnek a kisizamolasa egyszerre

TEG-es otlet alapbol, ezen also felso becsles (kesz)

n rogzitett pont, ezen minel jobb beosztas

itt az also becslesnel lefele, felso becslesnel felfele (kesobbi idopontbol korabbiba kell mennie), az nem eleg, hogy egy szinten van
1-2 el konstans 3, 2-3 el konstans 2, ha az 1-2 el tobb reszre van osztva es a 2-3 el a legelejen atmegy, akkor elromlik

folyamot kiszamoljuk megnezzuk melyik eleken van a legnagyobb elteres, ezeket erdemes tobb reszre osztani
ebbol valahogy ki lehet szamolni az uj folyamot



ut visszafejtese,
szakaszosan konstans fuggvenyek vannak,
ha ezeket vegigmegyek, akkor mi tortenik a vegen


Két egymás követõ élre mi történik? (klasszikus becslések, minél finomabb felosztások)

Nem akkor osztom fel, ha az integrál nagy.
Veszem az alsó és felsõ becslésre a folyamot, ha egy élen nagy az eltérés, akkor körülötte kell felosztani

Minden felosztás legyen az elõzõ finomítása

*/
