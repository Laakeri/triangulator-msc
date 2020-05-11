#include <bits/stdc++.h>
using namespace std;

int getrand(int a, int b, mt19937& gen) {
	return uniform_int_distribution<int>(a,b)(gen);
}

std::vector<std::vector<int>> readPmcs(ifstream& in, int k) {
	std::vector<std::vector<int>> ret;
	string tmp;
	while (getline(in, tmp)) {
		stringstream ss;
		ss<<tmp;
		ss>>tmp;
		if (tmp == "Pmc:") {
			ret.push_back({});
			int v;
			while (ss >> v) {
				assert(v>=0);
				ret.back().push_back(v);
			}
			if ((int)ret.back().size() > k) {
				ret.pop_back();
			}
		} else {
			break;
		}
	}
	return ret;
}

vector<int> g[101010];
int u[101010];

int dfs(int x) {
	if (u[x]) {
		return 0;
	}
	u[x] = 1;
	int r =1;
	for (int nx : g[x]) {
		r += dfs(nx);
	}
	return r;
}

int main(int argc, char** argv) {
	assert(argc == 4);
	string tt(argv[1]);
	int max_n = stoi(argv[2]);
	mt19937 gen(stoi(argv[3]));

	assert(tt == "cpmc" || tt == "brute" || tt == "basic" || tt == "ws");
	assert(max_n > 2 && max_n < 200);

	cout<<"testing: "<<tt<<" "<<max_n<<" "<<stoi(argv[3])<<endl;

	int ti=0;
	while (1) {
		ti++;
		int n = getrand(0, max_n, gen);
		int m = getrand(0, n*n, gen);
		int k = getrand(2, max(2, n), gen);
		if (tt != "ws") {
			k = n+1;
		}
		for (int i=0;i<=n;i++){
			g[i].clear();
			u[i] = 0;
		}
		set<pair<int, int>> es;
		for (int i=0;i<m;i++){
			int a = getrand(0, n-1, gen);
			int b = getrand(0, n-1, gen);
			if (a != b) {
				es.insert(minmax(a, b));
				g[a].push_back(b);
				g[b].push_back(a);
			}
		}
		if (dfs(0) != n) continue;
		ofstream out("tinput");
		vector<pair<int, int>> el;
		out<<"p "<<n<<" "<<es.size()<<endl;
		for (auto e : es) {
			el.push_back({e.first, e.second});
		}
		shuffle(el.begin(), el.end(), gen);
		for (auto e : el) {
			assert(e.first<e.second);
			assert(e.first >= 0 && e.second < n);
			if (getrand(0, 1, gen) == 0) {
				swap(e.first, e.second);
			}
			out<<"e "<<e.first<<" "<<e.second<<endl;
		}
		out.close();
		if (tt == "brute") {
			assert(system("./triangulator brute -pa < tinput > bout 2>/dev/null") == 0);
			assert(system("./triangulator basic -pa < tinput > out 2>/dev/null") == 0);
		} else if (tt == "basic") {
			assert(system("./triangulator count_pmc_basic < tinput > bout 2>/dev/null") == 0);
			assert(system("./triangulator count_pmc < tinput > out 2>/dev/null") == 0);
		} else if (tt == "ws") {
			assert(system("./triangulator basic -pa < tinput > bout 2>/dev/null") == 0);
			assert(system(("./triangulator ws -pa -k=" + to_string(k) + " < tinput > out 2>/dev/null").c_str()) == 0);
			cout<<"k="<<k<<endl;
		} else if (tt == "cpmc") {
			assert(system("./triangulator count_pmc < tinput > bout 2>/dev/null") == 0);
			cout<<"ok "<<ti<<endl;
			continue;
		} else {
			assert(false);
		}
		ifstream in1("out");
		ifstream in2("bout");
		int pmc1, pmc2;
		in1>>pmc1;
		in2>>pmc2;
		assert(pmc1 == pmc2);
		cout<<"ok "<<ti<<" "<<pmc1<<endl;
		/*
		ifstream in1("bout");
		auto pmcs1 = readPmcs(in1, k);
		in1.close();
		ifstream in2("out");
		auto pmcs2 = readPmcs(in2, k);
		in2.close();

		assert(pmcs1 == pmcs2);
		cout<<"ok "<<ti<<" "<<pmcs1.size()<<endl;
		*/
	}
}