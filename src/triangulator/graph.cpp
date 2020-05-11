#include "graph.hpp"

#include <vector>
#include <algorithm>
#include <set>
#include <cassert>
#include <ostream>
#include <iostream>
#include <queue>

#include "utils.hpp"
#include "matrix.hpp"
#include "bitset.hpp"

namespace triangulator {

Graph::Graph(int n)
  : n_(n), m_(0), adj_list_(n) {
  adj_mat_.Resize(n_, n_);
  adj_mat2_.resize(n_);
  std::vector<int> identity(n);
  for (int i = 0; i < n; i++) {
    identity[i] = i;
    adj_mat2_[i] = Bitset(n);
    adj_mat2_[i].SetTrue(i);
  }
  vertex_map_.Init(identity);
}

Graph::Graph(std::vector<Edge> edges) : vertex_map_(edges) {
  n_ = vertex_map_.Size();
  m_ = 0;
  adj_list_.resize(n_);
  adj_mat_.Resize(n_, n_);
  adj_mat2_.resize(n_);
  for (int i = 0; i < n_; i++) {
    adj_mat2_[i] = Bitset(n_);
    adj_mat2_[i].SetTrue(i);
  }
  for (auto edge : edges) {
    AddEdge(vertex_map_.Rank(edge.first), vertex_map_.Rank(edge.second));
  }
}

int Graph::n() const {
  return n_;
}

int Graph::m() const {
  return m_;
}

bool Graph::HasEdge(int v, int u) const {
  return adj_mat_[v][u];
}

bool Graph::HasEdge(Edge e) const {
  return HasEdge(e.first, e.second);
}

std::vector<Edge> Graph::Edges() const {
  std::vector<Edge> ret;
  for (int i = 0; i < n_; i++) {
    for (int a : adj_list_[i]) {
      if (a > i) ret.push_back({i, a});
    }
  }
  return ret;
}

std::vector<int> Graph::Vertices() const {
  std::vector<int> ret(n_);
  for (int i=0;i<n_;i++){
    ret[i] = i;
  }
  return ret;
}

const std::vector<int>& Graph::Neighbors(int v) const {
  return adj_list_[v];
}

Bitset Graph::Neighbors(const Bitset& vs) const {
  Bitset nbs(n_);
  for (int v : vs) {
    nbs |= adj_mat2_[v];
  }
  nbs.TurnOff(vs);
  return nbs;
}

bool Graph::IsConnected() const {
  auto cs = Components({});
  return (cs.size() == 1) && ((int)cs[0].size() == n_);
}

bool Graph::IsConnectedOrIsolated() const {
  auto cs = Components({});
  int f = 0;
  for (const auto& c : cs) {
    if ((int)c.size() > 1) f++;
  }
  return f <= 1;
}

void Graph::AddEdge(int v, int u) {
  if (HasEdge(v, u)) return;
  assert(v != u);
  m_++;
  adj_list_[v].push_back(u);
  adj_list_[u].push_back(v);
  adj_mat_[v][u] = true;
  adj_mat_[u][v] = true;
  adj_mat2_[v].SetTrue(u);
  adj_mat2_[u].SetTrue(v);
}

void Graph::AddEdge(Edge e) {
  AddEdge(e.first, e.second);
}

void Graph::AddEdges(const std::vector<Edge>& edges) {
  for (auto& edge : edges) AddEdge(edge);
}

void Graph::RemoveEdge(int v, int u) {
  assert(HasEdge(v, u) && HasEdge(u, v));
  m_--;
  adj_mat_[v][u] = false;
  adj_mat_[u][v] = false;
  adj_mat2_[v].SetFalse(u);
  adj_mat2_[u].SetFalse(v);
  int fo = 0;
  for (int i = 0; i < (int)adj_list_[v].size(); i++) {
    if (adj_list_[v][i] == u) {
      std::swap(adj_list_[v][i], adj_list_[v].back());
      adj_list_[v].pop_back();
      fo++;
      break;
    }
  }
  for (int i = 0; i < (int)adj_list_[u].size(); i++) {
    if (adj_list_[u][i] == v) {
      std::swap(adj_list_[u][i], adj_list_[u].back());
      adj_list_[u].pop_back();
      fo++;
      break;
    }
  }
  assert(fo == 2);
}

int Graph::Degeneracy() const {
  std::vector<std::vector<int>> q(n_);
  std::vector<int> dg(n_);
  int vs = 0;
  for (int i=0;i<n_;i++) {
    dg[i] = adj_list_[i].size();
    if (dg[i] > 0) {
      q[dg[i]].push_back(i);
      vs++;
    }
  }
  int mt = 0;
  int t = 0;
  for (int it=0;it<vs;it++) {
    int x = -1;
    while (x == -1) {
      if (q[t].empty()) {
        t++;
      } else {
        x = q[t].back();
        q[t].pop_back();
        if (dg[x] == -1) {
          x = -1;
        } else {
          assert(dg[x] == t);
        }
      }
    }
    mt = std::max(mt, t);
    assert(x>=0&&x<n_&&dg[x]==t&&t>=0);
    for (int nx : adj_list_[x]) {
      if (dg[nx] >= 0) {
        assert(dg[nx] >= 1);
        dg[nx]--;
        dg[x]--;
        q[dg[nx]].push_back(nx);
      }
    }
    assert(dg[x] == 0);
    dg[x] = -1;
    t = std::max(0, t-1);
  }
  return mt;
}

void Graph::Dfs(int v, std::vector<char>& block, std::vector<int>& component) const {
  block[v] = true;
  component.push_back(v);
  for (int nv : adj_list_[v]) {
    if (!block[nv]) {
      Dfs(nv, block, component);
    }
  }
}

std::vector<int> Graph::FindComponentAndMark(int v, std::vector<char>& block) const {
  std::vector<int> component;
  Dfs(v, block, component);
  return component;
}

std::vector<std::vector<int> > Graph::Components(const std::vector<int>& separator) const {
  std::vector<char> blocked(n_);
  for (int v : separator) {
    blocked[v] = true;
  }
  std::vector<std::vector<int> > components;
  for (int i = 0; i < n_; i++) {
    if (!blocked[i]) {
      components.push_back(FindComponentAndMark(i, blocked));
    }
  }
  return components;
}

std::vector<std::vector<int> > Graph::NComponents(const std::vector<int>& separator) const {
  std::vector<char> blocked(n_);
  for (int v : separator) {
    blocked[v] = true;
  }
  std::vector<std::vector<int> > components;
  for (int v : separator) {
    for (int nv : adj_list_[v]) {
      if (!blocked[nv]) {
        components.push_back(FindComponentAndMark(nv, blocked));
      }
    }
  }
  return components;
}

// Returns the vector in sorted order
std::vector<int> Graph::Neighbors(const std::vector<int>& vs) const {
  std::vector<int> neighbors;
  neighbors.reserve(vs.size());
  int sum = 0;
  for (int v : vs) {
    sum += adj_list_[v].size();
  }
  if (sum >= n_/4) { // Two cases for optimization
    std::vector<char> nbs(n_);
    for (int v : vs) {
      for (int nv : adj_list_[v]) {
        nbs[nv] = true;
      }
    }
    for (int v : vs) {
      nbs[v] = false;
    }
    for (int i = 0; i < n_; i++) {
      if (nbs[i]) neighbors.push_back(i);
    }
  } else {
    std::set<int> nbs;
    for (int v : vs) {
      for (int nv : adj_list_[v]) {
        nbs.insert(nv);
      }
    }
    for (int v : vs) {
      nbs.erase(v);
    }
    for (int v : nbs) {
      neighbors.push_back(v);
    }
  }
  return neighbors;
}

std::vector<Edge> Graph::EdgesIn(const std::vector<int>& vs) const {
  static std::vector<char> is;
  utils::InitZero(is, n_);
  for (int v : vs) {
    is[v] = true;
  }
  std::vector<Edge> edges;
  for (int v : vs) {
    if (adj_list_[v].size() <= vs.size()) { // Two cases for optimization
      for (int nv : adj_list_[v]) {
        if (is[nv] && nv > v) edges.push_back({v, nv});
      }
    }
    else {
      for (int nv : vs) {
        if (adj_mat_[v][nv] && nv > v) edges.push_back({v, nv});
      }
    }
  }
  return edges;
}

bool Graph::IsClique(const std::vector<int>& clique) const {
  for (int i = 0; i < (int)clique.size(); i++) {
    for (int ii = i + 1; ii < (int)clique.size(); ii++) {
      if (!HasEdge(clique[i], clique[ii])) return false;
    }
  }
  return true;
}

bool Graph::IsFull(int v, Bitset sep, Bitset vis) const {
  static std::vector<int> q;
  q.resize(n_);
  q[0] = v;
  vis.SetFalse(v);
  int i = 0;
  int s = 1;
  int chunks = vis.Chunks();
  while (i < s) {
    int x = q[i++];
    for (int j = 0; j < chunks; j++) {
      uint64_t go = vis.data_[j] & adj_mat2_[x].data_[j];
      while (go) {
        vis.data_[j] &= (~(go&-go));
        q[s++] = __builtin_ctzll(go) + j*BITS;
        go &= ~-go;
      }
    }
    bool has = false;
    for (int j = 0; j < chunks; j++) {
      sep.data_[j] &= (~adj_mat2_[x].data_[j]);
      if (sep.data_[j]) has = true;
    }
    if (!has) return true;
  }
  return false;
}

bool Graph::IsFull2(int v, Bitset sep, Bitset& vis) const {
  static std::vector<int> q;
  q.resize(n_);
  q[0] = v;
  vis.SetFalse(v);
  int i = 0;
  int s = 1;
  int chunks = vis.Chunks();
  while (i < s) {
    int x = q[i++];
    for (int j = 0; j < chunks; j++) {
      uint64_t go = vis.data_[j] & adj_mat2_[x].data_[j];
      while (go) {
        vis.data_[j] &= (~(go&-go));
        q[s++] = __builtin_ctzll(go) + j*BITS;
        go &= ~-go;
      }
    }
    bool has = false;
    for (int j = 0; j < chunks; j++) {
      sep.data_[j] &= (~adj_mat2_[x].data_[j]);
      if (sep.data_[j]) has = true;
    }
    if (!has) return true;
  }
  return false;
}

void Graph::Dfs22(int v, Bitset& sep, Bitset& vis, std::vector<int>& f, const Bitset& good) const {
  vis.SetTrue(v);
  int chunks = vis.Chunks();
  Bitset ne(n_);
  ne.SetTrue(v);
  bool fo = true;
  while (fo) {
    fo = false;
    for (int j = 0; j < chunks; j++) {
      uint64_t gv = vis.data_[j] & ne.data_[j];
      while (gv) {
        fo = true;
        vis.data_[j] &= (~(gv&-gv));
        int x = __builtin_ctzll(gv) + j*BITS;
        ne |= adj_mat2_[x];
        gv &= ~-gv;
      }
      uint64_t gs = sep.data_[j] & ne.data_[j];
      while (gs) {
        sep.data_[j] &= (~(gs&-gs));
        if (good.data_[j] & (gs&-gs)) {
          f.push_back(__builtin_ctzll(gs) + j*BITS);
        }
        gs &= ~-gs;
      }
    }
  }
}

void Graph::Dfs2(int v, Bitset& sep, Bitset& vis, std::vector<int>& f) const {
  vis.SetTrue(v);
  int chunks = vis.Chunks();
  Bitset ne(n_);
  ne.SetTrue(v);
  bool fo = true;
  while (fo) {
    fo = false;
    for (int j = 0; j < chunks; j++) {
      uint64_t gv = vis.data_[j] & ne.data_[j];
      while (gv) {
        fo = true;
        vis.data_[j] &= (~(gv&-gv));
        int x = __builtin_ctzll(gv) + j*BITS;
        ne |= adj_mat2_[x];
        gv &= ~-gv;
      }
      uint64_t gs = sep.data_[j] & ne.data_[j];
      while (gs) {
        sep.data_[j] &= (~(gs&-gs));
        f.push_back(__builtin_ctzll(gs) + j*BITS);
        gs &= ~-gs;
      }
    }
  }
}

std::vector<Bitset> Graph::NComponents(const Bitset& separator) const {
  Bitset vis(n_);
  vis.FillUpTo(n_);
  vis.TurnOff(separator);
  Bitset nbs = Neighbors(separator);
  std::vector<Bitset> ret;
  for (Bitset comp : BitComps(vis)) {
    if (comp.Intersects(nbs)) {
      ret.push_back(comp);
    }
  }
  return ret;
}

std::vector<Bitset> Graph::BitComps(Bitset vis) const {
  Bitset ne(n_);
  int chunks = vis.Chunks();
  std::vector<Bitset> ret;
  bool fo = false;
  while (1) {
    if (!fo) {
      for (int j = 0; j < chunks; j++) {
        if (vis.data_[j]) {
          int x = __builtin_ctzll(vis.data_[j]) + j*BITS;
          ne.SetTrue(x);
          ret.push_back(Bitset(n_));
          fo = true;
          break;
        }
      }
      if (!fo) return ret;
    }
    fo = false;
    for (int j = 0; j < chunks; j++) {
      uint64_t gv = vis.data_[j] & ne.data_[j];
      while (gv) {
        fo = true;
        vis.data_[j] &= (~(gv&-gv));
        int x = __builtin_ctzll(gv) + j*BITS;
        ne |= adj_mat2_[x];
        ret.back().SetTrue(x);
        gv &= ~-gv;
      }
    }
  }
}

std::vector<std::vector<int>> Graph::CompNeighs(const std::vector<int>& block) const {
  Bitset sep(n_);
  Bitset vis(n_);
  for (int i = 0; i < n_; i++) {
    if (adj_list_[i].size() > 0) {
      vis.SetTrue(i);
    }
  }
  for (int x : block) {
    sep.SetTrue(x);
    vis.SetFalse(x);
  }
  Bitset sb = sep;
  std::vector<std::vector<int>> ret;
  std::vector<int> f;
  f.reserve(block.size());
  for (int i = 0; i < n_; i++) {
    if (vis.Get(i)) {
      Dfs2(i, sep, vis, f);
      if (!f.empty()) {
        ret.push_back(f);
        f.clear();
      }
      sep = sb;
    }
  }
  return ret;
}

std::vector<int> Graph::CompNeigh(const std::vector<int>& block, int v) const {
  Bitset sep = utils::ToBitset(block, n_);
  assert(!sep.Get(v));
  Bitset vis = ~sep;
  std::vector<int> f;
  f.reserve(block.size());
  Dfs2(v, sep, vis, f);
  return f;
}

bool Graph::HasNFullComponents(const std::vector<int>& separator, int n) const {
  return HasNFullComponents(utils::ToBitset(separator, n_), n);
}

bool Graph::HasNFullComponents(const Bitset& separator, int n) const {
  Bitset vis(n_);
  vis.FillTrue();
  vis.TurnOff(separator);
  std::vector<int> f;
  f.reserve(separator.Popcount());
  int cnt = 0;
  Bitset ne(n_);
  for (int i = 0; i < n_; i++) {
    if (vis.Get(i) && adj_list_[i].size() > 0) {
      ne.SetTrue(i);
      Dfs2Bit(vis, ne);
      if (ne.Subsumes(separator)) {
        cnt++;
      }
      if (cnt >= n) return true;
      ne.Clear();
    }
  }
  return false;
}

bool Graph::IsMinsep(const Bitset& separator) const {
  return HasNFullComponents(separator, 2);
}

bool Graph::IsMinsep(const std::vector<int>& separator) const {
  return HasNFullComponents(utils::ToBitset(separator, n_), 2);
}

bool Graph::IsMinsep(const std::vector<int>& separator, int a, int b) const {
  if (a == b) return false;
  Bitset sep = utils::ToBitset(separator, n_);
  if (sep.Get(a) || sep.Get(b)) return false;
  Bitset vis = ~sep;
  std::vector<int> f;
  f.reserve(separator.size());
  Dfs2(a, sep, vis, f);
  if (f.size() != separator.size()) return false;
  if (!vis.Get(b)) return false;
  sep.SetTrue(separator);
  f.clear();
  Dfs2(b, sep, vis, f);
  return f.size() == separator.size();
}

Bitset Graph::AnotherComp(int x, const Bitset& minsep) const {
  // could be optimized
  int fulls = 0;
  bool ffx = false;
  bool fnx = false;
  Bitset ret;
  Bitset vis(n_);
  vis.FillUpTo(n_);
  vis.TurnOff(minsep);
  for (auto comp : BitComps(vis)) {
    if (Neighbors(comp).Popcount() == minsep.Popcount()) {
      fulls++;
      if (!comp.Get(x)) {
        fnx = true;
        ret = comp;
      } else {
        ffx = true;
      }
    }
  }
  assert(fulls == 2);
  assert(ffx && fnx);
  return ret;
}

std::vector<Edge> Graph::FillEdges(const std::vector<int>& clq) const {
  std::vector<Edge> ret;
  for (int i=0;i<(int)clq.size();i++){
    for (int ii=i+1;ii<(int)clq.size();ii++){
      if (!HasEdge(clq[i], clq[ii])) {
        ret.push_back(std::minmax(clq[i], clq[ii]));
      }
    }
  }
  return ret;
}

std::vector<Edge> Graph::FillEdges(const Graph& other) const {
  std::vector<Edge> ret;
  for (auto e : other.Edges()) {
    if (!HasEdge(e)) {
      ret.push_back(e);
    }
  }
  return ret;
}

std::vector<Edge> Graph::FillEdges(Bitset bs) const {
  int chunks = bs.chunks_;
  std::vector<Edge> ret;
  for (int i=0;i<chunks;i++){
    while (bs.data_[i]) {
      int v = i*BITS + __builtin_ctzll(bs.data_[i]);
      bs.data_[i] &= ~-bs.data_[i];
      for (int j=i;j<chunks;j++){
        uint64_t td = bs.data_[j] & (~adj_mat2_[v].data_[j]);
        while (td) {
          int u = j*BITS + __builtin_ctzll(td);
          td &= ~-td;
          ret.push_back({v, u});
        }
      }
    }
  }
  return ret;
}

int Graph::FillSize(Bitset bs) const {
  int chunks = bs.chunks_;
  int ans = 0;
  for (int i=0;i<chunks;i++){
    while (bs.data_[i]) {
      int v = i*BITS + __builtin_ctzll(bs.data_[i]);
      bs.data_[i] &= ~-bs.data_[i];
      for (int j=i;j<chunks;j++){
        ans += __builtin_popcountll(bs.data_[j] & (~adj_mat2_[v].data_[j]));
      }
    }
  }
  return ans;
}

void Graph::FillBS(Bitset bs) {
  int chunks = bs.chunks_;
  for (int i=0;i<chunks;i++){
    while (bs.data_[i]) {
      int v = i*BITS + __builtin_ctzll(bs.data_[i]);
      bs.data_[i] &= ~-bs.data_[i];
      for (int j=i;j<chunks;j++){
        uint64_t td = bs.data_[j] & (~adj_mat2_[v].data_[j]);
        while (td) {
          int u = j*BITS + __builtin_ctzll(td);
          td &= ~-td;
          AddEdge(v, u);
        }
      }
    }
  }
}

bool Graph::IsClique(Bitset bs) const {
  int chunks = bs.chunks_;
  for (int i=0;i<chunks;i++){
    while (bs.data_[i]) {
      int v = i*BITS + __builtin_ctzll(bs.data_[i]);
      bs.data_[i] &= ~-bs.data_[i];
      for (int j=i;j<chunks;j++){
        if (bs.data_[j] & (~adj_mat2_[v].data_[j])) {
          return false;
        }
      }
    }
  }
  return true;
}

void Graph::RemoveEdgesBetween(int v, const std::vector<int>& vs) {
  for (int u : vs) {
    assert(u != v);
    RemoveEdge(u, v);
  }
}

bool Graph::IsAlmostClique(const std::vector<int>& clq) const {
  std::vector<int> rm;
  for (int i=0;i<(int)clq.size();i++){
    for (int ii=i+1;ii<(int)clq.size();ii++){
      if (!HasEdge(clq[i], clq[ii])) {
        if (rm.size() == 0) {
          rm = {clq[i], clq[ii]};
        } else if (rm.size() == 1) {
          if (clq[i] != rm[0] && clq[ii] != rm[0]) {
            return false;
          }
        } else if (rm.size() == 2) {
          if (clq[i] == rm[0]) {
            assert(clq[ii] != rm[1]);
            rm = {rm[0]};
          } else if(clq[i] == rm[1]) {
            assert(clq[ii] != rm[0]);
            rm = {rm[1]};
          } else if(clq[ii] == rm[0]) {
            assert(clq[i] != rm[1]);
            rm = {rm[0]};
          } else if(clq[ii] == rm[1]) {
            assert(clq[i] != rm[0]);
            rm = {rm[1]};
          } else {
            return false;
          }
        } else {
          assert(0);
        }
      }
    }
  }
  return true;
}

std::vector<int> Graph::Distances(const std::vector<int>& start) const {
  assert(start.size() > 0);
  std::vector<int> d(n_);
  for (int i=0;i<n_;i++){
    d[i] = n_+1;
  }
  std::queue<int> q;
  for (int v : start) {
    d[v] = 0;
    q.push(v);
  }
  while (!q.empty()) {
    int x = q.front();
    q.pop();
    for (int nx : adj_list_[x]){
      if (d[nx] == n_+1) {
        d[nx] = d[x]+1;
        q.push(nx);
      }
    }
  }
  return d;
}

std::vector<std::vector<int>> Graph::DistanceMatrix() const {
  std::vector<std::vector<int>> ret;
  for (int i=0;i<n_;i++){
    ret.push_back(Distances({i}));
    assert((int)ret.back().size() == n_);
  }
  return ret;
}

int Graph::MapBack(int v) const {
  return vertex_map_.Kth(v);
}
std::vector<int> Graph::MapBack(std::vector<int> vs) const {
  for (int& v : vs) {
    v = MapBack(v);
  }
  return vs;
}
std::pair<int, int> Graph::MapBack(int v, int u) const {
  return {MapBack(v), MapBack(u)};
}
int Graph::MapInto(int v) const {
  return vertex_map_.Rank(v);
}
std::vector<int> Graph::MapInto(std::vector<int> vs) const {
  for (int& v : vs) {
    v = MapInto(v);
  }
  return vs;
}
std::set<std::pair<int, int>> Graph::MapInto(std::set<std::pair<int, int>> vs) const {
  std::set<std::pair<int, int>> ret;
  for (auto v : vs) {
    ret.insert(std::minmax(MapInto(v.first), MapInto(v.second)));
  }
  return ret;
}
Edge Graph::MapBack(Edge e) const {
  return {MapBack(e.first), MapBack(e.second)};
}
std::vector<Edge> Graph::MapBack(std::vector<Edge> es) const {
  for (Edge& e : es) {
    e = MapBack(e);
  }
  return es;
}
void Graph::InheritMap(const Graph& parent) {
  vertex_map_ = StaticSet<int>(parent.MapBack(vertex_map_.Values()));
}

std::vector<Bitset> Graph::CompNeighsBit(const Bitset& block) const {
  Bitset vis = ~block;
  std::vector<Bitset> ret;
  Bitset ne(n_);
  Bitset sep(n_);
  for (int i=0;i<n_;i++){
    if (vis.Get(i)) {
      ne.CopyFrom(adj_mat2_[i]);
      Dfs2Bit(vis, ne);
      sep.SetAnd(block, ne);
      if (sep.Popcount() > 0) {
        ret.push_back(sep);
      }
    }
  }
  return ret;
}


void Graph::Dfs2Bit(Bitset& vis, Bitset& ne) const {
  int chunks = vis.Chunks();
  bool fo = true;
  while (fo) {
    fo = false;
    for (int j = 0; j < chunks; j++) {
      uint64_t gv = vis.data_[j] & ne.data_[j];
      while (gv) {
        fo = true;
        vis.data_[j] &= (~(gv&-gv));
        int x = __builtin_ctzll(gv) + j*BITS;
        gv &= ~-gv;
        ne |= adj_mat2_[x];
      }
    }
  }
}

std::vector<Bitset> Graph::AllMinseps() const {
  assert(IsConnectedOrIsolated());
  Timer mstimer;
  mstimer.start();
  std::vector<Bitset> minseps;
  BitsetSet ff(n_, n_, 2);
  for (int i = 0; i < n_; i++) {
    if (Neighbors(i).empty()) continue;
    for (const Bitset& nbs : CompNeighsBit(adj_mat2_[i])) {
      if (ff.Insert(nbs)) {
        minseps.push_back(nbs);
      }
    }
  }
  Bitset block(n_);
  Bitset vis(n_);
  Bitset sep(n_);
  Bitset ne(n_);
  Bitset mask(n_);
  for (int i=0;i<n_;i++){
    if (!Neighbors(i).empty()) mask.SetTrue(i);
  }
  int chunks = vis.Chunks();
  for (int i = 0;i<(int)minseps.size(); i++) {
    if (i%10000 == 0 || (i < 10000 && i%1000 ==0) || (i < 1000 && i%100==0) || (i<100 && i%10==0)) Log::Write(3, "allminseps ", i, " ", minseps.size(), " ", (double)i/(double)minseps.size());
    for (int j=0;j<n_;j++){
      if (!minseps[i].Get(j)) continue;
      block.CopyFrom(minseps[i]);
      block |= adj_mat2_[j];
      vis.SetNegAnd(block, mask);
      for (int ch = 0; ch < chunks; ch++) {
        while (vis.data_[ch] > 0) {
          int k = __builtin_ctzll(vis.data_[ch]) + ch*BITS;
          sep.CopyFrom(block);
          ne.CopyFrom(adj_mat2_[k]);
          Dfs2Bit(vis, ne);
          sep.SetAnd(ne, block);
          if (ff.Insert(sep)) {
            minseps.push_back(sep);
          }
        }
      }
    }
  }
  Log::Write(3, "MSTIME ", mstimer.getTime().count(), " ", ff.ContainerSize());
  return minseps;
}

int Graph::MaximalIS(const Bitset& vs) const {
  Bitset is(n_);
  int ans = 0;
  for (int v : vs) {
    if (!is.Intersects(adj_mat2_[v])) {
      is.SetTrue(v);
      ans++;
    }
  }
  return ans;
}
} // namespace triangulator                                           
