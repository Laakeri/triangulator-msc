#include "mcs.hpp"

#include <vector>
#include <cassert>
#include <set>
#include <queue>

#include "utils.hpp"
#include "graph.hpp"

namespace triangulator {
namespace mcs {

std::vector<int> Mcs(const Graph& graph) {
  std::vector<int> order(graph.n());
  static std::vector<int> label;
  static std::vector<char> rm;
  static std::vector<std::vector<int> > labels;
  utils::InitZero(label, graph.n());
  utils::InitZero(rm, graph.n());
  if ((int)labels.size() < graph.n()) labels.resize(graph.n());
  for (int i = 0; i < graph.n(); i++) labels[i].clear();
  for (int i = 0; i < graph.n(); i++) labels[0].push_back(i);
  int max_label = 0;
  for (int it = graph.n() - 1; it >= 0; it--) {
    if (labels[max_label].size() == 0) {
      max_label--;
      it++;
      continue;
    }
    int x = labels[max_label].back();
    labels[max_label].pop_back();
    if (rm[x]) {
      it++;
      continue;
    }
    order[it] = x;
    for (int nx : graph.Neighbors(x)) {
      if (!rm[nx]) {
        label[nx]++;
        labels[label[nx]].push_back(nx);
        max_label = std::max(max_label, label[nx]);
      }
    }
    rm[x] = true;
  }
  return order;
}

McsMOutput McsM(const Graph& graph) {
  // these vectors are static to prevent reallocating memory every time this function is called
  Timer mcst;
  mcst.start();
  static std::vector<int> label;
  static std::vector<std::vector<int> > reach;
  Bitset rm(graph.n());
  Bitset rc(graph.n());
  utils::InitZero(label, graph.n());
  if ((int)reach.size() < graph.n()) reach.resize(graph.n());
  for (int i = 0; i < graph.n(); i++) reach[i].clear();
  std::vector<Edge> fill;
  std::vector<int> order(graph.n());
  std::vector<char> is_maximal_point(graph.n());
  // TODO: maybe better variable names?
  int chunks = rm.chunks_;
  int prev_label = -1;
  for (int it = graph.n() - 1; it >= 0; it--) {
    int x = 0;
    int max_label = 0;
    for (int i = 0; i < graph.n(); i++) {
      if (!rm.Get(i) && label[i] >= max_label) {
        x = i;
        max_label = label[x];
      }
    }
    assert(!rm.Get(x) && label[x] < graph.n());
    order[it] = x;
    is_maximal_point[it] = (label[x] <= prev_label);
    prev_label = label[x];
    rc.Clear();
    rc.SetTrue(x);
    rm.SetTrue(x);
    for (int y : graph.Neighbors(x)) {
      if (!rm.Get(y)) {
        rc.SetTrue(y);
        reach[label[y]].push_back(y);
      }
    }
    for (int i = 0; i < graph.n(); i++) {
      while (!reach[i].empty()) {
        int y = reach[i].back();
        reach[i].pop_back();
        for (int j=0;j<chunks;j++){
          uint64_t td = graph.adj_mat2_[y].data_[j] & (~rm.data_[j]) & (~rc.data_[j]);
          while (td) {
            int z = j*BITS + __builtin_ctzll(td);
            td &= ~-td;
            rc.SetTrue(z);
            if (label[z] > i) {
              reach[label[z]].push_back(z);
              label[z]++;
              fill.push_back({x, z});
            } else {
              reach[i].push_back(z);
            }
          }
        }
      }
    }
    for (int y : graph.Neighbors(x)) {
      if (!rm.Get(y)) label[y]++;
    }
  }
  return McsMOutput({fill, order, is_maximal_point});
}

std::vector<Graph> Atoms(const Graph& graph, const McsMOutput& mcs_m_output) {
  if (graph.n() == 0 || graph.m() == 0) {
    return {};
  }
  Graph filled_graph(graph);
  filled_graph.AddEdges(mcs_m_output.fill_edges);
  // Use static to not allocate memory every time this is called
  static std::vector<char> rm, block;
  utils::InitZero(rm, graph.n());
  utils::InitZero(block, graph.n());
  std::vector<Graph> atoms;
  for (int it = 0; it < graph.n(); it++) {
    int x = mcs_m_output.elimination_order[it];
    if (mcs_m_output.is_maximal_clique_point[it]) {
      std::vector<int> cand_clique;
      for (int nx : filled_graph.Neighbors(x)) {
        if (!rm[nx]) cand_clique.push_back(nx);
      }
      if (graph.IsClique(cand_clique)) {
        for (int y : cand_clique) block[y] = 1;
        std::vector<int> component = graph.FindComponentAndMark(x, block);
        for (int y : cand_clique) block[y] = 0;
        component.insert(component.end(), cand_clique.begin(), cand_clique.end());
        atoms.push_back(Graph(graph.EdgesIn(component)));
      }
    }
    rm[x] = true;
  }
  int found = 0; // TODO
  for (int x = 0; x < graph.n(); x++) {
    if (!block[x]) {
      std::vector<int> component = graph.FindComponentAndMark(x, block);
      atoms.push_back(Graph(graph.EdgesIn(component)));
      found++;
    }
  }
  assert(found == 1); // TODO
  return atoms;
}

int Heur(const Graph& graph, int v) {
  std::set<Edge> fes;
  for (const Bitset& cn : graph.CompNeighsBit(graph.adj_mat2_[v])) {
    for (auto fe : graph.FillEdges(cn.Elements())) {
      fes.insert(fe);
    }
  }
  return fes.size();
}

void LbTriang(Graph& graph) {
  Timer lbt;
  lbt.start();
  std::priority_queue<std::pair<int, int>> q;
  std::vector<int> hs(graph.n());
  for (int i=0;i<graph.n();i++) {
    hs[i] = Heur(graph, i);
    q.push({-hs[i], i});
  }
  while (!q.empty()) {
    int fi = -q.top().first;
    int v = q.top().second;
    q.pop();
    if (hs[v] != fi) continue;
    hs[v] = -1;
    std::set<int> upd1;
    for (const Bitset& cn : graph.CompNeighsBit(graph.adj_mat2_[v])) {
      for (auto fe : graph.FillEdges(cn.Elements())) {
        graph.AddEdge(fe);
        upd1.insert(fe.first);
        upd1.insert(fe.second);
      }
    }
    std::set<int> upd2;
    for (int u : upd1) {
      for (int nb : graph.Neighbors(u)) {
        if (hs[nb] != -1) {
          upd2.insert(nb);
        }
      }
    }
    for (int u : upd2) {
      assert(hs[u] != -1);
      hs[u] = Heur(graph, u);
      q.push({-hs[u], u});
    }
  }
}

void FastMT(Graph& graph) {
  Timer fmtt;
  fmtt.start();
  for (int i = 0; i < graph.n(); i++) {
    for (const Bitset& cn : graph.CompNeighsBit(graph.adj_mat2_[i])) {
      graph.FillBS(cn);
    }
  }
  fmtt.stop();
}

int Treewidth(const Graph& graph) {
  if (graph.m() == 0) return 0;
  std::vector<int> order = Mcs(graph);
  std::vector<int> inv_order = utils::PermInverse(order);
  int treewidth = 0;
  for (int i = 0; i < graph.n(); i++) {
    int x = order[i];
    int nb = 0;
    std::vector<int> clq;
    for (int nx : graph.Neighbors(x)) {
      if (inv_order[nx] > i) {
        nb++;
        clq.push_back(nx);
      }
    }
    treewidth = std::max(treewidth, nb);
  }
  return treewidth;
}

int64_t TotalTableSize(const Graph& graph, const std::vector<int64_t>& domains) {
  if (graph.m() == 0) return 0;
  std::vector<int> order = Mcs(graph);
  std::vector<int> inv_order = utils::PermInverse(order);
  int64_t ans = 0;
  std::vector<int> nbs(graph.n());
  std::vector<int> no_max(graph.n());
  for (int i=0;i<graph.n();i++){
    nbs[i] = graph.Neighbors(i).size();
  }
  for (int i = 0; i < graph.n(); i++) {
    int x = order[i];
    int nb = 0;
    std::vector<int> clq;
    int64_t cost = domains[graph.MapBack(x)];
    for (int nx : graph.Neighbors(x)) {
      if (inv_order[nx] > i) {
        cost *= domains[graph.MapBack(nx)];
        assert(nbs[nx] >= nbs[x]);
        if (nbs[nx] == nbs[x]) {
          no_max[nx] = 1;
        }
        nbs[nx]--;
        nb++;
        clq.push_back(nx);
      }
    }
    assert(nb == nbs[x]);
    if (!no_max[x]) ans += cost;
  }
  return ans;
}

int Treelength(const Graph& graph, const std::vector<std::vector<int>>& dists) {
  if (graph.m() == 0) return 0;
  std::vector<int> order = Mcs(graph);
  std::vector<int> inv_order = utils::PermInverse(order);
  std::vector<int> nbs(graph.n());
  int nn = 0;
  for (int i=0;i<graph.n();i++){
    nbs[i] = graph.Neighbors(i).size();
    nn = std::max(nn, graph.MapBack(i));
  }
  assert((int)dists.size() > nn);
  for (int i=0;i<(int)dists.size();i++){
    assert(dists[i].size() == dists.size());
  }
  int ans = 0;
  for (int i = 0; i < graph.n(); i++) {
    int x = order[i];
    int nb = 0;
    std::vector<int> clq;
    for (int nx : graph.Neighbors(x)) {
      if (inv_order[nx] > i) {
        assert(nbs[nx] >= nbs[x]);
        nbs[nx]--;
        nb++;
        clq.push_back(nx);
      }
    }
    assert(nb == nbs[x]);
    clq.push_back(x);
    for (int j=0;j<(int)clq.size();j++) {
      for (int jj=j+1;jj<(int)clq.size();jj++) {
        ans = std::max(ans, dists[graph.MapBack(clq[j])][graph.MapBack(clq[jj])]);
      }
    }
  }
  return ans;
}
} // namespace mcs
} // namespace triangulator
