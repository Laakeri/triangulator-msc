#include "ghtw.hpp"

#include <vector>
#include <cassert>
#include <map>
#include <iostream>

#include "enumerator.hpp"
#include "graph.hpp"
#include "mcs.hpp"
#include "utils.hpp"
#include "bt_algorithm.hpp"
#include "hypercost.hpp"

namespace triangulator {
bool GHTW::Reduce(Graph graph, const std::vector<std::vector<int>>& hyper_edges) {
  if (graph.m() == 0) return true;
  mcs::McsMOutput minimal_triangulation = mcs::McsM(graph);
  Graph fill_graph = graph;
  fill_graph.AddEdges(minimal_triangulation.fill_edges);
  int ghtw = ChordGHTW(fill_graph, hyper_edges);
  if (minimal_triangulation.fill_edges.size() == 0 || ghtw <= lb_) {
	  lb_ = std::max(lb_, ghtw);
  	utils::Append(fill_, graph.MapBack(minimal_triangulation.fill_edges));
  	return true;
  }
  lb_ = std::max(lb_, 2);
  std::vector<Graph> atoms = mcs::Atoms(graph, minimal_triangulation);
  Log::Write(3, "GHTW atoms ", atoms.size());
  if (atoms.size() == 0) return true;
  if (atoms.size() > 1) {
    for (auto& atom : atoms) {
      atom.InheritMap(graph);
      kernel_.push(atom);
    }
    return true;
  } else {
    atoms[0].InheritMap(graph);
    graph = atoms[0];
  }
  Graph lb_graph = graph;
  mcs::LbTriang(lb_graph);
  int ub = ChordGHTW(lb_graph, hyper_edges);
  if (ub <= lb_) {
    utils::Append(fill_, graph.MapBack(graph.FillEdges(lb_graph)));
    return true;
  }
  kernel_.push(graph);
  return false;
}

Graph GHTW::CompGHTW(Graph graph, const std::vector<std::vector<int>>& hyper_edges) {
  assert(!bad_);
  bad_ = true;
  assert(kernel_.size() == 0);
  assert(lb_ == 0);
  kernel_.push(graph);
  int failed = 0;
  while (failed < (int)kernel_.size()) {
  	Log::Write(3, "GHTWPP reducing ", kernel_.front().n(), " ", failed, " ", (int)kernel_.size(), " ", lb_);
  	if (Reduce(kernel_.front(), hyper_edges)) {
  		failed = 0;
  	} else {
  		failed++;
  	}
  	kernel_.pop();
  }
  std::vector<Graph> insts;
  while (!kernel_.empty()) {
  	insts.push_back(kernel_.front());
  	kernel_.pop();
  }
  auto cmp = [](const Graph& a, const Graph& b) {
  	return a.n() < b.n();
  };
  std::sort(insts.begin(), insts.end(), cmp);
  Log::Write(3, "GHTWPP finished ", insts.size());
  for (const Graph& inst : insts) {
  	Graph lb_graph = inst;
  	mcs::LbTriang(lb_graph);
  	int ub = ChordGHTW(lb_graph, hyper_edges);
		int cost = HardSolve(inst, hyper_edges, ub);
		if (cost == -1) {
      lb_ = std::max(lb_, ub);
			utils::Append(fill_, inst.MapBack(inst.FillEdges(lb_graph)));
		} else {
      lb_ = std::max(lb_, cost);
		}
  }
  graph.AddEdges(fill_);
  Log::Write(3, "lol ", lb_);
  assert(ChordGHTW(graph, hyper_edges) == lb_);
  return graph;
}

int GHTW::HardSolve(const Graph& graph, const std::vector<std::vector<int>>& hyper_edges, int ub) {
	Log::Write(3, "GHTW HARDSOLVE ", graph.n(), " ", graph.m(), " ", lb_, " ", ub);
	
  Enumerator pe;
  pe.SetGHTWFilter(ub-1);
	auto pmcs = pe.AdvancedPmcEnum(graph);
	Log::Write(3, "PMCS ", pmcs.size());

  Timer hct;
  hct.start();
  HyperCost hc(graph, hyper_edges, lb_, ub-1);

  std::vector<std::vector<Bitset>> ord_pmcs(ub);
  for (int i=0;i<(int)pmcs.size();i++){
    if (graph.MaximalIS(pmcs[i]) >= ub) continue;
    int c = hc.Cover(pmcs[i]);
    if (c < ub) {
      ord_pmcs[c].push_back({pmcs[i]});
    }
  }
  pmcs.clear();
  assert(ord_pmcs[0].size() == 0);
  std::vector<int> trs(ub);
  for (int i=1;i<ub;i++){
    for (const auto& pmc : ord_pmcs[i]) {
      pmcs.push_back(pmc);
    }
    ord_pmcs[i].clear();
    trs[i] = pmcs.size();
  }
  hct.stop();
  Log::Write(3, "hypercost time ", hct.getTime().count());

	Timer bt_time;
	bt_time.start();
	BtAlgorithm bt;
	bt.SetId();
	auto sol = bt.Solve(graph, pmcs);
	Log::Write(3, "bt time ", bt_time.getTime().count());
	Log::Write(3, "bt ", sol.first, " ", sol.second.size());

	if (sol.first == -1) {
		return -1;
	}
  bool fo = false;
  for (int i=1;i<ub;i++){
    if (sol.first < trs[i]) {
      sol.first = i;
      fo = true;
      break;
    }
  }
  assert(fo);
  assert(sol.first < ub);
	utils::Append(fill_, graph.MapBack(sol.second));
	return sol.first;
}

int ChordGHTW(const Graph& graph, const std::vector<std::vector<int>>& hyper_edges) {
  if (graph.m() == 0) return 0;
  std::vector<int> order = mcs::Mcs(graph);
  std::vector<int> inv_order = utils::PermInverse(order);
  std::vector<int> nbs(graph.n());
  std::vector<int> no_max(graph.n());
  for (int i=0;i<graph.n();i++){
    nbs[i] = graph.Neighbors(i).size();
  }
  int ans = 0;
  HyperCost hc(graph, hyper_edges, 0, graph.n());
  for (int i = 0; i < graph.n(); i++) {
    int x = order[i];
    int nb = 0;
    std::vector<int> clq;
    for (int nx : graph.Neighbors(x)) {
      if (inv_order[nx] > i) {
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
    if (!no_max[x]) {
      clq.push_back(x);
      ans = std::max(ans, hc.Cover(utils::ToBitset(clq, graph.n())));
    }
  }
  return ans;
}
} // namespace triangulator
