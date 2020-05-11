#include "ftw.hpp"

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
bool FHTW::Reduce(Graph graph, const std::vector<std::vector<int>>& hyper_edges) {
  if (graph.m() == 0) return true;
  mcs::McsMOutput minimal_triangulation = mcs::McsM(graph);
  Graph fill_graph = graph;
  fill_graph.AddEdges(minimal_triangulation.fill_edges);
  double fhtw = ChordFHTW(fill_graph, hyper_edges);
  if (minimal_triangulation.fill_edges.size() == 0 || fhtw <= lb_) {
	  lb_ = std::max(lb_, fhtw);
  	utils::Append(fill_, graph.MapBack(minimal_triangulation.fill_edges));
  	return true;
  }
  std::vector<Graph> atoms = mcs::Atoms(graph, minimal_triangulation);
  Log::Write(3, "FHTW atoms ", atoms.size());
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
  double ub = ChordFHTW(lb_graph, hyper_edges);
  if (ub <= lb_) {
    utils::Append(fill_, graph.MapBack(graph.FillEdges(lb_graph)));
    return true;
  }
  kernel_.push(graph);
  return false;
}

Graph FHTW::CompFHTW(Graph graph, const std::vector<std::vector<int>>& hyper_edges) {
  assert(!bad_);
  bad_ = true;
  assert(kernel_.size() == 0);
  assert(lb_ == 0);
  kernel_.push(graph);
  int failed = 0;
  while (failed < (int)kernel_.size()) {
  	Log::Write(3, "FHTWPP reducing ", kernel_.front().n(), " ", failed, " ", (int)kernel_.size(), " ", lb_);
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
  Log::Write(3, "FHTWPP finished ", insts.size());
  for (const Graph& inst : insts) {
  	Graph lb_graph = inst;
  	mcs::LbTriang(lb_graph);
  	double ub = ChordFHTW(lb_graph, hyper_edges);
		double cost = HardSolve(inst, hyper_edges, ub);
		if (cost <= -0.5) {
      lb_ = std::max(lb_, ub);
			utils::Append(fill_, inst.MapBack(inst.FillEdges(lb_graph)));
		} else {
      lb_ = std::max(lb_, cost);
		}
  }
  graph.AddEdges(fill_);
  Log::Write(3, "lol ", lb_);
  assert(std::abs(ChordFHTW(graph, hyper_edges) - lb_) < 0.01);
  return graph;
}

double FHTW::HardSolve(const Graph& graph, const std::vector<std::vector<int>>& hyper_edges, double ub) {
	Log::Write(3, "FHTW HARDSOLVE ", graph.n(), " ", graph.m(), " ", lb_, " ", ub);
	
  Enumerator pe;
	auto pmcs = pe.AdvancedPmcEnum(graph);
	Log::Write(3, "PMCS ", pmcs.size());

  Timer hct;
  hct.start();
  FracHyperCost fhc(graph, hyper_edges, lb_, ub);

  std::vector<std::pair<double, Bitset>> ord_pmcs;
  for (int i=0;i<(int)pmcs.size();i++){
    double c = fhc.Cover(pmcs[i]);
    if (c < ub) {
      ord_pmcs.push_back({c, pmcs[i]});
    }
  }
  pmcs.clear();
  std::sort(ord_pmcs.begin(), ord_pmcs.end());
  for (const auto& op : ord_pmcs) {
    pmcs.push_back(op.second);
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
  assert(sol.first < (int)ord_pmcs.size());
  double ret = ord_pmcs[sol.first].first;
	utils::Append(fill_, graph.MapBack(sol.second));
	return ret;
}

double ChordFHTW(const Graph& graph, const std::vector<std::vector<int>>& hyper_edges) {
  if (graph.m() == 0) return 0;
  std::vector<int> order = mcs::Mcs(graph);
  std::vector<int> inv_order = utils::PermInverse(order);
  std::vector<int> nbs(graph.n());
  std::vector<int> no_max(graph.n());
  for (int i=0;i<graph.n();i++){
    nbs[i] = graph.Neighbors(i).size();
  }
  double ans = 0;
  FracHyperCost fhc(graph, hyper_edges, 0, graph.n());
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
      ans = std::max(ans, fhc.Cover(utils::ToBitset(clq, graph.n())));
    }
  }
  return ans;
}
} // namespace triangulator
