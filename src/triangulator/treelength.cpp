#include "treelength.hpp"

#include <vector>
#include <cassert>
#include <map>
#include <iostream>

#include "enumerator.hpp"
#include "graph.hpp"
#include "mcs.hpp"
#include "utils.hpp"
#include "bt_algorithm.hpp"

namespace triangulator {
bool Treelength::Reduce(Graph graph) {
	if (graph.m() == 0) return true;
  mcs::McsMOutput minimal_triangulation = mcs::McsM(graph);
  Graph fill_graph = graph;
  fill_graph.AddEdges(minimal_triangulation.fill_edges);
  int treelength = mcs::Treelength(fill_graph, dists_);
  if (minimal_triangulation.fill_edges.size() == 0 || treelength <= lb_) {
  	lb_ = std::max(lb_, treelength);
    utils::Append(fill_, graph.MapBack(minimal_triangulation.fill_edges));
    return true;
  }
  lb_ = std::max(lb_, 2);
  std::vector<Graph> atoms = mcs::Atoms(graph, minimal_triangulation);
  Log::Write(3, "Treelength atoms ", atoms.size());
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
  int ub = mcs::Treelength(lb_graph, dists_);
  if (ub <= lb_) {
    utils::Append(fill_, graph.MapBack(graph.FillEdges(lb_graph)));
    return true;
  }
  kernel_.push(graph);
  return false;
}

int Treelength::TreeLength(Graph graph) {
  assert(!bad_);
  bad_ = true;
  assert(kernel_.size() == 0);
  assert(lb_ == 0);
  kernel_.push(graph);
  int failed = 0;
  dists_ = graph.DistanceMatrix();
  for (int i=0;i<graph.n();i++){
    for (int ii=0;ii<graph.n();ii++){
      assert(dists_[i][ii] == dists_[ii][i]);
      if (i != ii) {
        assert(dists_[i][ii] >= 1);
        assert(dists_[i][ii] < graph.n() || dists_[i][ii] == graph.n()+1);
      } else {
        assert(dists_[i][ii] == 0);
      }
    }
  }
  while (failed < (int)kernel_.size()) {
  	Log::Write(3, "Treelength reducing ", kernel_.front().n(), " ", failed, " ", (int)kernel_.size(), " ", lb_);
  	if (Reduce(kernel_.front())) {
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
  Log::Write(3, "Treelength finished ", lb_, " ", insts.size());
  for (const Graph& inst : insts) {
  	Graph lb_graph = inst;
  	mcs::LbTriang(lb_graph);
    int ub = mcs::Treelength(lb_graph, dists_);
		int ans = HardSolve(inst, ub);
		if (ans == -1) {
      lb_ = std::max(lb_, ub);
			utils::Append(fill_, inst.MapBack(inst.FillEdges(lb_graph)));
		} else {
      lb_ = std::max(lb_, ans);
		}
  }
  graph.AddEdges(fill_);
  assert(mcs::Treelength(graph, dists_) == lb_);
  return lb_;
}

int Treelength::HardSolve(const Graph& graph, int ub) {
	assert(ub < graph.n());
  Log::Write(3, "Treelength HARDSOLVE ", graph.n(), " ", graph.m(), " ", ub);
	Enumerator pe;
	pe.SetTLFilter(ub-1, dists_);
	auto pmcs = pe.AdvancedPmcEnum(graph);
	Log::Write(3, "PMCS ", pmcs.size());

  Timer dt_time;
  dt_time.start();
  std::vector<std::vector<Bitset>> ord_pmcs(ub);
  for (int i=0;i<(int)pmcs.size();i++){
    int ttl = utils::MaxInSub(dists_, graph.MapBack(pmcs[i].Elements()));
    assert(ttl >= 1 && ttl < graph.n());
    assert(ttl < ub);
    if (ttl < ub) {
      ord_pmcs[ttl].push_back(pmcs[i]);
    }
  }
  pmcs.clear();
  std::vector<int> trs(ub);
  for (int i=1;i<ub;i++) {
    for (const auto& pmc : ord_pmcs[i]) {
      pmcs.push_back(pmc);
    }
    ord_pmcs[i].clear();
    trs[i] = pmcs.size();
  }
  Log::Write(3, "dt time ", dt_time.getTime().count());

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
} // namespace triangulator
