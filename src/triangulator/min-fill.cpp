#include "min-fill.hpp"

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
// THIS SHOULD BE CALLED WITH GRAPH THAT HAS BEEN BROKE INTO ATOMS
bool MinFill::GreedyElim(Graph graph) {
	UniqQue q(graph.n());
	q.Add(graph.Vertices());
	bool reduced = false;
	while (!q.Empty()) {
		int v = q.Pop();
		if (graph.IsClique(graph.Neighbors(v))) {
			reduced = true;
			auto nbs = graph.Neighbors(v);
			graph.RemoveEdgesBetween(v, nbs);
			q.Add(nbs);
		} else {
			for (auto cn : graph.CompNeighsBit(graph.adj_mat2_[v])) {
				assert(graph.IsMinsep(cn));
				auto fes = graph.FillEdges(cn);
				if (fes.size() == 1) {
					graph.AddEdges(fes);
					utils::Append(fill_, graph.MapBack(fes));
					ans_++;
					reduced = true;
					for (auto e : fes) {
						q.Add(graph.Neighbors(e.first));
						q.Add(graph.Neighbors(e.second));
					}
					break;
				}
			}
		}
	}
	if (reduced) {
		Log::Write(3, "Minfill GreedyElim success");
		Graph t_graph(graph.Edges());
		if (t_graph.m() > 0) {
			t_graph.InheritMap(graph);
			kernel_.push(t_graph);
		}
		return true;
	} else {
		return false;
	}
}

bool MinFill::Reduce(Graph graph) {
	if (graph.m() == 0) return true;
  mcs::McsMOutput minimal_triangulation = mcs::McsM(graph);
  Graph fill_graph = graph;
  fill_graph.AddEdges(minimal_triangulation.fill_edges);
  if (minimal_triangulation.fill_edges.size() <= 1) {
  	utils::Append(fill_, graph.MapBack(minimal_triangulation.fill_edges));
  	ans_ += (int)minimal_triangulation.fill_edges.size();
  	return true;
  }
  std::vector<Graph> atoms = mcs::Atoms(graph, minimal_triangulation);
  Log::Write(3, "Minfill atoms ", atoms.size());
  if (atoms.size() == 0) return true;
  if (atoms.size() > 1) {
  	for (auto& atom : atoms) {
  		atom.InheritMap(graph);
  		if (atom.n() >= 4) {
  			kernel_.push(atom);
  		}
  	}
  	return true;
  } else {
  	atoms[0].InheritMap(graph);
  	graph = atoms[0];
  }
  if (GreedyElim(graph)) {
  	return true;
  }
  kernel_.push(graph);
  return false;
}

Graph MinFill::MinimumFill(Graph graph) {
  assert(!bad_);
  bad_ = true;
  assert(kernel_.size() == 0);
  assert(ans_ == 0);
  kernel_.push(graph);
  int failed = 0;
  while (failed < (int)kernel_.size()) {
  	Log::Write(3, "Minfill reducing ", kernel_.front().n(), " ", failed, " ", (int)kernel_.size(), " ", ans_);
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
  Log::Write(3, "Minfill finished ", ans_, " ", insts.size());
  for (const Graph& inst : insts) {
  	Graph lb_graph = inst;
  	mcs::LbTriang(lb_graph);
  	int ub = inst.FillEdges(lb_graph).size();
		int ans = HardSolve(inst, ub);
		if (ans == -1) {
			ans_ += ub;
			utils::Append(fill_, inst.MapBack(inst.FillEdges(lb_graph)));
		} else {
			assert(ans >= 1);
			ans_ += ans;
		}
  }
  assert((int)fill_.size() == ans_);
  graph.AddEdges(fill_);
  return graph;
}

int MinFill::HardSolve(const Graph& graph, int ub) {
	Log::Write(3, "Minfill HARDSOLVE ", graph.n(), " ", graph.m(), " ", ub);
	Enumerator pe;
	pe.SetMFFilter(ub-1);
	auto pmcs = pe.AdvancedPmcEnum(graph);
	Log::Write(3, "PMCS ", pmcs.size());

	Timer bt_time;
	bt_time.start();
	BtAlgorithm bt;
	bt.SetMF();
	auto sol = bt.Solve(graph, pmcs);
	Log::Write(3, "bt time ", bt_time.getTime().count());
	Log::Write(3, "bt ", sol.first, " ", sol.second.size());

	if (sol.first == -1 || sol.first >= ub) {
		return -1;
	}
	assert(sol.first == (int)sol.second.size());
	utils::Append(fill_, graph.MapBack(sol.second));
	return sol.first;
}
} // namespace triangulator
