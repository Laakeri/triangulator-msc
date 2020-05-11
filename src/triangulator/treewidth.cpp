#include "treewidth.hpp"

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
void TW::UpdLB(const Graph& graph) {
	lb_ = std::max(lb_, graph.Degeneracy());
}

// THIS SHOULD BE CALLED WITH GRAPH THAT HAS BEEN BROKE INTO ATOMS
bool TW::GreedyElim(Graph graph) {
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
			lb_ = std::max(lb_, (int)nbs.size());
		} else {
			for (auto cn : graph.CompNeighsBit(graph.adj_mat2_[v])) {
				assert(graph.IsMinsep(cn));
				if (graph.IsAlmostClique(cn.Elements())) {
					auto fes = graph.FillEdges(cn);
					graph.AddEdges(fes);
					utils::Append(fill_, graph.MapBack(fes));
					lb_ = std::max(lb_, cn.Popcount());
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
		Log::Write(3, "TWPP GreedyElim success");
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

bool TW::AlmostClique(Graph graph) {
	assert(graph.IsConnected());
	std::vector<int> found;
	int maxc = graph.n() + 1;
	for (int i=0;i<graph.n();i++) {
		Graph t_graph = graph;
		auto nbs = t_graph.Neighbors(i);
		assert((int)nbs.size() > 2);
		t_graph.RemoveEdgesBetween(i, nbs);
		mcs::McsMOutput minimal_triangulation = mcs::McsM(t_graph);
		t_graph.AddEdges(minimal_triangulation.fill_edges);
  	Bitset not_rm(graph.n());
  	not_rm.FillTrue();
		for (int it = 0; it < graph.n(); it++) {
	    int x = minimal_triangulation.elimination_order[it];
	    if (minimal_triangulation.is_maximal_clique_point[it]) {
	      Bitset cand_clique_b = t_graph.adj_mat2_[x] & not_rm;
	      cand_clique_b.SetFalse(x);
	      if (cand_clique_b.Popcount()>0 && graph.IsClique(cand_clique_b)) {
	      	auto cand_clique = cand_clique_b.Elements();
	      	Log::Write(3, "cc ", cand_clique);
	      	assert(graph.IsClique(cand_clique));
	      	assert(t_graph.IsMinsep(cand_clique));
	      	assert(!graph.IsMinsep(cand_clique));
	      	cand_clique.push_back(i);
		  		utils::SortAndDedup(cand_clique);
		  		assert(graph.IsMinsep(cand_clique));
		  		assert(!graph.IsClique(cand_clique));
		  		assert(graph.IsAlmostClique(cand_clique));
		  		auto comps = graph.NComponents(cand_clique);
		  		lb_ = std::max(lb_, (int)cand_clique.size());
		  		int tmaxc = 0;
		  		for (const auto& comp : graph.NComponents(cand_clique)) {
		  			tmaxc = std::max(tmaxc, (int)comp.size());
		  		}
		  		Log::Write(3, "Almostclique ", cand_clique.size(), " ", graph.n(), " ", tmaxc);
		  		if (tmaxc < maxc) {
		  			maxc = tmaxc;
		  			found = cand_clique;
		  		}
	      }
	    }
	    not_rm.SetFalse(x);
	  }
	}
	if (!found.empty()) {
		assert(maxc < graph.n());
		assert(graph.IsMinsep(found));
		assert(!graph.IsClique(found));
		assert(graph.IsAlmostClique(found));
		auto fes = graph.FillEdges(found);
		assert(fes.size()>0);
		graph.AddEdges(fes);
		utils::Append(fill_, graph.MapBack(fes));
		kernel_.push(graph);
		return true;
	} else {
		assert(maxc == graph.n() + 1);
		return false;
	}
}

bool TW::Reduce(Graph graph) {
	if (graph.m() == 0) return true;
  UpdLB(graph);
  mcs::McsMOutput minimal_triangulation = mcs::McsM(graph);
  Graph fill_graph = graph;
  fill_graph.AddEdges(minimal_triangulation.fill_edges);
  if (minimal_triangulation.fill_edges.size() <= 1 || mcs::Treewidth(fill_graph) <= lb_) {
  	utils::Append(fill_, graph.MapBack(minimal_triangulation.fill_edges));
  	lb_ = std::max(lb_, mcs::Treewidth(fill_graph));
  	return true;
  }
  std::vector<Graph> atoms = mcs::Atoms(graph, minimal_triangulation);
  Log::Write(3, "TWPP atoms ", atoms.size());
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
  Graph lb_graph = graph;
 	mcs::LbTriang(lb_graph);
 	int ub = mcs::Treewidth(lb_graph);
  if (ub <= lb_) {
  	utils::Append(fill_, graph.MapBack(graph.FillEdges(lb_graph)));
  	return true;
  }
  if (AlmostClique(graph)) {
  	return true;
  }
  kernel_.push(graph);
  return false;
}

Graph TW::Treewidth(Graph graph) {
  assert(!bad_);
  bad_ = true;
  assert(kernel_.size() == 0);
  assert(lb_ == 0);
  kernel_.push(graph);
  int failed = 0;
  while (failed < (int)kernel_.size()) {
  	Log::Write(3, "TWPP reducing ", kernel_.front().n(), " ", failed, " ", (int)kernel_.size(), " ", lb_);
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
  Log::Write(3, "TWPP finished ", lb_, " ", insts.size());
  for (const Graph& inst : insts) {
  	Graph lb_graph = inst;
  	mcs::LbTriang(lb_graph);
  	int ub = mcs::Treewidth(lb_graph);
  	if (ub <= lb_) {
  		utils::Append(fill_, inst.MapBack(inst.FillEdges(lb_graph)));
  		continue;
  	}
		int ans = HardSolve(inst, ub);
		if (ans == -1) {
			lb_ = ub;
			utils::Append(fill_, inst.MapBack(inst.FillEdges(lb_graph)));
		} else {
			assert(ans > 2); // Treewidth <= 2 will be decided by preprocessing
			lb_ = std::max(lb_, ans);
		}
  }
  graph.AddEdges(fill_);
  assert(mcs::Treewidth(graph) == lb_);
  return graph;
}

int TW::HardSolve(const Graph& graph, int ub) {
	Log::Write(3, "TW HARDSOLVE ", graph.n(), " ", graph.m(), " ", lb_, " ", ub);


	Enumerator pe;
	pe.SetTWFilter(ub-1);
	auto pmcs = pe.AdvancedPmcEnum(graph);
	Log::Write(3, "PMCS ", pmcs.size());

	Timer bt_time;
	bt_time.start();
	BtAlgorithm bt;
	auto sol = bt.Solve(graph, pmcs);
	Log::Write(3, "bt time ", bt_time.getTime().count());
	Log::Write(3, "bt ", sol.first, " ", sol.second.size());

	assert(sol.first < ub);
	if (sol.first == -1) {
		return -1;
	}
	utils::Append(fill_, graph.MapBack(sol.second));
	return sol.first;
}
} // namespace triangulator
