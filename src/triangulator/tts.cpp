#include "tts.hpp"

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
bool TTS::Reduce(const Graph& graph, const std::vector<int64_t>& domains) {
  if (graph.m() == 0) return true;
  mcs::McsMOutput minimal_triangulation = mcs::McsM(graph);
  if (minimal_triangulation.fill_edges.size() == 0) {
	  Graph fill_graph = graph;
	  fill_graph.AddEdges(minimal_triangulation.fill_edges);
	  ans_ += mcs::TotalTableSize(fill_graph, domains);
  	utils::Append(fill_, graph.MapBack(minimal_triangulation.fill_edges));
  	return true;
  }
  std::vector<Graph> atoms = mcs::Atoms(graph, minimal_triangulation);
  Log::Write(3, "TTSPP atoms ", atoms.size());
  if (atoms.size() == 0) return true;
	for (auto& atom : atoms) {
		atom.InheritMap(graph);
		kernel_.push(atom);
	}
	return atoms.size() > 1 || atoms[0].n() < graph.n();
}

Graph TTS::TotalTableSize(Graph graph, std::vector<int64_t> domains) {
  assert(!bad_);
  bad_ = true;
  assert(kernel_.size() == 0);
  assert(ans_ == 0);
  kernel_.push(graph);
  int failed = 0;
  while (failed < (int)kernel_.size()) {
  	Log::Write(3, "TTSPP reducing ", kernel_.front().n(), " ", failed, " ", (int)kernel_.size());
  	if (Reduce(kernel_.front(), domains)) {
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
  Log::Write(3, "TTSPP finished ", insts.size());
  for (const Graph& inst : insts) {
  	Graph lb_graph = inst;
  	mcs::LbTriang(lb_graph);
  	int64_t ub = mcs::TotalTableSize(lb_graph, domains);
		int64_t cost = HardSolve(inst, domains, ub);
		if (cost == -1) {
			ans_ += ub;
			utils::Append(fill_, inst.MapBack(inst.FillEdges(lb_graph)));
		} else {
			ans_ += cost;
		}
  }
  graph.AddEdges(fill_);
  Log::Write(3, "lol ", ans_);
  assert(mcs::TotalTableSize(graph, domains) == ans_);
  return graph;
}

int64_t TTS::HardSolve(const Graph& graph, const std::vector<int64_t>& domains, int64_t ub) {
	Log::Write(3, "TTS HARDSOLVE ", graph.n(), " ", graph.m(), " ", ub);
	Enumerator pe;
	pe.SetTTSFilter(ub-1, domains);
	auto pmcs = pe.AdvancedPmcEnum(graph);
	Log::Write(3, "PMCS ", pmcs.size());

	Timer bt_time;
	bt_time.start();
	BtAlgorithm bt;
	bt.SetTTS(domains);
	auto sol = bt.Solve(graph, pmcs);
	Log::Write(3, "bt time ", bt_time.getTime().count());
	Log::Write(3, "bt ", sol.first, " ", sol.second.size());

	assert(sol.first <= ub);
	if (sol.first == -1) {
		return -1;
	}
	utils::Append(fill_, graph.MapBack(sol.second));
	return sol.first;
}
} // namespace triangulator
