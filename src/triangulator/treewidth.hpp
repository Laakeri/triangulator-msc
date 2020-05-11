#pragma once

#include <vector>
#include <queue>

#include "graph.hpp"

namespace triangulator {
class TW {
 public:
  Graph Treewidth(Graph graph);
 private:
  int HardSolve(const Graph& graph, int ub);
  void UpdLB(const Graph& graph);
  bool GreedyElim(Graph graph);
  bool Reduce(Graph graph);
  bool AlmostClique(Graph graph);
  bool bad_ = false;
  std::vector<Edge> fill_;
  int lb_ = 0;
  std::queue<Graph> kernel_;
};
} // namespace triangulator
