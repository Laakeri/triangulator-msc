#pragma once

#include <vector>
#include <queue>

#include "graph.hpp"

namespace triangulator {
class MinFill {
 public:
  Graph MinimumFill(Graph graph);
 private:
  int HardSolve(const Graph& graph, int ub);
  void UpdLB(const Graph& graph);
  bool GreedyElim(Graph graph);
  bool Reduce(Graph graph);
  bool bad_ = false;
  int ans_ = 0;
  std::vector<Edge> fill_;
  std::queue<Graph> kernel_;
};
} // namespace triangulator
