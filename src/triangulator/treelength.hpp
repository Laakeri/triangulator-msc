#pragma once

#include <vector>
#include <queue>

#include "graph.hpp"

namespace triangulator {
class Treelength {
 public:
  int TreeLength(Graph graph);
 private:
  int HardSolve(const Graph& graph, int ub);
  bool Reduce(Graph graph);
  bool bad_ = false;
  int lb_ = 0;
  std::vector<Edge> fill_;
  std::queue<Graph> kernel_;
  std::vector<std::vector<int>> dists_;
};
} // namespace triangulator
