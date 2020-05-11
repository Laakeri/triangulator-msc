#pragma once

#include <vector>
#include <queue>

#include "graph.hpp"

namespace triangulator {
class GHTW {
 public:
  Graph CompGHTW(Graph graph, const std::vector<std::vector<int>>& hyperedges);
 private:
  int HardSolve(const Graph& graph, const std::vector<std::vector<int>>& hyperedges, int ub);
  bool Reduce(Graph graph, const std::vector<std::vector<int>>& hyperedges);
  std::vector<Edge> fill_;
  std::queue<Graph> kernel_;
  bool bad_ = false;
  int lb_ = 0;
};

int ChordGHTW(const Graph& graph, const std::vector<std::vector<int>>& hyper_edges);
} // namespace triangulator
