#pragma once

#include <vector>
#include <queue>

#include "graph.hpp"

namespace triangulator {
class FHTW {
 public:
  Graph CompFHTW(Graph graph, const std::vector<std::vector<int>>& hyperedges);
 private:
  double HardSolve(const Graph& graph, const std::vector<std::vector<int>>& hyperedges, double ub);
  bool Reduce(Graph graph, const std::vector<std::vector<int>>& hyperedges);
  std::vector<Edge> fill_;
  std::queue<Graph> kernel_;
  bool bad_ = false;
  double lb_ = 0;
};

double ChordFHTW(const Graph& graph, const std::vector<std::vector<int>>& hyper_edges);
} // namespace triangulator
