#pragma once

#include <vector>
#include <queue>

#include "graph.hpp"

namespace triangulator {
class TTS {
 public:
  Graph TotalTableSize(Graph graph, std::vector<int64_t> domains);
 private:
  int64_t HardSolve(const Graph& graph, const std::vector<int64_t>& domains, int64_t ub);
  bool Reduce(const Graph& graph, const std::vector<int64_t>& domains);
  std::vector<Edge> fill_;
  std::queue<Graph> kernel_;
  bool bad_ = false;
  int64_t ans_ = 0;
};
} // namespace triangulator
