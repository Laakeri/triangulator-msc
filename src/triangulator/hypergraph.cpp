#include "hypergraph.hpp"

#include <vector>
#include <algorithm>
#include <iostream>

#include "utils.hpp"

namespace triangulator {
namespace {
std::vector<std::pair<int, int>> PrimalEdges(const std::vector<std::vector<int>>& edges) {
  std::vector<std::pair<int, int>> es;
  for (auto& e : edges) {
    for (int i = 0; i < (int)e.size(); i++) {
      for (int ii = i + 1; ii < (int)e.size(); ii++) {
        es.push_back({e[i], e[ii]});
      }
    }
  }
  return es;
}
} // namespace

HyperGraph::HyperGraph(int n) : primal_(n) { }

HyperGraph::HyperGraph(std::vector<std::vector<int>> edges) : primal_(PrimalEdges(edges)), edges_(edges) {
  for (auto& e : edges_) {
    for (int& v : e) {
      v = primal_.MapInto(v);
    }
    utils::SortAndDedup(e);
  }
}
Graph HyperGraph::PrimalGraph() const {
  return primal_;
}
const std::vector<std::vector<int> >& HyperGraph::Edges() const {
  return edges_;
}
void HyperGraph::AddEdge(std::vector<int> edge) {
  utils::SortAndDedup(edge);
  for (int i = 0; i < (int)edge.size(); i++) {
    for (int ii = i + 1; ii < (int)edge.size(); ii++) {
      primal_.AddEdge(edge[i], edge[ii]);
    }
  }
  edges_.push_back(edge);
}
int HyperGraph::n() const {
  return primal_.n();
}
int HyperGraph::m() const {
  return edges_.size();
}
} // namespace triangulator
