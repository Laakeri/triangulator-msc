#pragma once

#include <vector>

#include "graph.hpp"

namespace triangulator {
namespace mcs {
struct McsMOutput {
  std::vector<Edge> fill_edges;
  std::vector<int> elimination_order;
  std::vector<char> is_maximal_clique_point;
};

// Computes an elimination order that is perfect iff the graph is chordal in O(n + m) [1]
std::vector<int> Mcs(const Graph& graph);

// Computes minimal triangulation of a graph in O(n m) using the MCS-M algorithm introduced in [1].
// The return value contains the fill edges, the elimination order and additional information for computing the atoms of the graph.
McsMOutput McsM(const Graph& graph);

void LbTriang(Graph& graph);

void FastMT(Graph& graph);

// Returns the atoms of a graph in O(n m).
// Takes the return value of MCS-M as input.
// The atoms are defined and the algorithm is given in [2]
std::vector<Graph> Atoms(const Graph& graph, const McsMOutput& mcs_m_output);

// Returns the treewidth of a chordal graph in O(n + m)
int Treewidth(const Graph& graph);

int64_t TotalTableSize(const Graph& graph, const std::vector<int64_t>& domains);

int Treelength(const Graph& graph, const std::vector<std::vector<int>>& dists);
// [1] Anne Berry, Jean R. S. Blair, Pinar Heggernes, and Barry W. Peyton. Maximum Cardinality Search for Computing Minimal Triangulations of Graphs. Algorithmica 39, (2004), 287-298
// [2] Robert E. Tarjan. Decomposition by clique separators. Discrete Mathematics 55, (1985), 221-232

} // namespace mcs
} // namespace triangulator
