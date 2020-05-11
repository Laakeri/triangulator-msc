#pragma once

#include <vector>
#include <string>
#include <istream>
#include <map>

#include "graph.hpp"
#include "hypergraph.hpp"
#include "utils.hpp"
#include "phyl_mat.hpp"

namespace triangulator {

class Io {
public:
  Graph ReadGraph(std::istream& in);
  HyperGraph ReadHyperGraph(std::istream& in);
  PhylMat ReadPhylMat(std::istream& in);
  std::pair<std::vector<int64_t>, Graph> ReadBayesNet(std::istream& in);
  std::string MapBack(int v);
private:
  bool dimacs_;
  StaticSet<std::string> vertex_map_;
};
} // namespace triangulator
