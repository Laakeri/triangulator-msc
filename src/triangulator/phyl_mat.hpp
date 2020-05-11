#pragma once

#include "matrix.hpp"
#include "graph.hpp"

namespace triangulator {
class PhylMat {
 public:
  explicit PhylMat(Matrix<int> mat);
  int Arity() const;
  bool HasMissing() const;
  const Graph& GetGraph() const;
  int Color(int i) const;
  int GetMat(int i, int j) const;
  int n() const;
  int m() const;
 private:
  Matrix<int> mat_;
  int arity_;
  bool has_missing_;
  Graph graph_;
  std::vector<int> color_;
};

} // namespace triangulator