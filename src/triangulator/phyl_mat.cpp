#include "phyl_mat.hpp"

#include <cassert>
#include <map>
#include <iostream>

#include "matrix.hpp"
#include "graph.hpp"
#include "staticset.hpp"

namespace triangulator {
PhylMat::PhylMat(Matrix<int> mat) : mat_(mat), arity_(0), has_missing_(false), graph_(0), color_(0) {
  for (int ii = 0; ii < (int)mat_.Columns(); ii++) {
    std::map<int, int> cnts;
    for (int i = 0; i < (int)mat_.Rows(); i++) {
      assert(mat_[i][ii] >= -1);
      if (mat_[i][ii] >= 0) {
        cnts[mat_[i][ii]]++;
      }
    }
    int goods = 0;
    for (auto c : cnts) {
      if (c.second >= 2) {
        goods++;
      }
    }
    for (int i = 0; i < (int)mat_.Rows(); i++) {
      if (mat_[i][ii] >= 0 && (goods < 2 || cnts[mat_[i][ii]] < 2)) {
        mat_[i][ii] = -1;
      }
    }
  }
  for (int i = 0; i < (int)mat_.Rows(); i++) {
    for (int ii = 0; ii < (int)mat_.Columns(); ii++) {
      if (mat_[i][ii] == -1) {
        has_missing_ = true;
      } else {
        arity_ = std::max(arity_, mat_[i][ii] + 1);
      }
    }
  }
  std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>> edges;
  for (int i = 0; i < (int)mat_.Rows(); i++) {
    for (int ii = 0; ii < (int)mat_.Columns(); ii++) {
      for (int j = ii + 1; j < (int)mat_.Columns(); j++) {
        if (mat_[i][ii] >= 0 && mat_[i][j] >= 0) {
          edges.push_back({{ii, mat_[i][ii]}, {j, mat_[i][j]}});
        }
      }
    }
  }
  StaticSet<std::pair<int, int>> ss(edges);
  graph_ = Graph(ss.Size());
  for (const auto& e : edges) {
    graph_.AddEdge(ss.Rank(e.first), ss.Rank(e.second));
  }
  color_.resize(graph_.n());
  for (int i = 0; i < graph_.n(); i++) {
    color_[i] = ss.Kth(i).first;
  }
  assert(arity_ >= 0);
}
int PhylMat::Arity() const {
  return arity_;
}
bool PhylMat::HasMissing() const {
  return has_missing_;
}
const Graph& PhylMat::GetGraph() const {
  return graph_;
}
int PhylMat::Color(int i) const {
  return color_[i];
}
int PhylMat::n() const {
  return mat_.Rows();
}
int PhylMat::m() const {
  return mat_.Columns();
}
int PhylMat::GetMat(int i, int j) const {
  return mat_[i][j];
}
} // namespace triangulator