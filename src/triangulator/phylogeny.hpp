#pragma once

#include <vector>

#include "graph.hpp"
#include "phyl_mat.hpp"

namespace triangulator {
bool PerfectPhylogeny(const PhylMat& pm);
std::vector<int> BinaryCharRemove(const PhylMat& pm);
std::vector<int> CharRemoveMaxSat(const PhylMat& pm);
} // namespace