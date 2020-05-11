#pragma once

#include <vector>

#include "graph.hpp"
#include "maxhs_iface.hpp"

namespace triangulator {
class BtAlgorithm {
public:
  enum class Problem {
  	tw,
  	mf,
  	tts,
    id,
    phylogeny
  };
  std::pair<int64_t, std::vector<Edge>> Solve(const Graph& graph, const std::vector<Bitset>& pmcs);
  void BuildColorMaxSAT(const Graph& graph, const std::vector<Bitset>& pmcs, const std::vector<int>& colors, const std::vector<int>& color_vars, MaxhsInterface& maxsat);
  void SetTW();
  void SetMF();
  void SetTTS(const std::vector<int64_t>& bn_domains);
  void SetId();
  void SetPhylogeny(const std::vector<int>& phyl_colors);
private:
	Problem problem_ = Problem::tw;
  std::vector<int64_t> bn_domains_;
  std::vector<int> phyl_colors_;
  int64_t MergeCost(int64_t c1, int64_t c2) const;
  int64_t CliqueCost(const Graph& graph, const Bitset& pmc, const Bitset& parent_sep, int64_t id) const;
};
} // namespace triangulator