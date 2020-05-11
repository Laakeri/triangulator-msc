#pragma once

#include <vector>
#include <tuple>

#include "graph.hpp"
#include "bitset.hpp"

namespace triangulator {
class Enumerator {
 private:
 	enum class FilterType {
 		tw,
 		mf,
 		tts,
 		ghtw,
    perf_phy,
    tl
 	};
 	bool GoodMinsep(const Bitset& minsep, const Graph& graph) const;
 	bool GoodPMC(const Bitset& pmc, const Graph& graph) const;
 	FilterType filter_ = FilterType::tw;
 	int twub_ = (int)1e9;
 	int mfub_ = (int)1e9;
 	int64_t ttsub_ = (int64_t)1e18;
 	int ghtwub_ = (int)1e9;
  int tlub_ = (int)1e9;

 	uint64_t crosses_ = 0;
	uint64_t dfss_ = 0;
	uint64_t gos_ = 0;
	std::vector<std::tuple<int, int>> impl_stack_;

	std::vector<int64_t> bn_domains_;
  std::vector<int> phyl_colors_;
  std::vector<std::vector<int>> tl_dists_;
 
	void Gogo(int ind, int height, const Graph& graph, const std::vector<std::tuple<Bitset, Bitset>>& s_index, const std::vector<std::vector<int>>& childs, const Bitset& t_sep, Bitset musthave, const std::vector<Bitset>& nohave, std::vector<Bitset>& con, std::vector<Bitset>& tpmcs);
 	void Combine(const Graph& graph, int x, 
             const std::vector<Bitset>& s_minseps, 
             const std::vector<Bitset>& t_minseps, 
             std::vector<Bitset>& tpmcs,
             Timer& pre_timer,
             Timer& rec_timer);
 public:
 	void SetTWFilter(int tw);
 	void SetMFFilter(int mf);
 	void SetTTSFilter(int64_t tts, std::vector<int64_t> bn_domains);
 	void SetGHTWFilter(int ghtw);
  void SetPerfPhyFilter(std::vector<int> phyl_colors);
  void SetTLFilter(int tl, std::vector<std::vector<int>> dists);
	std::vector<Bitset> AdvancedPmcEnum(const Graph& graph);
};
} // namespace triangulator
