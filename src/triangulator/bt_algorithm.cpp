#include "bt_algorithm.hpp"

#include <vector>
#include <map>
#include <cassert>
#include <queue>

#include "graph.hpp"
#include "bitset.hpp"
#include "utils.hpp"
#include "maxhs_iface.hpp"

using std::pair;
using std::vector;
using std::tuple;
using std::make_tuple;
using std::get;

namespace triangulator {
namespace {
struct Triplets {
  const vector<pair<int, int>> child_missing = {{-1, -1}};
  vector<tuple<int, int, int>> trips;
  const Graph& graph;
  const vector<Bitset>& pmcs;
  vector<Bitset> separators;
  vector<Bitset> components;
  pair<int, int> root;
  Triplets(const Graph& graph_, const vector<Bitset>& pmcs_) : graph(graph_), pmcs(pmcs_) {
    assert(graph.IsConnected());
    Bitset all_vertices(graph.n());
    all_vertices.FillUpTo(graph.n());
    Bitset no_vertices(graph.n());
    separators.push_back(no_vertices);
    components.push_back(all_vertices);
    vector<tuple<int, int, int>> triplets;
    for (const auto& pmc : pmcs) {
      for (const Bitset& sep : graph.CompNeighsBit(pmc)) {
        separators.push_back(sep);
      }
    }
    utils::SortAndDedup(separators);
    for (const auto& sep : separators) {
      all_vertices.TurnOff(sep);
      for (const auto& comp : graph.BitComps(all_vertices)) {
        components.push_back(comp);
      }
      all_vertices |= sep;
    }
    utils::SortAndDedup(components);
    root = {0, (int)components.size() - 1};
    assert(separators[root.first] == no_vertices);
    assert(components[root.second] == all_vertices);
    for (int i=0;i<(int)pmcs.size();i++) {
      const auto& pmc = pmcs[i];
      for (const Bitset& sep : graph.CompNeighsBit(pmc)) {
        assert(sep.Popcount() < pmc.Popcount());
        Bitset vis = all_vertices;
        vis.TurnOff(sep);
        Bitset comp = pmc;
        comp.TurnOff(sep);
        graph.Dfs2Bit(vis, comp);
        comp.TurnOff(sep);
        int sep_id = utils::BSFind(separators, sep);
        int comp_id = utils::BSFind(components, comp);
        assert(sep_id >= 0 && comp_id >= 0);
        trips.push_back(make_tuple(i, sep_id, comp_id));
      }
      trips.push_back(make_tuple(i, root.first, root.second));
    }
    auto cmp = [&](tuple<int, int, int> a, tuple<int, int, int> b) {
      int size_a = separators[get<1>(a)].Popcount() + components[get<2>(a)].Popcount();
      int size_b = separators[get<1>(b)].Popcount() + components[get<2>(b)].Popcount();
      return size_a < size_b;
    };
    std::sort(trips.begin(), trips.end(), cmp);
  }
  vector<pair<int, int>> Children(tuple<int, int, int> trip) {
    assert(get<0>(trip) >= 0 && get<1>(trip) >= 0 && get<2>(trip) >= 0);
    assert(get<0>(trip) < (int)pmcs.size() && get<1>(trip) < (int)separators.size() && get<2>(trip) < (int)components.size());
    Bitset c_comps = components[get<2>(trip)];
    c_comps.TurnOff(pmcs[get<0>(trip)]);
    vector<pair<int, int>> ret;
    for (const Bitset& c_comp : graph.BitComps(c_comps)) {
      Bitset c_sep = graph.Neighbors(c_comp);
      int sep_id = utils::BSFind(separators, c_sep);
      int comp_id = utils::BSFind(components, c_comp);
      if (sep_id < 0 || comp_id < 0) {
        return child_missing;
      }
      assert(sep_id < (int)separators.size() && comp_id < (int)components.size());
      ret.push_back({sep_id, comp_id});
    }
    return ret;
  }
  const Bitset& GetPmc(tuple<int, int, int> trip) {
    return pmcs[get<0>(trip)];
  }
  const Bitset& GetSep(tuple<int, int, int> trip) {
    return separators[get<1>(trip)];
  }
  const Bitset& GetComp(tuple<int, int, int> trip) {
    return components[get<2>(trip)];
  }
};
} // namespace

int64_t BtAlgorithm::MergeCost(int64_t c1, int64_t c2) const {
  if (problem_ == Problem::tw) {
    return std::max(c1, c2);
  } else if (problem_ == Problem::tts) {
    return c1+c2;
  } else if (problem_ == Problem::mf) {
    return c1+c2;
  } else if (problem_ == Problem::id) {
  	return std::max(c1, c2);
  } else if (problem_ == Problem::phylogeny) {
    return c1+c2;
  } else {
    assert(0);
  }
}

int64_t BtAlgorithm::CliqueCost(const Graph& graph, const Bitset& pmc, const Bitset& parent_sep, int64_t id) const {
  if (problem_ == Problem::tw) {
    return pmc.Popcount() - 1;
  } else if (problem_ == Problem::tts) {
    int64_t cost = 1;
    for (int v : pmc) {
      cost *= bn_domains_[graph.MapBack(v)];
    }
    return cost;
  } else if (problem_ == Problem::mf) {
    return graph.FillSize(pmc) - graph.FillSize(parent_sep);
  } else if (problem_ == Problem::id) {
  	return id;
  } else if (problem_ == Problem::phylogeny) {
    std::vector<std::pair<int, int>> cols;
    cols.reserve(pmc.Popcount());
    assert((int)phyl_colors_.size() == graph.n());
    for (int v : pmc) {
      cols.push_back({phyl_colors_[v], parent_sep.Get(v)});
    }
    std::sort(cols.begin(), cols.end());
    int cost = 0;
    for (int i = 1; i < (int)cols.size(); i++) {
      if (cols[i].first == cols[i-1].first && (!cols[i].second || !cols[i-1].second)) {
        cost++;
      }
    }
    return cost;
  } else {
    assert(0);
  }
}

pair<int64_t, vector<Edge>> BtAlgorithm::Solve(const Graph& graph, const vector<Bitset>& pmcs) {
  assert(graph.IsConnected());
  Triplets triplets(graph, pmcs);
  std::map<pair<int, int>, int64_t> dp;
  std::map<pair<int, int>, int> opt_choice;
  for (const auto& triplet : triplets.trips) {
    int64_t cost  = CliqueCost(graph, triplets.GetPmc(triplet), triplets.GetSep(triplet), get<0>(triplet));
    bool child_missing = false;
    for (const pair<int, int>& child : triplets.Children(triplet)) {
      if (child.first == -1) {
        child_missing = true;
        break;
      } else if (dp.count(child)) {
        cost = MergeCost(cost, dp[child]);
      } else {
        child_missing = true;
        break;
      }
    }
    if (child_missing) continue;
    pair<int, int> this_state = {get<1>(triplet), get<2>(triplet)};
    if (!dp.count(this_state) || cost < dp[this_state]) {
      dp[this_state] = cost;
      opt_choice[this_state] = get<0>(triplet);
    }
  }
  if (!dp.count(triplets.root)) return {-1, {}};
  else {
    vector<Edge> fill_edges;
    std::queue<pair<int, int>> reconstruct;
    reconstruct.push(triplets.root);
    while (!reconstruct.empty()) {
      pair<int, int> state = reconstruct.front();
      reconstruct.pop();
      assert(dp.count(state));

      tuple<int, int, int> triplet = make_tuple(opt_choice[state], state.first, state.second);

      for (Edge fe : graph.FillEdges(triplets.GetPmc(triplet))) {
        if (triplets.GetSep(triplet).Get(fe.first) && triplets.GetSep(triplet).Get(fe.second)) continue;
        fill_edges.push_back(fe);
      }
      for (const pair<int, int>& child : triplets.Children(make_tuple(opt_choice[state], state.first, state.second))) {
        assert(child.first >= 0 && child.second >= 0);
        assert(dp.count(child));
        reconstruct.push(child);
      }
    }
    return {dp[triplets.root], fill_edges};
  }
}

void BtAlgorithm::BuildColorMaxSAT(const Graph& graph, const vector<Bitset>& pmcs, const vector<int>& colors, const vector<int>& color_vars, MaxhsInterface& maxsat) {
  assert(graph.IsConnected());
  assert((int)colors.size() == graph.n());
  for (int c : colors) {
    assert(c >= 0 && c < (int)color_vars.size() && color_vars[c] > 0);
  }
  Triplets triplets(graph, pmcs);

  vector<int> pmc_vars;
  for (const auto& pmc : pmcs) {
    std::map<int, int> color_cs;
    for (int v : pmc) {
      color_cs[colors[v]]++;
    }
    vector<int> bad_colors;
    for (auto c : color_cs) {
      if (c.second >= 2) {
        bad_colors.push_back(c.first);
      }
    }
    if (bad_colors.size() == 0) {
      pmc_vars.push_back(maxsat.TrueLit);
    } else {
      pmc_vars.push_back(maxsat.NewVar());
      for (int c : bad_colors) {
        maxsat.AddClause({-pmc_vars.back(), color_vars[c]});
      }
    }
  }

  std::map<pair<int, int>, int> block_vars;
  std::map<pair<int, int>, vector<int>> child_trips;
  for (const auto& triplet : triplets.trips) {
    vector<int> cbs;
    bool child_missing = false;
    for (const pair<int, int>& child : triplets.Children(triplet)) {
      if (block_vars.count(child)) {
        cbs.push_back(block_vars[child]);
      } else {
        child_missing = true;
        break;
      }
    }
    if (child_missing) continue;
    int trip_var = maxsat.NewVar();
    for (int cb : cbs) {
      maxsat.AddClause({-trip_var, cb});
    }
    maxsat.AddClause({-trip_var, pmc_vars[get<0>(triplet)]});

    pair<int, int> this_state = {get<1>(triplet), get<2>(triplet)};
    if (!block_vars.count(this_state)) {
      block_vars[this_state] = maxsat.NewVar();
    }
    child_trips[this_state].push_back(trip_var);
  }
  for (const auto& bv : block_vars) {
    auto this_state = bv.first;
    int var = bv.second;
    vector<int> clause = child_trips[this_state];
    clause.push_back(-var);
    maxsat.AddClause(clause);
  }
  assert(block_vars.count(triplets.root));
  maxsat.AddClause({block_vars[triplets.root]});
}


void BtAlgorithm::SetTW() {
  problem_ = Problem::tw;
}

void BtAlgorithm::SetMF() {
  problem_ = Problem::mf;
}

void BtAlgorithm::SetTTS(const vector<int64_t>& bn_domains) {
  problem_ = Problem::tts;
  bn_domains_ = bn_domains;
}

void BtAlgorithm::SetId() {
	problem_ = Problem::id;
}

void BtAlgorithm::SetPhylogeny(const std::vector<int>& phyl_colors) {
  problem_ = Problem::phylogeny;
  phyl_colors_ = phyl_colors;
}
} // namespace triangulator

