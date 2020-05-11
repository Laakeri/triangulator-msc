#include "enumerator.hpp"

#include <vector>
#include <cassert>
#include <set>
#include <iostream>
#include <unordered_set>
#include <map>
#include <queue>

#include "graph.hpp"
#include "mcs.hpp"
#include "utils.hpp"

namespace triangulator {
void Enumerator::Gogo(int ind, int height, const Graph& graph, const std::vector<std::tuple<Bitset, Bitset>>& s_index, const std::vector<std::vector<int>>& childs, const Bitset& t_sep, Bitset musthave, const std::vector<Bitset>& nohave, std::vector<Bitset>& con, std::vector<Bitset>& tpmcs) {
  gos_++;
  const Bitset& ac = std::get<0>(s_index[ind]);
  const Bitset& sb_sep = std::get<1>(s_index[ind]);
  if (sb_sep.Subsumes(t_sep)) return;
  int chunks = t_sep.Chunks();
  for (int i=0;i<chunks;i++) {
    if (((ac.data_[i] | sb_sep.data_[i]) & t_sep.data_[i]) != t_sep.data_[i]) return;
  }
  for (int i=0;i<chunks;i++){
    if (musthave.data_[i] & nohave[ind].data_[i]) {
      return;
    }
  }
  bool no_pmc = false;
  if (!sb_sep.Subsumes(musthave)) {
    no_pmc = true;
  }
  if (!no_pmc) {
    for (auto im : impl_stack_) {
      if (sb_sep.Get(std::get<0>(im)) && !sb_sep.Get(std::get<1>(im))) {
        no_pmc = true;
        break;
      }
    }
  }
  int impl_stack_add = 0;
  if (!no_pmc){
    crosses_++;
    Bitset pmc = (sb_sep | t_sep);
    Bitset vis = ac;
    for (int i=0;i<chunks;i++){
      vis.data_[i] &= (~t_sep.data_[i]);
    }
    bool fail = false;
    std::vector<int> fcomp;
    fcomp.reserve(graph.n());
    for (int t : t_sep) {
      if (sb_sep.Get(t)) continue;
      bool noco = false;
      for (int j=0;j<chunks;j++){
        if ((pmc.data_[j] & con[t].data_[j]) != pmc.data_[j]) {
          noco = true;
          break;
        }
      }
      if (!noco) continue;
      for (int nt : graph.Neighbors(t)) {
        if (vis.Get(nt)) {
          Bitset go_pmc = pmc;
          graph.Dfs22(nt, go_pmc, vis, fcomp, t_sep);
          dfss_++;
          for (int j=0;j<chunks;j++){
            go_pmc.data_[j] = pmc.data_[j] & (~go_pmc.data_[j]);
          }
          for (int j=0;j<(int)fcomp.size();j++){
            for (int jj=0;jj<chunks;jj++){
              con[fcomp[j]].data_[jj] |= go_pmc.data_[jj];
            }
          }
          fcomp.clear();
        }
      }
      for (int j=0;j<chunks;j++){
        if ((pmc.data_[j] & con[t].data_[j]) != pmc.data_[j]) {
          fail = true;
          if (childs[ind].size() > 0) {
            uint64_t bads = pmc.data_[j] & (~con[t].data_[j]);
            while (bads) {
              int a = j*BITS + __builtin_ctzll(bads);
              bads &= ~-bads;
              if (!sb_sep.Get(a)) {
                musthave.SetTrue(a);
                musthave.SetTrue(t);
              } else {
                if (!musthave.Get(t) && !nohave[ind].Get(a)) {
                  impl_stack_.push_back(std::make_tuple(a, t));
                  impl_stack_add++;
                }
              }
            }
          } else {
            break;
          }
        }
      }
      if (fail && childs[ind].size() < 2) break;
    }
    for (int t : t_sep) {
      con[t].CopyFrom(graph.adj_mat2_[t]);
    }
    if (!fail) {
      if (GoodPMC(pmc, graph)) {
        tpmcs.push_back(pmc);
      }
    }
  }
  for (int ni : childs[ind]) {
    Gogo(ni, height+1, graph, s_index, childs, t_sep, musthave, nohave, con, tpmcs);
  }
  for (int i=0;i<impl_stack_add;i++){
    impl_stack_.pop_back();
  }
}

void Enumerator::Combine(const Graph& graph, int x, 
             const std::vector<Bitset>& s_minseps, 
             const std::vector<Bitset>& t_minseps, 
             std::vector<Bitset>& tpmcs,
             Timer& pre_timer,
             Timer& rec_timer) {
  if (s_minseps.size() == 0) return;
  if (t_minseps.size() == 0) return;
  pre_timer.start();
  std::vector<Bitset> con(graph.n());
  for (int j=0;j<graph.n();j++){
    con[j] = graph.adj_mat2_[j];
  }
  std::vector<std::tuple<Bitset, Bitset>> s_index;
  s_index.reserve(s_minseps.size());
  for (int i = 0; i < (int)s_minseps.size(); i++) {
    s_index.push_back(std::make_tuple(graph.AnotherComp(x, s_minseps[i]), s_minseps[i]));
  }
  auto cmp = [](const std::tuple<Bitset, Bitset>& a, const std::tuple<Bitset, Bitset>& b) {
    return std::get<0>(b)<std::get<0>(a);
  };
  std::sort(s_index.begin(), s_index.end(), cmp);
  std::vector<std::vector<int>> childs(s_index.size());
  std::vector<int> parent(s_index.size());
  std::vector<int> isgood(s_index.size());
  for (int i=0;i<(int)s_index.size();i++){
    parent[i]=-1;
  }
  Bitset dummy(graph.n());
  for (int i=0;i<(int)s_index.size();i++) {
    if (i>0) {
      assert(parent[i]>=0);
      assert(std::get<0>(s_index[parent[i]]).Subsumes(std::get<0>(s_index[i])));
    }
    for (int rm : std::get<1>(s_index[i])) {
      Bitset ac = std::get<0>(s_index[i]);
      for (int nx : graph.Neighbors(rm)) {
        ac.SetFalse(nx);
      }
      for (auto comp : graph.BitComps(ac)) {
        int pos = std::lower_bound(s_index.begin(), s_index.end(), std::make_tuple(comp, dummy), cmp)-s_index.begin();
        if (pos >= (int)s_index.size()) continue;
        assert(pos > i);
        assert(comp.chunks_ == dummy.chunks_);
        if (std::get<0>(s_index[pos]) == comp) {
          parent[pos] = i;
        }
      }
    }
  }
  isgood[0] = true; // in order to avoid akward special case
  for (int i=1;i<(int)s_index.size();i++){
    isgood[i] = GoodMinsep(std::get<1>(s_index[i]), graph);
    if (!isgood[parent[i]]) {
      parent[i] = parent[parent[i]];
    }
    assert(isgood[parent[i]]);
    if (isgood[i]) {
      childs[parent[i]].push_back(i);
    }
  }
  for (int i=1;i<(int)s_index.size();i++){
    assert(isgood[parent[i]]);
    assert(parent[i]>=0);
    assert(std::get<0>(s_index[parent[i]]).Subsumes(std::get<0>(s_index[i])));
  }
  std::vector<Bitset> nohave(s_index.size());
  for (int i=(int)s_index.size()-1;i>=0;i--){
    nohave[i] = (~std::get<1>(s_index[i]));
    for (int ni : childs[i]) {
      nohave[i] &= nohave[ni];
    }
  }
  pre_timer.stop();
  rec_timer.start();
  for (const Bitset& t_sep : t_minseps) {
    if (!GoodMinsep(t_sep, graph)) continue;
    impl_stack_.clear();
    Gogo(0, 0, graph, s_index, childs, t_sep, dummy, nohave, con, tpmcs);
  }
  rec_timer.stop();
  return;
}

std::vector<Bitset> Enumerator::AdvancedPmcEnum(const Graph& graph) {
  if (graph.n() == 0) return {};
  assert(graph.IsConnectedOrIsolated());
  Log::Write(2, "Advanced ", graph.n(), " ", graph.m());
  std::vector<int> order = mcs::Mcs(graph);
  std::reverse(order.begin(), order.end());
  std::vector<int> ord(graph.n());
  for (int i = 0; i < graph.n(); i++) {
    ord[order[i]] = i;
  }
  std::vector<Bitset> pmcs;

  Timer ims_timer;
  ims_timer.start();
  std::vector<Bitset> minseps = graph.AllMinseps();
  ims_timer.stop();
  Log::Write(3, "Enum minseps ", minseps.size());
  Log::Write(3, "IMS timer ", ims_timer.getTime().count());

  Graph t_graph = graph;
  crosses_ = 0;
  dfss_ = 0;
  gos_ = 0;
  uint64_t minsep_sq = (uint64_t)minseps.size()*(uint64_t)minseps.size();

  Timer ms_timer;
  Timer o_timer;
  Timer grow_timer;
  Timer pre_timer;
  Timer rec_timer;

  for (int i = graph.n()-1;i>=0;i--) {
    int x = order[i];
    Graph n_graph = t_graph;
    for (int y = 0; y < graph.n(); y++) {
      if (n_graph.HasEdge(x, y)) {
        n_graph.RemoveEdge(x, y);
      }
    }
    std::vector<Bitset> t_minseps;
    std::vector<Bitset> s_minseps;
    std::vector<Bitset> new_minseps;
    std::vector<Bitset> tpmcs;

    if (i == 0) {
      Bitset npmc(graph.n());
      npmc.SetTrue(x);
      tpmcs.push_back(npmc);
    }

    ms_timer.start();
    Bitset mvis(graph.n());
    for (int j=0;j<t_graph.n();j++){
      if (t_graph.Neighbors(j).size()>0) {
        mvis.SetTrue(j);
      }
    }
    for (Bitset ms : minseps) {
      if (ms.Get(x)) {
        ms.SetFalse(x);
        t_minseps.push_back(ms);
        continue;
      }
      Bitset vis = mvis;
      vis.TurnOff(ms);
      bool is_full = t_graph.IsFull(x, ms, vis);
      if (n_graph.IsMinsep(ms)) {
        new_minseps.push_back(ms);
        if (is_full && GoodMinsep(ms, graph)) {
          bool fail = false;
          ms.SetTrue(x);
          vis.SetFalse(x);
          for (int nx : t_graph.Neighbors(x)) {
            if (vis.Get(nx)) {
              if (t_graph.IsFull2(nx, ms, vis)) {
                fail = true;
                break;
              }
            }
          }
          if (!fail) {
            Bitset npmc = ms;
            if (GoodPMC(npmc, graph)) {
              tpmcs.push_back(npmc);
            }
          }
        }
      } else {
        s_minseps.push_back(ms);
        assert(is_full);
        if (is_full && GoodMinsep(ms, graph)) {
          Bitset npmc = ms;
          npmc.SetTrue(x);
          if (GoodPMC(npmc, graph)) {
            tpmcs.push_back(npmc);
          }
        }
      }
    }
    ms_timer.stop();

    Combine(t_graph, x, s_minseps, t_minseps, tpmcs, pre_timer, rec_timer);

    grow_timer.start();
    for (int j=i+1;j<graph.n();j++){
      for (int nx : graph.Neighbors(order[j])) {
        if (ord[nx] < j) {
          t_graph.AddEdge(nx, order[j]);
        }
      }
      Bitset vis(graph.n());
      for (int jj=0;jj<t_graph.n();jj++){
        if (t_graph.Neighbors(jj).size()>0) {
          vis.SetTrue(jj);
        }
      }
      std::vector<int> fcomp;
      fcomp.reserve(graph.n());
      int chunks = vis.Chunks();
      for (auto& pmc : tpmcs) {
        Bitset tvis = vis;
        for (int jj=0;jj<chunks;jj++){
          tvis.data_[jj] &= (~pmc.data_[jj]);
        }
        if (t_graph.IsFull(order[j], pmc, tvis)) {
          pmc.SetTrue(order[j]);
        }
      }
    }
    grow_timer.stop();

    Log::Write(3, "Enum.. ", i, " ", minseps.size(), " ", tpmcs.size(), " ", t_minseps.size(), " ", s_minseps.size());

    o_timer.start();
    for (const Bitset& pmc : tpmcs) {
      if (GoodPMC(pmc, graph)) {
        pmcs.push_back(pmc);
      }
    }

    t_graph = n_graph;
    minseps = new_minseps;
    for (const Bitset& ms : t_minseps) {
      minseps.push_back(ms);
    }
    utils::SortAndDedup(minseps);
    o_timer.stop();
  }
  utils::SortAndDedup(pmcs);
  for (const auto& pmc : pmcs) {
    assert(GoodPMC(pmc, graph));
  }
  Log::Write(3, "Dfss ", dfss_, "/", crosses_);
  Log::Write(3, "Crosses ", crosses_, "/", gos_);
  Log::Write(3, "Gos ", gos_, "/", minsep_sq);
  Log::Write(3, "Pre timer ", pre_timer.getTime().count());
  Log::Write(3, "Rec timer ", rec_timer.getTime().count());
  Log::Write(3, "MS timer ", ms_timer.getTime().count());
  Log::Write(3, "O timer ", o_timer.getTime().count());
  Log::Write(3, "Grow timer ", grow_timer.getTime().count());
  Log::Write(3, "IMS timer ", ims_timer.getTime().count());
  return pmcs;
}

void Enumerator::SetTWFilter(int tw) {
  filter_ = FilterType::tw;
  twub_ = tw;
}

void Enumerator::SetMFFilter(int mf) {
  filter_ = FilterType::mf;
  mfub_ = mf;
}

void Enumerator::SetTTSFilter(int64_t tts, std::vector<int64_t> bn_domains) {
  filter_ = FilterType::tts;
  ttsub_ = tts;
  bn_domains_ = bn_domains;
}

void Enumerator::SetGHTWFilter(int ghtw) {
  filter_ = FilterType::ghtw;
  ghtwub_ = ghtw;
}

void Enumerator::SetPerfPhyFilter(std::vector<int> phyl_colors) {
  filter_ = FilterType::perf_phy;
  phyl_colors_ = phyl_colors;
}

void Enumerator::SetTLFilter(int tl, std::vector<std::vector<int>> dists) {
  filter_ = FilterType::tl;
  tlub_ = tl;
  tl_dists_ = dists;
}

bool Enumerator::GoodMinsep(const Bitset& minsep, const Graph& graph) const {
  if (filter_ == FilterType::tw) {
    return minsep.Popcount() <= twub_;
  } else if (filter_ == FilterType::mf) {
    return graph.FillSize(minsep) <= mfub_;
  } else if (filter_ == FilterType::tts) {
    int64_t cost = 1;
    for (int v : minsep) {
      cost *= bn_domains_[graph.MapBack(v)];
      if (cost > ttsub_) return false;
    }
    return true;
  } else if (filter_ == FilterType::ghtw) {
    return graph.MaximalIS(minsep) <= ghtwub_;
  } else if (filter_ == FilterType::perf_phy) {
    assert((int)phyl_colors_.size() == graph.n());
    std::set<int> cols;
    for (int v : minsep) {
      if (cols.count(phyl_colors_[v])) return false;
      cols.insert(phyl_colors_[v]);
    }
    return true;
  } else if (filter_ == FilterType::tl) {
    int ttl = utils::MaxInSub(tl_dists_, graph.MapBack(minsep.Elements()));
    return ttl <= tlub_;
  } else {
    assert(0);
  }
}

bool Enumerator::GoodPMC(const Bitset& pmc, const Graph& graph) const {
  if (filter_ == FilterType::tw) {
    return pmc.Popcount() - 1 <= twub_;
  } else if (filter_ == FilterType::mf) {
    return graph.FillSize(pmc) <= mfub_;
  } else if (filter_ == FilterType::tts) {
    int64_t cost = 1;
    for (int v : pmc) {
      cost *= bn_domains_[graph.MapBack(v)];
      if (cost > ttsub_) return false;
    }
    return true;
  } else if (filter_ == FilterType::ghtw) {
    return graph.MaximalIS(pmc) <= ghtwub_;
  } else if (filter_ == FilterType::perf_phy) {
    assert((int)phyl_colors_.size() == graph.n());
    std::set<int> cols;
    for (int v : pmc) {
      if (cols.count(phyl_colors_[v])) return false;
      cols.insert(phyl_colors_[v]);
    }
    return true;
  } else if (filter_ == FilterType::tl) {
    int ttl = utils::MaxInSub(tl_dists_, graph.MapBack(pmc.Elements()));
    return ttl <= tlub_;
  } else {
    assert(0);
  }
}

} // namespace triangulator
