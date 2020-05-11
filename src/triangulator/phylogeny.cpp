#include "phylogeny.hpp"

#include <vector>

#include "graph.hpp"
#include "maxhs_iface.hpp"
#include "phyl_mat.hpp"
#include "enumerator.hpp"
#include "graph.hpp"
#include "mcs.hpp"
#include "utils.hpp"
#include "bt_algorithm.hpp"

using std::vector;
using std::pair;

namespace triangulator {
namespace {
class PhylogenyPreprocessor {
 public:
  // Returns pairs of (graph, colors)
  std::vector<std::pair<Graph, std::vector<int>>> Preprocess(const PhylMat& pm);
 private:
  std::vector<std::pair<Graph, std::vector<int>>> instances_;
  void Preprocess1(const Graph& graph, const PhylMat& opm);
  void AddInstance(const Graph& graph, const PhylMat& opm);
};

vector<pair<Graph, vector<int>>> PhylogenyPreprocessor::Preprocess(const PhylMat& pm) {
  instances_.clear();
  Graph graph = pm.GetGraph();
  Preprocess1(graph, pm);
  return instances_;
}

void PhylogenyPreprocessor::Preprocess1(const Graph& graph, const PhylMat& opm) {
  if (graph.n() == 0) return;
  mcs::McsMOutput minimal_triangulation = mcs::McsM(graph);
  if (minimal_triangulation.fill_edges.size() == 0) {
    // Does not contribute to CR.
    return;
  }
  bool bad_fill = false;
  for (auto e : minimal_triangulation.fill_edges) {
    if (opm.Color(graph.MapBack(e.first)) == opm.Color(graph.MapBack(e.second))) {
      bad_fill = true;
      break;
    }
  }
  if (!bad_fill) {
    // Found perfect fill for this component.
    Log::Write(3, "Eliminated component");
    return;
  }
  vector<Graph> atoms = mcs::Atoms(graph, minimal_triangulation);
  if (atoms.size() == 1) {
    AddInstance(graph, opm);
  } else {
    for (Graph& atom : atoms) {
      atom.InheritMap(graph);
      Preprocess1(atom, opm);
    }
  }
}

void PhylogenyPreprocessor::AddInstance(const Graph& graph, const PhylMat& opm) {
  vector<int> color(graph.n());
  for (int i = 0; i < graph.n(); i++) {
    int og = graph.MapBack(i);
    assert(og >= 0 && og < opm.GetGraph().n());
    color[i] = opm.Color(og);
  }
  instances_.push_back({graph, color});
}
} // namespace

vector<int> CharRemoveMaxSat(const PhylMat& pm) {
  Timer time_phase1;
  time_phase1.start();
  Log::Write(3, "Graph size ", pm.GetGraph().n());
  PhylogenyPreprocessor ppp;
  auto instances = ppp.Preprocess(pm);

  MaxhsInterface maxsat;
  vector<int> color_vars;
  for (int i=0;i<pm.m();i++){
  	color_vars.push_back(maxsat.NewVar());
  	maxsat.AddSoftClause({-color_vars.back()}, 1);
  }

  for (const pair<Graph, vector<int>>& instance : instances) {
  	Log::Write(3, "Inst ", instance.first.n(), " ", instance.first.m());
  	Enumerator pe;
  	auto pmcs = pe.AdvancedPmcEnum(instance.first);
  	Log::Write(3, "PMCS ", pmcs.size());

  	BtAlgorithm bt;
  	bt.BuildColorMaxSAT(instance.first, pmcs, instance.second, color_vars, maxsat);
  }
  time_phase1.stop();
  Log::Write(3, "Phase1 ", time_phase1.getTime().count());
  Timer time_phase2;
  time_phase2.start();
  assert(maxsat.Solve());
  time_phase2.stop();
  Log::Write(3, "Phase2 ", time_phase2.getTime().count());
  vector<int> sol;
  for (int i=0;i<pm.m();i++){
  	if (maxsat.SolutionValue(color_vars[i])) {
  		sol.push_back(i);
  	}
  }
  Log::Write(3, "remc ", sol.size());
  for (int c : sol) {
  	std::cout<<c<<" ";
  }
  std::cout<<std::endl;
  return sol;
}

bool PerfectPhylogeny(const PhylMat& pm) {
  Timer time_phase1;
  time_phase1.start();
  Log::Write(3, "Graph size ", pm.GetGraph().n());
  int deg = pm.GetGraph().Degeneracy();
  Log::Write(3, "degen ", deg, " ", pm.m());
  if (deg >= pm.m()) {
    return false;
  }
  PhylogenyPreprocessor ppp;
  auto instances = ppp.Preprocess(pm);
  auto cmp = [](const pair<Graph, vector<int>>& a, const pair<Graph, vector<int>>& b) {
    return a.second.size() < b.second.size();
  };
  std::sort(instances.begin(), instances.end(), cmp);
  bool perf = true;

  for (const pair<Graph, vector<int>>& instance : instances) {
    Log::Write(3, "Inst ", instance.first.n(), " ", instance.first.m());
    Enumerator pe;
    pe.SetPerfPhyFilter(instance.second);
    auto pmcs = pe.AdvancedPmcEnum(instance.first);
    Log::Write(3, "PMCS ", pmcs.size());

    BtAlgorithm bt;
    bt.SetPhylogeny(instance.second);
    auto sol = bt.Solve(instance.first, pmcs);
    Log::Write(3, "Sol ", sol.first);
    if (sol.first == -1) {
      perf = false;
      break;
    } else {
      assert(sol.first == 0);
    }
  }
  time_phase1.stop();
  Log::Write(3, "Phase1 ", time_phase1.getTime().count());
  return perf;
}

std::vector<int> BinaryCharRemove(const PhylMat& pm) {
	assert(pm.Arity() <= 2);
  Timer time_phase1;
  time_phase1.start();
  Log::Write(3, "Graph size ", pm.GetGraph().n());
  PhylogenyPreprocessor ppp;
  auto instances = ppp.Preprocess(pm);
  std::vector<int> remc;
  for (const pair<Graph, vector<int>>& instance : instances) {
    Log::Write(3, "Inst ", instance.first.n(), " ", instance.first.m());
    Enumerator pe;
    auto pmcs = pe.AdvancedPmcEnum(instance.first);
    Log::Write(3, "PMCS ", pmcs.size());

    BtAlgorithm bt;
    bt.SetPhylogeny(instance.second);
    auto sol = bt.Solve(instance.first, pmcs);
    Log::Write(3, "Sol ", sol.first);
    assert(sol.first >= 0);
    int f = 0;
    for (auto e : sol.second) {
      if (instance.second[e.first] == instance.second[e.second]) {
        remc.push_back(instance.second[e.first]);
        f++;
      }
    }
    Log::Write(3, "f ", f);
    assert(f == sol.first);
  }
  time_phase1.stop();
  Log::Write(3, "Phase1 ", time_phase1.getTime().count());
  std::sort(remc.begin(), remc.end());
  for (int i=1;i<(int)remc.size();i++){
    assert(remc[i] > remc[i-1]);
  }
  Log::Write(3, "remc ", remc.size());
  for (int c : remc) {
    std::cout<<c<<" ";
  }
  std::cout<<std::endl;
  return remc;
}

} // namespace