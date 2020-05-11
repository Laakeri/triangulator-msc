#include <iostream>
#include <vector>
#include <memory>
#include <fstream>
#include <iomanip>
#include <set>
#include <cassert>

#include "graph.hpp"
#include "hypergraph.hpp"
#include "io.hpp"
#include "enumerator.hpp"
#include "utils.hpp"
#include "mcs.hpp"
#include "treewidth.hpp"
#include "staticset.hpp"
#include "tts.hpp"
#include "min-fill.hpp"
#include "ghtw.hpp"
#include "ftw.hpp"
#include "phylogeny.hpp"
#include "treelength.hpp"

using namespace triangulator;

int main(int argc, char** argv) {
  assert(argc == 2);
  std::string ag(argv[1]);
  // Solving:
  if (ag == "tw") {
    Io io;
    Graph graph = io.ReadGraph(std::cin);
    TW tw;
    auto sol_graph = tw.Treewidth(graph);
    Log::Write(1, "Treewidth: ", mcs::Treewidth(sol_graph));
  } else if (ag == "minfill") {
    Io io;
    Graph graph = io.ReadGraph(std::cin);
    MinFill mf;
    auto sol_graph = mf.MinimumFill(graph);
    Log::Write(1, "Minfill: ", graph.FillEdges(sol_graph).size());
  } else if (ag == "ghtw") {
    Io io;
    HyperGraph hg = io.ReadHyperGraph(std::cin);
    GHTW ghtw;
    auto sol_graph = ghtw.CompGHTW(hg.PrimalGraph(), hg.Edges());
    Log::Write(1, "GHTW: ", ChordGHTW(sol_graph, hg.Edges()));
  } else if (ag == "fhtw") {
    Io io;
    HyperGraph hg = io.ReadHyperGraph(std::cin);
    FHTW fhtw;
    auto sol_graph = fhtw.CompFHTW(hg.PrimalGraph(), hg.Edges());
    Log::Write(1, "FHTW: ", ChordFHTW(sol_graph, hg.Edges()));
  } else if (ag == "tts") {
    Io io;
    auto inst = io.ReadBayesNet(std::cin);
    TTS tts;
    auto sol_graph = tts.TotalTableSize(inst.second, inst.first);
    Log::Write(1, "Treewidth: ", mcs::Treewidth(sol_graph));
    Log::Write(1, "TTS: ", mcs::TotalTableSize(sol_graph, inst.first));
  } else if (ag == "phyl_maxcomp") {
    Io io;
    triangulator::PhylMat pm = io.ReadPhylMat(std::cin);
    std::cout<<CharRemoveMaxSat(pm).size()<<std::endl;
  } else if (ag == "phyl_bin_maxcomp") {
    Io io;
    triangulator::PhylMat pm = io.ReadPhylMat(std::cin);
    std::cout<<BinaryCharRemove(pm).size()<<std::endl;
  } else if (ag == "phyl_perf") {
    Io io;
    triangulator::PhylMat pm = io.ReadPhylMat(std::cin);
    std::cout<<PerfectPhylogeny(pm)<<std::endl;
  } else if (ag == "tl") {
    Io io;
    Graph graph = io.ReadGraph(std::cin);
    Treelength tl;
    int treelength = tl.TreeLength(graph);
    Log::Write(1, "Treelength: ", treelength);
  }

  // Utilities with the BT algorithm:
  if (ag == "count_ms") {
    Io io;
    Graph graph = io.ReadGraph(std::cin);
    std::cout<<graph.AllMinseps().size()<<std::endl;
  } else if (ag == "count_pmc") {
    Io io;
    Graph graph = io.ReadGraph(std::cin);
    Enumerator e;
    std::cout<<e.AdvancedPmcEnum(graph).size()<<std::endl;
  } else if (ag == "bayes_tw") {
    Io io;
    auto inst = io.ReadBayesNet(std::cin);
    TW tw;
    auto sol_graph = tw.Treewidth(inst.second);
    Log::Write(1, "Treewidth: ", mcs::Treewidth(sol_graph));
    Log::Write(1, "TTS: ", mcs::TotalTableSize(sol_graph, inst.first));
  }

  // Converting graph formats:
  if (ag == "topace") {
    Io io;
    Graph graph = io.ReadGraph(std::cin);
    std::cout<<"p tw "<<graph.n()<<" "<<graph.m()<<std::endl;
    for (auto edge : graph.Edges()) {
      std::cout<<edge.first+1<<" "<<edge.second+1<<std::endl;
    }
  } else if (ag == "bayes_topace") {
    Io io;
    auto inst = io.ReadBayesNet(std::cin);
    std::cout<<"p tw "<<inst.second.n()<<" "<<inst.second.m()<<std::endl;
    for (auto edge : inst.second.Edges()) {
      std::cout<<edge.first<<" "<<edge.second<<std::endl;
    }
  } else if (ag == "toprimal") {
    Io io;
    HyperGraph hg = io.ReadHyperGraph(std::cin);
    Graph primal = hg.PrimalGraph();
    assert(hg.m() >= 1);
    std::cout<<"p tw "<<primal.n()<<" "<<primal.m()<<std::endl;
    for (auto edge : primal.Edges()) {
      std::cout<<edge.first<<" "<<edge.second<<std::endl;
    }
  } else if (ag == "tocnf") {
    Io io;
    Graph graph = io.ReadGraph(std::cin);
    std::cout<<"p cnf "<<graph.n()<<" "<<graph.m()<<std::endl;
    for (auto edge : graph.Edges()) {
      std::cout<<edge.first+1<<" "<<edge.second+1<<" 0"<<std::endl;
    }
  } else if (ag == "todimacs") {
    Io io;
    Graph graph = io.ReadGraph(std::cin);
    std::cout<<"p edge "<<graph.n()<<" "<<graph.m()<<std::endl;
    for (auto edge : graph.Edges()) {
      std::cout<<"e "<<edge.first+1<<" "<<edge.second+1<<std::endl;
    }
  } 
}
