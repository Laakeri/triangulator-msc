#pragma once

#include <vector>

#include <ilcplex/ilocplex.h>

#include "bitset.hpp"
#include "graph.hpp"

namespace triangulator {
class HyperCost {
 private:
 	const Graph graph_;
 	const std::vector<std::vector<int>> hyper_edges_;
 	IloEnv env_;
  IloModel model_;
  IloCplex cplex_;
  IloNumVarArray hes_;
  std::vector<IloConstraint> cover_;
  int ub_;
 public:
 	HyperCost(const Graph& graph, const std::vector<std::vector<int>>& hyper_edges, int lb, int ub);
 	int Cover(const Bitset& vs);
 	~HyperCost();
};

class FracHyperCost {
 private:
 	const Graph graph_;
 	const std::vector<std::vector<int>> hyper_edges_;
 	IloEnv env_;
  IloModel model_;
  IloCplex cplex_;
  IloNumVarArray hes_;
  std::vector<IloConstraint> cover_;
  double ub_;
 public:
 	FracHyperCost(const Graph& graph, const std::vector<std::vector<int>>& hyper_edges, double lb, double ub);
 	double Cover(const Bitset& vs);
 	~FracHyperCost();
};
} // namespace triangulator