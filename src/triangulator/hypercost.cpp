#include "hypercost.hpp"

#include <vector>

#include <ilcplex/ilocplex.h>

#include "bitset.hpp"
#include "graph.hpp"


namespace triangulator {
HyperCost::HyperCost(const Graph& graph, const std::vector<std::vector<int>>& hyper_edges, int lb, int ub) :
	graph_(graph), hyper_edges_(hyper_edges), model_(env_), cplex_(model_), hes_(env_), ub_(ub) {
	cplex_.setParam(IloCplex::Threads, 1);
  cplex_.setParam(IloCplex::MIPDisplay, 0);
  cplex_.setOut(env_.getNullStream());
  cplex_.setParam(IloCplex::Param::MIP::Tolerances::UpperCutoff, ub);
  cplex_.setParam(IloCplex::Param::MIP::Tolerances::LowerCutoff, lb);
  int n = 0;
  for (int i=0;i<graph_.n();i++){
  	n = std::max(n, graph_.MapBack(i)+1);
  }
  std::vector<std::vector<int>> covers;
  covers.resize(n);
  cover_.resize(n);
  IloExpr obj(env_);
  for (int i=0;i<(int)hyper_edges_.size();i++){
  	hes_.add(IloNumVar(env_, 0, 1, ILOBOOL));
  	obj += hes_[i];
  	for (int v : hyper_edges[i]) {
  		if (v < n) {
  			covers[v].push_back(i);
  		}
  	}
  }
  model_.add(IloMinimize(env_, obj));
  for (int i=0;i<n;i++){
  	IloExpr vcov(env_);
  	for (int he : covers[i]) {
  		vcov += hes_[he];
  	}
  	cover_[i] = (vcov >= 1);
  }
}

int HyperCost::Cover(const Bitset& vs) {
	for (int v : vs) {
		model_.add(cover_[graph_.MapBack(v)]);
	}
	cplex_.solve();
	int ans = 0;
	if (cplex_.getStatus() == IloAlgorithm::Infeasible) {
		ans = ub_ + 1;
	} else {
		assert(cplex_.getStatus() == IloAlgorithm::Optimal);
	  IloNumArray sol(env_);
	  cplex_.getValues(sol, hes_);
	  for (int i=0;i<(int)hyper_edges_.size();i++){
	  	if (sol[i] > 0.5) {
	  		ans++;
	  	}
	  }
	}
	for (int v : vs) {
		model_.remove(cover_[graph_.MapBack(v)]);
	}
	return ans;
}

HyperCost::~HyperCost() {
	env_.end();
}


FracHyperCost::FracHyperCost(const Graph& graph, const std::vector<std::vector<int>>& hyper_edges, double lb, double ub) :
	graph_(graph), hyper_edges_(hyper_edges), model_(env_), cplex_(model_), hes_(env_), ub_(ub) {
	cplex_.setParam(IloCplex::Threads, 1);
  cplex_.setParam(IloCplex::MIPDisplay, 0);
  cplex_.setParam(IloCplex::SimDisplay, 0);
  cplex_.setOut(env_.getNullStream());
  cplex_.setParam(IloCplex::Param::MIP::Tolerances::UpperCutoff, ub);
  cplex_.setParam(IloCplex::Param::MIP::Tolerances::LowerCutoff, lb);
  int n = 0;
  for (int i=0;i<graph_.n();i++){
  	n = std::max(n, graph_.MapBack(i)+1);
  }
  std::vector<std::vector<int>> covers;
  covers.resize(n);
  cover_.resize(n);
  IloExpr obj(env_);
  for (int i=0;i<(int)hyper_edges_.size();i++){
  	hes_.add(IloNumVar(env_, 0, 1, ILOFLOAT));
  	obj += hes_[i];
  	for (int v : hyper_edges[i]) {
  		if (v < n) {
  			covers[v].push_back(i);
  		}
  	}
  }
  model_.add(IloMinimize(env_, obj));
  for (int i=0;i<n;i++){
  	IloExpr vcov(env_);
  	for (int he : covers[i]) {
  		vcov += hes_[he];
  	}
  	cover_[i] = (vcov >= 1);
  }
}

double FracHyperCost::Cover(const Bitset& vs) {
	for (int v : vs) {
		model_.add(cover_[graph_.MapBack(v)]);
	}
	cplex_.solve();
	double ans = 0;
	if (cplex_.getStatus() == IloAlgorithm::Infeasible) {
		ans = ub_ + 1;
	} else {
		assert(cplex_.getStatus() == IloAlgorithm::Optimal);
	  IloNumArray sol(env_);
	  cplex_.getValues(sol, hes_);
	  for (int i=0;i<(int)hyper_edges_.size();i++){
	  	ans += sol[i];
	  }
	}
	for (int v : vs) {
		model_.remove(cover_[graph_.MapBack(v)]);
	}
	return ans;
}

FracHyperCost::~FracHyperCost() {
	env_.end();
}

} // namespace triangulator