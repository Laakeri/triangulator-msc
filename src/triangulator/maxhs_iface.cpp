#include "maxhs_iface.hpp"

#include <vector>
#include <cassert>

#include "utils.hpp"

#include "maxhs/core/MaxSolver.h"
#include "maxhs/core/Wcnf.h"

namespace triangulator {
namespace {
	vector<Lit> ConvertClause(const vector<int>& clause, int vars) {
		vector<Lit> ret;
		for (int i=0;i<(int)clause.size();i++){
			assert(clause[i] != 0);
			assert(abs(clause[i]) <= vars);
			ret.push_back(Minisat::mkLit(abs(clause[i]), clause[i] > 0));
		}
		return ret;
	}
	const int64_t top_weight = 1e9;
} // namespace

int MaxhsInterface::NewVar() {
	return ++vars_;
}

void MaxhsInterface::AddClause(std::vector<int> clause) {
	std::vector<int> cls;
	for (int lit : clause) {
		if (lit == TrueLit) {
			return;
		}
		if (lit != FalseLit) {
			cls.push_back(lit);
		}
	}
	clauses_.push_back(cls);
	weights_.push_back(top_weight);
}

void MaxhsInterface::AddSoftClause(std::vector<int> clause, int64_t weight) {
	for (int lit : clause) {
		assert(lit != TrueLit && lit != FalseLit);
	}
	clauses_.push_back(clause);
	weights_.push_back(weight);
}

bool MaxhsInterface::SolutionValue(int lit) {
	int var = abs(lit);
	assert(var > 0);
	assert(var < (int)solution_value_.size());
	if (lit > 0) {
		return solution_value_[var];
	} else {
		return !solution_value_[var];
	}
}

bool MaxhsInterface::Solve() {
	params.readOptions();
  params.verbosity = 0;
  params.mverbosity = 0;
	Wcnf formula;
	formula.set_dimacs_params(vars_, clauses_.size(), top_weight);
	assert(clauses_.size() == weights_.size());
	for (int i=0;i<(int)clauses_.size();i++){
		if (weights_[i] == top_weight) {
			auto mscl = ConvertClause(clauses_[i], vars_);
			formula.addHardClause(mscl);
		} else {
			assert(weights_[i] < top_weight);
			auto mscl = ConvertClause(clauses_[i], vars_);
			formula.addSoftClause(mscl, weights_[i]);
		}
	}
	MaxHS::MaxSolver maxhs(&formula);
	Log::Write(3, "MaxHS solving ", vars_, " ", clauses_.size());
	maxhs.solve();
	assert(maxhs.isSolved());
	assert(!maxhs.isUnsat());
	Log::Write(3, "MaxHS Solved ", maxhs.LB(), " ", maxhs.UB());
	for (auto val : maxhs.getBestModel()) {
		solution_value_.push_back(Minisat::toInt(val));
	}
	assert((int)solution_value_.size() == vars_ + 1);
	return true;
}

MaxhsInterface::MaxhsInterface() : vars_(0) { }

} // namespace