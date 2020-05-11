#pragma once

#include <vector>
#include <cstdint>

namespace triangulator {

class MaxhsInterface {
 public:
 	const int TrueLit = 1e9;
 	const int FalseLit = -TrueLit;
 	int NewVar();
 	void AddClause(std::vector<int> clause);
 	void AddSoftClause(std::vector<int> clause, int64_t weight);
 	bool SolutionValue(int lit);
 	bool Solve();
 	MaxhsInterface();
 private:
 	int vars_;
 	std::vector<char> solution_value_;
 	std::vector<std::vector<int>> clauses_;
 	std::vector<int64_t> weights_;
};

} // namespace