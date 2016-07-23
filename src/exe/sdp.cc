
#include "common.h"
#include "SdpProblem.h"
#include "SdpSolver.h"
#include "test_cholesky.h"
#include "Timer.h"

void read_graph(std::string const & file_name, SparseMatrix & output, bool complete_graph) {
	Timer timer;
	std::ifstream file(file_name.c_str());
	int i;
	int j;
	int n(0);
	std::map<int, IntSet> network_graph;
	size_t e(0);

	while (file >> i && file >> j) {
		n = std::max(n, i);
		n = std::max(n, j);
		if (network_graph[i - 1].insert(j - 1).second) {
			++e;
			network_graph[j - 1].insert(i - 1);
		}
	}
	//std::map<int, IntSet> sparsity_pattern;
	std::vector<IntSet> sparsity_pattern(n);
	size_t total_e(e);
	if (complete_graph) {
		for (auto const & kvp : network_graph) {
			//triplets.insert({ kvp.first, kvp.first, kvp.second.size()+1});
			// adding power flow balances 
			for (auto i : kvp.second) {
				for (auto j : kvp.second) {
					if (i < j) {
						++total_e;
						sparsity_pattern[i].insert(j);
						sparsity_pattern[j].insert(i);
					}
				}
			}
		}
	}
	else {
		for (auto const & i : network_graph)
			sparsity_pattern[i.first] = i.second;
	}
	//std::set<Triplet, TripletPredicate> triplets;
	std::vector<Triplet> triplets;
	triplets.reserve(2 * total_e + n);
	for (int i(0); i < n; ++i) {
		for (auto const j : sparsity_pattern[i]) {
			if (i < j) {
				triplets.push_back({ i, j, 1.0 });
				triplets.push_back({ j, i, 1.0 });
			}
		}
		triplets.push_back({ i, i, 1.0 + sparsity_pattern[i].size() });
	}
	size_t t = triplets.size();
	std::cout << "Read    a graph with " << e << " links" << std::endl;
	std::cout << "Created a graph with " << (t - n) / 2 << " links" << std::endl;
	std::cout << "Artificial links     " << ((t - n) / 2 - e) << ", " << timer.elapsed() << " s" << std::endl;
	output = SparseMatrix(n, n);
	output.setFromTriplets(triplets.begin(), triplets.end());
}

int main(int argc, char**argv) {
	SdpProblem sdp1;
	//get_sdp_1(sdp1);
	std::string const file_name(argv[1]);
	std::cout << "reading " << file_name << std::endl;
	bool complete_graph(false);
	if (argc > 2) {
		std::stringstream buffer(argv[2]);
		buffer >> complete_graph;
	}

	//sdp1.read(file_name);
	//IntSetPtrSet output;
	//sdp1.matrix_completion(output);

	SparseMatrix m;
	read_graph(file_name, m, complete_graph);
	Timer timer;
	IntSetPtrSet cliques;
	work_on(m, cliques);
	std::cout << "clique decomposition took  " << timer.elapsed() << std::endl;

	size_t max_size(0);
	std::map<size_t, size_t> clique_distribution;
	for (auto c : cliques) {
		max_size = std::max(c->size(), max_size);
		++clique_distribution[c->size()];

	}
	std::cout << std::setw(8) << "size" << "(" << std::setw(8) << "number" << ")" << std::endl;
	for (auto const & kvp : clique_distribution) {
		std::cout << std::setw(8) << kvp.first << "(" << std::setw(8) << kvp.second << ")" << std::endl;
	}
	return 0;
	//sdp1.print("my_sdp.dat");
	//SdpSolver solver(sdp1);
	////solver.launch_mosek();

	//solver.launch_xpress();
}