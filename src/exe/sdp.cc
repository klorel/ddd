
#include "common.h"
#include "SdpProblem.h"
#include "SdpSolver.h"
#include "test_cholesky.h"
#include "Timer.h"

void read_graph(std::string const & file_name, SparseMatrix & output, bool complete_graph) {
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
	std::map<int, IntSet> sparsity_pattern;
	if (complete_graph) {
		for (auto const & kvp : network_graph) {
			//triplets.insert({ kvp.first, kvp.first, kvp.second.size()+1});
			// adding power flow balances 
			for (auto i : kvp.second) {
				for (auto j : kvp.second) {
					if (i != j) {
						sparsity_pattern[i].insert(j);
						sparsity_pattern[j].insert(i);
					}
				}
			}
		}
	}
	else
		sparsity_pattern = network_graph;
	std::set<Triplet, TripletPredicate> triplets;
	for (auto const & i : sparsity_pattern) {
		for (auto const j : i.second) {
			triplets.insert({ i.first, j, 1.0 });
			triplets.insert({ j, i.first, 1.0 });
		}
		triplets.insert({ i.first, i.first, 1.0 + i.second.size() });
	}
	size_t t = triplets.size();
	std::cout << "Read    a graph with " << e << " links" << std::endl;
	std::cout << "Created a graph with " << (t - n) / 2 << " links" << std::endl;
	std::cout << "Artificial links     " << ((t - n)/2 - e) << std::endl;
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