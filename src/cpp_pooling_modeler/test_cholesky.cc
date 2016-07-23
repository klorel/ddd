
#include "test_cholesky.h"
#include "Timer.h"

std::ostream & operator<<(std::ostream & stream, IntSet const & rhs) {
	for (auto & v : rhs)
		stream << v << " ";
	return stream;
}

std::ostream & operator<<(std::ostream & stream, Labels const & labels) {
	for (auto label : labels) {
		stream << "[" << label.first << ", " << label.second << "] ";
	}
	return stream;
}
std::ostream & operator<<(std::ostream & stream, std::vector<int> const & labels) {
	for (auto label : labels) {
		stream << label << " ";
	}
	return stream;
}
void print_it_label(std::ostream & stream, ItLabels const & itLabels, Labels const & labels) {
	for (size_t i(0); i < itLabels.size(); ++i) {
		if (itLabels[i] != labels.end()) {
			std::cout << "label[" << i << "]->" << itLabels[i]->first << " (" << itLabels[i]->second << ")" << std::endl;
		}
		else {
			std::cout << "label[" << i << "]-> NULL " << std::endl;
		}
	}
}

void build_chordal_extension(SparseMatrix const & input, SparsityPattern & output) {

#if __MY_DEBUG__==2
	std::cout << "Building chrodal extension of "<<std::endl << input << std::endl;
#endif
	Timer timer;
	Eigen::SimplicialLLT< SparseMatrix, Eigen::Lower, Eigen::AMDOrdering<int> > lltof(input);
	//std::cout <<std::left<<std::setw(50)<< "SimplicialLLT : " << timer.elapsed() << std::endl;
	//Eigen::SimplicialLLT< SparseMatrix, Eigen::Lower, Eigen::COLAMDOrdering<int> > lltof(input);
	//timer.restart();
	SparseMatrix L = lltof.matrixL();
	//std::cout << std::left << std::setw(50) << "lltof.matrixL() : " << timer.elapsed() << std::endl;
#if __MY_DEBUG__==2

	std::cout << "The Cholesky factor L is" << std::endl << L << std::endl;
	std::cout << "To check this, let us compute L * L.transpose()" << std::endl;
	std::cout << L * L.transpose() << std::endl;	
#endif
	timer.restart();
	//std::cout << std::endl;
	int const n(static_cast<int>(input.cols()));
	output.assign(n, IntSet());
	size_t newLinks(0);
	for (int k = 0; k < L.outerSize(); ++k) {
		for (SparseMatrix::InnerIterator it(L, k); it; ++it)
		{
			double value = it.value();
			SparseMatrix::Index i = it.row();   // row index
			SparseMatrix::Index j = it.col();   // col index (here it is equal to k)
			SparseMatrix::Index index = it.index(); // inner index, here it is equal to it.row()

			if (std::fabs(value) > 1e-10) {
				output[i].insert((int)j);
				output[j].insert((int)i);
				if (std::fabs(input.coeff(i, j)) > 1e-10) {
					++newLinks;
				}
#if __MY_DEBUG__
				std::cout << std::setw(6) << i;
				std::cout << std::setw(6) << j;
				std::cout << std::endl;
#endif
			}

		}
	}
	std::cout << "chordal extension created " << newLinks << " new links" << std::endl;
//	for (int i(0); i < n; ++i) {
//		for (int j(0); j < i; ++j) {
//			if (std::fabs(L.coeff(i, j)) > 1e-10) {
//				++nz;
//				Timer local_timer;
//				output[i].insert(j);
//				output[j].insert(i);
//				insertTime += local_timer.elapsed();
//#if __MY_DEBUG__
//				std::cout << std::setw(6) << i;
//				std::cout << std::setw(6) << j;
//				std::cout << std::endl;
//#endif
//			}
//		}
//	}

	//std::cout << std::left << std::setw(50) << "chordal extraction : " << timer.elapsed() << std::endl;

}

void build_perfect_elimination_order(SparsityPattern const & input, IntVector & output) {
	//Donnees: Un graphe G = (V, E) et un sommet source s
	//Resultat : Un ordre total σ de V
	//Affecter l’etiquette ∅ a chaque sommet
	//label(s) ←{ n }
	//pour i ← n to 1 faire
	//	Choisir un sommet v d’etiquette maximal pour l’inclusion.
	//	σ(i) ← v
	//	pour chaque sommet non - numerote w ∈ N(v) faire
	//		label(w) ←{ i } ∪ label(w)
	//	fin
	//fin
	int const n((int)input.size());
	Labels labels;
	ItLabels itLabels(n);
	int source = -1;
	for (int i(0); i < n; ++i) {
		int const n_size(static_cast<int>(input[i].size()));
		if (source < 0 && n_size <= 2) {
			itLabels[i] = labels.insert(std::make_pair(n, i));
			source = i;
#if __MY_DEBUG__==2
			std::cout << "source is " << i << std::endl;
#endif
		}
		else {
			itLabels[i] = labels.insert(std::make_pair(0, i));
		}
		//itLabels[i] = labels.insert(std::make_pair(i == 2 ? n : 0, i));
	}

	output.assign(n,n);
	for (int i(n - 1); i >= 0; --i) {
#if __MY_DEBUG__==2
		std::cout << "Labels   : " << labels << std::endl;
		std::cout << "sigma    : " << output << std::endl;
		print_it_label(std::cout << "itLabels : " << std::endl, itLabels, labels);
		std::cout << std::endl;
#endif
		int v = labels.begin()->second;
#if __MY_DEBUG__==2
		std::cout << "v is " << v << std::endl;
#endif
		labels.erase(labels.begin());
		itLabels[v] = labels.end();
		output[i] = v;
		for (auto neighbor : input[v]) {
			if (itLabels[neighbor] != labels.end()) {
				int old_label = itLabels[neighbor]->first;
				labels.erase(itLabels[neighbor]);
				int new_label = old_label + 1;
				itLabels[neighbor] = labels.insert(std::make_pair(new_label, neighbor));
			}
		}
	}
#if __MY_DEBUG__
	std::cout << "Labels   : " << labels << std::endl;
	std::cout << "sigma    : " << output << std::endl;
#endif
#if __MY_DEBUG__==2
	print_it_label(std::cout << "itLabels : " << std::endl, itLabels, labels);
#endif
#if __MY_DEBUG__
	std::cout << std::endl;
#endif
}

void build_clique_decomposition(IntVector const & sigma, SparsityPattern const & chordal_extension, IntSetPtrSet & output) {
	output.clear();
	int const n((int)sigma.size());
	std::vector<bool> inGraph(n, true);
	std::vector<bool> inClique(n, false);
	// clique decomposition
	for (int i(0); i < n; ++i) {
		// simplicial node
		int const v(sigma[i]);
#if __MY_DEBUG__==2
		std::cout << "v  : " << v << std::endl;
#endif
		IntSetPtr cliquePtr(new IntSet);
		IntSet & clique(*cliquePtr);
		bool new_node(!inClique[v]);
		// adjacent nodes still in the graph forms a clique
#if __MY_DEBUG__==2
		std::cout << "neighbors : ";
#endif
		for (auto neighbor : chordal_extension[v]) {
			if (inGraph[neighbor]) {
#if __MY_DEBUG__==2
				std::cout << neighbor << " ";
#endif
				clique.insert(neighbor);
				if (!inClique[neighbor]) {
					new_node = true;
					//std::cout << neighbor << " was not in clique" << std::endl;
				}
			}
		}
#if __MY_DEBUG__==2
		std::cout << ", new_node is " << new_node << std::endl;
#endif
		if (!new_node)
			clique.clear();
		if (!clique.empty() || chordal_extension[v].empty()) {
			clique.insert(v);
			for (auto u : clique) {
				inClique[u] = true;
			}

			if (!output.insert(cliquePtr).second) {
#if __MY_DEBUG__==2
				std::cout << "failure " << clique << std::endl;
#endif
			}
			else {
#if __MY_DEBUG__==2
				std::cout << "success " << clique << std::endl;
#endif
			}
		}
		else {
			// empty clique ? 
		}
		inGraph[v].flip();
	}

}

void work_on(SparseMatrix const & sm, IntSetPtrSet & output) {
	Timer timer;
	SparsityPattern chordal_extension;	
	build_chordal_extension(sm, chordal_extension);
	//std::cout << std::left << std::setw(50) << "build_chordal_extension : "<<timer.elapsed() << std::endl;
	IntVector perfect_elimination_order;
	timer.restart();
	build_perfect_elimination_order(chordal_extension, perfect_elimination_order);
	//std::cout << std::left << std::setw(50) << "build_perfect_elimination_order : " << timer.elapsed() << std::endl;

	timer.restart();
	build_clique_decomposition(perfect_elimination_order, chordal_extension, output);
	//std::cout << std::left << std::setw(50) <<"build_clique_decomposition : " << timer.elapsed() << std::endl;
}



void build(SparsityPattern & input, SparseMatrix & output) {
	Triplets triplets;
	for (int i(0); i < input.size(); ++i) {
		for (auto const & j : input[i]) {
			if (i < j) {
				triplets.push_back({ i, j, 1 });
				triplets.push_back({ j, i, 1 });
			}
		}
		//triplets.push_back({ i, i, 1 + 1000.0*i });
		triplets.push_back({ i, i, 1 + 1.0*input[i].size() });
	}
	int const n(static_cast<int>(input.size()));
	output.resize(n, n);
	output.setZero();
	output.setFromTriplets(triplets.begin(), triplets.end());
}