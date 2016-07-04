
#include "test_cholesky.h"


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


void work_on(SparseMatrix const & sm, int & n_cliques, int & max_cliques) {

#if __MY_DEBUG__
	std::cout << "Sm is " << sm << std::endl;
#endif

	// compute the Cholesky decomposition of A
	//Eigen::SimplicialLLT< SparseMatrix, Eigen::Lower, Eigen::NaturalOrdering<int> > lltof(pSm);
	//Eigen::SimplicialLLT< SparseMatrix, Eigen::Lower, Eigen::COLAMDOrdering<int> > lltof(sm);
	Eigen::SimplicialLLT< SparseMatrix, Eigen::Lower, Eigen::AMDOrdering<int> > lltof(sm);
	// retrieve factor L  in the decomposition

	Eigen::MatrixXd L = lltof.matrixL();
	// The previous two lines can also be written as "L = A.llt().matrixL()"
#if __MY_DEBUG__
	std::cout << "The Cholesky factor L is" << std::endl << L << std::endl;
	std::cout << "To check this, let us compute L * L.transpose()" << std::endl;
	std::cout << L * L.transpose() << std::endl;
	std::cout << "This should equal the matrix pSm" << std::endl;
#endif
	//std::cout << std::endl;
	int const n(static_cast<int>(sm.cols()));
	SparsityPattern g(n);
	for (int i(0); i < n; ++i) {
		for (int j(0); j < i; ++j) {
			if (std::fabs(L.coeff(i, j)) > 1e-6) {
				g[i].insert(j);
				g[j].insert(i);
#if __MY_DEBUG__
				std::cout << std::setw(6) << i;
				std::cout << std::setw(6) << j;
				std::cout << std::endl;
#endif
			}
		}
	}

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
	Labels labels;
	ItLabels itLabels(n);
	int source = -1;
	for (int i(0); i < n; ++i) {
		int const n_size(static_cast<int>(g[i].size()));
		if (source < 0 && n_size <= 2) {
			itLabels[i] = labels.insert(std::make_pair(n, i));
			source = i;
#if __MY_DEBUG__
			std::cout << "source is " << i << std::endl;
#endif
		}
		else {
			itLabels[i] = labels.insert(std::make_pair(0, i));
		}
		//itLabels[i] = labels.insert(std::make_pair(i == 2 ? n : 0, i));
	}
	std::vector<int> sigma(n, n);

	for (int i(n - 1); i >= 0; --i) {
#if __MY_DEBUG__
		std::cout << "Labels   : " << labels << std::endl;
		std::cout << "sigma    : " << sigma << std::endl;
		print_it_label(std::cout << "itLabels : " << std::endl, itLabels, labels);
		std::cout << std::endl;
#endif
		int v = labels.begin()->second;
#if __MY_DEBUG__
		std::cout << "v is " << v << std::endl;
#endif
		labels.erase(labels.begin());
		itLabels[v] = labels.end();
		sigma[i] = v;
		for (auto neighbor : g[v]) {
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
	std::cout << "sigma    : " << sigma << std::endl;
	print_it_label(std::cout << "itLabels : " << std::endl, itLabels, labels);
	std::cout << std::endl;
#endif
	std::vector<bool> inGraph(n, true);
	std::vector<bool> inClique(n, false);
	std::set<IntSet> cliques;
	// clique decomposition
	for (int i(0); i < n; ++i) {
		// simplicial node
		int const v(sigma[i]);
#if __MY_DEBUG__
		std::cout << "v  : " << v << std::endl;
#endif
		IntSet clique;
		bool new_node(!inClique[v]);
		// adjacent nodes still in the graph forms a clique
#if __MY_DEBUG__
		std::cout << "neighbors : ";
#endif
		for (auto neighbor : g[v]) {
			if (inGraph[neighbor]) {
#if __MY_DEBUG__
				std::cout << neighbor << " ";
#endif
				clique.insert(neighbor);
				if (!inClique[neighbor]) {
					new_node = true;
					//std::cout << neighbor << " was not in clique" << std::endl;
				}
			}
		}
#if __MY_DEBUG__
		std::cout << ", new_node is " << new_node << std::endl;
#endif
		if (!new_node)
			clique.clear();
		if (!clique.empty()) {
			clique.insert(v);
			for (auto u : clique) {
				inClique[u] = true;
			}

			if (!cliques.insert(clique).second) {
#if __MY_DEBUG__
				std::cout << "failure " << clique << std::endl;
#endif
			}
			else {
#if __MY_DEBUG__
				std::cout << "success " << clique << std::endl;
#endif
			}
		}
		inGraph[v].flip();
	}
	size_t max_size(0);
	for (auto c : cliques) {
		max_size = std::max(c.size(), max_size);
#if __MY_DEBUG__+1
		std::cout << c.size() << " | " << c << std::endl;
#endif
	}
	n_cliques = static_cast<int>(cliques.size());
	max_cliques = static_cast<int>(max_size);
}


//int test_cholesky(int argc, char** argv){
//	std::srand(std::time(0));
//	int seed(std::rand());
//	std::cout << "rand " << seed << std::endl;
//	std::srand(seed);
//	std::srand(20354);
//	int m(0);
//	int p(0);
//	int a(0);
//	if (argc >= 2){
//		std::stringstream buffer(argv[1]);
//		buffer >> m;
//	}
//	if (argc >= 3) {
//		std::stringstream buffer(argv[2]);
//		buffer >> p;
//	}
//	if (argc >= 4) {
//		std::stringstream buffer(argv[3]);
//		buffer >> a;
//	}
//	SparseMatrix sm;
//	get_pooling(m, p, a, sm);
//	if(sm.cols()<50)
//	std::cout << "Sm is " << std::endl << sm << std::endl;
//
//	std::cout << std::setw(25) << "#";
//	std::cout << std::setw(25) << "[C_k|";
//	std::cout << std::setw(25) << "max_k card(C_k)";
//	std::cout << std::setw(25) << "Best max_k card(C_k)";
//	std::cout << std::setw(25) << "% reduction";
//	std::cout << std::endl;
//	int n_cliques(sm.cols());
//	int max_cliques(sm.cols());
//	int min_max_cliques(sm.cols());
//	int npermutations(1);
//	std::cout << "sm.cols() " << sm.cols() << std::endl;
//	for (size_t i(0); i < npermutations; ++i){
//		work_on(sm, n_cliques, max_cliques);
//		if(min_max_cliques>max_cliques)
//		{
//			min_max_cliques = std::min(min_max_cliques, max_cliques);
//			std::cout << std::setw(25) << i;
//			std::cout << std::setw(25) << n_cliques;
//			std::cout << std::setw(25) << max_cliques;
//			std::cout << std::setw(25) << min_max_cliques;
//			std::cout << std::setw(25) << std::floor((sm.cols()-min_max_cliques)*1.0/(sm.cols())*100);
//			std::cout << std::endl;
//		}
//
//	}
//	return 0;
//}

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

void build(SparseMatrix & input, IntPairSet & chordalExtension, IntSetPtrSet & output) {
	output.clear();
	//std::cout << input << std::endl;
	// compute the Cholesky decomposition of A
	//Eigen::SimplicialLLT< SparseMatrix, Eigen::Lower, Eigen::NaturalOrdering<int> > lltof(input);
	//Eigen::SimplicialLLT< SparseMatrix, Eigen::Lower, Eigen::COLAMDOrdering<int> > lltof(input);
	Eigen::SimplicialLLT< SparseMatrix, Eigen::Lower, Eigen::AMDOrdering<int> > lltof(input);
	// retrieve factor L  in the decomposition

	Eigen::MatrixXd L = lltof.matrixL();
	auto p = lltof.permutationP();
	IntVector oldOrder(input.cols());
	if (p.cols() > 0) {
		SparseMatrix f(p.cols(), 1);
		Triplets sequence;
		for (int i(0); i < p.cols(); ++i)
			sequence.push_back(Triplet(i, 0, 1.0*i));
		f.setFromTriplets(sequence.begin(), sequence.end());

		std::cout << "p.cols() is " << p.cols() << std::endl;
		std::cout << "p.rows() is " << p.rows() << std::endl;
		std::cout << "f.cols() is " << f.cols() << std::endl;
		std::cout << "f.rows() is " << f.rows() << std::endl;

		Eigen::SparseMatrix<double> pf = p*f;
		//std::cout << "pf.outerSize() is " << pf.outerSize() << std::endl;
		for (int k = 0; k < pf.outerSize(); ++k) {
			for (Eigen::SparseMatrix<double>::InnerIterator it(pf, k); it; ++it)
			{
				//std::cout << "---" << std::endl;
				//it.row();
				//it.col();
				//it.value();
				//it.index();
				int i = static_cast<int>(it.row());
				int j = static_cast<int>(it.value());
				//int const i(static_cast<int>(it.row()));
				//int const j(static_cast<int>(it.value()));
				//std::cout << i << " - " << j<< std::endl;
				//std::cout << it.row()<<" - "<<it.col()<<" : "<<it.value() << std::endl;
				oldOrder[i] = j;
			}
		}
		//std::cout << "end" << std::endl;
	}
	else {
		for (int i(0); i < oldOrder.size(); ++i)
			oldOrder[i] = i;
	
}
#if		__MY_DEBUG__
	std::cout << "input is " << std::endl << input << std::endl;
	std::cout << "permutation is " << std::endl << p.toDenseMatrix() << std::endl;
	std::cout << "p.cols() is " << p.cols() << std::endl;
	std::cout << "newOrder is " << std::endl << p*f << std::endl;
	std::cout << "permuted matrix  is " << std::endl << p*input*p.inverse() << std::endl;
#endif
	// The previous two lines can also be written as "L = A.llt().matrixL()"
#if __MY_DEBUG__
	std::cout << "The Cholesky factor L is" << std::endl << L << std::endl;
	std::cout << "To check this, let us compute L * L.transpose()" << std::endl;
	std::cout << L * L.transpose() << std::endl;
	std::cout << "This should equal the matrix pSm" << std::endl;
#endif
	//std::cout << std::endl;
	int const n(static_cast<int>(input.cols()));
	SparsityPattern g(n);
	chordalExtension.clear();
	for (int i(0); i < n; ++i) {
		for (int j(0); j < i; ++j) {
			if (std::fabs(L.coeff(i, j)) > 1e-6) {
				g[i].insert(j);
				g[j].insert(i);
				chordalExtension.insert({ oldOrder[i], oldOrder[j] });
#if __MY_DEBUG__
				std::cout << std::setw(6) << i;
				std::cout << std::setw(6) << j;
				std::cout << std::endl;
#endif
			}
		}
	}

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
	Labels labels;
	ItLabels itLabels(n);
	int source = -1;
	for (int i(0); i < n; ++i) {
		int const n_size(static_cast<int>(g[i].size()));
		if (source < 0 && n_size <= 2) {
			itLabels[i] = labels.insert(std::make_pair(n, i));
			source = i;
#if __MY_DEBUG__
			std::cout << "source is " << i << std::endl;
#endif
		}
		else {
			itLabels[i] = labels.insert(std::make_pair(0, i));
		}
		//itLabels[i] = labels.insert(std::make_pair(i == 2 ? n : 0, i));
	}
	std::vector<int> sigma(n, n);

	for (int i(n - 1); i >= 0; --i) {
#if __MY_DEBUG__
		std::cout << "Labels   : " << labels << std::endl;
		std::cout << "sigma    : " << sigma << std::endl;
		print_it_label(std::cout << "itLabels : " << std::endl, itLabels, labels);
		std::cout << std::endl;
#endif
		int v = labels.begin()->second;
#if __MY_DEBUG__
		std::cout << "v is " << v << std::endl;
#endif
		labels.erase(labels.begin());
		itLabels[v] = labels.end();
		sigma[i] = v;
		for (auto neighbor : g[v]) {
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
	std::cout << "sigma    : " << sigma << std::endl;
	print_it_label(std::cout << "itLabels : " << std::endl, itLabels, labels);
	std::cout << std::endl;
#endif
	std::vector<bool> inGraph(n, true);
	std::vector<bool> inClique(n, false);
	std::set<IntSet> cliques;
	// clique decomposition
	for (int i(0); i < n; ++i) {
		// simplicial node
		int const v(sigma[i]);
#if __MY_DEBUG__
		std::cout << "v  : " << v << std::endl;
#endif
		IntSet clique;
		bool new_node(!inClique[v]);
		// adjacent nodes still in the graph forms a clique
#if __MY_DEBUG__
		std::cout << "neighbors : ";
#endif
		for (auto neighbor : g[v]) {
			if (inGraph[neighbor]) {
#if __MY_DEBUG__
				std::cout << neighbor << " ";
#endif
				clique.insert(neighbor);
				if (!inClique[neighbor]) {
					new_node = true;
					//std::cout << neighbor << " was not in clique" << std::endl;
				}
			}
		}
#if __MY_DEBUG__
		std::cout << ", new_node is " << new_node << std::endl;
#endif
		if (!new_node)
			clique.clear();
		if (!clique.empty()) {
			clique.insert(v);
			for (auto u : clique) {
				inClique[u] = true;
			}

			if (!cliques.insert(clique).second) {
#if __MY_DEBUG__
				std::cout << "failure " << clique << std::endl;
#endif
			}
			else {
#if __MY_DEBUG__
				std::cout << "success " << clique << std::endl;
#endif
			}
		}
		inGraph[v].flip();
	}
	size_t max_size(0);
	for (auto c : cliques) {
		max_size = std::max(c.size(), max_size);
		IntSetPtr reordered(new IntSet);
		for (auto i : c) {
			reordered->insert(oldOrder[i]);
		}
		output.insert(reordered);
#if __MY_DEBUG__
		//std::cout << c.size() << " | " << c << std::endl;
		std::cout << c.size() << " | " << *reordered << std::endl;
#endif
	}
	for (int i(0); i < n; ++i) {
		if (g[i].empty()) {
			output.insert(IntSetPtr(new IntSet({ oldOrder[i] })));
		}
	}
}