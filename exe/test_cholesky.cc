
#include "test_cholesky.h"


std::mt19937 rng;
std::uniform_int_distribution<std::mt19937::result_type> dist;

void  get_pooling(int M, int P, int A, SparseMatrix & result) {

	Triplets triplets;
	std::map<std::pair<int, int>, int> ymp;
	std::map<std::pair<int, int>, int> zpa;
	int var(0);
	for (int m(0); m < M; ++m) {
		for (int p(0); p < P; ++p, ++var) {
			ymp[{m, p}] = var;
		}
	}
	for (int p(0); p < P; ++p) {
		for (int a(0); a < A; ++a, ++var) {
			zpa[{p, a}] = var;
		}
	}
	std::cout << "var is " << var << " | " << P*(M + A) << std::endl;

	for (int m(0); m < M; ++m) {
		for (int p(0); p < P; ++p) {
			for (int a(0); a < A; ++a) {
				triplets.push_back(Triplet(ymp[{m, p}], zpa[{p, a}], 1));
				triplets.push_back(Triplet(zpa[{p, a}], ymp[{m, p}], 1));
			}
		}
	}
	for (auto const & pa : zpa)
		triplets.push_back(Triplet(pa.second, pa.second, 1 + M));

	for (auto const & mp : ymp)
		triplets.push_back(Triplet(mp.second, mp.second, 1+A));

	result.resize(var, var);
	result.setFromTriplets(triplets.begin(), triplets.end());

}

std::ostream & operator<<(std::ostream & stream, IntSet const & rhs){
	for (auto & v : rhs)
		stream << v << " ";
	return stream;
}

std::ostream & operator<<(std::ostream & stream, Labels const & labels){
	for (auto label : labels){
		stream << "[" << label.first << ", " << label.second << "] ";
	}
	return stream;
}
std::ostream & operator<<(std::ostream & stream, std::vector<int> const & labels){
	for (auto label : labels){
		stream << label << " ";
	}
	return stream;
}
void print_it_label(std::ostream & stream, ItLabels const & itLabels, Labels const & labels){
	for (size_t i(0); i < itLabels.size(); ++i){
		if (itLabels[i] != labels.end()){
			std::cout << "label[" << i << "]->" << itLabels[i]->first << " (" << itLabels[i]->second << ")" << std::endl;
		}
		else{
			std::cout << "label[" << i << "]-> NULL " << std::endl;
		}
	}
}


void permutation(SparseMatrix const & input, SparseMatrix & result){
	int n = input.cols();

	std::vector<int> rowSwap(n);
	for (int i(0); i < n; ++i)
		rowSwap[i] = i;
	std::vector<int> colSwap(n);
	for (int i(0); i < n; ++i)
		colSwap[i] = i;

	Triplets triplets;
	for (int i(0); i < n; ++i){
		size_t r_pos = std::rand() % rowSwap.size();
		size_t c_pos = std::rand() % colSwap.size();

		triplets.push_back(Triplet(rowSwap[r_pos], colSwap[c_pos], 1));

		std::swap(colSwap.back(), colSwap[c_pos]);
		std::swap(rowSwap.back(), rowSwap[r_pos]);

		colSwap.pop_back();
		rowSwap.pop_back();
	}

	result.resize(n, n);
	result.setZero();
	result.setFromTriplets(triplets.begin(), triplets.end());
}


void work_on(SparseMatrix const & sm, int & n_cliques, int & max_cliques){

#if __MY_DEBUG__
	std::cout << "Sm is " << sm << std::endl;
#endif
	SparseMatrix pMatrix;
	permutation(sm, pMatrix);

#if __MY_DEBUG__
	std::cout << "permutation is " <<std::endl<< pMatrix << std::endl;
#endif
	//std::cout << "pSm is " << std::endl << pMatrix*sm*pMatrix.transpose() << std::endl;
	SparseMatrix pSm = pMatrix*sm*pMatrix.transpose();
	//return 0;

#if __MY_DEBUG__
	std::cout << "Matrix pSm " << std::endl << pSm << std::endl;
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
	int const n(sm.cols());
	SparsityPattern g(n);
	for (int i(0); i < n; ++i){
		for (int j(0); j < i; ++j){
			if (std::fabs(L.coeff(i, j))>1e-6){
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
	for (int i(0); i < n; ++i){
		int const n_size(g[i].size());
		if (source < 0 && n_size <= 2){
			itLabels[i] = labels.insert(std::make_pair(n, i));
			source = i;
#if __MY_DEBUG__
			std::cout << "source is " << i << std::endl;
#endif
		}
		else{
			itLabels[i] = labels.insert(std::make_pair(0, i));
		}
		//itLabels[i] = labels.insert(std::make_pair(i == 2 ? n : 0, i));
	}
	std::vector<int> sigma(n, n);

	for (int i(n - 1); i >= 0; --i){
#if __MY_DEBUG__
		std::cout << "Labels   : " << labels << std::endl;
		std::cout << "sigma    : " << sigma << std::endl;
		print_it_label(std::cout << "itLabels : " << std::endl, itLabels, labels);
		std::cout<< std::endl;
#endif
		int v = labels.begin()->second;
#if __MY_DEBUG__
		std::cout << "v is " << v << std::endl;
#endif
		labels.erase(labels.begin());
		itLabels[v] = labels.end();
		sigma[i] = v;
		for (auto neighbor : g[v]){
			if (itLabels[neighbor] != labels.end()){
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
	for (int i(0); i < n; ++i){
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
		for (auto neighbor : g[v]){
			if (inGraph[neighbor]){
#if __MY_DEBUG__
				std::cout << neighbor << " ";
#endif
				clique.insert(neighbor);
				if (!inClique[neighbor]){
					new_node = true;
					//std::cout << neighbor << " was not in clique" << std::endl;
				}
			}
		}
#if __MY_DEBUG__
		std::cout <<", new_node is "<<new_node << std::endl;
#endif
		if (!new_node)
			clique.clear();
		if (!clique.empty()){
			clique.insert(v);
			for (auto u : clique){
				inClique[u] = true;
			}

			if (!cliques.insert(clique).second){
#if __MY_DEBUG__
				std::cout << "failure " << clique << std::endl;
#endif
			}
			else{
#if __MY_DEBUG__
				std::cout << "success " << clique << std::endl;
#endif
			}
		}
		inGraph[v].flip();
	}
	size_t max_size(0);
	for (auto c : cliques){
		max_size = std::max(c.size(), max_size);
#if __MY_DEBUG__+1
		std::cout << c.size() << " | " << c << std::endl;
#endif
	}
	n_cliques = cliques.size();
	max_cliques =  max_size;
}
int test_cholesky(int argc, char** argv){
	std::srand(std::time(0));
	int seed(std::rand());
	std::cout << "rand " << seed << std::endl;
	std::srand(seed);
	std::srand(20354);
	int m(0);
	int p(0);
	int a(0);
	if (argc >= 2){
		std::stringstream buffer(argv[1]);
		buffer >> m;
	}
	if (argc >= 3) {
		std::stringstream buffer(argv[2]);
		buffer >> p;
	}
	if (argc >= 4) {
		std::stringstream buffer(argv[3]);
		buffer >> a;
	}
	SparseMatrix sm;
	get_pooling(m, p, a, sm);
	if(sm.cols()<50)
	std::cout << "Sm is " << std::endl << sm << std::endl;

	std::cout << std::setw(25) << "#";
	std::cout << std::setw(25) << "[C_k|";
	std::cout << std::setw(25) << "max_k card(C_k)";
	std::cout << std::setw(25) << "Best max_k card(C_k)";
	std::cout << std::setw(25) << "% reduction";
	std::cout << std::endl;
	int n_cliques(sm.cols());
	int max_cliques(sm.cols());
	int min_max_cliques(sm.cols());
	int npermutations(1);
	std::cout << "sm.cols() " << sm.cols() << std::endl;
	for (size_t i(0); i < npermutations; ++i){
		work_on(sm, n_cliques, max_cliques);
		if(min_max_cliques>max_cliques)
		{
			min_max_cliques = std::min(min_max_cliques, max_cliques);
			std::cout << std::setw(25) << i;
			std::cout << std::setw(25) << n_cliques;
			std::cout << std::setw(25) << max_cliques;
			std::cout << std::setw(25) << min_max_cliques;
			std::cout << std::setw(25) << std::floor((sm.cols()-min_max_cliques)*1.0/(sm.cols())*100);
			std::cout << std::endl;
		}

	}
	return 0;
}