#pragma once

#include "common.h"


#include <SparseCore>
#include <SparseCholesky>
#include <Eigenvalues>


typedef Eigen::SparseMatrix<double, 0, int> SparseMatrix;
typedef Eigen::Triplet<double> Triplet;
typedef std::vector<Triplet> Triplets;


class TripletPredicate {
public:
	bool operator()(Triplet const & t1, Triplet const & t2)const {
		std::pair<IntPair, double> v1({ {t1.col(), t1.row()}, t1.value() });
		std::pair<IntPair, double> v2({ {t2.col(), t2.row()}, t2.value() });
		return v1 < v2;
	}
};


#define __MY_DEBUG__ 0


std::ostream & operator<<(std::ostream & stream, IntSet const & rhs);
std::ostream & operator<<(std::ostream & stream, Labels const & labels);
std::ostream & operator<<(std::ostream & stream, std::vector<int> const & labels);

void print_it_label(std::ostream & stream, ItLabels const & itLabels, Labels const & labels);


void work_on(SparseMatrix const & sm, IntSetPtrSet &);

void build_chordal_extension(SparseMatrix const & input, SparsityPattern & output);

void build_perfect_elimination_order(SparsityPattern const & input, IntVector & output);

void build_clique_decomposition(IntVector const & sigma, SparsityPattern const & chordal_extension, IntSetPtrSet & output);

void display_info(IntSetPtrSet const &);

size_t get_info(IntSetPtrSet const & cliques, PosInt2PosInt& clique_distribution);