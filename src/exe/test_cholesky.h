#pragma once

#include "common.h"


#define __MY_DEBUG__ 0


std::ostream & operator<<(std::ostream & stream, IntSet const & rhs);
std::ostream & operator<<(std::ostream & stream, Labels const & labels);
std::ostream & operator<<(std::ostream & stream, std::vector<int> const & labels);

void print_it_label(std::ostream & stream, ItLabels const & itLabels, Labels const & labels);

void work_on(SparseMatrix const & sm, int & n_cliques, int & max_cliques);

void build(SparsityPattern & input, SparseMatrix & output);

void build(SparseMatrix & input, IntPairSet & chordalExtension, IntSetPtrSet & output);