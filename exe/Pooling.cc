/*
* Pooling.cpp
*
*  Created on: 4 sept. 2012
*      Author: manuel
*/

#include "Pooling.h"
#include "Problem.h"

Pooling::Pooling(size_t n, size_t k) {
	allocate(n, k);
}
void Pooling::allocate(size_t n, size_t k) {
	_n = n;
	_k = k;
	_nPlus.assign(n, IntSet());
	_nMoins.assign(n, IntSet());
	_b.assign(n, NumberVector(2, 0));
	_q.assign(n, NumberVector(k, 0));
	_c.assign(n, NumberVector(n, 0));
}
Pooling::~Pooling() {
}

Number Pooling::c(IntPair const & ij) const{
	return c(ij.first,ij.second);
}
Number Pooling::q(IntPair const &ij) const{
	return q(ij.first, ij.second);
}
Number Pooling::c(int i, int j) const {
	return _c[i][j];
}
Number Pooling::q(int i, int j) const {
	return _q[i][j];
}

Number Pooling::bL(int i) const {
	return _b[i][0];
}
Number Pooling::bU(int i) const {
	return _b[i][1];

}
void Pooling::newArc(int i, int j) {
	_nPlus[i].insert(j);
	_nMoins[j].insert(i);
}

void Pooling::Structure::build(Pooling const & pooling){
	S.clear();
	I.clear();
	T.clear();
	SI.clear();
	ST.clear();
	IT.clear();
	SIT.clear();
	IK.clear();
	TK.clear();
	bL.clear();
	bU.clear();
	for (int n(0); n < pooling._n; ++n) {
		if (pooling._nPlus[n].empty())
			T.insert(n);
		else if (pooling._nMoins[n].empty())
			S.insert(n);
		else
			I.insert(n);
	}

	std::vector<bool> isI(pooling._n, false);
	for (auto const& i : I){
		isI[i].flip();
		for (auto const & s : pooling._nMoins[i])
			SI.insert({ s, i });
		for (auto const &t : pooling._nPlus[i])
			IT.insert({ i, t });
		for (auto const & s : pooling._nMoins[i]){
			for (auto const &t : pooling._nPlus[i]){
				SIT.insert({ { s, i, }, t });
			}
		}
	}
	for (auto const & t : T){
		for (auto const & s : pooling._nMoins[t])
			if (!isI[s])
				ST.insert({ s, t });
	}
	for (int k(0); k < pooling._k; ++k) {
		for (auto const & i : I){
			if (pooling.q(i,k) != 0)
				IK.insert({ i, k });
		}
		for (auto const & t : T){
			if (pooling.q(t, k) != 0)
				TK.insert({ t, k });
		}
	}
	for (int n(0); n < pooling._n; ++n){
		if (pooling.bL(n) != 0){
			bL.insert(n);
		}
		if (pooling.bU(n) != 0){
			bU.insert(n);
		}
	}
}
void Pooling::get(Structure & result) const{
	result.build(*this);
}

void Pooling::read(std::string const & fileName) {
	std::ifstream file(fileName.c_str());
	if (file.good()) {
		int i(-1), j(-1);
		file >> i;
		file >> j;
		std::string header;
		allocate(i, j);
		while (file >> header) {
			if (header == "c") {
				file >> i;
				file >> j;
				file >> _c[i][j];
				newArc(i, j);
			}
			else if (header == "q") {
				file >> i;
				file >> j;
				file >> _q[i][j];
			}
			else if (header == "bu") {
				file >> i;
				file >> _b[i][1];
			}
			else if (header == "bl") {
				file >> i;
				file >> _b[i][0];
			}
		}
	}
	else {
		std::cout << "unable to open " << fileName << "\n";
	}
}
void Pooling::out(std::ostream & stream) const {
	stream << _n << " " << _k << "\n";
	for (size_t i(0); i < _n; ++i) {
		for (auto const & j : _nPlus[i])
			stream << "c " << i << " " << j << " " << _c[i][j] << "\n";
	}
	for (size_t i(0); i < _n; ++i)
		for (size_t k(0); k < _k; ++k)
			//			if (!isZero(_q.get(i, k)))
			stream << "q " << i << " " << k << " " << _q[i][k] << "\n";

	for (size_t i(0); i < _n; ++i)
		if (!isZero(_b[i][0]))
			stream << "bl " << i << " " << _b[i][0] << "\n";
	for (size_t i(0); i < _n; ++i)
		if (!isZero(_b[i][1]))
			stream << "bu " << i << " " << _b[i][1] << "\n";
}
void Pooling::outAmpl(std::ostream & stream) const {
	stream << "data;" << std::endl;
	stream << "param nNodes      := " << _n << ";" << std::endl;
	stream << "param nAttributes := " << _k << ";" << std::endl;
	stream << "param: Edges: c := " << std::endl;
	for (size_t i(0); i < _n; ++i) {
		for (auto const & j : _nPlus[i]) {
			stream << std::setw(4) << i + 1;
			stream << std::setw(4) << j + 1;
			stream << std::setw(25) << _c[i][j];
			stream << std::endl;
		}
	}
	stream << ";" << std::endl;
	stream << "param q := " << std::endl;
	for (size_t i(0); i < _n; ++i) {
		for (size_t k(0); k < _k; ++k) {
			if (!isZero(_q[i][k])) {
				stream << std::setw(4) << i + 1;
				stream << std::setw(4) << k + 1;
				stream << std::setw(25) << _q[i][k];
				stream << std::endl;
			}
		}
	}
	stream << ";" << std::endl;

	stream << "param bL := " << std::endl;
	for (size_t i(0); i < _n; ++i) {
		if (!isZero(_b[i][0])) {
			stream << std::setw(4) << i + 1;
			stream << std::setw(25) << _b[i][0];
			stream << std::endl;
		}
	}
	stream << ";" << std::endl;
	stream << "param bU := " << std::endl;
	for (size_t i(0); i < _n; ++i) {
		if (!isZero(_b[i][1])) {
			stream << std::setw(4) << i + 1;
			stream << std::setw(25) << _b[i][1];
			stream << std::endl;
		}
	}
	stream << ";" << std::endl;

}

size_t Pooling::n() const {
	return _n;
}

size_t Pooling::k() const {
	return _k;
}

IntSet const & Pooling::nPlus(int i) const {
	return _nPlus[i];
}
IntSet const & Pooling::nMoins(int i) const {
	return _nMoins[i];
}


void Pooling::pqFormulation(Problem & result)const{
	Structure structure;
	get(structure);

	result.newvarpool("x", structure.ST);
	result.newvarpool("y", structure.SI);
	result.newvarpool("z", structure.IT);


	// balance
	result.newctrpool("balance", structure.I);
	for (auto const & i : structure.I){
		result.ctr("balance", i).lb() = 1;
		result.ctr("balance", i).ub() = 1;
	}
	for (auto const & si : structure.SI){
		result.ctr("balance", si.second).f()+=result.variable("y", si);
	}
	// attributes
	result.newctrpool("attribute", structure.TK);	
	for (int k(0); k < _k; ++k){
		// <=0
		for (auto const & t : structure.T)
			result.ctr("attribute", t, k).ub() = 0;
		// (q[s,k]-q[t,k])*x[s,t]
		for (auto const & st : structure.ST){
			int const s(st.first);
			int const t(st.second);
			result.ctr("attribute", t, k).f() += (q(s, k) - q(t, k))* result.variable("x", st);
		}
		// (q[s,k]-q[t,k])*y[s,i]*z[i,t]
		for (auto const & sit : structure.SIT){
			int const s(sit.first.first);
			int const i(sit.first.second);
			int const t(sit.second);
			result.ctr("attribute", t, k).f() += (q(s, k) - q(t, k)) * result.variable("y", s, i) * result.variable("z", i, t);
		}
	}
	// bL
	result.newctrpool("bL", structure.bL);
	for (auto const & n : structure.bL){
		result.ctr("bL", n).lb() = bL(n);
	}
	// bU
	result.newctrpool("bU", structure.bU);
	for (auto const & n : structure.bU){
		result.ctr("bU", n).ub() = bU(n);
	}
	// 
	for (auto const & st : structure.ST){
		int const s(st.first);
		int const t(st.second);
		FunctionReal x(result.variable("x", s, t));
		if (bL(s) != 0){
			result.ctr("bL", s).f() += x;
		}
		if (bU(s) != 0){
			result.ctr("bU", s).f() += x;
		}
		if (bL(t) != 0){
			result.ctr("bL", t).f() += x;
		}
		if (bU(t) != 0){
			result.ctr("bU", t).f() += x;
		}
	}
	for (auto const & sit : structure.SIT){
		int const s(sit.first.first);
		int const i(sit.first.second);
		int const t(sit.second);
		FunctionReal y(result.variable("y", s, i));
		FunctionReal z(result.variable("z", i, t));
		if (bL(s) != 0){
			result.ctr("bL", s).f() += y*z;
		}
		if (bU(s) != 0){
			result.ctr("bU", s).f() += y*z;
		}
		if (bL(i) != 0){
			result.ctr("bL", i).f() += y*z;
		}
		if (bU(i) != 0){
			result.ctr("bU", i).f() += y*z;
		}
		if (bL(t) != 0){
			result.ctr("bL", t).f() += y*z;
		}
		if (bU(t) != 0){
			result.ctr("bU", t).f() += y*z;
		}
	}
	// positivity
	result.newctrpool("x_pos", structure.ST);
	for (auto const & st : structure.ST){
		result.ctr("x_pos", st).lb() = 0;
		FunctionReal x(result.variable("x", st));
		result.ctr("x_pos", st).f() += x;
	}
	result.newctrpool("y_pos", structure.SI);
	for (auto const & si : structure.SI){
		result.ctr("y_pos", si).lb() = 0;
		FunctionReal y(result.variable("y", si));
		result.ctr("y_pos", si).f() += y;
	}
	result.newctrpool("z_pos", structure.IT);
	for (auto const & it : structure.IT){
		result.ctr("z_pos", it).lb() = 0;
		FunctionReal y(result.variable("z", it));
		result.ctr("z_pos", it).f() += y;
	}
	// cost

	for (auto const & st : structure.ST){
		int const s(st.first);
		int const t(st.second);
		FunctionReal x(result.variable("x", s, t));
		result.minimize() += c(st)*x;
	}
	for (auto const & sit : structure.SIT){
		int const s(sit.first.first);
		int const i(sit.first.second);
		int const t(sit.second);
		FunctionReal y(result.variable("y", s, i));
		FunctionReal z(result.variable("z", i, t));
		result.minimize() +=(c(s, i)+c(i, t))*y*z;
	}	
}